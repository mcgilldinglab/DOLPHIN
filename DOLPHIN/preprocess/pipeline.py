import csv
import datetime
import json
import os
import shlex
import subprocess
import sys
import time
from concurrent import futures
from itertools import zip_longest
from typing import Dict, Iterable, List, Optional, Tuple


FEATURECOUNTS_MODE_CONFIG = {
    "gene": {
        "count_column_start": 6,
        "output_suffix": ".exongene.count.txt",
        "drop_zero_rows": False,
    },
    "exon": {
        "count_column_start": 6,
        "output_suffix": ".exon.count.txt",
        "drop_zero_rows": False,
    },
    "junction": {
        "count_column_start": 14,
        "output_suffix": ".exon.count.txt.jcounts",
        "drop_zero_rows": True,
    },
}


TENX_CHEMISTRY_PRESETS = {
    "10xv2": {
        "solo_type": "CB_UMI_Simple",
        "solo_cb_start": 1,
        "solo_cb_len": 16,
        "solo_umi_start": 17,
        "solo_umi_len": 10,
        "solo_barcode_read_length": 0,
        "solo_barcode_mate": 0,
        "read_order": "cdna_barcode",
    },
    "10xv3": {
        "solo_type": "CB_UMI_Simple",
        "solo_cb_start": 1,
        "solo_cb_len": 16,
        "solo_umi_start": 17,
        "solo_umi_len": 12,
        "solo_barcode_read_length": 0,
        "solo_barcode_mate": 0,
        "read_order": "cdna_barcode",
    },
}


def _quote(arg: str) -> str:
    return shlex.quote(str(arg))


def _shell_join(argv: Iterable[str]) -> str:
    return " ".join(_quote(arg) for arg in argv)


def _ensure_directory(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return os.path.abspath(path)


def _output_path_is_complete(path: str) -> bool:
    absolute_path = os.path.abspath(path)
    if not os.path.exists(absolute_path):
        return False
    if os.path.isdir(absolute_path):
        for entry in os.scandir(absolute_path):
            if not entry.is_file():
                continue
            if os.path.getsize(entry.path) > 0:
                return True
        return False
    return os.path.getsize(absolute_path) > 0


def _project_root() -> str:
    return os.path.abspath(
        os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)
    )


def _read_meminfo_kb(field_name: str) -> Optional[int]:
    meminfo_path = "/proc/meminfo"
    if not os.path.isfile(meminfo_path):
        return None

    with open(meminfo_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(field_name):
                continue
            parts = line.split()
            if len(parts) >= 2 and parts[1].isdigit():
                return int(parts[1])
    return None


def get_system_resource_snapshot() -> Dict[str, object]:
    cpu_count = os.cpu_count() or 1
    mem_total_kb = _read_meminfo_kb("MemTotal:")
    mem_available_kb = _read_meminfo_kb("MemAvailable:")
    load_average = None
    if hasattr(os, "getloadavg"):
        try:
            load_average = list(os.getloadavg())
        except OSError:
            load_average = None

    return {
        "cpu_count": int(cpu_count),
        "mem_total_kb": mem_total_kb,
        "mem_available_kb": mem_available_kb,
        "mem_total_gb": (
            round(mem_total_kb / (1024.0 * 1024.0), 2)
            if mem_total_kb is not None else None
        ),
        "mem_available_gb": (
            round(mem_available_kb / (1024.0 * 1024.0), 2)
            if mem_available_kb is not None else None
        ),
        "load_average": load_average,
    }


def recommend_full_length_parallelism(
    sample_count: int,
    star_threads: int,
    featurecounts_threads: int,
) -> Dict[str, object]:
    snapshot = get_system_resource_snapshot()
    cpu_count = max(1, int(snapshot["cpu_count"]))
    cpu_reserve = max(4, cpu_count // 10)
    usable_cpu = max(1, cpu_count - cpu_reserve)
    cpu_limited_parallel_cells = max(1, usable_cpu // max(1, int(star_threads)))

    mem_available_gb = snapshot.get("mem_available_gb")
    estimated_memory_per_cell_gb = round(max(8.0, 4.0 + 1.25 * int(star_threads)), 2)
    memory_reserve_gb = None
    if mem_available_gb is None:
        memory_limited_parallel_cells = cpu_limited_parallel_cells
    else:
        memory_reserve_gb = round(max(16.0, mem_available_gb * 0.10), 2)
        memory_budget_gb = max(
            estimated_memory_per_cell_gb,
            mem_available_gb - memory_reserve_gb,
        )
        memory_limited_parallel_cells = max(
            1,
            int(memory_budget_gb // estimated_memory_per_cell_gb),
        )

    safe_max_parallel_cells = max(
        1,
        min(
            int(sample_count),
            cpu_limited_parallel_cells,
            memory_limited_parallel_cells,
        ),
    )
    conservative_parallel_cells = max(
        1,
        min(safe_max_parallel_cells, max(1, safe_max_parallel_cells // 3)),
    )
    balanced_parallel_cells = max(
        conservative_parallel_cells,
        min(
            safe_max_parallel_cells,
            max(1, round(safe_max_parallel_cells * 0.60)),
        ),
    )
    aggressive_parallel_cells = safe_max_parallel_cells

    return {
        "system_snapshot": snapshot,
        "star_threads_per_cell": int(star_threads),
        "featurecounts_threads_per_cell": int(featurecounts_threads),
        "cpu_reserve": int(cpu_reserve),
        "usable_cpu": int(usable_cpu),
        "estimated_memory_per_cell_gb": estimated_memory_per_cell_gb,
        "memory_reserve_gb": memory_reserve_gb,
        "cpu_limited_parallel_cells": int(cpu_limited_parallel_cells),
        "memory_limited_parallel_cells": int(memory_limited_parallel_cells),
        "safe_max_parallel_cells": int(safe_max_parallel_cells),
        "default_profile": "balanced",
        "profiles": {
            "conservative": {
                "max_parallel_cells": int(conservative_parallel_cells),
                "description": (
                    "Lower risk. Leaves extra headroom for other workloads and "
                    "filesystem pressure."
                ),
            },
            "balanced": {
                "max_parallel_cells": int(balanced_parallel_cells),
                "description": (
                    "Recommended default. Uses a moderate fraction of the "
                    "available CPU budget for STAR-heavy cells."
                ),
            },
            "aggressive": {
                "max_parallel_cells": int(aggressive_parallel_cells),
                "description": (
                    "Highest throughput. Best when this machine is dedicated to "
                    "the preprocessing job."
                ),
            },
        },
        "rationale": (
            "Recommendations are based on CPU core count, currently available "
            "memory, and the fact that full-length preprocessing is dominated "
            "by the STAR alignment stage."
        ),
    }


def _infer_threads_from_stage_command(
    stages: List[Dict[str, object]],
    stage_prefix: str,
    flag: str,
    default: int,
) -> int:
    for stage in stages:
        stage_name = str(stage.get("stage_name", ""))
        if not stage_name.startswith(stage_prefix):
            continue
        try:
            argv = shlex.split(str(stage.get("command", "")))
        except ValueError:
            continue
        if flag not in argv:
            continue
        flag_index = argv.index(flag)
        if flag_index + 1 >= len(argv):
            continue
        try:
            return int(argv[flag_index + 1])
        except ValueError:
            continue
    return int(default)


def get_full_length_parallel_recommendations_from_plan(
    plan: Dict[str, object],
) -> Dict[str, object]:
    existing_recommendations = plan.get("parallel_recommendations")
    if existing_recommendations:
        return existing_recommendations

    stages = plan.get("stages", [])
    sample_count = int(
        plan.get("sample_count")
        or len(plan.get("cell_ids") or [])
        or sum(
            1
            for stage in stages
            if str(stage.get("stage_name", "")).startswith("align_exon.")
        )
    )
    star_threads = int(
        plan.get("star_threads_per_cell")
        or _infer_threads_from_stage_command(
            stages=stages,
            stage_prefix="align_exon.",
            flag="--runThreadN",
            default=4,
        )
    )
    featurecounts_threads = int(
        plan.get("featurecounts_threads_per_cell")
        or _infer_threads_from_stage_command(
            stages=stages,
            stage_prefix="count_gene.",
            flag="-T",
            default=4,
        )
    )
    return recommend_full_length_parallelism(
        sample_count=sample_count,
        star_threads=star_threads,
        featurecounts_threads=featurecounts_threads,
    )


def _guess_separator(path: str) -> str:
    with open(path, "r", encoding="utf-8") as handle:
        sample = handle.read(4096)

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        return dialect.delimiter
    except csv.Error:
        lower_path = path.lower()
        if lower_path.endswith(".tsv") or lower_path.endswith(".txt"):
            return "\t"
        if lower_path.endswith(".csv"):
            return ","
        return "\t" if "\t" in sample else ","


def read_metadata_table(metadata_path: str) -> Dict[str, object]:
    metadata_path = os.path.abspath(metadata_path)
    if not os.path.isfile(metadata_path):
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

    sep = _guess_separator(metadata_path)
    with open(metadata_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=sep)
        fieldnames = reader.fieldnames or []
        if "CB" not in fieldnames:
            raise ValueError(
                f"Metadata file must contain a 'CB' column: {metadata_path}"
            )

        rows = []
        cell_ids = []
        for row in reader:
            cb_value = row.get("CB")
            if cb_value is None or str(cb_value).strip() == "":
                raise ValueError(
                    f"Metadata file contains empty CB values: {metadata_path}"
                )
            normalized_row = dict(row)
            normalized_row["CB"] = str(cb_value).strip()
            rows.append(normalized_row)
            cell_ids.append(normalized_row["CB"])

    if not rows:
        raise ValueError(
            f"Metadata file does not contain any data rows: {metadata_path}"
        )
    return {
        "path": metadata_path,
        "columns": fieldnames,
        "rows": rows,
        "cell_ids": cell_ids,
    }


def load_reference_manifest(reference_manifest_path: str) -> Dict[str, object]:
    reference_manifest_path = os.path.abspath(reference_manifest_path)
    if not os.path.isfile(reference_manifest_path):
        raise FileNotFoundError(
            f"Reference manifest not found: {reference_manifest_path}"
        )

    with open(reference_manifest_path, "r", encoding="utf-8") as handle:
        manifest = json.load(handle)

    required_keys = [
        "reference_dir",
        "alignment_gtf_path",
        "adj_index_path",
        "adj_metadata_path",
        "exon_pkl_path",
    ]
    missing_keys = [key for key in required_keys if key not in manifest]
    if missing_keys:
        raise ValueError(
            f"Reference manifest is missing required keys: {missing_keys}"
        )
    return manifest


def resolve_star_index_dir(
    reference_manifest: Dict[str, object],
    explicit_star_index_dir: Optional[str] = None,
) -> Tuple[str, List[str]]:
    warnings = []
    if explicit_star_index_dir:
        star_index_dir = os.path.abspath(explicit_star_index_dir)
    else:
        star_index = reference_manifest.get("star_index") or {}
        star_index_dir = star_index.get("star_index_dir")
        if star_index_dir is None:
            star_index_dir = os.path.join(
                os.path.abspath(reference_manifest["reference_dir"]),
                "star_index",
            )
            warnings.append(
                "Reference manifest does not contain a STAR index entry. "
                "Using the default path under the reference directory."
            )

    if not os.path.isdir(star_index_dir):
        warnings.append(
            f"STAR index directory does not exist yet: {star_index_dir}"
        )
    return star_index_dir, warnings


def discover_full_length_fastqs(
    cell_ids: List[str],
    raw_fastq_dir: str,
    read_layout: str = "single",
    fastq_suffix: str = ".fastq.gz",
    skip_missing_fastqs: bool = True,
) -> Tuple[List[Dict[str, object]], List[Dict[str, str]]]:
    raw_fastq_dir = os.path.abspath(raw_fastq_dir)
    samples = []
    missing = []

    for cb in cell_ids:
        if read_layout == "single":
            read1_path = os.path.join(raw_fastq_dir, f"{cb}{fastq_suffix}")
            if not os.path.isfile(read1_path):
                missing.append({"cb": cb, "expected_path": read1_path})
                if skip_missing_fastqs:
                    continue
                raise FileNotFoundError(f"FASTQ file not found: {read1_path}")

            samples.append(
                {
                    "cell_id": cb,
                    "read_layout": "single",
                    "read1_path": read1_path,
                    "read2_path": None,
                }
            )
        else:
            read1_path = os.path.join(raw_fastq_dir, f"{cb}_1{fastq_suffix}")
            read2_path = os.path.join(raw_fastq_dir, f"{cb}_2{fastq_suffix}")
            if not os.path.isfile(read1_path) or not os.path.isfile(read2_path):
                missing.extend(
                    [
                        {"cb": cb, "expected_path": read1_path},
                        {"cb": cb, "expected_path": read2_path},
                    ]
                )
                if skip_missing_fastqs:
                    continue
                raise FileNotFoundError(
                    f"Paired FASTQ files not found: {read1_path}, {read2_path}"
                )

            samples.append(
                {
                    "cell_id": cb,
                    "read_layout": "paired",
                    "read1_path": read1_path,
                    "read2_path": read2_path,
                }
            )

    return samples, missing


def discover_full_length_fastqs_from_dir(
    raw_fastq_dir: str,
    read_layout: str = "single",
    fastq_suffix: str = ".fastq.gz",
) -> List[Dict[str, object]]:
    raw_fastq_dir = os.path.abspath(raw_fastq_dir)
    fastq_names = sorted(
        name for name in os.listdir(raw_fastq_dir)
        if name.endswith(fastq_suffix)
    )

    samples = []
    if read_layout == "single":
        for fastq_name in fastq_names:
            cell_id = fastq_name[: -len(fastq_suffix)]
            samples.append(
                {
                    "cell_id": cell_id,
                    "read_layout": "single",
                    "read1_path": os.path.join(raw_fastq_dir, fastq_name),
                    "read2_path": None,
                }
            )
        return samples

    paired_prefixes = []
    for fastq_name in fastq_names:
        if fastq_name.endswith(f"_1{fastq_suffix}"):
            prefix = fastq_name[: -len(f"_1{fastq_suffix}")]
            paired_prefixes.append(prefix)

    paired_prefixes = sorted(set(paired_prefixes))
    for prefix in paired_prefixes:
        read1_path = os.path.join(raw_fastq_dir, f"{prefix}_1{fastq_suffix}")
        read2_path = os.path.join(raw_fastq_dir, f"{prefix}_2{fastq_suffix}")
        if not os.path.isfile(read2_path):
            raise FileNotFoundError(
                f"Paired FASTQ files not found: {read1_path}, {read2_path}"
            )
        samples.append(
            {
                "cell_id": prefix,
                "read_layout": "paired",
                "read1_path": read1_path,
                "read2_path": read2_path,
            }
        )

    return samples


def discover_tenx_fastq_pair(
    raw_fastq_dir: str,
    sample_name: Optional[str] = None,
) -> Dict[str, str]:
    raw_fastq_dir = os.path.abspath(raw_fastq_dir)
    fastq_names = sorted(os.listdir(raw_fastq_dir))

    read1_candidates = [
        name for name in fastq_names
        if name.endswith(".fastq.gz") and "_1.fastq.gz" in name
    ]

    pairs = []
    for read1_name in read1_candidates:
        prefix = read1_name.rsplit("_1.fastq.gz", 1)[0]
        read2_name = f"{prefix}_2.fastq.gz"
        read2_path = os.path.join(raw_fastq_dir, read2_name)
        if os.path.isfile(read2_path):
            pairs.append(
                {
                    "sample_name": prefix,
                    "read1_path": os.path.join(raw_fastq_dir, read1_name),
                    "read2_path": read2_path,
                }
            )

    if sample_name is not None:
        for pair in pairs:
            if pair["sample_name"] == sample_name:
                return pair
        raise FileNotFoundError(
            f"No 10x FASTQ pair matched sample name '{sample_name}' in {raw_fastq_dir}"
        )

    if len(pairs) != 1:
        raise ValueError(
            "Could not infer a unique 10x FASTQ pair. "
            "Please provide --sample-name explicitly."
        )
    return pairs[0]


def infer_common_cb_suffix(cell_ids: Iterable[str]) -> Optional[str]:
    suffixes = []
    for cell_id in cell_ids:
        if "_" not in cell_id:
            return None
        suffixes.append(cell_id[cell_id.index("_"):])

    unique_suffixes = sorted(set(suffixes))
    if len(unique_suffixes) == 1:
        return unique_suffixes[0]
    return None


def collect_cb_suffix_counts(cell_ids: Iterable[str]) -> Dict[str, int]:
    suffix_counts: Dict[str, int] = {}
    for cell_id in cell_ids:
        if "_" not in cell_id:
            continue
        suffix = cell_id[cell_id.index("_"):]
        suffix_counts[suffix] = suffix_counts.get(suffix, 0) + 1
    return suffix_counts


def filter_metadata_by_cb_suffix(
    metadata_table: Dict[str, object],
    cb_suffix_filter: Optional[str] = None,
) -> Dict[str, object]:
    if cb_suffix_filter is None:
        return metadata_table

    filtered_rows = [
        row for row in metadata_table["rows"]
        if row["CB"].endswith(cb_suffix_filter)
    ]
    if not filtered_rows:
        raise ValueError(
            f"No metadata rows matched CB suffix filter: {cb_suffix_filter}"
        )

    return {
        "path": metadata_table["path"],
        "columns": metadata_table["columns"],
        "rows": filtered_rows,
        "cell_ids": [row["CB"] for row in filtered_rows],
    }


def _plan_record(
    stage_name: str,
    description: str,
    command: str,
    outputs: List[str],
    mode: str,
    sample_id: Optional[str] = None,
) -> Dict[str, object]:
    record = {
        "stage_name": stage_name,
        "description": description,
        "mode": mode,
        "command": command,
        "outputs": outputs,
    }
    if sample_id is not None:
        record["sample_id"] = sample_id
    return record


def read_cell_id_list(cell_id_list_path: str) -> List[str]:
    cell_id_list_path = os.path.abspath(cell_id_list_path)
    if not os.path.isfile(cell_id_list_path):
        raise FileNotFoundError(f"Cell ID list not found: {cell_id_list_path}")

    cell_ids = []
    with open(cell_id_list_path, "r", encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if not value:
                continue
            cell_ids.append(value.split("\t", 1)[0])

    if not cell_ids:
        raise ValueError(f"Cell ID list is empty: {cell_id_list_path}")
    return cell_ids


def write_cell_id_list(cell_ids: Iterable[str], output_path: str) -> str:
    output_path = os.path.abspath(output_path)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    unique_cell_ids = []
    seen = set()
    for cell_id in cell_ids:
        normalized = str(cell_id).strip()
        if not normalized or normalized in seen:
            continue
        seen.add(normalized)
        unique_cell_ids.append(normalized)

    if not unique_cell_ids:
        raise ValueError("Cannot write an empty cell ID list.")

    with open(output_path, "w", encoding="utf-8") as handle:
        for cell_id in unique_cell_ids:
            handle.write(f"{cell_id}\n")
    return output_path


def write_cell_manifest(cell_ids: Iterable[str], output_path: str) -> str:
    output_path = os.path.abspath(output_path)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    unique_cell_ids = []
    seen = set()
    for cell_id in cell_ids:
        normalized = str(cell_id).strip()
        if not normalized or normalized in seen:
            continue
        seen.add(normalized)
        unique_cell_ids.append(normalized)

    if not unique_cell_ids:
        raise ValueError("Cannot write an empty cell manifest.")

    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["CB"])
        for cell_id in unique_cell_ids:
            writer.writerow([cell_id])
    return output_path


def _resolve_featurecounts_input_files(
    input_dir: str,
    mode: str,
    requested_cell_ids: Optional[List[str]] = None,
) -> List[Tuple[str, str]]:
    if mode not in FEATURECOUNTS_MODE_CONFIG:
        raise ValueError(
            f"Unsupported featureCounts mode: {mode}. "
            f"Supported modes: {sorted(FEATURECOUNTS_MODE_CONFIG)}"
        )

    input_dir = os.path.abspath(input_dir)
    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f"featureCounts directory not found: {input_dir}")

    suffix = FEATURECOUNTS_MODE_CONFIG[mode]["output_suffix"]
    file_by_cell_id = {}
    for entry in os.scandir(input_dir):
        if not entry.is_file():
            continue
        if not entry.name.endswith(suffix):
            continue
        cell_id = entry.name[: -len(suffix)]
        file_by_cell_id[cell_id] = entry.path

    if requested_cell_ids is None:
        requested_cell_ids = sorted(file_by_cell_id)

    missing_cells = [
        cell_id for cell_id in requested_cell_ids if cell_id not in file_by_cell_id
    ]
    if missing_cells:
        raise FileNotFoundError(
            "featureCounts directory is missing requested cells. "
            f"Examples: {missing_cells[:10]}"
        )

    return [
        (cell_id, file_by_cell_id[cell_id])
        for cell_id in requested_cell_ids
    ]


def _read_featurecounts_header(handle, mode: str):
    comment_lines = []
    header_line = None
    if mode in {"gene", "exon"}:
        while True:
            line = handle.readline()
            if line == "":
                break
            if line.startswith("#"):
                comment_lines.append(line)
                continue
            header_line = line
            break
    else:
        header_line = handle.readline()

    if not header_line:
        raise ValueError("featureCounts file is missing its header line.")
    header = header_line.rstrip("\n").split("\t")
    return comment_lines, header


def _combine_dense_featurecounts_files(
    file_records: List[Tuple[str, str]],
    output_path: str,
    mode: str,
) -> Dict[str, object]:
    config = FEATURECOUNTS_MODE_CONFIG[mode]
    handles = []
    first_comment_lines = []
    first_annotation_columns = None
    row_count = 0

    try:
        for index, (cell_id, file_path) in enumerate(file_records):
            handle = open(file_path, "r", encoding="utf-8")
            comment_lines, header = _read_featurecounts_header(handle, mode=mode)
            annotation_columns = header[: config["count_column_start"]]
            if first_annotation_columns is None:
                first_annotation_columns = annotation_columns
                first_comment_lines = comment_lines
            elif annotation_columns != first_annotation_columns:
                raise ValueError(
                    f"Annotation columns differ in {file_path}. "
                    "Cannot combine full-length featureCounts tables safely."
                )
            handles.append((cell_id, file_path, handle))

        with open(output_path, "w", encoding="utf-8", newline="") as output_handle:
            for line in first_comment_lines:
                output_handle.write(line)
            output_handle.write(
                "\t".join(first_annotation_columns + [cell_id for cell_id, _, _ in handles])
                + "\n"
            )

            while True:
                first_line = handles[0][2].readline()
                if first_line == "":
                    break
                first_line = first_line.rstrip("\n")
                first_prefix, first_count = first_line.rsplit("\t", 1)
                counts = [first_count]

                for cell_id, file_path, handle in handles[1:]:
                    line = handle.readline()
                    if line == "":
                        raise ValueError(
                            "featureCounts tables ended at different row counts. "
                            f"First mismatch encountered in {file_path}"
                        )
                    prefix, count = line.rstrip("\n").rsplit("\t", 1)
                    if prefix != first_prefix:
                        raise ValueError(
                            "featureCounts annotation rows differ across cells. "
                            f"First mismatch encountered in {file_path}"
                        )
                    counts.append(count)

                output_handle.write(first_prefix + "\t" + "\t".join(counts) + "\n")
                row_count += 1

            for _, file_path, handle in handles[1:]:
                extra_line = handle.readline()
                if extra_line != "":
                    raise ValueError(
                        "featureCounts tables ended at different row counts. "
                        f"Extra rows found in {file_path}"
                    )
    finally:
        for _, _, handle in handles:
            handle.close()

    return {
        "mode": mode,
        "output_path": os.path.abspath(output_path),
        "cell_count": len(file_records),
        "row_count": row_count,
        "cell_ids": [cell_id for cell_id, _ in file_records],
        "source_files": [file_path for _, file_path in file_records],
    }


def _combine_sparse_junction_featurecounts_files(
    file_records: List[Tuple[str, str]],
    output_path: str,
) -> Dict[str, object]:
    config = FEATURECOUNTS_MODE_CONFIG["junction"]
    annotation_columns = None
    row_index_by_key = {}
    row_keys: List[Tuple[str, ...]] = []
    row_sparse_counts: List[Dict[int, str]] = []

    for cell_index, (cell_id, file_path) in enumerate(file_records):
        with open(file_path, "r", encoding="utf-8") as handle:
            _, header = _read_featurecounts_header(handle, mode="junction")
            current_annotation_columns = header[: config["count_column_start"]]
            if annotation_columns is None:
                annotation_columns = current_annotation_columns
            elif current_annotation_columns != annotation_columns:
                raise ValueError(
                    f"Junction annotation columns differ in {file_path}. "
                    "Cannot combine full-length junction tables safely."
                )

            for line in handle:
                stripped = line.rstrip("\n")
                if not stripped:
                    continue
                row = stripped.split("\t")
                if len(row) < config["count_column_start"] + 1:
                    raise ValueError(
                        f"Malformed junction featureCounts row in {file_path}: {stripped[:200]}"
                    )
                count_value = row[-1]
                try:
                    if float(count_value) <= 0:
                        continue
                except ValueError:
                    pass
                key = tuple(row[: config["count_column_start"]])
                row_index = row_index_by_key.get(key)
                if row_index is None:
                    row_index = len(row_keys)
                    row_index_by_key[key] = row_index
                    row_keys.append(key)
                    row_sparse_counts.append({})
                row_sparse_counts[row_index][cell_index] = count_value

    if annotation_columns is None:
        raise ValueError("No junction featureCounts files were combined.")

    with open(output_path, "w", encoding="utf-8", newline="") as output_handle:
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(annotation_columns + [cell_id for cell_id, _ in file_records])
        for key, sparse_counts in zip(row_keys, row_sparse_counts):
            counts = ["0"] * len(file_records)
            for cell_index, count_value in sparse_counts.items():
                counts[cell_index] = count_value
            writer.writerow(list(key) + counts)

    return {
        "mode": "junction",
        "output_path": os.path.abspath(output_path),
        "cell_count": len(file_records),
        "row_count": len(row_keys),
        "cell_ids": [cell_id for cell_id, _ in file_records],
        "source_files": [file_path for _, file_path in file_records],
    }


def combine_featurecounts_directory(
    input_dir: str,
    output_path: str,
    mode: str,
    cell_id_list_path: Optional[str] = None,
    summary_path: Optional[str] = None,
) -> Dict[str, object]:
    requested_cell_ids = (
        read_cell_id_list(cell_id_list_path)
        if cell_id_list_path is not None else None
    )
    file_records = _resolve_featurecounts_input_files(
        input_dir=input_dir,
        mode=mode,
        requested_cell_ids=requested_cell_ids,
    )
    output_path = os.path.abspath(output_path)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    if mode in {"gene", "exon"}:
        summary = _combine_dense_featurecounts_files(
            file_records=file_records,
            output_path=output_path,
            mode=mode,
        )
    else:
        summary = _combine_sparse_junction_featurecounts_files(
            file_records=file_records,
            output_path=output_path,
        )

    summary["input_dir"] = os.path.abspath(input_dir)
    summary["cell_id_list_path"] = (
        os.path.abspath(cell_id_list_path)
        if cell_id_list_path is not None else None
    )
    if summary_path is not None:
        summary_path = os.path.abspath(summary_path)
        os.makedirs(os.path.dirname(summary_path), exist_ok=True)
        with open(summary_path, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")
        summary["summary_path"] = summary_path
    return summary


def _normalize_merge_sample_labels(
    input_paths: List[str],
    sample_labels: Optional[List[str]] = None,
) -> List[Optional[str]]:
    if sample_labels is None:
        return [None] * len(input_paths)
    if len(sample_labels) != len(input_paths):
        raise ValueError(
            "sample_labels must either be omitted or provided once per input_path."
        )

    normalized_labels = []
    seen = set()
    for label in sample_labels:
        normalized = str(label).strip()
        if not normalized:
            raise ValueError("sample_labels may not be empty.")
        if any(character not in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-" for character in normalized):
            raise ValueError(
                "sample_labels may contain only letters, digits, '_' or '-'."
            )
        if normalized in seen:
            raise ValueError(f"Duplicate sample_label detected: {normalized}")
        seen.add(normalized)
        normalized_labels.append(normalized)
    return normalized_labels


def _prefixed_group_id(
    column_name: str,
    sample_label: Optional[str],
) -> str:
    group_id = derive_group_id_from_featurecounts_column(column_name)
    if sample_label:
        return f"{sample_label}__{group_id}"
    return group_id


def _merge_grouped_dense_featurecounts_files(
    input_paths: List[str],
    output_path: str,
    mode: str,
    sample_labels: List[Optional[str]],
) -> Dict[str, object]:
    config = FEATURECOUNTS_MODE_CONFIG[mode]
    handles = []
    first_comment_lines = []
    first_annotation_columns = None
    merged_group_ids = []
    seen_group_ids = set()
    row_count = 0

    try:
        for input_path, sample_label in zip(input_paths, sample_labels):
            handle = open(input_path, "r", encoding="utf-8")
            comment_lines, header = _read_featurecounts_header(handle, mode=mode)
            annotation_columns = header[: config["count_column_start"]]
            if first_annotation_columns is None:
                first_annotation_columns = annotation_columns
                first_comment_lines = comment_lines
            elif annotation_columns != first_annotation_columns:
                raise ValueError(
                    f"Annotation columns differ in grouped featureCounts file {input_path}."
                )

            grouped_columns = [
                _prefixed_group_id(column_name, sample_label)
                for column_name in header[config["count_column_start"]:]
            ]
            duplicate_ids = sorted(
                {
                    group_id
                    for group_id in grouped_columns
                    if group_id in seen_group_ids
                }
            )
            if duplicate_ids:
                raise ValueError(
                    "Merged grouped featureCounts would contain duplicate cell IDs. "
                    "Provide unique --sample-label values when merging multi-sample data. "
                    f"Examples: {duplicate_ids[:10]}"
                )
            seen_group_ids.update(grouped_columns)
            merged_group_ids.extend(grouped_columns)
            handles.append((input_path, handle, len(grouped_columns)))

        with open(output_path, "w", encoding="utf-8", newline="") as output_handle:
            for line in first_comment_lines:
                output_handle.write(line)
            output_handle.write(
                "\t".join(first_annotation_columns + merged_group_ids) + "\n"
            )

            while True:
                first_line = handles[0][1].readline()
                if first_line == "":
                    break
                first_fields = first_line.rstrip("\n").split("\t")
                if len(first_fields) < config["count_column_start"] + 1:
                    raise ValueError(
                        f"Malformed grouped featureCounts row in {handles[0][0]}: "
                        f"{first_line[:200]}"
                    )
                annotation_values = first_fields[: config["count_column_start"]]
                merged_counts = first_fields[config["count_column_start"]:]

                for input_path, handle, _ in handles[1:]:
                    line = handle.readline()
                    if line == "":
                        raise ValueError(
                            "Grouped featureCounts files ended at different row counts. "
                            f"First mismatch encountered in {input_path}"
                        )
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < config["count_column_start"] + 1:
                        raise ValueError(
                            f"Malformed grouped featureCounts row in {input_path}: {line[:200]}"
                        )
                    if fields[: config["count_column_start"]] != annotation_values:
                        raise ValueError(
                            "Grouped featureCounts annotation rows differ across inputs. "
                            f"First mismatch encountered in {input_path}"
                        )
                    merged_counts.extend(fields[config["count_column_start"]:])

                output_handle.write(
                    "\t".join(annotation_values + merged_counts) + "\n"
                )
                row_count += 1

            for input_path, handle, _ in handles[1:]:
                extra_line = handle.readline()
                if extra_line != "":
                    raise ValueError(
                        "Grouped featureCounts files ended at different row counts. "
                        f"Extra rows found in {input_path}"
                    )
    finally:
        for _, handle, _ in handles:
            handle.close()

    return {
        "mode": mode,
        "merge_strategy": "dense_row_aligned",
        "output_path": os.path.abspath(output_path),
        "group_count": len(merged_group_ids),
        "group_ids": merged_group_ids,
        "row_count": row_count,
        "input_paths": [os.path.abspath(path) for path in input_paths],
        "sample_labels": sample_labels,
    }


def _merge_grouped_junction_featurecounts_files(
    input_paths: List[str],
    output_path: str,
    sample_labels: List[Optional[str]],
) -> Dict[str, object]:
    config = FEATURECOUNTS_MODE_CONFIG["junction"]
    annotation_columns = None
    merged_group_ids = []
    seen_group_ids = set()
    handles = []
    row_count = 0

    try:
        for input_path, sample_label in zip(input_paths, sample_labels):
            handle = open(input_path, "r", encoding="utf-8")
            _, header = _read_featurecounts_header(handle, mode="junction")
            current_annotation_columns = header[: config["count_column_start"]]
            if annotation_columns is None:
                annotation_columns = current_annotation_columns
            elif current_annotation_columns != annotation_columns:
                raise ValueError(
                    f"Junction annotation columns differ in grouped featureCounts file {input_path}."
                )

            grouped_columns = [
                _prefixed_group_id(column_name, sample_label)
                for column_name in header[config["count_column_start"]:]
            ]
            duplicate_ids = sorted(
                {
                    group_id
                    for group_id in grouped_columns
                    if group_id in seen_group_ids
                }
            )
            if duplicate_ids:
                raise ValueError(
                    "Merged grouped junction featureCounts would contain duplicate cell IDs. "
                    "Provide unique --sample-label values when merging multi-sample data. "
                    f"Examples: {duplicate_ids[:10]}"
                )
            seen_group_ids.update(grouped_columns)
            handles.append((input_path, handle, grouped_columns))
            merged_group_ids.extend(grouped_columns)

        if annotation_columns is None:
            raise ValueError("No grouped junction featureCounts files were merged.")

        total_group_count = len(merged_group_ids)
        with open(output_path, "w", encoding="utf-8", newline="") as output_handle:
            writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
            writer.writerow(annotation_columns + merged_group_ids)

            offset = 0
            for input_path, handle, grouped_columns in handles:
                left_pad = ["0"] * offset
                right_pad = ["0"] * (total_group_count - offset - len(grouped_columns))
                for line in handle:
                    stripped = line.rstrip("\n")
                    if not stripped:
                        continue
                    row = stripped.split("\t")
                    if len(row) < config["count_column_start"] + len(grouped_columns):
                        raise ValueError(
                            f"Malformed grouped junction row in {input_path}: {stripped[:200]}"
                        )
                    writer.writerow(
                        row[: config["count_column_start"]]
                        + left_pad
                        + row[config["count_column_start"]:]
                        + right_pad
                    )
                    row_count += 1
                offset += len(grouped_columns)
    finally:
        for _, handle, _ in handles:
            handle.close()

    return {
        "mode": "junction",
        "merge_strategy": "sample_stacked_sparse",
        "output_path": os.path.abspath(output_path),
        "group_count": len(merged_group_ids),
        "group_ids": merged_group_ids,
        "row_count": row_count,
        "input_paths": [os.path.abspath(path) for path in input_paths],
        "sample_labels": sample_labels,
    }


def merge_grouped_featurecounts_files(
    input_paths: List[str],
    output_path: str,
    mode: str,
    sample_labels: Optional[List[str]] = None,
    summary_path: Optional[str] = None,
    cell_id_list_output_path: Optional[str] = None,
    cell_manifest_output_path: Optional[str] = None,
) -> Dict[str, object]:
    if mode not in FEATURECOUNTS_MODE_CONFIG:
        raise ValueError(
            f"Unsupported grouped merge mode: {mode}. "
            f"Supported modes: {sorted(FEATURECOUNTS_MODE_CONFIG)}"
        )
    if not input_paths:
        raise ValueError("At least one grouped featureCounts input_path is required.")

    normalized_input_paths = [os.path.abspath(path) for path in input_paths]
    for path in normalized_input_paths:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Grouped featureCounts input not found: {path}")

    normalized_labels = _normalize_merge_sample_labels(
        normalized_input_paths,
        sample_labels=sample_labels,
    )
    output_path = os.path.abspath(output_path)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    if mode in {"gene", "exon"}:
        summary = _merge_grouped_dense_featurecounts_files(
            input_paths=normalized_input_paths,
            output_path=output_path,
            mode=mode,
            sample_labels=normalized_labels,
        )
    else:
        summary = _merge_grouped_junction_featurecounts_files(
            input_paths=normalized_input_paths,
            output_path=output_path,
            sample_labels=normalized_labels,
        )

    if cell_id_list_output_path is not None:
        summary["cell_id_list_output_path"] = write_cell_id_list(
            summary["group_ids"],
            cell_id_list_output_path,
        )
    if cell_manifest_output_path is not None:
        summary["cell_manifest_output_path"] = write_cell_manifest(
            summary["group_ids"],
            cell_manifest_output_path,
        )
    if summary_path is not None:
        summary_path = os.path.abspath(summary_path)
        os.makedirs(os.path.dirname(summary_path), exist_ok=True)
        with open(summary_path, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")
        summary["summary_path"] = summary_path
    return summary


def build_full_length_preprocess_plan(
    reference_manifest_path: str,
    metadata_path: Optional[str],
    raw_fastq_dir: str,
    output_dir: str,
    read_layout: str = "single",
    fastq_suffix: str = ".fastq.gz",
    star_binary: str = "STAR",
    featurecounts_binary: str = "featureCounts",
    star_threads: int = 4,
    featurecounts_threads: int = 4,
    trim: bool = False,
    java_binary: str = "java",
    trimmomatic_jar: Optional[str] = None,
    adapter_fasta: Optional[str] = None,
    skip_missing_fastqs: bool = True,
    explicit_star_index_dir: Optional[str] = None,
) -> Dict[str, object]:
    if read_layout not in {"single", "paired"}:
        raise ValueError("read_layout must be 'single' or 'paired'")
    if trim and (trimmomatic_jar is None or adapter_fasta is None):
        raise ValueError(
            "Trimmomatic trimming requires both trimmomatic_jar and adapter_fasta."
        )

    reference_manifest = load_reference_manifest(reference_manifest_path)
    star_index_dir, warnings = resolve_star_index_dir(
        reference_manifest,
        explicit_star_index_dir=explicit_star_index_dir,
    )
    if metadata_path is not None:
        metadata_df = read_metadata_table(metadata_path)
        samples, missing_fastqs = discover_full_length_fastqs(
            cell_ids=metadata_df["cell_ids"],
            raw_fastq_dir=raw_fastq_dir,
            read_layout=read_layout,
            fastq_suffix=fastq_suffix,
            skip_missing_fastqs=skip_missing_fastqs,
        )
        input_cell_ids = metadata_df["cell_ids"]
        input_discovery = "metadata"
    else:
        metadata_df = None
        samples = discover_full_length_fastqs_from_dir(
            raw_fastq_dir=raw_fastq_dir,
            read_layout=read_layout,
            fastq_suffix=fastq_suffix,
        )
        missing_fastqs = []
        input_cell_ids = [sample["cell_id"] for sample in samples]
        input_discovery = "raw_fastq_dir"
        warnings.append(
            "No metadata provided. Full-length samples were inferred directly "
            "from FASTQ filenames."
        )

    if not samples:
        raise ValueError(
            f"No full-length FASTQ inputs were discovered in {os.path.abspath(raw_fastq_dir)}"
        )

    output_dir = _ensure_directory(output_dir)
    trim_dir = _ensure_directory(os.path.join(output_dir, "01_trim"))
    exon_star_dir = _ensure_directory(os.path.join(output_dir, "03_exon_star"))
    exon_gene_dir = _ensure_directory(os.path.join(output_dir, "04_exon_gene_cnt"))
    exon_gene_grouped_dir = _ensure_directory(os.path.join(output_dir, "04_exon_gene_cnt_grouped"))
    exon_junction_dir = _ensure_directory(os.path.join(output_dir, "05_exon_junct_cnt"))
    exon_junction_grouped_dir = _ensure_directory(os.path.join(output_dir, "05_exon_junct_cnt_grouped"))
    logs_dir = _ensure_directory(os.path.join(output_dir, "logs"))
    cell_id_list_path = write_cell_id_list(
        [sample["cell_id"] for sample in samples],
        os.path.join(output_dir, "cell_ids.txt"),
    )
    cell_manifest_path = write_cell_manifest(
        [sample["cell_id"] for sample in samples],
        os.path.join(output_dir, "cell_manifest.tsv"),
    )

    stages = []
    aligned_bam_paths = []
    for sample in samples:
        cell_id = sample["cell_id"]
        read_inputs = [sample["read1_path"]]
        if sample["read_layout"] == "paired":
            read_inputs.append(sample["read2_path"])

        working_reads = list(read_inputs)
        if trim:
            if sample["read_layout"] == "single":
                trimmed_paths = [
                    os.path.join(trim_dir, f"{cell_id}.trim.fastq.gz")
                ]
                working_reads_after_trim = trimmed_paths
            else:
                trimmed_paths = [
                    os.path.join(trim_dir, f"{cell_id}_1.paired.trim.fastq.gz"),
                    os.path.join(trim_dir, f"{cell_id}_1.unpaired.trim.fastq.gz"),
                    os.path.join(trim_dir, f"{cell_id}_2.paired.trim.fastq.gz"),
                    os.path.join(trim_dir, f"{cell_id}_2.unpaired.trim.fastq.gz"),
                ]
                working_reads_after_trim = [trimmed_paths[0], trimmed_paths[2]]
            trimmomatic_args = [
                java_binary,
                "-jar",
                trimmomatic_jar,
                "SE" if sample["read_layout"] == "single" else "PE",
                *working_reads,
                *trimmed_paths,
                f"ILLUMINACLIP:{adapter_fasta}:2:30:10",
                "LEADING:3",
                "TRAILING:3",
                "SLIDINGWINDOW:4:15",
                "MINLEN:36",
            ]
            stages.append(
                _plan_record(
                    stage_name=f"trim.{cell_id}",
                    description=f"Trim FASTQ reads for full-length cell {cell_id}.",
                    command=_shell_join(trimmomatic_args),
                    outputs=trimmed_paths,
                    mode="full_length",
                    sample_id=cell_id,
                )
            )
            working_reads = working_reads_after_trim

        sample_star_dir = _ensure_directory(os.path.join(exon_star_dir, cell_id))
        star_prefix = os.path.join(sample_star_dir, f"{cell_id}.")
        star_args = [
            star_binary,
            "--runThreadN",
            str(int(star_threads)),
            "--genomeDir",
            star_index_dir,
            "--readFilesIn",
            *working_reads,
            "--readFilesCommand",
            "gunzip",
            "-c",
            "--outSAMtype",
            "BAM",
            "SortedByCoordinate",
            "--outFileNamePrefix",
            star_prefix,
        ]
        stages.append(
            _plan_record(
                stage_name=f"align_exon.{cell_id}",
                description=f"Align full-length cell {cell_id} against the DOLPHIN STAR index.",
                command=_shell_join(star_args),
                outputs=[
                    os.path.join(sample_star_dir, f"{cell_id}.Aligned.sortedByCoord.out.bam"),
                ],
                mode="full_length",
                sample_id=cell_id,
            )
        )
        aligned_bam_paths.append(
            os.path.join(sample_star_dir, f"{cell_id}.Aligned.sortedByCoord.out.bam")
        )

        gene_count_path = os.path.join(exon_gene_dir, f"{cell_id}.exongene.count.txt")
        gene_count_args = [
            featurecounts_binary,
            "-T",
            str(int(featurecounts_threads)),
            "-t",
            "exon",
            "-O",
            "-M",
            "-a",
            reference_manifest["alignment_gtf_path"],
            "-o",
            gene_count_path,
            os.path.join(sample_star_dir, f"{cell_id}.Aligned.sortedByCoord.out.bam"),
        ]
        stages.append(
            _plan_record(
                stage_name=f"count_gene.{cell_id}",
                description=f"Generate exon-gene counts for full-length cell {cell_id}.",
                command=_shell_join(gene_count_args),
                outputs=[gene_count_path],
                mode="full_length",
                sample_id=cell_id,
            )
        )

        exon_count_path = os.path.join(exon_junction_dir, f"{cell_id}.exon.count.txt")
        exon_count_args = [
            featurecounts_binary,
            "-T",
            str(int(featurecounts_threads)),
            "-t",
            "exon",
            "-f",
            "-O",
            "-J",
            "-M",
            "-a",
            reference_manifest["alignment_gtf_path"],
            "-o",
            exon_count_path,
            os.path.join(sample_star_dir, f"{cell_id}.Aligned.sortedByCoord.out.bam"),
        ]
        stages.append(
            _plan_record(
                stage_name=f"count_exon_junction.{cell_id}",
                description=f"Generate exon and junction counts for full-length cell {cell_id}.",
                command=_shell_join(exon_count_args),
                outputs=[exon_count_path, f"{exon_count_path}.jcounts"],
                mode="full_length",
                sample_id=cell_id,
            )
        )

    grouped_gene_count_path = os.path.join(
        exon_gene_grouped_dir,
        "full_length.exongene.count.txt",
    )
    grouped_gene_summary_path = os.path.join(
        exon_gene_grouped_dir,
        "full_length.exongene.count.txt.summary.json",
    )
    grouped_gene_args = [
        sys.executable,
        "-m",
        "DOLPHIN.preprocess.cli",
        "combine-featurecounts",
        "--mode",
        "gene",
        "--input-dir",
        exon_gene_dir,
        "--output-path",
        grouped_gene_count_path,
        "--cell-id-list-path",
        cell_id_list_path,
        "--summary-path",
        grouped_gene_summary_path,
    ]
    grouped_gene_validate_args = [
        sys.executable,
        "-m",
        "DOLPHIN.preprocess.cli",
        "validate-grouped-featurecounts",
        "--mode",
        "gene",
        "--input-path",
        grouped_gene_count_path,
    ]
    stages.append(
        _plan_record(
            stage_name="grouped_gene_count.full_length",
            description=(
                "Combine per-cell exon-gene featureCounts tables into one grouped matrix "
                "without changing the original per-cell results."
            ),
            command=(
                f"{_shell_join(grouped_gene_args)} && "
                f"{_shell_join(grouped_gene_validate_args)}"
            ),
            outputs=[grouped_gene_count_path],
            mode="full_length",
        )
    )

    grouped_exon_count_path = os.path.join(
        exon_junction_grouped_dir,
        "full_length.exon.count.txt",
    )
    grouped_exon_summary_path = os.path.join(
        exon_junction_grouped_dir,
        "full_length.exon.count.txt.summary.json",
    )
    grouped_exon_args = [
        sys.executable,
        "-m",
        "DOLPHIN.preprocess.cli",
        "combine-featurecounts",
        "--mode",
        "exon",
        "--input-dir",
        exon_junction_dir,
        "--output-path",
        grouped_exon_count_path,
        "--cell-id-list-path",
        cell_id_list_path,
        "--summary-path",
        grouped_exon_summary_path,
    ]
    grouped_junction_summary_path = os.path.join(
        exon_junction_grouped_dir,
        "full_length.exon.count.txt.jcounts.summary.json",
    )
    grouped_junction_args = [
        sys.executable,
        "-m",
        "DOLPHIN.preprocess.cli",
        "combine-featurecounts",
        "--mode",
        "junction",
        "--input-dir",
        exon_junction_dir,
        "--output-path",
        f"{grouped_exon_count_path}.jcounts",
        "--cell-id-list-path",
        cell_id_list_path,
        "--summary-path",
        grouped_junction_summary_path,
    ]
    grouped_exon_validate_args = [
        sys.executable,
        "-m",
        "DOLPHIN.preprocess.cli",
        "validate-grouped-featurecounts",
        "--mode",
        "exon",
        "--input-path",
        grouped_exon_count_path,
    ]
    grouped_junction_validate_args = [
        sys.executable,
        "-m",
        "DOLPHIN.preprocess.cli",
        "validate-grouped-featurecounts",
        "--mode",
        "junction",
        "--input-path",
        f"{grouped_exon_count_path}.jcounts",
    ]
    stages.append(
        _plan_record(
            stage_name="grouped_exon_junction_count.full_length",
            description=(
                "Combine per-cell exon and junction featureCounts tables into grouped "
                "matrices while preserving the original per-cell outputs."
            ),
            command=(
                f"{_shell_join(grouped_exon_args)} && "
                f"{_shell_join(grouped_junction_args)} && "
                f"{_shell_join(grouped_exon_validate_args)} && "
                f"{_shell_join(grouped_junction_validate_args)}"
            ),
            outputs=[grouped_exon_count_path, f"{grouped_exon_count_path}.jcounts"],
            mode="full_length",
        )
    )

    parallel_recommendations = recommend_full_length_parallelism(
        sample_count=len(samples),
        star_threads=int(star_threads),
        featurecounts_threads=int(featurecounts_threads),
    )
    return {
        "mode": "full_length",
        "project_root": _project_root(),
        "reference_manifest_path": os.path.abspath(reference_manifest_path),
        "reference_manifest": reference_manifest,
        "metadata_path": (
            os.path.abspath(metadata_path) if metadata_path is not None else None
        ),
        "input_discovery": input_discovery,
        "input_cell_count": len(input_cell_ids),
        "raw_fastq_dir": os.path.abspath(raw_fastq_dir),
        "output_dir": output_dir,
        "star_index_dir": star_index_dir,
        "warnings": warnings,
        "missing_fastqs": missing_fastqs,
        "sample_count": len(samples),
        "star_threads_per_cell": int(star_threads),
        "featurecounts_threads_per_cell": int(featurecounts_threads),
        "parallel_recommendations": parallel_recommendations,
        "cell_ids": [sample["cell_id"] for sample in samples],
        "cell_id_list_path": cell_id_list_path,
        "cell_manifest_path": cell_manifest_path,
        "logs_dir": logs_dir,
        "stages": stages,
    }


def build_tenx_preprocess_plan(
    reference_manifest_path: str,
    metadata_path: Optional[str],
    raw_fastq_dir: str,
    output_dir: str,
    solo_whitelist_path: str,
    sample_name: Optional[str] = None,
    chemistry: str = "10xv2",
    star_binary: str = "STAR",
    featurecounts_binary: str = "featureCounts",
    samtools_binary: str = "samtools",
    star_threads: int = 16,
    featurecounts_threads: int = 8,
    solo_cell_filter_args: Optional[List[str]] = None,
    barcode_tag: str = "CB",
    rg_tag: str = "RG",
    append_suffix: Optional[str] = None,
    metadata_cb_suffix: Optional[str] = None,
    explicit_star_index_dir: Optional[str] = None,
) -> Dict[str, object]:
    chemistry = chemistry.lower()
    if chemistry not in TENX_CHEMISTRY_PRESETS:
        raise ValueError(
            f"Unsupported 10x chemistry preset: {chemistry}. "
            f"Supported presets: {sorted(TENX_CHEMISTRY_PRESETS)}"
        )

    reference_manifest = load_reference_manifest(reference_manifest_path)
    star_index_dir, warnings = resolve_star_index_dir(
        reference_manifest,
        explicit_star_index_dir=explicit_star_index_dir,
    )
    if metadata_path is not None:
        metadata_df = filter_metadata_by_cb_suffix(
            read_metadata_table(metadata_path),
            cb_suffix_filter=metadata_cb_suffix,
        )
    else:
        metadata_df = None
        if metadata_cb_suffix is not None:
            raise ValueError(
                "--metadata-cb-suffix requires --metadata-path."
            )
    fastq_pair = discover_tenx_fastq_pair(raw_fastq_dir=raw_fastq_dir, sample_name=sample_name)
    chemistry_preset = TENX_CHEMISTRY_PRESETS[chemistry]

    suffix_counts = {}
    if metadata_df is not None:
        suffix_counts = collect_cb_suffix_counts(metadata_df["cell_ids"])
        if append_suffix is None:
            append_suffix = infer_common_cb_suffix(metadata_df["cell_ids"])
        if append_suffix is None and suffix_counts:
            raise ValueError(
                "10x metadata CB values include multiple suffix groups. "
                "Please provide metadata for one sample or pass "
                "--metadata-cb-suffix so the CB->RG conversion can match "
                "the metadata cell IDs."
            )
        if metadata_cb_suffix is not None:
            warnings.append(
                f"Filtered 10x metadata rows by CB suffix: {metadata_cb_suffix}"
            )
    else:
        warnings.append(
            "No metadata provided. 10x cell selection will rely on STARsolo "
            "barcode outputs instead of annotation metadata."
        )

    output_dir = _ensure_directory(output_dir)
    sample_output_dir = _ensure_directory(os.path.join(output_dir, fastq_pair["sample_name"]))
    exon_star_dir = _ensure_directory(os.path.join(sample_output_dir, "03_exon_star"))
    rg_bam_dir = _ensure_directory(os.path.join(sample_output_dir, "04_rg_tagged_bam"))
    gene_grouped_dir = _ensure_directory(os.path.join(sample_output_dir, "05_exon_gene_cnt_grouped"))
    gene_split_dir = _ensure_directory(os.path.join(sample_output_dir, "05_exon_gene_cnt"))
    exon_grouped_dir = _ensure_directory(os.path.join(sample_output_dir, "06_exon_junct_cnt_grouped"))
    exon_split_dir = _ensure_directory(os.path.join(sample_output_dir, "06_exon_junct_cnt"))
    logs_dir = _ensure_directory(os.path.join(sample_output_dir, "logs"))

    sample_star_dir = _ensure_directory(os.path.join(exon_star_dir, fastq_pair["sample_name"]))
    star_prefix = os.path.join(sample_star_dir, f"{fastq_pair['sample_name']}.")
    aligned_bam_path = os.path.join(
        sample_star_dir,
        f"{fastq_pair['sample_name']}.Aligned.sortedByCoord.out.bam",
    )
    solo_out_dir = os.path.join(
        sample_star_dir,
        f"{fastq_pair['sample_name']}.Solo.out",
    )
    filtered_barcode_list_path = os.path.join(
        solo_out_dir,
        "Gene",
        "filtered",
        "barcodes.tsv",
    )
    raw_barcode_list_path = os.path.join(
        solo_out_dir,
        "Gene",
        "raw",
        "barcodes.tsv",
    )
    rg_bam_path = os.path.join(
        rg_bam_dir,
        f"{fastq_pair['sample_name']}.cb_rg.bam",
    )
    if solo_cell_filter_args is None:
        solo_cell_filter_args = ["CellRanger2.2", "3000", "0.99", "10"]
    use_filtered_barcodes = not (
        len(solo_cell_filter_args) == 1 and solo_cell_filter_args[0].lower() == "none"
    )

    starsolo_args = [
        star_binary,
        "--runThreadN",
        str(int(star_threads)),
        "--genomeDir",
        star_index_dir,
        "--readFilesIn",
        fastq_pair["read2_path"],
        fastq_pair["read1_path"],
        "--readFilesCommand",
        "gunzip",
        "-c",
        "--soloType",
        chemistry_preset["solo_type"],
        "--soloCBwhitelist",
        os.path.abspath(solo_whitelist_path),
        "--soloCBstart",
        str(int(chemistry_preset["solo_cb_start"])),
        "--soloCBlen",
        str(int(chemistry_preset["solo_cb_len"])),
        "--soloUMIstart",
        str(int(chemistry_preset["solo_umi_start"])),
        "--soloUMIlen",
        str(int(chemistry_preset["solo_umi_len"])),
        "--soloBarcodeReadLength",
        str(int(chemistry_preset["solo_barcode_read_length"])),
        "--soloCellFilter",
        *solo_cell_filter_args,
        "--soloFeatures",
        "Gene",
        "GeneFull",
        "SJ",
        "--outSAMtype",
        "BAM",
        "SortedByCoordinate",
        "--outSAMattributes",
        "NH",
        "HI",
        "nM",
        "AS",
        "CR",
        "UR",
        "CB",
        "UB",
        "GX",
        "GN",
        "--outFileNamePrefix",
        star_prefix,
    ]

    stages = [
        _plan_record(
            stage_name=f"align_starsolo.{fastq_pair['sample_name']}",
            description=(
                "Align the 10x library once with STARsolo against the shared "
                "DOLPHIN STAR index."
            ),
            command=_shell_join(starsolo_args),
            outputs=(
                [aligned_bam_path, filtered_barcode_list_path]
                if use_filtered_barcodes
                else [aligned_bam_path, raw_barcode_list_path]
            ),
            mode="tenx",
            sample_id=fastq_pair["sample_name"],
        )
    ]

    cb_to_rg_filter_args = [
        sys.executable,
        "-m",
        "DOLPHIN.preprocess.cli",
        "cb-to-rg",
        "--barcode-tag",
        barcode_tag,
        "--rg-tag",
        rg_tag,
    ]
    if use_filtered_barcodes:
        cb_to_rg_filter_args.extend(["--allow-list-path", filtered_barcode_list_path])
    if append_suffix:
        cb_to_rg_filter_args.extend(["--append-suffix", append_suffix])
    samtools_view_args = [
        samtools_binary,
        "view",
        "-h",
    ]
    if use_filtered_barcodes:
        samtools_view_args.extend(
            ["-D", f"{barcode_tag}:{filtered_barcode_list_path}"]
        )
    samtools_view_args.append(aligned_bam_path)
    rg_bam_command = (
        f"{_shell_join(samtools_view_args)} | "
        f"{_shell_join(cb_to_rg_filter_args)} | "
        f"{_quote(samtools_binary)} view -b -o {_quote(rg_bam_path)} -"
    )
    stages.append(
        _plan_record(
            stage_name=f"cb_to_rg.{fastq_pair['sample_name']}",
            description=(
                "Copy the corrected cell barcode tag into RG so featureCounts "
                "can group counts per cell in one pooled pass."
            ),
            command=rg_bam_command,
            outputs=[rg_bam_path],
            mode="tenx",
            sample_id=fastq_pair["sample_name"],
        )
    )

    grouped_gene_count_path = os.path.join(
        gene_grouped_dir,
        f"{fastq_pair['sample_name']}.exongene.count.txt",
    )
    grouped_gene_args = [
        featurecounts_binary,
        "-T",
        str(int(featurecounts_threads)),
        "-t",
        "exon",
        "-O",
        "-M",
        "--byReadGroup",
        "-a",
        reference_manifest["alignment_gtf_path"],
        "-o",
        grouped_gene_count_path,
        rg_bam_path,
    ]
    grouped_gene_validate_args = [
        sys.executable,
        "-m",
        "DOLPHIN.preprocess.cli",
        "validate-grouped-featurecounts",
        "--mode",
        "gene",
        "--input-path",
        grouped_gene_count_path,
    ]
    stages.append(
        _plan_record(
            stage_name=f"grouped_gene_count.{fastq_pair['sample_name']}",
            description="Generate pooled per-cell exon-gene counts using featureCounts --byReadGroup.",
            command=(
                f"{_shell_join(grouped_gene_args)} && "
                f"{_shell_join(grouped_gene_validate_args)}"
            ),
            outputs=[grouped_gene_count_path],
            mode="tenx",
            sample_id=fastq_pair["sample_name"],
        )
    )

    grouped_exon_count_path = os.path.join(
        exon_grouped_dir,
        f"{fastq_pair['sample_name']}.exon.count.txt",
    )
    grouped_exon_args = [
        featurecounts_binary,
        "-T",
        str(int(featurecounts_threads)),
        "-t",
        "exon",
        "-f",
        "-O",
        "-J",
        "-M",
        "--byReadGroup",
        "-a",
        reference_manifest["alignment_gtf_path"],
        "-o",
        grouped_exon_count_path,
        rg_bam_path,
    ]
    grouped_exon_validate_args = [
        sys.executable,
        "-m",
        "DOLPHIN.preprocess.cli",
        "validate-grouped-featurecounts",
        "--mode",
        "exon",
        "--input-path",
        grouped_exon_count_path,
    ]
    grouped_junction_validate_args = [
        sys.executable,
        "-m",
        "DOLPHIN.preprocess.cli",
        "validate-grouped-featurecounts",
        "--mode",
        "junction",
        "--input-path",
        f"{grouped_exon_count_path}.jcounts",
    ]
    stages.append(
        _plan_record(
            stage_name=f"grouped_exon_junction_count.{fastq_pair['sample_name']}",
            description="Generate pooled per-cell exon and junction counts using featureCounts --byReadGroup.",
            command=(
                f"{_shell_join(grouped_exon_args)} && "
                f"{_shell_join(grouped_exon_validate_args)} && "
                f"{_shell_join(grouped_junction_validate_args)}"
            ),
            outputs=[grouped_exon_count_path, f"{grouped_exon_count_path}.jcounts"],
            mode="tenx",
            sample_id=fastq_pair["sample_name"],
        )
    )

    split_gene_summary_path = os.path.join(
        gene_split_dir,
        "_split_gene_count.summary.json",
    )

    split_gene_args = [
        sys.executable,
        "-m",
        "DOLPHIN.preprocess.cli",
        "split-featurecounts",
        "--mode",
        "gene",
        "--input-path",
        grouped_gene_count_path,
        "--output-dir",
        gene_split_dir,
        "--summary-path",
        split_gene_summary_path,
    ]
    if metadata_path is not None:
        split_gene_args.extend(["--metadata-path", os.path.abspath(metadata_path)])
    elif use_filtered_barcodes:
        split_gene_args.extend(["--cell-id-list-path", filtered_barcode_list_path])
    stages.append(
        _plan_record(
            stage_name=f"split_gene_count.{fastq_pair['sample_name']}",
            description="Split grouped exon-gene counts back into one file per cell for DOLPHIN compatibility.",
            command=_shell_join(split_gene_args),
            outputs=[split_gene_summary_path],
            mode="tenx",
            sample_id=fastq_pair["sample_name"],
        )
    )

    split_exon_summary_path = os.path.join(
        exon_split_dir,
        "_split_exon_count.summary.json",
    )
    split_exon_args = [
        sys.executable,
        "-m",
        "DOLPHIN.preprocess.cli",
        "split-featurecounts",
        "--mode",
        "exon",
        "--input-path",
        grouped_exon_count_path,
        "--output-dir",
        exon_split_dir,
        "--summary-path",
        split_exon_summary_path,
    ]
    split_junction_summary_path = os.path.join(
        exon_split_dir,
        "_split_junction_count.summary.json",
    )
    split_junction_args = [
        sys.executable,
        "-m",
        "DOLPHIN.preprocess.cli",
        "split-featurecounts",
        "--mode",
        "junction",
        "--input-path",
        f"{grouped_exon_count_path}.jcounts",
        "--output-dir",
        exon_split_dir,
        "--summary-path",
        split_junction_summary_path,
    ]
    if metadata_path is not None:
        split_exon_args.extend(["--metadata-path", os.path.abspath(metadata_path)])
        split_junction_args.extend(["--metadata-path", os.path.abspath(metadata_path)])
    elif use_filtered_barcodes:
        split_exon_args.extend(["--cell-id-list-path", filtered_barcode_list_path])
        split_junction_args.extend(["--cell-id-list-path", filtered_barcode_list_path])
    else:
        warnings.append(
            "STARsolo cell filtering is disabled and no metadata was provided. "
            "Grouped featureCounts outputs will be split for every grouped barcode."
        )
    stages.append(
        _plan_record(
            stage_name=f"split_exon_count.{fastq_pair['sample_name']}",
            description="Split grouped exon counts into one file per cell for DOLPHIN compatibility.",
            command=_shell_join(split_exon_args),
            outputs=[split_exon_summary_path],
            mode="tenx",
            sample_id=fastq_pair["sample_name"],
        )
    )
    stages.append(
        _plan_record(
            stage_name=f"split_junction_count.{fastq_pair['sample_name']}",
            description="Split grouped junction counts into one file per cell for DOLPHIN compatibility.",
            command=_shell_join(split_junction_args),
            outputs=[split_junction_summary_path],
            mode="tenx",
            sample_id=fastq_pair["sample_name"],
        )
    )

    return {
        "mode": "tenx",
        "project_root": _project_root(),
        "reference_manifest_path": os.path.abspath(reference_manifest_path),
        "reference_manifest": reference_manifest,
        "metadata_path": (
            os.path.abspath(metadata_path) if metadata_path is not None else None
        ),
        "raw_fastq_dir": os.path.abspath(raw_fastq_dir),
        "output_dir": sample_output_dir,
        "sample_name": fastq_pair["sample_name"],
        "fastq_pair": fastq_pair,
        "star_index_dir": star_index_dir,
        "solo_whitelist_path": os.path.abspath(solo_whitelist_path),
        "solo_cell_filter_args": solo_cell_filter_args,
        "filtered_barcode_list_path": (
            filtered_barcode_list_path if use_filtered_barcodes else None
        ),
        "chemistry": chemistry,
        "append_suffix": append_suffix,
        "metadata_cb_suffix": metadata_cb_suffix,
        "cell_count": (
            len(metadata_df["cell_ids"]) if metadata_df is not None else None
        ),
        "cell_ids": (
            metadata_df["cell_ids"] if metadata_df is not None else []
        ),
        "metadata_suffix_counts": suffix_counts,
        "warnings": warnings,
        "logs_dir": logs_dir,
        "stages": stages,
    }


def write_preprocess_plan(
    plan: Dict[str, object],
    output_dir: str,
    script_name: str = "run_preprocess.sh",
    manifest_name: str = "preprocess_plan.json",
) -> Dict[str, str]:
    output_dir = _ensure_directory(output_dir)
    manifest_path = os.path.join(output_dir, manifest_name)
    script_path = os.path.join(output_dir, script_name)

    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(plan, handle, indent=2, sort_keys=True)

    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        f"# mode: {plan['mode']}",
        f"# reference_manifest: {plan['reference_manifest_path']}",
        f"export PYTHONPATH={_quote(plan.get('project_root', _project_root()))}:${{PYTHONPATH:-}}",
        "",
    ]
    for stage in plan["stages"]:
        lines.append(f"# [{stage['stage_name']}] {stage['description']}")
        lines.append(stage["command"])
        lines.append("")

    with open(script_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
        handle.write("\n")
    os.chmod(script_path, 0o755)

    return {
        "manifest_path": manifest_path,
        "script_path": script_path,
    }


def _sanitize_stage_name(stage_name: str) -> str:
    return "".join(
        character if character.isalnum() or character in {"-", "_", "."} else "_"
        for character in stage_name
    )


def _write_execution_records(
    records: List[Dict[str, object]],
    timing_path: str,
    summary_path: str,
) -> None:
    with open(timing_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "stage_name",
                "status",
                "started_at",
                "finished_at",
                "elapsed_seconds",
                "log_path",
            ]
        )
        for record in records:
            writer.writerow(
                [
                    record["stage_name"],
                    record["status"],
                    record["started_at"],
                    record["finished_at"],
                    f"{float(record['elapsed_seconds']):.6f}",
                    record["log_path"],
                ]
            )

    summary = {
        "stage_count": len(records),
        "completed_stage_count": sum(
            1 for record in records if record["status"] == "completed"
        ),
        "skipped_stage_count": sum(
            1 for record in records if record["status"] == "skipped_existing"
        ),
        "failed_stage_count": sum(
            1 for record in records if record["status"] == "failed"
        ),
        "stages": records,
    }
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")


def run_preprocess_plan(
    plan_manifest_path: str,
    log_dir: Optional[str] = None,
    timing_path: Optional[str] = None,
    summary_path: Optional[str] = None,
    resume: bool = False,
    continue_on_error: bool = False,
    max_parallel_cells: Optional[int] = None,
    parallel_profile: Optional[str] = None,
) -> Dict[str, object]:
    plan_manifest_path = os.path.abspath(plan_manifest_path)
    if not os.path.isfile(plan_manifest_path):
        raise FileNotFoundError(f"Preprocess plan manifest not found: {plan_manifest_path}")

    with open(plan_manifest_path, "r", encoding="utf-8") as handle:
        plan = json.load(handle)

    plan_output_dir = plan.get("output_dir") or os.path.dirname(plan_manifest_path)
    default_log_dir = plan.get("logs_dir") or os.path.join(plan_output_dir, "logs")
    log_dir = _ensure_directory(log_dir or default_log_dir)
    timing_path = os.path.abspath(
        timing_path or os.path.join(log_dir, "execution_timing.tsv")
    )
    summary_path = os.path.abspath(
        summary_path or os.path.join(log_dir, "execution_summary.json")
    )
    project_root = os.path.abspath(plan.get("project_root", _project_root()))

    env = os.environ.copy()
    existing_pythonpath = env.get("PYTHONPATH")
    env["PYTHONPATH"] = (
        f"{project_root}:{existing_pythonpath}"
        if existing_pythonpath else project_root
    )

    if max_parallel_cells is not None and max_parallel_cells < 1:
        raise ValueError("max_parallel_cells must be >= 1")
    if max_parallel_cells is not None and parallel_profile is not None:
        raise ValueError(
            "Use either max_parallel_cells or parallel_profile, not both."
        )

    resolved_parallel_profile = parallel_profile
    resolved_max_parallel_cells = (
        int(max_parallel_cells) if max_parallel_cells is not None else None
    )
    if parallel_profile is not None:
        if plan.get("mode") == "full_length":
            recommendations = get_full_length_parallel_recommendations_from_plan(plan)
        else:
            recommendations = plan.get("parallel_recommendations") or {}
        profiles = recommendations.get("profiles") or {}
        selected_profile = profiles.get(parallel_profile)
        if selected_profile is None:
            raise ValueError(
                f"Parallel profile '{parallel_profile}' is not available in the plan."
            )
        resolved_max_parallel_cells = int(selected_profile["max_parallel_cells"])

    prepared_stages = []
    for index, stage in enumerate(plan["stages"], start=1):
        prepared_stage = dict(stage)
        prepared_stage["stage_index"] = index
        prepared_stage["outputs"] = [
            os.path.abspath(path) for path in stage.get("outputs", [])
        ]
        prepared_stage["log_path"] = os.path.join(
            log_dir,
            f"{index:04d}_{_sanitize_stage_name(stage['stage_name'])}.log",
        )
        prepared_stages.append(prepared_stage)

    if (
        plan.get("mode") == "full_length"
        and resolved_max_parallel_cells is not None
        and resolved_max_parallel_cells > 1
    ):
        records = _run_full_length_plan_in_parallel(
            prepared_stages=prepared_stages,
            project_root=project_root,
            env=env,
            timing_path=timing_path,
            summary_path=summary_path,
            resume=resume,
            continue_on_error=continue_on_error,
            max_parallel_cells=resolved_max_parallel_cells,
        )
    else:
        records = []
        for stage in prepared_stages:
            record = _execute_preprocess_stage(
                stage=stage,
                project_root=project_root,
                env=env,
                resume=resume,
            )
            records.append(record)
            _write_execution_records(
                records,
                timing_path=timing_path,
                summary_path=summary_path,
            )

            if record["status"] == "failed" and not continue_on_error:
                break

    final_summary = {
        "plan_manifest_path": plan_manifest_path,
        "mode": plan.get("mode"),
        "project_root": project_root,
        "log_dir": log_dir,
        "timing_path": timing_path,
        "summary_path": summary_path,
        "parallel_profile": resolved_parallel_profile,
        "max_parallel_cells": resolved_max_parallel_cells,
        "stage_count": len(plan["stages"]),
        "completed_stage_count": sum(
            1 for record in records if record["status"] == "completed"
        ),
        "skipped_stage_count": sum(
            1 for record in records if record["status"] == "skipped_existing"
        ),
        "failed_stage_count": sum(
            1 for record in records if record["status"] == "failed"
        ),
        "records": records,
    }
    return final_summary


def _execute_preprocess_stage(
    stage: Dict[str, object],
    project_root: str,
    env: Dict[str, str],
    resume: bool = False,
) -> Dict[str, object]:
    stage_name = stage["stage_name"]
    outputs = [os.path.abspath(path) for path in stage.get("outputs", [])]
    log_path = os.path.abspath(stage["log_path"])

    started_at = datetime.datetime.now(datetime.timezone.utc).isoformat()
    start_time = time.time()
    status = "completed"
    return_code = 0

    if resume and outputs and all(_output_path_is_complete(path) for path in outputs):
        status = "skipped_existing"
        with open(log_path, "w", encoding="utf-8") as handle:
            handle.write("Skipped because all declared outputs already exist.\n")
    else:
        with open(log_path, "w", encoding="utf-8") as handle:
            handle.write(f"# stage_name: {stage_name}\n")
            handle.write(f"# description: {stage.get('description', '')}\n")
            handle.write(f"# command: {stage['command']}\n\n")
            handle.flush()
            completed = subprocess.run(
                stage["command"],
                shell=True,
                executable="/bin/bash",
                cwd=project_root,
                env=env,
                stdout=handle,
                stderr=subprocess.STDOUT,
                check=False,
                text=True,
            )
            return_code = completed.returncode
            if return_code != 0:
                status = "failed"

    finished_at = datetime.datetime.now(datetime.timezone.utc).isoformat()
    elapsed_seconds = time.time() - start_time
    return {
        "stage_index": int(stage["stage_index"]),
        "stage_name": stage_name,
        "status": status,
        "started_at": started_at,
        "finished_at": finished_at,
        "elapsed_seconds": elapsed_seconds,
        "log_path": log_path,
        "outputs": outputs,
        "return_code": return_code,
    }


def _group_stages_by_sample(
    prepared_stages: List[Dict[str, object]],
) -> List[List[Dict[str, object]]]:
    grouped: Dict[str, List[Dict[str, object]]] = {}
    for stage in prepared_stages:
        sample_id = stage.get("sample_id")
        if sample_id is None:
            stage_name = stage["stage_name"]
            sample_id = stage_name.split(".", 1)[-1] if "." in stage_name else stage_name
        grouped.setdefault(str(sample_id), []).append(stage)
    return [grouped[key] for key in grouped]


def _run_sample_stage_group(
    sample_stages: List[Dict[str, object]],
    project_root: str,
    env: Dict[str, str],
    resume: bool,
    continue_on_error: bool,
) -> List[Dict[str, object]]:
    group_records = []
    for stage in sample_stages:
        record = _execute_preprocess_stage(
            stage=stage,
            project_root=project_root,
            env=env,
            resume=resume,
        )
        group_records.append(record)
        if record["status"] == "failed" and not continue_on_error:
            break
    return group_records


def _run_full_length_plan_in_parallel(
    prepared_stages: List[Dict[str, object]],
    project_root: str,
    env: Dict[str, str],
    timing_path: str,
    summary_path: str,
    resume: bool,
    continue_on_error: bool,
    max_parallel_cells: int,
) -> List[Dict[str, object]]:
    sample_stages = [
        stage
        for stage in prepared_stages
        if stage.get("sample_id") is not None
    ]
    global_stages = [
        stage
        for stage in prepared_stages
        if stage.get("sample_id") is None
    ]
    sample_groups = _group_stages_by_sample(sample_stages)
    records_by_index: Dict[int, Dict[str, object]] = {}
    next_group_index = 0
    stop_submitting = False

    with futures.ThreadPoolExecutor(max_workers=max_parallel_cells) as executor:
        running: Dict[futures.Future, List[Dict[str, object]]] = {}

        while (
            next_group_index < len(sample_groups)
            and len(running) < max_parallel_cells
        ):
            sample_group = sample_groups[next_group_index]
            future = executor.submit(
                _run_sample_stage_group,
                sample_stages=sample_group,
                project_root=project_root,
                env=env,
                resume=resume,
                continue_on_error=continue_on_error,
            )
            running[future] = sample_group
            next_group_index += 1

        while running:
            done, _ = futures.wait(
                running.keys(),
                return_when=futures.FIRST_COMPLETED,
            )
            for completed_future in done:
                running.pop(completed_future)
                group_records = completed_future.result()
                for record in group_records:
                    records_by_index[int(record["stage_index"])] = record
                ordered_records = [
                    records_by_index[index]
                    for index in sorted(records_by_index)
                ]
                _write_execution_records(
                    ordered_records,
                    timing_path=timing_path,
                    summary_path=summary_path,
                )
                if (
                    not continue_on_error
                    and any(record["status"] == "failed" for record in group_records)
                ):
                    stop_submitting = True

            while (
                not stop_submitting
                and next_group_index < len(sample_groups)
                and len(running) < max_parallel_cells
            ):
                sample_group = sample_groups[next_group_index]
                future = executor.submit(
                    _run_sample_stage_group,
                    sample_stages=sample_group,
                    project_root=project_root,
                    env=env,
                    resume=resume,
                    continue_on_error=continue_on_error,
                )
                running[future] = sample_group
                next_group_index += 1

    ordered_records = [records_by_index[index] for index in sorted(records_by_index)]
    if stop_submitting and not continue_on_error:
        return ordered_records

    for stage in global_stages:
        record = _execute_preprocess_stage(
            stage=stage,
            project_root=project_root,
            env=env,
            resume=resume,
        )
        records_by_index[int(stage["stage_index"])] = record
        ordered_records = [records_by_index[index] for index in sorted(records_by_index)]
        _write_execution_records(
            ordered_records,
            timing_path=timing_path,
            summary_path=summary_path,
        )
        if record["status"] == "failed" and not continue_on_error:
            break

    return [records_by_index[index] for index in sorted(records_by_index)]


def _parse_featurecounts_file(input_path: str, count_column_start: int):
    with open(input_path, "r", encoding="utf-8") as handle:
        comment_lines = []
        header = None
        data_rows = []
        for line in handle:
            if header is None and line.startswith("#"):
                comment_lines.append(line)
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                continue
            data_rows.append(line.rstrip("\n").split("\t"))

    if header is None:
        raise ValueError(f"No header found in featureCounts file: {input_path}")
    if len(header) <= count_column_start:
        raise ValueError(
            f"featureCounts file has too few columns for grouped output: {input_path}"
        )

    annotation_columns = header[:count_column_start]
    count_columns = header[count_column_start:]
    return comment_lines, header, data_rows, annotation_columns, count_columns


def validate_grouped_featurecounts(
    input_path: str,
    mode: str,
) -> Dict[str, object]:
    if mode not in FEATURECOUNTS_MODE_CONFIG:
        raise ValueError(
            f"Unsupported grouped featureCounts mode: {mode}. "
            f"Supported modes: {sorted(FEATURECOUNTS_MODE_CONFIG)}"
        )

    input_path = os.path.abspath(input_path)
    config = FEATURECOUNTS_MODE_CONFIG[mode]
    _, header, data_rows, annotation_columns, count_columns = _parse_featurecounts_file(
        input_path=input_path,
        count_column_start=config["count_column_start"],
    )

    group_ids = [derive_group_id_from_featurecounts_column(column) for column in count_columns]
    duplicate_group_ids = sorted(
        {
            group_id
            for group_id in group_ids
            if group_ids.count(group_id) > 1
        }
    )
    if duplicate_group_ids:
        raise ValueError(
            f"Grouped featureCounts file contains duplicate group IDs: {duplicate_group_ids}"
        )

    return {
        "mode": mode,
        "input_path": input_path,
        "annotation_column_count": len(annotation_columns),
        "group_count": len(count_columns),
        "group_ids": group_ids,
        "row_count": len(data_rows),
        "header": header,
    }


def derive_group_id_from_featurecounts_column(column_name: str) -> str:
    basename = os.path.basename(column_name)
    if ":" in basename:
        return basename.rsplit(":", 1)[-1]
    return basename


def split_grouped_featurecounts(
    input_path: str,
    output_dir: str,
    mode: str,
    metadata_path: Optional[str] = None,
    cell_id_list_path: Optional[str] = None,
    summary_path: Optional[str] = None,
) -> Dict[str, object]:
    if mode not in FEATURECOUNTS_MODE_CONFIG:
        raise ValueError(
            f"Unsupported split mode: {mode}. "
            f"Supported modes: {sorted(FEATURECOUNTS_MODE_CONFIG)}"
        )
    if metadata_path is not None and cell_id_list_path is not None:
        raise ValueError(
            "Use either metadata_path or cell_id_list_path, not both."
        )

    input_path = os.path.abspath(input_path)
    output_dir = _ensure_directory(output_dir)
    config = FEATURECOUNTS_MODE_CONFIG[mode]
    comment_lines, header, data_rows, annotation_columns, count_columns = _parse_featurecounts_file(
        input_path=input_path,
        count_column_start=config["count_column_start"],
    )

    requested_cells = None
    if metadata_path is not None:
        metadata_df = read_metadata_table(metadata_path)
        requested_cells = metadata_df["cell_ids"]
    elif cell_id_list_path is not None:
        requested_cells = read_cell_id_list(cell_id_list_path)

    column_by_cell_id = {}
    for count_column in count_columns:
        group_id = derive_group_id_from_featurecounts_column(count_column)
        if group_id in column_by_cell_id:
            raise ValueError(
                f"Duplicate grouped featureCounts column for cell '{group_id}': "
                f"{column_by_cell_id[group_id]} and {count_column}"
            )
        column_by_cell_id[group_id] = count_column

    if requested_cells is None:
        requested_cells = sorted(column_by_cell_id)

    written_files = []
    missing_cells = []
    for cell_id in requested_cells:
        count_column = column_by_cell_id.get(cell_id)
        if count_column is None:
            missing_cells.append(cell_id)
            continue

        count_column_index = header.index(count_column)
        output_path = os.path.join(output_dir, f"{cell_id}{config['output_suffix']}")
        with open(output_path, "w", encoding="utf-8") as handle:
            for line in comment_lines:
                handle.write(line)
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(annotation_columns + [count_column])
            for row in data_rows:
                count_value = row[count_column_index]
                if config["drop_zero_rows"]:
                    try:
                        if float(count_value) <= 0:
                            continue
                    except ValueError:
                        pass
                writer.writerow(row[:config["count_column_start"]] + [count_value])
        written_files.append(output_path)

    summary = {
        "mode": mode,
        "input_path": input_path,
        "output_dir": output_dir,
        "metadata_path": (
            os.path.abspath(metadata_path) if metadata_path is not None else None
        ),
        "cell_id_list_path": (
            os.path.abspath(cell_id_list_path)
            if cell_id_list_path is not None else None
        ),
        "requested_cell_count": len(requested_cells),
        "written_file_count": len(written_files),
        "written_files": written_files,
        "missing_cells": missing_cells,
    }
    if summary_path is not None:
        summary_path = os.path.abspath(summary_path)
        os.makedirs(os.path.dirname(summary_path), exist_ok=True)
        with open(summary_path, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")
        summary["summary_path"] = summary_path
    return summary


def stream_cb_to_rg(
    input_stream,
    output_stream,
    barcode_tag: str = "CB",
    rg_tag: str = "RG",
    append_suffix: Optional[str] = None,
    allow_list_path: Optional[str] = None,
) -> Dict[str, int]:
    alignments_seen = 0
    alignments_tagged = 0
    alignments_written = 0
    alignments_filtered_out = 0

    allowed_barcodes = None
    rg_ids = None
    if allow_list_path is not None:
        allowed_barcodes = set(read_cell_id_list(allow_list_path))
        rg_ids = sorted(
            f"{barcode}{append_suffix}" if append_suffix else barcode
            for barcode in allowed_barcodes
        )

    barcode_prefix = f"{barcode_tag}:"
    rg_prefix = f"{rg_tag}:"
    header_lines = []
    existing_rg_ids = set()
    header_flushed = False

    for raw_line in input_stream:
        if raw_line.startswith("@"):
            header_lines.append(raw_line)
            if raw_line.startswith("@RG\t"):
                for field in raw_line.rstrip("\n").split("\t")[1:]:
                    if field.startswith("ID:"):
                        existing_rg_ids.add(field.split(":", 1)[1])
                        break
            continue

        if not header_flushed:
            for header_line in header_lines:
                output_stream.write(header_line)
            if rg_ids is not None:
                for rg_id in rg_ids:
                    if rg_id in existing_rg_ids:
                        continue
                    output_stream.write(
                        f"@RG\tID:{rg_id}\tSM:{rg_id}\tLB:DOLPHIN\tPL:ILLUMINA\n"
                    )
            header_flushed = True

        alignments_seen += 1
        stripped_line = raw_line.rstrip("\n")
        fields = stripped_line.split("\t")

        barcode_value = None
        new_optional_fields = []
        for field in fields[11:]:
            if field.startswith(barcode_prefix):
                parts = field.split(":", 2)
                if len(parts) == 3:
                    barcode_value = parts[2]
                new_optional_fields.append(field)
                continue
            if field.startswith(rg_prefix):
                continue
            new_optional_fields.append(field)

        if allowed_barcodes is not None and barcode_value is None:
            alignments_filtered_out += 1
            continue

        if barcode_value is not None:
            if allowed_barcodes is not None and barcode_value not in allowed_barcodes:
                alignments_filtered_out += 1
                continue
            if append_suffix:
                barcode_value = f"{barcode_value}{append_suffix}"
            new_optional_fields.append(f"{rg_tag}:Z:{barcode_value}")
            alignments_tagged += 1

        new_fields = fields[:11] + new_optional_fields
        output_stream.write("\t".join(new_fields) + "\n")
        alignments_written += 1

    if not header_flushed:
        for header_line in header_lines:
            output_stream.write(header_line)
        if rg_ids is not None:
            for rg_id in rg_ids:
                if rg_id in existing_rg_ids:
                    continue
                output_stream.write(
                    f"@RG\tID:{rg_id}\tSM:{rg_id}\tLB:DOLPHIN\tPL:ILLUMINA\n"
                )

    return {
        "alignments_seen": alignments_seen,
        "alignments_tagged": alignments_tagged,
        "alignments_written": alignments_written,
        "alignments_filtered_out": alignments_filtered_out,
    }


def rewrite_bam_cb_to_rg(
    input_bam_path: str,
    output_bam_path: str,
    barcode_tag: str = "CB",
    rg_tag: str = "RG",
    append_suffix: Optional[str] = None,
    allow_list_path: Optional[str] = None,
    threads: int = 1,
    write_index: bool = False,
) -> Dict[str, object]:
    import pysam

    input_bam_path = os.path.abspath(input_bam_path)
    output_bam_path = os.path.abspath(output_bam_path)
    os.makedirs(os.path.dirname(output_bam_path), exist_ok=True)

    allowed_barcodes = None
    rg_ids = None
    if allow_list_path is not None:
        allowed_barcodes = set(read_cell_id_list(allow_list_path))
        rg_ids = sorted(
            f"{barcode}{append_suffix}" if append_suffix else barcode
            for barcode in allowed_barcodes
        )

    alignments_seen = 0
    alignments_tagged = 0
    alignments_written = 0
    alignments_filtered_out = 0

    with pysam.AlignmentFile(
        input_bam_path,
        "rb",
        threads=max(1, int(threads)),
        check_sq=False,
    ) as bam_in:
        header = bam_in.header.to_dict()
        if rg_ids is not None:
            rg_records = list(header.get("RG", []))
            existing_ids = {
                record.get("ID")
                for record in rg_records
                if record.get("ID") is not None
            }
            for rg_id in rg_ids:
                if rg_id in existing_ids:
                    continue
                rg_records.append(
                    {
                        "ID": rg_id,
                        "SM": rg_id,
                        "LB": "DOLPHIN",
                        "PL": "ILLUMINA",
                    }
                )
            header["RG"] = rg_records

        with pysam.AlignmentFile(
            output_bam_path,
            "wb",
            header=header,
            threads=max(1, int(threads)),
        ) as bam_out:
            for alignment in bam_in.fetch(until_eof=True):
                alignments_seen += 1
                if not alignment.has_tag(barcode_tag):
                    if allowed_barcodes is not None:
                        alignments_filtered_out += 1
                        continue
                    bam_out.write(alignment)
                    alignments_written += 1
                    continue

                barcode_value = alignment.get_tag(barcode_tag)
                if allowed_barcodes is not None and barcode_value not in allowed_barcodes:
                    alignments_filtered_out += 1
                    continue

                rg_value = (
                    f"{barcode_value}{append_suffix}"
                    if append_suffix else barcode_value
                )
                alignment.set_tag(rg_tag, rg_value, value_type="Z")
                bam_out.write(alignment)
                alignments_tagged += 1
                alignments_written += 1

    if write_index:
        pysam.index(output_bam_path)

    return {
        "input_bam_path": input_bam_path,
        "output_bam_path": output_bam_path,
        "alignments_seen": alignments_seen,
        "alignments_tagged": alignments_tagged,
        "alignments_written": alignments_written,
        "alignments_filtered_out": alignments_filtered_out,
        "rg_header_count": len(rg_ids or []),
        "threads": int(threads),
        "write_index": bool(write_index),
    }

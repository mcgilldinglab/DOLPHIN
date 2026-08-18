import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


MODULE_DIR = Path(__file__).resolve().parent
PACKAGE_ROOT = MODULE_DIR.parents[2]
ROOT = Path(os.environ.get("DOLPHIN_AS_ROOT", PACKAGE_ROOT / "alternative_splicing_runs"))
VENDORED_SITE = Path(
    os.environ.get(
        "DOLPHIN_AS_VENDOR_SITE",
        "/mnt/md0/kailu/DOLPHIN_codex/runtime_envs/as_site_manual",
    )
)
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

from find_cell_neighbor import run_find_neighbor
from get_single_bam_reads import run_reads_count
from process_junction_aggregation import run_junction_aggregation
from process_reads_aggregation import run_reads_aggregation
from convert_psi_to_h5ad import run_convert_psi
from convert_random_psi import run_psi_random
from generate_differential_as import run_differential_as


def _cpu_count():
    return os.cpu_count() or 1


def _default_tool(env_name, bundled_path, tool_name):
    env_value = os.environ.get(env_name)
    if env_value:
        return env_value
    if bundled_path and os.path.exists(bundled_path):
        return bundled_path
    resolved = shutil.which(tool_name)
    if resolved:
        return resolved
    return tool_name


def _default_results_root():
    return os.environ.get("DOLPHIN_AS_RESULTS_ROOT", str(ROOT / "results"))


def _default_logs_root():
    return os.environ.get("DOLPHIN_AS_LOGS_ROOT", str(ROOT / "logs"))


def _default_outrigger_work_root():
    return os.environ.get(
        "DOLPHIN_AS_OUTRIGGER_WORK_ROOT",
        os.path.join(tempfile.gettempdir(), "dolphin_as_work"),
    )


def _default_prepared_inputs_root():
    return os.environ.get(
        "DOLPHIN_AS_PREPARED_INPUTS_ROOT",
        os.path.join(tempfile.gettempdir(), "dolphin_as_prepared_inputs"),
    )


def _default_read_count_workers():
    return min(8, max(1, _cpu_count() // 2))


def _default_aggregation_workers():
    return min(4, max(1, _cpu_count() // 8))


def _default_star_threads():
    return min(8, max(1, _cpu_count() // 4))


def _default_star_jobs():
    return min(2, max(1, _cpu_count() // 16))


def _default_outrigger_python():
    return os.environ.get("DOLPHIN_AS_OUTRIGGER_PYTHON", sys.executable)


def _default_outrigger_pythonpath():
    parts = [
        str(MODULE_DIR / "runtime_support" / "pandas_compat"),
        str(MODULE_DIR / "runtime_support" / "outrigger_patched"),
    ]
    if VENDORED_SITE.exists():
        parts.append(str(VENDORED_SITE))
    extra = os.environ.get("DOLPHIN_AS_EXTRA_PYTHONPATH")
    if extra:
        parts.append(extra)
    return os.pathsep.join(parts)


DEFAULTS = {
    "embedding_h5ad": "/mnt/md0/kailu/DOLPHIN_codex/execution_runs/full_length_dolphin_model_v2/data/DOLPHIN_Z.h5ad",
    "metadata_path": "/mnt/md0/kailu/DOLPHIN_codex/data_example/flashseq/metadata.csv",
    "bam_root": "/mnt/md0/kailu/DOLPHIN_codex/execution_runs/full_length_flashseq_clean_v1/outputs/03_exon_star",
    "junction_root": "/mnt/md0/kailu/DOLPHIN_codex/execution_runs/full_length_flashseq_clean_v1/outputs/03_exon_star",
    "bam_file_extension": ".Aligned.sortedByCoord.out.bam",
    "junction_file_extension": ".SJ.out.tab",
    "results_root": _default_results_root(),
    "logs_root": _default_logs_root(),
    "star_binary": _default_tool(
        "DOLPHIN_AS_STAR",
        "/mnt/md0/kailu/DOLPHIN_codex/runtime_envs/dolphin_preprocess_env/bin/STAR",
        "STAR",
    ),
    "samtools_binary": _default_tool(
        "DOLPHIN_AS_SAMTOOLS",
        "/mnt/md0/kailu/DOLPHIN_codex/runtime_envs/dolphin_preprocess_env/bin/samtools",
        "samtools",
    ),
    "gtf_path": "/mnt/md0/kailu/microbiota_fungi_codex/04_human_gtf_star_index/inputs/Homo_sapiens.GRCh38.107.gtf",
    "gffutils_db": "/mnt/md0/kailu/DOLPHIN_codex/AS_run/reference_cache/Homo_sapiens.GRCh38.107.gtf.db",
    "genome_sizes_path": "/mnt/md0/kailu/DOLPHIN_codex/AS_run/reference_cache/grch38_ensembl107.chrom.sizes",
    "fasta_path": "/mnt/md0/kailu/DOLPHIN_codex/reference_inputs/grch38_release107/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
    "star_index_dir": "/mnt/md0/kailu/microbiota_fungi_codex/04_human_gtf_star_index/outputs/star_index_ensembl107",
    "out_name": "full_length_flashseq",
    "neighbor_k": 10,
    "cluster_name": "celltype1",
    "min_event_cells": 10,
    "star_threads": _default_star_threads(),
    "star_jobs": _default_star_jobs(),
    "read_count_workers": _default_read_count_workers(),
    "aggregation_workers": _default_aggregation_workers(),
    "aggregation_mode": "bam",
    "prepared_inputs_root": _default_prepared_inputs_root(),
    "rg_bam_path": None,
    "metadata_normalization": "none",
    "outrigger_python": _default_outrigger_python(),
    "outrigger_pythonpath": _default_outrigger_pythonpath(),
    "outrigger_work_root": _default_outrigger_work_root(),
    "use_gffutils_db": False,
}


def _run(cmd, *, cwd=None, env=None):
    print("$", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, env=env, check=True)


def _write_status(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))


def _resolve_input_file(root: str, sample: str, extension: str):
    root_path = Path(root)
    candidates = [
        root_path / f"{sample}{extension}",
        root_path / sample / f"{sample}{extension}",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Missing source file for {sample}: tried {candidates}")


def _resolve_cell_ids(metadata_path, max_cells=None):
    cell_ids = []
    with open(metadata_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cell_ids.append(row["CB"])
    if max_cells is not None:
        cell_ids = cell_ids[:max_cells]
    return cell_ids


def _normalize_tenx_barcode(value):
    if value is None:
        return value
    text = str(value)
    text = text.split("_", 1)[0]
    text = text.split("-", 1)[0]
    return text


def _prepare_metadata_for_run(metadata_path, embedding_h5ad, output_path, summary_path, normalization_mode="none"):
    import anndata
    import pandas as pd

    df = pd.read_csv(metadata_path, sep="\t")
    if "CB" not in df.columns:
        raise ValueError(f"Metadata file must contain a CB column: {metadata_path}")

    source_rows = int(df.shape[0])
    if normalization_mode == "tenx_barcode":
        df["CB"] = df["CB"].map(_normalize_tenx_barcode)
    elif normalization_mode != "none":
        raise ValueError(f"Unsupported metadata normalization mode: {normalization_mode}")

    value_cols = [col for col in df.columns if col != "CB"]
    dedup_rows = []
    conflicts = []
    for cb, group in df.groupby("CB", sort=False):
        row = {"CB": cb}
        for col in value_cols:
            values = [value for value in group[col].dropna().astype(str).unique().tolist() if value and value != "nan"]
            if len(values) <= 1:
                row[col] = values[0] if values else pd.NA
            else:
                row[col] = values[0]
                conflicts.append({"CB": cb, "column": col, "values": values})
        dedup_rows.append(row)

    prepared = pd.DataFrame(dedup_rows)
    adata = anndata.read_h5ad(embedding_h5ad, backed="r")
    obs_names = list(adata.obs_names)
    aligned = pd.DataFrame({"CB": obs_names}).merge(prepared, on="CB", how="left")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    aligned.to_csv(output_path, sep="\t", index=False)

    summary = {
        "source_metadata_path": str(metadata_path),
        "embedding_h5ad": str(embedding_h5ad),
        "normalization_mode": normalization_mode,
        "source_rows": source_rows,
        "prepared_rows": int(prepared.shape[0]),
        "aligned_rows": int(aligned.shape[0]),
        "conflict_count": len(conflicts),
        "missing_cb_after_alignment": int(aligned["CB"].isna().sum()),
        "missing_values_by_column": {
            col: int(aligned[col].isna().sum())
            for col in value_cols
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2))
    if conflicts:
        pd.DataFrame(conflicts).to_csv(
            output_path.with_suffix(".conflicts.tsv"),
            sep="\t",
            index=False,
        )


def _link_sj_files(star_per_cell_dir, sj_dir):
    if sj_dir.exists() or sj_dir.is_symlink():
        if sj_dir.is_symlink() or sj_dir.is_file():
            sj_dir.unlink()
        else:
            shutil.rmtree(sj_dir)
    sj_dir.mkdir(parents=True, exist_ok=True)
    linked = 0
    for cell_dir in sorted(star_per_cell_dir.iterdir()):
        if not cell_dir.is_dir():
            continue
        if cell_dir.name == sj_dir.name:
            continue
        matches = list(cell_dir.glob("*.aggr.SJ.out.tab"))
        for src in matches:
            dst = sj_dir / src.name
            os.symlink(src.resolve(), dst)
            linked += 1
    return linked


def _sync_tree(src, dst):
    src = Path(src)
    dst = Path(dst)
    if dst.exists():
        shutil.rmtree(dst)
    if src.exists():
        shutil.copytree(src, dst)


def _run_star_single(cell, cell_agg_dir, star_aggr_dir, args):
    cell_out = star_aggr_dir / cell
    cell_out.mkdir(parents=True, exist_ok=True)
    bam_path = cell_agg_dir / f"{cell}.aggr.final.bam"
    prefix = cell_out / f"{cell}.aggr."
    cmd = [
        args.star_binary,
        "--runThreadN",
        str(args.star_threads),
        "--genomeDir",
        args.star_index_dir,
        "--readFilesIn",
        str(bam_path),
        "--readFilesCommand",
        args.samtools_binary,
        "view",
        "-F",
        "0x100",
        "--outSAMtype",
        "BAM",
        "Unsorted",
        "--readFilesType",
        "SAM",
        "SE",
        "--outFileNamePrefix",
        str(prefix),
    ]
    _run(cmd)


def _prepare_outrigger_workdir(final_dir, local_root, out_name):
    final_dir = Path(final_dir)
    local_dir = Path(local_root) / out_name / "outrigger_output"
    local_dir.parent.mkdir(parents=True, exist_ok=True)
    if local_dir.exists():
        shutil.rmtree(local_dir)
    if final_dir.exists():
        shutil.copytree(final_dir, local_dir)
    return local_dir


def _normalize_read_counts_csv(path):
    import pandas as pd

    df = pd.read_csv(path)
    if "sample" in df.columns:
        df["sample"] = (
            df["sample"]
            .str.replace(".Aligned.sortedByCoord.out", "", regex=False)
            .str.replace(".aggr.final", "", regex=False)
        )
    df.to_csv(path, index=False)


def _validate_read_count_samples(metadata_path, read_counts_csv, sample_col="sample", metadata_col="CB", preview=10):
    import pandas as pd

    metadata = pd.read_csv(metadata_path, sep="\t")
    read_counts = pd.read_csv(read_counts_csv)

    metadata_ids = set(metadata[metadata_col].astype(str))
    read_count_ids = set(read_counts[sample_col].astype(str))

    if metadata_ids == read_count_ids:
        return

    overlap = metadata_ids & read_count_ids
    if overlap:
        return

    metadata_only = sorted(metadata_ids - read_count_ids)[:preview]
    read_count_only = sorted(read_count_ids - metadata_ids)[:preview]
    raise ValueError(
        "read_counts sample IDs do not overlap metadata cell IDs. "
        f"metadata_only_example={metadata_only}; "
        f"read_counts_only_example={read_count_only}"
    )


def _write_subset_metadata(metadata_path, cell_ids, out_path):
    keep = set(cell_ids)
    with open(metadata_path, newline="") as src, open(out_path, "w", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        writer = csv.DictWriter(dst, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()
        for row in reader:
            if row["CB"] in keep:
                writer.writerow(row)


def _stage_subset_inputs(cell_ids, bam_root, junction_root, bam_extension, junction_extension, staging_root):
    bam_stage = staging_root / "bam"
    junction_stage = staging_root / "junction"
    bam_stage.mkdir(parents=True, exist_ok=True)
    junction_stage.mkdir(parents=True, exist_ok=True)

    for cell in cell_ids:
        bam_src = _resolve_input_file(bam_root, cell, bam_extension)
        sj_src = _resolve_input_file(junction_root, cell, junction_extension)
        bam_dst = bam_stage / f"{cell}{bam_extension}"
        sj_dst = junction_stage / f"{cell}{junction_extension}"
        if bam_dst.exists() or bam_dst.is_symlink():
            bam_dst.unlink()
        if sj_dst.exists() or sj_dst.is_symlink():
            sj_dst.unlink()
        os.symlink(bam_src, bam_dst)
        os.symlink(sj_src, sj_dst)
    return bam_stage, junction_stage


def _ensure_rg_split_bams(rg_bam_path, output_root, samtools_binary, expected_cells):
    output_root = Path(output_root)
    expected_paths = [output_root / f"{cell}.Aligned.sortedByCoord.out.bam" for cell in expected_cells]
    if expected_paths and all(path.exists() for path in expected_paths):
        return output_root

    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    format_string = str(output_root / "%!.Aligned.sortedByCoord.out.bam")
    _run(
        [
            samtools_binary,
            "split",
            "-d",
            "RG",
            "-M",
            str(max(5000, len(expected_cells) + 100)),
            "-f",
            format_string,
            str(rg_bam_path),
        ]
    )
    missing = [path.name for path in expected_paths if not path.exists()]
    if missing:
        raise RuntimeError(
            f"samtools split finished but {len(missing)} expected per-cell BAMs are missing; "
            f"examples: {missing[:5]}"
        )
    return output_root


def _collect_required_cells_from_neighbors(neighbor_csv, main_cells):
    import pandas as pd

    main_cells = set(main_cells)
    df = pd.read_csv(neighbor_csv)
    df = df[df["main_name"].isin(main_cells)]
    required = set(df["main_name"]).union(set(df["neighbor"]))
    return sorted(required)


def run_pipeline(args):
    root = Path(args.results_root) / args.out_name
    logs_dir = Path(args.logs_root)
    root.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    status_path = logs_dir / f"{args.out_name}.status.json"
    summary_path = logs_dir / f"{args.out_name}.summary.json"

    metadata_prepared = root / "metadata_prepared.tsv"
    metadata_prepare_summary = root / "metadata_preparation_summary.json"
    _prepare_metadata_for_run(
        args.metadata_path,
        args.embedding_h5ad,
        metadata_prepared,
        metadata_prepare_summary,
        normalization_mode=args.metadata_normalization,
    )
    cell_ids = _resolve_cell_ids(metadata_prepared, args.max_cells)
    started = time.time()
    _write_status(
        status_path,
        {
            "current_step": "starting",
            "started_at": started,
            "n_cells_requested": len(cell_ids),
            "out_name": args.out_name,
        },
    )

    # Ensure legacy shell-outs find STAR/samtools by name.
    tool_bin = str(Path(args.samtools_binary).parent)
    os.environ["PATH"] = tool_bin + os.pathsep + os.environ.get("PATH", "")

    metadata_for_run = metadata_prepared
    bam_root_for_run = Path(args.bam_root)
    junction_root_for_run = Path(args.junction_root)
    staging_root = None

    if args.rg_bam_path:
        prepared_bam_key = Path(args.rg_bam_path).name.replace(".bam", "")
        bam_root_for_run = _ensure_rg_split_bams(
            args.rg_bam_path,
            Path(args.prepared_inputs_root) / prepared_bam_key / "bam",
            args.samtools_binary,
            cell_ids,
        )

    if args.max_cells is not None:
        staging_root = root / "staging_inputs"
        staging_root.mkdir(parents=True, exist_ok=True)
        metadata_for_run = staging_root / "metadata_subset.tsv"
        _write_subset_metadata(metadata_prepared, cell_ids, metadata_for_run)

    neighbors_csv = root / f"N_{args.out_name}_{args.neighbor_k}.csv"
    read_counts_csv = root / f"{args.out_name}_read_counts.csv"
    cell_agg_dir = root / "cell_aggregation"
    star_aggr_dir = root / "aggregated_star"
    star_sj_dir = star_aggr_dir / "SJ"
    outrigger_dir = root / "outrigger_output"
    outrigger_work_dir = _prepare_outrigger_workdir(
        outrigger_dir,
        args.outrigger_work_root,
        args.out_name,
    )
    alt_splice_dir = root / "alternative_splicing"

    timings = {}

    # Step 1: neighbors
    if not neighbors_csv.exists():
        t0 = time.time()
        _write_status(status_path, {"current_step": "find_neighbors", "started_at": started})
        run_find_neighbor(
            embedding_data=args.embedding_h5ad,
            out_name=args.out_name,
            N_neighbor=args.neighbor_k,
            out_directory=str(root),
        )
        timings["find_neighbors"] = time.time() - t0

    if args.max_cells is not None:
        required_cells = _collect_required_cells_from_neighbors(neighbors_csv, cell_ids)
        bam_root_for_run, junction_root_for_run = _stage_subset_inputs(
            required_cells,
            str(bam_root_for_run),
            str(junction_root_for_run),
            args.bam_file_extension,
            args.junction_file_extension,
            staging_root,
        )

    # Step 2: read counts
    if not read_counts_csv.exists():
        t0 = time.time()
        _write_status(status_path, {"current_step": "read_counts", "started_at": started})
        run_reads_count(
            out_name=args.out_name,
            bam_file_path=str(bam_root_for_run),
            out_directory=str(root),
            samtools_binary=args.samtools_binary,
            jobs=args.read_count_workers,
        )
        _normalize_read_counts_csv(read_counts_csv)
        timings["read_counts"] = time.time() - t0

    _validate_read_count_samples(metadata_for_run, read_counts_csv)

    # Step 3: aggregation
    if args.aggregation_mode == "bam":
        expected_aggr = [cell_agg_dir / f"{cell}.aggr.final.bam" for cell in cell_ids]
        if not all(path.exists() for path in expected_aggr):
            t0 = time.time()
            _write_status(status_path, {"current_step": "cell_aggregation", "started_at": started})
            run_reads_aggregation(
                metadata_path=str(metadata_for_run),
                bam_file_path=str(bam_root_for_run),
                bam_file_extension=args.bam_file_extension,
                junction_file_path=str(junction_root_for_run),
                junction_file_extension=args.junction_file_extension,
                neighbor_file=str(neighbors_csv),
                read_count_path=str(read_counts_csv),
                N_neighbor=args.neighbor_k,
                out_directory=str(root),
                samtools_binary=args.samtools_binary,
                neighbor_workers=args.aggregation_workers,
            )
            timings["cell_aggregation"] = time.time() - t0

        # Step 4: STAR on aggregated BAM
        missing_star = []
        for cell in cell_ids:
            out_prefix_dir = star_aggr_dir / cell
            sj_path = out_prefix_dir / f"{cell}.aggr.SJ.out.tab"
            if not sj_path.exists():
                missing_star.append(cell)

        if missing_star:
            t0 = time.time()
            _write_status(status_path, {"current_step": "star_aggregate_alignment", "started_at": started})
            star_aggr_dir.mkdir(parents=True, exist_ok=True)
            if args.star_jobs > 1:
                with ThreadPoolExecutor(max_workers=args.star_jobs) as executor:
                    futures = [
                        executor.submit(_run_star_single, cell, cell_agg_dir, star_aggr_dir, args)
                        for cell in missing_star
                    ]
                    for future in futures:
                        future.result()
            else:
                for cell in missing_star:
                    _run_star_single(cell, cell_agg_dir, star_aggr_dir, args)
            linked = _link_sj_files(star_aggr_dir, star_sj_dir)
            timings["star_aggregate_alignment"] = time.time() - t0
            timings["linked_sj_files"] = linked
        else:
            _link_sj_files(star_aggr_dir, star_sj_dir)
    elif args.aggregation_mode == "junction":
        missing_sj = []
        for cell in cell_ids:
            out_prefix_dir = star_aggr_dir / cell
            sj_path = out_prefix_dir / f"{cell}.aggr.SJ.out.tab"
            if not sj_path.exists():
                missing_sj.append(cell)
        if missing_sj:
            t0 = time.time()
            _write_status(status_path, {"current_step": "junction_aggregation", "started_at": started})
            run_junction_aggregation(
                metadata_path=str(metadata_for_run),
                junction_file_path=str(junction_root_for_run),
                junction_file_extension=args.junction_file_extension,
                neighbor_file=str(neighbors_csv),
                read_count_path=str(read_counts_csv),
                N_neighbor=args.neighbor_k,
                out_directory=str(root),
                gtf_path=args.gtf_path,
                fasta_path=args.fasta_path,
            )
            linked = _link_sj_files(star_aggr_dir, star_sj_dir)
            timings["junction_aggregation"] = time.time() - t0
            timings["linked_sj_files"] = linked
        else:
            _link_sj_files(star_aggr_dir, star_sj_dir)
    else:
        raise ValueError(f"Unsupported aggregation mode: {args.aggregation_mode}")

    # Step 5: Outrigger
    outrigger_index_complete = (
        (outrigger_dir / "index" / "se" / "events.csv").exists()
        and (outrigger_dir / "index" / "mxe" / "events.csv").exists()
    )
    outrigger_validate_complete = (
        (outrigger_dir / "index" / "se" / "validated" / "events.csv").exists()
        and (outrigger_dir / "index" / "mxe" / "validated" / "events.csv").exists()
    )
    outrigger_summary = outrigger_dir / "psi" / "outrigger_summary.csv"
    outrigger_psi_complete = outrigger_summary.exists()

    outrigger_work_dir.mkdir(parents=True, exist_ok=True)
    sj_files = sorted(str(path) for path in star_sj_dir.glob("*.SJ.out.tab"))
    if not sj_files:
        raise RuntimeError("No aggregated SJ.out.tab files found for Outrigger.")

    env = os.environ.copy()
    env["PYTHONPATH"] = args.outrigger_pythonpath + os.pathsep + env.get("PYTHONPATH", "")

    if not outrigger_index_complete:
        t0 = time.time()
        _write_status(status_path, {"current_step": "outrigger_index", "started_at": started})
        _run(
            [
                args.outrigger_python,
                "-m",
                "outrigger.commandline",
                "index",
                "--output",
                str(outrigger_work_dir),
                "--sj-out-tab",
                *sj_files,
                *(
                    ["--gffutils-db", args.gffutils_db]
                    if args.use_gffutils_db and args.gffutils_db and os.path.exists(args.gffutils_db)
                    else ["--gtf-filename", args.gtf_path]
                ),
            ],
            cwd=str(outrigger_work_dir),
            env=env,
        )
        _sync_tree(outrigger_work_dir, outrigger_dir)
        timings["outrigger_index"] = time.time() - t0

    if not outrigger_validate_complete:
        t0 = time.time()
        _write_status(status_path, {"current_step": "outrigger_validate", "started_at": started})
        _run(
            [
                args.outrigger_python,
                "-m",
                "outrigger.commandline",
                "validate",
                "--output",
                str(outrigger_work_dir),
                "--genome",
                args.genome_sizes_path,
                "--fasta",
                args.fasta_path,
            ],
            cwd=str(outrigger_work_dir),
            env=env,
        )
        _sync_tree(outrigger_work_dir, outrigger_dir)
        timings["outrigger_validate"] = time.time() - t0

    if not outrigger_psi_complete:
        t0 = time.time()
        _write_status(status_path, {"current_step": "outrigger_psi", "started_at": started})
        _run(
            [
                args.outrigger_python,
                "-m",
                "outrigger.commandline",
                "psi",
                "--output",
                str(outrigger_work_dir),
            ],
            cwd=str(outrigger_work_dir),
            env=env,
        )
        _sync_tree(outrigger_work_dir, outrigger_dir)
        timings["outrigger_psi"] = time.time() - t0

    # Step 6: PSI h5ad
    psi_h5ad = alt_splice_dir / f"{args.out_name}_PSI.h5ad"
    if not psi_h5ad.exists():
        t0 = time.time()
        _write_status(status_path, {"current_step": "convert_psi", "started_at": started})
        run_convert_psi(
            metadata_path=str(metadata_for_run),
            outrigger_path=str(outrigger_dir),
            out_name=args.out_name,
            out_directory=str(root),
        )
        timings["convert_psi"] = time.time() - t0

    # Step 7: random PSI
    psi_random_h5ad = alt_splice_dir / f"{args.out_name}_PSI_random.h5ad"
    if not psi_random_h5ad.exists():
        t0 = time.time()
        _write_status(status_path, {"current_step": "psi_random", "started_at": started})
        run_psi_random(
            outrigger_psi_data=str(psi_h5ad),
            out_name=args.out_name,
            out_directory=str(root),
            seed_num=args.seed_num,
        )
        timings["psi_random"] = time.time() - t0

    # Step 8: differential AS
    psi_das_h5ad = alt_splice_dir / f"{args.out_name}_PSI_DAS.h5ad"
    if not psi_das_h5ad.exists():
        t0 = time.time()
        _write_status(status_path, {"current_step": "differential_as", "started_at": started})
        run_differential_as(
            outrigger_psi_data=str(psi_h5ad),
            out_name=args.out_name,
            cluster_name=args.cluster_name,
            out_directory=str(root),
            n_cell=args.min_event_cells,
        )
        timings["differential_as"] = time.time() - t0

    finished = time.time()
    summary = {
        "status": "completed",
        "started_at": started,
        "finished_at": finished,
        "duration_seconds": finished - started,
        "out_root": str(root),
        "n_cells_requested": len(cell_ids),
        "aggregation_mode": args.aggregation_mode,
        "timings": timings,
        "outputs": {
            "neighbors_csv": str(neighbors_csv),
            "metadata_prepared": str(metadata_prepared),
            "metadata_prepare_summary": str(metadata_prepare_summary),
            "read_counts_csv": str(read_counts_csv),
            "prepared_bam_root": str(bam_root_for_run),
            "cell_aggregation_dir": str(cell_agg_dir),
            "aggregated_star_dir": str(star_aggr_dir),
            "outrigger_dir": str(outrigger_dir),
            "psi_h5ad": str(psi_h5ad),
            "psi_random_h5ad": str(psi_random_h5ad),
            "psi_das_h5ad": str(psi_das_h5ad),
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2))
    _write_status(status_path, {"current_step": "completed", "started_at": started, "finished_at": finished})
    print(json.dumps(summary, indent=2))


def build_parser(default_overrides=None):
    defaults = dict(DEFAULTS)
    if default_overrides:
        defaults.update(default_overrides)

    parser = argparse.ArgumentParser(description="Run the finalized DOLPHIN alternative splicing pipeline.")
    parser.add_argument("--embedding-h5ad", default=defaults["embedding_h5ad"])
    parser.add_argument("--metadata-path", default=defaults["metadata_path"])
    parser.add_argument("--bam-root", default=defaults["bam_root"])
    parser.add_argument("--junction-root", default=defaults["junction_root"])
    parser.add_argument("--bam-file-extension", default=defaults["bam_file_extension"])
    parser.add_argument("--junction-file-extension", default=defaults["junction_file_extension"])
    parser.add_argument("--results-root", default=defaults["results_root"])
    parser.add_argument("--logs-root", default=defaults["logs_root"])
    parser.add_argument("--star-binary", default=defaults["star_binary"])
    parser.add_argument("--samtools-binary", default=defaults["samtools_binary"])
    parser.add_argument("--gtf-path", default=defaults["gtf_path"])
    parser.add_argument("--gffutils-db", default=defaults["gffutils_db"])
    parser.add_argument("--genome-sizes-path", default=defaults["genome_sizes_path"])
    parser.add_argument("--fasta-path", default=defaults["fasta_path"])
    parser.add_argument("--star-index-dir", default=defaults["star_index_dir"])
    parser.add_argument("--out-name", default=defaults["out_name"])
    parser.add_argument("--neighbor-k", type=int, default=defaults["neighbor_k"])
    parser.add_argument("--cluster-name", default=defaults["cluster_name"])
    parser.add_argument("--min-event-cells", type=int, default=defaults["min_event_cells"])
    parser.add_argument("--star-threads", type=int, default=defaults["star_threads"])
    parser.add_argument("--star-jobs", type=int, default=defaults["star_jobs"])
    parser.add_argument("--read-count-workers", type=int, default=defaults["read_count_workers"])
    parser.add_argument("--aggregation-workers", type=int, default=defaults["aggregation_workers"])
    parser.add_argument("--aggregation-mode", choices=["bam", "junction"], default=defaults["aggregation_mode"])
    parser.add_argument("--prepared-inputs-root", default=defaults["prepared_inputs_root"])
    parser.add_argument("--rg-bam-path", default=defaults["rg_bam_path"])
    parser.add_argument(
        "--metadata-normalization",
        choices=["none", "tenx_barcode"],
        default=defaults["metadata_normalization"],
    )
    parser.add_argument("--outrigger-python", default=defaults["outrigger_python"])
    parser.add_argument("--outrigger-pythonpath", default=defaults["outrigger_pythonpath"])
    parser.add_argument("--outrigger-work-root", default=defaults["outrigger_work_root"])
    parser.add_argument("--use-gffutils-db", action="store_true", default=defaults["use_gffutils_db"])
    parser.add_argument("--seed-num", type=int, default=0)
    parser.add_argument("--max-cells", type=int)
    return parser


def main(default_overrides=None):
    run_pipeline(build_parser(default_overrides=default_overrides).parse_args())


if __name__ == "__main__":
    main()

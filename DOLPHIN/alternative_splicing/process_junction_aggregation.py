import csv
import os
from collections import Counter
from functools import lru_cache
from pathlib import Path

import pandas as pd
import pysam
from tqdm import tqdm


def _resolve_source_file(root: str, sample: str, extension: str) -> Path:
    root_path = Path(root)
    candidates = [
        root_path / f"{sample}{extension}",
        root_path / sample / f"{sample}{extension}",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Missing source file for {sample}: tried {candidates}")


def _strand_to_code(value):
    if value == "+":
        return 1
    if value == "-":
        return 2
    return 0


def _code_to_strand(value):
    if value == 1:
        return "+"
    if value == 2:
        return "-"
    return None


def _parse_gtf_attributes(raw):
    values = {}
    for chunk in raw.strip().split(";"):
        chunk = chunk.strip()
        if not chunk or " " not in chunk:
            continue
        key, value = chunk.split(" ", 1)
        values[key] = value.strip().strip('"')
    return values


def _iter_gtf_introns(gtf_path):
    transcript_exons = {}
    with open(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            chrom = parts[0]
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue
            strand = parts[6]
            attrs = _parse_gtf_attributes(parts[8])
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue
            transcript_exons.setdefault((chrom, strand, transcript_id), []).append((start, end))

    for (chrom, strand, _transcript_id), exons in transcript_exons.items():
        exons = sorted(exons)
        for left, right in zip(exons, exons[1:]):
            intron_start = left[1] + 1
            intron_end = right[0] - 1
            if intron_start <= intron_end:
                yield chrom, intron_start, intron_end, strand


@lru_cache(maxsize=4)
def _load_annotated_introns(gtf_path):
    introns = {}
    for chrom, start, end, strand in _iter_gtf_introns(gtf_path):
        key = (str(chrom), int(start), int(end))
        introns.setdefault(key, set()).add(_strand_to_code(strand))
    return introns


@lru_cache(maxsize=2)
def _open_fasta(fasta_path):
    return pysam.FastaFile(fasta_path)


def _compute_motif_code(chrom, start, end, fasta):
    donor = fasta.fetch(str(chrom), int(start) - 1, int(start) + 1).upper()
    acceptor = fasta.fetch(str(chrom), int(end) - 2, int(end)).upper()
    motif = f"{donor}/{acceptor}"
    if motif == "GT/AG":
        return 1
    if motif == "CT/AC":
        return 2
    if motif == "GC/AG":
        return 3
    if motif == "CT/GC":
        return 4
    if motif == "AT/AC":
        return 5
    if motif == "GT/AT":
        return 6
    return 0


def _strand_from_motif_code(motif_code):
    if motif_code in {1, 3, 5}:
        return 1
    if motif_code in {2, 4, 6}:
        return 2
    return 0


def _backfill_junction_annotation(chrom, start, end, fallback_annotated, gtf_introns, fasta):
    key = (str(chrom), int(start), int(end))
    motif_code = _compute_motif_code(chrom, start, end, fasta)
    motif_strand = _strand_from_motif_code(motif_code)

    intron_strands = gtf_introns.get(key)
    if intron_strands:
        annotated = 1
        if len(intron_strands) == 1:
            strand_code = next(iter(intron_strands))
        else:
            strand_code = motif_strand
    else:
        annotated = int(fallback_annotated)
        strand_code = motif_strand

    return {
        "strand": int(strand_code),
        "motif": int(motif_code),
        "annotated": int(annotated),
    }


def _init_star_row(chrom, start, end, strand=0, motif=0, annotated=0, uniq=0.0, multi=0.0, max_overhang=0):
    return {
        "chrom": str(chrom),
        "start": int(start),
        "end": int(end),
        "strand": int(strand),
        "motif": int(motif),
        "annotated": int(annotated),
        "uniq": float(uniq),
        "multi": float(multi),
        "max_overhang": int(max_overhang),
    }


def _load_sample_junction_rows(
    sample,
    junction_file_path,
    junction_file_extension,
    junction_cache,
    gtf_introns=None,
    fasta=None,
):
    cached = junction_cache.get(sample)
    if cached is not None:
        return cached

    src = _resolve_source_file(junction_file_path, sample, junction_file_extension)
    rows = {}

    if junction_file_extension.endswith(".jcounts"):
        with open(src, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                junction_cache[sample] = rows
                return rows
            count_field = reader.fieldnames[-1]
            for row in reader:
                chrom = row.get("Chr_SP1")
                chrom2 = row.get("Chr_SP2")
                start = row.get("Location_SP1")
                end = row.get("Location_SP2")
                if not chrom or not start or not end:
                    continue
                if chrom2 and chrom2 != chrom:
                    continue
                try:
                    count = int(float(row[count_field]))
                    start = int(float(start))
                    end = int(float(end))
                except (TypeError, ValueError):
                    continue
                if count <= 0:
                    continue
                # featureCounts .jcounts stores exon-flanking coordinates,
                # while STAR SJ.out.tab expects the intron coordinates.
                # Align the 10x junction route with the BAM/STAR route by
                # shifting to STAR-compatible intron boundaries.
                start += 1
                end -= 1
                if start > end:
                    continue
                key = (str(chrom), start, end)
                if key in rows:
                    rows[key]["uniq"] += count
                    continue
                status = row.get("Status")
                strand_code = _strand_to_code(row.get("Strand_SP1")) or _strand_to_code(row.get("Strand_SP2"))
                annotated = 0 if str(status).upper() == "NOVEL" else 1
                motif = 0
                if gtf_introns is not None and fasta is not None:
                    annotation = _backfill_junction_annotation(
                        chrom=chrom,
                        start=start,
                        end=end,
                        fallback_annotated=annotated,
                        gtf_introns=gtf_introns,
                        fasta=fasta,
                    )
                    strand_code = annotation["strand"]
                    motif = annotation["motif"]
                    annotated = annotation["annotated"]
                rows[key] = _init_star_row(
                    chrom=chrom,
                    start=start,
                    end=end,
                    strand=strand_code,
                    motif=motif,
                    annotated=annotated,
                    uniq=count,
                    multi=0,
                    max_overhang=0,
                )
    else:
        with open(src, newline="") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                try:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    strand = int(parts[3])
                    motif = int(parts[4])
                    annotated = int(parts[5])
                    uniq = int(parts[6])
                    multi = int(parts[7])
                    max_overhang = int(parts[8])
                except ValueError:
                    continue
                if uniq <= 0 and multi <= 0:
                    continue
                key = (chrom, start, end)
                rows[key] = _init_star_row(
                    chrom=chrom,
                    start=start,
                    end=end,
                    strand=strand,
                    motif=motif,
                    annotated=annotated,
                    uniq=uniq,
                    multi=multi,
                    max_overhang=max_overhang,
                )

    junction_cache[sample] = rows
    return rows


def _load_keep_junctions(
    neighbors,
    junction_file_path,
    junction_file_extension,
    threshold,
    junction_cache,
    gtf_introns=None,
    fasta=None,
):
    seen_order = []
    seen_once = set()
    counts = Counter()

    for neighbor in neighbors:
        per_file = _load_sample_junction_rows(
            neighbor,
            junction_file_path,
            junction_file_extension,
            junction_cache,
            gtf_introns=gtf_introns,
            fasta=fasta,
        )
        for key in per_file:
            counts[key] += 1
            if key not in seen_once:
                seen_once.add(key)
                seen_order.append(key)

    return [key for key in seen_order if counts[key] >= threshold]


def _clone_row(row):
    return {
        "chrom": row["chrom"],
        "start": row["start"],
        "end": row["end"],
        "strand": row["strand"],
        "motif": row["motif"],
        "annotated": row["annotated"],
        "uniq": float(row["uniq"]),
        "multi": float(row["multi"]),
        "max_overhang": row["max_overhang"],
    }


def _merge_scaled_row(target_row, source_row, scale):
    target_row["uniq"] += source_row["uniq"] * scale
    target_row["multi"] += source_row["multi"] * scale
    target_row["max_overhang"] = max(target_row["max_overhang"], source_row["max_overhang"])
    if target_row["strand"] == 0 and source_row["strand"] != 0:
        target_row["strand"] = source_row["strand"]
    if target_row["motif"] == 0 and source_row["motif"] != 0:
        target_row["motif"] = source_row["motif"]
    target_row["annotated"] = max(target_row["annotated"], source_row["annotated"])


def _write_sj_out_tab(rows, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        for key in sorted(rows):
            row = rows[key]
            uniq = int(round(row["uniq"]))
            multi = int(round(row["multi"]))
            if uniq <= 0 and multi <= 0:
                continue
            handle.write(
                "\t".join(
                    [
                        row["chrom"],
                        str(row["start"]),
                        str(row["end"]),
                        str(row["strand"]),
                        str(row["motif"]),
                        str(row["annotated"]),
                        str(uniq),
                        str(multi),
                        str(row["max_overhang"]),
                    ]
                )
                + "\n"
            )


def run_junction_aggregation(
    metadata_path: str,
    junction_file_path: str,
    junction_file_extension: str,
    neighbor_file: str,
    read_count_path: str,
    N_neighbor: int = 10,
    out_directory: str = "./",
    gtf_path: str | None = None,
    fasta_path: str | None = None,
):
    """
    Aggregate per-cell junction summaries directly into STAR `SJ.out.tab`-compatible files.

    This keeps the neighborhood / majority-voting logic aligned with the BAM-based pipeline:
    - target cell contributes all of its own junction counts without filtering
    - neighbor cells contribute only junctions passing majority voting
    - neighbor contributions are scaled to the target read depth

    Supported inputs:
    - STAR `SJ.out.tab`
    - DOLPHIN per-cell `.jcounts`
    """

    pd_aggr = pd.read_csv(neighbor_file)
    pd_single_size = pd.read_csv(read_count_path)
    sample_list = list(pd.read_csv(metadata_path, sep="\t")["CB"])
    read_count_map = dict(zip(pd_single_size["sample"], pd_single_size["num_seqs"]))
    neighbor_map = {
        key: list(group["neighbor"])
        for key, group in pd_aggr.groupby("main_name", sort=False)
    }

    star_aggr_dir = Path(out_directory) / "aggregated_star"
    star_aggr_dir.mkdir(parents=True, exist_ok=True)
    junction_cache = {}
    gtf_introns = _load_annotated_introns(gtf_path) if gtf_path and junction_file_extension.endswith(".jcounts") else None
    fasta = _open_fasta(fasta_path) if fasta_path and junction_file_extension.endswith(".jcounts") else None

    for target in tqdm(sample_list):
        target_size = read_count_map[target]
        neighbors = neighbor_map[target]
        keep_keys = set(
            _load_keep_junctions(
                neighbors,
                junction_file_path,
                junction_file_extension,
                int(N_neighbor / 2),
                junction_cache,
                gtf_introns=gtf_introns,
                fasta=fasta,
            )
        )

        target_rows = _load_sample_junction_rows(
            target,
            junction_file_path,
            junction_file_extension,
            junction_cache,
            gtf_introns=gtf_introns,
            fasta=fasta,
        )
        aggregated = {key: _clone_row(row) for key, row in target_rows.items()}

        for neighbor in neighbors:
            if neighbor == target:
                continue
            neighbor_rows = _load_sample_junction_rows(
                neighbor,
                junction_file_path,
                junction_file_extension,
                junction_cache,
                gtf_introns=gtf_introns,
                fasta=fasta,
            )
            neighbor_size = read_count_map[neighbor]
            if neighbor_size <= 0:
                continue
            scale = target_size / neighbor_size
            for key, row in neighbor_rows.items():
                if keep_keys and key not in keep_keys:
                    continue
                if key not in aggregated:
                    aggregated[key] = _init_star_row(
                        chrom=row["chrom"],
                        start=row["start"],
                        end=row["end"],
                        strand=row["strand"],
                        motif=row["motif"],
                        annotated=row["annotated"],
                        uniq=0,
                        multi=0,
                        max_overhang=row["max_overhang"],
                    )
                _merge_scaled_row(aggregated[key], row, scale)

        cell_out = star_aggr_dir / target
        output_path = cell_out / f"{target}.aggr.SJ.out.tab"
        _write_sj_out_tab(aggregated, output_path)

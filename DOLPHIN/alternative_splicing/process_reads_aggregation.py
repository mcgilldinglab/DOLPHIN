import csv
import os
import shlex
import shutil
import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
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


def _run_shell(command: str):
    subprocess.run(command, shell=True, executable="/bin/sh", check=True)


def _load_sample_junctions(sample, junction_file_path, junction_file_extension, junction_cache):
    cached = junction_cache.get(sample)
    if cached is not None:
        return cached

    src = _resolve_source_file(junction_file_path, sample, junction_file_extension)
    per_file = set()
    if junction_file_extension.endswith(".jcounts"):
        with open(src, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                junction_cache[sample] = per_file
                return per_file
            count_field = reader.fieldnames[-1]
            for row in reader:
                chrom = row.get("Chr_SP1")
                chrom2 = row.get("Chr_SP2")
                first_base = row.get("Location_SP1")
                last_base = row.get("Location_SP2")
                if not chrom or not first_base or not last_base:
                    continue
                if chrom2 and chrom2 != chrom:
                    continue
                try:
                    count = int(float(row[count_field]))
                    first_base = str(int(float(first_base)))
                    last_base = str(int(float(last_base)))
                except (TypeError, ValueError):
                    continue
                if count <= 0:
                    continue
                per_file.add((chrom, first_base, last_base))
    else:
        with open(src, newline="") as handle:
            for line in handle:
                row = line.rstrip("\n").split("\t")
                if len(row) < 3:
                    continue
                per_file.add((row[0], row[1], row[2]))

    junction_cache[sample] = per_file
    return per_file


def _load_keep_junctions(neighbors, junction_file_path, junction_file_extension, threshold, junction_cache):
    seen_order = []
    seen_once = set()
    counts = Counter()

    for neighbor in neighbors:
        per_file = _load_sample_junctions(
            neighbor,
            junction_file_path,
            junction_file_extension,
            junction_cache,
        )
        for key in per_file:
            counts[key] += 1
            if key not in seen_once:
                seen_once.add(key)
                seen_order.append(key)

    return [key for key in seen_order if counts[key] >= threshold]


def _write_keep_junction_bed(keys, path: Path):
    with open(path, "w", newline="") as handle:
        for chrom, first_base, last_base in keys:
            handle.write(f"{chrom}\t{first_base}\t{last_base}\n")


def _write_majority_spliced_bam(source_bam: Path, keep_bed: Path, out_bam: Path, samtools_binary: str):
    source_q = shlex.quote(str(source_bam))
    keep_q = shlex.quote(str(keep_bed))
    out_q = shlex.quote(str(out_bam))
    sam_q = shlex.quote(str(samtools_binary))
    _run_shell(
        f"{sam_q} view -h {source_q} | "
        "awk '$0 ~ /^@/ || $6 ~ /N/' | "
        f"{sam_q} view -b -L {keep_q} - > {out_q}"
    )


def _prepare_neighbor_bam(
    sample,
    target,
    target_size,
    read_count_map,
    bam_file_path,
    bam_file_extension,
    keep_any,
    keep_bed,
    target_temp_dir,
    samtools_binary,
):
    _n_seq = read_count_map[sample]
    source_bam = _resolve_source_file(bam_file_path, sample, bam_file_extension)
    norm_bam = target_temp_dir / f"{sample}.norm.bam"
    if norm_bam.exists() or norm_bam.is_symlink():
        norm_bam.unlink()
    if _n_seq == target_size:
        os.symlink(source_bam, norm_bam)
    elif _n_seq < target_size:
        _cat_self_n = int(target_size / _n_seq)
        if _cat_self_n == 1:
            _add_seq_perct = (target_size - _n_seq) / _n_seq
        else:
            _add_seq_perct = (target_size - _n_seq * _cat_self_n) / _n_seq
        sample_bam = target_temp_dir / f"{sample}.sample.bam"
        with open(sample_bam, "wb") as out_f:
            subprocess.run(
                [samtools_binary, "view", "-b", "-s", str(_add_seq_perct), str(source_bam)],
                stdout=out_f,
                check=True,
            )
        merge_inputs = [str(source_bam)] * _cat_self_n + [str(sample_bam)]
        subprocess.run([samtools_binary, "merge", str(norm_bam), *merge_inputs], check=True)
        sample_bam.unlink(missing_ok=True)
    else:
        _keep_seq_perct = target_size / _n_seq
        with open(norm_bam, "wb") as out_f:
            subprocess.run(
                [samtools_binary, "view", "-b", "-s", str(_keep_seq_perct), str(source_bam)],
                stdout=out_f,
                check=True,
            )

    if sample != target and keep_any:
        mj_bam = target_temp_dir / f"{sample}.mj.junction.norm.bam"
        _write_majority_spliced_bam(
            norm_bam,
            keep_bed,
            mj_bam,
            samtools_binary,
        )

def run_reads_aggregation(
    metadata_path: str,
    bam_file_path: str,
    bam_file_extension: str,
    junction_file_path: str,
    junction_file_extension: str,
    neighbor_file: str,
    read_count_path: str,
    N_neighbor: int = 10,
    out_directory: str = "./",
    samtools_binary: str = "samtools",
    neighbor_workers: int = 1,
):  
    """
    Aggregate single-cell BAM files by incorporating reads from neighbor cells, 
    applying junction majority voting and read count normalization.

    This function performs the following operations:
    1. Reads metadata, neighbor, and read count files.
    2. For each target cell:
    - Identifies frequent junctions using majority voting across neighbors (Junctions are retained only if they appear in at least half of the neighbors). 
    - Normalizes neighbor BAM files to match the target read count (via up/downsampling).
    - Filters junction reads to retain only those appearing frequently.
    - Aggregates filtered reads from neighbors and unfiltered reads from the target cell.
    3. Outputs a final merged BAM file per target cell.

    Parameters
    ----------
    metadata_path : str
        Path to the metadata file.
    bam_file_path : str
        Directory containing the BAM files generated by STAR.
    bam_file_extension : str
        Suffix of BAM files after the sample name. 
        Example: for "SRR18379095.std.Aligned.sortedByCoord.out.bam", use ".std.Aligned.sortedByCoord.out.bam".
    junction_file_path : str
        Directory containing junction files. Supported formats:
        - STAR `SJ.out.tab`
        - DOLPHIN per-cell `.jcounts`
    junction_file_extension : str
        Suffix of junction files after the sample name.
        Example: for "SRR18379095.std.SJ.out.tab", use ".std.SJ.out.tab".
    neighbor_file : str
        CSV file specifying neighbor relationships. Must include 'main_name' and 'neighbor' columns.
    read_count_path : str
        Path to a CSV file with two columns: 'sample' and 'num_seqs', representing read counts for each cell.
    N_neighbor : int, optional
        Number of neighbors per target cell. Default is 10.
    out_directory : str, optional
        Output directory to save results. Default is the current directory.
    samtools_binary : str, optional
        samtools executable path. Default is "samtools".
    neighbor_workers : int, optional
        Number of concurrent neighbor normalization/filter workers per target cell.
        Default is 1.

    Returns
    -------
    None
        For each cell in the metadata save final aggregated BAM
        <out_directory>/cell_aggregation/<cell_id>.aggr.final.bam
    """
    pd_aggr = pd.read_csv(neighbor_file)
    pd_single_size = pd.read_csv(read_count_path)
    sample_list = list(pd.read_csv(metadata_path, sep='\t')["CB"])
    read_count_map = dict(zip(pd_single_size["sample"], pd_single_size["num_seqs"]))
    neighbor_map = {
        key: list(group["neighbor"])
        for key, group in pd_aggr.groupby("main_name", sort=False)
    }

    final_out_dir = os.path.join(out_directory, "cell_aggregation")
    os.makedirs(final_out_dir, exist_ok=True)
    temp_out_dir = os.path.join(final_out_dir, "temp")
    os.makedirs(temp_out_dir, exist_ok=True)
    junction_cache = {}
    
    for target in tqdm(sample_list):
        target_size = read_count_map[target]
        _neighbor = neighbor_map[target]
        target_temp_dir = Path(temp_out_dir) / target
        target_temp_dir.mkdir(parents=True, exist_ok=True)
        '''
        Majority voting: find the frequent junction reads
        '''
        keep_keys = _load_keep_junctions(
            _neighbor,
            junction_file_path,
            junction_file_extension,
            int(N_neighbor / 2),
            junction_cache,
        )
        keep_bed = target_temp_dir / "keep_junction.bed"
        _write_keep_junction_bed(keep_keys, keep_bed)
        keep_any = bool(keep_keys)
        '''
        Bam file batch size normalization
        '''
        if neighbor_workers > 1:
            with ThreadPoolExecutor(max_workers=neighbor_workers) as executor:
                futures = [
                    executor.submit(
                        _prepare_neighbor_bam,
                        _n,
                        target,
                        target_size,
                        read_count_map,
                        bam_file_path,
                        bam_file_extension,
                        keep_any,
                        keep_bed,
                        target_temp_dir,
                        samtools_binary,
                    )
                    for _n in _neighbor
                ]
                for future in futures:
                    future.result()
        else:
            for _n in _neighbor:
                _prepare_neighbor_bam(
                    _n,
                    target,
                    target_size,
                    read_count_map,
                    bam_file_path,
                    bam_file_extension,
                    keep_any,
                    keep_bed,
                    target_temp_dir,
                    samtools_binary,
                )
        '''
        Final concate all normalized bam files    
        '''     
        # Final concate all normalized fastq files
        mj_bams = sorted(str(path) for path in target_temp_dir.glob("*.mj.junction.norm.bam"))
        final_bam = Path(final_out_dir) / f"{target}.aggr.final.bam"
        merge_cmd = [samtools_binary, "merge", "-f", str(final_bam), *mj_bams, str(target_temp_dir / f"{target}.norm.bam")]
        subprocess.run(merge_cmd, check=True)
        shutil.rmtree(target_temp_dir)

import glob
import json
import math
import os
import shutil
import subprocess
from time import perf_counter

import numpy as np
import pandas as pd
from intervaltree import Interval, IntervalTree
from tqdm import tqdm

from DOLPHIN.preprocess import gtfpy

pd.set_option('display.max_columns',500)
pd.set_option('display.max_rows',100)


def _init_timing_log(timing_log_path):
    if timing_log_path is None:
        return

    os.makedirs(os.path.dirname(timing_log_path), exist_ok=True)
    with open(timing_log_path, "w", encoding="utf-8") as handle:
        handle.write("step\tseconds\tminutes\textra_json\n")


def _append_timing_record(timing_log_path, step, seconds, extra=None):
    if timing_log_path is None:
        return

    extra = extra or {}
    with open(timing_log_path, "a", encoding="utf-8") as handle:
        handle.write(
            f"{step}\t{seconds:.6f}\t{seconds / 60.0:.6f}\t"
            f"{json.dumps(extra, sort_keys=True)}\n"
        )
    print(f"[Timing] {step}: {seconds:.2f}s")


def _run_timed_step(step_name, timing_log_path, func, *args, extra=None, **kwargs):
    start_time = perf_counter()
    result = func(*args, **kwargs)
    elapsed = perf_counter() - start_time
    _append_timing_record(timing_log_path, step_name, elapsed, extra=extra)
    return result


def _resolve_executable(executable_name):
    if os.path.sep in executable_name:
        executable_path = os.path.abspath(executable_name)
        if os.path.isfile(executable_path) and os.access(executable_path, os.X_OK):
            return executable_path
        raise FileNotFoundError(f"Executable not found or not executable: {executable_name}")

    resolved_path = shutil.which(executable_name)
    if resolved_path is None:
        raise FileNotFoundError(f"Could not find executable on PATH: {executable_name}")
    return resolved_path


def _extract_gtf_attribute(attribute_series, field_name):
    return attribute_series.str.extract(fr'{field_name} "([^"]+)"', expand=False)

def prepare_exon_gtf(input_gtf_path, output_dir="./"):
    """
    Load an Ensembl GTF file and extract exon-level annotations with unique start/end per gene.

    Parameters
    ----------
    input_gtf_path : str
        Path to the original Ensembl .gtf file.
    output_dir : str, optional
        Directory to save intermediate results (default: './dolphin_exon_gtf/').

    Returns
    -------
    df_exon_nodup : pandas.DataFrame
        Filtered exon annotation table with duplicates (same gene_id, start, end) removed.
    """

    # 1. Ensure output directory exists
    os.makedirs(os.path.join(output_dir, "dolphin_exon_gtf"), exist_ok=True)

    print(f"[Step] Reading GTF file from: {input_gtf_path}")

    # 2. Load GTF and only parse exon rows and needed attributes.
    df_gtf = gtfpy.readGTF(input_gtf_path)
    df_exon = df_gtf[df_gtf["feature"] == "exon"].copy()
    print(
        f"[Status] Loaded {df_gtf.shape[0]} total entries and kept "
        f"{df_exon.shape[0]} exon entries."
    )

    if df_exon.empty:
        raise ValueError("No exon entries were found in the input GTF.")

    attribute_cols = [
        "gene_id",
        "gene_name",
        "exon_number",
        "gene_version",
        "gene_source",
        "gene_biotype",
        "gene_type",
    ]
    for col in attribute_cols:
        df_exon[col] = _extract_gtf_attribute(df_exon["attribute"], col)

    missing_required = [
        col for col in ("gene_id", "exon_number")
        if df_exon[col].isna().all()
    ]
    if missing_required:
        raise ValueError(
            f"Input GTF is missing required attributes: {missing_required}."
        )

    final_cols = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "gene_id",
        "gene_name",
        "exon_number",
        "gene_version",
        "gene_source",
        "gene_biotype",
        "gene_type",
    ]
    df_exon = df_exon[final_cols]

    # 4. Sort and convert coordinates
    df_exon[["start", "end"]] = df_exon[["start", "end"]].apply(pd.to_numeric)
    df_exon = df_exon.sort_values(
        by=["gene_id", "start", "end"],
        ascending=[True, True, True],
    )

    # 5. Remove duplicate exons (same gene, start, end)
    df_exon_nodup = df_exon.drop_duplicates(subset=["gene_id", "start", "end"], keep="first")

    print(f"[Status] Removed duplicates: {df_exon_nodup.shape[0]} unique exon entries remain.")
    
    return df_exon_nodup

def exon_uniq(
    df_exon_nodup: pd.DataFrame,
    gene: str = None,
    gene_df: pd.DataFrame = None,
) -> pd.DataFrame:
    """
    Merge overlapping exons for a single gene using interval trees.

    Parameters
    ----------
    df_exon_nodup : pandas.DataFrame
        DataFrame containing all exons (from `prepare_exon_gtf`), including gene IDs and coordinates.
    gene : str
        The gene ID whose exons will be processed.

    Returns
    -------
    pandas.DataFrame
        A cleaned exon DataFrame for the given gene, where overlapping exons are merged,
        exon coordinates are updated, and exon numbers are reindexed. Exons that are
        invalid or cannot be matched to any merged region are excluded.
    """

    if gene_df is not None:
        _df_exon = gene_df.reset_index(drop=True).copy()
    else:
        _df_exon = df_exon_nodup.loc[df_exon_nodup["gene_id"] == gene].reset_index(drop=True).copy()

    if gene is None and not _df_exon.empty:
        gene = _df_exon["gene_id"].iloc[0]

    # Return empty DataFrame if the gene has no exons
    if _df_exon.empty:
        print(f"[Warning] Gene {gene} has no exon entries.")
        return pd.DataFrame(columns=df_exon_nodup.columns)

    # Filter out invalid exons where start >= end
    _df_exon = _df_exon[_df_exon["start"] < _df_exon["end"]].reset_index(drop=True)
    if _df_exon.empty:
        print(f"[Warning] Gene {gene} has only invalid exon coordinates (start >= end).")
        return pd.DataFrame(columns=df_exon_nodup.columns)

    # Create IntervalTree from valid exon intervals
    try:
        valid_intervals = [
            Interval(row["start"], row["end"])
            for _, row in _df_exon.iterrows()
            if row["start"] < row["end"]
        ]
        tree = IntervalTree(valid_intervals)
        tree.merge_overlaps()
    except Exception as e:
        print(f"[Error] IntervalTree failed for gene {gene}: {e}")
        return pd.DataFrame(columns=df_exon_nodup.columns)

    # Extract merged, sorted non-overlapping intervals
    merged_intervals = sorted([(iv.begin, iv.end) for iv in tree])

    if not merged_intervals:
        print(f"[Warning] Gene {gene} has no merged exon intervals after processing.")
        return pd.DataFrame(columns=df_exon_nodup.columns)

    template = _df_exon.iloc[0].copy()
    rows = []
    for exon_number, (start, end) in enumerate(merged_intervals, start=1):
        row = template.copy()
        row["start"] = int(start)
        row["end"] = int(end)
        row["exon_number"] = exon_number
        rows.append(row)

    return pd.DataFrame(rows, columns=_df_exon.columns)


def process_exons_in_memory(df_exon_nodup, flush_every=5000, timing_log_path=None):
    """
    Process gene groups fully in memory and only write the final outputs.

    Compared with the legacy temp-pickle workflow, this avoids repeated disk
    writes and repeated full-dataframe scans for every gene.
    """
    grouped = df_exon_nodup.groupby("gene_id", sort=False)
    total_genes = grouped.ngroups

    completed_chunks = []
    current_chunk = []
    chunk_start = perf_counter()
    chunk_gene_count = 0
    chunk_index = 0

    for idx, (gene, gene_df) in enumerate(
        tqdm(grouped, total=total_genes, desc="Processing all genes"),
        start=1,
    ):
        processed = exon_uniq(df_exon_nodup, gene=gene, gene_df=gene_df)
        if not processed.empty:
            current_chunk.append(processed)
        chunk_gene_count += 1

        if idx % flush_every == 0 or idx == total_genes:
            chunk_df = (
                pd.concat(current_chunk, ignore_index=True)
                if current_chunk
                else pd.DataFrame(columns=df_exon_nodup.columns)
            )
            completed_chunks.append(chunk_df)
            _append_timing_record(
                timing_log_path,
                f"process_exons_in_memory.chunk_{chunk_index}",
                perf_counter() - chunk_start,
                extra={
                    "genes_in_chunk": chunk_gene_count,
                    "processed_genes": int(idx),
                    "rows_written": int(chunk_df.shape[0]),
                    "flush_every": int(flush_every),
                    "total_genes": int(total_genes),
                },
            )
            chunk_index += 1
            current_chunk = []
            chunk_start = perf_counter()
            chunk_gene_count = 0

    if not completed_chunks:
        return pd.DataFrame(columns=df_exon_nodup.columns)

    return pd.concat(completed_chunks, ignore_index=True)

def save_by_batch(df_exon_nodup, save_num=10000, output_dir="./", timing_log_path=None):
    """
    Process exon annotations for each gene in batches and save results as serialized .pkl files.

    This function applies `exon_uniq()` to each gene in the input DataFrame and saves the processed
    exon data in batches. Each batch contains up to `save_num` genes and is written to a pickle file.
    A log file is generated to record processing status and potential errors.

    Parameters
    ----------
    df_exon_nodup : pandas.DataFrame
        DataFrame containing filtered exon annotations (typically from `prepare_exon_gtf`).
    save_num : int, optional
        Number of genes to include per output batch file (default is 10,000).
    output_dir : str, optional
        Path to the output directory where batch `.pkl` files and the log file will be stored
        (default is "./").

    Returns
    -------
    None
        This function writes intermediate results to disk but does not return any object.
    """

    print(f"[Step] Start processing and saving exons by batch...")

    # 1. Create temp folder
    temp_dir = os.path.join(output_dir, "dolphin_exon_gtf", "temp")
    os.makedirs(temp_dir, exist_ok=True)

    # delete all .pkl 文件
    for f in glob.glob(os.path.join(temp_dir, "df_exon_gtf_*.pkl")):
        os.remove(f)

    # delete all log 文件
    if os.path.exists(os.path.join(temp_dir,"process_log.txt")):
        os.remove("temp/process_log.txt")
    
    # 2. Gene list
    gene_list = df_exon_nodup["gene_id"].unique().tolist()
    total_genes = len(gene_list)

    # 3. Initialize list for all batch DataFrames (optional, for further processing)
    list_of_df = []
    log_path = os.path.join(temp_dir, "process_log.txt")

    # 4. Track current batch number
    current_batch = 0
    df_out = pd.DataFrame()

    # 5. Process all genes with a unified tqdm progress bar
    batch_start = perf_counter()
    batch_gene_count = 0
    for i, gene in enumerate(tqdm(gene_list, total=total_genes, desc="Processing all genes")):
        try:
            _temp = exon_uniq(df_exon_nodup, gene)
            df_out = pd.concat([df_out, _temp], ignore_index=True)
            batch_gene_count += 1

            with open(log_path, "a") as f:
                f.write(f"Batch {current_batch}, Processed gene: {gene}\n")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Batch {current_batch}, Error processing gene {gene}: {e}\n")
            continue  # skip gene with error

        # Save this batch every save_num genes or at the end
        if (i + 1) % save_num == 0 or (i + 1 == total_genes):
            rows_in_batch = int(df_out.shape[0])
            output_path = os.path.join(temp_dir, f"df_exon_gtf_{current_batch}.pkl")
            df_out.to_pickle(output_path)
            list_of_df.append(df_out)
            _append_timing_record(
                timing_log_path,
                f"save_by_batch.batch_{current_batch}",
                perf_counter() - batch_start,
                extra={
                    "genes_in_batch": batch_gene_count,
                    "processed_genes": int(i + 1),
                    "rows_written": rows_in_batch,
                    "save_num": int(save_num),
                    "total_genes": int(total_genes),
                },
            )
            df_out = pd.DataFrame()  # reset for next batch
            current_batch += 1
            batch_start = perf_counter()
            batch_gene_count = 0
    
    print(f"[Done] Finished saving all exon batches.")


def combine_saved_batches(folder="./", prefix="df_exon_gtf_"):
    """
    Combine multiple saved exon batch files into a single concatenated DataFrame.

    This function reads all `.pkl` files in the specified folder that start with the given prefix,
    concatenates them in order, and returns a single DataFrame containing all exon records.

    Parameters
    ----------
    folder : str, optional
        Directory where batch `.pkl` files are stored.
        Default is "./", which typically points to the parent of "dolphin_exon_gtf/temp".
    prefix : str, optional
        Filename prefix used to identify batch `.pkl` files.
        Default is ``df_exon_gtf_``.
        
    Returns
    -------
    pandas.DataFrame
        A single DataFrame containing concatenated exon entries from all batch files.
        The rows are ordered according to batch and file sorting.
    """

    temp_dir = os.path.join(folder, "dolphin_exon_gtf", "temp")
    # List all matching .pkl files in the folder
    pkl_files = [os.path.join(temp_dir, f) for f in os.listdir(temp_dir)
                 if f.startswith(prefix) and f.endswith(".pkl")]

    # Sort to ensure batch order is respected
    pkl_files.sort()

    if not pkl_files:
        raise FileNotFoundError(f"No .pkl files found in '{temp_dir}' with prefix '{prefix}'.")

    # Read and concatenate
    dataframes = [pd.read_pickle(f) for f in pkl_files]
    combined_df = pd.concat(dataframes, ignore_index=True)

    print(f"Successfully combined {len(pkl_files)} files into a single DataFrame with {combined_df.shape[0]} rows.")
    return combined_df

def check_exon_overlap(gtf_df, expected_gene_list=None):
    """
    Check for overlapping adjacent exon intervals within each gene.

    This function checks whether any exons within the same gene have overlapping intervals,
    based on their `start` and `end` positions. Optionally, it compares the set of gene IDs
    in the provided DataFrame with an expected list to detect any missing or extra genes.

    Parameters
    ----------
    gtf_df : pandas.DataFrame
        A DataFrame containing exon annotations with at least the columns: 'gene_id', 'start', and 'end'.
    expected_gene_list : list of str, optional
        A list of expected gene IDs used to validate that all genes were processed and included in `gtf_df`.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing exon entries that overlap with their adjacent exons within the same gene.
        The result may be empty if no overlaps are detected.
    """
    
    df_check = gtf_df.copy()

    # Sort by gene and start for correct order
    df_check = df_check.sort_values(by=["gene_id", "start"]).reset_index(drop=True)

    # Get the start of the next exon within each gene
    df_check['_next_start'] = df_check.groupby('gene_id')['start'].shift(-1)

    df_check['_overlap'] = df_check['_next_start'].isna() | (
        df_check['_next_start'] >= df_check['end']
    )

    # Filter exons that overlap with the next one
    overlap_issues = df_check[df_check['_overlap'] == False]
    print(f"Found {overlap_issues.shape[0]} overlapping exon entries.")

    # If expected gene list is provided, validate gene count
    if expected_gene_list is not None:
        unique_genes_in_df = set(gtf_df["gene_id"].unique())
        expected_genes = set(expected_gene_list)
        missing_genes = expected_genes - unique_genes_in_df
        extra_genes = unique_genes_in_df - expected_genes

        if not missing_genes and not extra_genes:
            print(f"All {len(expected_genes)} expected genes are present in the merged DataFrame.")
        else:
            print(f"Gene count mismatch:")
            if missing_genes:
                print(f"  - {len(missing_genes)} gene(s) missing: {list(missing_genes)[:5]} ...")
            if extra_genes:
                print(f"  - {len(extra_genes)} unexpected gene(s): {list(extra_genes)[:5]} ...")

    return overlap_issues

def save_gtf_outputs(gtf_df, output_dir="./", base_name="dolphin.exon"):
    """
    Save the final exon DataFrame to both GTF and Pickle formats.

    This function writes the given exon annotation table to two output files:
    one in standard GTF format, and the other as a serialized Python pickle (.pkl).

    Parameters
    ----------
    gtf_df : pandas.DataFrame
        The exon annotation DataFrame to be saved.
    output_dir : str, optional
        Directory where the output files will be saved (default is the current directory).
    base_name : str, optional
        Filename prefix used for both output files (default is "dolphin.exon").

    Returns
    -------
        <output_dir>/dolphin_exon_gtf/<base_name>.gtf : GTF-format annotation file
        <output_dir>/dolphin_exon_gtf/<base_name>.pkl : Pickle-serialized DataFrame
    """
    reference_dir = os.path.join(os.path.abspath(output_dir), "dolphin_exon_gtf")
    os.makedirs(reference_dir, exist_ok=True)

    # Keep the original step0 naming and also emit the alignment-facing name used in step1.
    gtf_path = os.path.join(reference_dir, f"{base_name}.gtf")
    pkl_path = os.path.join(reference_dir, f"{base_name}.pkl")
    alignment_gtf_path = os.path.join(reference_dir, "dolphin_exon_gtf.gtf")

    gtf_output_df = gtf_df.copy()
    if "gene_id" in gtf_output_df.columns and "transcript_id" not in gtf_output_df.columns:
        # The DOLPHIN merged-exon reference is gene-centric. Provide a stable
        # transcript_id so featureCounts -J can parse the GTF contract.
        gtf_output_df["transcript_id"] = gtf_output_df["gene_id"]

    fixed_cols = [col for col in gtf_output_df.columns[:8]]
    preferred_attr_order = [
        "gene_id",
        "transcript_id",
        "gene_name",
        "exon_number",
        "gene_version",
        "gene_source",
        "gene_biotype",
        "gene_type",
    ]
    remaining_attr_cols = [
        col for col in gtf_output_df.columns[8:]
        if col not in preferred_attr_order
    ]
    ordered_attr_cols = [
        col for col in preferred_attr_order if col in gtf_output_df.columns
    ] + remaining_attr_cols
    gtf_output_df = gtf_output_df[fixed_cols + ordered_attr_cols]

    # Write to GTF
    gtfpy.writeGTF(gtf_output_df, gtf_path)

    # Write to Pickle
    gtf_df.to_pickle(pkl_path)

    if os.path.abspath(alignment_gtf_path) != os.path.abspath(gtf_path):
        shutil.copyfile(gtf_path, alignment_gtf_path)

    print(f"GTF file saved to: {gtf_path}")
    print(f"Pickle file saved to: {pkl_path}")
    print(f"Alignment GTF file saved to: {alignment_gtf_path}")

    return gtf_path, pkl_path, alignment_gtf_path

def generate_nonoverlapping_exons(
    input_gtf_path: str,
    output_dir: str = "./",
    batch_size: int = 10000,
    timing_log_path=None,
):
    """
    End-to-end pipeline to process an Ensembl GTF file and generate non-overlapping exons per gene.

    This function performs the following steps:
    1. Load and filter exon features from a GTF file.
    2. Remove duplicate exons (by gene_id, start, end).
    3. Process each gene to merge overlapping exons using IntervalTree.
    4. Save intermediate results in batches.
    5. Combine all batches into a final exon DataFrame.
    6. Optionally check for residual overlaps.
    7. Save the final results in GTF and Pickle formats.

    Parameters
    ----------
    input_gtf_path : str
        Path to the input Ensembl-format GTF file.
    output_dir : str
        Directory to save intermediate and final output files.
    batch_size : int
        Number of genes to process and save per batch (default: 10000).

    Returns
    -------
    gtf_all : pd.DataFrame
        Final merged and cleaned exon annotation table.
    overlap_issues : pd.DataFrame
        DataFrame of overlapping exons detected post-processing (if any).
    """
    
    total_start = perf_counter()

    # Step 1: Load exon entries from GTF and remove duplicates
    df_exon_nodup = _run_timed_step(
        "generate_nonoverlapping_exons.prepare_exon_gtf",
        timing_log_path,
        prepare_exon_gtf,
        input_gtf_path,
        output_dir=output_dir,
        extra={"input_gtf_path": os.path.abspath(input_gtf_path)},
    )

    # Step 2: Process exons in memory
    def _process_exons():
        return process_exons_in_memory(
            df_exon_nodup,
            flush_every=batch_size,
            timing_log_path=timing_log_path,
        )

    gtf_all = _run_timed_step(
        "generate_nonoverlapping_exons.process_exons_in_memory",
        timing_log_path,
        _process_exons,
        extra={"flush_every": int(batch_size)},
    )

    # Step 3: Check for residual overlaps
    overlap_issues = _run_timed_step(
        "generate_nonoverlapping_exons.check_exon_overlap",
        timing_log_path,
        check_exon_overlap,
        gtf_all,
        expected_gene_list=df_exon_nodup["gene_id"].unique().tolist(),
    )

    # Step 4: Save final GTF and pickle files
    _run_timed_step(
        "generate_nonoverlapping_exons.save_gtf_outputs",
        timing_log_path,
        save_gtf_outputs,
        gtf_all,
        output_dir=output_dir,
    )

    _append_timing_record(
        timing_log_path,
        "generate_nonoverlapping_exons.total",
        perf_counter() - total_start,
        extra={
            "batch_size": int(batch_size),
            "gene_count": int(df_exon_nodup["gene_id"].nunique()),
            "input_gtf_path": os.path.abspath(input_gtf_path),
        },
    )
    
    print(f"[Success] Exon GTF processing pipeline completed.")

    return gtf_all, overlap_issues


def build_star_genome_index(
    genome_fasta_path,
    alignment_gtf_path,
    output_dir,
    star_executable="STAR",
    run_thread_n=16,
    genome_sa_sparse_d=None,
    sjdb_overhang=None,
):
    """
    Build a STAR genome index using the DOLPHIN alignment GTF.

    Parameters
    ----------
    genome_fasta_path : str
        Path to the genome FASTA file.
    alignment_gtf_path : str
        Path to the DOLPHIN alignment GTF file.
    output_dir : str
        Directory where STAR index files will be written.
    star_executable : str, optional
        STAR executable name or absolute path.
    run_thread_n : int, optional
        Number of threads to pass to STAR.
    genome_sa_sparse_d : int, optional
        Optional STAR ``--genomeSAsparseD`` value.
    sjdb_overhang : int, optional
        Optional STAR ``--sjdbOverhang`` value.

    Returns
    -------
    dict
        Manifest entries for the STAR index.
    """
    genome_fasta_path = os.path.abspath(genome_fasta_path)
    alignment_gtf_path = os.path.abspath(alignment_gtf_path)
    output_dir = os.path.abspath(output_dir)

    if not os.path.isfile(genome_fasta_path):
        raise FileNotFoundError(f"Genome FASTA not found: {genome_fasta_path}")
    if not os.path.isfile(alignment_gtf_path):
        raise FileNotFoundError(f"Alignment GTF not found: {alignment_gtf_path}")

    star_path = _resolve_executable(star_executable)
    os.makedirs(output_dir, exist_ok=True)

    command = [
        star_path,
        "--runMode",
        "genomeGenerate",
        "--runThreadN",
        str(int(run_thread_n)),
        "--genomeDir",
        output_dir,
        "--genomeFastaFiles",
        genome_fasta_path,
        "--sjdbGTFfile",
        alignment_gtf_path,
    ]
    if genome_sa_sparse_d is not None:
        command.extend(["--genomeSAsparseD", str(int(genome_sa_sparse_d))])
    if sjdb_overhang is not None:
        command.extend(["--sjdbOverhang", str(int(sjdb_overhang))])

    log_path = os.path.join(output_dir, "star_genome_generate.log")
    print(f"[Step] Building STAR genome index in: {output_dir}")

    with open(log_path, "w", encoding="utf-8") as log_handle:
        log_handle.write("Command:\n")
        log_handle.write(" ".join(command) + "\n\n")
        try:
            subprocess.run(
                command,
                check=True,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
            )
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"STAR genomeGenerate failed. See log: {log_path}"
            ) from exc

    return {
        "star_executable": star_path,
        "star_index_dir": output_dir,
        "star_genome_generate_log_path": log_path,
        "run_thread_n": int(run_thread_n),
        "genome_sa_sparse_d": (
            int(genome_sa_sparse_d) if genome_sa_sparse_d is not None else None
        ),
        "sjdb_overhang": int(sjdb_overhang) if sjdb_overhang is not None else None,
    }


def prepare_reference_bundle(
    input_gtf_path: str,
    output_dir: str = "./",
    batch_size: int = 10000,
    genome_fasta_path: str = None,
    star_executable: str = "STAR",
    star_threads: int = 16,
    star_index_dir: str = None,
    genome_sa_sparse_d: int = None,
    sjdb_overhang: int = None,
):
    """
    Generate all reference files needed by downstream DOLPHIN preprocessing.

    This wraps the exon-level GTF generation and adjacency metadata generation
    into a single step so users can prepare the alignment and graph references
    from one input Ensembl GTF.

    Parameters
    ----------
    input_gtf_path : str
        Path to the input Ensembl-format GTF file.
    output_dir : str, optional
        Parent directory that will contain the ``dolphin_exon_gtf`` folder.
    batch_size : int, optional
        Number of genes to process per intermediate batch.
    genome_fasta_path : str, optional
        Path to the genome FASTA file. If provided, STAR genomeGenerate will
        be run after the DOLPHIN GTF bundle is created.
    star_executable : str, optional
        STAR executable name or absolute path.
    star_threads : int, optional
        Number of threads to use for STAR genomeGenerate.
    star_index_dir : str, optional
        Output directory for the STAR index. Defaults to
        ``<reference_dir>/star_index``.
    genome_sa_sparse_d : int, optional
        Optional STAR ``--genomeSAsparseD`` value.
    sjdb_overhang : int, optional
        Optional STAR ``--sjdbOverhang`` value.

    Returns
    -------
    dict
        A manifest describing the generated reference files.
    """
    from .generate_adj_index import generate_adj_index_table, generate_adj_metadata_table

    reference_dir = os.path.join(os.path.abspath(output_dir), "dolphin_exon_gtf")
    timing_log_path = os.path.join(reference_dir, "step_timing.tsv")
    _init_timing_log(timing_log_path)

    total_start = perf_counter()

    gtf_all, overlap_issues = generate_nonoverlapping_exons(
        input_gtf_path=input_gtf_path,
        output_dir=output_dir,
        batch_size=batch_size,
        timing_log_path=timing_log_path,
    )

    exon_gtf = os.path.join(reference_dir, "dolphin.exon.gtf")
    alignment_gtf = os.path.join(reference_dir, "dolphin_exon_gtf.gtf")
    exon_pkl = os.path.join(reference_dir, "dolphin.exon.pkl")

    _run_timed_step(
        "prepare_reference_bundle.generate_adj_index_table",
        timing_log_path,
        generate_adj_index_table,
        exon_pkl,
        output_dir=reference_dir,
    )
    _run_timed_step(
        "prepare_reference_bundle.generate_adj_metadata_table",
        timing_log_path,
        generate_adj_metadata_table,
        exon_pkl,
        output_dir=reference_dir,
    )

    star_index_manifest = None
    if genome_fasta_path is not None:
        resolved_star_index_dir = (
            os.path.abspath(star_index_dir)
            if star_index_dir is not None
            else os.path.join(reference_dir, "star_index")
        )
        star_index_manifest = _run_timed_step(
            "prepare_reference_bundle.build_star_genome_index",
            timing_log_path,
            build_star_genome_index,
            genome_fasta_path=genome_fasta_path,
            alignment_gtf_path=alignment_gtf,
            output_dir=resolved_star_index_dir,
            star_executable=star_executable,
            run_thread_n=star_threads,
            genome_sa_sparse_d=genome_sa_sparse_d,
            sjdb_overhang=sjdb_overhang,
        )

    total_seconds = perf_counter() - total_start
    _append_timing_record(
        timing_log_path,
        "prepare_reference_bundle.total",
        total_seconds,
        extra={
            "batch_size": int(batch_size),
            "input_gtf_path": os.path.abspath(input_gtf_path),
        },
    )

    manifest = {
        "input_gtf_path": os.path.abspath(input_gtf_path),
        "reference_dir": reference_dir,
        "exon_gtf_path": exon_gtf,
        "alignment_gtf_path": alignment_gtf,
        "exon_pkl_path": exon_pkl,
        "adj_index_path": os.path.join(reference_dir, "dolphin_adj_index.csv"),
        "adj_metadata_path": os.path.join(reference_dir, "dolphin_adj_metadata_table.csv"),
        "gene_metadata_path": os.path.join(reference_dir, "dolphin_gene_meta.csv"),
        "overlap_issue_count": int(overlap_issues.shape[0]),
        "gene_count": int(gtf_all["gene_id"].nunique()),
        "exon_count": int(gtf_all.shape[0]),
        "timing_log_path": timing_log_path,
        "total_seconds": total_seconds,
        "genome_fasta_path": (
            os.path.abspath(genome_fasta_path)
            if genome_fasta_path is not None
            else None
        ),
        "star_index": star_index_manifest,
    }

    manifest_path = os.path.join(reference_dir, "reference_manifest.json")
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)

    print(f"Reference manifest saved to: {manifest_path}")

    return manifest

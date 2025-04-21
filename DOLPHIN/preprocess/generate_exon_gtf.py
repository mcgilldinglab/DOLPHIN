import pandas as pd
from DOLPHIN.preprocess import gtfpy
import numpy as np
import math
from tqdm import tqdm 
from intervaltree import Interval, IntervalTree
import glob, os

pd.set_option('display.max_columns',500)
pd.set_option('display.max_rows',100)

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

    # 2. Load GTF and parse
    df_gtf = gtfpy.readGTF(input_gtf_path)
    GTFpa = gtfpy.parseGTF(df_gtf)
    
    print(f"[Status] GTF loaded and parsed with {GTFpa.shape[0]} total entries.")

    # 3. Filter exon features and keep necessary columns
    df_exon = GTFpa[GTFpa["feature"] == "exon"][[
        'seqname', 'source', 'feature', 'start', 'end', 'score',
        'strand', 'frame', 'gene_id', 'gene_version', 'gene_name',
        'gene_source', 'gene_biotype', 'exon_number'
    ]]

    # 4. Sort and convert coordinates
    df_exon[["start", "end"]] = df_exon[["start", "end"]].apply(pd.to_numeric)
    df_exon = df_exon.sort_values(by=["seqname", "start", "end"], ascending=[True, True, True])

    # 5. Remove duplicate exons (same gene, start, end)
    df_exon_nodup = df_exon.drop_duplicates(subset=["gene_id", "start", "end"], keep="first")

    print(f"[Status] Removed duplicates: {df_exon_nodup.shape[0]} unique exon entries remain.")
    
    return df_exon_nodup

# def exon_uniq(df_exon_nodup, gene):
#     # Subset exons for the given gene
#     _df_exon = df_exon_nodup.loc[df_exon_nodup["gene_id"] == gene].reset_index(drop=True).copy()

#     # Create interval tree from exon start/end
#     # tree = IntervalTree(Interval(row["start"], row["end"]) for _, row in _df_exon.iterrows())
#     tree = IntervalTree(
#         Interval(row["start"], row["end"])
#         for _, row in _df_exon.iterrows()
#         if row["start"] < row["end"]
#     )
    
#     # Merge overlapping exons
#     tree.merge_overlaps()
    
#     # Get sorted list of unique intervals
#     merged_intervals = sorted([(iv.begin, iv.end) for iv in tree])

#     # Assign merged coordinates to exons
#     _df_exon["_start"] = math.nan
#     _df_exon["_end"] = math.nan

#     for k in range(_df_exon.shape[0]):
#         for start, end in merged_intervals:
#             if start <= _df_exon.at[k, "start"] and end >= _df_exon.at[k, "end"]:
#                 _df_exon.at[k, "_start"] = start
#                 _df_exon.at[k, "_end"] = end
#                 break  # stop at first match

#     # Warn if any exon did not get matched
#     for i in range(_df_exon.shape[0]):
#         if math.isnan(_df_exon.at[i, '_start']):
#             print(f"Attention: {gene}, start={_df_exon.at[i,'start']}, end={_df_exon.at[i,'end']} was not assigned.")

#     # Keep only unique exons based on merged coordinates
#     _df_exon_out = _df_exon.drop_duplicates(subset=["_start", "_end"], keep="first").copy()
#     _df_exon_out["start"] = _df_exon_out["_start"].astype(int)
#     _df_exon_out["end"] = _df_exon_out["_end"].astype(int)
#     _df_exon_out.drop(columns=["_start", "_end"], inplace=True)

#     # Re-index exon numbers
#     _df_exon_out = _df_exon_out.reset_index(drop=True)
#     _df_exon_out["exon_number"] = _df_exon_out.index + 1

#     return _df_exon_out

def exon_uniq(df_exon_nodup: pd.DataFrame, gene: str) -> pd.DataFrame:
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

    # Subset exons for the given gene
    _df_exon = df_exon_nodup.loc[df_exon_nodup["gene_id"] == gene].reset_index(drop=True).copy()

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

    # Initialize columns for matched merged coordinates
    _df_exon["_start"] = math.nan
    _df_exon["_end"] = math.nan

    # Assign each exon to the first merged interval that fully contains it
    for k in range(_df_exon.shape[0]):
        for start, end in merged_intervals:
            if start <= _df_exon.at[k, "start"] and end >= _df_exon.at[k, "end"]:
                _df_exon.at[k, "_start"] = start
                _df_exon.at[k, "_end"] = end
                break  # stop after the first match

    # Warn about unmatched exons
    unmatched = _df_exon[_df_exon["_start"].isna()]
    for _, row in unmatched.iterrows():
        print(f"[Warning] Gene {gene}: exon (start={row['start']}, end={row['end']}) was not assigned to any merged region.")

    # Drop exons that failed to match any merged interval
    _df_exon_out = _df_exon.dropna(subset=["_start", "_end"]).copy()

    # If no exons remain after matching, return empty
    if _df_exon_out.empty:
        print(f"[Warning] Gene {gene} has no exons remaining after merging.")
        return pd.DataFrame(columns=df_exon_nodup.columns)

    # Remove duplicate merged intervals, keeping only one exon per interval
    _df_exon_out = _df_exon_out.drop_duplicates(subset=["_start", "_end"], keep="first").copy()

    # Convert merged coordinates to integer and replace original coordinates
    _df_exon_out["start"] = _df_exon_out["_start"].astype(int)
    _df_exon_out["end"] = _df_exon_out["_end"].astype(int)
    _df_exon_out.drop(columns=["_start", "_end"], inplace=True)

    # Re-index exon numbers
    _df_exon_out = _df_exon_out.reset_index(drop=True)
    _df_exon_out["exon_number"] = _df_exon_out.index + 1

    return _df_exon_out

def save_by_batch(df_exon_nodup, save_num=10000, output_dir="./"):
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
    for i, gene in enumerate(tqdm(gene_list, total=total_genes, desc="Processing all genes")):
        try:
            _temp = exon_uniq(df_exon_nodup, gene)
            df_out = pd.concat([df_out, _temp], ignore_index=True)

            with open(log_path, "a") as f:
                f.write(f"Batch {current_batch}, Processed gene: {gene}\n")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Batch {current_batch}, Error processing gene {gene}: {e}\n")
            continue  # skip gene with error

        # Save this batch every save_num genes or at the end
        if (i + 1) % save_num == 0 or (i + 1 == total_genes):
            output_path = os.path.join(temp_dir, f"df_exon_gtf_{current_batch}.pkl")
            df_out.to_pickle(output_path)
            list_of_df.append(df_out)
            df_out = pd.DataFrame()  # reset for next batch
            current_batch += 1
    
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

    # Initialize overlap check column
    df_check['_overlap'] = True

    # Check if next exon starts after or at the end of current exon
    for i in range(df_check.shape[0]):
        if math.isnan(df_check.loc[i, '_next_start']):
            df_check.loc[i, '_overlap'] = True  # Last exon of gene
        else:
            df_check.loc[i, '_overlap'] = df_check.loc[i, '_next_start'] >= df_check.loc[i, 'end']

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
    # Build full output paths
    gtf_path = os.path.join(os.path.join(output_dir, "dolphin_exon_gtf", f"{base_name}.gtf"))
    pkl_path = os.path.join(os.path.join(output_dir, "dolphin_exon_gtf", output_dir, f"{base_name}.pkl"))

    # Write to GTF
    gtfpy.writeGTF(gtf_df, gtf_path)

    # Write to Pickle
    gtf_df.to_pickle(pkl_path)

    print(f"GTF file saved to: {gtf_path}")
    print(f"Pickle file saved to: {pkl_path}")

def generate_nonoverlapping_exons(input_gtf_path: str, output_dir: str = "./", batch_size: int = 10000):
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
    
    # Step 1: Load exon entries from GTF and remove duplicates
    df_exon_nodup = prepare_exon_gtf(input_gtf_path, output_dir=output_dir)

    # Step 2: Process and save exons by gene in batches
    save_by_batch(df_exon_nodup, save_num=batch_size, output_dir=output_dir)

    # Step 3: Combine saved batches
    gtf_all = combine_saved_batches(folder=output_dir)

    # Step 4: Check for residual overlaps
    overlap_issues = check_exon_overlap(gtf_all, expected_gene_list=df_exon_nodup["gene_id"].unique().tolist())

    # Step 5: Save final GTF and pickle files
    save_gtf_outputs(gtf_all, output_dir=output_dir)
    
    print(f"[Success] Exon GTF processing pipeline completed.")

    return gtf_all, overlap_issues

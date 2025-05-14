import os
import math
import shutil
import pandas as pd
import anndata
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from functools import partial
from .func_step02_adj_mat_main_part1_main_1 import combine_adj

def _combine_adj_wrapper(args):
    return combine_adj(*args)

def run_adjacency_combination(
    metadata_path: str,
    graph_directory: str,
    adj_meta_file: str,
    out_name: str,
    out_directory="./",
    adj_run_num=50,
    clean_temp: bool = True,
    parallel: bool = True,
):
    """
    Run adjacency matrix combination in batches and merge results into a final AnnData object.

    Parameters
    ----------
    metadata_path : str
        Path to the metadata file with cell barcodes.
    graph_directory : str
        Path to directory containing cell-level _adj.csv files.
    adj_meta_file : str
        Path adjacency metatable dolphin_adj_metadata_table.csv.
    out_name : str
        Output name prefix.
    out_directory : str
        Output folder to save results.
    adj_run_num : int
        Number of cells to combine per batch. Default is 25.
    clean_temp : bool
        Whether to delete temporary intermediate batch files.
    parallel : bool
        If True, run batches in parallel. Default is True.
    """
    print("Start Combining Adjacency Matrix...")
    df_label = pd.read_csv(metadata_path, sep='\t')
    total_sample_size = len(df_label)

    final_out_dir = os.path.join(out_directory, "data")
    os.makedirs(final_out_dir, exist_ok=True)
    temp_out_dir = os.path.join(final_out_dir, "temp")
    os.makedirs(temp_out_dir, exist_ok=True)

    # 1. Prepare all batch arguments
    args_list = []
    for i in range(0, total_sample_size, adj_run_num):
        args_list.append((
            df_label,
            graph_directory,
            adj_meta_file,
            i,
            adj_run_num,
            temp_out_dir,
            out_name
        ))

    # 2. Run combine_adj for each batch
    if parallel:
        print(f"Running in parallel with batch size = {adj_run_num} ...")
        with Pool(processes=cpu_count()) as pool:
            for idx, _ in enumerate(pool.imap_unordered(_combine_adj_wrapper, args_list), start=1):
                print(f"[{idx}/{len(args_list)}] Finished batch")
    else:
        print(f"Running sequentially with batch size = {adj_run_num} ...")
        for idx, args in enumerate(args_list, start=1):
            _ = _combine_adj_wrapper(args)
            print(f"[{idx}/{len(args_list)}] Finished batch")

    # 3. Merge batch .h5ad files into one final file
    print("Merging .h5ad batches...")
    total_batches = math.ceil(total_sample_size / adj_run_num)
    adata_list = [
        anndata.read_h5ad(os.path.join(temp_out_dir, f"Adjacency_{out_name}_{i}.h5ad"))
        for i in range(total_batches)
    ]
    combined_adata = adata_list[0]
    for ad in adata_list[1:]:
        combined_adata = combined_adata.concatenate(ad, index_unique=None, batch_key=None)
    combined_adata.write(os.path.join(final_out_dir, f"Adjacency_{out_name}.h5ad"))

    if clean_temp:
        print("Cleaning up temporary files...")
        shutil.rmtree(temp_out_dir)

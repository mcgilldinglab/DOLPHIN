import os
import math
import shutil
import pandas as pd
import anndata
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from functools import partial
from .func_step02_adj_mat_main_part2_main_3_combine import combine_adj_comp

def _combine_adj_comp_wrapper(args):
    return combine_adj_comp(*args)

def run_adjacency_compress_combination(
    metadata_path: str,
    out_name: str,
    out_directory: str = "./",
    adj_run_num: int = 50,
    clean_temp: bool = True,
    parallel: bool = True,
):
    """
    Combine compressed adjacency matrices in batches and merge into a final AnnData object.

    Parameters
    ----------
    metadata_path : str
        Path to the metadata file with cell barcodes.
    out_name : str
        Output name prefix.
    out_directory : str
        Output folder to save results.
    adj_run_num : int
        Number of cells to combine per batch. Default is 50.
    clean_temp : bool
        Whether to delete temporary intermediate batch files.
    parallel : bool
        If True, run batches in parallel.
        
    Returns
    -------
    None
        Saves the compressed adjacency matrix to the output directory as `AdjacencyComp_<out_name>.h5ad`.

    """
    print("Start Combining Compressed Adjacency Matrix...")

    df_label = pd.read_csv(metadata_path, sep='\t')
    total_sample_size = len(df_label)

    final_out_dir = os.path.join(out_directory, "data")
    temp_out_dir = os.path.join(final_out_dir, "temp")
    os.makedirs(temp_out_dir, exist_ok=True)

    # 1. Prepare batch arguments
    args_list = [
        (df_label, i, adj_run_num, temp_out_dir, out_name)
        for i in range(0, total_sample_size, adj_run_num)
    ]

    # 2. Run in parallel or sequential
    if parallel:
        print(f"Running in parallel with batch size = {adj_run_num} ...")
        with Pool(processes=cpu_count()) as pool:
            for idx, _ in enumerate(tqdm(pool.imap_unordered(_combine_adj_comp_wrapper, args_list), total=len(args_list))):
                pass
    else:
        print(f"Running sequentially with batch size = {adj_run_num} ...")
        for idx, args in enumerate(tqdm(args_list), start=1):
            _combine_adj_comp_wrapper(args)

    # 3. Merge all temporary .h5ad files
    print("Merging .h5ad batches...")
    total_batches = math.ceil(total_sample_size / adj_run_num)
    adata_list = [
        anndata.read_h5ad(os.path.join(temp_out_dir, f"AdjacencyComp_{out_name}_{i}.h5ad"))
        for i in range(total_batches)
    ]
    combined_adata = adata_list[0]
    for ad in adata_list[1:]:
        combined_adata = combined_adata.concatenate(ad, index_unique=None, batch_key=None)

    final_output_path = os.path.join(final_out_dir, f"AdjacencyComp_{out_name}.h5ad")
    combined_adata.write(final_output_path)

    # 4. Clean up temporary files
    if clean_temp:
        print("Cleaning up temporary files...")
        shutil.rmtree(temp_out_dir)

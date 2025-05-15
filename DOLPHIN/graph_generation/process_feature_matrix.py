from tqdm import tqdm
import math
import anndata
import os
import pandas as pd
import torch
from .func_step01_fea_mat_main_part1 import combine_fea
from .func_step01_fea_mat_main_part2 import fea_comp

def run_feature_combination(
    metadata_path: str,
    graph_directory: str,
    gene_annotation,
    gtf_pkl_path: str,
    out_name: str,
    out_directory="./",
    fea_run_num=100,
    clean_temp: bool = True
):
    """
    Run feature matrix combination in batches and merge the results into a final AnnData object.

    This function reads cell metadata and processes each cell's feature vectors in batches.
    It combines features and saves intermediate `.h5ad` files for each batch,
    and finally concatenates them into one unified `.h5ad` file. This is useful for large-scale
    datasets where memory-efficient batch processing is necessary.

    Parameters
    ----------
    fea_run_num : int
        Number of cells to process per batch, default is 100.
    metadata_path : str
        Path to the metadata file (e.g., a csv file with cell information).
    graph_directory : str
        Path to the directory containing graph input files.
    gene_annotation : Any
        Gene annotation data (can be a list, dict, or DataFrame depending on context).
    gtf_pkl_path : str
        Path to the GTF pickle file.
    out_directory : str
        Output directory to save the combined feature matrix, default save to ./data/ folder.
    out_name : str
        Output filename for the feature matrix CSV.
    clean_temp : bool
        Whether to delete the temporary folder after processing. Default is True.

    Returns
    -------
    None
        Saves the following output files:
        - Batch-wise `.h5ad` files for each group of samples.
        - A final merged `.h5ad` file: `Feature_<out_name>.h5ad`.

    """
    
    df_label = pd.read_csv(metadata_path, sep='\t')
    total_sample_size = len(df_label)
    
    # Set output and temp directories
    final_out_dir = os.path.join(out_directory, "data")
    os.makedirs(final_out_dir, exist_ok=True)
    temp_out_dir = os.path.join(final_out_dir, "temp")
    os.makedirs(temp_out_dir, exist_ok=True)

    print("Start Combining Feature Matrix...")
    with tqdm(total=total_sample_size, desc="Combining Features") as pbar_fea:
        for i in range(0, total_sample_size, fea_run_num):
            pbar_fea = combine_fea(
                pbar_fea,
                df_label,
                graph_directory,
                gene_annotation,
                gtf_pkl_path,
                start_idx=i,
                sample_num=fea_run_num,
                output_path=temp_out_dir,
                output_name=out_name
            )
    
    #### combine all feature .h5ad files
    total_number_fea_anndata = math.ceil(total_sample_size/fea_run_num)
    for _idx, _fea_idx in enumerate(range(0, total_number_fea_anndata)):
        _temp_ad = anndata.read_h5ad(os.path.join(temp_out_dir, "Feature_"+out_name+"_"+str(_fea_idx)+".h5ad"))
        if _idx ==0:
            combine_anndata_fea = _temp_ad
        else:
            combine_anndata_fea = combine_anndata_fea.concatenate(_temp_ad, index_unique = None, batch_key = None)
    combine_anndata_fea.write(os.path.join(final_out_dir, "Feature_"+out_name+".h5ad"))

    fea_comp(final_out_dir, out_name)
    
    if clean_temp:
        import shutil
        print("Cleaning up temporary files...")
        shutil.rmtree(temp_out_dir)

from .func_step03_GNN_main import get_graph_input
import torch
from tqdm import tqdm
import math
import anndata
import os
import pandas as pd

def run_model_input(
    metadata_path: str,
    out_name: str,
    out_directory: str = "./",
    gnn_run_num: int = 100,
    celltypename: str = None
    ):
    
    """
    Combines feature matrix and adjacency matrix and generates input for the DOLPHIN model.
    
    Parameters
    ----------
    metadata_path : str
        Path to the metadata file (e.g., a csv file with cell information).
    out_name : str
        Output filename for the feature matrix CSV.
    out_directory : str
        Output directory to save the combined feature matrix, default save to ./data/ folder.
    gnn_run_num : int
        Number of samples per GNN batch.
    celltypename : str, optional
        Column name in metadata indicating cell types. Default is None.

    Returns
    -------
    None
        Saves the final torch tensor file as `model_<out_name>.pt` in the output directory.

        This file contains a list of PyTorch Geometric `Data` objects, one per cell. Each object includes:

        - x : Feature matrix of the cell (normalized exon counts, shaped `[num_features, 1]`)
        - edge_index : Graph connectivity (exon-exon connection indices)
        - edge_attr : Edge weights for the exon graph
        - y : Label for the cell (optional; set to numerical index if `celltypename` is not provided)
        - x_fea : Original feature vector for the cell
        - x_adj : Raw adjacency matrix for the cell
        - sample_name : The ID of the cell

    """
    
    df_label = pd.read_csv(metadata_path, sep='\t')
    total_sample_size = len(df_label)
    
    mapper = None
    if celltypename and celltypename in df_label.columns:
        unique_celltypes = sorted(df_label[celltypename].dropna().unique())
        mapper = {celltype: idx for idx, celltype in enumerate(unique_celltypes)}
    
    final_out_dir = os.path.join(out_directory, "data")
    os.makedirs(final_out_dir, exist_ok=True)
    temp_out_dir = os.path.join(final_out_dir, "temp")
    os.makedirs(temp_out_dir, exist_ok=True)
        
    print("Start Construct Data Input for model input")
    with tqdm(total=total_sample_size) as pbar_gnn:
        for i in range(0, total_sample_size):
            if i%gnn_run_num==0:
                pbar_gnn= get_graph_input(pbar_gnn, 
                                          i, 
                                          gnn_run_num, 
                                          temp_out_dir, 
                                          final_out_dir,
                                          out_name, 
                                          celltypename=celltypename, 
                                          mapper=mapper)

    ##### combine all  geometric .pt files
    total_number_gnn_anndata = math.ceil(total_sample_size/gnn_run_num)
    for _idx, _gnn_idx in enumerate(range(0, total_number_gnn_anndata)):
        _temp_gnn = torch.load(os.path.join(temp_out_dir, "model_"+out_name+"_"+str(_gnn_idx)+".pt"))
        if _idx ==0:
            combine_gnn = _temp_gnn
        else:
            combine_gnn += _temp_gnn
    torch.save(combine_gnn, os.path.join(final_out_dir, "model_"+out_name+".pt"))

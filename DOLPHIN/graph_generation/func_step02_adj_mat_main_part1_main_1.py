import pandas as pd
import numpy as np
import os
import anndata
from scipy.sparse import csr_matrix
import gc
from ._anndata_compat import enable_nullable_string_writes

def build_adjacency_var_annotation(gtf_jun_pkl_path):
    df_jun_gtf = pd.read_csv(gtf_jun_pkl_path)
    var = pd.DataFrame(index=df_jun_gtf["Gene_Junc_name"].astype(str).values)
    var["gene_id"] = df_jun_gtf["Geneid"].astype(str).values
    var["gene_name"] = df_jun_gtf["GeneName"].astype(str).values
    return var


def combine_adj(
    pd_gt,
    graph_path,
    gtf_jun_pkl_path,
    start_idx,
    sample_num,
    output_path,
    output_name,
    adj_var=None,
):
    """
    Combine a batch of adjacency matrices into an AnnData object.
    """
    
    print(f"Processing batch: {start_idx} → {min(start_idx+sample_num, len(pd_gt))}")

    cell_list = list(pd_gt["CB"])
    end_idx = min(start_idx + sample_num, len(cell_list))
    batch_cells = cell_list[start_idx:end_idx]
    cell_adj_list = []

    for _cb in batch_cells:
        cell_adj_list.append(np.loadtxt(os.path.join(graph_path, _cb + "_adj.csv")))

    X = np.vstack(cell_adj_list).astype(np.float32, copy=False)
    batch_meta = pd_gt.set_index("CB").loc[batch_cells].reset_index()

    # Create AnnData
    obs_names = batch_meta["CB"].astype(str).values
    obs = pd.DataFrame(index=obs_names)
    for _col_name in pd_gt.columns:
        obs[_col_name] = batch_meta[_col_name].values
    if adj_var is None:
        adj_var = build_adjacency_var_annotation(gtf_jun_pkl_path)

    adata = anndata.AnnData(X=csr_matrix(X), obs=obs, var=adj_var.copy())

    # Save
    out_path = os.path.join(output_path, f"Adjacency_{output_name}_{int(start_idx/sample_num)}.h5ad")
    enable_nullable_string_writes()
    adata.write(out_path)
    

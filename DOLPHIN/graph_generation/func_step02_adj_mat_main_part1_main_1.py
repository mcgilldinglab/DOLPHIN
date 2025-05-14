import pandas as pd
import numpy as np
import os
import anndata
from scipy.sparse import csr_matrix
import gc

def combine_adj(pd_gt, graph_path, gtf_jun_pkl_path, start_idx, sample_num, output_path, output_name):
    """
    Combine a batch of adjacency matrices into an AnnData object.
    """
    
    print(f"Processing batch: {start_idx} → {min(start_idx+sample_num, len(pd_gt))}")

    cell_list = list(pd_gt["CB"])
    end_idx = min(start_idx + sample_num, len(cell_list))
    cell_adj_list = []

    for _cb in cell_list[start_idx:end_idx]:
        # Load adjacency vector
        _temp_adj = pd.read_csv(os.path.join(graph_path, _cb + "_adj.csv"), header=None).values.flatten()
        # Get label row
        _temp_label = pd_gt[pd_gt["CB"] == _cb].to_numpy().flatten()
        # Concatenate
        _temp_all = np.concatenate((_temp_label, _temp_adj), axis=0)
        cell_adj_list.append(_temp_all)

    # Stack all rows
    cell_adj = np.vstack(cell_adj_list)
    df_adj = pd.DataFrame(cell_adj)

    # Load gene-junction annotation
    df_jun_gtf = pd.read_csv(gtf_jun_pkl_path)

    # Create AnnData
    obs = df_adj.iloc[:, :pd_gt.shape[1]]
    obs.columns = pd_gt.columns  # assign correct column names
    obs.index = obs["CB"].astype(str)  # make sure index is string

    var = pd.DataFrame(index=df_jun_gtf["Gene_Junc_name"].astype(str).values)
    var["gene_id"] = df_jun_gtf["Geneid"].astype(str).values
    var["gene_name"] = df_jun_gtf["GeneName"].astype(str).values

    X = df_adj.iloc[:, pd_gt.shape[1]:].astype(np.float32).values
    adata = anndata.AnnData(X=csr_matrix(X), obs=obs, var=var)

    # Save
    out_path = os.path.join(output_path, f"Adjacency_{output_name}_{int(start_idx/sample_num)}.h5ad")
    adata.write(out_path)
    
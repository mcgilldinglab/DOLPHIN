import pandas as pd
import numpy as np
import os
import anndata
from scipy.sparse import csr_matrix
import gc

def combine_adj(pbar, pd_gt, graph_path, gtf_jun_pkl_path, start_idx, sample_num, output_path, output_name):
    '''
    combine all adjacency matrix together
    '''
    cell_adj = np.array([[]])
    # _df_temp = pd_gt
    cell_list = list(pd_gt["CB"])

    if start_idx+sample_num > len(cell_list):
        end_idx = len(cell_list)
    else:
        end_idx = start_idx+sample_num 

    for i, _cb in enumerate(cell_list[start_idx:end_idx]):
        _temp_adj = np.array([np.loadtxt(os.path.join(graph_path, _cb + "_adj.csv"))])
        _temp_lable = pd_gt[pd_gt["CB"] == _cb].to_numpy()
        _temp_all = np.concatenate((_temp_lable, _temp_adj), axis=1)
        if i == 0:
            cell_adj = _temp_all
        else: 
            cell_adj = np.concatenate((cell_adj, _temp_all), axis=0)
        pbar.update(1)

    df_adj = pd.DataFrame(cell_adj)

    df_jun_gtf = pd.read_pickle(gtf_jun_pkl_path)

    ## adata for adjacency matrix
    obs_names = df_adj[0].values
    obs = pd.DataFrame(index=obs_names)
    for _i, _col_name in enumerate(pd_gt.columns):
        obs[_col_name] = df_adj[_i].values

    ## dataframe for annotating the variables = geneid
    var_names = df_jun_gtf["Gene_Junc_name"].values
    var = pd.DataFrame(index=var_names)
    var["gene_id"] = df_jun_gtf["Geneid"].values
    var["gene_name"] = df_jun_gtf["GeneName"].values

    # ##the data matrix 
    X = df_adj.iloc[:,pd_gt.shape[1]::].values
    adata = anndata.AnnData(X=X, obs=obs, var=var, dtype=np.float32)
    adata.X = csr_matrix(adata.X)

    adata.write(os.path.join(output_path, "Adjacency_"+output_name+"_"+str(int(start_idx/sample_num))+".h5ad"))
    
    del adata
    del cell_adj
    gc.collect()
    
    return pbar
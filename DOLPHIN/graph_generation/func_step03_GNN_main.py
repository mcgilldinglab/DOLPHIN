from torch_geometric.data import Data
import anndata
import numpy as np
import torch
import pandas as pd
import os
import gc
import scipy.sparse as sp

def get_graph_input(pbar, start_idx, sample_num, temp_folder, final_output_path, output_name, celltypename=None, mapper=None):
    adata_adj_raw = anndata.read_h5ad(os.path.join(final_output_path, "AdjacencyCompReHvg_" + output_name+".h5ad"))
    adata_adj = anndata.read_h5ad(os.path.join(final_output_path, "AdjacencyCompReHvgEdge_" + output_name+".h5ad"))
    adata_fea = anndata.read_h5ad(os.path.join(final_output_path, "FeatureCompHvg_" + output_name+".h5ad"))

    adata_fea.X = adata_fea.X.toarray()
    adata_adj_raw.X = adata_adj_raw.X.toarray()
    adata_adj.X = adata_adj.X.toarray()

    gene_list = list(sorted(set(adata_fea.var["gene_id"])))
    sample_list = list(adata_fea.obs_names)

    #get the start and end dimension for adj
    df_edge_idx = pd.DataFrame(adata_fea.var).groupby("gene_id").count().reset_index().rename(columns={"gene_name":"count"})
    df_edge_idx["count"] = df_edge_idx["count"].shift(1).fillna(0)
    df_edge_idx["_idx"] = df_edge_idx["count"].cumsum()
    #mapped sample type
    if celltypename is not None and mapper is not None and celltypename in adata_fea.obs.columns:
        mapped_cell_type = pd.DataFrame(adata_fea.obs)[celltypename].map(mapper)
    else:
        mapped_cell_type = pd.Series(np.zeros(len(adata_fea)), index=adata_fea.obs.index, dtype=int)
        
    if (start_idx+sample_num) > len(sample_list):
        end_idx = len(sample_list)
    else:
        end_idx = start_idx+sample_num

    cell_data_all = []
    # for sample_idx in range(0, len(sample_list)):
    for sample_idx in range(start_idx, end_idx):
        # print(sample_idx)
        #### feature
        _fea = adata_fea[sample_list[sample_idx],:].X
        cell_x = torch.tensor(np.reshape(_fea, (-1,1)), dtype=torch.float32)

        _adj = adata_adj_raw[sample_list[sample_idx],:].X

        #### edge index and weight
        cell_edge = np.array([[],[]])
        cell_edge_weight = np.array([]).reshape(-1,1)
        for i in range(0,len(gene_list)):
            _orig_adj = adata_adj[sample_list[sample_idx],adata_adj.var["gene_id"] == gene_list[i]].X
            _size = int(_orig_adj.shape[1] ** 0.5) 
            _adj_2d = _orig_adj.reshape(_size,_size)
            _edge = np.nonzero(_adj_2d)
            _val = _adj_2d[_edge]
            cell_edge_weight = np.concatenate((cell_edge_weight, _val.reshape(-1, 1)), axis=0)
            np_edge = np.stack(_edge, axis=0) + df_edge_idx[df_edge_idx["gene_id"] == gene_list[i]]["_idx"][i]
            cell_edge = np.concatenate((cell_edge, np_edge), axis=1)

        #### create graph data for one cell
        _data = Data(x = cell_x, edge_index = torch.tensor(cell_edge, dtype=torch.long) , edge_attr = torch.tensor(cell_edge_weight, dtype=torch.float), y = torch.tensor(mapped_cell_type[sample_idx]), x_fea = torch.tensor(_fea, dtype=torch.float32), x_adj = torch.tensor(_adj, dtype=torch.float32), sample_name = sample_list[sample_idx])
        cell_data_all.append(_data)
        pbar.update(1)
        
    torch.save(cell_data_all, os.path.join(temp_folder, "model_" + output_name + "_" + str(int(start_idx/sample_num))+".pt"))
    
    del cell_data_all
    gc.collect()
    
    return pbar
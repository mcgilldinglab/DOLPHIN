from torch_geometric.data import Data
import anndata
import numpy as np
import torch
import pandas as pd
import os
import gc

def get_graph_input(pbar, start_idx, sample_num, output_path, output_name, celltypename=None, mapper=None):
    """
    Converts per-cell exon feature and adjacency data into PyTorch Geometric Data objects.

    This function processes slices of precomputed feature and adjacency matrices stored as AnnData `.h5ad` files,
    and transforms them into `torch_geometric.data.Data` objects. It assembles per-cell graphs with node features,
    edge indices, edge weights, and optional cell type labels. The output is saved as a `.pt` file.

    Parameters
    ----------
    pbar : tqdm.tqdm
        A progress bar object to track batch graph generation progress.
    start_idx : int
        Starting index of the cells to be processed in the batch.
    sample_num : int
        Number of cells to process in this batch.
    output_path : str
        Path where the required `.h5ad` inputs are located and `.pt` file will be saved.
    output_name : str
        Base name used to locate input files and name the output.
    celltypename : str, optional
        The column name in `adata_fea.obs` that holds the cell type label for each cell.
        If not provided, all cells will be assigned label 0.
    mapper : dict, optional
        A dictionary mapping cell type names to integer class labels.
        Required if `celltypename` is specified.

    Returns
    -------
    pbar : tqdm.tqdm
        The updated progress bar object after graph creation.

    Notes
    -----
    - This function reads three files: feature matrix (`FeatureCompHvg_*.h5ad`), 
      processed adjacency (`AdjacencyCompReHvg_*.h5ad`), and raw adjacency (`*_raw.h5ad`).
    - The adjacency matrices are reshaped into square matrices per gene, and edge lists are built from non-zero entries.
    - Graph data includes:
        - `x`: per-exon feature vector (as node feature),
        - `edge_index`: graph connectivity,
        - `edge_attr`: edge weights from adjacency values,
        - `y`: cell type label (optional)
        - `x_fea`, `x_adj`: raw flat vectors used as reconstruction targets for the decoders.

    - Output is saved as `geometric_<output_name>.pt`.

    Example
    -------
    >>> from tqdm import tqdm
    >>> pbar = tqdm(total=50)
    >>> get_graph_input(pbar, 0, 50, './output/', 'batch1', celltypename='celltype', mapper={'tumor': 0, 'immune': 1})
    """
    
    #convert adjacency matrix to row and columns
    adata_adj_raw = anndata.read_h5ad(os.path.join(output_path, "AdjacencyCompReHvg_" + output_name+"_raw.h5ad"))
    adata_adj = anndata.read_h5ad(os.path.join(output_path, "AdjacencyCompReHvg_" + output_name+".h5ad"))
    adata_fea = anndata.read_h5ad(os.path.join(output_path, "FeatureCompHvg_" + output_name+".h5ad"))

    adata_fea.X = adata_fea.X.toarray()

    gene_list = list(sorted(set(adata_fea.var["gene_id"])))
    sample_list = list(adata_fea.obs_names)
    #get the start and end dimension for adj
    df_edge_idx = pd.DataFrame(adata_fea.var).groupby("gene_id").count().reset_index().rename(columns={"gene_name":"count"})
    df_edge_idx["count"] = df_edge_idx["count"].shift(1)
    df_edge_idx["_idx"] = df_edge_idx["count"].cumsum()
    df_edge_idx = df_edge_idx.fillna(0)

    #mapped sample type
    # mapper = {'tumor': 0, 'endothelial': 1, 'fibroblast': 2, 'follicular': 3, 'myeloid': 4, 'plasma': 5, 'T & NK': 6}
    '''HERE need modification, 这里需要增加一个celltype的名字,如果有的话assing mapper,没有的话就都设成0'''
    if celltypename is not None:
        mapped_cell_type = pd.DataFrame(adata_fea.obs)[celltypename].map(mapper)
    else:
        mapped_cell_type = pd.Series(0, index=adata_fea.obs.index)
        
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
        
    torch.save(cell_data_all, os.path.join(output_path, "geometric_" + output_name + "_" + str(int(start_idx/sample_num))+".pt"))
    
    del cell_data_all
    gc.collect()
    
    return pbar
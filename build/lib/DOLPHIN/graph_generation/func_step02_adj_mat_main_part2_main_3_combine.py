import anndata
from scipy.sparse import csr_matrix
import os
import gc

def combine_adj_comp(pd_gt, start_idx, sample_num, output_path, output_name):
    sample_list = list(pd_gt["CB"])

    #combine all the single h5ad file together, process in compute canada
    _temp_a1 = anndata.read_h5ad(os.path.join(output_path, "adj_comp_matrix", sample_list[start_idx]+".h5ad"))
    adata_adj_up = _temp_a1.copy()

    if (start_idx+sample_num) > len(sample_list):
        end_idx = len(sample_list)
    else:
        end_idx = start_idx+sample_num
    
    for i in range(start_idx+1, end_idx):
        _temp_adj = anndata.read_h5ad(os.path.join(output_path, "adj_comp_matrix", sample_list[i]+".h5ad"))
        adata_adj_up = adata_adj_up.concatenate(_temp_adj, index_unique = None, batch_key = None)

    adata_adj_up.X = csr_matrix(adata_adj_up.X)
    adata_adj_up.write(os.path.join(output_path, "AdjacencyComp_"+output_name+"_"+str(int(start_idx/sample_num))+".h5ad"))

    del adata_adj_up
    gc.collect()
    
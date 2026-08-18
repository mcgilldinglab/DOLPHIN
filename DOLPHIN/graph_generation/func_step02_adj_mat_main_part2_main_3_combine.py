import anndata
from scipy.sparse import csr_matrix
import os
import gc
from ._anndata_compat import enable_nullable_string_writes

def combine_adj_comp(pd_gt, start_idx, sample_num, output_path, output_name):
    sample_list = list(pd_gt["CB"])

    end_idx = min(start_idx + sample_num, len(sample_list))
    batch_sample_list = sample_list[start_idx:end_idx]
    adata_list = [
        anndata.read_h5ad(os.path.join(output_path, "adj_comp_matrix", sample_id + ".h5ad"))
        for sample_id in batch_sample_list
    ]
    adata_adj_up = anndata.concat(
        adata_list,
        index_unique=None,
        merge="same",
    )

    adata_adj_up.X = csr_matrix(adata_adj_up.X)
    enable_nullable_string_writes()
    adata_adj_up.write(
        os.path.join(
            output_path,
            "AdjacencyComp_" + output_name + "_" + str(int(start_idx / sample_num)) + ".h5ad",
        )
    )

    del adata_adj_up
    del adata_list
    gc.collect()
    

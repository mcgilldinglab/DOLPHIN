import os

import anndata
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse

from ._anndata_compat import enable_nullable_string_writes


##################################################################################################################
## Convert Full Feature Matrix to Compact Feature Matrix
## Only the exon whose expression level is zero across all the samples needs to be removed.
##################################################################################################################
def fea_comp(output_path, output_name):
    adata_fea_orig = anndata.read_h5ad(os.path.join(output_path, "Feature_" + output_name + ".h5ad"))

    matrix = adata_fea_orig.X
    if issparse(matrix):
        keep_mask = np.asarray(matrix.getnnz(axis=0)).ravel() > 0
    else:
        keep_mask = np.asarray(matrix).sum(axis=0) != 0

    adata_fea_comp = adata_fea_orig[:, keep_mask].copy()

    df_var = pd.DataFrame(adata_fea_comp.var).copy().reset_index()
    df_var["orig_idx_order"] = df_var["index"].apply(lambda x: int(str(x).split("-")[-1]))
    df_var["gene_order"] = df_var["gene_id"].apply(lambda x: int(str(x)[4:]))
    df_var = df_var.sort_values(by=["gene_order", "orig_idx_order"]).reset_index(drop=True)
    df_var["new_index"] = df_var.groupby(["gene_id"]).cumcount() + 1
    df_var["var_new_index"] = df_var["gene_id"].astype(str) + "-" + df_var["new_index"].astype(str)

    reordered_matrix = adata_fea_comp.X[:, df_var.index.to_numpy()]
    if not issparse(reordered_matrix):
        reordered_matrix = csr_matrix(reordered_matrix)
    else:
        reordered_matrix = reordered_matrix.tocsr()

    obs = pd.DataFrame(adata_fea_comp.obs)
    var = pd.DataFrame(index=df_var["var_new_index"].values)
    var["gene_id"] = df_var["gene_id"].values
    var["gene_name"] = df_var["gene_name"].values

    adata = anndata.AnnData(
        X=reordered_matrix.astype(np.float32),
        obs=obs,
        var=var,
        dtype=np.float32,
    )

    enable_nullable_string_writes()
    adata.write(os.path.join(output_path, "FeatureComp_" + output_name + ".h5ad"))

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import math
import os
from scipy import sparse

def run_feature_hvg(
    out_name: str,
    out_directory="./"):
    """
    Filters the feature matrix to retain only highly variable genes (HVGs), normalizes the data,
    and prepares it as input for the DOLPHIN model.   
     
    Parameters
    ----------
    out_name : str
        Output filename for the feature matrix CSV.
    out_directory : str
        Output directory to save the combined feature matrix, default save to ./data/ folder.
    Returns
    -------
    None
        Saves the final feature matrix as `FeatureCompHvg_<out_name>.h5ad` in the specified output directory.
    """
    final_out_dir = os.path.join(out_directory, "data")
    os.makedirs(final_out_dir, exist_ok=True)
    
    hvg_path = os.path.join(final_out_dir, "ExonGene_hvg_"+out_name+".h5ad")
    feature_anndata = os.path.join(final_out_dir, "FeatureComp_"+out_name+".h5ad")

    ## load feature matrix
    adata = anndata.read_h5ad(feature_anndata)

    hvg_adata = anndata.read_h5ad(hvg_path)
    cell_keep = list(hvg_adata.obs.index)
    hvg_list = set(hvg_adata.var["Geneid"])

    adata = adata[adata.obs.index.isin(cell_keep), :]

    #normalize the feature matrix and round to nearest integer 
    sc.pp.normalize_total(adata)

    #only keep highly varaible genes
    adata_hvg = adata[:,adata.var["gene_id"].isin(hvg_list)].copy()
    
    # if sparse.issparse(adata_hvg.X):
    #     X_rounded = sparse.csr_matrix(np.round(adata_hvg.X.toarray()).astype(np.float32))  # .A converts to dense
    # else:
    #     X_rounded = np.round(adata_hvg.X).astype(np.float32)
    
    # X_rounded_sparse = sparse.csr_matrix(np.round(X_rounded).astype(np.float32))
    # adata_final = anndata.AnnData(X=X_rounded_sparse, obs=adata_hvg.obs.copy(), var=adata_hvg.var.copy())

    if sparse.issparse(adata_hvg.X):
        X_dense = adata_hvg.X.toarray()
    else:
        X_dense = adata_hvg.X

    X_rounded_sparse = sparse.csr_matrix(np.round(X_dense).astype(np.float32))

    adata_final = anndata.AnnData(X=X_rounded_sparse, obs=adata_hvg.obs.copy(), var=adata_hvg.var.copy())

    print("Keep "+ str(len(set(adata_final.var["gene_id"]))) + " genes")
    print("The Final Feature Matrix Size is " + str(adata_final.shape[0]) + " Cells and " +str(adata_final.shape[1])+ " exons")

    #remove unneccary var names
    adata_final.write(os.path.join(final_out_dir, "FeatureCompHvg_"+out_name+".h5ad"))
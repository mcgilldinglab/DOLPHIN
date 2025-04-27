import anndata
import os
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

##################################################################################################################
## Convert Full Feature Matrix to Compact Feature Matrix
## Only the exon whose expression level is zero across all the samples needs to be removed.
##################################################################################################################
def fea_comp(output_path, output_name):
    #check the exon is zero across all the samples
    adata_fea_orig = anndata.read_h5ad(os.path.join(output_path, "Feature_"+output_name+".h5ad"))

    #### remove exon is empty across all the cells
    sc.pp.filter_genes(adata_fea_orig, min_cells=1)
    #check if any exon has zero values across all samples == False mean no zeros: https://github.com/scverse/scanpy/issues/1083
    np.any(adata_fea_orig.X.sum(axis=0) == 0)
    adata_fea_orig

    #### sort var based on gene id
    df_fea_orig_var = pd.DataFrame(adata_fea_orig.var)
    df_fea_comp = adata_fea_orig.to_df() #cleaned feature matrix
    df_fea_comp_add_var = pd.merge(df_fea_comp.T, df_fea_orig_var, how="left",left_index=True, right_index=True)
    df_fea_comp_add_var = df_fea_comp_add_var.reset_index()
    df_fea_comp_add_var["orig_idx_order"] = df_fea_comp_add_var["index"].apply(lambda x: int(x.split('-')[-1]))
    df_fea_comp_add_var["gene_order"] = df_fea_comp_add_var["gene_id"].apply(lambda x: int(x[4:]))
    df_fea_comp_reorder = df_fea_comp_add_var.sort_values(by=["gene_order", "orig_idx_order"])
    df_fea_comp_reorder["new_index"] = df_fea_comp_reorder.groupby(["gene_id"]).cumcount()+1
    df_fea_comp_reorder["var_new_index"] = df_fea_comp_reorder['gene_id'].astype(str) +"-"+ df_fea_comp_reorder["new_index"].astype(str)
    df_fea_comp_reorder = df_fea_comp_reorder.reset_index(drop=True)
    df_fea_comp_reorder

    ## adata for feature matrix
    obs = pd.DataFrame(adata_fea_orig.obs)

    ## dataframe for annotating the variables = geneid
    var_names = df_fea_comp_reorder["var_new_index"].values
    var = pd.DataFrame(index=var_names)
    var["gene_id"] = df_fea_comp_reorder["gene_id"].values
    var["gene_name"] = df_fea_comp_reorder["gene_name"].values

    # # ##the data matrix 
    X = df_fea_comp_reorder.drop(columns=["index"]).T.iloc[:-7,:].values
    adata = anndata.AnnData(X=X, obs=obs, var=var, dtype=np.float32)
    adata.X = csr_matrix(adata.X)
    
    adata.write(os.path.join(output_path, "FeatureComp_"+output_name+".h5ad"))

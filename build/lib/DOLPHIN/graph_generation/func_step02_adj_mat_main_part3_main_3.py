import pandas as pd
import anndata
import pickle
import os
from scipy.sparse import csr_matrix

def adj_comp_re_part3(output_path, output_name):
    df_adj_comp_re = pd.read_pickle(os.path.join(output_path, "AdjacencyCompRe_part2_"+output_name+".pkl"))

    df_adj_comp_re["new_index"] = df_adj_comp_re.groupby(["gene_id"]).cumcount()+1
    df_adj_comp_re["var_new_index"] = df_adj_comp_re['gene_id'].astype(str) +"-"+ df_adj_comp_re["new_index"].astype(str)
    df_adj_comp_re = df_adj_comp_re.reset_index(drop=True)

    ## == step6: adata for adjacency matrix
    adata_adj_comp = anndata.read_h5ad(os.path.join(output_path, "AdjacencyComp_"+output_name+".h5ad"))

    obs = pd.DataFrame(adata_adj_comp.obs)
    # obs = obs.drop(columns=['batch'])
    ## dataframe for annotating the variables = geneid
    var_names = df_adj_comp_re["var_new_index"].values
    var = pd.DataFrame(index=var_names)
    var["gene_id"] = df_adj_comp_re["gene_id"].values
    # var["gene_name"] = df_gtf_an["GeneName"].values
    var["gene_name"] = df_adj_comp_re["gene_name"].values

    # # ##the data matrix 
    X = df_adj_comp_re.drop(columns=["gene_id", "gene_name","_ck_empty","_idx","adj_vec_size","_flag", "gene_order", "exon_order","exon_name","new_index","var_new_index"]).T.iloc[:,:].values
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    adata.X = csr_matrix(adata.X)
    
    adata.write(os.path.join(output_path, "AdjacencyCompRe_" + output_name+".h5ad"))

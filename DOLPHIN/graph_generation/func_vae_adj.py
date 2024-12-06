import anndata
import os
import scanpy as sc
import pandas as pd
import numpy as np
from tqdm import tqdm


def get_adj_hvg(output_path, output_name):

    hvg_path = os.path.join(output_path, "ExonGene_hvg_"+output_name+".h5ad")

    ## load h5ad file
    adj_anndata = anndata.read_h5ad(os.path.join(output_path, "AdjacencyCompRe_"+output_name+".h5ad"))

    hvg_adata = anndata.read_h5ad(hvg_path)
    hvg_list = set(hvg_adata.var["Geneid"])
    cell_keep = list(hvg_adata.obs.index)

    #only keep highly varaible genes = 5000 genes, where highly_variable == True
    adj_anndata = adj_anndata[:,adj_anndata.var["gene_id"].isin(hvg_list)]

    print("Keep "+ str(len(set(adj_anndata.var["gene_id"]))) + " genes")
    print("The Final Adjacency Matrix Size is " + str(adj_anndata.shape[0]) + " Cells with dimension of " +str(adj_anndata.shape[1]))

    ## Normalization
    # #convert h5ad to dataframe
    df_adj = adj_anndata.to_df().transpose()
    df_adj["gene_id"] = df_adj.index
    df_adj["gene_id"] = df_adj["gene_id"].apply(lambda x: x.rsplit('-',1)[0])
    # #sum the count per gene
    df_adj_sum = df_adj.groupby("gene_id").sum().reset_index()

    srr = adj_anndata.obs.index.to_list()

    for i in tqdm(range(0,len(srr))):
        # print(i)
        _sub_df_adj_sum = df_adj_sum[[srr[i], "gene_id"]]
        _sub_df_adj_sum = _sub_df_adj_sum.rename(columns={srr[i]: "sum_count"})
        _sub_df_adj = df_adj[[srr[i], "gene_id"]].reset_index()
        _sub_merge = pd.merge(_sub_df_adj, _sub_df_adj_sum, how="left",left_on="gene_id",right_on="gene_id")
        _sub_merge[srr[i]] = np.where(_sub_merge["sum_count"]!= 0, _sub_merge[srr[i]]/_sub_merge["sum_count"], 0)
        _sub_merge.set_index("index", inplace=True)
        _sub_merge.index.name = None
        _sub_merge = _sub_merge[[srr[i]]]
        if i == 0:
            df_adj_norm = _sub_merge
        else:
            df_adj_norm = pd.merge(df_adj_norm, _sub_merge, how="outer", left_index=True, right_index=True)

    #conver to h5ad file
    # # # ##the data matrix 
    X = df_adj_norm.transpose().iloc[:,:].values
    new_adata = anndata.AnnData(X, obs=pd.DataFrame(adj_anndata.obs), var=pd.DataFrame(adj_anndata.var))

    results_file = os.path.join(output_path, "AdjacencyCompReHvg_"+ output_name+".h5ad")

    new_adata.write(results_file)

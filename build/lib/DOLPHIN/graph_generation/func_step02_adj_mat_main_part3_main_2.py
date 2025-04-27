import pandas as pd
import numpy as np
import anndata
# from matplotlib import pyplot as plt
# from collections import Counter
# from itertools import combinations
# from itertools import tee
import pickle
#doing convolution as pre-process method
# from scipy.linalg import fractional_matrix_power
import os
# from scipy.sparse import csr_matrix

def adj_comp_re_part2(output_path, output_name):
    adata_adj_comp = anndata.read_h5ad(os.path.join(output_path, "AdjacencyComp_"+output_name+".h5ad"))
    #convert sparse matrix back to dense matrix
    adata_adj_comp.X = adata_adj_comp.X.toarray()

    ## == step2: get compact adjacency h5ad, and convert to dataframe for groupby analysis
    #get geneid and gene name dictionary
    df_adata_adj_comp_var = pd.DataFrame(adata_adj_comp.var)
    dict_gene_name_id = dict(zip(df_adata_adj_comp_var.gene_name, df_adata_adj_comp_var.gene_id))
    #2A. label the exon who is empty across all the samples
    df_adata_adj_comp = adata_adj_comp.to_df().T
    df_adata_adj_comp["_ck_empty"] = np.amax(df_adata_adj_comp.to_numpy(), axis = 1)
    df_adata_adj_comp["_idx"] = df_adata_adj_comp.index
    df_adata_adj_comp["gene_name"] =  df_adata_adj_comp["_idx"].apply(lambda x: x[:x.rfind("-")])
    df_adata_adj_comp["gene_id"] = df_adata_adj_comp['gene_name'].map(dict_gene_name_id)

    #2B. get the adj vector size per gene
    # _df_adj_size = pd.DataFrame(df_adata_adj_comp.groupby("gene_name")["gene_name"].count()).rename(columns={"gene_name":"adj_vec_size"}).reset_index()
    _dict_adj_size = df_adata_adj_comp.groupby("gene_id")["gene_id"].count().to_dict()
    df_adata_adj_comp["adj_vec_size"] = df_adata_adj_comp['gene_id'].map(_dict_adj_size)

    ### step2C: convert type
    df_adata_adj_comp["gene_id"] = df_adata_adj_comp["gene_id"].astype("object")
    df_adata_adj_comp["gene_order"] = df_adata_adj_comp["gene_id"].apply(lambda x: int(x[4:]))
    df_adata_adj_comp["exon_order"] = df_adata_adj_comp["_idx"].apply(lambda x: int(x.split("-")[-1]))
    df_adata_adj_comp["exon_name"] =  df_adata_adj_comp['gene_id'].astype(str) +"-"+ df_adata_adj_comp["exon_order"].astype(str)

    ##== step3A: find zero edges across all the samples and flag the correspoding exon
    df_add_flag = pd.read_pickle(os.path.join(output_path, "AdjacencyCompRe_part1_"+output_name+".pkl"))
    df_add_flag["exon_name"] =  df_add_flag['gene_id'].astype(str) +"-"+ df_add_flag["exon_idx"].astype(str)

    # #sort based on origal order
    df_add_flag = df_add_flag.set_index('exon_name')
    df_add_flag = df_add_flag.reindex(index=df_adata_adj_comp['exon_name'])
    df_add_flag = df_add_flag.reset_index()

    ##== step3C: Merge flag column with adj_comp dataframe
    df_adata_adj_comp["_flag"] = df_add_flag["_flag"].values

    print("Start")
    #filter since large dataframe, seperate filter
    num_list = [652842, 652778, 651448, 654284, 652881, 652770, 653024, 652559, 653107, 653788]
    cum = 0
    for n_idx in range(len(num_list)):
        print(n_idx)
        if n_idx == 0:
            temp_df = df_adata_adj_comp.head(num_list[n_idx])
        else:
            df_sub = df_adata_adj_comp.iloc[cum:]
            temp_df = df_sub.head(num_list[n_idx])
        cum = cum + num_list[n_idx]
        out_df = temp_df[temp_df["_flag"] == 0]
        out_df.to_pickle(os.path.join(output_path, "AdjacencyCompRe_part2_"+output_name+"_"+str(n_idx)+".pkl"))

    for n_idx in range(len(num_list)):
        df_temp = pd.read_pickle(os.path.join(output_path, "AdjacencyCompRe_part2_"+output_name+"_"+str(n_idx)+".pkl"))
        if n_idx == 0:
            df_adj_comp_re = df_temp
        else:
            df_adj_comp_re = pd.concat([df_adj_comp_re, df_temp], ignore_index=True)

    df_adj_comp_re.to_pickle(os.path.join(output_path, "AdjacencyCompRe_part2_"+output_name+".pkl"))
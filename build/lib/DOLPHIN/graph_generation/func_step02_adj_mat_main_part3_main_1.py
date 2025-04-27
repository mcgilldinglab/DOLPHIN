import pandas as pd
import numpy as np
import anndata
# from collections import Counter
# from itertools import combinations
# from itertools import tee
import pickle
#doing convolution as pre-process method
# from scipy.linalg import fractional_matrix_power
import os
from scipy.sparse import csr_matrix

def adj_comp_re_part1(output_path, output_name):
    adata_adj_comp = anndata.read_h5ad(os.path.join(output_path, "AdjacencyComp_"+output_name+".h5ad"))
    #convert sparse matrix back to dense matrix
    adata_adj_comp.X = adata_adj_comp.X.toarray()

    adata_fea_gtf = anndata.read_h5ad(os.path.join(output_path, "Feature_"+output_name+".h5ad"))
    sample_list = list(adata_fea_gtf.obs_names)

    ###===Step 1: get the exon index from feature matrix
    #1A. Rename the feature gtf matrix, since the post-processing is groupby gene id not gene name
    df_fea_orig_var = pd.DataFrame(adata_fea_gtf.var)
    df_fea_gtf = pd.merge(adata_fea_gtf.to_df().T, df_fea_orig_var, how="left",left_index=True, right_index=True)
    df_fea_gtf["gene_order"] = df_fea_gtf["gene_id"].apply(lambda x: int(x[4:]))
    df_fea_gtf["gtf_index"] = df_fea_gtf.index
    df_fea_gtf["exon_order"] = df_fea_gtf["gtf_index"].apply(lambda x: int(x.split("-")[-1]))
    df_fea_gtf_order = df_fea_gtf.sort_values(by=["gene_order", "exon_order"])
    df_fea_gtf_order = df_fea_gtf_order.reset_index(drop=True)
    df_fea_gtf_order["new_gtf_exon_index"] = df_fea_gtf_order.groupby(["gene_id"]).cumcount()+1
    df_fea_gtf_order["new_gtf_index"] = df_fea_gtf_order['gene_id'].astype(str) +"-"+ df_fea_gtf_order["new_gtf_exon_index"].astype(str)
    df_fea_gtf_new_index = df_fea_gtf_order.set_index("new_gtf_index")
    df_fea_gtf_new_index = df_fea_gtf_new_index[sample_list]
    df_fea_gtf_new_index = df_fea_gtf_new_index.T

    #1B. find the exon who is empty across all the samples
    df_fea_exon_keep = df_fea_gtf_new_index.loc[:, (df_fea_gtf_new_index > 0).any()]
    exon_keep_list = df_fea_exon_keep.columns

    #1C. seperate exon_keep_list to gene_id and it's correspoding idx
    fea_gene_keep_list = []
    fea_ex_id_keep_list = []
    for exon_name in exon_keep_list:
        us_idx = [i for i, ltr in enumerate(exon_name) if ltr == "-"]
        fea_ex_id_keep_list.append(exon_name.split('-')[-1])
        fea_gene_keep_list.append(exon_name[:us_idx[-1]])

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

    ##== step3A: find zero edges across all the samples and flag the correspoding exon
    def flag_zero(x):
        print(x["gene_id"].iloc[0])
        _size = int(x["adj_vec_size"].iloc[0])
        _2d_size = int(x["adj_vec_size"].iloc[0] ** 0.5)

        fea_exon_to_keep = np.array(fea_ex_id_keep_list)[np.where(np.array(fea_gene_keep_list) == x["gene_id"].iloc[0])]

        _orig_list = np.zeros(_size)
        _df_out = pd.DataFrame(np.zeros(_size).reshape(_2d_size,_2d_size))
        for i in range(0, _size):
            _orig_list[i] = x["_ck_empty"].iloc[i]
        _df_orig = pd.DataFrame(_orig_list.reshape(_2d_size,_2d_size))
        for _idx in range(0, _2d_size):
            if ((_df_orig[_idx] == 0).all()) & ((_df_orig == 0).all(axis=1)[_idx]) & (str(_idx+1) not in fea_exon_to_keep):
                #fill all column with 1
                _df_out[_idx] = int(1)
                #fill all row with 1
                _df_out.loc[_idx] = int(1)
        print("Done")
        return(_df_out.values.flatten())

    print("start")

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
        series_add_ck_zero = temp_df.groupby("gene_id").apply(lambda x:flag_zero(x))

        # #set the series as column back in the original dataframe
        df_add_ck_zero = pd.DataFrame(series_add_ck_zero).rename(columns={0:"_flag"})
        #use numpy array instead of concatenate pandas together
        gene_id_list_all = []
        exon_id_list_all = []
        flag_list_all = []
        for gid in df_add_ck_zero.index.values:
            _current_val = df_add_ck_zero["_flag"][gid]
            #create pandas series
            for i in range(0, len(_current_val)):
                gene_id_list_all.append(gid)
                exon_id_list_all.append(i+1)
            for j in _current_val:
                flag_list_all.append(j)
        df_add_flag = pd.DataFrame({"gene_id":gene_id_list_all, "_flag": flag_list_all, "exon_idx": exon_id_list_all})
        df_add_flag.to_pickle(os.path.join(output_path, "AdjacencyCompRe_part1_" + output_name+"_"+str(n_idx)+".pkl"))
    #combine all temp pkl dataframe
    for n_idx in range(len(num_list)):
        df_temp = pd.read_pickle(os.path.join(output_path, "AdjacencyCompRe_part1_"+output_name+"_"+str(n_idx)+".pkl"))
        if n_idx == 0:
            df_add_flag = df_temp
        else:
            df_add_flag = pd.concat([df_add_flag, df_temp], ignore_index=True)

    df_add_flag.to_pickle(os.path.join(output_path, "AdjacencyCompRe_part1_"+output_name+".pkl"))
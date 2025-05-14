import pandas as pd
import numpy as np
import anndata
import os
from scipy.sparse import csr_matrix
import gc

def combine_fea(pbar, pd_gt, graph_path, gene_annotation, gtf_pkl_path, start_idx, sample_num, output_path, output_name):
    #step1:load gene annoataion file, gene id and gene name (one to one)
    df_an = pd.read_csv(gene_annotation, sep = ',')
    
    #step2: annotate gtf to get gene names:gtf file order is the same as feature/adj orders
    gtf = pd.read_pickle(gtf_pkl_path)
    gtf = gtf[["seqname","gene_id","start","end"]]
    gtf = gtf.rename(columns={"seqname":"Chr","gene_id":"Geneid","start":"Start","end":"End"})
    gtf["Start"] = gtf["Start"].astype(int)
    gtf["End"] = gtf["End"].astype(int)
    
    ## step3: combine gtf and annotation files to get gene-exon name
    
    df_gtf_an = pd.merge(gtf, df_an, how='left', left_on=["Geneid"], right_on=["gene_id"])
    df_gtf_an.gene_name.fillna(df_gtf_an.Geneid, inplace=True)
    df_gtf_an["Gene_Exon_name"] = df_gtf_an.groupby(['gene_name']).cumcount()+1
    df_gtf_an["Gene_Exon_name"] = df_gtf_an['gene_name'].astype(str) +"-"+ df_gtf_an["Gene_Exon_name"].astype(str)

    ## Step4: combine all feature matrix together per sample id
    
    cell_fea = np.array([[]])
    _df_temp = pd_gt
    cell_list = list(pd_gt["CB"])

    if start_idx+sample_num > len(cell_list):
        end_idx = len(cell_list)
    else:
        end_idx = start_idx+sample_num 
    for i, _cb in enumerate(cell_list[start_idx:end_idx]):
        _temp_fea = np.array([np.loadtxt(os.path.join(graph_path, _cb + "_fea.csv"))])
        _temp_lable = pd_gt[pd_gt["CB"] == _cb].to_numpy()
        _temp_all = np.concatenate((_temp_lable, _temp_fea), axis=1)
        if i == 0:
            cell_fea = _temp_all
        else: 
            cell_fea = np.concatenate((cell_fea, _temp_all), axis=0)
        pbar.update(1)

    df_fea = pd.DataFrame(cell_fea)

    ## adata for feature matrix
    obs_names = df_fea[0].values
    obs = pd.DataFrame(index=obs_names)
    for _i, _col_name in enumerate(pd_gt.columns):
        obs[_col_name] = df_fea[_i].values

    ## dataframe for annotating the variables = geneid
    var_names = df_gtf_an["Gene_Exon_name"].values
    var = pd.DataFrame(index=var_names)
    var["gene_id"] = df_gtf_an["Geneid"].values
    var["gene_name"] = df_gtf_an["gene_name"].values

    # # ##the data matrix 
    X = df_fea.iloc[:,pd_gt.shape[1]::]
    adata = anndata.AnnData(X=X, obs=obs, var=var, dtype=np.float32)
    adata.X = csr_matrix(adata.X)

    adata.write(os.path.join(output_path, "Feature_"+output_name+"_"+str(int(start_idx/sample_num))+".h5ad"))

    del adata
    del cell_fea
    gc.collect()
    
    return pbar
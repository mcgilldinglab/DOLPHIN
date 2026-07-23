import pandas as pd
import numpy as np
import anndata
import os
from scipy.sparse import csr_matrix
import gc
from ._anndata_compat import enable_nullable_string_writes

def build_feature_var_annotation(gene_annotation, gtf_pkl_path):
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
    df_gtf_an["gene_name"] = df_gtf_an["gene_name"].fillna(df_gtf_an["Geneid"])
    df_gtf_an["Gene_Exon_name"] = df_gtf_an.groupby(['gene_name']).cumcount()+1
    df_gtf_an["Gene_Exon_name"] = df_gtf_an['gene_name'].astype(str) +"-"+ df_gtf_an["Gene_Exon_name"].astype(str)

    var_names = df_gtf_an["Gene_Exon_name"].values
    var = pd.DataFrame(index=var_names)
    var["gene_id"] = df_gtf_an["Geneid"].values
    var["gene_name"] = df_gtf_an["gene_name"].values
    return var


def combine_fea(
    pbar,
    pd_gt,
    graph_path,
    gene_annotation,
    gtf_pkl_path,
    start_idx,
    sample_num,
    output_path,
    output_name,
    feature_var=None,
):
    ## Step4: combine all feature matrix together per sample id
    cell_list = list(pd_gt["CB"])

    if start_idx+sample_num > len(cell_list):
        end_idx = len(cell_list)
    else:
        end_idx = start_idx+sample_num 
    batch_cells = cell_list[start_idx:end_idx]
    feature_rows = []
    for _cb in batch_cells:
        feature_rows.append(np.loadtxt(os.path.join(graph_path, _cb + "_fea.csv")))
        pbar.update(1)

    X = np.vstack(feature_rows).astype(np.float32, copy=False)
    batch_meta = pd_gt.set_index("CB").loc[batch_cells].reset_index()

    ## adata for feature matrix
    obs_names = batch_meta["CB"].astype(str).values
    obs = pd.DataFrame(index=obs_names)
    for _i, _col_name in enumerate(pd_gt.columns):
        obs[_col_name] = batch_meta[_col_name].values
    if feature_var is None:
        feature_var = build_feature_var_annotation(gene_annotation, gtf_pkl_path)

    # # ##the data matrix 
    adata = anndata.AnnData(X=X, obs=obs, var=feature_var.copy(), dtype=np.float32)
    adata.X = csr_matrix(adata.X)

    enable_nullable_string_writes()
    adata.write(os.path.join(output_path, "Feature_"+output_name+"_"+str(int(start_idx/sample_num))+".h5ad"))

    del adata
    del feature_rows
    gc.collect()
    
    return pbar

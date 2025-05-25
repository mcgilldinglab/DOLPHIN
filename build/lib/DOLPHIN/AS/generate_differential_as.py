import os
import anndata
import pandas as pd

def run_differential_as(
    outrigger_psi_data: str,
    out_name: str,
    cluster_name: str,
    out_directory: str = "./",
    n_cell: int = 10,
):  
    """
    This function imputes missing PSI values in preparation for downstream differential alternative splicing analysis.

    This function imputes missing PSI values by calculating the average PSI 
    across all splicing events within each cell cluster. The cluster-wise mean is used to fill in 
    NaN-psi values, ensuring that sparsely detected events still receive a representative value. 
    Only splicing events detected in at least n_cell cells are retained for analysis.

    Parameters
    ----------
    outrigger_psi_data : str
        Path to the input <out_name>_PSI.h5ad file containing PSI values with NaNs.
    cluster_name : str
        Column name in adata.obs of <out_name>_PSI.h5ad that contains cell cluster labels.
    out_directory : str, optional
        Directory where the output .h5ad file will be saved. Default is the current directory ("./").
    n_cell : int, optional
        Minimum number of cells in which a splicing event must be detected to be retained. Default is 10.
    out_name : str
        Prefix for the output file. The result will be saved as "<out_name>_PSI_DAS.h5ad".

    Returns
    -------
    adata : anndata.AnnData
        AnnData object with cluster-mean-imputed PSI values, ready for differential splicing analysis.
        Saved to: `<out_directory>/alternative_splicing/<out_name>_PSI_DAS.h5ad`.
    """

    final_out_dir = os.path.join(out_directory, "alternative_splicing")
    os.makedirs(final_out_dir, exist_ok=True)

    adata = anndata.read_h5ad(outrigger_psi_data)
    df_psi_raw = adata.to_df()

    ### filter event which only keep the ones who is exist in at lease in 10 cells
    #get number of cells per each event
    df_psi_raw_t = df_psi_raw.T
    print(f"Total number of splicing events before filtering: {df_psi_raw_t.shape[0]}")
    df_psi_raw_t["count"] = df_psi_raw_t.ge(0.0).sum(axis=1)
    df_psi_raw_t_filter = df_psi_raw_t[df_psi_raw_t["count"] >=n_cell].copy()
    df_psi_raw_t_filter.drop(columns=["count"], inplace=True)
    df_psi_raw_filter = df_psi_raw_t_filter.T
    print(f"Number of splicing events after filtering (>= {n_cell} cells with valid PSI): {df_psi_raw_filter.shape[1]}")

    #get sample mean
    df_psi_raw_filter['sample_mean'] = df_psi_raw_filter.mean(axis=1, skipna=True)
    dict_cluster = dict(zip(adata.obs.index, adata.obs[cluster_name]))
    df_psi_raw_filter[cluster_name] = df_psi_raw_filter.index.map(dict_cluster)
    
    df_value_mean = df_psi_raw_filter.groupby([cluster_name], as_index=False)["sample_mean"].mean()
    df_value_mean.sample_mean = df_value_mean.sample_mean.round(3).astype(str)
    dict_mean_cluster = dict(zip(df_value_mean[cluster_name], df_value_mean.sample_mean))
    # print("average psi values per cluster:", dict_mean_cluster)
    
    df_psi_raw_filter = df_psi_raw_filter.replace({cluster_name: dict_mean_cluster})
    df_psi_raw_filter = df_psi_raw_filter.rename(columns={cluster_name: "psi_mean_cluster"})
    
    df_psi_mod_mean_cluster = df_psi_raw_filter.apply(
        lambda row: row.fillna(row["psi_mean_cluster"]),
        axis=1
    )

    cols = df_psi_mod_mean_cluster.columns
    df_psi_mod_mean_cluster[cols] = df_psi_mod_mean_cluster[cols].apply(pd.to_numeric, errors='coerce')

    df_psi_mod_mean_cluster = df_psi_mod_mean_cluster.drop(columns=["sample_mean","psi_mean_cluster"])
    
    ### Create a new adata
    obs_sub = adata.obs.reindex(df_psi_mod_mean_cluster.index).copy()
    var_sub = adata.var.reindex(df_psi_mod_mean_cluster.columns).copy()

    adata_new = anndata.AnnData(
        X=df_psi_mod_mean_cluster.values,
        obs=obs_sub,
        var=var_sub
    )

    adata_new.write(os.path.join(final_out_dir, out_name+"_PSI_DAS.h5ad")) 
    
    return adata_new
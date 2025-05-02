import os
import pandas as pd
import numpy as np
import anndata
import matplotlib.pyplot as plt

def prepare_data_for_differential_AS_analysis(
    input_h5ad_path: str,
    outrigger_path: str,
    output_path: str,
    output_name: str,
    min_cells: int = 10,
    cluster_column: str = "predicted.celltype.l1"
):
    """
    Prepares PSI matrix and metadata for differential alternative splicing analysis.

    This function filters splicing events based on minimum cell coverage, fills missing
    PSI values with cluster-mean imputation, maps AS events to gene names, and writes
    the processed data to an AnnData `.h5ad` file.

    Parameters
    ----------
    input_h5ad_path : str
        Path to the input AnnData `.h5ad` file containing raw PSI matrix.
    outrigger_path : str
        Path to the Outrigger output directory (containing event metadata under `index/`).
    output_path : str
        Directory where the processed `.h5ad` file will be saved.
    output_name : str
        Output filename (without extension) for the resulting AnnData file.
    min_cells : int, optional
        Minimum number of cells in which a splicing event must be detected (default is 10).
    cluster_column : str, optional
        Column in `adata.obs` to use as cell cluster label for averaging (default: "predicted.celltype.l1").

    Returns
    -------
    adata_leiden : anndata.AnnData
        AnnData object containing the filtered and imputed PSI matrix with gene annotations.

    Notes
    -----
    - This function includes both SE and MXE events from Outrigger.
    - Missing PSI values are filled using cluster-level mean values.
    - The returned `.h5ad` file is suitable for downstream differential AS testing.

    Example
    -------
    >>> prepare_as_data_for_differential_analysis(
            input_h5ad_path="../data/fsla_PSI_N10_exon.h5ad",
            outrigger_path="./outrigger_output/",
            output_path="../data/",
            output_name="fsla_PSI_N10_exon_GO"
        )
    """
    adata = anndata.read_h5ad(input_h5ad_path)
    df_psi_raw = adata.to_df()

    # Step 1: Filter events present in at least `min_cells` cells
    df_psi_raw_t = df_psi_raw.T
    df_psi_raw_t["count"] = df_psi_raw_t.ge(0.0).sum(axis=1)
    df_psi_raw_t_filter = df_psi_raw_t[df_psi_raw_t["count"] >= min_cells].drop(columns=["count"])
    df_psi_raw_filter = df_psi_raw_t_filter.T

    # Step 2: Fill missing PSI values by cluster mean
    df_psi_raw_filter['sample_mean'] = df_psi_raw_filter.mean(axis=1, skipna=True)
    df_psi_raw_filter = pd.merge(df_psi_raw_filter, pd.DataFrame(adata.obs), left_index=True, right_index=True)
    df_value_mean = df_psi_raw_filter.groupby([cluster_column], as_index=False)["sample_mean"].mean()
    df_value_mean.sample_mean = df_value_mean.sample_mean.round(3).astype(str)
    dict_mean_cluster = dict(zip(df_value_mean[cluster_column], df_value_mean.sample_mean))

    df_psi_raw_filter = df_psi_raw_filter.replace({cluster_column: dict_mean_cluster})
    df_psi_raw_filter = df_psi_raw_filter.rename(columns={cluster_column: "psi_mean_cluster"})

    df_psi_mod_mean_cluster = df_psi_raw_filter.apply(lambda x: x.fillna(value=df_psi_raw_filter["psi_mean_cluster"]))
    df_psi_mod_mean_cluster[df_psi_mod_mean_cluster.columns] = df_psi_mod_mean_cluster[df_psi_mod_mean_cluster.columns].apply(pd.to_numeric, errors='coerce')

    df_psi_mod_mean_cluster = df_psi_mod_mean_cluster.drop(columns=["sample_mean", "psi_mean_cluster"])
    df_obs_go_mean_cluster = pd.merge(pd.DataFrame(adata.obs), df_psi_mod_mean_cluster, left_index=True, right_index=True)
    df_obs_go_mean_cluster["CB"] = df_obs_go_mean_cluster.index

    # Step 3: Load SE and MXE event annotation and map to gene names
    pd_mxe_event = pd.read_csv(os.path.join(outrigger_path, "index", "mxe", "validated", "events.csv"))
    pd_se_event = pd.read_csv(os.path.join(outrigger_path, "index", "se", "validated", "events.csv"))

    pd_mxe_event["AS_event_type"] = "MXE"
    pd_se_event["AS_event_type"] = "SE"
    pd_event = pd.concat([pd_mxe_event, pd_se_event], ignore_index=True)

    pd_event["isoform1_gene_name_mod"] = pd_event["isoform1_gene_name"].fillna(pd_event["isoform1_gene_id"])
    pd_event["isoform2_gene_name_mod"] = pd_event["isoform2_gene_name"].fillna(pd_event["isoform2_gene_id"])

    pd_event_iso1_freq = pd_event[["event_id", "isoform1_gene_name_mod"]].groupby(
        ["event_id", "isoform1_gene_name_mod"], dropna=False
    ).size().reset_index(name="count1").sort_values(["event_id", "count1"], ascending=False).groupby("event_id").head(1)

    pd_event_iso2_freq = pd_event[["event_id", "isoform2_gene_name_mod"]].groupby(
        ["event_id", "isoform2_gene_name_mod"], dropna=False
    ).size().reset_index(name="count2").sort_values(["event_id", "count2"], ascending=False).groupby("event_id").head(1)

    pd_event_gene = pd.merge(pd_event_iso1_freq, pd_event_iso2_freq, on="event_id", how="outer")
    pd_event_gene["gene_name"] = np.select(
        [
            (pd_event_gene["isoform1_gene_name_mod"].notna()) & (pd_event_gene["isoform1_gene_name_mod"] == pd_event_gene["isoform2_gene_name_mod"]),
            (pd_event_gene["isoform1_gene_name_mod"].notna()) & (pd_event_gene["isoform2_gene_name_mod"].isna()),
            (pd_event_gene["isoform2_gene_name_mod"].notna()) & (pd_event_gene["isoform1_gene_name_mod"].isna()),
            (pd_event_gene["isoform1_gene_name_mod"].notna()) & (pd_event_gene["isoform2_gene_name_mod"].notna()) & (pd_event_gene["isoform1_gene_name_mod"] != pd_event_gene["isoform2_gene_name_mod"]),
            (pd_event_gene["isoform2_gene_name_mod"].isna()) & (pd_event_gene["isoform1_gene_name_mod"].isna())
        ],
        [
            pd_event_gene["isoform1_gene_name_mod"],
            pd_event_gene["isoform1_gene_name_mod"],
            pd_event_gene["isoform2_gene_name_mod"],
            pd_event_gene["isoform1_gene_name_mod"] + "," + pd_event_gene["isoform2_gene_name_mod"],
            "Empty"
        ]
    )
    pd_event_gene["gene_name"] = pd_event_gene["gene_name"].apply(lambda x: ",".join(sorted(set(x.split(",")))) if "," in x else x)
    dict_event_gene = dict(zip(pd_event_gene.event_id, pd_event_gene.gene_name))

    # Step 4: Construct AnnData
    obs = df_obs_go_mean_cluster[["CB", cluster_column]]
    var_names = df_obs_go_mean_cluster.columns[3:-3]
    var = pd.DataFrame(index=var_names)
    var["gene_name"] = var.index
    var.replace({"gene_name": dict_event_gene}, inplace=True)

    X = df_obs_go_mean_cluster.iloc[:, 3:-3].values
    adata_leiden = anndata.AnnData(X=X, obs=obs, var=var, dtype=np.float32)
    adata_leiden.write(os.path.join(output_path, output_name + ".h5ad"))

    return adata_leiden

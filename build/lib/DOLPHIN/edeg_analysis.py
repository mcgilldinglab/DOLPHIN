import pandas as pd
import numpy as np
import statsmodels.stats.multitest as smm
from DOLPHIN.preprocess import gtfpy
import scanpy as sc
from scipy.stats import norm
import os

def generate_mast_weighted_pvals(
    mast_input: str,
    gtf_path: str,
    feature_matrix: str,
    leiden_res: str = "leiden_0_4_8",
    leiden_cluster_num: str = "2",
    pval_column_prefix: str = "MAST"
) -> pd.DataFrame:
    """
    Processes MAST results and computes weighted gene-level p-values using the Stouffer method.

    This function loads exon-level MAST results, maps exons to genomic coordinates using a GTF file and a
    feature matrix (AnnData), calculates exon lengths and weights, and applies the weighted Stouffer method
    to derive gene-level p-values. Adjusted p-values are computed using FDR and Bonferroni correction.

    Parameters
    ----------
    mast_input : str
        Path to the CSV file containing exon-level MAST results from Seurat.
    gtf_path : str
        Path to the GTF annotation file.
    feature_matrix : str
        Path to the `.h5ad` file containing the exon-level feature matrix with gene annotations.
    leiden_res : str, optional
        Name of the Leiden clustering resolution used (not used in calculation, for record only).
    leiden_cluster_num : str, optional
        Cluster ID used in the MAST result (not used in calculation, for record only).
    pval_column_prefix : str, optional
        Prefix used for naming final p-value columns (default: "MAST").

    Returns
    -------
    pd_MAST : pandas.DataFrame
        Processed DataFrame with added columns for exon length, weights, weighted Stouffer p-values,
        and multiple testing corrected adjusted p-values.

    Notes
    -----
    - Assumes the exon ID format is `geneID-exonNumber`.
    - Missing exon mappings are silently ignored (will be NaN in output).
    """

    # Load MAST result
    pd_MAST = pd.read_csv(mast_input)
    pd_MAST = pd_MAST.rename(columns={"Unnamed: 0": "Exon_names"})
    pd_MAST['Exon_names'] = pd_MAST['Exon_names'].str.replace('/', '_')
    pd_MAST['Gene_names'] = pd_MAST['Exon_names'].apply(lambda x: x[:x.rfind('-')] if '-' in x else x)

    pd_MAST = pd_MAST.rename(columns={
        "Gene_names": "Cancer_gene_names",
        "Exon_names": "Cancer_names",
        "p_val": f"{pval_column_prefix}_pvals",
        "p_val_adj": f"{pval_column_prefix}_pvals_adj",
        "avg_log2FC": f"{pval_column_prefix}_logfoldchanges",
        "pct.1": f"{pval_column_prefix}_pct.1",
        "pct.2": f"{pval_column_prefix}_pct.2"
    })

    # Load GTF and create exon key
    df_gtf_orig = gtfpy.readGTF(gtf_path)
    attribute_split = df_gtf_orig['attribute'].str.split(';', expand=True)
    df_gtf = pd.concat([df_gtf_orig, attribute_split], axis=1)[[0, 2, 5, "start", "end"]]
    df_gtf["gene_id"] = df_gtf[0].str.extract(r'"(.*?)"')
    df_gtf["exon_num"] = df_gtf[5].str.extract(r'"(.*?)"')
    df_gtf["_new"] = df_gtf["gene_id"] + "-" + df_gtf["exon_num"]

    # Load feature matrix to build index mapping
    adata = sc.read_h5ad(feature_matrix)
    adata.var["_new"] = adata.var.groupby('gene_id').cumcount() + 1
    adata.var['_new'] = adata.var['gene_id'].astype(str) + '-' + adata.var['_new'].astype(str)

    dict_start = dict(zip(df_gtf['_new'], df_gtf['start']))
    dict_end = dict(zip(df_gtf['_new'], df_gtf['end']))
    adata.var['start'] = adata.var['_new'].map(dict_start)
    adata.var['end'] = adata.var['_new'].map(dict_end)

    dict_map_start = dict(zip(adata.var.index, adata.var["start"]))
    dict_map_end = dict(zip(adata.var.index, adata.var["end"]))

    # Map exon coordinates to MAST
    pd_MAST["start"] = pd_MAST["Cancer_names"].map(dict_map_start)
    pd_MAST["end"] = pd_MAST["Cancer_names"].map(dict_map_end)
    pd_MAST['start'] = pd.to_numeric(pd_MAST['start'], errors='coerce')
    pd_MAST['end'] = pd.to_numeric(pd_MAST['end'], errors='coerce')
    pd_MAST['exon_length'] = pd_MAST['end'] - pd_MAST['start'] + 1

    # Weight calculation: exon length normalized per gene
    pd_MAST[f'{pval_column_prefix}_exon_length'] = np.where(
        pd_MAST[f"{pval_column_prefix}_pvals"].isna(), np.nan, pd_MAST['exon_length']
    )
    pd_MAST[f'{pval_column_prefix}_exon_weight'] = pd_MAST.groupby('Cancer_gene_names')[
        f'{pval_column_prefix}_exon_length'
    ].transform(lambda x: x / x.sum())

    # Weighted Stouffer method per gene
    def weighted_stouffer_method(p_values, weights):
        non_nan_pvals = p_values[~np.isnan(p_values)]
        non_nan_weights = weights[~np.isnan(weights)]

        if len(non_nan_pvals) == 0 or len(non_nan_pvals) != len(non_nan_weights):
            return np.nan
        z_scores = norm.isf(non_nan_pvals)
        combined_z = np.sum(non_nan_weights * z_scores) / np.sqrt(np.sum(non_nan_weights ** 2))
        return norm.sf(combined_z)

    gene_pvals = pd_MAST.groupby('Cancer_gene_names').apply(
        lambda x: weighted_stouffer_method(
            x[f'{pval_column_prefix}_pvals'].values,
            x[f'{pval_column_prefix}_exon_weight'].values
        )
    ).reset_index(name=f'{pval_column_prefix}_weighted_stouffer_pval')

    dict_weighted_p = dict(zip(gene_pvals['Cancer_gene_names'], gene_pvals[f'{pval_column_prefix}_weighted_stouffer_pval']))
    pd_MAST[f'{pval_column_prefix}_weighted_stouffer_pval'] = pd_MAST['Cancer_gene_names'].map(dict_weighted_p)

    # Adjust with FDR and Bonferroni
    df_gene = pd_MAST.dropna(subset=[f'{pval_column_prefix}_weighted_stouffer_pval']).drop_duplicates('Cancer_gene_names')
    _, adj_fdr, _, _ = smm.multipletests(df_gene[f'{pval_column_prefix}_weighted_stouffer_pval'], method='fdr_bh')
    _, adj_bonf, _, _ = smm.multipletests(df_gene[f'{pval_column_prefix}_weighted_stouffer_pval'], method='bonferroni')

    dict_fdr = dict(zip(df_gene["Cancer_gene_names"], adj_fdr))
    dict_bonf = dict(zip(df_gene["Cancer_gene_names"], adj_bonf))

    pd_MAST[f"{pval_column_prefix}_weighted_stouffer_pval_adj"] = pd_MAST["Cancer_gene_names"].map(dict_fdr)
    pd_MAST[f"{pval_column_prefix}_weighted_stouffer_pval_adj_bonf"] = pd_MAST["Cancer_gene_names"].map(dict_bonf)

    return pd_MAST

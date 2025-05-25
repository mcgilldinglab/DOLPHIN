import pandas as pd
import numpy as np
import statsmodels.stats.multitest as smm
import scanpy as sc
from scipy.stats import norm

def stouffer_method(pvals):
    non_nan_pvals = pvals[~np.isnan(pvals)]
    
    if (len(non_nan_pvals) == 0):
        return np.nan
    
    # Convert p-values to z-scores
    z_scores = norm.isf(non_nan_pvals)  # Inverse survival function (equivalent to 1 - CDF)
    
    # Calculate the unweighted z-score
    combined_z = np.sum(z_scores) / np.sqrt(len(non_nan_pvals))
    
    # Convert the combined z-score back to a p-value
    combined_p_value = norm.sf(combined_z)  # Survival function (1 - CDF)
    
    return combined_p_value

def weighted_avg_log2fc(log2fc, weights):
    non_nan_log2fc = log2fc[~np.isnan(log2fc)]
    non_nan_weights = weights[~np.isnan(weights)]
    
    if (len(non_nan_log2fc) == 0):
        return np.nan
    if (len(non_nan_log2fc) != len(non_nan_weights)):
        return np.nan
    
    return np.sum(np.abs(non_nan_log2fc) * non_nan_weights)

def run_jdeg(seurat_output, output):
    """
    Aggregate junction-level marker results into gene-level statistics using stouffer method.

    This function converts junction-level marker results from Seurat using MAST into gene-level statistical
    insights. It applies the Stouffer method to combine junction-level p-values, and computes average absolute log2 fold changes.

    Parameters
    ----------
    seurat_output : str
        Path to the CSV file containing exon-level differential expression results from Seurat.
        
    output : str
        Path where the final gene-level results will be saved as a CSV file.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing gene-level statistics including:
        - MAST_abs_avg_log2FC: Average of absolute log2 fold changes per gene.
        - MAST_stouffer_pval: Combined p-value using the Stouffer method.
        - MAST_stouffer_pval_adj_bonf: Bonferroni-adjusted Stouffer p-value.
        
    """
    
    pd_MAST = pd.read_csv(seurat_output)
    pd_MAST = pd_MAST.rename(columns={"Unnamed: 0":"Exon_names"})
    pd_MAST['Exon_names'] = pd_MAST['Exon_names'].str.replace('/', '_')
    pd_MAST['Gene_names'] = pd_MAST['Exon_names'].apply(lambda x: x[:x.rfind('-')] if '-' in x else x)

    pd_MAST['MAST_abs_avg_log2FC'] = pd_MAST.groupby('Gene_names')['avg_log2FC'].transform(lambda x: x.dropna().abs().mean())
    
    temp_stouffer_p = pd_MAST.groupby('Gene_names').apply(
        lambda x: stouffer_method(x["p_val"].values)
    ).reset_index(name='temp_p_value')

    dict_temp_p_value = dict(zip(temp_stouffer_p['Gene_names'], temp_stouffer_p['temp_p_value']))

    pd_MAST['MAST_stouffer_pval'] = pd_MAST['Gene_names'].map(dict_temp_p_value)
    
    ## adjusted weighted fisher's p-value    
    df_sub_temp = pd_MAST.dropna(subset=['MAST_stouffer_pval']).copy()
    df_sub_temp = df_sub_temp.drop_duplicates(subset=['Gene_names'])
    _, adjusted_pvals_bonf, _, _ = smm.multipletests(df_sub_temp['MAST_stouffer_pval'], method='bonferroni')
    
    df_sub_temp["temp_adj_p_bonf"] = adjusted_pvals_bonf
    
    dict_p_adj_bonf = dict(zip(df_sub_temp["Gene_names"], df_sub_temp["temp_adj_p_bonf"]))

    pd_MAST["MAST_stouffer_pval_adj_bonf"] = pd_MAST["Gene_names"].map(dict_p_adj_bonf)
                     
    pd_MAST.to_csv(output, index=False)
    
    return pd_MAST
    
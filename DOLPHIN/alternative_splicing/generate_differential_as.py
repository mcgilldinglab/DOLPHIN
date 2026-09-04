import os

import anndata
import numpy as np
import pandas as pd
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests

try:
    from ._anndata_compat import enable_nullable_string_writes, write_h5ad_preserve_strings
except ImportError:  # Support direct execution from the module directory.
    from _anndata_compat import enable_nullable_string_writes, write_h5ad_preserve_strings


def impute_psi_by_cell_type(psi, cell_types):
    """Impute each event from its mean among cells of the same cell type."""
    cell_types = pd.Series(cell_types, index=psi.index)
    if cell_types.isna().any():
        missing = cell_types.index[cell_types.isna()].tolist()[:10]
        raise ValueError(
            "Cell-type labels are required for differential-AS imputation. "
            f"Missing labels for cells: {missing}"
        )
    event_means = psi.groupby(cell_types, observed=True, dropna=False).transform("mean")
    return psi.fillna(event_means)


def calculate_differential_as(psi, groups, group1=None, group2=None):
    """Run two-sided Wilcoxon rank-sum tests and BH correction per event."""
    groups = pd.Series(groups, index=psi.index).astype("string")
    if groups.isna().any():
        missing = groups.index[groups.isna()].tolist()[:10]
        raise ValueError(
            "Biological group labels are required for differential AS. "
            f"Missing labels for cells: {missing}"
        )
    available_groups = [str(value) for value in pd.unique(groups.dropna())]
    if group1 is None and group2 is None:
        if len(available_groups) != 2:
            raise ValueError(
                "Differential AS requires exactly two groups when --das-group1 and "
                "--das-group2 are not supplied. "
                f"Observed groups: {available_groups}"
            )
        group1, group2 = available_groups
    elif group1 is None or group2 is None:
        raise ValueError("Both DAS comparison groups must be supplied together.")

    group1 = str(group1)
    group2 = str(group2)
    if group1 == group2:
        raise ValueError("DAS comparison groups must be different.")
    missing_groups = [
        group for group in (group1, group2) if group not in available_groups
    ]
    if missing_groups:
        raise ValueError(
            f"DAS groups not found in metadata: {missing_groups}. "
            f"Observed groups: {available_groups}"
        )

    group1_mask = groups == group1
    group2_mask = groups == group2
    rows = []
    for event_id in psi.columns:
        values1 = psi.loc[group1_mask, event_id].dropna().to_numpy(dtype=float)
        values2 = psi.loc[group2_mask, event_id].dropna().to_numpy(dtype=float)
        if values1.size and values2.size:
            statistic, p_value = ranksums(values1, values2, alternative="two-sided")
        else:
            statistic, p_value = np.nan, np.nan
        mean1 = float(np.mean(values1)) if values1.size else np.nan
        mean2 = float(np.mean(values2)) if values2.size else np.nan
        rows.append(
            {
                "event_id": event_id,
                "group1": group1,
                "group2": group2,
                "n_group1": int(values1.size),
                "n_group2": int(values2.size),
                "mean_psi_group1": mean1,
                "mean_psi_group2": mean2,
                "delta_psi": mean1 - mean2,
                "wilcoxon_statistic": statistic,
                "p_value": p_value,
            }
        )

    results = pd.DataFrame(rows)
    results["p_value_adj_bh"] = np.nan
    finite = np.isfinite(results["p_value"].to_numpy(dtype=float))
    if finite.any():
        results.loc[finite, "p_value_adj_bh"] = multipletests(
            results.loc[finite, "p_value"], method="fdr_bh"
        )[1]
    results["significant_bh_0_05"] = results["p_value_adj_bh"] < 0.05
    return results


def run_differential_as(
    outrigger_psi_data: str,
    out_name: str,
    cluster_name: str,
    out_directory: str = "./",
    n_cell: int = 10,
    group_column: str = "Condition",
    group1: str = None,
    group2: str = None,
):
    """Impute PSI by event/cell type and run differential AS testing."""
    final_out_dir = os.path.join(out_directory, "alternative_splicing")
    os.makedirs(final_out_dir, exist_ok=True)

    adata = anndata.read_h5ad(outrigger_psi_data)
    for column, purpose in (
        (cluster_name, "cell-type imputation"),
        (group_column, "differential-AS comparison"),
    ):
        if column not in adata.obs:
            raise ValueError(
                f"Metadata column '{column}' required for {purpose} is missing "
                f"from {outrigger_psi_data}."
            )

    psi_raw = adata.to_df().apply(pd.to_numeric, errors="coerce")
    invalid = ((psi_raw < 0) | (psi_raw > 1)).any(axis=None)
    if invalid:
        raise ValueError("PSI values must be within [0, 1] or NaN.")

    print(f"Total number of splicing events before filtering: {psi_raw.shape[1]}")
    keep_events = psi_raw.notna().sum(axis=0) >= int(n_cell)
    psi_filtered = psi_raw.loc[:, keep_events].copy()
    print(
        f"Number of splicing events after filtering (>= {n_cell} cells with valid PSI): "
        f"{psi_filtered.shape[1]}"
    )
    if psi_filtered.shape[1] == 0:
        raise ValueError(
            f"No splicing events were detected in at least {n_cell} cells."
        )

    psi_imputed = impute_psi_by_cell_type(
        psi_filtered,
        adata.obs.loc[psi_filtered.index, cluster_name],
    )
    results = calculate_differential_as(
        psi_imputed,
        adata.obs.loc[psi_imputed.index, group_column],
        group1=group1,
        group2=group2,
    )

    var_sub = adata.var.reindex(psi_imputed.columns).copy()
    annotations = var_sub.copy()
    annotations["event_id"] = annotations.index
    annotations = annotations[["event_id"] + [
        column for column in annotations.columns if column != "event_id"
    ]]
    results = results.merge(annotations.reset_index(drop=True), on="event_id", how="left")

    das_results_path = os.path.join(final_out_dir, f"{out_name}_DAS.csv")
    results.to_csv(das_results_path, index=False)

    adata_new = anndata.AnnData(
        X=psi_imputed.to_numpy(dtype=np.float32),
        obs=adata.obs.reindex(psi_imputed.index).copy(),
        var=var_sub,
    )
    adata_new.uns["differential_as"] = {
        "cell_type_column": cluster_name,
        "group_column": group_column,
        "group1": str(results["group1"].iloc[0]),
        "group2": str(results["group2"].iloc[0]),
        "method": "two-sided Wilcoxon rank-sum with Benjamini-Hochberg correction",
        "results_file": os.path.abspath(das_results_path),
    }

    enable_nullable_string_writes()
    write_h5ad_preserve_strings(
        adata_new,
        os.path.join(final_out_dir, f"{out_name}_PSI_DAS.h5ad"),
    )
    return adata_new

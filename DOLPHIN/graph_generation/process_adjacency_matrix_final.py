import os

import anndata
import numpy as np
import pandas as pd
from scipy import sparse

from ._anndata_compat import enable_nullable_string_writes, write_h5ad_preserve_strings


def _to_csr_matrix(matrix):
    if sparse.issparse(matrix):
        return matrix.tocsr()
    return sparse.csr_matrix(np.asarray(matrix))


def _contiguous_slices(values):
    values = np.asarray(values, dtype=object)
    if values.size == 0:
        return []

    slices = []
    start = 0
    current = values[0]
    for idx in range(1, values.size + 1):
        if idx == values.size or values[idx] != current:
            slices.append((current, start, idx))
            if idx < values.size:
                start = idx
                current = values[idx]
    return slices


def _gene_numeric_order(gene_id: str) -> int:
    return int(str(gene_id)[4:])


def _feature_exon_order(var_name: str) -> int:
    return int(str(var_name).rsplit("-", 1)[-1])


def _build_exon_keep_lookup(adata_fea_gtf):
    feature_matrix = _to_csr_matrix(adata_fea_gtf.X)
    feature_var = pd.DataFrame(adata_fea_gtf.var)
    feature_var = feature_var.copy()
    feature_var["var_name"] = list(adata_fea_gtf.var_names)
    feature_var["orig_pos"] = np.arange(feature_var.shape[0])
    feature_var["gene_order"] = feature_var["gene_id"].map(_gene_numeric_order)
    feature_var["exon_order"] = feature_var["var_name"].map(_feature_exon_order)
    feature_var = feature_var.sort_values(
        by=["gene_order", "exon_order"],
        kind="stable",
    ).reset_index(drop=True)
    feature_nonzero = np.asarray((feature_matrix != 0).sum(axis=0)).ravel() > 0
    feature_var["is_expressed"] = feature_nonzero[feature_var["orig_pos"].to_numpy()]

    exon_keep_lookup = {}
    for gene_id, group in feature_var.groupby("gene_id", sort=False):
        positions = np.flatnonzero(group["is_expressed"].to_numpy()) + 1
        exon_keep_lookup[gene_id] = {str(position) for position in positions}
    return exon_keep_lookup


def _build_keep_mask(adata_adj, exon_keep_lookup):
    adj_matrix = _to_csr_matrix(adata_adj.X)
    adj_var = pd.DataFrame(adata_adj.var)
    adj_var = adj_var.copy()
    adj_var["var_name"] = list(adata_adj.var_names)

    adj_nonzero = np.asarray(adj_matrix.max(axis=0).toarray()).ravel()
    keep_mask = np.ones(adj_var.shape[0], dtype=bool)

    for gene_id, start, stop in _contiguous_slices(adj_var["gene_id"].to_numpy()):
        gene_var = adj_var.iloc[start:stop]
        gene_size = int(np.sqrt(stop - start))
        if gene_size * gene_size != (stop - start):
            raise ValueError(
                f"Adjacency vector length {stop - start} for gene {gene_id} "
                "is not a perfect square."
            )

        lex_order = np.argsort(gene_var["var_name"].to_numpy(), kind="stable")
        sorted_values = adj_nonzero[start:stop][lex_order].reshape((gene_size, gene_size))

        zero_rows = np.all(sorted_values == 0, axis=1)
        zero_cols = np.all(sorted_values == 0, axis=0)
        keep_ids = exon_keep_lookup.get(gene_id, set())

        drop_exons = np.array(
            [
                zero_rows[idx] and zero_cols[idx] and str(idx + 1) not in keep_ids
                for idx in range(gene_size)
            ],
            dtype=bool,
        )
        if np.any(drop_exons):
            flag_matrix = np.zeros((gene_size, gene_size), dtype=bool)
            flag_matrix[drop_exons, :] = True
            flag_matrix[:, drop_exons] = True
            keep_mask[start:stop] = ~flag_matrix.reshape(-1)

    return keep_mask


def _build_output_var(adj_var: pd.DataFrame, keep_mask: np.ndarray) -> pd.DataFrame:
    filtered_var = adj_var.loc[keep_mask, ["gene_id", "gene_name"]].copy()
    filtered_var["gene_id"] = filtered_var["gene_id"].astype(str)
    filtered_var["gene_name"] = filtered_var["gene_name"].astype(str)
    new_indices = []
    current_gene = None
    counter = 0
    for gene_id in filtered_var["gene_id"].to_numpy():
        if gene_id != current_gene:
            current_gene = gene_id
            counter = 1
        else:
            counter += 1
        new_indices.append(f"{gene_id}-{counter}")
    filtered_var.index = new_indices
    filtered_var.index.name = None
    return filtered_var


def run_adjacency_matrix_final(
    out_name: str,
    out_directory: str = "./",
    batch_size=1000):
    """
    Generates the final adjacency matrix by filtering out invalid edges based on FeatureComp data.
    """
    print("Start Generating Final Adjacency Matrix...")

    final_out_dir = os.path.join(out_directory, "data")

    adata_adj = anndata.read_h5ad(os.path.join(final_out_dir, f"AdjacencyComp_{out_name}.h5ad"))
    adata_fea_gtf = anndata.read_h5ad(os.path.join(final_out_dir, f"Feature_{out_name}.h5ad"))

    exon_keep_lookup = _build_exon_keep_lookup(adata_fea_gtf)
    adj_var = pd.DataFrame(adata_adj.var).copy()
    keep_mask = _build_keep_mask(adata_adj, exon_keep_lookup)

    adj_matrix = _to_csr_matrix(adata_adj.X)
    filtered_matrix = adj_matrix[:, keep_mask].tocsr()
    filtered_var = _build_output_var(adj_var, keep_mask)
    obs = pd.DataFrame(adata_adj.obs).copy()

    adata = anndata.AnnData(X=filtered_matrix, obs=obs, var=filtered_var)
    enable_nullable_string_writes()
    write_h5ad_preserve_strings(
        adata,
        os.path.join(final_out_dir, "AdjacencyCompRe_" + out_name + ".h5ad"),
    )

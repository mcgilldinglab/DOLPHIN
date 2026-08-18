from multiprocessing import Process
import os
import shutil
import warnings

import anndata
import numpy as np
import pandas as pd
from anndata._core.views import ImplicitModificationWarning
from scipy import sparse

from ._anndata_compat import enable_nullable_string_writes

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=ImplicitModificationWarning)


_ADJ_COMP_CONTEXT = {}


def _as_csr_matrix(matrix):
    if sparse.issparse(matrix):
        return matrix.tocsr()
    return sparse.csr_matrix(np.asarray(matrix))


def _build_contiguous_gene_slices(gene_ids):
    gene_ids = np.asarray(gene_ids, dtype=object)
    if gene_ids.size == 0:
        return []

    slices = []
    start = 0
    current_gene = gene_ids[0]
    for idx in range(1, gene_ids.size + 1):
        if idx == gene_ids.size or gene_ids[idx] != current_gene:
            slices.append((current_gene, start, idx))
            if idx < gene_ids.size:
                start = idx
                current_gene = gene_ids[idx]
    return slices


def _build_gene_tasks(feature_var, adj_var, feature_matrix):
    feature_gene_slices = _build_contiguous_gene_slices(feature_var["gene_id"].to_numpy())
    adj_gene_slices = {
        gene_id: (start, stop)
        for gene_id, start, stop in _build_contiguous_gene_slices(adj_var["gene_id"].to_numpy())
    }
    exon_nonzero_mask = np.asarray(feature_matrix.max(axis=0).toarray()).ravel() != 0

    gene_tasks = []
    for gene_id, fea_start, fea_stop in feature_gene_slices:
        exon_count = fea_stop - fea_start
        if exon_count <= 1:
            continue
        if not np.any(exon_nonzero_mask[fea_start:fea_stop]):
            continue
        adj_bounds = adj_gene_slices.get(gene_id)
        if adj_bounds is None:
            continue
        adj_start, adj_stop = adj_bounds
        adj_vector_size = adj_stop - adj_start
        gene_size = int(np.sqrt(adj_vector_size))
        if gene_size * gene_size != adj_vector_size:
            raise ValueError(
                f"Adjacency vector length {adj_vector_size} for gene {gene_id} "
                "is not a perfect square."
            )
        gene_tasks.append((fea_start, fea_stop, adj_start, adj_stop, gene_size))
    return gene_tasks


def _compress_one_cell(cell_index):
    context = _ADJ_COMP_CONTEXT
    sample_id = context["sample_ids"][cell_index]
    print(f"Starting processing cell: {sample_id}")

    feature_row = context["feature_matrix"].getrow(cell_index)
    adj_row = context["adj_matrix"].getrow(cell_index)
    feature_indices = feature_row.indices
    feature_values = feature_row.data
    adj_indices = adj_row.indices
    adj_values = adj_row.data
    compressed_indices = []
    compressed_values = []

    for fea_start, fea_stop, adj_start, adj_stop, gene_size in context["gene_tasks"]:
        feature_left = np.searchsorted(feature_indices, fea_start, side="left")
        feature_right = np.searchsorted(feature_indices, fea_stop, side="left")
        if feature_right <= feature_left:
            continue
        local_feature_indices = feature_indices[feature_left:feature_right] - fea_start
        local_feature_values = feature_values[feature_left:feature_right]
        non_zero_index = local_feature_indices[local_feature_values > 0]
        if non_zero_index.size < 2:
            continue

        local_adj_flat = np.zeros(adj_stop - adj_start, dtype=context["adj_dtype"])
        adj_left = np.searchsorted(adj_indices, adj_start, side="left")
        adj_right = np.searchsorted(adj_indices, adj_stop, side="left")
        if adj_right > adj_left:
            local_adj_indices = adj_indices[adj_left:adj_right] - adj_start
            local_adj_flat[local_adj_indices] = adj_values[adj_left:adj_right]
        old_adj = local_adj_flat.reshape((gene_size, gene_size))
        new_adj = np.zeros((gene_size, gene_size), dtype=context["adj_dtype"])

        pair_row, pair_col = np.triu_indices(non_zero_index.size, k=1)
        if pair_row.size:
            ii = non_zero_index[pair_row]
            jj = non_zero_index[pair_col]
            new_adj[ii, jj] = old_adj[ii, jj]

        consecutive_i = non_zero_index[:-1]
        consecutive_j = non_zero_index[1:]
        if consecutive_i.size:
            zero_mask = old_adj[consecutive_i, consecutive_j] == 0
            if np.any(zero_mask):
                new_adj[consecutive_i[zero_mask], consecutive_j[zero_mask]] = 1

        local_nonzero = np.flatnonzero(new_adj)
        if local_nonzero.size:
            compressed_indices.append(adj_start + local_nonzero)
            compressed_values.append(new_adj.reshape(-1)[local_nonzero])

    if compressed_indices:
        final_indices = np.concatenate(compressed_indices).astype(np.int32, copy=False)
        final_values = np.concatenate(compressed_values).astype(context["adj_dtype"], copy=False)
        indptr = np.array([0, final_indices.size], dtype=np.int32)
        row_matrix = sparse.csr_matrix(
            (final_values, final_indices, indptr),
            shape=(1, context["adj_num_vars"]),
        )
    else:
        row_matrix = sparse.csr_matrix((1, context["adj_num_vars"]), dtype=context["adj_dtype"])

    print(f"Finished processing cell: {sample_id}")
    return row_matrix


def _shard_filename(shard_index, start_idx, end_idx):
    return f"shard_{shard_index:04d}_{start_idx:06d}_{end_idx - 1:06d}.h5ad"


def _process_range(start_idx, end_idx, shard_index):
    total_cells = len(_ADJ_COMP_CONTEXT["sample_ids"])
    end_idx = min(end_idx, total_cells)
    rows = []
    for cell_index in range(start_idx, end_idx):
        rows.append(_compress_one_cell(cell_index))

    if rows:
        shard_matrix = sparse.vstack(rows, format="csr")
    else:
        shard_matrix = sparse.csr_matrix(
            (0, _ADJ_COMP_CONTEXT["adj_num_vars"]),
            dtype=_ADJ_COMP_CONTEXT["adj_dtype"],
        )

    shard_obs = _ADJ_COMP_CONTEXT["obs"].iloc[start_idx:end_idx].copy()
    shard_var = _ADJ_COMP_CONTEXT["var"]
    out_file = os.path.join(
        _ADJ_COMP_CONTEXT["out_path_dir"],
        _shard_filename(shard_index, start_idx, end_idx),
    )
    adata_to_write = anndata.AnnData(
        X=shard_matrix,
        obs=shard_obs,
        var=shard_var,
    )
    enable_nullable_string_writes()
    adata_to_write.write_h5ad(out_file)


def run_adjacency_compression(
    metadata_path,
    out_name,
    out_directory,
    num_processes=25,
    min_cells_per_shard=8,
):
    """
    Compute and update exon-level adjacency matrices per gene for each cell.

    This function loads gene feature and adjacency data, selects genes that have
    more than one exon and non-zero expression, and reconstructs the adjacency
    matrices for each selected gene based on the positions of exons with non-zero values.

    If both exons in an adjacency pair have zero expression, the corresponding
    adjacency value will be set to 0.

    The updated adjacency matrix is saved for each cell in `.h5ad` format.
    """

    adata_adj_orig = anndata.read_h5ad(
        os.path.join(out_directory, "data", f"Adjacency_{out_name}.h5ad")
    )
    adata_fea_orig = anndata.read_h5ad(
        os.path.join(out_directory, "data", f"Feature_{out_name}.h5ad")
    )

    feature_matrix = _as_csr_matrix(adata_fea_orig.X)
    adj_matrix = _as_csr_matrix(adata_adj_orig.X)
    sample_ids = list(pd.read_csv(metadata_path, sep="\t")["CB"])
    cell_index = pd.Index(adata_adj_orig.obs_names)
    sample_indices = [cell_index.get_loc(sample_id) for sample_id in sample_ids]

    feature_var = pd.DataFrame(adata_fea_orig.var)
    adj_var = pd.DataFrame(adata_adj_orig.var)
    gene_tasks = _build_gene_tasks(
        feature_var=feature_var,
        adj_var=adj_var,
        feature_matrix=feature_matrix,
    )

    temp_out_dir = os.path.join(out_directory, "data", "temp")
    old_cell_dir = os.path.join(temp_out_dir, "adj_comp_matrix")
    out_path_dir = os.path.join(temp_out_dir, "adj_comp_shards")
    if os.path.isdir(old_cell_dir):
        shutil.rmtree(old_cell_dir)
    if os.path.isdir(out_path_dir):
        shutil.rmtree(out_path_dir)
    os.makedirs(out_path_dir, exist_ok=True)

    _ADJ_COMP_CONTEXT.clear()
    _ADJ_COMP_CONTEXT.update(
        {
            "feature_matrix": feature_matrix[sample_indices, :].tocsr(),
            "adj_matrix": adj_matrix[sample_indices, :].tocsr(),
            "adj_num_vars": adj_matrix.shape[1],
            "adj_dtype": adj_matrix.dtype,
            "gene_tasks": gene_tasks,
            "sample_ids": sample_ids,
            "obs": pd.DataFrame(adata_adj_orig.obs).loc[sample_ids].copy(),
            "var": pd.DataFrame(adata_adj_orig.var),
            "out_path_dir": out_path_dir,
        }
    )

    num_processes = max(1, min(int(num_processes), len(sample_ids)))
    cells_per_shard = max(
        int(min_cells_per_shard),
        int(np.ceil(len(sample_ids) / num_processes)),
    )
    shard_ranges = []
    for start_idx in range(0, len(sample_ids), cells_per_shard):
        shard_ranges.append((start_idx, min(start_idx + cells_per_shard, len(sample_ids))))

    processes = []
    manifest_rows = []
    for process_idx, (start_idx, end_idx) in enumerate(shard_ranges):
        manifest_rows.append(
            {
                "shard_index": process_idx,
                "start_idx": start_idx,
                "end_idx": end_idx,
                "file_name": _shard_filename(process_idx, start_idx, end_idx),
            }
        )
        proc = Process(target=_process_range, args=(start_idx, end_idx, process_idx))
        processes.append(proc)
        proc.start()

    for proc in processes:
        proc.join()
        if proc.exitcode != 0:
            raise RuntimeError(
                f"Adjacency compression worker exited with code {proc.exitcode}."
            )

    pd.DataFrame(manifest_rows).to_csv(
        os.path.join(out_path_dir, "manifest.tsv"),
        sep="\t",
        index=False,
    )

import gc
import os

import anndata
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
from torch_geometric.data import Data


def _to_csr_matrix(matrix):
    if sp.issparse(matrix):
        return matrix.tocsr()
    return sp.csr_matrix(np.asarray(matrix))


def _row_to_numpy(matrix, row_idx):
    if sp.issparse(matrix):
        return matrix.getrow(row_idx).toarray()
    return np.asarray(matrix[row_idx : row_idx + 1])


def _build_graph_context(final_output_path, output_name, celltypename=None, mapper=None):
    adata_adj_raw = anndata.read_h5ad(
        os.path.join(final_output_path, "AdjacencyCompReHvg_" + output_name + ".h5ad")
    )
    adata_adj = anndata.read_h5ad(
        os.path.join(final_output_path, "AdjacencyCompReHvgEdge_" + output_name + ".h5ad")
    )
    adata_fea = anndata.read_h5ad(
        os.path.join(final_output_path, "FeatureCompHvg_" + output_name + ".h5ad")
    )

    feature_matrix = _to_csr_matrix(adata_fea.X)
    adj_raw_matrix = _to_csr_matrix(adata_adj_raw.X)
    adj_matrix = _to_csr_matrix(adata_adj.X)

    fea_var = pd.DataFrame(adata_fea.var)
    gene_list = list(sorted(set(fea_var["gene_id"])))
    feature_counts = fea_var.groupby("gene_id").size()
    feature_offsets = {}
    offset = 0
    for gene_id in gene_list:
        feature_offsets[gene_id] = offset
        offset += int(feature_counts[gene_id])

    adj_gene_ids = np.asarray(adata_adj.var["gene_id"], dtype=object)
    gene_adj_indices = {}
    gene_adj_sizes = {}
    for gene_id in gene_list:
        indices = np.flatnonzero(adj_gene_ids == gene_id)
        gene_adj_indices[gene_id] = indices
        gene_adj_sizes[gene_id] = int(np.sqrt(indices.size))

    sample_list = list(adata_fea.obs_names)
    if celltypename is not None and mapper is not None and celltypename in adata_fea.obs.columns:
        labels = (
            pd.DataFrame(adata_fea.obs)[celltypename]
            .map(mapper)
            .fillna(0)
            .astype(int)
            .to_numpy()
        )
    else:
        labels = np.zeros(len(adata_fea), dtype=int)

    return {
        "feature_matrix": feature_matrix,
        "adj_raw_matrix": adj_raw_matrix,
        "adj_matrix": adj_matrix,
        "feature_h5ad_path": os.path.join(final_output_path, "FeatureCompHvg_" + output_name + ".h5ad"),
        "adj_raw_h5ad_path": os.path.join(final_output_path, "AdjacencyCompReHvg_" + output_name + ".h5ad"),
        "adj_edge_h5ad_path": os.path.join(final_output_path, "AdjacencyCompReHvgEdge_" + output_name + ".h5ad"),
        "feature_offsets": feature_offsets,
        "gene_adj_indices": gene_adj_indices,
        "gene_adj_sizes": gene_adj_sizes,
        "gene_list": gene_list,
        "sample_list": sample_list,
        "labels": labels,
    }


def _build_cell_edge_arrays(context, sample_idx):
    sample_name = context["sample_list"][sample_idx]
    label = int(context["labels"][sample_idx])
    cell_edge_parts = []
    cell_weight_parts = []

    adj_row = context["adj_matrix"].getrow(sample_idx)
    for gene_id in context["gene_list"]:
        gene_indices = context["gene_adj_indices"][gene_id]
        if gene_indices.size == 0:
            continue

        gene_values = adj_row[:, gene_indices].toarray().ravel()
        if gene_values.size == 0:
            continue

        nz = np.flatnonzero(gene_values)
        if nz.size == 0:
            continue

        size = context["gene_adj_sizes"][gene_id]
        rows = nz // size
        cols = nz % size
        offset = context["feature_offsets"][gene_id]
        cell_edge_parts.append(np.vstack((rows + offset, cols + offset)))
        cell_weight_parts.append(gene_values[nz].reshape(-1, 1))

    if cell_edge_parts:
        cell_edge = np.concatenate(cell_edge_parts, axis=1)
        cell_edge_weight = np.concatenate(cell_weight_parts, axis=0)
    else:
        cell_edge = np.empty((2, 0), dtype=np.int64)
        cell_edge_weight = np.empty((0, 1), dtype=np.float32)

    return {
        "sample_name": sample_name,
        "label": label,
        "edge_index": cell_edge,
        "edge_attr": cell_edge_weight,
    }


def _build_cell_data(context, sample_idx):
    feature_row = _row_to_numpy(context["feature_matrix"], sample_idx).astype(np.float32, copy=False)
    feature_vec = feature_row.ravel()
    adj_raw_row = _row_to_numpy(context["adj_raw_matrix"], sample_idx).astype(np.float32, copy=False)
    edge_data = _build_cell_edge_arrays(context, sample_idx)

    return Data(
        x=torch.tensor(feature_vec.reshape(-1, 1), dtype=torch.float32),
        edge_index=torch.tensor(edge_data["edge_index"], dtype=torch.long),
        edge_attr=torch.tensor(edge_data["edge_attr"], dtype=torch.float32),
        y=torch.tensor(edge_data["label"]),
        x_fea=torch.tensor(feature_row, dtype=torch.float32),
        x_adj=torch.tensor(adj_raw_row, dtype=torch.float32),
        sample_name=edge_data["sample_name"],
    )


def write_graph_store(edge_store_path, context, pbar=None):
    total_samples = len(context["sample_list"])
    string_dtype = h5py.string_dtype(encoding="utf-8")

    with h5py.File(edge_store_path, "w") as h5:
        edge_ptr_ds = h5.create_dataset("edge_ptr", shape=(total_samples + 1,), dtype=np.int64)
        sample_name_ds = h5.create_dataset("sample_name", shape=(total_samples,), dtype=string_dtype)
        label_ds = h5.create_dataset("label", shape=(total_samples,), dtype=np.int32)
        edge_row_ds = h5.create_dataset(
            "edge_row",
            shape=(0,),
            maxshape=(None,),
            dtype=np.int32,
            chunks=True,
            compression="lzf",
        )
        edge_col_ds = h5.create_dataset(
            "edge_col",
            shape=(0,),
            maxshape=(None,),
            dtype=np.int32,
            chunks=True,
            compression="lzf",
        )
        edge_weight_ds = h5.create_dataset(
            "edge_weight",
            shape=(0,),
            maxshape=(None,),
            dtype=np.float32,
            chunks=True,
            compression="lzf",
        )

        current_edge_offset = 0
        edge_ptr_ds[0] = 0
        for sample_idx in range(total_samples):
            edge_data = _build_cell_edge_arrays(context, sample_idx)
            edge_index = edge_data["edge_index"]
            edge_attr = edge_data["edge_attr"].ravel()
            sample_name_ds[sample_idx] = edge_data["sample_name"]
            label_ds[sample_idx] = edge_data["label"]

            edge_count = int(edge_index.shape[1])
            if edge_count:
                next_edge_offset = current_edge_offset + edge_count
                edge_row_ds.resize((next_edge_offset,))
                edge_col_ds.resize((next_edge_offset,))
                edge_weight_ds.resize((next_edge_offset,))
                edge_row_ds[current_edge_offset:next_edge_offset] = edge_index[0].astype(np.int32, copy=False)
                edge_col_ds[current_edge_offset:next_edge_offset] = edge_index[1].astype(np.int32, copy=False)
                edge_weight_ds[current_edge_offset:next_edge_offset] = edge_attr.astype(np.float32, copy=False)
                current_edge_offset = next_edge_offset

            edge_ptr_ds[sample_idx + 1] = current_edge_offset
            if pbar is not None:
                pbar.update(1)

    return {
        "n_cells": total_samples,
        "n_edges": int(current_edge_offset),
    }


def get_graph_input(
    pbar,
    start_idx,
    sample_num,
    temp_folder,
    final_output_path,
    output_name,
    celltypename=None,
    mapper=None,
    context=None,
):
    if context is None:
        context = _build_graph_context(final_output_path, output_name, celltypename=celltypename, mapper=mapper)

    sample_list = context["sample_list"]
    end_idx = min(start_idx + sample_num, len(sample_list))

    cell_data_all = []
    for sample_idx in range(start_idx, end_idx):
        cell_data_all.append(_build_cell_data(context, sample_idx))
        pbar.update(1)

    torch.save(
        cell_data_all,
        os.path.join(temp_folder, "model_" + output_name + "_" + str(int(start_idx / sample_num)) + ".pt"),
    )

    del cell_data_all
    gc.collect()

    return pbar

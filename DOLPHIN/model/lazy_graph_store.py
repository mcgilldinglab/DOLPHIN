import json
import os

import anndata
import h5py
import numpy as np
import scipy.sparse as sp
import torch
from torch.utils.data import Dataset
from torch_geometric.data import Data


def _to_csr_matrix(matrix):
    if sp.issparse(matrix):
        return matrix.tocsr()
    return sp.csr_matrix(np.asarray(matrix))


def _row_to_numpy(matrix, row_idx):
    if sp.issparse(matrix):
        return matrix.getrow(row_idx).toarray()
    return np.asarray(matrix[row_idx : row_idx + 1])


def is_graph_store_input(path):
    return os.path.isdir(path) or path.endswith(".graph.json") or path.endswith("manifest.json")


def _resolve_manifest_path(path):
    if os.path.isdir(path):
        candidate_paths = [
            os.path.join(path, "manifest.json"),
            os.path.join(path, "model.graph.json"),
        ]
        for candidate_path in candidate_paths:
            if os.path.exists(candidate_path):
                return candidate_path
        raise FileNotFoundError(f"No graph-store manifest found under {path}")
    return path


class LazyGraphDataset(Dataset):
    def __init__(self, manifest_path):
        resolved_manifest = _resolve_manifest_path(manifest_path)
        with open(resolved_manifest, "r", encoding="utf-8") as handle:
            self.manifest = json.load(handle)

        self.feature_adata = anndata.read_h5ad(self.manifest["feature_h5ad"])
        self.adj_raw_adata = anndata.read_h5ad(self.manifest["adjacency_raw_h5ad"])
        self.feature_matrix = _to_csr_matrix(self.feature_adata.X)
        self.adj_raw_matrix = _to_csr_matrix(self.adj_raw_adata.X)

        self.edge_store = h5py.File(self.manifest["edge_store_h5"], "r")
        self.edge_ptr = self.edge_store["edge_ptr"]
        self.edge_row = self.edge_store["edge_row"]
        self.edge_col = self.edge_store["edge_col"]
        self.edge_weight = self.edge_store["edge_weight"]
        self.sample_name = self.edge_store["sample_name"]
        self.label = self.edge_store["label"]

    def __len__(self):
        return int(self.manifest["n_cells"])

    def __getitem__(self, idx):
        feature_row = _row_to_numpy(self.feature_matrix, idx).astype(np.float32, copy=False)
        feature_vec = feature_row.ravel()
        adj_raw_row = _row_to_numpy(self.adj_raw_matrix, idx).astype(np.float32, copy=False)

        start = int(self.edge_ptr[idx])
        end = int(self.edge_ptr[idx + 1])
        if end > start:
            edge_index = np.vstack(
                (
                    np.asarray(self.edge_row[start:end], dtype=np.int64),
                    np.asarray(self.edge_col[start:end], dtype=np.int64),
                )
            )
            edge_attr = np.asarray(self.edge_weight[start:end], dtype=np.float32).reshape(-1, 1)
        else:
            edge_index = np.empty((2, 0), dtype=np.int64)
            edge_attr = np.empty((0, 1), dtype=np.float32)

        sample_name = self.sample_name[idx]
        if isinstance(sample_name, bytes):
            sample_name = sample_name.decode("utf-8")

        return Data(
            x=torch.tensor(feature_vec.reshape(-1, 1), dtype=torch.float32),
            edge_index=torch.tensor(edge_index, dtype=torch.long),
            edge_attr=torch.tensor(edge_attr, dtype=torch.float32),
            y=torch.tensor(int(self.label[idx])),
            x_fea=torch.tensor(feature_row, dtype=torch.float32),
            x_adj=torch.tensor(adj_raw_row, dtype=torch.float32),
            sample_name=sample_name,
        )

    def close(self):
        if getattr(self, "edge_store", None) is not None:
            self.edge_store.close()
            self.edge_store = None

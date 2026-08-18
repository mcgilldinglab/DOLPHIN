import json
from .func_step03_GNN_main import _build_graph_context, write_graph_store
from tqdm import tqdm
import os
import pandas as pd

def run_model_input(
    metadata_path: str,
    out_name: str,
    out_directory: str = "./",
    gnn_run_num: int = 100,
    celltypename: str = None
    ):
    
    """
    Combines feature matrix and adjacency matrix and generates input for the DOLPHIN model.
    
    Parameters
    ----------
    metadata_path : str
        Path to the metadata file (e.g., a csv file with cell information).
    out_name : str
        Output filename for the feature matrix CSV.
    out_directory : str
        Output directory to save the combined feature matrix, default save to ./data/ folder.
    gnn_run_num : int
        Number of samples per GNN batch.
    celltypename : str, optional
        Column name in metadata indicating cell types. Default is None.

    Returns
    -------
    None
        Saves a graph-store manifest as `model_<out_name>.graph.json` in the output directory.

        The graph store keeps the canonical pooled h5ad matrices on disk and stores
        only edge_index / edge_attr / labels in a compact edge store. During model
        training, a lazy loader reconstructs the same PyTorch Geometric `Data`
        objects that the legacy `.pt` path produced.

    """
    
    df_label = pd.read_csv(metadata_path, sep='\t')
    total_sample_size = len(df_label)
    
    mapper = None
    if celltypename and celltypename in df_label.columns:
        unique_celltypes = sorted(df_label[celltypename].dropna().unique())
        mapper = {celltype: idx for idx, celltype in enumerate(unique_celltypes)}
    
    final_out_dir = os.path.join(out_directory, "data")
    os.makedirs(final_out_dir, exist_ok=True)
        
    print("Start Construct Data Input for model input")
    graph_context = _build_graph_context(
        final_out_dir,
        out_name,
        celltypename=celltypename,
        mapper=mapper,
    )
    total_sample_size = len(graph_context["sample_list"])
    edge_store_path = os.path.join(final_out_dir, "model_" + out_name + ".edges.h5")
    with tqdm(total=total_sample_size) as pbar_gnn:
        graph_store_summary = write_graph_store(edge_store_path, graph_context, pbar=pbar_gnn)

    manifest = {
        "format": "dolphin_graph_store_v1",
        "out_name": out_name,
        "n_cells": int(total_sample_size),
        "n_edges": int(graph_store_summary["n_edges"]),
        "feature_h5ad": graph_context["feature_h5ad_path"],
        "adjacency_raw_h5ad": graph_context["adj_raw_h5ad_path"],
        "adjacency_edge_h5ad": graph_context["adj_edge_h5ad_path"],
        "edge_store_h5": edge_store_path,
        "celltypename": celltypename,
        "gnn_run_num_legacy": int(gnn_run_num),
    }
    manifest_path = os.path.join(final_out_dir, "model_" + out_name + ".graph.json")
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    print(f"Saved graph-store manifest to {manifest_path}")

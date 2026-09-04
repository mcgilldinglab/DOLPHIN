import anndata
import os
import pandas as pd
import numpy as np
from scipy import sparse
from ._anndata_compat import enable_nullable_string_writes


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


def _normalize_within_gene(adata):
    matrix = _to_csr_matrix(adata.X)
    gene_ids = np.asarray(adata.var["gene_id"], dtype=object)
    gene_slices = _contiguous_slices(gene_ids)
    slice_gene_ids = [gene_id for gene_id, _, _ in gene_slices]
    if len(slice_gene_ids) != len(set(slice_gene_ids)):
        raise ValueError(
            "Adjacency columns for each gene must be contiguous before within-gene "
            "normalization."
        )

    normalized_blocks = []
    for _, start, stop in gene_slices:
        # Preserve the flattened adjacency order recorded in adata.var exactly.
        block = matrix[:, start:stop].astype(np.float64).tocsr(copy=True)
        gene_sum = np.asarray(block.sum(axis=1)).ravel()
        inv = np.divide(
            1.0,
            gene_sum,
            out=np.zeros_like(gene_sum, dtype=np.float64),
            where=gene_sum != 0,
        )
        normalized_blocks.append(block.multiply(inv[:, None]))

    if normalized_blocks:
        return sparse.hstack(normalized_blocks, format="csr").astype(np.float32)
    return sparse.csr_matrix(matrix.shape, dtype=np.float32)

def run_adjacency_hvg(
    out_name: str,
    out_directory="./"):
    """
    Processes the adjacency matrix by retaining only highly variable genes (HVGs)
    and performing within-gene normalization for graph construction in downstream models.
    
    Parameters
    ----------
    out_name : str
        Output filename for the feature matrix CSV.
    out_directory : str
        Output directory to save the combined feature matrix, default save to ./data/ folder.
    Returns
    -------
    None
        Saves two `.h5ad` files to the specified output directory:
        - `AdjacencyCompReHvg_<out_name>.h5ad`: HVG-filtered adjacency matrix.
        - `AdjacencyCompReHvgEdge_<out_name>.h5ad`: HVG-filtered and within-gene normalized adjacency matrix.
    """
    
    final_out_dir = os.path.join(out_directory, "data")
    os.makedirs(final_out_dir, exist_ok=True)
    
    hvg_path = os.path.join(final_out_dir, "ExonGene_hvg_"+out_name+".h5ad")

    ## load h5ad file
    adj_anndata = anndata.read_h5ad(os.path.join(final_out_dir, "AdjacencyCompRe_"+out_name+".h5ad"))

    hvg_adata = anndata.read_h5ad(hvg_path)
    hvg_list = set(hvg_adata.var["Geneid"])
    #only keep highly varaible genes
    adj_anndata = adj_anndata[:,adj_anndata.var["gene_id"].isin(hvg_list)].copy()

    print("Keep "+ str(len(set(adj_anndata.var["gene_id"]))) + " genes")
    print("The Final Adjacency Matrix Size is " + str(adj_anndata.shape[0]) + " Cells with dimension of " +str(adj_anndata.shape[1]))

    enable_nullable_string_writes()
    adj_anndata.write(os.path.join(final_out_dir, "AdjacencyCompReHvg_"+ out_name+".h5ad"))
    
    ## Normalization
    X_sparse = _normalize_within_gene(adj_anndata).astype(np.float32)

    # Create AnnData with sparse matrix
    new_adata = anndata.AnnData(X=X_sparse, obs=pd.DataFrame(adj_anndata.obs), var=pd.DataFrame(adj_anndata.var))

    # Save to h5ad
    results_file = os.path.join(final_out_dir, "AdjacencyCompReHvgEdge_" + out_name + ".h5ad")
    enable_nullable_string_writes()
    new_adata.write(results_file)

    # #conver to h5ad file
    # # # # ##the data matrix 
    # X = df_adj_norm.transpose().iloc[:,:].values
    # new_adata = anndata.AnnData(X, obs=pd.DataFrame(adj_anndata.obs), var=pd.DataFrame(adj_anndata.var))

    # results_file = os.path.join(final_out_dir, "AdjacencyCompReHvgEdge_"+ out_name+".h5ad")

    # new_adata.write(results_file)

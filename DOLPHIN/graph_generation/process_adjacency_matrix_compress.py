from multiprocessing import Process
import numpy as np
from itertools import combinations, tee
import pandas as pd
import anndata
import os
from scipy import sparse
import warnings
from anndata._core.views import ImplicitModificationWarning

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=ImplicitModificationWarning)

def run_adjacency_compression(
    metadata_path,
    out_name,
    out_directory,
    num_processes=25
):
    """
    Compute and update exon-level adjacency matrices per gene for each cell.

    This function loads gene feature and adjacency data, selects genes that have 
    more than one exon and non-zero expression, and reconstructs the adjacency 
    matrices for each selected gene based on the positions of exons with non-zero values.

    If both exons in an adjacency pair have zero expression, the corresponding 
    adjacency value will be set to 0.

    The updated adjacency matrix is saved for each cell in `.h5ad` format.

    Parameters
    ----------
    out_name : str
        Prefix for loading input H5AD files, e.g. "LUAD".
        Expects files named "Adjacency_<out_name>.h5ad" and "Feature_<out_name>.h5ad".
    metadata_path : str
        Path to a metadata file containing cell IDs under column "CB".
    out_directory : str
        Root directory for reading data and saving outputs. Final outputs will go to:
        <out_directory>/data/temp/adj_comp_matrix/
    num_processes : int, optional
        Number of parallel processes to run. Default is 25.

    Notes
    -----
    - Requires `anndata`, `numpy`, `pandas`, and `multiprocessing`.
    - Only genes with more than one exon AND non-zero expression in at least
      one cell are considered.
    - For each cell, the updated adjacency matrix is saved as:
      <cell>.h5ad in the output directory.

    """

    def pairwise(iterable):
        a, b = tee(iterable)
        next(b, None)
        return zip(a, b)

    # Load input AnnData files
    adata_adj_orig = anndata.read_h5ad(os.path.join(out_directory, "data", f"Adjacency_{out_name}.h5ad"))
    adata_adj_up = anndata.AnnData.copy(adata_adj_orig)
    adata_adj_up.X = adata_adj_orig.X.toarray()

    adata_fea_orig = anndata.read_h5ad(os.path.join(out_directory, "data", f"Feature_{out_name}.h5ad"))
    adata_fea_orig.X = adata_fea_orig.X.toarray()

    # Determine valid gene list
    df_fea = adata_fea_orig[:, :].to_df().T
    df_fea['gene_id'] = adata_fea_orig[0, :].var["gene_id"]
    df_fea_sum = df_fea.groupby("gene_id", observed=False).sum().reset_index().set_index("gene_id")
    nonzero_gene_list = df_fea_sum[~(df_fea_sum == 0).all(axis=1)].index.tolist()

    exon_count = pd.DataFrame(adata_fea_orig.var).reset_index().groupby("gene_id", observed=False).last().reset_index()
    exon_count["Exon_Number"] = exon_count["index"].apply(lambda x: int(x.split('-')[-1]))
    exon_gene_list = exon_count[exon_count["Exon_Number"] > 1]["gene_id"].tolist()

    gene_list = list(set(nonzero_gene_list) & set(exon_gene_list))

    pd_gt = pd.read_csv(metadata_path, sep='\t')
    sample_list = list(pd_gt["CB"])

    # Create output directory
    out_path_dir = os.path.join(out_directory, "data", "temp", "adj_comp_matrix")
    os.makedirs(out_path_dir, exist_ok=True)

    def process_one_sample(j):
        print(f"Starting processing cell: {sample_list[j]}")  

        for i in range(len(gene_list)):
            # print(f"Processing: Sample {sample_list[j]} - Gene {gene_list[i]}")
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]

            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape((_size, _size))
            _new_adj = np.zeros((_size, _size))

            _adj_pair = list(pairwise(non_zero_index))

            for v in combinations(non_zero_index, 2):
                if (_old_adj[v] == 0) and (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]

            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()

        # out_file = os.path.join(out_path_dir, f"{sample_list[j]}.h5ad")
        # adata_adj_up[sample_list[j], :].X = sparse.csr_matrix(adata_adj_up[sample_list[j], :].X)
        # adata_adj_up[sample_list[j], :].write_h5ad(out_file)
        
        out_file = os.path.join(out_path_dir, f"{sample_list[j]}.h5ad")
        adata_to_write = adata_adj_up[sample_list[j], :].copy()
        adata_to_write.X = sparse.csr_matrix(adata_to_write.X)
        adata_to_write.write_h5ad(out_file)

        print(f"Finished processing cell: {sample_list[j]}")

    def process_range(start_idx, end_idx):
        for j in range(start_idx, min(end_idx, len(sample_list))):
            process_one_sample(j)

    # Launch multiprocessing
    processes = []
    start = 0
    step = int(np.ceil(len(sample_list) / num_processes))

    for p in range(num_processes):
        start_idx = start + p * step
        end_idx = start + (p + 1) * step
        proc = Process(target=process_range, args=(start_idx, end_idx))
        processes.append(proc)
        proc.start()

    for proc in processes:
        proc.join()
        
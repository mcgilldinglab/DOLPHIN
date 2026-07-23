import os
import pandas as pd
import numpy as np
import anndata
from tqdm import tqdm
import sys

from .grouped_featurecounts import load_grouped_featurecounts_table
from ._anndata_compat import enable_nullable_string_writes

def run_raw_gene(
    metadata_path: str,
    gtf_path: str,
    out_name: str,
    featurecount_path: str = None,
    n_hvg = 2000,
    out_directory="./",
    grouped_featurecount_path: str = None):
    
    """
    Combines featureCounts results into a gene count matrix and selects highly variable genes (HVGs).

    This function reads featureCounts gene-level count files and sample metadata,
    constructs a combined gene count matrix, and identifies the top highly variable genes
    for downstream analysis.
    
    Parameters
    ----------
    metadata_path : str
        Path to the metadata file (e.g., a csv file with cell information).
    featurecount_path : str
        Path to the directory containing gene-level featureCounts output files.
        Kept for backward compatibility with the one-file-per-cell layout.
    gtf_path : str
        Path to the GTF file used for gene annotation.
    out_name : str
        Output filename for the feature matrix CSV.
    out_directory : str
        Output directory to save the combined feature matrix, default save to ./data/ folder.
    n_hvg: int
        Number of highly variable genes to select. Defaults to 2000.
    grouped_featurecount_path : str, optional
        Optional grouped gene count matrix produced by featureCounts. If provided,
        this becomes the preferred input instead of `featurecount_path`.

    Returns
    -------
    None
        Saves the final annotated AnnData object as `ExonGene_<out_name>.h5ad` in the specified output directory.
    """
    
    final_out_dir = os.path.join(out_directory, "data")
    os.makedirs(final_out_dir, exist_ok=True)
    
    pd_gt = pd.read_csv(metadata_path, sep='\t')
    cell_ids = list(pd_gt["CB"])
    pd_gt = pd_gt.set_index("CB", drop=False)
    pd_gt.index.name = None
    if grouped_featurecount_path is not None:
        if featurecount_path is not None:
            raise ValueError(
                "Use either grouped_featurecount_path or featurecount_path, not both."
            )
        pd_count = load_grouped_featurecounts_table(
            grouped_featurecount_path,
            mode="gene",
            requested_cell_ids=cell_ids,
        )
        missing_cells = [cell_id for cell_id in cell_ids if cell_id not in pd_count.columns]
        if missing_cells:
            sys.exit(
                "There is a mismatch between the metadata and grouped FeatureCounts "
                f"results. Missing cells include: {missing_cells[:10]}"
            )
        pd_count = pd_count[["Geneid"] + cell_ids].copy()
    else:
        if featurecount_path is None:
            raise ValueError(
                "featurecount_path is required when grouped_featurecount_path is not provided."
            )

        cnt_files = []
        for f in os.listdir(featurecount_path):
            if f.endswith("count.txt"):
                if f.split(".")[0] in cell_ids:
                    cnt_files.append(f)
        if len(cnt_files) != pd_gt.shape[0]:
            sys.exit("There is a mismatch between the metadata and FeatureCounts results. Please check.")
            
        pd_count = pd.DataFrame([])
        for i, f in enumerate(tqdm(cnt_files)):
            _cb = f.split(".")[0]
            pd_cb = pd.read_csv(os.path.join(featurecount_path,f), sep="\t", skiprows=1)
            pd_cb.columns = [*pd_cb.columns[:-1], _cb]
            pd_cb = pd_cb[["Geneid", _cb]]
            if i == 0:
                pd_count = pd_cb
            else:
                pd_count= pd.merge(pd_count, pd_cb, left_on=["Geneid"], right_on=["Geneid"], how='outer')

    ####complete function, may need slight modification based on your gtf format####
    def get_ens_dict(file_path):
        gtf_dict = {}
        found_annotation = False
        with open(file_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                if 'gene_id "' not in line or 'gene_name "' not in line:
                    continue
                found_annotation = True
                gene_id = line.split('gene_id "')[1].split('"')[0]
                gene_name = line.split('gene_name "')[1].split('"')[0]
                gtf_dict[gene_id] = gene_name
        if not found_annotation:
            print('you need to change gene_id " and gene_name " formats')
        return gtf_dict

    gtf_dict = get_ens_dict(gtf_path)

    pd_count = pd_count.set_index("Geneid", drop=False)
    pd_count.index.name = None
    gene_name = pd_count["Geneid"].map(gtf_dict).fillna(pd_count["Geneid"])
    pd_count["GeneName"] = gene_name
    pd_count.set_index("GeneName", drop=False, inplace=True)
    pd_count.index.name = None

    #conver to h5ad file
    ## dataframe for annotating the observations = sample name
    obs = pd_gt.loc[cell_ids, pd_gt.columns].copy()
    obs.index.name = None

    ## dataframe for annotating the variables = geneid
    var = pd_count[["Geneid", "GeneName"]]

    # # # ##the data matrix 
    X = pd_count.loc[:, cell_ids].T.to_numpy(dtype=np.float32, copy=True)
    adata = anndata.AnnData(X, obs=obs, var=var, dtype=np.float32)
    
    enable_nullable_string_writes()
    adata.write(os.path.join(final_out_dir, "ExonGene_"+out_name+".h5ad"))
    
    ## get hvg
    import scanpy as sc
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    enable_nullable_string_writes()
    adata.write(os.path.join(final_out_dir, "ExonGene_hvg_"+out_name+".h5ad"))
    

import os
import pandas as pd
import numpy as np
import anndata
from tqdm import tqdm
import sys
import scanpy as sc

def run_raw_gene(
    metadata_path: str,
    featurecount_path: str,
    gtf_path: str,
    out_name: str,
    n_hvg = 2000,
    out_directory="./"):
    
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
    gtf_path : str
        Path to the GTF file used for gene annotation.
    out_name : str
        Output filename for the feature matrix CSV.
    out_directory : str
        Output directory to save the combined feature matrix, default save to ./data/ folder.
    n_hvg: int
        Number of highly variable genes to select. Defaults to 2000.

    Returns
    -------
    None
        Saves the final annotated AnnData object as `ExonGene_<out_name>.h5ad` in the specified output directory.
    """
    
    final_out_dir = os.path.join(out_directory, "data")
    os.makedirs(final_out_dir, exist_ok=True)
    
    pd_gt = pd.read_csv(metadata_path, sep='\t')

    cnt_files = []
    for f in os.listdir(featurecount_path):
        if f.endswith("count.txt"):
            if f.split(".")[0] in list(pd_gt["CB"]):
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

    pd_count = pd_count.set_index("Geneid", drop=False)
    pd_count.index.name=None

    pd_count_t = pd_count.drop('Geneid', axis=1).T
    pd_count_t = pd.merge(pd_gt, pd_count_t, left_on="CB",right_index=True)
    pd_count_t.set_index('CB', drop=False, inplace=True)
    pd_count_t.index.name = None

    ####complete function, may need slight modification based on your gtf format####
    def get_ens_dict(file_path):
        with open(file_path) as f:
            gtf = list(f)

        gtf = [x for x in gtf if not x.startswith('#')]
        gtf = [x for x in gtf if 'gene_id "' in x and 'gene_name "' in x]
        if len(gtf) == 0:
            print('you need to change gene_id " and gene_name " formats')
        
        gtf = list(map(lambda x: (x.split('gene_id "')[1].split('"')[0], x.split('gene_name "')[1].split('"')[0]), gtf))
        gtf = dict(set(gtf))
        return gtf

    gtf_dict = get_ens_dict(gtf_path)

    pd_count["GeneName"] = pd_count["Geneid"]
    pd_count = pd_count.replace({"GeneName": gtf_dict})

    pd_count.set_index("GeneName", drop=False, inplace=True)
    pd_count.index.name=None

    #conver to h5ad file
    ## dataframe for annotating the observations = sample name
    obs = pd_count_t[pd_gt.columns]

    ## dataframe for annotating the variables = geneid
    var = pd_count[["Geneid", "GeneName"]]

    # # # ##the data matrix 
    X = pd_count_t.iloc[:,pd_gt.shape[1]:].values
    adata = anndata.AnnData(X, obs=obs, var=var, dtype=np.float32)
    
    adata.write(os.path.join(final_out_dir, "ExonGene_"+out_name+".h5ad"))
    
    ## get hvg
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    adata.write(os.path.join(final_out_dir, "ExonGene_hvg_"+out_name+".h5ad"))
    
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import math
import os

def get_fea_hvg(output_path, output_name):
    hvg_path = os.path.join(output_path, "ExonGene_hvg_"+output_name+".h5ad")
    feature_anndata = os.path.join(output_path, "FeatureComp_"+output_name+".h5ad")

    ## load feature matrix
    adata = anndata.read_h5ad(feature_anndata)

    hvg_adata = anndata.read_h5ad(hvg_path)
    cell_keep = list(hvg_adata.obs.index)
    hvg_list = set(hvg_adata.var["Geneid"])

    adata = adata[adata.obs.index.isin(cell_keep), :]

    #based on the total counts plot, I will normalized to median value
    sc.pp.normalize_total(adata)

    #only keep highly varaible genes = 5000 genes, where highly_variable == True
    adata = adata[:,adata.var["gene_id"].isin(hvg_list)]

    print("Keep "+ str(len(set(adata.var["gene_id"]))) + " genes")
    print("The Final Feature Matrix Size is " + str(adata.shape[0]) + " Cells and " +str(adata.shape[1])+ " exons")

    #remove unneccary var names
    adata.write(os.path.join(output_path, "FeatureCompHvg_"+output_name+".h5ad"))
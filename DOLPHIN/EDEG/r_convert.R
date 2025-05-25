#!/usr/bin/env Rscript

library(dplyr)
library(Seurat)
library(patchwork)
library(reticulate)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]
output_path <- args[2]

# Step 1: Automatically bind to the current Conda Python environment
try({
  conda_env <- Sys.getenv("CONDA_DEFAULT_ENV")
  use_condaenv(conda_env, required = TRUE)
}, silent = TRUE)

# Step 2: Import Python 'scanpy' module
scanpy <- import("scanpy")

# Step 3: Convert h5ad file to Seurat object
data_preprocess <- function(data_file) {
  adata <- scanpy$read(data_file)
  meta <- adata$obs
  gene <- adata$var
  adata2 <- as.matrix(adata$X)
  adata_raw <- t(adata2)
  rownames(adata_raw) <- rownames(gene)
  colnames(adata_raw) <- rownames(meta)
  merge <- CreateSeuratObject(adata_raw)
  merge <- AddMetaData(merge, meta)
  return(merge)
}

# Run the conversion
seurat_obj <- data_preprocess(data_file)

# Save the Seurat object to output path
saveRDS(seurat_obj, file = output_path)

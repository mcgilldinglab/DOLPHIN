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
  meta_dict <- py_to_r(adata$obs$to_dict("list"))
  gene_dict <- py_to_r(adata$var$to_dict("list"))
  adata2 <- as.matrix(adata$X)
  adata_raw <- t(adata2)
  gene_names <- tryCatch(py_to_r(adata$var_names$to_list()), error = function(e) NULL)
  cell_names <- tryCatch(py_to_r(adata$obs_names$to_list()), error = function(e) NULL)
  if (is.null(gene_names) || length(gene_names) != nrow(adata_raw)) {
    gene_names <- names(gene_dict[[1]])
  }
  if (is.null(cell_names) || length(cell_names) != ncol(adata_raw)) {
    stop("Failed to recover cell names from AnnData obs_names.")
  }
  rownames(adata_raw) <- gene_names
  colnames(adata_raw) <- cell_names
  meta <- as.data.frame(meta_dict, stringsAsFactors = FALSE)
  rownames(meta) <- cell_names
  meta[] <- lapply(meta, function(x) {
    if (is.factor(x)) {
      as.character(x)
    } else if (is.list(x)) {
      vapply(x, as.character, character(1))
    } else {
      as.vector(x)
    }
  })
  merge <- CreateSeuratObject(adata_raw)
  merge <- AddMetaData(merge, meta)
  return(merge)
}

# Run the conversion
seurat_obj <- data_preprocess(data_file)

# Save the Seurat object to output path
saveRDS(seurat_obj, file = output_path)

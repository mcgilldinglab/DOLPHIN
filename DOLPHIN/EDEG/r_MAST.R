#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)

input_rds <- args[[1]]
output_dir <- args[[2]]
ident_column <- ifelse(length(args) >= 3, args[[3]], "cluster")
cluster_values_raw <- ifelse(length(args) >= 4, args[[4]], "")
logfc_threshold <- ifelse(length(args) >= 5, as.numeric(args[[5]]), 0.0)
normalize_flag <- ifelse(length(args) >= 6, toupper(args[[6]]), "TRUE")

cluster_values <- character(0)
if (nzchar(cluster_values_raw)) {
  cluster_values <- strsplit(cluster_values_raw, ",", fixed = TRUE)[[1]]
  cluster_values <- trimws(cluster_values)
  cluster_values <- cluster_values[nzchar(cluster_values)]
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seurat_obj <- readRDS(file = input_rds)

if (normalize_flag == "TRUE") {
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
}

if (!(ident_column %in% colnames(seurat_obj@meta.data))) {
  stop(sprintf("ident_column '%s' missing from Seurat metadata.", ident_column))
}

Idents(seurat_obj) <- ident_column

if (length(cluster_values) == 0) {
  cluster_values <- levels(Idents(seurat_obj))
}

for (cluster_name in cluster_values) {
  de_mast <- FindMarkers(
    seurat_obj,
    ident.1 = cluster_name,
    test.use = "MAST",
    logfc.threshold = logfc_threshold
  )
  output_csv <- file.path(output_dir, paste0(cluster_name, "_MAST_DE.csv"))
  write.csv(de_mast, file = output_csv, row.names = TRUE)
}

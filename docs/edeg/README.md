# EDEG / JDEG

This section covers exon-level and junction-level differential gene analysis.

## Preserved Route

The retained public route is the historical one-cluster-vs-rest workflow:

- convert `h5ad` to Seurat `rds`
- run Seurat `FindMarkers(..., test.use = "MAST")` for one cluster against the rest
- write one `*_MAST_DE.csv` per cluster
- aggregate exon-level markers to gene-level EDEG
- aggregate junction-level markers to gene-level JDEG

Pairwise MAST is not preserved as a public route in this clean package.

## Input Contract

EDEG historical input:

- use `FeatureComp_<out_name>.h5ad`
- convert to `FeatureComp_<out_name>.rds`
- run historical one-vs-rest MAST
- feed the resulting `<cluster>_MAST_DE.csv` into `run_edeg(...)`

JDEG historical input:

- use the junction-analysis matrix prepared for Seurat/MAST
- run historical one-vs-rest MAST
- feed the resulting `<cluster>_MAST_DE.csv` into `run_jdeg(...)`

## What The Statistics Mean

- MAST first identifies differential exon-level or junction-level features
- EDEG then aggregates exon-level evidence across all exons belonging to the same gene
- JDEG then aggregates junction-level evidence across all junction features belonging to the same gene

## Preserved Components

- explicit Python wrapper for the R conversion step
- explicit Python wrapper for the historical Seurat/MAST step
- reusable exon-level aggregation function
- reusable junction-level aggregation function
- small CLI for repeatable runs

## Package Location

- [DOLPHIN/EDEG](../../DOLPHIN/EDEG)

## Important Note

The original repository did not contain a standalone `r_MAST.R` file.
The tutorial and historical run records performed the MAST step either inside an R notebook or through ad hoc cluster-vs-rest R scripts.

This clean package restores a formal `r_MAST.R` that matches the historical one-cluster-vs-rest workflow.

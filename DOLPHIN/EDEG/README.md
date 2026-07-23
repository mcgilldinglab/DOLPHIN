# EDEG / JDEG

This module preserves two layers:

1. tutorial-compatible workflow pieces
   - `call_convert.py`
   - `call_MAST.py`
   - `r_convert.R`
   - `r_MAST.R`
   - `generate_EDEG.py`
   - `generate_JDEG.py`

2. cleaner reusable analysis function
   - `edeg_analysis.py`

The preserved public route is the historical one-cluster-vs-rest workflow.

## Historical MAST Contract

1. convert `.h5ad` to Seurat `.rds`
2. run Seurat `FindMarkers(..., test.use = "MAST")` in one-cluster-vs-rest mode
3. write one marker CSV per cluster
4. aggregate exon-level markers to gene-level EDEG
5. aggregate junction-level markers to gene-level JDEG

This module does not expose pairwise MAST as a public workflow.

## Input / Output Contract

EDEG historical input:

- exon matrix input is `FeatureComp_<out_name>.h5ad`
- converted Seurat object is `FeatureComp_<out_name>.rds`
- MAST output is `<cluster>_MAST_DE.csv`
- final exon-to-gene aggregation output is `<cluster>_EDEG.csv`

JDEG historical input:

- junction matrix input is the junction-analysis matrix used for Seurat/MAST
- MAST output is `<cluster>_MAST_DE.csv`
- final junction-to-gene aggregation output is `<cluster>_JDEG.csv`

## Statistical Meaning

The sequence is:

- MAST first identifies differential exons or differential junction features
- `generate_EDEG.py` then combines exon-level evidence for all exons under the same gene
- `generate_JDEG.py` then combines junction-level evidence for all junction features under the same gene

The aggregation logic comes from the original DOLPHIN EDEG/JDEG code. The missing piece that was restored here is the standalone MAST runner.

# Pipeline Outputs

This document defines the retained outputs for the finalized DOLPHIN pipeline.

## Step 0: Reference Preparation

Input:
- source GTF
- optional genome FASTA

Main outputs:
- `reference_manifest.json`
- `dolphin.exon.gtf`
- `dolphin.exon.pkl`
- `dolphin_adj_index.csv`
- `dolphin_adj_metadata_table.csv`
- `dolphin_gene_meta.csv`
- optional `star_index/`

Purpose:
- the exon/junction reference bundle used by all later preprocessing and graph-generation steps

## Step 1: Raw Preprocessing

### Full-length

Main retained outputs:
- `cell_manifest.tsv`
- grouped gene count table
- grouped exon count table
- grouped junction count table

Compatibility outputs:
- per-cell gene counts
- per-cell exon/junction counts

Purpose:
- convert raw FASTQ into stable grouped count matrices and a cell manifest

### 10x

Main retained outputs:
- `cell_manifest.tsv`
- grouped gene count table
- grouped exon count table
- grouped junction count table

Compatibility outputs:
- split per-cell gene counts
- split per-cell exon/junction counts
- RG-tagged BAM when requested by the pipeline

Purpose:
- convert STARsolo-aligned 10x data into grouped counts and per-cell split counts

## Step 2: Graph Generation

Main retained outputs:
- `Feature_<out_name>.h5ad`
- `FeatureComp_<out_name>.h5ad`
- `Adjacency_<out_name>.h5ad`
- `AdjacencyComp_<out_name>.h5ad`
- `AdjacencyCompRe_<out_name>.h5ad`
- `ExonGene_<out_name>.h5ad`
- `ExonGene_hvg_<out_name>.h5ad`
- `FeatureCompHvg_<out_name>.h5ad`
- `AdjacencyCompReHvg_<out_name>.h5ad`
- `AdjacencyCompReHvgEdge_<out_name>.h5ad`

Not retained in the final efficient path:
- per-cell `*_fea.csv`
- per-cell `*_adj.csv`
- per-cell graph H5AD fragments as user-facing outputs

Purpose:
- transform grouped or split featureCounts outputs into pooled graph matrices and model-ready HVG matrices

## Step 3: Model Input

Main retained output:
- `model_<out_name>.graph.json`
- `model_<out_name>.edges.h5`

Legacy compatibility output:
- `model_<out_name>.pt` only if explicitly generated for backward compatibility

Purpose:
- lightweight graph-store input for DOLPHIN model training

## Step 4: Model Training

Input:
- `FeatureCompHvg_<out_name>.h5ad`
- `model_<out_name>.graph.json`
- `model_<out_name>.edges.h5`

Main retained output:
- `DOLPHIN_Z.h5ad`

Purpose:
- final latent embedding for downstream clustering and biological analysis

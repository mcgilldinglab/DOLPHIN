# DOLPHIN Preprocessing Status

All paths below are relative to `/mnt/md0/kailu/DOLPHIN_codex`.

## Scope

This document records the current stable status of DOLPHIN preprocessing after optimization.

Here, "raw preprocessing" means:

1. Reference preparation (`step0`)
2. Raw FASTQ to grouped count matrices (`step1`)

The canonical downstream handoff from preprocessing is:

- one `cell_manifest.tsv`
- one grouped gene count table
- one grouped exon count table
- one grouped junction count table

The compatibility one-cell-per-file outputs are still available, but they are no longer the recommended main input for downstream graph generation.

## Current Status

Raw preprocessing is complete for both data types currently tested:

- `full-length / flashseq`: complete
- `10x / SRR8513796`: complete

The stable full-length raw-preprocessing run is:

- `execution_runs/full_length_flashseq_clean_v1`

The stable 10x raw-preprocessing run is:

- `execution_runs/tenx_srr8513796_clean_v2`

## Step0: Reference Preparation

Stable reference bundle:

- `reference_inputs/grch38_release107/dolphin_exon_gtf`

Main outputs:

- `dolphin.exon.gtf`
- `dolphin_exon_gtf.gtf`
- `dolphin.exon.pkl`
- `dolphin_adj_index.csv`
- `dolphin_adj_metadata_table.csv`
- `dolphin_gene_meta.csv`
- `reference_manifest.json`
- `step_timing.tsv`
- `star_index/`

Timing summary from `reference_inputs/grch38_release107/dolphin_exon_gtf/step_timing.tsv`:

- `generate_nonoverlapping_exons.total`: `103.104687 s` (`1.72 min`)
- `prepare_reference_bundle.build_star_genome_index`: `708.971279 s` (`11.82 min`)
- `prepare_reference_bundle.total`: `822.451221 s` (`13.71 min`)

Interpretation:

- If only DOLPHIN exon reference files are needed, the heavy Python reference-generation part is about `1.7 min`.
- If STAR genome index is also built, the full step0 wall time is about `13.7 min`.

## Step1: Full-Length Raw Preprocessing

Stable run:

- `execution_runs/full_length_flashseq_clean_v1`

Inputs used:

- raw FASTQ directory for all cells
- prepared reference bundle from step0

Stable canonical outputs:

- `execution_runs/full_length_flashseq_clean_v1/cell_manifest.tsv`
- `execution_runs/full_length_flashseq_clean_v1/cell_ids.txt`
- `execution_runs/full_length_flashseq_clean_v1/outputs/04_exon_gene_cnt_grouped/full_length.exongene.count.txt`
- `execution_runs/full_length_flashseq_clean_v1/outputs/05_exon_junct_cnt_grouped/full_length.exon.count.txt`
- `execution_runs/full_length_flashseq_clean_v1/outputs/05_exon_junct_cnt_grouped/full_length.exon.count.txt.jcounts`

Compatibility outputs still preserved:

- `execution_runs/full_length_flashseq_clean_v1/outputs/04_exon_gene_cnt/`
- `execution_runs/full_length_flashseq_clean_v1/outputs/05_exon_junct_cnt/`

Timing summary from `execution_runs/full_length_flashseq_clean_v1/logs/execution_timing.tsv`:

- `align_exon` total: `23062.537399 s` (`384.38 min`, `6.41 h`)
- `count_gene` total: `552.942609 s` (`9.22 min`)
- `count_exon_junction` total: `1515.656326 s` (`25.26 min`)
- total retained step1 time: `25131.136334 s` (`418.85 min`, `6.98 h`)

Average per-cell timing over `795` cells:

- alignment: about `29.01 s / cell`
- gene count: about `0.70 s / cell`
- exon + junction count: about `1.91 s / cell`

## Step1: 10x Raw Preprocessing

Stable run:

- `execution_runs/tenx_srr8513796_clean_v2`

Inputs used:

- one 10x raw FASTQ sample
- prepared reference bundle from step0
- STARsolo barcode handling

Stable canonical outputs:

- `execution_runs/tenx_srr8513796_clean_v2/cell_manifest.tsv`
- `execution_runs/tenx_srr8513796_clean_v2/cell_ids.txt`
- `execution_runs/tenx_srr8513796_clean_v2/outputs/SRR8513796/05_exon_gene_cnt_grouped/SRR8513796.exongene.count.txt`
- `execution_runs/tenx_srr8513796_clean_v2/outputs/SRR8513796/06_exon_junct_cnt_grouped/SRR8513796.exon.count.txt`
- `execution_runs/tenx_srr8513796_clean_v2/outputs/SRR8513796/06_exon_junct_cnt_grouped/SRR8513796.exon.count.txt.jcounts`

Compatibility outputs still preserved:

- `execution_runs/tenx_srr8513796_clean_v2/outputs/SRR8513796/05_exon_gene_cnt/`
- `execution_runs/tenx_srr8513796_clean_v2/outputs/SRR8513796/06_exon_junct_cnt/`

Retained intermediate needed by the current stable 10x pipeline:

- `execution_runs/tenx_srr8513796_clean_v2/outputs/SRR8513796/04_rg_tagged_bam/SRR8513796.cb_rg.bam`

Timing summary from `execution_runs/tenx_srr8513796_clean_v2/outputs/SRR8513796/logs/execution_timing.tsv`:

- `cb_to_rg`: `2107.146539 s` (`35.12 min`)
- `grouped_gene_count`: `61.779264 s` (`1.03 min`)
- `grouped_exon_junction_count`: `511.111892 s` (`8.52 min`)
- `split_gene_count`: `226.320669 s` (`3.77 min`)
- `split_exon_count`: `498.388960 s` (`8.31 min`)
- `split_junction_count`: `217.108160 s` (`3.62 min`)
- total retained post-alignment time: `3621.855484 s` (`60.36 min`)

Important note:

- the retained `v2` run reused an already completed STARsolo alignment, so the current timing file does not preserve a standalone full alignment wall time
- the stable grouped outputs themselves are valid and fully generated

## Canonical Preprocessing Output Contract

After raw preprocessing, downstream graph generation should consume:

- `cell_manifest.tsv`
- grouped gene count table
- grouped exon count table
- grouped junction count table

Current canonical format by data type:

### Full-Length

- `cell_manifest.tsv`
- `full_length.exongene.count.txt`
- `full_length.exon.count.txt`
- `full_length.exon.count.txt.jcounts`

### 10x

- `cell_manifest.tsv`
- `SRR8513796.exongene.count.txt`
- `SRR8513796.exon.count.txt`
- `SRR8513796.exon.count.txt.jcounts`

## What Is Considered Stable

The following are considered stable, validated preprocessing baselines:

- step0 reference bundle generation
- full-length raw FASTQ to grouped counts
- 10x raw FASTQ to grouped counts
- grouped count files as the main downstream handoff

The following are preserved only for compatibility/debugging and are not the preferred primary handoff:

- one-cell-per-file featureCounts outputs
- 10x split compatibility count directories

## Downstream Handoff

The current stable downstream handoff for graph generation is:

- step0 reference bundle from `reference_inputs/grch38_release107/dolphin_exon_gtf`
- grouped raw-count outputs from one of the stable step1 runs

The downstream graph-generation baseline is additionally documented in:

- `DOLPHIN_optimized/GRAPH_GENERATION_BASELINE.md`

## Additional Full-Length Graph Matrices Already Generated

Although this document focuses on raw preprocessing, the following full-length pooled graph matrices have already been generated and retained:

- `execution_runs/graph_step2_full_length_allcells_opt_v1/data/Feature_full_length_allcells.h5ad`
- `execution_runs/graph_step2_full_length_allcells_opt_v1/data/FeatureComp_full_length_allcells.h5ad`
- `execution_runs/graph_step2_full_length_allcells_opt_v1/data/Adjacency_full_length_allcells.h5ad`
- `execution_runs/graph_step45_full_length_allcells_opt_v2/data/AdjacencyComp_full_length_allcells.h5ad`

These belong to graph generation rather than raw preprocessing, but they confirm that the current full-length preprocessing outputs are usable downstream.

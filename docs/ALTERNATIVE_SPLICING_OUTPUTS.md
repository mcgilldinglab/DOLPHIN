# Alternative Splicing Outputs

This document defines the retained outputs for the finalized DOLPHIN alternative splicing module.

## Shared Outputs

Each AS run writes:

- `metadata_prepared.tsv`
- `metadata_preparation_summary.json`
- `N_<out_name>_<k>.csv`
- `<out_name>_read_counts.csv`
- `aggregated_star/`
- `outrigger_output/`
- `alternative_splicing/<out_name>_PSI.h5ad`
- `alternative_splicing/<out_name>_PSI_random.h5ad`
- `alternative_splicing/<out_name>_PSI_DAS.h5ad`
- `<logs_root>/<out_name>.status.json`
- `<logs_root>/<out_name>.summary.json`

## BAM Route Extra Outputs

The BAM route additionally retains:

- `cell_aggregation/`
  - per-cell aggregated BAM
- `aggregated_star/<cell>/<cell>.aggr.SJ.out.tab`
  - STAR realignment junction summaries

## Junction Route Extra Outputs

The junction route retains:

- `aggregated_star/<cell>/<cell>.aggr.SJ.out.tab`
  - direct aggregated junction summaries

It does not retain:

- aggregated BAM
- STAR realignment BAM outputs

## Modality Differences

full-length:

- consumes existing per-cell BAM and per-cell STAR `SJ.out.tab`

10x BAM:

- consumes pooled RG-tagged BAM
- prepares reusable per-cell BAM cache before aggregation

10x junction:

- consumes per-cell `.jcounts`
- backfills `annotated`, `motif`, and `strand` from reference/GTF before Outrigger

## Comparison Outputs

When BAM vs junction comparison is run, keep:

- `bam_vs_junction_comparison.json`

This file records:

- shared cells
- shared events
- Pearson correlation
- absolute difference summary
- shape differences between routes

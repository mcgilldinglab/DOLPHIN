# Alternative Splicing

This document describes the finalized DOLPHIN alternative splicing code organization.

## Package Location

- [DOLPHIN/alternative_splicing](../../DOLPHIN/alternative_splicing)

## Finalized Routes

Four validated routes are preserved:

1. `full-length` + `bam`
2. `full-length` + `junction`
3. `10x` + `bam`
4. `10x` + `junction`

## What Is Shared

- neighbor finding from DOLPHIN embeddings
- per-cell read counting
- Outrigger `index`, `validate`, and `psi`
- conversion to `PSI.h5ad`, `PSI_random.h5ad`, and `PSI_DAS.h5ad`
- event-level differential AS testing in `<out_name>_DAS.csv`
- comparison helper for BAM vs junction outputs

## What Differs

By route:

- BAM route performs:
  - neighbor-aware BAM aggregation
  - STAR realignment on aggregated BAM
- junction route performs:
  - neighbor-aware direct junction aggregation
  - no BAM aggregation
  - no STAR realignment

By modality:

- full-length starts from per-cell BAM + per-cell STAR `SJ.out.tab`
- 10x BAM starts from pooled RG-tagged BAM and splits to per-cell BAM cache
- 10x junction starts from per-cell `.jcounts` and backfills STAR-style splice annotation fields from reference/GTF

## Current Runtime Assumptions

- `STAR`, `samtools`, `bedtools` available on `PATH` or via env vars
- direct-junction routes require reference FASTA + GTF
- Outrigger compatibility layer is loaded via:
  - `runtime_support/pandas_compat`
  - `runtime_support/outrigger_patched`

## Recommended Entry Points

- full-length BAM:
  - `DOLPHIN/alternative_splicing/run_full_length.py`
- full-length junction:
  - `DOLPHIN/alternative_splicing/run_full_length_junction.py`
- 10x BAM:
  - `DOLPHIN/alternative_splicing/run_tenx.py`
- 10x junction:
  - `DOLPHIN/alternative_splicing/run_tenx_junction.py`

## Differential AS Metadata

Differential AS uses two distinct metadata roles:

- `--cluster-name` identifies the cell-type column used to impute each missing PSI value from the mean of that same event in the same cell type.
- `--das-group-column` identifies the biological condition used for the two-sided Wilcoxon rank-sum test.
- `--das-group1` and `--das-group2` select the comparison. They may be omitted only when the group column contains exactly two values.

P values are adjusted across tested events with the Benjamini-Hochberg procedure. The imputed matrix is retained as `<out_name>_PSI_DAS.h5ad`; event statistics, `delta_psi`, raw P values, and BH-adjusted P values are retained as `<out_name>_DAS.csv`.

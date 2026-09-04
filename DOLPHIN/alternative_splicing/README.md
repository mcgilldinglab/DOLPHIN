# Alternative Splicing Module

This module preserves the finalized DOLPHIN alternative splicing workflows used in validation:

- `full-length` with BAM aggregation
- `full-length` with direct junction aggregation
- `10x` with BAM aggregation
- `10x` with direct junction aggregation

## Code Layout

Shared across full-length and 10x:

- `pipeline.py`
  - common orchestration for neighbor finding, read counting, aggregation, Outrigger stages, and post-processing
- `find_cell_neighbor.py`
  - shared nearest-neighbor discovery from DOLPHIN embeddings
- `get_single_bam_reads.py`
  - shared per-cell BAM read counting
- `convert_psi_to_h5ad.py`
  - shared `outrigger_summary.csv -> PSI.h5ad` conversion
- `convert_random_psi.py`
  - shared randomized PSI export
- `generate_differential_as.py`
  - event-by-cell-type PSI imputation and Wilcoxon/BH differential testing
- `build_gffutils_db.py`
  - shared one-time reference DB builder
- `compare_bam_vs_junction_as.py`
  - shared BAM-vs-junction comparison helper
- `runtime_support/`
  - shared runtime compatibility layer for pandas / Outrigger

Route-specific:

- `process_reads_aggregation.py`
  - BAM aggregation route
  - used by both `full-length` and `10x`
- `process_junction_aggregation.py`
  - direct junction aggregation route
  - used by both `full-length` and `10x`

Modality-specific:

- `presets.py`
  - explicit modality defaults
  - `FULL_LENGTH_*`
  - `TENX_*`
- `run_full_length.py`
  - full-length BAM route launcher
- `run_full_length_junction.py`
  - full-length direct-junction route launcher
- `run_tenx.py`
  - 10x BAM route launcher
- `run_tenx_junction.py`
  - 10x direct-junction route launcher

## Shared vs Different Behavior

Common logic for both modalities:

1. Find neighbors in embedding space.
2. Count per-cell reads.
3. Aggregate evidence using one of two routes:
   - BAM route
   - junction route
4. Run Outrigger:
   - `index`
   - `validate`
   - `psi`
5. Convert outputs to:
   - `PSI.h5ad`
   - `PSI_random.h5ad`
   - `PSI_DAS.h5ad`
   - `DAS.csv`

full-length-specific input assumptions:

- input already includes per-cell BAM
- input already includes per-cell STAR `SJ.out.tab`
- metadata does not require barcode normalization

10x-specific input assumptions:

- BAM route requires RG-tagged pooled BAM split into per-cell BAM cache
- direct junction route consumes per-cell `.jcounts`
- metadata requires barcode normalization (`tenx_barcode`)
- direct junction route backfills `annotated`, `motif`, and `strand` from reference/GTF because `.jcounts` do not carry STAR splice semantics

## Important Boundaries

- Outrigger event logic is not changed.
- Outrigger runtime compatibility and DB-generation behavior are patched only to make the tool stable on modern environments and large runs.
- `max_overhang` is preserved as `0` in the direct-junction 10x route; it is not inferred from coordinates.

Differential AS uses `--cluster-name` for cell-type-specific, event-specific PSI imputation and `--das-group-column` for the biological comparison. Use `--das-group1` and `--das-group2` when the comparison column contains more than two values.

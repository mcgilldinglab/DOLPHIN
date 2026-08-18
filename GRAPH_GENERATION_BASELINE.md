# Graph Generation Baseline

This file records the current clean handoff from preprocessing into graph generation.

## Canonical Step2 Inputs

Step2 should now start from grouped featureCounts outputs for both full-length and 10x:

- `grouped gene count`
- `grouped exon count`
- `grouped junction count`
- `cell manifest` with a `CB` column
- shared reference bundle generated in step0:
  - `dolphin.exon.pkl`
  - `dolphin_adj_index.csv`
  - `dolphin_gene_meta.csv`
  - `dolphin_adj_metadata_table.csv`

Legacy per-cell count files are kept only for backward compatibility and validation.

## Preferred Complete-Matrix Entry

If the user only needs the final all-cell graph matrices, the preferred stable entry
point is now:

- `DOLPHIN.graph_generation.run_graph_matrix_generation`

This wrapper:

- starts from grouped exon/junction counts plus `cell_manifest.tsv`
- produces the retained outputs:
  - `Feature_<out_name>.h5ad`
  - `FeatureComp_<out_name>.h5ad`
  - `Adjacency_<out_name>.h5ad`
- removes `06_graph_mtx/` after step3 when `retain_graph_csv=False`
- keeps the regression-validated internal logic unchanged

The more aggressive direct-path prototype that skips `06_graph_mtx/` entirely is
not yet promoted, because grouped TSV loading is still too memory-heavy.

## Stable Baselines

Reference bundle:

- `/mnt/md0/kailu/DOLPHIN_codex/reference_inputs/grch38_release107/dolphin_exon_gtf`

Full-length grouped counts:

- `/mnt/md0/kailu/DOLPHIN_codex/execution_runs/full_length_flashseq_clean_v1/cell_manifest.tsv`
- `/mnt/md0/kailu/DOLPHIN_codex/execution_runs/full_length_flashseq_clean_v1/outputs/04_exon_gene_cnt_grouped/full_length.exongene.count.txt`
- `/mnt/md0/kailu/DOLPHIN_codex/execution_runs/full_length_flashseq_clean_v1/outputs/05_exon_junct_cnt_grouped/full_length.exon.count.txt`
- `/mnt/md0/kailu/DOLPHIN_codex/execution_runs/full_length_flashseq_clean_v1/outputs/05_exon_junct_cnt_grouped/full_length.exon.count.txt.jcounts`

10x grouped counts:

- `/mnt/md0/kailu/DOLPHIN_codex/execution_runs/tenx_srr8513796_clean_v2/cell_manifest.tsv`
- `/mnt/md0/kailu/DOLPHIN_codex/execution_runs/tenx_srr8513796_clean_v2/cell_ids.txt`
- `/mnt/md0/kailu/DOLPHIN_codex/execution_runs/tenx_srr8513796_clean_v2/outputs/SRR8513796/05_exon_gene_cnt_grouped/SRR8513796.exongene.count.txt`
- `/mnt/md0/kailu/DOLPHIN_codex/execution_runs/tenx_srr8513796_clean_v2/outputs/SRR8513796/06_exon_junct_cnt_grouped/SRR8513796.exon.count.txt`
- `/mnt/md0/kailu/DOLPHIN_codex/execution_runs/tenx_srr8513796_clean_v2/outputs/SRR8513796/06_exon_junct_cnt_grouped/SRR8513796.exon.count.txt.jcounts`

## Junction Zero Rows

Grouped junction tables can contain many zero values across cells after merging. That is acceptable.

The important rule is:

- when a single cell is extracted from a grouped junction table, rows with `Count <= 0` must be dropped before graph construction

This is already handled in:

- `/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/graph_generation/grouped_featurecounts.py`

## Validation Already Preserved

Validation summaries are kept in:

- `/mnt/md0/kailu/DOLPHIN_codex/tmp_validation/summaries`

Key checks already passed:

- full-length grouped input reproduces legacy graph outputs exactly on the retained validation subset
- 10x grouped counting reproduces split-per-cell featureCounts on the retained smoke-test subset
- multi-sample grouped merge preserves exact per-cell extraction on the retained smoke-test subset

## Step2 Execution Order

The current downstream order remains:

1. `run_parallel_gene_processing`
2. `run_feature_combination`
3. `run_adjacency_combination`
4. `run_adjacency_compression`
5. `run_adjacency_compress_combination`
6. `run_adjacency_matrix_final`
7. `run_raw_gene`
8. `run_feature_hvg`
9. `run_adjacency_hvg`
10. `run_model_input`

## Multi-Sample Rule

For multiple 10x samples:

- preprocess each sample independently
- merge grouped counts after preprocessing
- ensure cell IDs are globally unique before entering step2

The merged grouped matrices are then the canonical step2 input.

# Environment And Validation

Use this reference before running DOLPHIN commands, debugging failures, optimizing performance, or writing installation instructions.

## Environment Contract

The package includes `environment.yaml`. Treat it as the authoritative Python environment unless `setup.py` is later expanded with complete `install_requires`.

Important Python dependencies include:

- `anndata`
- `pandas`
- `scanpy`
- `igraph`
- `leidenalg`
- `numpy`
- `scipy`
- `scikit-learn`
- `torch`
- `torch-geometric`
- `pyro-ppl`
- `h5py`
- `pybedtools`

Important external binaries:

- `STAR`
- `samtools`
- `bedtools`
- `featureCounts` from Subread
- `Rscript`

Important R packages for EDEG/JDEG:

- `Seurat`
- `MAST`
- `reticulate`
- `dplyr`
- `patchwork`

Do not assume `/usr/bin/python3` has the required scientific stack. If imports fail with missing `anndata`, `pandas`, `pyro`, `torch_geometric`, `igraph`, or `leidenalg`, check whether the conda environment from `environment.yaml` is active.

## Runtime Roots And Scratch

Prefer configurable local scratch for heavy intermediate I/O, especially Outrigger work directories and AS aggregation. Do not hard-code one server path as the only option.

Robust behavior:

- Allow users to set scratch/work roots explicitly.
- Use local SSD or `/tmp` when available and sufficiently large.
- Copy or sync final retained outputs back to the durable results directory.
- Record where temporary work was performed in status or summary JSON when possible.

Before long runs, check free space for both scratch and final output roots.

### Alternative-Splicing Path Configuration

Do not place developer workstation or server paths in AS presets. Accept user paths through CLI flags or their environment-variable equivalents:

- `--embedding-h5ad` or `DOLPHIN_AS_EMBEDDING_H5AD`
- `--metadata-path` or `DOLPHIN_AS_METADATA_PATH`
- `--bam-root` or `DOLPHIN_AS_BAM_ROOT`
- `--junction-root` or `DOLPHIN_AS_JUNCTION_ROOT`
- `--gtf-path` or `DOLPHIN_AS_GTF_PATH`
- `--gffutils-db` or `DOLPHIN_AS_GFFUTILS_DB`
- `--genome-sizes-path` or `DOLPHIN_AS_GENOME_SIZES_PATH`
- `--fasta-path` or `DOLPHIN_AS_FASTA_PATH`
- `--star-index-dir` or `DOLPHIN_AS_STAR_INDEX_DIR`
- `--star-binary`, `--samtools-binary`, and `--bedtools-bin-dir`, or their documented environment variables/active `PATH`
- `--results-root`, `--logs-root`, `--prepared-inputs-root`, and `--outrigger-work-root`, or their documented environment variables

The BAM route requires a standard-GTF STAR index. The junction route does not require a STAR index. Both routes retain the same downstream Outrigger event/index/PSI semantics.

## Running And Debugging

For planned pipeline execution:

- Record command, input paths, output paths, route, cell count, worker/thread settings, and start/end time.
- Keep status JSON and summary JSON when the pipeline supports them.
- Monitor stderr/stdout directly for long external tools.

For failed external-tool stages:

- Do not rely on wrappers that hide errors.
- Re-run the failing external command directly with visible stdout/stderr.
- Preserve or inspect partial output before deleting it.
- Delete only the failed stage output unless broader cleanup is required and approved.

Outrigger-specific debugging:

- `outrigger index` is external Outrigger logic; do not alter its event construction semantics.
- It is acceptable to use packaged compatibility support for pandas/gffutils/SQLite stability.
- If a run repeatedly stops without useful logs, run the exact Outrigger command directly rather than through a high-level wrapper.

## Smoke Tests

Smoke tests must use the same pipeline logic as full runs. Use cell limits or small manifests to reduce scope, but do not switch from BAM to junction, skip steps, or change aggregation semantics just to make a smoke test run.

Useful smoke validation:

- Confirm expected files are produced.
- Confirm matrices are non-empty when the selected subset should contain signal.
- Confirm status/summary JSON says completed.
- Confirm no hidden traceback or silent failure occurred.

## Performance Optimization

Optimization is allowed in DOLPHIN orchestration and aggregation code when outputs are preserved.

Rules:

- Keep the old route as a baseline until the optimized route is validated.
- Use the same input, same cell subset, same parameters, and same output contract.
- Record per-step timings before and after optimization.
- Compare output shapes, IDs, nonzero counts, and numeric values.
- Use exact equality where expected; use correlation only for intentionally different routes such as BAM aggregation vs direct junction aggregation.

Do not optimize by changing Outrigger event/index/PSI semantics.

## Common Validation Targets

Preprocessing/graph/model:

- grouped count files exist
- pooled H5AD outputs exist
- graph-store manifest and edge store exist
- model training writes `DOLPHIN_Z.h5ad`
- downstream clustering/ARI can use the embedding and ground-truth cell type labels

Alternative splicing:

- `N_<out_name>_<k>.csv`
- `<out_name>_read_counts.csv`
- `aggregated_star/`
- `outrigger_output/psi/outrigger_summary.csv`
- `<out_name>_PSI.h5ad`
- `<out_name>_PSI_random.h5ad`
- `<out_name>_PSI_DAS.h5ad`
- `<out_name>_DAS.csv`
- `bam_vs_junction_comparison.json` when comparing routes

EDEG/JDEG:

- converted `.rds`
- one `<cluster>_MAST_DE.csv` per requested cluster
- `<cluster>_EDEG.csv`
- `<cluster>_JDEG.csv`

GitHub readiness:

- run Python syntax validation
- verify `find_packages()` includes runtime support packages
- verify console entrypoints resolve in a dependency-complete environment
- remove `__pycache__`, `.pyc`, and large runtime result files before final packaging

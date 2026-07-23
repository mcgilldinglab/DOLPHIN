# DOLPHIN GitHub-Ready Clean Package

This directory is the curated final package for the current DOLPHIN pipeline.

It contains only the code paths that were used for the finalized end-to-end runs:

- `preprocessing`
  - reference preparation
  - full-length raw preprocessing
  - 10x raw preprocessing
  - graph generation and model-input preparation
- `alternative_splicing`
  - full-length BAM route
  - full-length direct-junction route
  - 10x BAM route
  - 10x direct-junction route
- `EDEG`
  - historical one-cluster-vs-rest exon/junction differential workflow
  - industrialized exon/junction aggregation entrypoints
- `model`
  - DOLPHIN training
  - lazy graph-store loading
  - final embedding export

It intentionally excludes legacy, duplicated, cached, compiled, and experimental files that were not part of the final retained pipeline.

## Directory Layout

- [`DOLPHIN/preprocess`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/preprocess)
- [`DOLPHIN/graph_generation`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/graph_generation)
- [`DOLPHIN/alternative_splicing`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/alternative_splicing)
- [`DOLPHIN/EDEG`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/EDEG)
- [`DOLPHIN/model`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/model)
- [`docs/preprocessing`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/docs/preprocessing)
- [`docs/alternative_splicing`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/docs/alternative_splicing)
- [`docs/edeg`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/docs/edeg)
- [`docs/model`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/docs/model)

## What This Package Preserves

- unified preprocessing for full-length and 10x
- unified graph generation route
- dual AS routes for both full-length and 10x
  - BAM aggregation
  - direct junction aggregation
- exon-level and junction-level differential gene analysis helpers
- sparse pooled graph-matrix outputs
- lightweight graph-store model input
- direct model training to `DOLPHIN_Z.h5ad`

## What This Package Does Not Preserve

- `build/`
- `__pycache__/`
- compiled Cython artifacts
- experimental Cython acceleration files
- unused legacy helper modules
- historical execution scripts and one-off validation code

## Reference Documents

- outputs and handoff contract: [docs/PIPELINE_OUTPUTS.md](/mnt/md0/kailu/DOLPHIN_codex/github_ready/docs/PIPELINE_OUTPUTS.md)
- original source coverage map: [docs/ORIGINAL_FUNCTION_COVERAGE.md](/mnt/md0/kailu/DOLPHIN_codex/github_ready/docs/ORIGINAL_FUNCTION_COVERAGE.md)
- alternative splicing outputs: [docs/ALTERNATIVE_SPLICING_OUTPUTS.md](/mnt/md0/kailu/DOLPHIN_codex/github_ready/docs/ALTERNATIVE_SPLICING_OUTPUTS.md)
- preprocessing notes: [docs/preprocessing/README.md](/mnt/md0/kailu/DOLPHIN_codex/github_ready/docs/preprocessing/README.md)
- alternative splicing notes: [docs/alternative_splicing/README.md](/mnt/md0/kailu/DOLPHIN_codex/github_ready/docs/alternative_splicing/README.md)
- EDEG notes: [docs/edeg/README.md](/mnt/md0/kailu/DOLPHIN_codex/github_ready/docs/edeg/README.md)
- model notes: [docs/model/README.md](/mnt/md0/kailu/DOLPHIN_codex/github_ready/docs/model/README.md)
- curation notes: [CURATION_NOTES.md](/mnt/md0/kailu/DOLPHIN_codex/github_ready/CURATION_NOTES.md)

# Curation Notes

This clean package keeps only the finalized pipeline that was used for the retained full-length and 10x runs.

## Retained Code

- [`DOLPHIN/preprocess`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/preprocess)
- [`DOLPHIN/graph_generation`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/graph_generation)
- [`DOLPHIN/model`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/model)

## Explicitly Excluded

- `build/`
- `DOLPHIN.egg-info/`
- `__pycache__/`
- `example_runs/`
- `AS/`
- `EDEG/`
- `cell_reads_aggregation/`
- [`DOLPHIN/preprocess/preprocess_raw.py`](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/preprocess_raw.py)
- experimental Cython files for `count_adj`
  - `_count_adj_cython.pyx`
  - `_count_adj_cython.c`
  - `_count_adj_cython*.so`

## Rationale

- `preprocess_raw.py` is not part of the retained final preprocessing entrypoints.
- the Cython path was tested and not retained as the default implementation.
- compiled artifacts should not be committed into the clean GitHub package.
- historical benchmark and validation helpers belong in execution/validation directories, not in the distributable package.

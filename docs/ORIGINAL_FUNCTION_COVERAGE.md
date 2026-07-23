# Original Function Coverage

This document maps the original source tree under
[`DOLPHIN_original_clone/DOLPHIN`](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_original_clone/DOLPHIN)
to the finalized engineering state.

## Coverage Table

| Original scope | Original location | Final status | Final location | Notes |
|---|---|---|---|---|
| Reference and raw preprocessing | `preprocess/` | finalized | [`github_ready/DOLPHIN/preprocess`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/preprocess) | Includes reference preparation, full-length preprocessing, 10x preprocessing, CLI, and output contract. |
| Graph generation | `graph_generation/` | finalized | [`github_ready/DOLPHIN/graph_generation`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/graph_generation) | Includes pooled feature/adjacency generation, compressed matrices, HVG filtering, and model-input preparation. |
| DOLPHIN model | `model/` | finalized | [`github_ready/DOLPHIN/model`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/model) | Includes lazy graph-store loading, training, and final embedding export. |
| Alternative splicing orchestration | `AS/` | superseded by finalized module | [`github_ready/DOLPHIN/alternative_splicing`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/alternative_splicing) | Historical `AS/` pieces are merged into a single module with explicit full-length/10x and BAM/junction routes. |
| Cell-neighbor read aggregation | `cell_reads_aggregation/` | superseded by finalized module | [`github_ready/DOLPHIN/alternative_splicing`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/alternative_splicing) | Historical read-aggregation utilities now live under the same alternative splicing package. |
| Exon-level differential gene analysis | `EDEG/` | finalized with historical MAST route preserved | [`github_ready/DOLPHIN/EDEG`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/EDEG) | Public route is historical one-cluster-vs-rest only. Historical EDEG input is `FeatureComp_<out_name>.h5ad -> FeatureComp_<out_name>.rds -> <cluster>_MAST_DE.csv -> <cluster>_EDEG.csv`. |
| Reusable gene-level aggregation helper | `edeg_analysis.py` | retained as helper | [`github_ready/DOLPHIN/EDEG/edeg_analysis.py`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/EDEG/edeg_analysis.py) | Kept as a reusable aggregation implementation, but not exposed as a separate public workflow. |

## Current Engineering State

- `preprocess`: finalized
- `graph_generation`: finalized
- `model`: finalized
- `alternative_splicing`: functionally finalized and curated as the retained replacement for historical `AS/` and `cell_reads_aggregation/`
- `EDEG/JDEG`: functionally finalized with historical one-cluster-vs-rest MAST preserved as the public contract

## Legacy Compatibility Notes

- Historical directories such as `AS/` and `cell_reads_aggregation/` may still exist in working trees for traceability, but the retained implementation is the unified
  [`alternative_splicing`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/alternative_splicing)
  package.
- Historical `EDEG` execution evidence showed multiple ad hoc MAST styles. The retained public route is the historical one-cluster-vs-rest workflow only.
- `build/`, `build/lib`, `__pycache__`, compiled artifacts, and one-off recovered files are not part of the retained public package.

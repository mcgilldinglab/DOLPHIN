# Preprocessing

This folder covers the finalized preprocessing path:

- reference preparation
- full-length raw preprocessing
- 10x raw preprocessing
- graph generation to model input

## Code Roots

- [`DOLPHIN/preprocess`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/preprocess)
- [`DOLPHIN/graph_generation`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/graph_generation)

## Final Entry Points

- reference / preprocessing CLI:
  - [`DOLPHIN/preprocess/cli.py`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/preprocess/cli.py)
- preprocessing plan / execution:
  - [`DOLPHIN/preprocess/pipeline.py`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/preprocess/pipeline.py)
- graph generation:
  - [`DOLPHIN/graph_generation/preprocess_raw_reads.py`](/mnt/md0/kailu/DOLPHIN_codex/github_ready/DOLPHIN/graph_generation/preprocess_raw_reads.py)

## Final Design Decisions

- full-length and 10x share the same retained graph-generation contract
- grouped count matrices are canonical preprocessing outputs
- graph generation keeps pooled H5AD outputs, not per-cell CSV outputs
- adjacency compression uses the sparse retained path, not the earlier dense temporary layout
- model input uses graph-store files, not a mandatory large serialized `.pt`

## Output Contract

See [docs/PIPELINE_OUTPUTS.md](/mnt/md0/kailu/DOLPHIN_codex/github_ready/docs/PIPELINE_OUTPUTS.md).

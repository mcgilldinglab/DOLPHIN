# Model

This folder covers the finalized DOLPHIN model path.

## Code Root

- [`DOLPHIN/model`](../../DOLPHIN/model)

## Final Entry Points

- main API:
  - [`DOLPHIN/model/run_model.py`](../../DOLPHIN/model/run_model.py)
- training loop:
  - [`DOLPHIN/model/train.py`](../../DOLPHIN/model/train.py)
- graph-store loader:
  - [`DOLPHIN/model/lazy_graph_store.py`](../../DOLPHIN/model/lazy_graph_store.py)
- model definition:
  - [`DOLPHIN/model/model.py`](../../DOLPHIN/model/model.py)

## Final Design Decisions

- model architecture is unchanged relative to the retained default pipeline
- training accepts lightweight graph-store input
- runtime still reconstructs standard `torch_geometric.data.Data`
- final embedding is written to `DOLPHIN_Z.h5ad`

## Required Inputs

- `FeatureCompHvg_<out_name>.h5ad`
- `model_<out_name>.graph.json`
- `model_<out_name>.edges.h5`

## Main Output

- `DOLPHIN_Z.h5ad`

## Output Contract

See [the pipeline output contract](../PIPELINE_OUTPUTS.md).

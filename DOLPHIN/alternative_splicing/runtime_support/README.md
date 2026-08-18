# Runtime Support

This directory contains the compatibility layer required to run the finalized AS pipeline reproducibly.

Contents:

- `pandas_compat/sitecustomize.py`
  - pandas 3 compatibility for older Outrigger code
  - bedtools path bootstrapping
  - sqlite timeout / `FeatureDB.update` pragma patch for `gffutils`
- `outrigger_patched/`
  - patched Outrigger runtime used by finalized runs

Notes:

- The finalized 10x junction route relied on `gffutils` timeout / pragma behavior.
- In this clean package, that behavior is made explicit in `sitecustomize.py` instead of remaining only as a hidden environment edit.
- The AS module should be run with `outrigger_pythonpath` pointing at this directory before any external site-package path.

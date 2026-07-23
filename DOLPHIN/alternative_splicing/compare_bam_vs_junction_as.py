import argparse
import json
from pathlib import Path

import anndata
import numpy as np


def _to_dense(matrix):
    if hasattr(matrix, "toarray"):
        return matrix.toarray()
    return np.asarray(matrix)


def _pearson_flat(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 2:
        return None
    xv = x[mask].astype(float)
    yv = y[mask].astype(float)
    if np.std(xv) == 0 or np.std(yv) == 0:
        return None
    return float(np.corrcoef(xv, yv)[0, 1])


def _compare_h5ad(bam_path, junction_path):
    bam = anndata.read_h5ad(bam_path)
    junction = anndata.read_h5ad(junction_path)

    shared_obs = [obs for obs in bam.obs_names if obs in junction.obs_names]
    shared_var = [var for var in bam.var_names if var in junction.var_names]

    bam_sub = bam[shared_obs, shared_var]
    junction_sub = junction[shared_obs, shared_var]

    x_bam = _to_dense(bam_sub.X)
    x_junction = _to_dense(junction_sub.X)

    finite_mask = np.isfinite(x_bam) & np.isfinite(x_junction)
    if finite_mask.any():
        abs_diff = np.abs(x_bam[finite_mask] - x_junction[finite_mask])
        max_abs_diff = float(abs_diff.max())
        mean_abs_diff = float(abs_diff.mean())
    else:
        max_abs_diff = 0.0
        mean_abs_diff = 0.0

    return {
        "bam_shape": [int(bam.n_obs), int(bam.n_vars)],
        "junction_shape": [int(junction.n_obs), int(junction.n_vars)],
        "shared_obs": len(shared_obs),
        "shared_var": len(shared_var),
        "exact_equal": bool(np.array_equal(x_bam, x_junction, equal_nan=True)),
        "max_abs_diff": max_abs_diff,
        "mean_abs_diff": mean_abs_diff,
        "pearson_flat": _pearson_flat(x_bam.ravel(), x_junction.ravel()),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam-root", required=True)
    parser.add_argument("--junction-root", required=True)
    parser.add_argument("--out-json", required=True)
    parser.add_argument("--bam-prefix", required=True)
    parser.add_argument("--junction-prefix", required=True)
    args = parser.parse_args()

    bam_root = Path(args.bam_root)
    junction_root = Path(args.junction_root)

    comparisons = {
        "PSI": _compare_h5ad(
            bam_root / "alternative_splicing" / f"{args.bam_prefix}_PSI.h5ad",
            junction_root / "alternative_splicing" / f"{args.junction_prefix}_PSI.h5ad",
        ),
        "PSI_random": _compare_h5ad(
            bam_root / "alternative_splicing" / f"{args.bam_prefix}_PSI_random.h5ad",
            junction_root / "alternative_splicing" / f"{args.junction_prefix}_PSI_random.h5ad",
        ),
        "PSI_DAS": _compare_h5ad(
            bam_root / "alternative_splicing" / f"{args.bam_prefix}_PSI_DAS.h5ad",
            junction_root / "alternative_splicing" / f"{args.junction_prefix}_PSI_DAS.h5ad",
        ),
    }

    out_json = Path(args.out_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(comparisons, indent=2))
    print(json.dumps(comparisons, indent=2))


if __name__ == "__main__":
    main()

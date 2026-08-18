import os
import shutil

import h5py
import pandas as pd
from anndata.experimental import concat_on_disk

from ._anndata_compat import enable_nullable_string_writes


def _ordered_adj_comp_paths(sample_ids, temp_out_dir):
    shard_manifest = os.path.join(temp_out_dir, "adj_comp_shards", "manifest.tsv")
    if os.path.exists(shard_manifest):
        shard_df = pd.read_csv(shard_manifest, sep="\t")
        shard_df = shard_df.sort_values("shard_index", kind="stable").reset_index(drop=True)
        return [
            os.path.join(temp_out_dir, "adj_comp_shards", file_name)
            for file_name in shard_df["file_name"].tolist()
        ]

    adj_comp_dir = os.path.join(temp_out_dir, "adj_comp_matrix")
    file_paths = []
    missing = []
    for sample_id in sample_ids:
        path = os.path.join(adj_comp_dir, sample_id + ".h5ad")
        if os.path.exists(path):
            file_paths.append(path)
        else:
            missing.append(sample_id)
    if missing:
        preview = ", ".join(missing[:10])
        raise FileNotFoundError(
            "Missing compressed adjacency files for "
            f"{len(missing)} cells. Examples: {preview}"
        )
    return file_paths


def _copy_var_metadata(source_path, destination_path):
    # `concat_on_disk(..., merge=None)` keeps the sparse matrix and obs on disk
    # efficiently, then we copy the reference var metadata from one input file.
    with h5py.File(destination_path, "a") as dst, h5py.File(source_path, "r") as src:
        for key in ("var", "varm", "varp"):
            if key in dst:
                del dst[key]
            if key in src:
                src.copy(key, dst)


def _remove_path(path):
    if os.path.isdir(path):
        shutil.rmtree(path)
    elif os.path.exists(path):
        os.remove(path)


def _temp_h5ad_path(final_output_path):
    if final_output_path.endswith(".h5ad"):
        return final_output_path[: -len(".h5ad")] + ".partial.h5ad"
    return final_output_path + ".partial.h5ad"

def run_adjacency_compress_combination(
    metadata_path: str,
    out_name: str,
    out_directory: str = "./",
    adj_run_num: int = 50,
    clean_temp: bool = True,
    parallel: bool = True,
    max_processes: int | None = None,
    max_loaded_elems: int = 100_000_000,
):
    """
    Combine compressed adjacency matrices into a final AnnData object.

    Parameters
    ----------
    metadata_path : str
        Path to the metadata file with cell barcodes.
    out_name : str
        Output name prefix.
    out_directory : str
        Output folder to save results.
    adj_run_num : int
        Retained for backward compatibility. The optimized implementation no
        longer materializes intermediate batch-level h5ad files.
    clean_temp : bool
        Whether to delete temporary intermediate batch files.
    parallel : bool
        Retained for backward compatibility.
    max_processes : int | None
        Retained for backward compatibility.
    max_loaded_elems : int
        Passed to `anndata.experimental.concat_on_disk` to bound sparse-array
        loading in memory during concatenation.
        
    Returns
    -------
    None
        Saves the compressed adjacency matrix to the output directory as `AdjacencyComp_<out_name>.h5ad`.

    """
    print("Start Combining Compressed Adjacency Matrix...")

    df_label = pd.read_csv(metadata_path, sep='\t')
    sample_ids = list(df_label["CB"])

    final_out_dir = os.path.join(out_directory, "data")
    temp_out_dir = os.path.join(final_out_dir, "temp")
    os.makedirs(temp_out_dir, exist_ok=True)
    final_output_path = os.path.join(final_out_dir, f"AdjacencyComp_{out_name}.h5ad")
    temp_output_path = _temp_h5ad_path(final_output_path)

    input_files = _ordered_adj_comp_paths(sample_ids, temp_out_dir)

    _remove_path(temp_output_path)
    _remove_path(final_output_path)

    enable_nullable_string_writes()
    concat_on_disk(
        input_files,
        temp_output_path,
        axis=0,
        join="inner",
        merge=None,
        index_unique=None,
        max_loaded_elems=max_loaded_elems,
    )
    _copy_var_metadata(input_files[0], temp_output_path)
    os.replace(temp_output_path, final_output_path)

    # 4. Clean up temporary files
    if clean_temp:
        print("Cleaning up temporary files...")
        shutil.rmtree(temp_out_dir)

from multiprocessing import Pool
import logging
import os
import shutil
import sys

import anndata
from anndata.experimental import concat_on_disk
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from ._anndata_compat import enable_nullable_string_writes
from .func_preprocess_raw_reads import (
    build_reference_lookups,
    configure_grouped_count_sources,
    gene,
    get_gtf,
)
from .func_step01_fea_mat_main_part1 import build_feature_var_annotation
from .func_step01_fea_mat_main_part2 import fea_comp
from .func_step02_adj_mat_main_part1_main_1 import build_adjacency_var_annotation
from .grouped_featurecounts import load_grouped_featurecounts_table


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

_WORKER_CONTEXT = {}


def _ensure_metadata(metadata_path: str) -> pd.DataFrame:
    metadata = pd.read_csv(metadata_path, sep="\t")
    if "CB" not in metadata.columns:
        raise ValueError(f"Metadata file must contain a 'CB' column: {metadata_path}")
    metadata["CB"] = metadata["CB"].astype(str)
    return metadata


def _candidate_count_directories(out_directory, grouped_exon_count_path, grouped_junction_count_path):
    candidates = []
    for candidate in (
        os.path.join(out_directory, "05_exon_junct_cnt"),
        os.path.join(out_directory, "06_exon_junct_cnt"),
    ):
        candidates.append(candidate)

    for grouped_path in (grouped_exon_count_path, grouped_junction_count_path):
        if not grouped_path:
            continue
        grouped_dir = os.path.dirname(os.path.abspath(grouped_path))
        parent_dir = os.path.dirname(grouped_dir)
        basename = os.path.basename(grouped_dir)
        if basename.endswith("_grouped"):
            candidates.append(os.path.join(parent_dir, basename[: -len("_grouped")]))
        candidates.append(os.path.join(parent_dir, "05_exon_junct_cnt"))
        candidates.append(os.path.join(parent_dir, "06_exon_junct_cnt"))

    deduped = []
    seen = set()
    for candidate in candidates:
        if not candidate:
            continue
        normalized = os.path.abspath(candidate)
        if normalized in seen:
            continue
        deduped.append(normalized)
        seen.add(normalized)
    return deduped


def _resolve_count_directory(
    out_directory: str,
    grouped_exon_count_path: str | None = None,
    grouped_junction_count_path: str | None = None,
    count_directory: str | None = None,
) -> str | None:
    candidates = []
    if count_directory:
        candidates.append(os.path.abspath(count_directory))
    candidates.extend(
        _candidate_count_directories(
            out_directory=out_directory,
            grouped_exon_count_path=grouped_exon_count_path,
            grouped_junction_count_path=grouped_junction_count_path,
        )
    )

    for candidate in candidates:
        if not os.path.isdir(candidate):
            continue
        try:
            if any(name.endswith(".exon.count.txt") for name in os.listdir(candidate)):
                return candidate
        except OSError:
            continue
    return None


def _init_worker(
    gtf,
    adj_ind,
    main_path,
    gene_index_map,
    adj_start_by_gene,
    grouped_exon_counts,
    grouped_junction_counts,
):
    _WORKER_CONTEXT.clear()
    _WORKER_CONTEXT.update(
        {
            "gtf": gtf,
            "adj_ind": adj_ind,
            "main_path": main_path,
            "gene_index_map": gene_index_map,
            "adj_start_by_gene": adj_start_by_gene,
        }
    )
    configure_grouped_count_sources(grouped_exon_counts, grouped_junction_counts)


def _worker(args):
    index, cb = args
    try:
        g = gene(
            _WORKER_CONTEXT["gtf"],
            _WORKER_CONTEXT["adj_ind"],
            cb,
            main_path=_WORKER_CONTEXT["main_path"],
            gene_index_map=_WORKER_CONTEXT["gene_index_map"],
            adj_start_by_gene=_WORKER_CONTEXT["adj_start_by_gene"],
        )
        g.get_all()
    except Exception as exc:
        logger.error("[Error] Index %s, CB %s: %s", index, cb, exc)
        raise


def _worker_sparse(args):
    index, cb = args
    try:
        g = gene(
            _WORKER_CONTEXT["gtf"],
            _WORKER_CONTEXT["adj_ind"],
            cb,
            main_path=_WORKER_CONTEXT["main_path"],
            gene_index_map=_WORKER_CONTEXT["gene_index_map"],
            adj_start_by_gene=_WORKER_CONTEXT["adj_start_by_gene"],
        )
        return g.get_all_sparse()
    except Exception as exc:
        logger.error("[Error] Index %s, CB %s: %s", index, cb, exc)
        raise


def run_parallel_gene_processing(
    metadata_path: str,
    gtf_path: str,
    adj_index_path: str,
    main_folder: str = ".",
    grouped_exon_count_path: str | None = None,
    grouped_junction_count_path: str | None = None,
    n_processes: int | None = None,
):
    print("Starting Raw Reads Processing...")

    subfolder_5 = os.path.join(main_folder, "05_exon_junct_cnt")
    use_grouped = grouped_exon_count_path and grouped_junction_count_path
    if (not os.path.isdir(subfolder_5)) and (not use_grouped):
        print(
            "Error: required input counts not found. Provide either "
            "'05_exon_junct_cnt' under main_folder or grouped exon/junction paths."
        )
        sys.exit(1)

    subfolder_6 = os.path.join(main_folder, "06_graph_mtx")
    os.makedirs(subfolder_6, exist_ok=True)

    metadata = _ensure_metadata(metadata_path)
    cell_ids = metadata["CB"].tolist()
    gtf, adj_ind = get_gtf(gtf_path, adj_index_path)
    gene_index_map, adj_start_by_gene = build_reference_lookups(gtf, adj_ind)

    grouped_exon_counts = None
    grouped_junction_counts = None
    if use_grouped and not os.path.isdir(subfolder_5):
        grouped_exon_counts = load_grouped_featurecounts_table(
            grouped_exon_count_path,
            mode="exon",
            requested_cell_ids=cell_ids,
        )
        grouped_junction_counts = load_grouped_featurecounts_table(
            grouped_junction_count_path,
            mode="junction",
            requested_cell_ids=cell_ids,
        )

    if n_processes is None:
        n_processes = os.cpu_count() or 1

    logger.info("Running gene processing using %s processes...", n_processes)
    args_list = [(idx, cb) for idx, cb in enumerate(cell_ids)]
    with Pool(
        processes=n_processes,
        initializer=_init_worker,
        initargs=(
            gtf,
            adj_ind,
            main_folder,
            gene_index_map,
            adj_start_by_gene,
            grouped_exon_counts,
            grouped_junction_counts,
        ),
    ) as pool:
        list(pool.imap(_worker, args_list))

    configure_grouped_count_sources(None, None)
    return {"graph_directory": subfolder_6}


def _copy_var_metadata(source_path: str, destination_path: str):
    with h5py.File(destination_path, "a") as dst, h5py.File(source_path, "r") as src:
        for key in ("var", "varm", "varp"):
            if key in dst:
                del dst[key]
            if key in src:
                src.copy(key, dst)


def _remove_path(path: str):
    if os.path.isdir(path):
        shutil.rmtree(path)
    elif os.path.exists(path):
        os.remove(path)


def _temp_h5ad_path(final_output_path: str) -> str:
    if final_output_path.endswith(".h5ad"):
        return final_output_path[: -len(".h5ad")] + ".partial.h5ad"
    return final_output_path + ".partial.h5ad"


def _build_sparse_matrix_from_rows(indices_list, values_list, matrix_size: int, dtype=np.float32):
    row_count = len(indices_list)
    if row_count == 0:
        return sparse.csr_matrix((0, matrix_size), dtype=dtype)

    indptr = np.zeros(row_count + 1, dtype=np.int64)
    total_nnz = sum(len(indices) for indices in indices_list)
    combined_indices = np.empty(total_nnz, dtype=np.int32)
    combined_values = np.empty(total_nnz, dtype=dtype)
    cursor = 0
    for row_idx, (indices, values) in enumerate(zip(indices_list, values_list), start=1):
        nnz = len(indices)
        if nnz:
            combined_indices[cursor : cursor + nnz] = indices
            combined_values[cursor : cursor + nnz] = values.astype(dtype, copy=False)
            cursor += nnz
        indptr[row_idx] = cursor
    return sparse.csr_matrix(
        (combined_values[:cursor], combined_indices[:cursor], indptr),
        shape=(row_count, matrix_size),
        dtype=dtype,
    )


def _build_obs_dataframe(metadata: pd.DataFrame, cell_ids):
    batch_meta = metadata.set_index("CB").loc[list(cell_ids)].reset_index()
    obs_names = batch_meta["CB"].astype(str).values
    obs = pd.DataFrame(index=obs_names)
    for column_name in metadata.columns:
        obs[column_name] = batch_meta[column_name].values
    return obs


def _write_direct_shard(
    shard_paths: list,
    rows: list,
    metadata: pd.DataFrame,
    shard_dir: str,
    prefix: str,
    out_name: str,
    shard_index: int,
    matrix_size: int,
    var: pd.DataFrame,
    index_key: str,
    value_key: str,
):
    if not rows:
        return

    cell_ids = [row["sample_id"] for row in rows]
    indices_list = [row[index_key] for row in rows]
    values_list = [row[value_key] for row in rows]
    matrix = _build_sparse_matrix_from_rows(indices_list, values_list, matrix_size)
    obs = _build_obs_dataframe(metadata, cell_ids)
    adata = anndata.AnnData(X=matrix, obs=obs, var=var.copy())
    enable_nullable_string_writes()
    output_path = os.path.join(shard_dir, f"{prefix}_{out_name}_{shard_index}.h5ad")
    adata.write_h5ad(output_path)
    shard_paths.append(output_path)


def _concat_shards_to_final(shard_paths: list[str], final_output_path: str):
    if not shard_paths:
        raise ValueError("No shard files were generated for concatenation.")

    temp_output_path = _temp_h5ad_path(final_output_path)
    _remove_path(temp_output_path)
    _remove_path(final_output_path)

    enable_nullable_string_writes()
    concat_on_disk(
        shard_paths,
        temp_output_path,
        axis=0,
        join="inner",
        merge=None,
        index_unique=None,
    )
    _copy_var_metadata(shard_paths[0], temp_output_path)
    os.replace(temp_output_path, final_output_path)


def run_direct_graph_matrix_construction(
    metadata_path: str,
    gtf_path: str,
    adj_index_path: str,
    gene_annotation_path: str,
    adj_meta_file: str,
    out_name: str,
    out_directory: str = ".",
    grouped_exon_count_path: str | None = None,
    grouped_junction_count_path: str | None = None,
    count_directory: str | None = None,
    n_processes: int | None = None,
    feature_batch_size: int = 100,
    clean_temp: bool = True,
):
    metadata = _ensure_metadata(metadata_path)
    cell_ids = metadata["CB"].tolist()
    gtf, adj_ind = get_gtf(gtf_path, adj_index_path)
    gene_index_map, adj_start_by_gene = build_reference_lookups(gtf, adj_ind)

    resolved_count_directory = _resolve_count_directory(
        out_directory=out_directory,
        grouped_exon_count_path=grouped_exon_count_path,
        grouped_junction_count_path=grouped_junction_count_path,
        count_directory=count_directory,
    )
    linked_count_dir = None
    grouped_exon_counts = None
    grouped_junction_counts = None
    count_source = "grouped"

    if resolved_count_directory is not None:
        legacy_link = os.path.join(out_directory, "05_exon_junct_cnt")
        if not os.path.exists(legacy_link):
            os.symlink(resolved_count_directory, legacy_link)
            linked_count_dir = legacy_link
        count_source = "legacy_per_cell"
    else:
        if not grouped_exon_count_path or not grouped_junction_count_path:
            raise ValueError(
                "Unable to resolve per-cell count directory and grouped exon/junction "
                "paths were not both provided."
            )
        grouped_exon_counts = load_grouped_featurecounts_table(
            grouped_exon_count_path,
            mode="exon",
            requested_cell_ids=cell_ids,
        )
        grouped_junction_counts = load_grouped_featurecounts_table(
            grouped_junction_count_path,
            mode="junction",
            requested_cell_ids=cell_ids,
        )

    data_dir = os.path.join(out_directory, "data")
    temp_dir = os.path.join(data_dir, "temp")
    feature_shard_dir = os.path.join(temp_dir, "feature_direct_shards")
    adjacency_shard_dir = os.path.join(temp_dir, "adjacency_direct_shards")
    os.makedirs(feature_shard_dir, exist_ok=True)
    os.makedirs(adjacency_shard_dir, exist_ok=True)

    feature_var = build_feature_var_annotation(
        gene_annotation=gene_annotation_path,
        gtf_pkl_path=gtf_path,
    )
    adjacency_var = build_adjacency_var_annotation(adj_meta_file)

    if n_processes is None:
        n_processes = os.cpu_count() or 1
    shard_size = max(1, int(feature_batch_size))

    feature_shard_paths = []
    adjacency_shard_paths = []
    args_list = [(idx, cb) for idx, cb in enumerate(cell_ids)]
    with Pool(
        processes=n_processes,
        initializer=_init_worker,
        initargs=(
            gtf,
            adj_ind,
            out_directory,
            gene_index_map,
            adj_start_by_gene,
            grouped_exon_counts,
            grouped_junction_counts,
        ),
    ) as pool:
        row_buffer = []
        shard_index = 0
        for row in pool.imap(_worker_sparse, args_list):
            row_buffer.append(row)
            if len(row_buffer) >= shard_size:
                _write_direct_shard(
                    feature_shard_paths,
                    row_buffer,
                    metadata,
                    feature_shard_dir,
                    "Feature",
                    out_name,
                    shard_index,
                    len(feature_var),
                    feature_var,
                    "feature_indices",
                    "feature_values",
                )
                _write_direct_shard(
                    adjacency_shard_paths,
                    row_buffer,
                    metadata,
                    adjacency_shard_dir,
                    "Adjacency",
                    out_name,
                    shard_index,
                    len(adjacency_var),
                    adjacency_var,
                    "adjacency_indices",
                    "adjacency_values",
                )
                row_buffer = []
                shard_index += 1
        if row_buffer:
            _write_direct_shard(
                feature_shard_paths,
                row_buffer,
                metadata,
                feature_shard_dir,
                "Feature",
                out_name,
                shard_index,
                len(feature_var),
                feature_var,
                "feature_indices",
                "feature_values",
            )
            _write_direct_shard(
                adjacency_shard_paths,
                row_buffer,
                metadata,
                adjacency_shard_dir,
                "Adjacency",
                out_name,
                shard_index,
                len(adjacency_var),
                adjacency_var,
                "adjacency_indices",
                "adjacency_values",
            )

    feature_output = os.path.join(data_dir, f"Feature_{out_name}.h5ad")
    adjacency_output = os.path.join(data_dir, f"Adjacency_{out_name}.h5ad")
    _concat_shards_to_final(feature_shard_paths, feature_output)
    _concat_shards_to_final(adjacency_shard_paths, adjacency_output)

    fea_comp(data_dir, out_name)
    feature_comp_output = os.path.join(data_dir, f"FeatureComp_{out_name}.h5ad")

    if clean_temp and os.path.isdir(temp_dir):
        shutil.rmtree(temp_dir)
    if linked_count_dir and os.path.islink(linked_count_dir):
        os.unlink(linked_count_dir)

    configure_grouped_count_sources(None, None)
    return {
        "count_source": count_source,
        "feature_output": feature_output,
        "feature_comp_output": feature_comp_output,
        "adjacency_output": adjacency_output,
        "feature_shard_count": len(feature_shard_paths),
        "adjacency_shard_count": len(adjacency_shard_paths),
    }


def run_graph_matrix_generation(
    metadata_path: str,
    gtf_path: str,
    adj_index_path: str,
    gene_annotation_path: str,
    adj_meta_file: str,
    out_name: str,
    out_directory: str = ".",
    grouped_exon_count_path: str | None = None,
    grouped_junction_count_path: str | None = None,
    count_directory: str | None = None,
    n_processes: int | None = None,
    feature_batch_size: int = 100,
    adjacency_batch_size: int = 50,
    adjacency_parallel: bool = True,
    retain_graph_csv: bool = False,
    clean_temp: bool = True,
):
    if not retain_graph_csv:
        return run_direct_graph_matrix_construction(
            metadata_path=metadata_path,
            gtf_path=gtf_path,
            adj_index_path=adj_index_path,
            gene_annotation_path=gene_annotation_path,
            adj_meta_file=adj_meta_file,
            out_name=out_name,
            out_directory=out_directory,
            grouped_exon_count_path=grouped_exon_count_path,
            grouped_junction_count_path=grouped_junction_count_path,
            count_directory=count_directory,
            n_processes=n_processes,
            feature_batch_size=feature_batch_size,
            clean_temp=clean_temp,
        )

    graph_directory = run_parallel_gene_processing(
        metadata_path=metadata_path,
        gtf_path=gtf_path,
        adj_index_path=adj_index_path,
        main_folder=out_directory,
        grouped_exon_count_path=grouped_exon_count_path,
        grouped_junction_count_path=grouped_junction_count_path,
        n_processes=n_processes,
    )["graph_directory"]

    from .process_feature_matrix import run_feature_combination
    from .process_adjacency_matrix import run_adjacency_combination

    run_feature_combination(
        metadata_path=metadata_path,
        graph_directory=graph_directory,
        gene_annotation=gene_annotation_path,
        gtf_pkl_path=gtf_path,
        out_name=out_name,
        out_directory=out_directory,
        fea_run_num=feature_batch_size,
        clean_temp=clean_temp,
    )
    run_adjacency_combination(
        metadata_path=metadata_path,
        graph_directory=graph_directory,
        adj_meta_file=adj_meta_file,
        out_name=out_name,
        out_directory=out_directory,
        adj_run_num=adjacency_batch_size,
        clean_temp=clean_temp,
        parallel=adjacency_parallel,
    )

    feature_output = os.path.join(out_directory, "data", f"Feature_{out_name}.h5ad")
    feature_comp_output = os.path.join(out_directory, "data", f"FeatureComp_{out_name}.h5ad")
    adjacency_output = os.path.join(out_directory, "data", f"Adjacency_{out_name}.h5ad")

    if not retain_graph_csv and os.path.isdir(graph_directory):
        shutil.rmtree(graph_directory)

    return {
        "count_source": "grouped" if grouped_exon_count_path else "legacy_per_cell",
        "feature_output": feature_output,
        "feature_comp_output": feature_comp_output,
        "adjacency_output": adjacency_output,
        "graph_directory": graph_directory,
    }

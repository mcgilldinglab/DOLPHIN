"""Lazy exports for the DOLPHIN graph_generation package."""

from importlib import import_module


_EXPORTS = {
    "run_parallel_gene_processing": (
        "DOLPHIN.graph_generation.preprocess_raw_reads",
        "run_parallel_gene_processing",
    ),
    "run_direct_graph_matrix_construction": (
        "DOLPHIN.graph_generation.preprocess_raw_reads",
        "run_direct_graph_matrix_construction",
    ),
    "run_graph_matrix_generation": (
        "DOLPHIN.graph_generation.preprocess_raw_reads",
        "run_graph_matrix_generation",
    ),
    "run_feature_combination": (
        "DOLPHIN.graph_generation.process_feature_matrix",
        "run_feature_combination",
    ),
    "run_adjacency_combination": (
        "DOLPHIN.graph_generation.process_adjacency_matrix",
        "run_adjacency_combination",
    ),
    "run_adjacency_compression": (
        "DOLPHIN.graph_generation.process_adjacency_matrix_compress",
        "run_adjacency_compression",
    ),
    "run_adjacency_compress_combination": (
        "DOLPHIN.graph_generation.process_adjacency_matrix_compress_combine",
        "run_adjacency_compress_combination",
    ),
    "run_adjacency_matrix_final": (
        "DOLPHIN.graph_generation.process_adjacency_matrix_final",
        "run_adjacency_matrix_final",
    ),
    "run_raw_gene": (
        "DOLPHIN.graph_generation.process_raw_gene",
        "run_raw_gene",
    ),
    "run_feature_hvg": (
        "DOLPHIN.graph_generation.process_feature_hvg",
        "run_feature_hvg",
    ),
    "run_adjacency_hvg": (
        "DOLPHIN.graph_generation.process_adjacency_hvg",
        "run_adjacency_hvg",
    ),
    "run_model_input": (
        "DOLPHIN.graph_generation.process_graph_final",
        "run_model_input",
    ),
}


__all__ = sorted(_EXPORTS)


def __getattr__(name):
    if name in _EXPORTS:
        module_name, attr_name = _EXPORTS[name]
        module = import_module(module_name)
        attr = getattr(module, attr_name)
        globals()[name] = attr
        return attr
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(__all__))

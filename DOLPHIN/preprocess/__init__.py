"""Lazy exports for the DOLPHIN preprocess package.

This keeps lightweight CLI helpers importable even when optional heavy
dependencies are not needed for the current command.
"""

from importlib import import_module


_SUBMODULES = {
    "gtfpy": "DOLPHIN.preprocess.gtfpy",
    "generate_exon_gtf": "DOLPHIN.preprocess.generate_exon_gtf",
    "generate_adj_index": "DOLPHIN.preprocess.generate_adj_index",
    "pipeline": "DOLPHIN.preprocess.pipeline",
}


_EXPORTS = {
    "prepare_reference_bundle": (
        "DOLPHIN.preprocess.generate_exon_gtf",
        "prepare_reference_bundle",
    ),
    "build_star_genome_index": (
        "DOLPHIN.preprocess.generate_exon_gtf",
        "build_star_genome_index",
    ),
    "generate_nonoverlapping_exons": (
        "DOLPHIN.preprocess.generate_exon_gtf",
        "generate_nonoverlapping_exons",
    ),
    "prepare_exon_gtf": (
        "DOLPHIN.preprocess.generate_exon_gtf",
        "prepare_exon_gtf",
    ),
    "exon_uniq": (
        "DOLPHIN.preprocess.generate_exon_gtf",
        "exon_uniq",
    ),
    "save_by_batch": (
        "DOLPHIN.preprocess.generate_exon_gtf",
        "save_by_batch",
    ),
    "combine_saved_batches": (
        "DOLPHIN.preprocess.generate_exon_gtf",
        "combine_saved_batches",
    ),
    "check_exon_overlap": (
        "DOLPHIN.preprocess.generate_exon_gtf",
        "check_exon_overlap",
    ),
    "save_gtf_outputs": (
        "DOLPHIN.preprocess.generate_exon_gtf",
        "save_gtf_outputs",
    ),
    "generate_adj_index_table": (
        "DOLPHIN.preprocess.generate_adj_index",
        "generate_adj_index_table",
    ),
    "generate_adj_metadata_table": (
        "DOLPHIN.preprocess.generate_adj_index",
        "generate_adj_metadata_table",
    ),
    "build_full_length_preprocess_plan": (
        "DOLPHIN.preprocess.pipeline",
        "build_full_length_preprocess_plan",
    ),
    "build_tenx_preprocess_plan": (
        "DOLPHIN.preprocess.pipeline",
        "build_tenx_preprocess_plan",
    ),
    "write_preprocess_plan": (
        "DOLPHIN.preprocess.pipeline",
        "write_preprocess_plan",
    ),
    "split_grouped_featurecounts": (
        "DOLPHIN.preprocess.pipeline",
        "split_grouped_featurecounts",
    ),
    "merge_grouped_featurecounts_files": (
        "DOLPHIN.preprocess.pipeline",
        "merge_grouped_featurecounts_files",
    ),
    "stream_cb_to_rg": (
        "DOLPHIN.preprocess.pipeline",
        "stream_cb_to_rg",
    ),
    "rewrite_bam_cb_to_rg": (
        "DOLPHIN.preprocess.pipeline",
        "rewrite_bam_cb_to_rg",
    ),
}


__all__ = sorted(set(_SUBMODULES) | set(_EXPORTS))


def __getattr__(name):
    if name in _SUBMODULES:
        module = import_module(_SUBMODULES[name])
        globals()[name] = module
        return module

    if name in _EXPORTS:
        module_name, attr_name = _EXPORTS[name]
        module = import_module(module_name)
        attr = getattr(module, attr_name)
        globals()[name] = attr
        return attr

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(__all__))

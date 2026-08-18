FULL_LENGTH_DEFAULTS = {
    "bam_file_extension": ".Aligned.sortedByCoord.out.bam",
    "junction_file_extension": ".SJ.out.tab",
    "out_name": "full_length",
    "cluster_name": "celltype1",
    "metadata_normalization": "none",
}


TENX_DEFAULTS = {
    "bam_file_extension": ".Aligned.sortedByCoord.out.bam",
    "junction_file_extension": ".exon.count.txt.jcounts",
    "out_name": "tenx",
    "cluster_name": "celltype",
    "metadata_normalization": "tenx_barcode",
}


FULL_LENGTH_BAM_DEFAULTS = {
    **FULL_LENGTH_DEFAULTS,
    "aggregation_mode": "bam",
    "out_name": "full_length_bam",
}


FULL_LENGTH_JUNCTION_DEFAULTS = {
    **FULL_LENGTH_DEFAULTS,
    "aggregation_mode": "junction",
    "out_name": "full_length_junction",
}


TENX_BAM_DEFAULTS = {
    **TENX_DEFAULTS,
    "aggregation_mode": "bam",
    "out_name": "tenx_bam",
}


TENX_JUNCTION_DEFAULTS = {
    **TENX_DEFAULTS,
    "aggregation_mode": "junction",
    "out_name": "tenx_junction",
}

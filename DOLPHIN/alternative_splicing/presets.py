FULL_LENGTH_DEFAULTS = {
    "embedding_h5ad": "/mnt/md0/kailu/DOLPHIN_codex/execution_runs/full_length_dolphin_model_v2/data/DOLPHIN_Z.h5ad",
    "metadata_path": "/mnt/md0/kailu/DOLPHIN_codex/data_example/flashseq/metadata.csv",
    "bam_root": "/mnt/md0/kailu/DOLPHIN_codex/execution_runs/full_length_flashseq_clean_v1/outputs/03_exon_star",
    "junction_root": "/mnt/md0/kailu/DOLPHIN_codex/execution_runs/full_length_flashseq_clean_v1/outputs/03_exon_star",
    "bam_file_extension": ".Aligned.sortedByCoord.out.bam",
    "junction_file_extension": ".SJ.out.tab",
    "out_name": "full_length_flashseq",
    "cluster_name": "celltype1",
    "metadata_normalization": "none",
    "rg_bam_path": None,
}


TENX_DEFAULTS = {
    "embedding_h5ad": "/mnt/md0/kailu/DOLPHIN_codex/execution_runs/tenx_srr8513796_dolphin_model_v1/data/DOLPHIN_Z.h5ad",
    "metadata_path": "/mnt/md0/kailu/DOLPHIN_codex/data_example/10x/10x_metadata.csv",
    "bam_root": "/mnt/md0/kailu/DOLPHIN_codex/execution_runs/tenx_srr8513796_clean_v2/outputs/SRR8513796/04_rg_tagged_bam",
    "junction_root": "/mnt/md0/kailu/DOLPHIN_codex/execution_runs/tenx_srr8513796_clean_v2/outputs/SRR8513796/06_exon_junct_cnt",
    "bam_file_extension": ".Aligned.sortedByCoord.out.bam",
    "junction_file_extension": ".exon.count.txt.jcounts",
    "out_name": "tenx_srr8513796",
    "cluster_name": "celltype",
    "metadata_normalization": "tenx_barcode",
    "rg_bam_path": "/mnt/md0/kailu/DOLPHIN_codex/execution_runs/tenx_srr8513796_clean_v2/outputs/SRR8513796/04_rg_tagged_bam/SRR8513796.cb_rg.bam",
}


FULL_LENGTH_BAM_DEFAULTS = {
    **FULL_LENGTH_DEFAULTS,
    "aggregation_mode": "bam",
    "out_name": "full_length_flashseq_final",
}


FULL_LENGTH_JUNCTION_DEFAULTS = {
    **FULL_LENGTH_DEFAULTS,
    "aggregation_mode": "junction",
    "out_name": "full_length_flashseq_junction_final",
}


TENX_BAM_DEFAULTS = {
    **TENX_DEFAULTS,
    "aggregation_mode": "bam",
    "out_name": "tenx_srr8513796_final_bam",
}


TENX_JUNCTION_DEFAULTS = {
    **TENX_DEFAULTS,
    "aggregation_mode": "junction",
    "out_name": "tenx_srr8513796_final_junction",
}

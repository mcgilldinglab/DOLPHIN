import pandas as pd
import os
import shutil

def aggregate_and_normalize_bam(
    sample_list_file,
    neighbor_csv,
    single_size_file,
    src_path,
    dist_path,
    n_neighbor=10
):
    """
    Aggregates BAM and junction data from neighboring cells using majority voting
    and generates normalized BAM files for downstream splicing analysis.

    This function performs cell-level aggregation of splice junctions using DOLPHIN model embedding with KNN,
    followed by normalization of BAM file sizes through upsampling/downsampling, and filters
    junction reads based on majority-voted junctions. Final merged BAM files are created for each target cell.

    Parameters
    ----------
    sample_list_file : str
        Path to the file containing the list of sample (cell) names.
    neighbor_csv : str
        Path to CSV file listing each cell's nearest neighbors (e.g. generated from KNN).
    single_size_file : str
        Path to file listing the number of sequences for each individual single cell BAM.
    src_path : str
        Directory containing original single-cell BAM and SJ (splice junction) files.
    dist_path : str
        Destination directory where final BAM outputs will be stored.
    n_neighbor : int, optional
        Number of neighbors used for majority voting (default is 10).

    Returns
    -------
    None
        The function runs in-place and writes output BAM files to `dist_path`.

    Notes
    -----
    - Uses `samtools` for BAM filtering, merging, and sampling.
    - Temporary intermediate files are deleted after final BAM file is created.
    - Majority junctions are defined as junctions occurring in at least half of the neighbor cells.

    Example
    -------
    >>> aggregate_and_normalize_bam(
            sample_list_file="samples.txt",
            neighbor_csv="dolphin_aggregation_KNN10.csv",
            single_size_file="simu_single_bam.txt",
            src_path="/path/to/source/",
            dist_path="/path/to/output/",
        )
    """

    pd_aggr = pd.read_csv(neighbor_csv)
    pd_single_size = pd.read_csv(single_size_file)
    sample_list = list(pd.read_csv(sample_list_file, names=["CB"])["CB"])

    for target in sample_list:
        os.makedirs(os.path.join(dist_path, 'final_bam'), exist_ok=True)
        print(target)
        target_size = pd_single_size[pd_single_size["sample"] == target].iloc[0]["num_seqs"]
        _neighbor = list(pd_aggr[pd_aggr["main_name"] == target]["neighbor"])
        os.makedirs(os.path.join(dist_path, target), exist_ok=True)

        # Majority voting for junctions
        for _i, _temp_n in enumerate(_neighbor):
            _df_junc = pd.read_csv(
                os.path.join(src_path, "single_SJ", _temp_n + ".SJ.out.tab"),
                sep="\t", usecols=[0, 1, 2, 7],
                names=["chr", "first_base", "last_base", "multi_map" + _temp_n]
            )
            if _i == 0:
                df_merge = _df_junc
            else:
                df_merge = pd.merge(
                    df_merge, _df_junc, how="outer",
                    on=["chr", "first_base", "last_base"]
                )

        df_merge["nont_na"] = n_neighbor - df_merge.drop(columns=["chr", "first_base", "last_base"]).isna().sum(axis=1)
        df_keep_junct = df_merge[df_merge["nont_na"] >= (n_neighbor // 2)]
        df_keep_junct[["chr", "first_base", "last_base"]].to_csv(
            os.path.join(dist_path, target, "keep_junction.bed"),
            sep="\t", index=False, header=False
        )

        # Normalize BAM files (up/down-sampling)
        for _n in _neighbor:
            _n_seq = pd_single_size[pd_single_size["sample"] == _n].iloc[0]["num_seqs"]
            original_bam = os.path.join(src_path, "single_bam", _n + ".Aligned.sortedByCoord.out.bam")
            target_dir = os.path.join(dist_path, target)
            shutil.copyfile(original_bam, os.path.join(target_dir, _n + ".bam"))

            if _n_seq == target_size:
                os.rename(os.path.join(target_dir, _n + ".bam"), os.path.join(target_dir, _n + ".norm.bam"))
            elif _n_seq < target_size:
                _cat_self_n = int(target_size / _n_seq)
                _add_seq_perct = (target_size - _n_seq * _cat_self_n) / _n_seq
                sample_bam = os.path.join(target_dir, _n + ".sample.bam")
                os.system(f"samtools view -b -s {_add_seq_perct} {target_dir}/{_n}.bam > {sample_bam}")
                combine_list = [os.path.join(target_dir, f"{_n}.bam")] * _cat_self_n + [sample_bam]
                os.system(f"samtools merge {target_dir}/{_n}.norm.bam {' '.join(combine_list)}")
                os.remove(sample_bam)
                os.remove(os.path.join(target_dir, _n + ".bam"))
            elif _n_seq > target_size:
                _keep_seq_perct = target_size / _n_seq
                os.system(f"samtools view -b -s {_keep_seq_perct} {target_dir}/{_n}.bam > {target_dir}/{_n}.norm.bam")
                os.remove(os.path.join(target_dir, _n + ".bam"))

            # Extract junction reads and filter
            if _n != target:
                norm_bam = os.path.join(target_dir, _n + ".norm.bam")
                junction_bam = os.path.join(target_dir, _n + ".junction.norm.bam")
                mj_bam = os.path.join(target_dir, _n + ".mj.junction.norm.bam")
                os.system(f"samtools view -h {norm_bam} | awk '$0 ~ /^@/ || $6 ~ /N/' | samtools view -b > {junction_bam}")
                os.system(f"samtools index {junction_bam}")
                os.system(f"samtools view -h -L {target_dir}/keep_junction.bed {junction_bam} > {mj_bam}")

        # Merge final BAM files
        mj_bams = f"{target_dir}/*.mj.junction.norm.bam {target_dir}/{target}.norm.bam"
        final_bam = os.path.join(dist_path, "final_bam", target + ".final.bam")
        os.system(f"samtools merge {final_bam} {mj_bams}")
        shutil.rmtree(target_dir)

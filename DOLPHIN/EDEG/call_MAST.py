import os
import subprocess


def call_MAST_func(
    input_rds,
    output_dir,
    ident_column="cluster",
    cluster_values=None,
    logfc_threshold=0.0,
    normalize=True,
    rscript_binary="Rscript",
):
    """
    Run historical Seurat/MAST in one-cluster-vs-rest mode.
    """

    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script = os.path.join(script_dir, "r_MAST.R")

    if not os.path.exists(r_script):
        raise FileNotFoundError(f"Cannot find R script at: {r_script}")

    if cluster_values is None:
        cluster_values_arg = ""
    elif isinstance(cluster_values, str):
        cluster_values_arg = cluster_values
    else:
        cluster_values_arg = ",".join(map(str, cluster_values))

    subprocess.run(
        [
            rscript_binary,
            r_script,
            input_rds,
            output_dir,
            ident_column,
            cluster_values_arg,
            str(logfc_threshold),
            "TRUE" if normalize else "FALSE",
        ],
        check=True,
    )

import subprocess

r_script = "r_convert.R"
# input_file = "/mnt/data/kailu9/DOLPHIN_run_input_output/DOLPHIN_tutorial/data/FeatureComp_fsla.h5ad"
# output_rds = "/mnt/data/kailu9/DOLPHIN_run_input_output/DOLPHIN_tutorial/data/FeatureComp_fsla.rds"
# # cluster_name = "celltype"  # 可选参数

def call_MAST_func(input_anndata, output_rds, cluster_name):
    subprocess.run(["Rscript", r_script, input_anndata, output_rds, cluster_name])

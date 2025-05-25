import subprocess

r_script = "r_MAST.R"

def run_h5ad_rds(input_anndata, output_rds):
    subprocess.run(["Rscript", r_script, input_anndata, output_rds])

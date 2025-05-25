import os
import subprocess

def run_h5ad_rds(input_anndata, output_rds):
    """
    Convert a .h5ad file to a Seurat RDS object by calling an R script.

    Parameters
    ----------
    input_anndata : str
        Path to the input .h5ad file to be converted.
    
    output_rds : str
        Path where the resulting .rds file will be saved.

    Returns
    -------
    None
        The Seurat object is saved to the location specified by output_rds.
    """
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script = os.path.join(script_dir, "r_convert.R")

    if not os.path.exists(r_script):
        raise FileNotFoundError(f"Cannot find R script at: {r_script}")

    subprocess.run(["Rscript", r_script, input_anndata, output_rds], check=True)

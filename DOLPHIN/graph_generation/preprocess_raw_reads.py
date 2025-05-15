from multiprocessing import Pool
import pandas as pd
from .func_preprocess_raw_reads import gene, get_gtf
import os
import logging
import sys


# 设置 logger（替代 print）
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def _worker(args):
    """
    Worker function that processes one cell barcode.
    Receives all dependencies as part of args tuple (cb, gtf, adj_ind).
    """
    index, cb, gtf, adj_ind, main_path = args
    
    try:
        g = gene(gtf, adj_ind, cb, main_path=main_path)
        g.get_all()
    except Exception as e:
        logger.error(f"[Error] Index {index}, CB {cb}: {e}")

def run_parallel_gene_processing(
    metadata_path: str,
    gtf_path: str,
    adj_index_path: str,
    main_folder: str = ".",
    n_processes: int = None,
):
    """
    Run gene.get_all() processing in parallel across multiple cell barcodes.

    This function processes exon count and junction raw count data for each cell and converts them into
    flattened feature and adjacency vectors. Each output vector corresponds to a single cell
    and follows a consistent ordering defined by the provided GTF `.pkl` file and adjacency
    index `.csv` file. This ensures the output matrices are aligned across all cells and can
    be directly used in downstream graph-based models or statistical analysis.

    It also performs parallelization using a thread for better performance.

    Parameters
    ----------
    metadata_path : str
        Path to a metadata file (e.g., `.csv` or `.txt`) containing a column of cell barcodes (CB).
    gtf_path : str
        Path to the pickled GTF file containing exon information. Should be generated ahead of time.
    adj_index_path : str
        Path to the adjacency index CSV file. This defines adjacency matrix layout per gene.
    main_folder : str, optional
        Path to the working directory. Must contain subfolder `05_exon_junct_cnt` with count files.
        Output will be written to `06_graph_mtx` under this folder. Default is current directory `"./"`.
    n_processes : int, optional
        Number of threads or processes to run in parallel. If None, uses all available CPU cores.

    Returns
    -------
    None
        Saves the following files to the `06_graph_mtx` subdirectory inside `main_folder`:
        
        - `<cell_id>_fea.csv`: Flattened feature vector (exon counts) for each cell.
        - `<cell_id>_adj.csv`: Flattened adjacency matrix vector for each cell.

    """
    
    print("Starting Raw Reads Processing...")
    # Check that the required input subfolder exists
    subfolder_5 = os.path.join(main_folder, "05_exon_junct_cnt")
    if not os.path.isdir(subfolder_5):
        print(f"Error: Required subfolder '05_exon_junct_cnt' not found in: {main_folder}")
        sys.exit(1)  # Exit with non-zero status (indicates failure)
        
    # Create output folder if it does not exist
    subfolder_6 = os.path.join(main_folder, "06_graph_mtx")
    os.makedirs(subfolder_6, exist_ok=True)
    
    # Load metadata
    pd_gt = pd.read_csv(metadata_path, sep='\t')
    mr_cb_list = list(pd_gt["CB"])

    # Load GTF and adjacency index
    gtf, adj_ind = get_gtf(gtf_path, adj_index_path)

    # Determine number of processes
    if n_processes is None:
        n_processes = os.cpu_count()

    logger.info(f"Running gene processing using {n_processes} processes...")

    # Prepare arguments for workers
    args_list = [(i, cb, gtf, adj_ind, main_folder) for i, cb in enumerate(mr_cb_list)]

    with Pool(processes=n_processes) as pool:
        pool.map(_worker, args_list)

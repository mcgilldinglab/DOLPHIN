import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import pandas as pd
import os
from tqdm import tqdm


def _run_flagstat(bam_file: Path, samtools_binary: str):
    result = subprocess.run(
        [samtools_binary, "flagstat", str(bam_file)],
        capture_output=True,
        text=True,
        check=True,
    )
    first_line = result.stdout.splitlines()[0]
    count = int(first_line.split(" ")[0])
    return bam_file, bam_file.stem, count, result.stdout

def run_reads_count(
    out_name: str,
    bam_file_path: str,
    out_directory: str = "./",
    samtools_binary: str = "samtools",
    jobs: int = 1,
):
    """
    Count the number of reads in each BAM file using `samtools flagstat`.

    Parameters
    ----------
    out_name : str
        Prefix for output files.
    bam_file_path : str
        Directory path to search for BAM files.
    out_directory : str, optional
        Output directory to save results. Default is current directory.

    Returns
    -------
    None
        Writes two files:
        - <out_name>_flagstat_raw.txt: raw output from samtools flagstat
        - <out_name>_read_counts.csv: table with sample name and read count
    """
    os.makedirs(out_directory, exist_ok=True)

    ### step1: get the single cell bam reads number
    # Paths for output files
    raw_flagstat_path = Path(out_directory) / f"{out_name}_flagstat_raw.txt"
    summary_csv_path = Path(out_directory) / f"{out_name}_read_counts.csv"

    bam_files = list(Path(bam_file_path).rglob("*.bam"))

    results = [None] * len(bam_files)

    if jobs > 1:
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            future_map = {
                executor.submit(_run_flagstat, bam_file, samtools_binary): idx
                for idx, bam_file in enumerate(bam_files)
            }
            for future in tqdm(as_completed(future_map), total=len(future_map), desc="Processing BAM files"):
                idx = future_map[future]
                results[idx] = future.result()
    else:
        for idx, bam_file in enumerate(tqdm(bam_files, desc="Processing BAM files")):
            results[idx] = _run_flagstat(bam_file, samtools_binary)

    sample_names = []
    read_counts = []

    with raw_flagstat_path.open("w") as raw_out:
        for bam_file, sample_name, count, stdout in results:
            raw_out.write(f"{bam_file}\n")
            raw_out.write(stdout + "\n")
            sample_names.append(sample_name)
            read_counts.append(count)

    # Save summary table
    df = pd.DataFrame({"sample": sample_names, "num_seqs": read_counts})
    df.to_csv(summary_csv_path, index=False)

    print(f"Saved summary: {summary_csv_path}")
    print(f"Saved raw flagstat log: {raw_flagstat_path}")
    

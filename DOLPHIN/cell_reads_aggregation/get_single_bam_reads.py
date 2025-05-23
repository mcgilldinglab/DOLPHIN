import subprocess
from pathlib import Path
import pandas as pd
import os
from tqdm import tqdm

def run_reads_count(
    out_name: str,
    bam_file_path: str,
    out_directory: str = "./"
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

    # bam_files = Path(bam_file_path).rglob("*.bam")
    bam_files = list(Path(bam_file_path).rglob("*.bam"))

    sample_names = []
    read_counts = []

    with raw_flagstat_path.open("w") as raw_out:
        for bam_file in tqdm(bam_files, desc="Processing BAM files"):
            sample_name = bam_file.stem  # Remove .bam
            # print(f"Processing: {sample_name}")

            # Run samtools flagstat
            result = subprocess.run(
                ["samtools", "flagstat", str(bam_file)],
                capture_output=True,
                text=True
            )

            # Write to raw output log
            raw_out.write(f"{bam_file}\n")
            raw_out.write(result.stdout + "\n")

            # Extract read count (1st line of output, before ' + ' or space)
            first_line = result.stdout.splitlines()[0]
            count = int(first_line.split(" ")[0])  # More robust than assuming 14-line blocks

            sample_names.append(sample_name)
            read_counts.append(count)

    # Save summary table
    df = pd.DataFrame({"sample": sample_names, "num_seqs": read_counts})
    df.to_csv(summary_csv_path, index=False)

    print(f"Saved summary: {summary_csv_path}")
    print(f"Saved raw flagstat log: {raw_flagstat_path}")
    
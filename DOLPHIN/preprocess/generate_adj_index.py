import pickle
import pandas as pd
import os

def generate_adj_index_table(exon_pkl_path: str, output_dir: str = "./dolphin_exon_gtf/") -> pd.DataFrame:
    """
    Generate and save an adjacency index table for gene-level exon graphs from a exon pickle file.

    This function reads a `.pkl` file containing exon annotations (as a pandas DataFrame),
    groups exons by `gene_id`, calculates the number of exons per gene, and computes the
    flattened adjacency matrix indices for each gene using the formula:

        ind = exon_count^2
        ind_st = cumulative sum of previous `ind` values

    The resulting table is saved as `dolphin_adj_index.csv` in the specified output directory.

    Parameters
    ----------
    exon_pkl_path : str
        Path to the pickle file (.pkl) containing the exon DataFrame. The DataFrame must include a 'gene_id' column.

    output_dir : str, optional
        Directory where the output `dolphin_adj_index.csv` will be saved. Default is './dolphin_exon_gtf/'.

    Returns
    -------
    adj_df : pandas.DataFrame
        A DataFrame with the following columns:
        - 'geneid': gene ID
        - 'ind_st': starting index in the concatenated adjacency matrix
        - 'ind': size of the flattened square adjacency matrix for that gene (exon_count^2)

    Raises
    ------
    AssertionError
        If the gene order in the output does not match the input DataFrame's gene appearance order.
    """

    # 1. Load the DataFrame
    with open(exon_pkl_path, "rb") as f:
        exon_df = pickle.load(f)

    # 2. Preserve gene order as they first appear
    grouped = exon_df.groupby("gene_id", sort=False)
    gene_ordered = exon_df["gene_id"].drop_duplicates()
    exon_counts = grouped.size().reindex(gene_ordered).dropna()

    # 3. Build index table
    rows = []
    current_ind = 0

    for geneid, exon_count in exon_counts.items():
        ind = exon_count * exon_count
        rows.append({
            "geneid": geneid,
            "ind_st": float(current_ind),
            "ind": float(ind)
        })
        current_ind += ind

    # 4. Create result DataFrame
    adj_df = pd.DataFrame(rows)

    # 5. Check gene order consistency
    assert adj_df["geneid"].tolist() == gene_ordered.tolist(), "Gene order mismatch!"

    # 6. Save to CSV
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "dolphin_adj_index.csv")
    adj_df.to_csv(output_path, index=False)
    print(f"[Saved] Adjacency index table saved to: {output_path}")

    return adj_df

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


def generate_adj_metadata_table(exon_pkl_path: str, output_dir: str = "./dolphin_exon_gtf/") -> pd.DataFrame:
    """
    Generate metadata table for flattened exon adjacency matrices per gene.

    Ensures unique and non-missing gene names:
    - If gene_name is missing or empty, fallback to gene_id.
    - If gene_name is duplicated across gene_ids, disambiguate using gene_name-gene_id.

    Parameters
    ----------
    exon_pkl_path : str
        Path to exon .pkl file.
    output_dir : str, optional
        Output directory to save CSV file. Default is './dolphin_exon_gtf/'.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: 'Geneid', 'GeneName', 'Gene_Junc_name'
    """

    # 1. Load exon DataFrame
    with open(exon_pkl_path, "rb") as f:
        exon_df = pickle.load(f)

    # 2. Get gene order and exon counts
    gene_ordered = exon_df["gene_id"].drop_duplicates()
    grouped = exon_df.groupby("gene_id", sort=False)
    exon_counts = grouped.size().reindex(gene_ordered).dropna()

    # 3. Extract gene_id to gene_name map
    gene_info = exon_df.drop_duplicates("gene_id")[["gene_id", "gene_name"]].set_index("gene_id")

    # -- Fill missing gene names and empty strings with gene_id
    gene_info = gene_info.reset_index()
    gene_info["gene_name"] = gene_info.apply(
        lambda row: row["gene_id"] if pd.isna(row["gene_name"]) or row["gene_name"] == "" else row["gene_name"],
        axis=1
    )

    # -- Disambiguate duplicated gene_names
    name_counts = gene_info["gene_name"].value_counts()
    duplicated_names = name_counts[name_counts > 1].index

    gene_info["gene_name"] = gene_info.apply(
        lambda row: f"{row['gene_name']}-{row['gene_id']}" if row["gene_name"] in duplicated_names else row["gene_name"],
        axis=1
    )

    gene_info = gene_info.set_index("gene_id")

    # 4. Build metadata
    rows = []
    for geneid in gene_ordered:
        gene_name = gene_info.loc[geneid, "gene_name"]
        exon_count = exon_counts[geneid]
        n_junc = exon_count * exon_count

        for i in range(1, n_junc + 1):
            rows.append({
                "Geneid": geneid,
                "GeneName": gene_name,
                "Gene_Junc_name": f"{gene_name}-{i}"
            })

    # 5. Construct and save
    meta_df = pd.DataFrame(rows)
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "dolphin_adj_metadata_table.csv")
    meta_df.to_csv(output_path, index=False)
    print(f"[Saved] Metadata table saved to: {output_path}")

    return meta_df

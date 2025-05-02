import os
import pandas as pd
import numpy as np
import anndata
from tqdm import tqdm
from functools import reduce

def convert_psi_to_h5ad(
    outrigger_path: str,
    output_path: str,
    output_name: str,
    splice_type: str = "se"
):
    """
    Converts Outrigger PSI output into an AnnData object and saves it in `.h5ad` format.

    This function reads PSI values (percent spliced-in) from Outrigger's summary output,
    constructs a matrix of samples (cells) by splicing events, maps events to gene names,
    and returns a complete AnnData object ready for alternative splicing analysis.

    Parameters
    ----------
    outrigger_path : str
        Path to the directory containing Outrigger output files (e.g., "psi/", "index/").
    output_path : str
        Directory where the `.h5ad` file will be saved.
    output_name : str
        Filename (without extension) for the output AnnData file.
    splice_type : str, optional
        Type of splicing event to include (default is "se" for skipped exon).

    Returns
    -------
    adata : anndata.AnnData
        An AnnData object containing PSI values with event-gene annotation.

    Notes
    -----
    - PSI values are extracted from `outrigger_summary.csv`.
    - Only events of type `splice_type` (e.g., "se") are included.
    - Event-to-gene mapping uses both isoform gene names and gene IDs, with fallback logic.
    - Output file will be written to: `<output_path>/<output_name>.h5ad`.

    Example
    -------
    >>> convert_psi_to_h5ad(
            outrigger_path="./outrigger_output/",
            output_path="../data/",
            output_name="DOLPHIN_PSI_N10"
        )
    """

    # Load PSI summary and filter by splice type
    pd_psi = pd.read_csv(os.path.join(outrigger_path, "psi", "outrigger_summary.csv"))
    pd_psi = pd_psi[pd_psi["splice_type"] == splice_type]
    pd_psi["sample_id"] = pd_psi["sample_id"].apply(lambda x: x.split(".")[0])
    all_sample = set(pd_psi["sample_id"])

    # Convert to sample x event matrix
    d = {}
    for _srr in tqdm(all_sample, desc="Processing samples"):
        _temp_df = pd_psi[pd_psi["sample_id"] == _srr].rename(columns={"psi": _srr})[["event_id", _srr]]
        d[_srr] = _temp_df

    df_merged = reduce(lambda left, right: pd.merge(left, right, on="event_id", how="outer"), d.values())
    df_merged = df_merged.set_index("event_id")
    df_recon = df_merged.transpose()

    # Observation table (samples)
    obs = pd.DataFrame(df_recon.index).rename(columns={0: "CB"})
    obs = obs.set_index("CB", drop=False)
    obs.index.name = None

    # Load event metadata and extract gene names
    pd_event = pd.read_csv(os.path.join(outrigger_path, "index", splice_type, "events.csv"))
    pd_event["AS_event_type"] = splice_type.upper()
    pd_event["isoform1_gene_name_mod"] = pd_event["isoform1_gene_name"].fillna(pd_event["isoform1_gene_id"])
    pd_event["isoform2_gene_name_mod"] = pd_event["isoform2_gene_name"].fillna(pd_event["isoform2_gene_id"])

    pd_event_iso1 = pd_event[["event_id", "isoform1_gene_name_mod"]]
    pd_event_iso2 = pd_event[["event_id", "isoform2_gene_name_mod"]]

    pd_event_iso1_freq = pd_event_iso1.groupby(["event_id", "isoform1_gene_name_mod"]).size().reset_index(name="count1")
    pd_event_iso2_freq = pd_event_iso2.groupby(["event_id", "isoform2_gene_name_mod"]).size().reset_index(name="count2")

    pd_event_iso1_freq = pd_event_iso1_freq.sort_values(["event_id", "count1"], ascending=False).groupby("event_id").head(1)
    pd_event_iso2_freq = pd_event_iso2_freq.sort_values(["event_id", "count2"], ascending=False).groupby("event_id").head(1)

    pd_event_gene = pd.merge(pd_event_iso1_freq, pd_event_iso2_freq, on="event_id", how="outer")
    pd_event_gene["gene_name"] = np.select(
        [
            (pd_event_gene["isoform1_gene_name_mod"].notna()) & (pd_event_gene["isoform1_gene_name_mod"] == pd_event_gene["isoform2_gene_name_mod"]),
            (pd_event_gene["isoform1_gene_name_mod"].notna()) & (pd_event_gene["isoform2_gene_name_mod"].isna()),
            (pd_event_gene["isoform2_gene_name_mod"].notna()) & (pd_event_gene["isoform1_gene_name_mod"].isna()),
            (pd_event_gene["isoform1_gene_name_mod"].notna()) & (pd_event_gene["isoform2_gene_name_mod"].notna()) & (pd_event_gene["isoform1_gene_name_mod"] != pd_event_gene["isoform2_gene_name_mod"]),
            (pd_event_gene["isoform2_gene_name_mod"].isna()) & (pd_event_gene["isoform1_gene_name_mod"].isna())
        ],
        [
            pd_event_gene["isoform1_gene_name_mod"],
            pd_event_gene["isoform1_gene_name_mod"],
            pd_event_gene["isoform2_gene_name_mod"],
            pd_event_gene["isoform1_gene_name_mod"] + "," + pd_event_gene["isoform2_gene_name_mod"],
            "Empty"
        ]
    )
    pd_event_gene["gene_name"] = pd_event_gene["gene_name"].apply(lambda x: ",".join(sorted(set(x.split(",")))) if "," in x else x)
    gene_name_dict = dict(zip(pd_event_gene["event_id"], pd_event_gene["gene_name"]))

    # Variable table (events)
    var = pd.DataFrame(index=df_recon.columns)
    var["gene_name"] = var.index
    var.replace({"gene_name": gene_name_dict}, inplace=True)

    # Build AnnData
    adata = anndata.AnnData(X=df_recon.values, obs=obs, var=var, dtype=np.float32)
    adata.write(os.path.join(output_path, output_name + ".h5ad"))

    return adata

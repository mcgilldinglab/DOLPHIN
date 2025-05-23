import pandas as pd
import numpy as np
import anndata
import os
from tqdm import tqdm
from scipy.sparse import csr_matrix

def run_adjacency_matrix_final(
    out_name: str,
    out_directory: str = "./",
    batch_size=1000):
    """
    Generates the final adjacency matrix by filtering out invalid edges based on FeatureComp data.

    This function processes a compressed adjacency matrix and removes edges corresponding to exons 
    that are consistently zero across all cells.
 
    Parameters
    ----------
    out_name : str
        Output name prefix.
    out_directory : str, optional
        Output folder to save results.
        
    Returns
    -------
    None
        Saves the final adjacency matrix to the output directory as `AdjacencyCompRe_<out_name>.h5ad`.

    """

    print("Start Generating Final Adjacency Matrix...")
    
    final_out_dir = os.path.join(out_directory, "data")
    # temp_out_dir = os.path.join(final_out_dir, "temp")
    # os.makedirs(temp_out_dir, exist_ok=True)
    
    # === Step 1: Load adjacency and feature matrices
    adata_adj = anndata.read_h5ad(os.path.join(final_out_dir, f"AdjacencyComp_{out_name}.h5ad"))
    adata_adj.X = adata_adj.X.toarray()

    adata_fea_gtf = anndata.read_h5ad(os.path.join(final_out_dir, f"Feature_{out_name}.h5ad"))
    sample_list = list(adata_fea_gtf.obs_names)

    # === Step 2: Identify expressed exons
    df_fea_orig_var = pd.DataFrame(adata_fea_gtf.var)
    df_fea_gtf = pd.merge(adata_fea_gtf.to_df().T, df_fea_orig_var, how="left",left_index=True, right_index=True)
    df_fea_gtf["gene_order"] = df_fea_gtf["gene_id"].apply(lambda x: int(x[4:]))
    df_fea_gtf["gtf_index"] = df_fea_gtf.index
    df_fea_gtf["exon_order"] = df_fea_gtf["gtf_index"].apply(lambda x: int(x.split("-")[-1]))
    df_fea_gtf_order = df_fea_gtf.sort_values(by=["gene_order", "exon_order"])
    df_fea_gtf_order = df_fea_gtf_order.reset_index(drop=True)
    df_fea_gtf_order["new_gtf_exon_index"] = df_fea_gtf_order.groupby(["gene_id"]).cumcount()+1
    df_fea_gtf_order["new_gtf_index"] = df_fea_gtf_order['gene_id'].astype(str) +"-"+ df_fea_gtf_order["new_gtf_exon_index"].astype(str)
    df_fea_gtf_new_index = df_fea_gtf_order.set_index("new_gtf_index")
    df_fea_gtf_new_index = df_fea_gtf_new_index[sample_list]
    df_fea_gtf_new_index = df_fea_gtf_new_index.T

    df_fea_exon_keep = df_fea_gtf_new_index.loc[:, (df_fea_gtf_new_index > 0).any()]
    expressed_exons = df_fea_exon_keep.columns
    
    exon_keep_dict = {
        gene: set(exon.split("-")[-1] for exon in expressed_exons if exon.startswith(gene))
        for gene in set(x.rsplit("-", 1)[0] for x in expressed_exons)
    }

    # === Step 3: Prepare adjacency DataFrame
    df_adj = adata_adj.to_df().T
    df_adj["_ck_empty"] = np.amax(df_adj.to_numpy(), axis=1)
    df_adj["_idx"] = df_adj.index

    var_df = pd.DataFrame(adata_adj.var)
    gene_name_to_id = dict(zip(var_df["gene_name"], var_df["gene_id"]))

    df_adj["gene_name"] = df_adj["_idx"].apply(lambda x: x[:x.rfind("-")])
    df_adj["gene_id"] = df_adj["gene_name"].map(gene_name_to_id)
    df_adj["adj_vec_size"] = df_adj.groupby("gene_id")["gene_id"].transform("count")
    df_adj["gene_order"] = df_adj["gene_id"].apply(lambda x: int(x[4:]))
    df_adj["exon_order"] = df_adj["_idx"].apply(lambda x: int(x.split("-")[-1]))
    df_adj["exon_name"] = df_adj["gene_id"].astype(str) + "-" + df_adj["exon_order"].astype(str)

    # === Step 4: Flag zero edges
    def flag_zero(group):
        gene_id = group["gene_id"].iloc[0]
        adj_size = int(group["adj_vec_size"].iloc[0])
        n = int(adj_size ** 0.5)
        matrix = group["_ck_empty"].to_numpy().reshape((n, n))
        mask = np.zeros((n, n), dtype=np.uint8)

        keep_ids = exon_keep_dict.get(gene_id, set())

        for i in range(n):
            if np.all(matrix[i, :] == 0) and np.all(matrix[:, i] == 0) and str(i + 1) not in keep_ids:
                mask[i, :] = 1
                mask[:, i] = 1
        return pd.Series(mask.flatten())

    # === Step 5: Process in gene_id batches
    all_flags = []
    ordered_gene_ids = df_adj["gene_id"].drop_duplicates().tolist()
    gene_batches = [ordered_gene_ids[i:i + batch_size] for i in range(0, len(ordered_gene_ids), batch_size)]

    for batch in tqdm(gene_batches, desc="Processing gene batches"):
        temp_df = df_adj[df_adj["gene_id"].isin(batch)].copy()
        temp_df = temp_df.sort_index()

        flags = temp_df.groupby("gene_id").apply(flag_zero).explode().astype(np.uint8).reset_index(level=0)
        flags.columns = ["gene_id", "_flag"]
        flags["exon_idx"] = flags.groupby("gene_id").cumcount() + 1
        flags["exon_name"] = flags["gene_id"] + "-" + flags["exon_idx"].astype(str)
        all_flags.append(flags)

    df_flags = pd.concat(all_flags, ignore_index=True)

    # === Step 6: Merge flags with original matrix
    df_flags = df_flags.set_index("exon_name").reindex(df_adj["exon_name"]).reset_index()
    df_adj["_flag"] = df_flags["_flag"].values

    # === Step 7: Filter flagged entries and save
    df_filtered = df_adj[df_adj["_flag"] == 0].copy()
    
    # === Step 8: Convert to h5ad
    df_filtered["new_index"] = df_filtered.groupby(["gene_id"]).cumcount()+1
    df_filtered["var_new_index"] = df_filtered['gene_id'].astype(str) +"-"+ df_filtered["new_index"].astype(str)
    df_filtered = df_filtered.reset_index(drop=True)

    ## == step6: adata for adjacency matrix
    adata_adj_comp = anndata.read_h5ad(os.path.join(final_out_dir, "AdjacencyComp_"+out_name+".h5ad"))

    obs = pd.DataFrame(adata_adj_comp.obs)
    ## dataframe for annotating the variables = geneid
    var_names = df_filtered["var_new_index"].values
    var = pd.DataFrame(index=var_names)
    var["gene_id"] = df_filtered["gene_id"].values
    var["gene_name"] = df_filtered["gene_name"].values

    # # ##the data matrix 
    X = df_filtered.drop(columns=["gene_id", "gene_name","_ck_empty","_idx","adj_vec_size","_flag", "gene_order", "exon_order","exon_name","new_index","var_new_index"]).T.iloc[:,:].values
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    adata.X = csr_matrix(adata.X)
    
    adata.write(os.path.join(final_out_dir, "AdjacencyCompRe_" + out_name +".h5ad"))

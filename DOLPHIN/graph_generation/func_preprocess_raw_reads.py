"""
This function generates standardized feature and adjacency matrices for all cells, 
to be used as input for the GVAE (Graph Variational Autoencoder) model.

The structure and dimensions of both matrices are determined entirely by the reference GTF file. 
No filtering is applied — all genes and exons listed in the GTF are included.

Feature Matrix:
- Each exon contributes a 1×1 scalar feature (e.g., raw count).
- Each gene forms an Ni×1 feature vector, where Ni is the number of exons for gene i.
- The full feature vector for a single cell is of size (N1 + N2 + ... + Nn) × 1, 
  where n is the total number of genes in the GTF.

Adjacency Matrix:
- Each gene contributes an Ni×Ni adjacency submatrix, where Ni is the number of exons for gene i.
- The full adjacency matrix for a single cell is a flattened vector of size (N1² + N2² + ... + Nn²).

To ensure compatibility across all cells:
- Genes and exons absent in a given cell are zero-filled.
- This guarantees that all output matrices have the same dimensions, regardless of actual gene/exon coverage.

Note:
- No normalization is applied to either the feature matrix or the adjacency matrix.
"""

import pandas as pd
import numpy as np
import os

from .grouped_featurecounts import (
    extract_grouped_cell_exon_counts,
    extract_grouped_cell_junction_counts,
)


GROUPED_EXON_COUNTS = None
GROUPED_JUNCTION_COUNTS = None
EXON_POSITION_COLUMN = "__dolphin_exon_pos"


def _locate_junction_positions(exon_start, exon_end, positions):
    """Map genomic splice-site coordinates to exon slots."""
    if positions.size == 0:
        return np.array([], dtype=np.float64)

    located = np.full(positions.shape, np.nan, dtype=np.float64)
    before_first = positions < exon_start[0]
    located[before_first] = -1.0

    candidate = np.searchsorted(exon_start, positions, side="right") - 1
    valid_candidate = candidate >= 0
    safe_candidate = np.clip(candidate, 0, len(exon_start) - 1)
    next_candidate = np.clip(safe_candidate + 1, 0, len(exon_start) - 1)

    inside_exon = valid_candidate & (positions <= exon_end[safe_candidate])
    located[inside_exon] = safe_candidate[inside_exon].astype(np.float64)

    between_exons = (
        np.isnan(located)
        & valid_candidate
        & (candidate < len(exon_start) - 1)
        & (positions > exon_end[safe_candidate])
        & (positions < exon_start[next_candidate])
    )
    located[between_exons] = safe_candidate[between_exons].astype(np.float64) + 0.5

    after_last = (
        np.isnan(located)
        & (safe_candidate == len(exon_start) - 1)
        & (positions > exon_end[-1])
    )
    located[after_last] = safe_candidate[after_last].astype(np.float64) + 0.5
    return located


def build_raw_adjacency_flat(
    exon_start,
    exon_end,
    exon_count,
    junction_start,
    junction_end,
    junction_weight,
):
    """Build the raw junction adjacency matrix for one gene."""
    exon_start = np.asarray(exon_start, dtype=np.int64)
    exon_end = np.asarray(exon_end, dtype=np.int64)
    exon_count = np.asarray(exon_count, dtype=np.float64).reshape(-1)
    junction_start = np.asarray(junction_start, dtype=np.int64)
    junction_end = np.asarray(junction_end, dtype=np.int64)
    junction_weight = np.asarray(junction_weight, dtype=np.float64)

    n_node = exon_start.size
    adjacency = np.zeros((n_node, n_node), dtype=np.float64)
    if n_node == 0 or junction_start.size == 0:
        return adjacency.ravel()
    if exon_end.size != n_node or exon_count.size != n_node:
        raise ValueError("Exon starts, ends, and counts must have the same length.")
    if not (junction_start.size == junction_end.size == junction_weight.size):
        raise ValueError("Junction starts, ends, and weights must have the same length.")

    start_slots = _locate_junction_positions(exon_start, exon_end, junction_start)
    end_slots = _locate_junction_positions(exon_start, exon_end, junction_end)

    start_half = np.isfinite(start_slots) & np.isclose(np.mod(start_slots, 1.0), 0.5)
    end_half = np.isfinite(end_slots) & np.isclose(np.mod(end_slots, 1.0), 0.5)
    start_slots[start_half] -= 0.5
    end_slots[end_half] += 0.5
    end_slots[end_slots >= n_node] = n_node - 1

    finite = np.isfinite(start_slots) & np.isfinite(end_slots)
    start_idx = np.zeros(start_slots.size, dtype=np.int64)
    end_idx = np.zeros(end_slots.size, dtype=np.int64)
    start_idx[finite] = start_slots[finite].astype(np.int64)
    end_idx[finite] = end_slots[finite].astype(np.int64)
    valid = (
        finite
        & (start_idx >= 0)
        & (end_idx >= 0)
        & (start_idx < n_node)
        & (end_idx < n_node)
        & (start_idx < end_idx)
        & (exon_count[start_idx] > 0)
        & (exon_count[end_idx] > 0)
        & np.isfinite(junction_weight)
        & (junction_weight > 0)
    )
    np.add.at(adjacency, (start_idx[valid], end_idx[valid]), junction_weight[valid])
    return adjacency.ravel()


def configure_grouped_count_sources(exon_counts=None, junction_counts=None):
    global GROUPED_EXON_COUNTS, GROUPED_JUNCTION_COUNTS
    GROUPED_EXON_COUNTS = exon_counts
    GROUPED_JUNCTION_COUNTS = junction_counts


def build_reference_lookups(gtf, adj_ind):
    gene_index_map = {
        gene_id: gene_indices.to_list()
        for gene_id, gene_indices in gtf.groupby("Geneid").groups.items()
    }
    adj_start_by_gene = (
        adj_ind[["geneid", "ind_st"]]
        .drop_duplicates(subset=["geneid"])
        .set_index("geneid")["ind_st"]
        .astype(int)
        .to_dict()
    )
    return gene_index_map, adj_start_by_gene

"""
Load all the necessary files
"""
def get_gtf(path_gtf, path_adj_index):
    full_gtf = pd.read_pickle(path_gtf)
    full_gtf = full_gtf[["seqname","gene_id","start","end"]]
    full_gtf = full_gtf.rename(columns={"seqname":"Chr","gene_id":"Geneid","start":"Start","end":"End"})
    full_gtf["Start"] = full_gtf["Start"].astype(int)
    full_gtf["End"] = full_gtf["End"].astype(int)

    #get cell level adj matrix, generate the index position table for each gene adjacency matrix
    df_adj_ind = pd.read_csv(path_adj_index, skiprows = 0, sep = ',',low_memory=False,index_col=False)

    return(full_gtf, df_adj_ind)

"""
gene class
"""
class gene(object):
    def __init__(
        self,
        gtf,
        adj_ind,
        srr,
        id = "",
        show = "N",
        ck = "N",
        main_path="./",
        gene_index_map=None,
        adj_start_by_gene=None,
        verbose=False,
    ):
        self.srr = srr
        self.id = id #if id == "", then will run the entire cell
        self.show = show #if show == "Y", will show all before and after adj table, show visualization,default is N
        self.ck = ck #if ck == "Y", will check the adjacency table
        self.exon_table = pd.DataFrame()
        self.junct_table = pd.DataFrame()
        self.feat_mat = []
        self.adj_mat = []
        self.gtf = gtf
        self.adj_ind = adj_ind
        self.main_path = main_path
        self._cell_exon_counts = None
        self._cell_junction_counts = None
        self._cell_exon_by_gene = None
        self._cell_junction_by_gene = None
        self._empty_exon_counts = None
        self._empty_junction_counts = None
        self._gene_exon_template_cache = {}
        self._current_exon_start = None
        self._current_exon_end = None
        self._current_exon_count = None
        self.gene_index_map = gene_index_map or {}
        self.adj_start_by_gene = adj_start_by_gene or {}
        self.verbose = bool(verbose)

    def _using_grouped_counts(self):
        return GROUPED_EXON_COUNTS is not None and GROUPED_JUNCTION_COUNTS is not None

    def _load_cell_exon_counts(self):
        if self._cell_exon_counts is not None:
            return self._cell_exon_counts

        if self._using_grouped_counts():
            df_temp_exon = extract_grouped_cell_exon_counts(
                grouped_exon_counts=GROUPED_EXON_COUNTS,
                cell_id=self.srr,
            )
        else:
            df_temp_exon = pd.read_csv(
                os.path.join(self.main_path, "05_exon_junct_cnt", self.srr + ".exon.count.txt"),
                skiprows=1,
                sep='\t',
                low_memory=False,
            )
            df_temp_exon.columns.values[6] = "Count"
            df_temp_exon = df_temp_exon[df_temp_exon["Count"] > 0]
            df_temp_exon["Type"] = "Exon"

        self._cell_exon_counts = df_temp_exon
        return self._cell_exon_counts

    def _load_cell_exon_by_gene(self):
        if self._cell_exon_by_gene is not None:
            return self._cell_exon_by_gene

        df_temp_exon = self._load_cell_exon_counts()
        self._empty_exon_counts = df_temp_exon.iloc[0:0].copy()
        if not df_temp_exon.empty:
            df_temp_exon = (
                df_temp_exon.groupby(
                    ["Geneid", "Chr", "Start", "End", "Type"],
                    sort=False,
                    as_index=False,
                )["Count"]
                .sum()
            )
        self._cell_exon_by_gene = {
            gene_id: self._annotate_exon_positions(gene_id, dataframe)
            for gene_id, dataframe in df_temp_exon.groupby("Geneid", sort=False)
        }
        return self._cell_exon_by_gene

    def _annotate_exon_positions(self, gene_id, dataframe):
        if dataframe.empty:
            return dataframe

        original_gene_id = self.id
        self.id = gene_id
        _, _, position_map, _, _ = self._get_gene_exon_template()
        self.id = original_gene_id

        annotated = dataframe.copy()
        annotated[EXON_POSITION_COLUMN] = [
            position_map[(row.Chr, row.Geneid, row.Start, row.End)]
            for row in annotated.itertuples(index=False)
        ]
        return annotated

    def _load_cell_junction_counts(self):
        if self._cell_junction_counts is not None:
            return self._cell_junction_counts

        if self._using_grouped_counts():
            df_temp_sj = extract_grouped_cell_junction_counts(
                grouped_junction_counts=GROUPED_JUNCTION_COUNTS,
                cell_id=self.srr,
            )
        else:
            df_temp_sj = pd.read_csv(
                os.path.join(self.main_path, "05_exon_junct_cnt", self.srr + ".exon.count.txt.jcounts"),
                skiprows=0,
                sep='\t',
                low_memory=False,
            )
            df_temp_sj = extract_grouped_cell_junction_counts(
                grouped_junction_counts=df_temp_sj,
                cell_id=df_temp_sj.columns[-1],
            )

        self._cell_junction_counts = df_temp_sj
        return self._cell_junction_counts

    def _load_cell_junction_by_gene(self):
        if self._cell_junction_by_gene is not None:
            return self._cell_junction_by_gene

        df_temp_sj = self._load_cell_junction_counts()
        self._empty_junction_counts = df_temp_sj.iloc[0:0].copy()
        self._cell_junction_by_gene = {
            gene_id: dataframe
            for gene_id, dataframe in df_temp_sj.groupby("Geneid", sort=False)
        }
        return self._cell_junction_by_gene

    def _get_gene_exon_template(self):
        cached_template = self._gene_exon_template_cache.get(self.id)
        if cached_template is not None:
            return cached_template

        gene_indices = self.gene_index_map.get(self.id)
        if gene_indices is None:
            template = self.gtf[(self.gtf["Geneid"] == self.id)][
                ["Chr", "Geneid", "Start", "End"]
            ].copy()
        else:
            template = self.gtf.loc[gene_indices, ["Chr", "Geneid", "Start", "End"]].copy()
        template["Count"] = 0.0
        template["Type"] = "Exon"
        template_index = pd.MultiIndex.from_frame(
            template[["Chr", "Geneid", "Start", "End"]]
        )
        position_map = {
            exon_key: index
            for index, exon_key in enumerate(template_index.tolist())
        }
        exon_start = template["Start"].to_numpy(dtype=np.int64, copy=False)
        exon_end = template["End"].to_numpy(dtype=np.int64, copy=False)
        cached_template = (
            template,
            template_index,
            position_map,
            exon_start,
            exon_end,
        )
        self._gene_exon_template_cache[self.id] = cached_template
        return cached_template

    """
    Get the gene list for each individual cell
    """
    def get_gene(self):
        exon_by_gene = self._load_cell_exon_by_gene()
        #convert gene to sorted list, only run exist gene per cell, not run all genes in gtf
        self.gene_list = sorted(exon_by_gene)

    """
    Read the exon count table and only keep count number > 0
    """
    def exon_read(self):
        exon_by_gene = self._load_cell_exon_by_gene()
        if self._empty_exon_counts is None:
            self._empty_exon_counts = self._load_cell_exon_counts().iloc[0:0].copy()
        df_exon = exon_by_gene.get(self.id, self._empty_exon_counts)
        (
            template,
            _template_index,
            _position_map,
            exon_start,
            exon_end,
        ) = self._get_gene_exon_template()
        exon_counts = np.zeros(len(template), dtype=np.float64)
        if not df_exon.empty:
            exon_positions = df_exon[EXON_POSITION_COLUMN].to_numpy(
                dtype=np.int64,
                copy=False,
            )
            exon_values = df_exon["Count"].to_numpy(dtype=np.float64, copy=False)
            np.add.at(exon_counts, exon_positions, exon_values)
        self._current_exon_start = exon_start
        self._current_exon_end = exon_end
        self._current_exon_count = exon_counts
        if (self.ck.upper() == "Y") or (self.show.upper() == "Y"):
            exon_table = template.copy()
            exon_table["Count"] = exon_counts
            self.exon_table = exon_table
        else:
            self.exon_table = pd.DataFrame()
    
    """
    Read the junction count table and clean
    """
    def jun_read(self):
        junction_by_gene = self._load_cell_junction_by_gene()
        if self._empty_junction_counts is None:
            self._empty_junction_counts = (
                self._load_cell_junction_counts().iloc[0:0].copy()
            )
        self.junct_table = junction_by_gene.get(self.id, self._empty_junction_counts)
    
    """
    Combine exon and junction table in order to check
    """
    def show_full_table(self):
        self.exon_read()
        self.jun_read()
        df_comb = pd.concat([self.exon_table, self.junct_table], ignore_index=True)
        df_comb = df_comb.sort_values(by=["Chr","Geneid","Start"],ascending=[True,True,True]).reset_index()
        print("Combined exon and junction table:")
        display(df_comb.head(200))

    """
    Define FEATURE MATRIX, one exon has feature matrix N * 1, one gene has N * N.
    """
    def count_fea(self):
        #node feature
        if self._current_exon_count is not None:
            node_feature = self._current_exon_count
        else:
            node_feature = (
                self.exon_table.sort_values("Start")["Count"]
                .to_numpy(dtype=np.float64, copy=False)
            )
        self.feat_mat = node_feature.reshape(-1, 1)

    """
    Define ADJACENCY MATRIX, one exon has adj matrix N * N
    """
    def count_adj(self):
        if self._current_exon_start is not None and self._current_exon_end is not None:
            exon_start = self._current_exon_start
            exon_end = self._current_exon_end
        else:
            exon_start = self.exon_table["Start"].to_numpy(dtype=np.int64, copy=False)
            exon_end = self.exon_table["End"].to_numpy(dtype=np.int64, copy=False)

        junction_start = self.junct_table["Start"].to_numpy(dtype=np.int64, copy=False)
        junction_end = self.junct_table["End"].to_numpy(dtype=np.int64, copy=False)
        junction_weight = self.junct_table["Count"].to_numpy(dtype=np.float64, copy=False)
        exon_count = np.asarray(self.feat_mat, dtype=np.float64).reshape(-1)
        self.adj_mat = build_raw_adjacency_flat(
            exon_start=exon_start,
            exon_end=exon_end,
            exon_count=exon_count,
            junction_start=junction_start,
            junction_end=junction_end,
            junction_weight=junction_weight,
        )

        if (self.ck.upper() == "Y") or (self.show.upper() == "Y"):
            start_slots = _locate_junction_positions(exon_start, exon_end, junction_start)
            end_slots = _locate_junction_positions(exon_start, exon_end, junction_end)
            original_edges = pd.DataFrame(
                {
                    "edge_orig": np.column_stack([start_slots, end_slots]).tolist(),
                    "weight": junction_weight,
                }
            )
            adjacency = self.adj_mat.reshape(len(exon_start), len(exon_start))
            kept_start, kept_end = np.nonzero(adjacency)
            kept_edges = pd.DataFrame(
                {
                    "start": kept_start,
                    "end": kept_end,
                    "weight": adjacency[kept_start, kept_end],
                }
            )
            print("Original junction mapping:")
            display(original_edges)
            print("Raw adjacency edges:")
            display(kept_edges)
            print("Number of Node:", len(exon_start))
            print("Number of Edge:", kept_edges.shape[0])
            print("ADJ MATRIX:", self.adj_mat)

    """
    Show ADJACENCY MATRIX
    """
    #initialize direct graph
    def adj_show(self):
        from matplotlib import pyplot as plt
        import networkx as nx

        #initialize direct graph
        G = nx.DiGraph()

        #add node
        for i in range(0,len(self.feat_mat)):
            G.add_node(i)
                
        #add junction
        for i in range(0,len(self.feat_mat)):
            for j in range(0,len(self.adj_mat)):
                if self.adj_mat[i][j] != 0:
                    G.add_weighted_edges_from([(i,j,self.adj_mat[i][j])])
        print("Node list:", G.nodes()) # returns a list
        print("Edge List:", G.edges()) # returns a list

        pos=nx.spring_layout(G)
        nx.draw(G, pos, with_labels=True, font_weight="bold")
        edge_weight = nx.get_edge_attributes(G,"weight")
        nx.draw_networkx_edge_labels(G, pos, edge_labels = edge_weight,label_pos=0.5)
        plt.show()
    
    """
    Run all the functions per gene/cell to update attributes
    """
    def get_all(self):
        if self.id == "":
            #get gene list per cell
            self.get_gene()
            #initialize feature matrix
            vec_f = np.zeros(shape=(len(self.gtf)))
            #initialize adjacency matrix
            adj_size = int(self.adj_ind["ind"].sum())
            vec_a = np.zeros(shape=(adj_size))
            for i in range(0,len(self.gene_list)):
            # for i in tqdm(range(len(self.gene_list)), desc=f"[{self.srr}] Processing genes", leave=False):
                self.id = self.gene_list[i]
                if self.verbose:
                    print("Sample = ",self.srr, ", Gene id = ",self.id, "is running.")
                self.exon_read()
                self.jun_read()
                self.count_fea()
                #combine all feature matrix of genes per cell 
                temp_indx = self.gene_index_map.get(self.id)
                if temp_indx is None:
                    temp_indx = self.gtf[(self.gtf['Geneid'] == self.id)].index.to_list()
                vec_f[np.asarray(temp_indx, dtype=int)] = self.feat_mat[:, 0]
                
                self.count_adj()
                #combine all adjacency matrix of genes per cell 
                start_indx = self.adj_start_by_gene.get(self.id)
                if start_indx is None:
                    start_indx = int(self.adj_ind[(self.adj_ind["geneid"] == self.id)]['ind_st'].values[0])
                vec_a[int(start_indx):int(start_indx) + len(self.adj_mat)] = self.adj_mat

            #save the final data to csv file
            np.savetxt(os.path.join(self.main_path, "06_graph_mtx", self.srr+ "_fea.csv"), vec_f, fmt='%10.4f', delimiter = ',')
            np.savetxt(os.path.join(self.main_path, "06_graph_mtx", self.srr+ "_adj.csv"), vec_a, fmt='%10.4f', delimiter = ',')
            print("========================Sample",self.srr,"is Done========================")
            return(vec_f,vec_a)
        else:
            #if only running one gene per sample, no need to produce tensor output, matrix is enough for debug
            if self.verbose:
                print("Sample = ",self.srr, ", Gene id = ",self.id, "is running.")
            self.exon_read()
            self.jun_read()
            self.count_fea()
            self.count_adj()
            if (self.show.upper() == "Y"):
                self.show_full_table()
                self.adj_show()
            print(self.feat_mat)
            print("Done")

    def get_all_sparse(self):
        if self.id != "":
            raise ValueError("get_all_sparse only supports whole-cell processing.")

        self.get_gene()
        feature_indices = []
        feature_values = []
        adjacency_indices = []
        adjacency_values = []
        feature_size = len(self.gtf)
        adjacency_size = int(self.adj_ind["ind"].sum())

        for gene_id in self.gene_list:
            self.id = gene_id
            if self.verbose:
                print("Sample = ", self.srr, ", Gene id = ", self.id, "is running.")

            self.exon_read()
            self.jun_read()
            self.count_fea()

            temp_indx = self.gene_index_map.get(self.id)
            if temp_indx is None:
                temp_indx = self.gtf[(self.gtf["Geneid"] == self.id)].index.to_list()
            temp_indx = np.asarray(temp_indx, dtype=np.int64)
            feat_vals = self.feat_mat[:, 0].astype(np.float32, copy=False)
            feat_nonzero = np.flatnonzero(feat_vals)
            if feat_nonzero.size:
                feature_indices.append(temp_indx[feat_nonzero].astype(np.int32, copy=False))
                feature_values.append(feat_vals[feat_nonzero])

            self.count_adj()
            start_indx = self.adj_start_by_gene.get(self.id)
            if start_indx is None:
                start_indx = int(
                    self.adj_ind[(self.adj_ind["geneid"] == self.id)]["ind_st"].values[0]
                )
            adj_vals = self.adj_mat.astype(np.float32, copy=False)
            adj_nonzero = np.flatnonzero(adj_vals)
            if adj_nonzero.size:
                adjacency_indices.append(
                    (adj_nonzero + int(start_indx)).astype(np.int32, copy=False)
                )
                adjacency_values.append(adj_vals[adj_nonzero])

        if feature_indices:
            feature_indices = np.concatenate(feature_indices)
            feature_values = np.concatenate(feature_values)
        else:
            feature_indices = np.empty(0, dtype=np.int32)
            feature_values = np.empty(0, dtype=np.float32)

        if adjacency_indices:
            adjacency_indices = np.concatenate(adjacency_indices)
            adjacency_values = np.concatenate(adjacency_values)
        else:
            adjacency_indices = np.empty(0, dtype=np.int32)
            adjacency_values = np.empty(0, dtype=np.float32)

        return {
            "sample_id": self.srr,
            "feature_size": int(feature_size),
            "feature_indices": feature_indices,
            "feature_values": feature_values,
            "adjacency_size": int(adjacency_size),
            "adjacency_indices": adjacency_indices,
            "adjacency_values": adjacency_values,
        }

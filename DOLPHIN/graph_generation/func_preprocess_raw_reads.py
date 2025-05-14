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
from matplotlib import pyplot as plt
import networkx as nx
import os
from tqdm import tqdm

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
    def __init__(self, gtf, adj_ind, srr, id = "", show = "N", ck = "N", main_path="./"):
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

    """
    Get the gene list for each individual cell
    """
    def get_gene(self):
        df_temp_exon = pd.read_csv(os.path.join(self.main_path, "05_exon_junct_cnt", self.srr + ".exon.count.txt"), skiprows = 1, sep = '\t',low_memory=False)
        df_temp_exon.columns.values[6] = "Count"
        df_temp_exon = df_temp_exon[df_temp_exon["Count"] > 0]
        df_count_gene = df_temp_exon.groupby(['Geneid'])['Geneid'].count().reset_index(name='count')
        #convert gene to list, only run exist gene per cell, not run all genes in gtf
        self.gene_list = df_count_gene['Geneid'].tolist()

    """
    Read the exon count table and only keep count number > 0
    """
    def exon_read(self):
        df_temp_exon = pd.read_csv(os.path.join(self.main_path, "05_exon_junct_cnt", self.srr + ".exon.count.txt"), skiprows = 1, sep = '\t',low_memory=False)
        df_temp_exon.columns.values[6] = "Count"
        df_temp_exon = df_temp_exon[df_temp_exon["Count"] > 0]
        df_temp_exon["Type"] = "Exon"
        df_exon = df_temp_exon[(df_temp_exon['Geneid']==self.id)]
        #modify exon if some exon is missing from annotation file, to make sure all the cell have the same amout of exons
        df_gtf = self.gtf[(self.gtf["Geneid"] == self.id)]
        self.exon_table = pd.merge(df_gtf,df_exon,on=["Chr","Geneid", "Start", "End"], how="outer")
        self.exon_table["Count"] = self.exon_table["Count"].fillna(0)
        self.exon_table["Type"] = self.exon_table["Type"].fillna("Exon")
    
    """
    Read the junction count table and clean
    """
    def jun_read(self):
        df_temp_sj = pd.read_csv(os.path.join(self.main_path, "05_exon_junct_cnt", self.srr + ".exon.count.txt.jcounts"), skiprows = 0, sep = '\t')
        df_temp_sj.columns.values[8] = "Count"

        #remove junction if primary gene is NaN
        df_temp_sj = df_temp_sj.dropna(subset=['PrimaryGene'])
        df_temp_sj = df_temp_sj[["PrimaryGene","Site1_chr","Site1_location","Site2_location","Count"]]
        df_temp_sj = df_temp_sj.rename(columns={"PrimaryGene":"Geneid","Site1_chr":"Chr","Site1_location":"Start","Site2_location":"End"})
        df_temp_sj["Type"] = "Junction"
        self.junct_table = df_temp_sj[(df_temp_sj['Geneid']==self.id)]
    
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
        #define number of node = number of exon
        n_node = len(self.exon_table)
        
        #node feature
        node_feature = self.exon_table.sort_values('Start').reset_index()
        node_feature.drop(["Geneid","Chr","Start","End","index","Strand","Length","Type"],axis=1,inplace=True)

        #covert to matrix
        self.feat_mat = np.zeros(shape=(n_node,1))
        for i in range(0,n_node):
            self.feat_mat[i][0] = node_feature["Count"][i]

    """
    Define ADJACENCY MATRIX, one exon has adj matrix N * N
    """
    def count_adj(self):
        n_node = len(self.exon_table)
        #1. initialize all list
        _exon_list = [] #original exon locations
        _jun_list = [] #original junction locations
        _wgt_list = [] #edge weight correspoding to jucntion read

        for i in range(self.exon_table.shape[0]):
            _exon_list.append([self.exon_table.iloc[i]['Start'],self.exon_table.iloc[i]['End']])
        for j in range(self.junct_table.shape[0]):
            _jun_list.append([self.junct_table.iloc[j]['Start'],self.junct_table.iloc[j]['End']])
            _wgt_list.append(self.junct_table.iloc[j]['Count'])

        _out_list = [[np.nan,np.nan] for _ in range(len(_jun_list))] #output adjacent list

        #label each junction about which slot they are in,len(_out_list) = len(_jun_list)
        #start and end check seperately
        for k in range(len(_jun_list)):
            for l in range(n_node):
                for m in (0,1): #0 for start, 1 for end:
                    #junction falls before first exon
                    if _jun_list[k][m] < _exon_list[0][0]:
                        _out_list[k][m] = -1
                    #junction falls into exon region:
                    elif (_exon_list[l][0] <= _jun_list[k][m] <= _exon_list[l][1]):
                        _out_list[k][m] = l
                    #junction falls into exon-exon region:
                    elif l < n_node-2: 
                        if (_exon_list[l][1] < _jun_list[k][m] < _exon_list[l+1][0]):
                            _out_list[k][m] = l+0.5
                    #junction greater then the last exon end:
                    elif (l == n_node-1) & (_exon_list[l][1] < _jun_list[k][m]):
                        _out_list[k][m] = l+0.5
        
        #postprocess of edge dataset
        _start = [i[0] for i in _out_list]
        _end = [i[1] for i in _out_list]
        _df_adj = pd.DataFrame(list(zip(_out_list,_start,_end,_wgt_list)),columns = ["edge_orig","start","end","weight"])
        _df_adj_orig = _df_adj
        _df_adj["edge_orig"] = _df_adj["edge_orig"].astype("str")
        _df_adj["_status"] = ""
        _df_adj["edge_mod"] = ""
        #only modify edge which may influency matrix
        for k in range(_df_adj.shape[0]):
            #1.entire junction before or within first exon- delete
            if _df_adj['start'][k] == _df_adj['end'][k] == -1:
                _df_adj.loc[k,'_status'] = "D"
            #2.junction after last exon - delete
            elif (_df_adj['start'][k] == _df_adj['end'][k] > n_node-0.5):
                _df_adj.loc[k,'_status'] = "D"
            #3.junction falls into beween exon.
            else:
                if _df_adj['start'][k] % 1 == 0.5:
                    _df_adj.loc[k, 'start'] = _df_adj['start'][k] - 0.5
                if _df_adj['end'][k] % 1 == 0.5:
                    #after the last exon will be end at last exon
                    if int(_df_adj['end'][k] +0.5) >= n_node:
                        _df_adj.loc[k, 'end'] = n_node-1
                    else:  
                    #else be end at the next exon 
                        _df_adj.loc[k, 'end'] = _df_adj['end'][k] + 0.5
             #4.junction with exon, will keep since not influcen adjency matrix
            _df_adj.loc[k,"edge_mod"] = "["+str(_df_adj["start"][k])+", "+str(_df_adj["end"][k])+"]"

        #add dummy junction if no junction data in between, with weight equal -99(1)!!!!!!!!!!
        _df_ck_adj = _df_adj[["start","end"]]

        #if more than 1 node but junction table is empty, need to add dummy junction between nodes
        if n_node > 1 and _df_adj.empty:
            for j in range(1,n_node+1):
                _df_adj = pd.concat([_df_adj, pd.DataFrame.from_records([{"edge_mod":"["+str(j-1)+", "+str(j)+"]","start":j-1, "end":j,"weight":1, "_status": "A" }])], ignore_index= True)
            _df_adj_sum = _df_adj
        #delete any edges if there is only one node
        elif(n_node == 1):
            _df_adj_sum = pd.DataFrame()
        else:
            #add dummy junction and assign weight with 1, 
            """
            if the feature matrix has count = 0 at this node, which cannot assign to 1, should be zero, 
            right now, this  patch is fixed by program "compact_matrix.ipynb" <= THIS NEED TO BE FIXED HERE THOUGH
            """
            for j in range(1,n_node):
                if not ([True,  True] in np.equal.outer(_df_ck_adj.to_numpy(copy=False),  [j-1,j]).any(axis=1).tolist()):
                    _df_adj = pd.concat([_df_adj, pd.DataFrame.from_records([{"edge_mod":"["+str(j-1)+", "+str(j)+"]","start":j-1, "end":j,"weight":1, "_status": "A" }])], ignore_index= True)
            #summary dataframe by edge
            _df_adj_sum = _df_adj.groupby( ['edge_mod','start','end','_status'], as_index=False).sum(numeric_only=True)

            #if A exist in the dataframe, add 1 on other junctions
            if "A" in _df_adj_sum['_status'].values:
                _df_adj_sum.loc[_df_adj_sum['_status'] != "A", "weight"] = _df_adj_sum["weight"] + 1
            
            #delete dataframe with "D" status
            for k in range(_df_adj_sum.shape[0]):
                #1.entire junction before or within first exon- delete
                if (_df_adj_sum['_status'][k] == "D"):
                    _df_adj_sum=_df_adj_sum.drop(_df_adj_sum.index[k])
            _df_adj_sum = _df_adj_sum.reset_index().drop(columns=["index"])

        # check adj matrix
        if ((self.ck.upper() == "Y") | (self.show.upper() == "Y")):
            print('Original edge table:')
            display(_df_adj_orig)
            print('Modified edge table:')
            display(_df_adj_sum)
            print("Number of Node:", n_node)
            print("Number of Edge:", _df_adj_sum.shape[0])
        if (self.ck.upper() == "Y"):
            #the size of the modified adj table show >= exon number 
            if (_df_adj_sum.shape[0] <n_node-1): 
                print("=============EDGE NUMBER < NODE NUMBER-1===================")
            #check if the modified table has and D row didn't delete
            #check if modified table has node is with x.5
            for m in range(_df_adj_sum.shape[0]):
                if (_df_adj_sum["start"][m] % 1 == 0.5) | (_df_adj_sum["end"][m] % 1 == 0.5) | (_df_adj_sum["_status"][m] == "D"):
                    print("=============SPECIAL CASE NEED MODIFICATION===================")
                    print("GENEID ====, ", self.id, "Start=====, ", _df_adj_sum["start"][m])
                #check junction span over more than 1 exon:
                if _df_adj_sum["end"][m] - _df_adj_sum["start"][m] > 1:
                    print("=============JUNCTION SPAN===================")
                    print("GENEID ====, ", self.id, "Start=====, ", _df_adj_sum["start"][m])


        #initialize an adjacent matrix 
        am = np.zeros(shape=(n_node,n_node))

        #update adjacency using edge dataframe weight
        for i in range(0,len(am)):
            for j in range(0,len(am)):
                for k in range(_df_adj_sum.shape[0]):
                    if (_df_adj_sum['start'][k] == i) & (_df_adj_sum['end'][k] == j):
                        am[i][j] = _df_adj_sum["weight"][k]
                        break
        
        self.adj_mat = am.flatten()
                    
        if ((self.ck.upper() == "Y") | (self.show.upper() == "Y")):
            print("ADJ MATRIX:", self.adj_mat)
    
    """
    Show ADJACENCY MATRIX
    """
    #initialize direct graph
    def adj_show(self):
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
                print("Sample = ",self.srr, ", Gene id = ",self.id, "is running.")
                self.exon_read()
                self.jun_read()
                self.count_fea()
                #combine all feature matrix of genes per cell 
                temp_indx = self.gtf[(self.gtf['Geneid'] == self.id)].index.to_list()
                for ind in range(0,len(temp_indx)):
                    vec_f[temp_indx[ind]] = self.feat_mat[ind]
                
                self.count_adj()
                #combine all adjacency matrix of genes per cell 
                start_indx = int(self.adj_ind[(self.adj_ind["geneid"] == self.id)]['ind_st'].values[0])
                for i in range(0,len(self.adj_mat)):
                    vec_a[i+start_indx] = self.adj_mat[i]

            #save the final data to csv file
            np.savetxt(os.path.join(self.main_path, "06_graph_mtx", self.srr+ "_fea.csv"), vec_f, fmt='%10.4f', delimiter = ',')
            np.savetxt(os.path.join(self.main_path, "06_graph_mtx", self.srr+ "_adj.csv"), vec_a, fmt='%10.4f', delimiter = ',')
            print("========================Sample",self.srr,"is Done========================")
            return(vec_f,vec_a)
        else:
            #if only running one gene per sample, no need to produce tensor output, matrix is enough for debug
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

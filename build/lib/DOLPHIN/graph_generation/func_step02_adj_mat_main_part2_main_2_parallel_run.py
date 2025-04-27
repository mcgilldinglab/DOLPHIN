from multiprocessing import Process

import pandas as pd
import numpy as np
import anndata
from matplotlib import pyplot as plt
from collections import Counter
from itertools import combinations
from itertools import tee
import pickle
#doing convolution as pre-process method
from scipy.linalg import fractional_matrix_power
import os

start = 0
incr = 32

out_name = "LUAD"
out_directory = os.path.join("..", "data")
metadata = "../final_cancer_cell_id_07_0.csv"

# os.mkdir(os.path.join(out_directory,"adj_comp_matrix"))

print("Start loading")
# adata_adj_orig = anndata.read_h5ad("../data/Adjacency_fslaSimu.h5ad")
adata_adj_orig=anndata.read_h5ad(os.path.join(out_directory, "Adjacency_"+out_name+".h5ad"))
adata_adj_up = anndata.AnnData.copy(adata_adj_orig)
adata_adj_up.X = adata_adj_orig.X.todense()

# adata_fea_orig = anndata.read_h5ad("../data/Feature_fslaSimu.h5ad")
adata_fea_orig = anndata.read_h5ad(os.path.join(out_directory, "Feature_"+out_name+".h5ad"))
adata_fea_orig.X = adata_fea_orig.X.todense()

print("loading is done")
#get gene list which has non-zero values in any of samples
df_fea = adata_fea_orig[:,:].to_df().T
df_fea['gene_id'] = adata_fea_orig[0,:].var["gene_id"]
df_fea_sum = df_fea.groupby("gene_id").sum().reset_index()
df_fea_sum = df_fea_sum.set_index("gene_id")
nonzero_gene_list = df_fea_sum[~(df_fea_sum == 0).all(axis=1)].index.tolist()
allzero_gene_list = df_fea_sum[(df_fea_sum == 0).all(axis=1)].index.tolist()

#get exon number per gene
exon_count = pd.DataFrame(adata_fea_orig.var).reset_index().groupby("gene_id").last().reset_index()
exon_count["Exon_Number"] = exon_count["index"].apply(lambda x: int(x.split('-')[-1]))
exon_gene_list = exon_count[exon_count["Exon_Number"] > 1]["gene_id"].tolist()
one_exon_gene_list = exon_count[exon_count["Exon_Number"] == 1]["gene_id"].tolist()

#only keep gene who has more than one exon and has value in any of samples
gene_list = list(set(nonzero_gene_list) & set(exon_gene_list))

pd_gt = pd.read_csv(metadata, sep='\t')
sample_list = list(pd_gt["CB"])

print("start process")
"""
Load all the necessary files
"""

def func1():
    for j in range(start,start+incr*1):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func2():
    for j in range(start+incr*1,start+incr*2):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func3():
    for j in range(start+incr*2,start+incr*3):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func4():
    for j in range(start+incr*3,start+incr*4):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func5():
    for j in range(start+incr*4,start+incr*5):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func6():
    for j in range(start+incr*5,start+incr*6):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func7():
    for j in range(start+incr*6,start+incr*7):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func8():
    for j in range(start+incr*7,start+incr*8):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func9():
    for j in range(start+incr*8,start+incr*9):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func10():
    for j in range(start+incr*9,start+incr*10):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func11():
    for j in range(start+incr*10,start+incr*11):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func12():
    for j in range(start+incr*11,start+incr*12):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func13():
    for j in range(start+incr*12,start+incr*13):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func14():
    for j in range(start+incr*13,start+incr*14):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func15():
    for j in range(start+incr*14,start+incr*15):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func16():
    for j in range(start+incr*15,start+incr*16):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func17():
    for j in range(start+incr*16,start+incr*17):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func18():
    for j in range(start+incr*17,start+incr*18):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func19():
    for j in range(start+incr*18,start+incr*19):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func20():
    for j in range(start+incr*19,start+incr*20):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func21():
    for j in range(start+incr*20,start+incr*21):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func22():
    for j in range(start+incr*21,start+incr*22):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func23():
    for j in range(start+incr*22,start+incr*23):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func24():
    for j in range(start+incr*23,start+incr*24):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

def func25():
    for j in range(start+incr*24,start+incr*25):
        for i in range(0, len(gene_list)):
            print(sample_list[j], i)
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            temp_fea = np.array(adata_fea_orig[sample_list[j], adata_fea_orig.var["gene_id"] == gene_list[i]].X[0])
            temp_adj = np.array(adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X[0])
            non_zero_index = np.where(temp_fea > 0)[0]
            _size = int(np.sqrt(len(temp_adj)))
            _old_adj = temp_adj.reshape(_size,_size)
            _new_adj = np.zeros((_size,_size))

            #get adjacency list
            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return zip(a, b)

            _adj_pair = list(map(lambda x: (x[0], x[1]), pairwise(non_zero_index)))

            for v in combinations(non_zero_index,2):
                if (_old_adj[v] == 0) & (v in _adj_pair):
                    _new_adj[v] = 1
                else:
                    _new_adj[v] = _old_adj[v]
                    
            #update adjacency data
            adata_adj_up[sample_list[j], adata_adj_up.var["gene_id"] == gene_list[i]].X = _new_adj.flatten()
        #save each individual sample h5ad file, it was runned on compute canda, saved under "/home/kailu/projects/def-junding/kailu/flash_seq_gvae/data/adj_comp_matrix/"
        adata_adj_up[sample_list[j],:].write_h5ad(os.path.join(out_directory,"adj_comp_matrix",sample_list[j]+".h5ad"))

# #run all the below functions at the same time
if __name__=='__main__':
    p1 = Process(target=func1)
    p1.start()
    p2 = Process(target=func2)
    p2.start()
    p3 = Process(target=func3)
    p3.start()
    p4 = Process(target=func4)
    p4.start()
    p5 = Process(target=func5)
    p5.start()
    p6 = Process(target=func6)
    p6.start()
    p7 = Process(target=func7)
    p7.start()
    p8 = Process(target=func8)
    p8.start()
    p9 = Process(target=func9)
    p9.start()
    p10 = Process(target=func10)
    p10.start()
    p11 = Process(target=func11)
    p11.start()
    p12 = Process(target=func12)
    p12.start()
    p13 = Process(target=func13)
    p13.start()
    p14 = Process(target=func14)
    p14.start()
    p15 = Process(target=func15)
    p15.start()
    p16 = Process(target=func16)
    p16.start()
    p17 = Process(target=func17)
    p17.start()
    p18 = Process(target=func18)
    p18.start()
    p19 = Process(target=func19)
    p19.start()
    p20 = Process(target=func20)
    p20.start()
    p21 = Process(target=func21)
    p21.start()
    p22 = Process(target=func22)
    p22.start()
    p23 = Process(target=func23)
    p23.start()
    p24 = Process(target=func24)
    p24.start()
    p25 = Process(target=func25)
    p25.start()
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

print('here1')
'''
step1:
load gene annoataion file, gene id and gene name (one to one)
'''
df_an = pd.read_csv('/mnt/data/kailu/10x/GSE125970/02_h5ad/ensembl_mod_annotation.csv', skiprows = 0, sep = ',')
print('here2')

''' 
step2: annotate gtf to get gene names
gtf file order is the same as feature/adj orders
'''
gtf = pd.read_pickle("/mnt/data/kailu/STAR_example/model/data/gene_annotation/gtf.pkl")
gtf = gtf[["seqname","gene_id","start","end"]]
gtf = gtf.rename(columns={"seqname":"Chr","gene_id":"Geneid","start":"Start","end":"End"})
gtf["Start"] = gtf["Start"].astype(int)
gtf["End"] = gtf["End"].astype(int)
print('here3')

''' 
step3: combine gtf and annotation files to get gene-exon name
'''
df_gtf_an = pd.merge(gtf, df_an, how='left', left_on=["Geneid"], right_on=["gene_id"])
df_gtf_an["Gene_Exon_name"] = df_gtf_an.groupby(['gene_name']).cumcount()+1
df_gtf_an["Gene_Exon_name"] = df_gtf_an['gene_name'].astype(str) +"-"+ df_gtf_an["Gene_Exon_name"].astype(str)
print('here4')

'''
Get gtf adjacency information
'''
adj_index = pd.read_csv('/mnt/data/kailu/10x/GSE125970/data/adj_index.csv', skiprows = 0, sep = ',',low_memory=False,index_col=False)

df_jun_gtf = pd.DataFrame()
exon_gene_list = list(dict.fromkeys(adj_index["geneid"].to_list())) #61860 genes
for i in range(0,len(exon_gene_list)):
# for i in range(0,2):
    print(i)
    _rep_num = adj_index.loc[i]['ind']
    _temp_genename = df_gtf_an[df_gtf_an["Geneid"] == exon_gene_list[i]]["gene_name"].iloc[0]
    _temp_exon = pd.DataFrame([{'Geneid':exon_gene_list[i], 'GeneName':_temp_genename}]) 
    _temp_exon = _temp_exon.loc[_temp_exon.index.repeat(_rep_num)]
    df_jun_gtf = pd.concat([df_jun_gtf, _temp_exon])
df_jun_gtf["Gene_Junc_name"] = df_jun_gtf.groupby(['GeneName']).cumcount()+1
df_jun_gtf["Gene_Junc_name"] = df_jun_gtf['GeneName'].astype(str) +"-"+ df_jun_gtf["Gene_Junc_name"].astype(str)
df_jun_gtf = df_jun_gtf.reset_index(drop=True)

df_jun_gtf.to_pickle("df_jun_gtf.pkl")
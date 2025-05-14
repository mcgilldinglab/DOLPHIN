from tqdm import tqdm
import math
import anndata
import os
import pandas as pd
import torch
from func_step01_fea_mat_main_part1 import combine_fea
from func_step01_fea_mat_main_part2 import fea_comp
from func_step02_adj_mat_main_part1_main_1 import combine_adj
# #miss one "step02_adj_mat_main_part2_main_2_parallel_run_X"
from func_step02_adj_mat_main_part2_main_3_combine import combine_adj_comp
from func_step02_adj_mat_main_part3_main_1 import adj_comp_re_part1
from func_step02_adj_mat_main_part3_main_2 import adj_comp_re_part2
from func_step02_adj_mat_main_part3_main_3 import adj_comp_re_part3
from func_exon_gene import get_exon_gene
from func_vae_fea import get_fea_adj
from func_vae_adj import get_vae_adj
from func_step03_GNN_main import get_graph_input

batch_name='batch_50k_s0'
metadata = "/mnt/data/kailu/00_scExon/07_simulation_3prime/00_simulation/samples.txt"
df_label = pd.read_csv(metadata, sep='\t', names=["CB"])
total_sample_size = len(df_label)
out_name = "Simu_"+batch_name
fea_run_num = 100
adj_run_num = 25
gnn_run_num = 100
out_directory = os.path.join("..", "data", batch_name)
os.makedirs(out_directory, exist_ok=True) # create an output folder is not exist
graph_directory = os.path.join("..", "06_graph_mtx", batch_name)
exon_gene_featurecount_directory=os.path.join("..", "04_exon_gene_cnt", batch_name)
gene_annotation = "/mnt/data/kailu/00_scExon/00_pre_process/gene_id_annotation.csv"
gtf_pkl_path = "/mnt/data/kailu/STAR_example/ensembl_mod/gtf.pkl"
gtf_jun_pkl_path = "/mnt/data/kailu/STAR_example/ensembl/df_jun_gtf.pkl"
gtf_path= "/mnt/data/kailu/STAR_example/ensembl/Homo_sapiens.GRCh38.107.gtf"
# celltype_mapper = {'B': 0, 'CD4 T': 1, 'CD8 T': 2, 'DC': 3, 'Mono': 4, 'NK': 5, 'other T': 6, 'other': 7}

# ==================================
# Adjacency Compression Version
# ==================================
# generate compressed version adjacency matrix
print("Start Generating Compressed version of Adjacency Matrix")
with tqdm(total=total_sample_size) as pbar:
    for i in range(0, total_sample_size):
        if i%adj_run_num==0:
            pbar=combine_adj_comp(pbar, df_label, start_idx = i, sample_num=adj_run_num, output_path=out_directory, output_name= out_name)

##combine the compressed adjacency matrix
_anndata_ca_dict = {}

total_number_anndata = math.ceil(total_sample_size/adj_run_num)
_ca_all_names_lst = []
print("Combining ", total_number_anndata, " compressed adjacency anndata")
for j in range(0, total_number_anndata):
    print(j)
    _temp_ad = anndata.read_h5ad(os.path.join(out_directory, "AdjacencyComp_"+out_name+"_"+str(j)+".h5ad"))
    _anndata_ca_dict["CA"+str(j%5)] = _temp_ad

    if j < (len(range(0, total_number_anndata)) // 5) * 5:
        if (j+1)%5 == 0:
            for _i, (_, _ad) in enumerate(_anndata_ca_dict.items()):
                if _i == 0:
                    combine_anndata = _ad
                else:
                    combine_anndata = combine_anndata.concatenate(_ad, index_unique = None, batch_key = None)
            _out_ca_name = "AdjacencyComp_"+out_name+"_"+str(j-4)+"_"+str(j)+".h5ad"
            combine_anndata.write(os.path.join(out_directory,_out_ca_name))
            _anndata_ca_dict = {}
            _ca_all_names_lst.append(_out_ca_name)
    elif j>=(len(range(0, total_number_anndata)) // 5) * 5 and j==max(range(0, total_number_anndata)):
        for _i, (_, _ad) in enumerate(_anndata_ca_dict.items()):
            if _i == 0:
                combine_anndata = _ad
            else:
                combine_anndata = combine_anndata.concatenate(_ad, index_unique = None, batch_key = None)
        _out_ca_name = "AdjacencyComp_"+out_name+"_"+str(total_number_anndata%5+1)+"_"+str(j)+".h5ad" 
        combine_anndata.write(os.path.join(out_directory, _out_ca_name))
        _anndata_ca_dict = {}
        _ca_all_names_lst.append(_out_ca_name)

_anndata_ca_all_dict = {}
for _i, _ca_name in enumerate(_ca_all_names_lst):
    _temp_ad = anndata.read_h5ad(os.path.join(out_directory, _ca_name))
    _anndata_ca_all_dict["CA_all"+str(_i)] = _temp_ad
    for _i, (_, _ad) in enumerate(_anndata_ca_all_dict.items()):
        if _i == 0:
            combine_anndata = _ad
        else:
            combine_anndata = combine_anndata.concatenate(_ad, index_unique = None, batch_key = None)
    combine_anndata.write(os.path.join(out_directory,"AdjacencyComp_"+out_name+".h5ad"))

# ==================================
# Adjacency Compression Remove Exon Version
# ==================================
adj_comp_re_part1(out_directory, out_name)
adj_comp_re_part2(out_directory, out_name)
adj_comp_re_part3(out_directory, out_name)

# ====================================
# Find HVG using ExonGene count table, modify the code below to make three more functions
# ====================================
get_exon_gene(df_label, exon_gene_featurecount_directory, gtf_path, out_directory, out_name)
get_fea_adj(out_directory, out_name)
get_vae_adj(out_directory, out_name)

# ==================================
# Construct data input for Geometric
# ==================================
print("Start Construct Data Input for Pytorch Geometric")
with tqdm(total=total_sample_size) as pbar_gnn:
    for i in range(0, total_sample_size):
        if i%gnn_run_num==0:
            pbar_gnn= get_graph_input(pbar_gnn, i, gnn_run_num, out_directory, out_name, celltypename=None, mapper=None)

##### combine all  geometric .pt files
total_number_gnn_anndata = math.ceil(total_sample_size/gnn_run_num)
for _idx, _gnn_idx in enumerate(range(0, total_number_gnn_anndata)):
    _temp_gnn = torch.load(os.path.join(out_directory, "geometric_"+out_name+"_"+str(_gnn_idx)+".pt"))
    if _idx ==0:
        combine_gnn = _temp_gnn
    else:
        combine_gnn += _temp_gnn
torch.save(combine_gnn, os.path.join(out_directory, "geometric_"+out_name+".pt"))

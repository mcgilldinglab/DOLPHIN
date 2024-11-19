#hyperparameter search
import yaml
import os
from train import run_train
import re
import itertools
import random
import numpy as np
import torch

"""
The main function, hyperparameter search
"""

import scanpy as sc

#load data
in_path_gra = "../DOLPHIN/test_datasets/geometric_smart_seq_example.pt"
in_path_fea = "../DOLPHIN/test_datasets/smart_seq_example_FeatureCompHvg5k.h5ad"

output_path = "./output"

#check number of run_x folder under the output_path, so rerun will not start from 0.
# n_folders = 0
# if len(os.listdir(output_path)) > 0:
#     for i in os.listdir(output_path):
#         if i.find("run_") > -1:
#             temp_folder = int(re.sub('run_','', i))
#             if temp_folder > n_folders:
#                 n_folders = temp_folder
#     n_folders = n_folders + 1

seed_number = 0

def seed_everything(seed: int):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = True
    torch.backends.cudnn.enabled = True
    
seed_everything(seed_number)

import argparse
from argparse import Namespace
import pickle

# https://stackoverflow.com/questions/50869939/how-to-load-argparse-from-namespace-string
parser = argparse.ArgumentParser()
parser.add_argument('--p')
args = parser.parse_args()

with open(args.p, 'rb') as fp:
    current_hyper = pickle.load(fp)
print(current_hyper)

# 1. Define an objective function to be maximized.
# def objective(trial):
device = 'cuda:0'
# define paths
# run_name = args.p.split("/")[-1].split(".pkl")[0]

current_out_path = output_path+"/run_seednumber_" + str(seed_number)

# 2. Suggest values of the hyperparameters using a trial object.
params = {"gat_channel":current_hyper["gat_channel"],
        "nhead": current_hyper["nhead"],
        "gat_dropout": current_hyper["gat_dropout"],
        "concat": False,
        "list_gra_enc_hid": current_hyper["list_gra_enc_hid"],
        "gra_p_dropout":current_hyper["gra_p_dropout"],
        "z_dim": current_hyper["z_dim"],
        "list_fea_dec_hid":current_hyper["list_fea_dec_hid"],
        "list_adj_dec_hid":current_hyper["list_adj_dec_hid"],
        "lr": 1.0e-3,
        "batch": 20,
        "epochs": 2000,
        "pre_train_load": False,
        "kl_beta": current_hyper["kl_beta"],
        "fea_lambda": current_hyper["fea_lambda"], 
        "adj_lambda": 1 - current_hyper["fea_lambda"],
        "early_stopping": False,
        "progress_checking":False,
        "draw_freq":50
}

## 3. create a folder to store all information per trial
os.makedirs(current_out_path)

#save hyperparameters
with open(current_out_path+'/params.yml', 'w') as outfile:
    yaml.dump(params, outfile, default_flow_style=False)

## 4. train the model and get evaluation score
score = run_train(in_path_gra, in_path_fea, current_out_path, params, device)

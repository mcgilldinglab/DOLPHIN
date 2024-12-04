#hyperparameter search
# import yaml
import os
from .train import run_train
# import re
# import itertools
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

device = 'cuda:0'

current_out_path = output_path+"/run_seednumber_" + str(seed_number)

current_hyper = {'gat_channel': [2],
 'nhead': 9,
 'gat_dropout': 0.3,
 'gra_p_dropout': 0.2,
 'list_gra_enc_hid': [256],
 'list_fea_dec_hid': [256],
 'list_adj_dec_hid': [128],
 'z_dim': 50,
 'kl_beta': 0.7,
 'fea_lambda': 0.5}

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
        "batch": 1,
        "epochs": 200,
        "pre_train_load": False,
        "kl_beta": current_hyper["kl_beta"],
        "fea_lambda": current_hyper["fea_lambda"], 
        "adj_lambda": 1 - current_hyper["fea_lambda"],
        "early_stopping": False,
        "progress_checking":False,
        "draw_freq":50
}

def run_DOLPHIN(data_type, fea_in, graph_in, current_out_path, params = params, device='cuda:0'):
    if data_type == "Full-length":
        run_train(graph_in, fea_in, current_out_path, params, device)
    else:
        run_train(graph_in, fea_in, current_out_path, params, device)
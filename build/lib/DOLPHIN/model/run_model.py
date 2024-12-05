import os
from .train import run_train
import numpy as np
import torch
import pickle
import random

"""
The main function, hyperparameter search
"""

def run_DOLPHIN(data_type, graph_in, fea_in, current_out_path='./', params=None, device='cuda:0', seed_num=0):
    # print(seed_num)
    random.seed(seed_num)
    os.environ['PYTHONHASHSEED'] = str(seed_num)
    np.random.seed(seed_num)
    torch.manual_seed(seed_num)
    torch.cuda.manual_seed(seed_num)
    torch.cuda.manual_seed_all(seed_num)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = True
    torch.backends.cudnn.enabled = True
        
    default_params_full_length = {"gat_channel":[2],
            "nhead": 9,
            "gat_dropout": 0.3,
            "concat": False,
            "list_gra_enc_hid": [256],
            "gra_p_dropout":0.2,
            "z_dim": 50,
            "list_fea_dec_hid":[256],
            "list_adj_dec_hid":[128],
            "lr": 1.0e-3,
            "batch": 20,
            "epochs": 200,
            "kl_beta": 0.7,
            "fea_lambda": 0.5, 
            "adj_lambda": 0.5,
    }
    
    default_params_10x = {"gat_channel":[9],
            "nhead": 1,
            "gat_dropout": 0.1,
            "concat": False,
            "list_gra_enc_hid": [512],
            "gra_p_dropout":0.3,
            "z_dim": 35,
            "list_fea_dec_hid":[512],
            "list_adj_dec_hid":[256],
            "lr": 1.0e-3,
            "batch": 20,
            "epochs": 200,
            "kl_beta": 0.7,
            "fea_lambda": 0.5, 
            "adj_lambda": 0.5,
    }
    
    # Select the default parameters based on data_type
    if data_type == "full-length":
        default_params = default_params_full_length
    else:
        default_params = default_params_10x

    # If params is None, use default parameters
    if params is None:
        params = default_params
    else:
        # Update the default parameters with user-provided parameters
        default_params.update(params)
        params = default_params

    # print(params)
    
    run_train(graph_in, fea_in, current_out_path, params, device='cuda')


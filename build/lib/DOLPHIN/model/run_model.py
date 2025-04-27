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
    
    """
    Run the DOLPHIN model on single-cell RNA-seq data to obtain latent cell embeddings.

    Parameters
    ----------
    data_type : str
        Specifies the type of input single-cell RNA-seq data.
        - "full-length": For full-length RNA-seq data.
        - "10x": For 10x Genomics RNA-seq data.

    graph_in : object
        The input graph structure (precomputed from exon-level data).

    fea_in : anndata.AnnData
        The input feature matrix, provided as an AnnData object.

    current_out_path : str, optional
        Output directory where the resulting cell embeddings will be saved.
        The embeddings will be written to `DOLPHIN_Z.h5ad`. Default is `'./'`.

    params : dict, optional
        A dictionary of model hyperparameters. If not provided, default parameters will be used
        depending on `data_type`. Customizable parameters include:

        - "gat_channel"       : Number of GAT output channels per head.
        - "nhead"             : Number of GAT attention heads.
        - "gat_dropout"       : Dropout rate in the GAT layer.
        - "list_gra_enc_hid"  : Encoder MLP hidden layer sizes.
        - "gra_p_dropout"     : Dropout rate in the encoder.
        - "z_dim"             : Dimensionality of the latent space.
        - "list_fea_dec_hid"  : Feature decoder MLP hidden layer sizes.
        - "list_adj_dec_hid"  : Adjacency decoder MLP hidden layer sizes.
        - "lr"                : Learning rate.
        - "batch"             : Mini-batch size.
        - "epochs"            : Number of training epochs.
        - "kl_beta"           : KL divergence loss weight.
        - "fea_lambda"        : Feature reconstruction loss weight.
        - "adj_lambda"        : Adjacency reconstruction loss weight.

    device : str, optional
        Device to run the model on. Default is `'cuda:0'` (recommended for GPU acceleration).

    seed_num : int, optional
        Random seed for reproducibility. Default is `0`.

    Returns
    -------
    None
        Saves the latent cell embedding matrix to `DOLPHIN_Z.h5ad` under `current_out_path`.
    """
    
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


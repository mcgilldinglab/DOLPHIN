import pyro
import torch
import matplotlib.pyplot as plt
import gc
import anndata
from torch_geometric.loader import DataLoader
from .model import define_svi
from .utils import z_leiden_cluster

"""
Define training process
"""

# train function
def train_step(svi, train_loader, device):
    epoch_loss = 0.
    for x_gra in train_loader:
        # if "cuda" in device:
        x_gra = x_gra.to(device)
        epoch_loss += svi.step(x_gra)
    normalizer_train = len(train_loader.dataset)
    total_epoch_loss_train = epoch_loss / normalizer_train
    return total_epoch_loss_train

def run_train(in_path_gp, in_path_fea, out_path, params, device, pretrain_fea=None, pretrain_adj=None):
    pyro.clear_param_store()
    gc.collect()
    torch.cuda.empty_cache()

    #### step 1: Input anndata dataset
    adata_fea = anndata.read_h5ad(in_path_fea)

    #load the graph data
    pg_celldata = torch.load(in_path_gp)

    #dataloader
    all_cell_loader = DataLoader(pg_celldata, batch_size=params["batch"], shuffle=False)
    #### step 2: Define model
    vae, svi = define_svi(1, pg_celldata[0].x_fea.shape[1], pg_celldata[0].x_adj.shape[1], params, device)

    #### step 3: Train model
    train_elbo = []
    val_elbo = []
    # save runned epoch value
    all_epoch = []
    
    for epoch in range(params["epochs"]):
        all_epoch.append(epoch)
        
        ##################
        ###Train
        ##################
        total_epoch_loss_train = train_step(svi, all_cell_loader, device)
        train_elbo.append(total_epoch_loss_train)
        
        print("[epoch %03d] training loss: %.4f" % (epoch, total_epoch_loss_train))

        gc.collect()
        torch.cuda.empty_cache()

    z_leiden_cluster(out_path, vae, adata_fea, all_cell_loader, epoch, device)

import pyro
import torch
import matplotlib.pyplot as plt
import gc
import anndata
from torch_geometric.loader import DataLoader
from .model import define_svi
from .utils import z_leiden_cluster
from .pytorchtools import EarlyStopping

"""
Define training process
"""

# train function
def train_step(svi, train_loader, device):
    epoch_loss = 0.
    for x_gra in train_loader:
        if "cuda" in device:
            x_gra = x_gra.cuda(device)
        epoch_loss += svi.step(x_gra)

    # return epoch loss by average loss per epoch
    normalizer_train = len(train_loader.dataset)
    total_epoch_loss_train = epoch_loss / normalizer_train
    return total_epoch_loss_train

# def validate_step(svi, val_loader, device):
#     val_epoch_loss = 0.

#     with torch.no_grad():
#         for x_gra_val in val_loader:
#             if "cuda" in device:
#                 x_gra_val = x_gra_val.cuda(device)
#             val_epoch_loss += svi.evaluate_loss(x_gra_val)

#     # return epoch loss by average loss per epoch
#     val_normalizer_train = len(val_loader.dataset)
#     total_epoch_loss_val = val_epoch_loss / val_normalizer_train
#     return total_epoch_loss_val

def run_train(in_path_gp, in_path_fea, out_path, params, device, pretrain_fea=None, pretrain_adj=None):
    pyro.clear_param_store()
    gc.collect()
    torch.cuda.empty_cache()

    #### step 1: Input anndata dataset
    adata_fea = anndata.read_h5ad(in_path_fea)

    pg_celldata = torch.load(in_path_gp)

    #split dataset for train and validate dataset
    train_idx, val_idx = torch.utils.data.random_split(pg_celldata, [len(pg_celldata) - int(len(pg_celldata)*0.2), int(len(pg_celldata)*0.2)])
    train_cell = [pg_celldata[i] for i in train_idx.indices]
    val_cell = [pg_celldata[i] for i in val_idx.indices]

    #dataloader
    # train_cell_loader = DataLoader(train_cell, batch_size=params["batch"], shuffle=False, drop_last=True)
    # val_cell_loader = DataLoader(val_cell, batch_size=params["batch"], shuffle=False, drop_last=True) 
    #for checking all ari score purpose
    all_cell_loader = DataLoader(pg_celldata, batch_size=params["batch"], shuffle=False)
    #### step 2: Define model
    vae, svi = define_svi(1, pg_celldata[0].x_fea.shape[1], pg_celldata[0].x_adj.shape[1], params, device)

    #### step 3: Train model
    train_elbo = []
    # val_elbo = []
    # save runned epoch value
    all_epoch = []
    # if params["early_stopping"]:
    #     final_epoch = []
        
     # initialize the early_stopping object
    if params["early_stopping"]:
        early_stopping = EarlyStopping(patience=45, verbose=True, path=out_path+"/checkpoint.pt")
    for epoch in range(params["epochs"]):
        all_epoch.append(epoch)
        
        ##################
        ###Train
        ##################
        total_epoch_loss_train = train_step(svi, all_cell_loader, device)
        train_elbo.append(total_epoch_loss_train)
        
        ##################
        ###Validate
        ##################
        # total_epoch_loss_val = validate_step(svi, val_cell_loader, device)
        # val_elbo.append(total_epoch_loss_val)
        # final_train_loss.append(total_epoch_loss_train)
        # print("[epoch %03d] training loss: %.4f; validation loss: %.4f; " % (epoch, total_epoch_loss_train, total_epoch_loss_val))
        # print("[epoch %03d] training loss: %.4f" % (epoch, total_epoch_loss_train))
        ##check latent space clustering as criteria for hyperparameter search, run for entire dataset
        ## create a folder to store each trail and output
            
        ##################
        ###Early Stopping
        ##################
        # early_stopping needs the validation loss to check if it has decresed, 
        # # and if it has, it will make a checkpoint of the current model
        # if params["early_stopping"]:
        #     early_stopping(total_epoch_loss_val, vae)

        #     if early_stopping.early_stop:
        #         print("!!!!!!!Early stopping")
        #         print("!!!!!!Early stopping at epoch: %03d" % epoch)
        #         break
        
        ##################
        ###ARI Checking
        ##################
        if params["progress_checking"]:
            if ((epoch+1) % params["draw_freq"] == 0) and epoch > 0:
                final_epoch.append(epoch)
                z_leiden_cluster(out_path, vae, adata_fea, all_cell_loader, epoch, device)
                # final_type1_ARI.append(type1_ARI)
                # final_type2_ARI.append(type2_ARI)

        gc.collect()
        # allow reusable memory to be free
        torch.cuda.empty_cache()

    z_leiden_cluster(out_path, vae, adata_fea, all_cell_loader, epoch, device)

    return 0

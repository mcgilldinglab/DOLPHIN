import pyro
import torch
import matplotlib.pyplot as plt
import gc
import anndata
from torch_geometric.loader import DataLoader
from model import define_svi
from utils import z_leiden_cluster
from pytorchtools import EarlyStopping

"""
Define training process
"""

# train function
def train_step(svi, train_loader, device):
    # initialize loss accumulator
    epoch_loss = 0.
    # do a training epoch over each mini-batch x returned
    # by the data loader
    for x_gra in train_loader:
        # if on GPU put mini-batch into CUDA memory
        if "cuda" in device:
            x_gra = x_gra.cuda(device)
        # do ELBO gradient and accumulate loss
        epoch_loss += svi.step(x_gra)

    # return epoch loss by average loss per epoch
    normalizer_train = len(train_loader.dataset)
    total_epoch_loss_train = epoch_loss / normalizer_train
    return total_epoch_loss_train

def validate_step(svi, val_loader, device):
    # initialize loss accumulator
    val_epoch_loss = 0.
    # do a training epoch over each mini-batch x returned
    # by the data loader
    with torch.no_grad():
        for x_gra_val in val_loader:
            # if on GPU put mini-batch into CUDA memory
            if "cuda" in device:
                x_gra_val = x_gra_val.cuda(device)
            # do ELBO gradient and accumulate loss
            val_epoch_loss += svi.evaluate_loss(x_gra_val)

    # return epoch loss by average loss per epoch
    val_normalizer_train = len(val_loader.dataset)
    total_epoch_loss_val = val_epoch_loss / val_normalizer_train
    return total_epoch_loss_val

def run_train(in_path_gp, in_path_fea, out_path, params, device, pretrain_fea=None, pretrain_adj=None):
    ## step 0: memory free and clean parameters
    #clear param store
    pyro.clear_param_store()
    #garbage collection: clean up objects, it's automatic but does not run when device is out of memory
    gc.collect()
    # allow reusable memory to be freed
    torch.cuda.empty_cache()

    #### step 1: Input anndata dataset
    adata_fea = anndata.read_h5ad(in_path_fea)

    #load the graph data
    pg_celldata = torch.load(in_path_gp)

    #split dataset for train and validate dataset
    train_idx, val_idx = torch.utils.data.random_split(pg_celldata, [len(pg_celldata) - int(len(pg_celldata)*0.2), int(len(pg_celldata)*0.2)])
    train_cell = [pg_celldata[i] for i in train_idx.indices]
    val_cell = [pg_celldata[i] for i in val_idx.indices]

    #dataloader
    train_cell_loader = DataLoader(train_cell, batch_size=params["batch"], shuffle=False, drop_last=True)
    val_cell_loader = DataLoader(val_cell, batch_size=params["batch"], shuffle=False, drop_last=True) #drop the last incomplete batch
    #for checking all ari score purpose
    all_cell_loader = DataLoader(pg_celldata, batch_size=params["batch"], shuffle=False)
    #### step 2: Define model
    vae, svi = define_svi(1, pg_celldata[0].x_fea.shape[1], pg_celldata[0].x_adj.shape[1], params, device)

    #### step 3: Train model
    train_elbo = []
    val_elbo = []
    # save runned epoch value
    all_epoch = []
    if params["early_stopping"]:
        # save epoch list
        final_epoch = []
        # save training loss
        # final_train_loss = []
        # save ARI score
        final_type_ARI = []
        
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
        print("[epoch %03d] training loss: %.4f" % (epoch, total_epoch_loss_train))
        ##check latent space clustering as criteria for hyperparameter search, run for entire dataset
        ## create a folder to store each trail and output
            
        ##################
        ###Early Stopping
        ##################
        # early_stopping needs the validation loss to check if it has decresed, 
        # and if it has, it will make a checkpoint of the current model
        if params["early_stopping"]:
            early_stopping(total_epoch_loss_val, vae)

            if early_stopping.early_stop:
                print("!!!!!!!Early stopping")
                print("!!!!!!Early stopping at epoch: %03d" % epoch)
                break
        
        ##################
        ###ARI Checking
        ##################
        if params["progress_checking"]:
            if ((epoch+1) % params["draw_freq"] == 0) and epoch > 0:
                final_epoch.append(epoch)
                type1_ARI, type2_ARI = z_leiden_cluster(out_path, vae, adata_fea, all_cell_loader, epoch, device)
                final_type1_ARI.append(type1_ARI)
                final_type2_ARI.append(type2_ARI)
        else:
            if (epoch == 1999) | (epoch == 999):
                type1_ARI, type2_ARI = z_leiden_cluster(out_path, vae, adata_fea, all_cell_loader, epoch, device)
        
        gc.collect()
        # allow reusable memory to be free
        torch.cuda.empty_cache()

    #return the highest ARI score per each trial
    criteria = type2_ARI
    
    #save final scores
    if params["early_stopping"]:
        f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (15,5))
        ax1.plot(all_epoch, train_elbo, label='Training Loss') 
        ax1.plot(all_epoch, val_elbo, label='Validation Loss') 
        ax1.legend()
        ax1.set_title("Loss")
        
        ax2.plot(final_epoch, final_type1_ARI)
        ax2.set_title("Cell tyep1 ARI score")
        
        ax3.plot(final_epoch, final_type2_ARI)
        ax3.set_title("Cell tyep2 ARI score")
        plt.savefig(out_path+"/train_all_loss.png", bbox_inches="tight")
    else:
        plt.plot(all_epoch, train_elbo, label='Training Loss') 
        # plt.plot(all_epoch, val_elbo, label='Validation Loss') 
        plt.legend()
        plt.title("Loss")
        
        plt.savefig(out_path+"/train_all_loss.png", bbox_inches="tight")

    return criteria
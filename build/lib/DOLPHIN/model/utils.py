import torch
import numpy as np
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import scanpy as sc
import matplotlib.pyplot as plt
import os

"""
Define dataset and latent space plot
"""

####function to produce leiden clustering of latent space
#get latent space
def z_leiden_cluster(out_path, model, adata, data_loader, epoch, device):
    mu_latent_all=[]
    #switch to evalution step
    model.eval()
    with torch.no_grad():
    #step1:get all cells latent representation
        for batch_gra in data_loader:
            local_gra=batch_gra.cuda(device)
            z_mu, z_mu_add_var = model.getZ(local_gra)
            np_z_mu = z_mu.cpu().detach().numpy()
            mu_latent_all.append(np_z_mu)
        mu_latent_all = np.concatenate(mu_latent_all)
    
    sample_type = "mu"
        
    ### step2: create annadata for latent space
    adata.obsm["X_z"] = mu_latent_all
    #save the latent representation results
    adata.write(os.path.join(out_path,"DOLPHIN_Z.h5ad"))
    
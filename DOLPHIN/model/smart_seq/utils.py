import torch
import numpy as np
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import scanpy as sc
import matplotlib.pyplot as plt

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
    #step1:get all cells latent representation (795 * latent_dim)
    #no need to shuffle here, since it's not training, it's the same order as celldata
        for batch_gra in data_loader:
            local_gra=batch_gra.cuda(device)
            #get the latent space
            z_mu, z_mu_add_var = model.getZ(local_gra)
            #remove gradient and convert to numpy
            np_z_mu = z_mu.cpu().detach().numpy()
            #combine all cell latent space together
            mu_latent_all.append(np_z_mu)
        mu_latent_all = np.concatenate(mu_latent_all)
    
    sample_type = "mu"
        
    ### step2: create annadata for latent space
    adata.obsm["X_z"] = mu_latent_all
    #save the latent representation results
    if (epoch == 1999) | (epoch ==999):
        adata.write(out_path+"/Epoch_"+str(epoch)+"_exon_results.h5ad")
    
    # ### step3: find leiden clustering of the latent space
    sc.pp.neighbors(adata, use_rep="X_z")
    sc.tl.umap(adata)
    
    ##########
    #celltype
    ##########
    celltype1_resolution_val_list = []
    celltype1_ari_val_list = []
    # ### step4: calculateARI score
    #Get best ARI by changing resolution
    for res in np.arange(0.0,2,0.1):
        celltype1_resolution_val_list.append(round(res,1))
        sc.tl.leiden(adata, resolution=round(res,1))
        ARI = adjusted_rand_score(adata.obs["predicted.celltype.l1"].tolist(), adata.obs["leiden"].tolist())
        #store ari for speicific sample type
        celltype1_ari_val_list.append(ARI)
        #store for all the ari
        # celltype1_all_ari_list.append(ARI)
    
    final_res_celltype1 = celltype1_resolution_val_list[np.where(np.array(celltype1_ari_val_list)==max(celltype1_ari_val_list))[0][0]]
    
    sc.tl.leiden(adata,resolution=round(final_res_celltype1,1))
    
    # scib.me.cluster_optimal_resolution(adata, cluster_key="leiden", label_key="predicted.celltype.l1") #score here is corresponding to NMI
    # celltype_l1_ari = scib.me.ari(adata, cluster_key="leiden", label_key="predicted.celltype.l1")
    with plt.rc_context():
        sc.pl.umap(adata, color=['leiden', "predicted.celltype.l1", "predicted.celltype.l2"], wspace=0.5, title='celltype1_ARI='+ str(round(max(celltype1_ari_val_list), 3)) + ", Epoch:" + str(epoch) + ", Sample_type: " + sample_type, show=False)
        plt.savefig(out_path+"/celltype_l1_" + sample_type + "_Epoch_" + str(epoch) + "_ARI_" + str(round(max(celltype1_ari_val_list), 3)) +".png", bbox_inches="tight")
    
    ##########
    #subcelltype
    ##########
    celltype2_resolution_val_list = []
    celltype2_ari_val_list = []
    # ### step4: calculateARI score
    #Get best ARI by changing resolution
    for res in np.arange(0.0,2,0.1):
        celltype2_resolution_val_list.append(round(res,1))
        sc.tl.leiden(adata, resolution=round(res,1))
        ARI = adjusted_rand_score(adata.obs["predicted.celltype.l2"].tolist(), adata.obs["leiden"].tolist())
        #store ari for speicific sample type
        celltype2_ari_val_list.append(ARI)
        #store for all the ari
        # celltype2_all_ari_list.append(ARI)
    
    final_res_celltype2 = celltype2_resolution_val_list[np.where(np.array(celltype2_ari_val_list)==max(celltype2_ari_val_list))[0][0]]
    
    sc.tl.leiden(adata,resolution=round(final_res_celltype2,1))
    
    #celltype
    # scib.me.cluster_optimal_resolution(adata, cluster_key="leiden", label_key="predicted.celltype.l2") #score here is corresponding to NMI
    # celltype_l1_ari = scib.me.ari(adata, cluster_key="leiden", label_key="predicted.celltype.l2")
    with plt.rc_context():
        sc.pl.umap(adata, color=['leiden', "predicted.celltype.l1", "predicted.celltype.l2"], wspace=0.5, title='celltype2_ARI='+ str(round(max(celltype2_ari_val_list), 3)) + ", Epoch:" + str(epoch) + ", Sample_type: " + sample_type, show=False)
        plt.savefig(out_path+"/celltype_l2_" + sample_type + "_Epoch_" + str(epoch) + "_ARI_" + str(round(max(celltype2_ari_val_list), 3)) +".png", bbox_inches="tight")
    
    #chanage back to train mode
    model.train()
    return max(celltype1_ari_val_list), max(celltype2_ari_val_list)
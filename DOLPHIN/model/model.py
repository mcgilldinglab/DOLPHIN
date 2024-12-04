import torch
import torch.nn as nn
import pyro
import pyro.distributions as dist
from pyro.infer import SVI, Trace_ELBO
from pyro.optim import Adam
from torch_geometric.nn import GATConv
import pyro.poutine as poutine

"""
Define VAE model
"""
pyro.distributions.enable_validation(False)
pyro.set_rng_seed(0)

# Gat layers +encoder
class Gra_Encoder(nn.Module):
    def __init__(self, in_node_fea, gat_channel, nhead, gat_dropout, concat, in_fea, hidden_dim, z_dim, p_dropout):
        super().__init__()

        #graph attention convolution layer
        self.hidden_channels = gat_channel
        self.convs = nn.ModuleList()
        for i, h_dim in enumerate(self.hidden_channels):
            if (i == 0) or (concat == False):
                self.convs.append(GATConv(in_node_fea, h_dim, heads = nhead, dropout = gat_dropout, concat = concat))
            else:
                self.convs.append(GATConv(in_node_fea * nhead, h_dim, heads = nhead, dropout = gat_dropout, concat = concat))
            in_node_fea = h_dim
        self.act = nn.ReLU(inplace=True)

        # Build Encoder
        # modules = []
        # if concat == True:
        #     for h_dim in hidden_dim:
        #         modules.append(
        #             nn.Sequential(
        #                 nn.Linear(in_fea * gat_channel[-1] * nhead, h_dim),
        #                 nn.BatchNorm1d(h_dim, momentum=0.01, eps=0.001), 
        #                 nn.LayerNorm(h_dim, elementwise_affine=False),
        #                 nn.ReLU(),
        #                 nn.Dropout(p=p_dropout))
        #         )
        #         in_fea = h_dim
        # else:
        #     for h_dim in hidden_dim:
        #         modules.append(
        #             nn.Sequential(
        #                 nn.Linear(in_fea * gat_channel[-1], h_dim),
        #                 nn.BatchNorm1d(h_dim, momentum=0.01, eps=0.001), 
        #                 nn.LayerNorm(h_dim, elementwise_affine=False),
        #                 nn.ReLU(),
        #                 nn.Dropout(p=p_dropout))
        #         )
        #         in_fea = h_dim
        
        # self.encoder = nn.Sequential(*modules)
        
        modules = []
        if concat == True:
            modules.append(nn.Sequential(
                        nn.Linear(in_fea * gat_channel[-1] * nhead, hidden_dim[0]),
                        nn.BatchNorm1d(hidden_dim[0], momentum=0.01, eps=0.001), 
                        nn.LayerNorm(hidden_dim[0], elementwise_affine=False),
                        nn.ReLU(),
                        nn.Dropout(p=p_dropout)))
        else:
            modules.append(nn.Sequential(
                        nn.Linear(in_fea * gat_channel[-1], hidden_dim[0]),
                        nn.BatchNorm1d(hidden_dim[0], momentum=0.01, eps=0.001), 
                        nn.LayerNorm(hidden_dim[0], elementwise_affine=False),
                        nn.ReLU(),
                        nn.Dropout(p=p_dropout)))
            
        if len(hidden_dim) > 1:
            in_val = hidden_dim[0]
            for h_dim in hidden_dim[1:]:
                modules.append(
                    nn.Sequential(
                        nn.Linear(hidden_dim[0] , h_dim),
                        nn.BatchNorm1d(h_dim, momentum=0.01, eps=0.001), 
                        nn.LayerNorm(h_dim, elementwise_affine=False),
                        nn.ReLU(),
                        nn.Dropout(p=p_dropout))
                )
                in_val = h_dim
        
        self.encoder = nn.Sequential(*modules)
        
        #last layer
        self.fc_mu = nn.Linear(hidden_dim[-1], z_dim)
        self.fc_var = nn.Linear(hidden_dim[-1], z_dim)
        
    def forward(self, data, batch):
        #graph attention layer
        x = data.x
        for i in range(len(self.hidden_channels)):
            x = self.convs[i](x, data.edge_index, data.edge_attr)
            x = self.act(x)
        x_gat_conv = x.reshape(batch, -1)
        
        #encoder
        x_gat_ln = torch.log(1 + x_gat_conv) #for data steability
        fea_out = self.encoder(x_gat_ln)

        #mu and variance
        z_loc = self.fc_mu(fea_out) # mean
        z_scale = torch.exp(self.fc_var(fea_out)) + 1e-4 #variance
        return z_loc, z_scale

class Fea_Decoder(nn.Module):
    def __init__(self, z_dim, hidden_dim, out_dim):
        super().__init__()
        # input layer
        hidden_dim_decoder = hidden_dim[:]
        hidden_dim_decoder.reverse()

        # layers
        modules = []
        for h_dim in hidden_dim_decoder:
            modules.append(
                nn.Sequential(
                    nn.Linear(z_dim, h_dim),
                    nn.BatchNorm1d(h_dim, momentum=0.01, eps=0.001), 
                    nn.LayerNorm(h_dim, elementwise_affine=False),
                    nn.ReLU(),
                    nn.Dropout(p=0))
            )
            z_dim = h_dim
        
        self.decoder = nn.Sequential(*modules)

        # final layer
        self.final = nn.Linear(hidden_dim_decoder[-1], out_dim)
        
    def forward(self, z):
        # define the forward computation on the latent z
        hidden = self.decoder(z)
        # return the parameter for the output distribution
        xrec= torch.sigmoid(self.final(hidden))
        return xrec


class Adj_Decoder(nn.Module):
    def __init__(self, z_dim, hidden_dim, out_dim):
        super().__init__()
        # input layer
        hidden_dim_decoder = hidden_dim[:]
        hidden_dim_decoder.reverse()

        # layers
        modules = []
        for h_dim in hidden_dim_decoder:
            modules.append(
                nn.Sequential(
                    nn.Linear(z_dim, h_dim),
                    nn.BatchNorm1d(h_dim, momentum=0.01, eps=0.001), 
                    nn.LayerNorm(h_dim, elementwise_affine=False),
                    nn.ReLU(),
                    nn.Dropout(p=0))
            )
            z_dim = h_dim
        
        self.decoder = nn.Sequential(*modules)

        # final layer
        self.final = nn.Linear(hidden_dim_decoder[-1], out_dim)
        
    def forward(self, z):
        # define the forward computation on the latent z
        hidden = self.decoder(z)
        # return the parameter for the output distribution
        xrec= torch.sigmoid(self.final(hidden))
        return xrec
    
############VAE
class VAE(nn.Module):
    def __init__(self, 
                 in_node_fea, gat_channel, nhead, gat_dropout, concat, 
                 in_fea, list_gra_enc_hid, z_dim, gra_p_dropout, 
                 list_fea_dec_hid, 
                 in_adj, list_adj_dec_hid,
                 kl_beta, fea_lambda, adj_lambda 
                 ):
        super().__init__()

        """
        === GAT Encoder
        in_node_fea: number of input feature per node
        gat_channel: number of feature per node after gat layer
        nhead: number of attention head
        gat_dropout -> float: gat layer dropout rate
        concat -> boolean: concatenate head or average
        in_fea: original feature vector size
        list_gra_enc_hid: feature encoder hidden dimension in list
        gra_p_dropout: feature encoder dropout

        === latent 
        z_dim: latent dismension

        === feature decoder
        list_fea_dec_hid: feature decoder hidden layer list

        === adjacency decoder
        in_adj: original adjacency vector size
        list_adj_dec_hid: feature decoder hidden layer list

        === lambdas
        kl_beta: factor scale the kl divergence
        fea_lambda: factor sacle feature reconstruction loss
        adj_lambda: factor sacle adjacency reconstruction loss
        """

        ### convolution graph encoder
        self.in_node_fea = in_node_fea
        self.gat_channel = gat_channel
        self.nhead = nhead
        self.gat_dropout = gat_dropout
        self.concat = concat

        self.in_fea = in_fea
        self.list_gra_enc_hid = list_gra_enc_hid
        self.gra_p_dropout = gra_p_dropout
        self.z_dim = z_dim
        self.gra_encoder = Gra_Encoder(self.in_node_fea, self.gat_channel, self.nhead, self.gat_dropout, self.concat, self.in_fea, self.list_gra_enc_hid, self.z_dim, self.gra_p_dropout)

        ### feature decoder
        self.list_fea_dec_hid = list_fea_dec_hid
        self.fea_decoder = Fea_Decoder(self.z_dim, self.list_fea_dec_hid, self.in_fea)
        #setup the parameters of the generative model
        self.fea_log_theta = torch.nn.Parameter(torch.randn(self.in_fea))
        self.fea_gate_logits = torch.nn.Parameter(torch.randn(self.in_fea))

        ### adjacency decoder
        self.in_adj = in_adj
        self.list_adj_dec_hid = list_adj_dec_hid
        self.adj_decoder = Adj_Decoder(self.z_dim, self.list_adj_dec_hid, self.in_adj)
        self.adj_log_theta = torch.nn.Parameter(torch.randn(self.in_adj))
        self.adj_gate_logits = torch.nn.Parameter(torch.randn(self.in_adj))

        ### lambdas
        self.kl_beta = kl_beta
        self.fea_lambda = fea_lambda
        self.adj_lambda = adj_lambda
        
    def model(self, x_gra): #the input for the model will be batch_size * x(x is the size of one cell sample)
        # register PyTorch module `decoder` with Pyro p(z) p(x|z), like nn.module
        pyro.module("fea_decoder", self) #give Decoder a name = "fea_decoder", let pyro knows all the parameters inside Decoder
        pyro.module("adj_decoder", self)
        with pyro.plate("data", x_gra.y.shape[0]): #use pyro.plate to declare the mini-batch size
            # setup hyperparameters(mu, variance) for prior p(z), it's actually constrained in guide
                # define size with torch.Size, but reture with same dtype as x and will be on the same device as x
                # x.shape[0] is the batch size
                #define batch size inside plate
            x_fea = x_gra.x_fea
            x_adj = x_gra.x_adj
            batch = x_fea.shape[0]
            z_loc = x_fea.new_zeros(torch.Size((batch, self.z_dim))) 
            z_scale = x_fea.new_ones(torch.Size((batch, self.z_dim)))

            # sample from prior (value will be sampled by guide when computing the ELBO)
            # pyro create a normal distribution as prior(q(z)) and then sample from this prior to get z
            # z_loc and z_scale is the set of parameters, each pair of parameter will produce a normal distribution,
            # z will sample from each of this distributions, so the size of z = batch_size * latent_space
            # .to_event(1) sample z from multivariate normal distribution instead from each univariate normal distribution
            with poutine.scale(scale = self.kl_beta):
                z = pyro.sample("latent", dist.Normal(z_loc, z_scale).to_event(1))

            ## Reconstruction Loss
            ## since it's single cell data, assume the parameters output from decoder is zero inflated negative binomial distribution
            ## we want to get the probability of p(x|z) (likelihood), measure the probability of seeing the output given the z that was sampled.
            
            #sample from z to get fea
            ## parameterization from zero inflated negative binomial
            # get the "normalized" mean of the negative binomial
            fea_px_scale = self.fea_decoder(z)
            # get the mean of the negative binomial
            fea_px_rate =  fea_px_scale
            # get the dispersion parameter
            fea_theta = torch.exp(self.fea_log_theta)
            fea_glog=self.fea_gate_logits

            # build count distribution
            fea_nb_logits = (fea_px_rate + 1e-4).log() - (fea_theta + 1e-4).log()
            fea_x_dist = dist.ZeroInflatedNegativeBinomial(total_count=fea_theta, logits=fea_nb_logits, gate_logits=fea_glog)
            with poutine.scale(scale = self.fea_lambda):
                # score against actual counts: reconstruction loss of feature matrix
                fea_rx=pyro.sample("obs_fea", fea_x_dist.to_event(1), obs=x_fea)

            #sample from z to get adj
            ## parameterization from zero inflated negative binomial
            # get the "normalized" mean of the negative binomial
            adj_px_scale = self.adj_decoder(z)
            # get the mean of the negative binomial
            adj_px_rate =  adj_px_scale
            # get the dispersion parameter
            adj_theta = torch.exp(self.adj_log_theta)
            adj_glog=self.adj_gate_logits

            # build count distribution
            adj_nb_logits = (adj_px_rate + 1e-4).log() - (adj_theta + 1e-4).log()
            adj_x_dist = dist.ZeroInflatedNegativeBinomial(total_count=adj_theta, logits=adj_nb_logits, gate_logits=adj_glog)
            with poutine.scale(scale = self.adj_lambda):
                # score against actual counts
                adj_rx=pyro.sample("obs2", adj_x_dist.to_event(1), obs=x_adj)

            return fea_rx, adj_rx
        
    def guide(self, x_gra):
        # define the guide (i.e. variational distribution) q(z|x)
        pyro.module("gra_encoder", self)
        with pyro.plate("data", x_gra.y.shape[0]):
            batch = x_gra.x_fea.shape[0]
            [qz_m, qz_v] = self.gra_encoder(x_gra, batch)

            # sample the latent code z given x
            with poutine.scale(scale = self.kl_beta):
                rz=pyro.sample("latent", dist.Normal(qz_m, qz_v.sqrt()).to_event(1))
            return rz
            
    def getZ(self, x_gra):
        batch = x_gra.x_fea.shape[0]
        [z_mu, z_var] = self.gra_encoder(x_gra, batch)
        # sample in latent space
        # z = dist.Normal(z_mu, z_var.sqrt()).sample()
        return z_mu, z_mu + z_var

def define_svi(in_node_fea, in_fea, in_adj, params, device, pretrain_path_fea = None, pretrain_path_adj = None):
    ### setup the VAE 
    vae = VAE(in_node_fea, params["gat_channel"], params["nhead"], params["gat_dropout"], params["concat"],
              in_fea, params["list_gra_enc_hid"], params["z_dim"], params["gra_p_dropout"], 
              params["list_fea_dec_hid"],
              in_adj, params["list_adj_dec_hid"],
              params["kl_beta"], params["fea_lambda"], params["adj_lambda"]
              )
    ###load pre_train model
    if params["pre_train_load"]:
        ### load adj pretrain model dict
        adj_pretrained_dict= torch.load(pretrain_path_adj)

        #### update current model
        ## use adj_pretrain to initialize gra_encoder and adj_deocder
        for name in vae.state_dict().keys():
            if name in adj_pretrained_dict["model_state_dict"].keys():
                print(name)
                vae.state_dict()[name].copy_(adj_pretrained_dict["model_state_dict"][name])

        ### load feature pretrain model dict_only load fea_decoder
        fea_pretrained_dict= torch.load(pretrain_path_fea)

        ### remove gra_encoder weight from feature
        new_fea_pretrained_dict = fea_pretrained_dict["model_state_dict"].copy()
        for key, value in fea_pretrained_dict["model_state_dict"].items():
            if "gra_encoder" in key:
                new_fea_pretrained_dict.pop(key)

        #### update current model
        ## use adj_pretrain to initialize gra_encoder and adj_deocder
        for name in vae.state_dict().keys():
            if name in new_fea_pretrained_dict.keys():
                print(name)
                vae.state_dict()[name].copy_(new_fea_pretrained_dict[name])

    # put the vae model in gpu
    vae.cuda(device)
    vae.train()

    # setup the optimizer
    adam_args = {"lr": params["lr"]}
    optimizer = Adam(adam_args) 

    # setup the inference algorithm
    elbo = Trace_ELBO(strict_enumeration_warning=False)
    svi = SVI(vae.model, vae.guide, optimizer, loss=elbo)

    return vae, svi
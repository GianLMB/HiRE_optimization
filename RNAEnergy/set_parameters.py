import numpy as np
import torch
import torch.nn as nn
from torch.nn import Parameter
import LocalEnergyVct as le
import os
import re
import pandas as pd
import pickle


class LocalEnergyOpt(nn.Module):

    def __init__(self, top_pars=None, dat_pars=None, device='cpu', set_to_one=False):
        super(LocalEnergyOpt, self).__init__()
        if top_pars is None:
            top_pars = {'bond_type': torch.zeros(15,2), 'angle_type': torch.zeros(13,2), 'torsion_type': torch.zeros(25,3)}
        if dat_pars is None:
            dat_pars = torch.zeros(47)
        if set_to_one:
            top_pars['bond_type'][:,0] = 1.
            top_pars['angle_type'][:,0] = 1.
            top_pars['torsion_type'][:,0] = 1.
        self.device = device
        self.opt_pars = Parameter(torch.tensor(dat_pars, dtype=torch.float, device=self.device, requires_grad=True))
        self.bond_type = Parameter(torch.tensor(top_pars['bond_type'], dtype=torch.float, device=self.device, requires_grad=True))
        self.angle_type = Parameter(torch.tensor(top_pars['angle_type'], dtype=torch.float, device=self.device, requires_grad=True))
        self.tor_type = Parameter(torch.tensor(top_pars['torsion_type'][:,(0,2)], dtype=torch.float, device=self.device, requires_grad=True))
        self.multiplicity = torch.tensor(top_pars['torsion_type'][:,1], dtype=torch.int64, device=self.device, requires_grad=False)


    def forward(self, X):

        X_lengths = X['lengths']
        X_features = X['features']

        if len(X_lengths.shape) == 1:
            X_lengths = X_lengths.unsqueeze(0)
            X_features = X_features.unsqueeze(0)

        energy = torch.zeros(X_lengths.shape[0],3,device=self.device)

        for i in range(X_lengths.shape[0]):
            lengths = X_lengths[i]
            features = X_features[i]
            atoms = features[:lengths[0],0].long()
            # res_labels
            # res_pointer
            # mass
            # charge
            coords = features[:lengths[5],5].view(-1,3)
            bonds = features[:lengths[6],6].long().view(-1,3)
            angles = features[:lengths[7],7].long().view(-1,4)
            tors = features[:lengths[8],8].long().view(-1,5)
            energy[i,0] = le.bonds_energy(coords,bonds,self.bond_type,self.opt_pars)
            energy[i,1] = le.angles_energy(atoms,coords,angles,self.angle_type,self.opt_pars)
            energy[i,2] = le.torsions_energy(atoms,coords,tors,self.tor_type,self.multiplicity,self.opt_pars)

        return energy


def make_dat(par_filepath,overwrite=False):   # model parameters saved as .pth file

    model = LocalEnergyOpt()
    model.load_state_dict((torch.load(par_filepath)))
    varied_par_indexes = [*range(0,10), *range(28,38), *[46]]
    model_parameters = [p.data for p in model.parameters()]
    new_values = np.zeros(47)
    new_values[varied_par_indexes] = model_parameters[0][varied_par_indexes]

    with open('scale_RNA.dat','r') as f:
        data = np.asarray([line.strip().split(maxsplit=2) for line in f])
    data[:,1] = [f'{v:.3f}' for v in new_values]

    if overwrite:
        out_file = 'scale_RNA.dat'
    else:
        out_file = 'scale_RNAnew.dat'

    # np.savetxt(out_file,data,fmt="%4s %9s %50s")

    col_format = "{:>4}" + "{:>9}" + "         " + "{:<60}" + "\n"
    with open(out_file, 'w') as of:
        for x in data:
            of.write(col_format.format(*x))

    return


def make_top(par_filepath,overwrite=False):   # model parameters saved as .pth file

    model = LocalEnergyOpt()
    model.load_state_dict((torch.load(par_filepath)))
    model_parameters = [p.detach().numpy() for p in model.parameters()]

    with open('parametres_RNA.top', 'r') as f:
        print(f.readlines())
    return



make_top("data/Results/400_same_order.pth")

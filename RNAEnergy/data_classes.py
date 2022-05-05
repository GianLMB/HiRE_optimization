from __future__ import print_function, division
import torch
import pickle
import os
import pandas as pd
import numpy as np
from torch.utils.data import Dataset
import torch.nn as nn
from torch.nn import Parameter
import LocalEnergyVct as le

# Ignore warnings
import warnings
warnings.filterwarnings("ignore")


class ToTensor(object):
    """Convert ndarrays in sample to Tensors."""

    def __call__(self, sample):
        sample['coords'] = torch.from_numpy(sample['coords'])
        sample['mass'] = torch.tensor(sample['mass'])
        sample['charge'] = torch.tensor(sample['charge'])
        sample['energy'] = torch.tensor(sample['energy'])
        return sample


class RNASeqDataset_local(Dataset):
    """RNA sequences dataset."""

    def __init__(self, pkl_file, root_dir, transform=None):
        self.seq_list = pickle.load(open(pkl_file, 'rb'))
        self.root_dir = root_dir
        self.transform = transform

    def __len__(self):
        return len(self.seq_list)

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()

        sample = self.seq_list[idx]

        if self.transform:
            sample = self.transform(sample)

        return sample




class RNASeqDataset(Dataset):
    """RNA sequences dataset."""

    def __init__(self, device='cpu', csv_file='data/SeqCSV/seq_frame.csv', root_dir='data/SeqCSV/', transform=None):
        self.seq_frame = pd.read_csv(csv_file)
        self.root_dir = root_dir
        # self.transform = transform
        size = len(self.seq_frame)
        lengths = self.seq_frame.iloc[:, 1:].astype('int64')
        lengths = torch.from_numpy(np.array(lengths)).to(device)

        # get features size
        seq_name = os.path.join(self.root_dir, self.seq_frame.iloc[0, 0] + '.csv')
        features = pd.read_csv(seq_name)
        row, col = np.array(features).shape

        features = torch.zeros(size,row,col)
        for i in range(size):
            seq_name = os.path.join(self.root_dir, self.seq_frame.iloc[i, 0] + '.csv')
            seq = pd.read_csv(seq_name)
            features[i,:,:] = torch.from_numpy(np.array(seq))
        features = features.to(device)
        self.dataset = {'lengths': lengths, 'features': features}

    def __len__(self):
        return len(self.seq_frame)

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
        lengths = self.dataset['lengths'][idx]
        features = self.dataset['features'][idx]
        sample = {'lengths': lengths, 'features': features}

        return sample




class LocalEnergyOpt(nn.Module):

    def __init__(self, top_pars=None, dat_pars=None, device='cpu', set_to_one=False):
        super(LocalEnergyOpt, self).__init__()
        if top_pars is None:
            top_pars = {'bond_type': torch.zeros(15,2), 'angle_type': torch.zeros(13,2), 'torsion_type': torch.zeros(25,3)}
            top_pars['torsion_type'][:,1] = torch.tensor([*([1]*19), *[3,1,1,2,3,2]])
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


    def forward(self,X):

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

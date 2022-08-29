# Contains classes for dataset and model with less degrees of freedom (i.e. with global particle-dependent 
# parameters extracted from scaleRNA.dat file set to 1)

from __future__ import print_function, division
import torch
import os
import pandas as pd
import numpy as np
from torch.utils.data import Dataset
import torch.nn as nn
from torch.nn import Parameter
import src.energy_functions as le

# Ignore warnings
import warnings
warnings.filterwarnings("ignore")


class RNADataset(Dataset):
    """
    RNA sequences dataset, class inherited from torch.utils.data.Dataset.
    Data from sequences are extracted and encoded in a dictionary that can be retrieved, after initialisation, in the dataset attribute. 
    It contains two tensors:
        - 'lengths', which contain the original length of the features for each RNA sequence;
        - 'features', a padded tensor with all the features (particles type, type of interactions, index of 
             particles involved ...) for each sequence.   
    """

    def __init__(self, device='cpu', csv_file='../data/CSV_minimized/seq_frame.csv', root_dir='../data/CSV_minimized/'):
        """
        device: 'cuda' or 'cpu', indicates where the dataset is allocated. The default is 'cpu'.
        csv_file (string): Path to the csv file with annotations. The default is '../data/CSV_minimized/seq_frame.csv'.
        root_dir (string): Directory with all the features of the RNA sequences. The default is '../data/CSV_minimized/'.
        """
        
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

        self.names = []
        features = torch.zeros(size, row, col)
        for i in range(size):
            self.names.append(self.seq_frame.iloc[i, 0])
            seq_name = os.path.join(self.root_dir, self.seq_frame.iloc[i, 0] + '.csv')
            seq = pd.read_csv(seq_name)
            features[i, :, :] = torch.from_numpy(np.array(seq))
        features = features.to(device)
        self.dataset = {'lengths': lengths, 'features': features}
        # print(lengths.shape, features.shape)

    def __len__(self):
        return len(self.seq_frame)

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
        lengths = self.dataset['lengths'][idx]
        features = self.dataset['features'][idx]
        name = self.names[idx]
        sample = {'lengths': lengths, 'features': features, 'name': name}

        return sample


class Model(nn.Module):
    """
    Creates the Model for the optimization process, setting some parameters to 1. The others are divided in
    - globals, for global coefficients
    - couplings, for particle dependent couplings
    - equil, containing geometric parameters
    The Model computes the energy for the three types of local interactions.
    """

    def __init__(self, top_pars=None, dat_pars=None, device='cpu'):
        """
        top_pars (string): dictionary containing parameters extracted from TOP files. The default is None.
        dat_pars (string): dictionary containing parameters extracted from DAT file. The default is None.
        device: 'cuda' or 'cpu', indicates where the Model is allocated. The default is 'cpu'.
        """
        
        super(Model, self).__init__()
        if top_pars is None:
            top_pars = {'bond_type': torch.zeros(15,2), 'angle_type': torch.zeros(13,2), 'torsion_type': torch.zeros(25,3)}
            top_pars['torsion_type'][:,1] = torch.tensor([*([1]*19), *[3,1,1,2,3,2]])
        if dat_pars is None:
            dat_pars = torch.zeros(3)
        else:
            dat_pars = [dat_pars[0], dat_pars[9], dat_pars[-1]]

        self.device = device
        self.couplings = self.Interactions()
        self.equil = self.Interactions()
        self.globals = Parameter(torch.tensor(dat_pars, dtype=torch.float, device=self.device, requires_grad=True))
        self.couplings.bonds = Parameter(torch.tensor(top_pars['bond_type'][:,0], dtype=torch.float, device=self.device, requires_grad=True))
        self.equil.bonds = Parameter(torch.tensor(top_pars['bond_type'][:,1], dtype=torch.float, device=self.device, requires_grad=False))
        self.couplings.angles = Parameter(torch.tensor(top_pars['angle_type'][:,0], dtype=torch.float, device=self.device, requires_grad=True))
        self.equil.angles = Parameter(torch.tensor(top_pars['angle_type'][:,1], dtype=torch.float, device=self.device, requires_grad=False))
        self.couplings.torsions = Parameter(torch.tensor(top_pars['torsion_type'][:,0], dtype=torch.float, device=self.device, requires_grad=True))
        self.equil.torsions = Parameter(torch.tensor(top_pars['torsion_type'][:,2], dtype=torch.float, device=self.device, requires_grad=False))
        self.multiplicity = torch.tensor(top_pars['torsion_type'][:,1], dtype=torch.int64, device=self.device, requires_grad=False)

    class Interactions(nn.Module):
        pass

    def forward(self, X):

        X_lengths = X['lengths']
        X_features = X['features']

        if len(X_lengths.shape) == 1:
            X_lengths = X_lengths.unsqueeze(0)
            X_features = X_features.unsqueeze(0)

        energy = torch.zeros(X_lengths.shape[0], 3, device=self.device)

        for i in range(X_lengths.shape[0]):
            lengths = X_lengths[i]
            features = X_features[i]
            # atoms = features[:lengths[0], 0].long()
            # res_labels
            # res_pointer
            # mass
            # charge
            coords = features[:lengths[5], 5].view(-1, 3)
            bonds = features[:lengths[6], 6].long().view(-1, 3)
            angles = features[:lengths[7], 7].long().view(-1, 4)
            tors = features[:lengths[8], 8].long().view(-1, 5)
            energy[i,0] = le.bonds_energy(coords, bonds, self.couplings, self.equil, self.globals)
            energy[i,1] = le.angles_energy(coords, angles, self.couplings, self.equil,  self.globals)
            energy[i,2] = le.torsions_energy(coords, tors, self.couplings, self.equil,  self.multiplicity, self.globals)

        return energy 

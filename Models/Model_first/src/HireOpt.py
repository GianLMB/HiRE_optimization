# Deprecated: old version of the optimization model. In newer versions, functions and classes 
# are organised in different files, and scripts to execute are in "../script/" folder

from __future__ import print_function, division
import os
import pickle
import pandas as pd
import torch
from torch.utils.data import Dataset, DataLoader, TensorDataset
import numpy as np
import matplotlib.pyplot as plt
import torch.nn as nn
from torch.nn import Parameter
from numba import cuda
print(torch.cuda.is_available())
import LocalEnergy as le

# ignore warnings
import warnings
warnings.filterwarnings("ignore")

device = torch.device("cuda" if torch.cuda.is_available() else 'cpu')


def get_target(X):
    if len(X['features'].shape) == 2:
        X['features'] = X['features'].unsqueeze(0)
    # print(torch.sum(X['features'][:,0:3,9],dim=1))
    target = (torch.sum(X['features'][:,0:3,9],dim=1)/X['features'].shape[0]).squeeze().to(device)
    return target


def loss_fn(energy,target,lc=1):
    batch_size = energy.shape[0]
    # jac = torch.autograd.functional.jacobian(energy, pars, create_graph=True)
    loss = (energy - target).pow(2).sum() / batch_size  #  + lc * grad.pow(2).sum())
    return loss


class ToTensor(object):
    """Convert ndarrays in sample to Tensors."""

    def __call__(self, sample):
        lengths, features = sample['lengths'], sample['features']
        return {'lengths': torch.from_numpy(lengths),
                'features': torch.from_numpy(features)}


class RNASeqDataset(Dataset):
    """RNA sequences dataset."""

    def __init__(self, csv_file='data/SeqCSV/seq_frame.csv', root_dir='data/SeqCSV/', transform=None):
        self.seq_frame = pd.read_csv(csv_file)
        self.root_dir = root_dir
        self.transform = transform
        # TODO: rewrite init and getitem methods to load all dataset into gpu

    def __len__(self):
        return len(self.seq_frame)

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()

        name = self.seq_frame.iloc[idx, 0]
        seq_name = os.path.join(self.root_dir, name+'.csv')
        lengths = self.seq_frame.iloc[idx, 1:].astype('int64')
        lengths = np.array(lengths)
        features = pd.read_csv(seq_name)
        features = np.array(features)
        sample = {'lengths': lengths, 'features': features}

        if self.transform:
            sample = self.transform(sample)

        return sample


class LocalEnergyOpt(nn.Module):

    def __init__(self,fixed_pars,opt_pars):
        super(LocalEnergyOpt, self).__init__()
        self.opt_pars = Parameter(torch.tensor(opt_pars, dtype=torch.float, device=device, requires_grad=True))
        self.bond_type = Parameter(torch.tensor(fixed_pars['bond_type'], dtype=torch.float, device=device, requires_grad=False))
        self.angle_type = Parameter(torch.tensor(fixed_pars['angle_type'], dtype=torch.float, device=device, requires_grad=False))
        self.tor_type = Parameter(torch.tensor(fixed_pars['torsion_type'], dtype=torch.float, device=device, requires_grad=False))

    def forward(self,X):

        X_lengths = X['lengths']
        X_features = X['features']

        if len(X_lengths.shape) == 1:
            X_lengths = X_lengths.unsqueeze(0)
            X_features = X_features.unsqueeze(0)

        energy = torch.zeros(X_lengths.shape[0]).to(device)
        for i in range(X_lengths.shape[0]):
            lengths = X_lengths[i]
            features = X_features[i]
            if torch.is_tensor(lengths):
                lengths = lengths.tolist()
            atoms = np.array(features[:lengths[0],0],dtype=int)
            # res_labels
            # res_pointer
            # mass
            # charge
            coords = features[:lengths[5],5].reshape(-1,3)
            bonds = np.array(features[:lengths[6],6],dtype=int).reshape(-1,3)
            angles = np.array(features[:lengths[7],7],dtype=int).reshape(-1,4)
            tors = np.array(features[:lengths[8],8],dtype=int).reshape(-1,5)  # all indexes: not necessary to convert to tensors
            ebonds = le.bonds_energy(coords,bonds,self.bond_type,self.opt_pars)
            eangles = le.angles_energy(atoms,coords,angles,self.angle_type,self.opt_pars)
            etors = le.torsions_energy(atoms,coords,tors,self.tor_type,self.opt_pars)
            energy[i] += ebonds + eangles + etors

        return energy


seq_data = RNASeqDataset(transform=ToTensor())
print('dataset created')

for X in seq_data:
    X['lengths'] = X['lengths'].to(device)
    X['features'] = X['features'].to(device)
print('data loaded to gpu')

batch_size = 32
seq_dataloader = DataLoader(seq_data,batch_size=batch_size,shuffle=True,num_workers=1,pin_memory=True)  # all fields must have same length

fixed_pars = pickle.load(open('data/SeqCSV/fixed_pars.p', 'rb'))
opt_pars = pickle.load(open('data/SeqCSV/pars.p', 'rb'))
model = LocalEnergyOpt(fixed_pars,opt_pars).to(device)

'''''
# check if initial energies were not NaN -> ok
for batch, X in enumerate(seq_dataloader):
    # print(X['name'])
    pred = model(X)
    print(pred)
    target = get_target(X)
    print(target)
'''''

lr = 1e-8
optimizer = torch.optim.SGD(model.parameters(), lr=lr)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', factor = 0.5, patience = 500, cooldown = 1000, threshold = 1e-8, verbose = True)

epochs = 5
num_batches = len(seq_dataloader)
loss_epoch = []
for index_epoch in range(epochs):
    # print(index_epoch)
    train_loss = 0

    for X in seq_dataloader:
        optimizer.zero_grad()
        energy = model(X)
        # der = energy.backward(create_graph=True)
        target = get_target(X)
        loss = loss_fn(energy,target)
        # loss += torch.sum(der)
        train_loss += loss.item()/num_batches
        print(loss.item())
        loss.backward()
        optimizer.step()

    print(index_epoch, train_loss)
    print('\n \n \n \n \n')
    loss_epoch.append(train_loss)




'''''
while 1:
    for i in range(len(dataset)):
        X = dataset[i]
        # print(X['name'])
        if len(X['energy']) == 0:
            continue
        if X['name'] == '2ho7':
            continue
        energy = model(X)
        target = sum(X['energy'][0:4])
        loss = loss_fn(energy,target)
        loss.backward()
        if i == 1:
            print(loss.item())
            print(model.opt_pars.grad)

        with torch.no_grad():
            model.opt_pars -= lr*model.opt_pars.grad
'''''

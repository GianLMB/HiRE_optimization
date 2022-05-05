from __future__ import print_function, division
import os
import pickle
import pandas as pd
import torch
import time
from torch.utils.data import Dataset, DataLoader
import numpy as np
import matplotlib.pyplot as plt
import torch.nn as nn
from torch.nn import Parameter
from numba import cuda
print(torch.cuda.is_available())
import LocalEnergyVct as le

# ignore warnings
import warnings
warnings.filterwarnings("ignore")

# import torch.multiprocessing as mp
# mp.set_start_method('spawn')


device = torch.device("cuda" if torch.cuda.is_available() else 'cpu')


def get_target(X):
    if len(X['features'].shape) == 2:
        X['features'] = X['features'].unsqueeze(0)
    # print(torch.sum(X['features'][:,0:3,9],dim=1))
    target = (X['features'][:,0:3,9]).to(device)  # /X['features'].shape[0]).squeeze().to(device)
    return target


def loss_fn(energy,target):
    batch_size = energy.shape[0]
    loss = (energy - target).pow(2).sum() / batch_size
    return loss


def loss_with_grad(pred,target,model,lc=0.1):
    batch_size = pred.shape[0]
    grad2 = 0.
    for en in pred.view(-1,):
        grad_list = torch.autograd.grad(en, model.parameters(), create_graph=True)
        for t in grad_list:
            grad2 += t.pow(2).sum().squeeze()
    # print((pred - target).pow(2).sum(), lc*grad2)
    loss = ((pred - target).pow(2).sum() + lc*grad2) / batch_size
    return loss


def train(dataloader, model, loss_fn, optimizer):
    # size = len(dataloader.dataset)
    # num_batches = len(dataloader)
    # model.train()
    num_batches = 0
    train_loss = 0

    for X in dataloader:
        # Compute prediction error
        pred = model(X)
        target = get_target(X)
        loss = loss_fn(pred, target)
        if torch.isnan(loss):
            continue
        num_batches += 1
        train_loss += loss.item()

        # Backpropagation
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    train_loss /= num_batches
    print(f'Avg loss = {train_loss:>0.4f}, valid batches = {num_batches}')

    return train_loss



class ToTensor(object):
    """Convert ndarrays in sample to Tensors."""

    def __call__(self, sample):
        lengths, features = sample['lengths'], sample['features']
        return {'lengths': torch.from_numpy(lengths),
                'features': torch.from_numpy(features)}


class RNASeqDataset(Dataset):
    """RNA sequences dataset."""

    def __init__(self, device, csv_file='data/SeqCSV/seq_frame.csv', root_dir='data/SeqCSV/', transform=None):
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

    def __init__(self,fixed_pars,opt_pars):
        super(LocalEnergyOpt, self).__init__()
        self.opt_pars = Parameter(torch.tensor(opt_pars, dtype=torch.float, device=device, requires_grad=True))
        self.bond_type = Parameter(torch.tensor(fixed_pars['bond_type'], dtype=torch.float, device=device, requires_grad=True))
        self.angle_type = Parameter(torch.tensor(fixed_pars['angle_type'], dtype=torch.float, device=device, requires_grad=True))
        self.tor_type = Parameter(torch.tensor(fixed_pars['torsion_type'], dtype=torch.float, device=device, requires_grad=True))

    def forward(self,X):

        X_lengths = X['lengths']
        X_features = X['features']

        if len(X_lengths.shape) == 1:
            X_lengths = X_lengths.unsqueeze(0)
            X_features = X_features.unsqueeze(0)

        energy = torch.zeros(X_lengths.shape[0],3).to(device)

        for i in range(X_lengths.shape[0]):
            lengths = X_lengths[i]
            features = X_features[i]
            if torch.is_tensor(lengths):
                lengths = lengths.tolist()
            atoms = features[:lengths[0],0].long()
            # res_labels
            # res_pointer
            # mass
            # charge
            coords = features[:lengths[5],5].view(-1,3)
            bonds = features[:lengths[6],6].long().view(-1,3)
            angles = features[:lengths[7],7].long().view(-1,4)
            tors = features[:lengths[8],8].long().view(-1,5)  # all indexes: not necessary to convert to tensors
            energy[i,0] = le.bonds_energy(coords,bonds,self.bond_type,self.opt_pars)
            energy[i,1] = le.angles_energy(atoms,coords,angles,self.angle_type,self.opt_pars)
            energy[i,2] = le.torsions_energy(atoms,coords,tors,self.tor_type,self.opt_pars)

        return energy


seq_data = RNASeqDataset(device=device)
print(f'dataset allocated on {device}')
print(len(seq_data))
print(seq_data[0])

batch_size = 32
seq_dataloader = DataLoader(seq_data,batch_size=batch_size,shuffle=True) #,num_workers=1,pin_memory=True)
# num_batches = len(seq_dataloader)
# print(num_batches)

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

lr = 1e-7
optimizer = torch.optim.SGD(model.parameters(), lr=lr)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', factor = 0.5, patience = 500, cooldown = 1000, threshold = 1e-12, verbose = True)
loss_fn = loss_fn


epochs = 100
loss_epoch = []
for index_epoch in range(epochs):
    print(f'epoch {index_epoch+1}/{epochs} \n-------------------------')
    t0 = time.time()
    train_loss = train(seq_dataloader, model, loss_fn, optimizer)
    tf = time.time()
    print(f'time for epoch: {tf-t0} \n')
    loss_epoch.append(train_loss)

for p in model.parameters():
    print(p.data)
torch.save(model.state_dict(), 'data/Results/try.npy')
plt.plot(loss_epoch)
plt.show()

torch.cuda.empty_cache()


# TODO: analyze speed of convergence varying batch size and learning rate
'''
batch_size_list = (8,16,32,64,128,256,512,len(seq_data))
lr_list = np.logspace(-9,-7,3)
epochs = 5
res = np.zeros((len(batch_size_list),len(lr_list),2))
for i, batch_size in enumerate(batch_size_list):
    seq_dataloader = DataLoader(seq_data,batch_size=batch_size,shuffle=True,num_workers=1,pin_memory=True)
    for j, lr in enumerate(lr_list):
        model = LocalEnergyOpt(fixed_pars,opt_pars).to(device)
        optimizer = torch.optim.SGD(model.parameters(), lr=lr)
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', factor = 0.5, patience = 500, cooldown = 1000, threshold = 1e-12, verbose = True)
        t0 = time.time()
        print(f'Batch size = {batch_size}, Learning rate = {lr}\n')
        for index_epoch in range(epochs):
            print(f'epoch {index_epoch+1}/{epochs} \n-------------------------')
            train_loss = train(seq_dataloader, model, loss_fn, optimizer)
        tf = time.time()
        res[i][j][0] = train_loss
        res[i][j][1] = tf - t0
        print(f'\nCurrent values for res: \n{res[i][j]} \n\n')
'''

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

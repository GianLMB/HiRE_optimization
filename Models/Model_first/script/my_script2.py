# script to perform the optimization of the prameters, second version:
# the model is defined through the class LocalEnergy, and the optimizer assigns different
# learning rates to different groups of parameters.
# All the functions for loss, train and test are also defined here

from __future__ import print_function, division
import sys
sys.path.append('..')

import time
import torch
import pickle
from tqdm import tqdm
from torch.utils.data import DataLoader, random_split
from src.data_classes import RNASeqDataset, LocalEnergy
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt


# ignore warnings
import warnings
warnings.filterwarnings("ignore")





# FUNCTIONS
# -------------------------------------------------------------------------------------------------

def get_target(X):
    if len(X['features'].shape) == 2:
        X['features'] = X['features'].unsqueeze(0)
    target = (X['features'][:,0:3,9]).squeeze()  # .sum(dim=0).squeeze() / X['features'].shape[0]
    return target


# Functions without gradients in the loss

def loss_fn(energy,target):
    batch_size = energy.shape[0]
    loss = (energy - target).pow(2).sum() / batch_size
    return loss

def bonds_loss_fn(energy,target):
    batch_size = energy.shape[0]
    loss = (energy - target).pow(2).sum(dim=0)[0] / batch_size
    return loss


def angles_loss_fn(energy,target):
    batch_size = energy.shape[0]
    loss = (energy - target).pow(2).sum(dim=0)[1] / batch_size
    return loss


def torsions_loss_fn(energy,target):
    batch_size = energy.shape[0]
    loss = (energy - target).pow(2).sum(dim=0)[2] / batch_size
    return loss
    
def train(dataloader, model, loss_fn, optimizer):
    num_batches = len(dataloader)
    model.train()
    train_loss = 0
    for X in tqdm(dataloader):
        pred = model(X)
        target = get_target(X)
        loss = loss_fn(pred, target)
        train_loss += loss.item()
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    train_loss /= num_batches
    print(f'Avg loss = {train_loss:>0.4f}, batches = {num_batches}')
    return train_loss


def test(dataloader, model, loss_fn):
    num_batches = len(dataloader)
    model.eval()
    test_loss = 0
    with torch.no_grad():
        for X in tqdm(dataloader):
            pred = model(X)
            target = get_target(X)
            loss = loss_fn(pred, target)
            test_loss += loss.item()
    test_loss /= num_batches
    print(f'Avg test_loss = {test_loss:>0.4f}, batches = {num_batches}')
    return test_loss


# Functions with gradients in the loss implemented

def loss_fn_with_grad(pred,target,model,lc=0.05):
    batch_size = pred.shape[0]
    grad2 = 0.
    pred_split = torch.split(pred.view(-1,),1)
    grad_list = torch.autograd.grad(pred_split, model.parameters(), create_graph=True)
    for t in grad_list:
        grad2 += t.pow(2).sum().squeeze()
    # print((pred - target).pow(2).sum(), lc*grad2)
    loss = ((pred - target).pow(2).sum() + lc*grad2) / batch_size
    return loss


def train_with_grad(dataloader, model, loss_fn, optimizer):
    num_batches = len(dataloader)
    model.train()
    train_loss = 0
    for X in tqdm(dataloader):
        pred = model(X)
        target = get_target(X)
        loss = loss_fn(pred, target, model)
        train_loss += loss.item()
        optimizer.zero_grad()
        loss.backward(retain_graph=True)
        optimizer.step()
    train_loss /= num_batches
    print(f'Avg loss = {train_loss:>0.4f}, batches = {num_batches}')
    return train_loss


def test_with_grad(dataloader, model, loss_fn):
    num_batches = len(dataloader)
    model.eval()
    test_loss = 0
    for X in tqdm(dataloader):
        pred = model(X)
        target = get_target(X)
        loss = loss_fn(pred, target, model)
        test_loss += loss
    test_loss /= num_batches
    print(f'Avg test_loss = {test_loss:>0.4f}, batches = {num_batches}')
    return test_loss


def bonds_loss_fn(energy,target):
    batch_size = energy.shape[0]
    loss = (energy - target).pow(2).sum(dim=0)[0] / batch_size
    return loss


def angles_loss_fn(energy,target):
    batch_size = energy.shape[0]
    loss = (energy - target).pow(2).sum(dim=0)[1] / batch_size
    return loss


def torsions_loss_fn(energy,target):
    batch_size = energy.shape[0]
    loss = (energy - target).pow(2).sum(dim=0)[2] / batch_size
    return loss
# -------------------------------------------------------------------------------------------------


if __name__ == '__main__':

    # CUDA for Pytorch
    # print(torch.cuda.is_available())
    # device = torch.device("cuda" if torch.cuda.is_available() else 'cpu')
    device = 'cpu'


    # Parameters
    params = {'batch_size': 4,
          'shuffle': True,
          'num_workers': 0,
          'drop_last': True,
          'pin_memory': False}
    epochs = 400


    # Datasets and Dataloaders
    seq_data = RNASeqDataset(device=device)
    print(f'dataset allocated on {device}')

    tot_length = len(seq_data)
    test_length = int(0.2*tot_length)
    train_set, test_set = random_split(seq_data, [tot_length - test_length, test_length], generator=torch.Generator().manual_seed(42))
    print(f'Training set: {len(train_set)} elements')
    print(f'Test set: {len(test_set)} elements')

    train_dataloader = DataLoader(train_set,**params)
    test_dataloader = DataLoader(test_set,**params)


    # Model and optimizer
    # top_pars = pickle.load(open('data/SeqCSV/fixed_pars.p', 'rb'))
    # dat_pars = pickle.load(open('data/SeqCSV/pars.p', 'rb'))
    model = LocalEnergy().to(device)
    # torch.save(model.state_dict(), 'data/NewResults/initial_values.pth')
    model.load_state_dict(torch.load('../results/NewResults/initial_values.pth'))
    # print(model.state_dict())

    # """
    sd = model.state_dict()
    sd['couplings.angles'][4] = 1200
    sd['couplings.angles'][5] = 1200
    # sd['couplings.torsions'][:] = 10
    # sd["opt_pars"][28:38] = 10
    # sd["opt_pars"][1:9] = 10
    print(sd)

    model.load_state_dict(sd)
    torch.save(model.state_dict(), '../results/NewResults/AnglesResults/initial_big_purines_new.pth')
    # """

    optimizer = torch.optim.Adam([
        {'params': model.opt_pars, 'lr': 1e-4},
        {'params': model.couplings.bonds, 'lr': 1e-3},
        {'params': model.couplings.angles, 'lr': 1e-3},
        {'params': model.couplings.torsions, 'lr': 1e-5},
        {'params': model.equil.bonds, 'lr': 1e-7},
        {'params': model.equil.angles, 'lr': 1e-5},
        {'params': model.equil.torsions, 'lr': 1e-7}
    ])
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', factor=0.5, patience=10, cooldown=1000, threshold=1e-9, verbose=True)
    my_loss = angles_loss_fn  # _with_grad
    my_train = train  # _with_grad
    my_test = test  # _with_grad


    # Loop over epochs
    train_loss = []
    test_loss = []
    for index_epoch in range(epochs):
        print(f'epoch {index_epoch+1}/{epochs} \n-------------------------')
        t0 = time.time()
        train_tmp = my_train(train_dataloader, model, my_loss, optimizer)
        test_tmp = my_test(test_dataloader, model, my_loss)
        tf = time.time()
        print(f'time for epoch: {tf-t0} \n')
        # for p in model.parameters():
        #     print(p.grad)
        train_loss.append(train_tmp)
        test_loss.append(test_tmp)
        print(model.state_dict())

    for p in model.parameters():
        print(p.data)

    torch.save(model.state_dict(), '../results/NewResults/AnglesResults/final_big_purines_new.pth')

    plt.plot(train_loss)
    plt.plot(test_loss)

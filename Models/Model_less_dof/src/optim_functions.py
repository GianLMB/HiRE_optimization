# Functions for optimization process: loss function, training and test
# Implemented with or without considering gradients in the loss function

from __future__ import print_function, division
import torch
from tqdm import tqdm

# ignore warnings
import warnings
warnings.filterwarnings("ignore")


# FUNCTIONS
# -------------------------------------------------------------------------------------------------

# extract AMBER energies for sequences in a minibatch X
def get_target(X):
    if len(X['features'].shape) == 2:
        X['features'] = X['features'].unsqueeze(0)
    target = (X['features'][:,0:3,9]).squeeze()  # .sum(dim=0).squeeze() / X['features'].shape[0]
    return target


# Functions without gradients in the loss

# Loss function: sum of squared differences between energies obtained with AMBER and HiRE model
def loss_fn(energy,target):
    batch_size = energy.shape[0]
    loss = (energy - target).pow(2).sum() / batch_size
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

def loss_with_grad(pred,target,model,lc=0.05):
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

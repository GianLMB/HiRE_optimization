from __future__ import print_function, division
import sys
sys.path.append('..')

import time
import torch
import pickle
from tqdm import tqdm
from torch.utils.data import DataLoader, random_split
from src.data_classes import RNASeqDataset, LocalEnergy
import src.preconditioned_stochastic_gradient_descent as psgd
from sklearn.metrics import r2_score


# ignore warnings
import warnings
warnings.filterwarnings("ignore")


# FUNCTIONS
# -------------------------------------------------------------------------------------------------

def get_target(X):
    if len(X['features'].shape) == 2:
        X['features'] = X['features'].unsqueeze(0)
    target = X['features'][:,0:3,9]  # .sum(dim=0).squeeze() / X['features'].shape[0]
    return target


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

# --------------------------------------------------------------------------------------


if __name__ == '__main__':

    # CUDA for Pytorch
    # print(torch.cuda.is_available())
    # device = torch.device("cuda" if torch.cuda.is_available() else 'cpu')
    device = 'cpu'


    # Parameters
    params = {'batch_size': 1,
          'shuffle': True,
          'num_workers': 0,
          'drop_last': True,
          'pin_memory': False}
    epochs = 100

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
    model = LocalEnergy(device=device).to(device)
    # torch.save(model.state_dict(), 'data/NewResults/initial_values.pth')
    model.load_state_dict(torch.load('data/NewResults/initial_values_samedist.pth'))


    """
    sd = model.state_dict()
    print(sd)
    sd['couplings.angles'] /= 10
    sd['couplings.torsions'] *= 10
    model.load_state_dict(sd)
    print(model.state_dict())
    torch.save(model.state_dict(), 'data/NewResults/initial_values_samedist.pth')
    """

    # scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', factor=0.5, patience=50, cooldown=1000, threshold=1e-12, verbose=True)
    my_loss = loss_fn


    # Loop over epochs

    train_loss = []
    test_loss = []
    train_batches = len(train_dataloader)
    test_batches = len(test_dataloader)
    num_para = sum([torch.numel(w) for w in model.parameters()])
    Q = 0.1 * torch.eye(num_para).to(device)
    lr_list = [1e-4, 0.1, 1e-2, 1e-3, 1e-6, 1e-6, 1e-7]

    for index_epoch in range(epochs):
        print(f'epoch {index_epoch+1}/{epochs} \n-------------------------')
        t0 = time.time()
        train_lossel = 0
        test_lossel = 0

        for X in tqdm(train_dataloader):
            pred = model(X)
            target = get_target(X)
            loss = my_loss(pred, target)
            train_lossel += loss.item()
            grads = torch.autograd.grad(loss, model.parameters(), create_graph=True)
            vs = [torch.randn_like(w) for w in model.parameters()]
            for (v,g) in zip(vs,grads):
                v[g==0] = 0
            Hvs = torch.autograd.grad(grads, model.parameters(), vs)
            with torch.no_grad():
                Q = psgd.update_precond_dense(Q, vs, Hvs, step=0.1)
                pre_grads = psgd.precond_grad_dense(Q, grads)
                [w.subtract_(lr * g) for (w, g, lr) in zip(model.parameters(), pre_grads, lr_list)]

        train_lossel /= train_batches
        train_loss.append(train_lossel)

        print(model.state_dict())
        print("Avg train loss = {}, batches = {}".format(train_lossel,train_batches))

        with torch.no_grad():
            for X in tqdm(test_dataloader):
                pred = model(X)
                target = get_target(X)
                loss = my_loss(pred, target)
                test_lossel += loss.item()

        test_lossel /= test_batches
        test_loss.append(test_lossel)
        print("Avg test loss = {}, batches = {}".format(test_lossel, test_batches))

        tf = time.time()
        print(f'time for epoch: {tf - t0} \n')

        """
        for p in model.parameters():
            print(p)
        """

    torch.save(model.state_dict(), 'data/NewResults/torsions_150_b1_e3e2e6_psgd_samedist.pth')

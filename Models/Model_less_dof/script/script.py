# script to perform optimization of the parameters
# all hyperparameters for training process can be modified, and new parameters
# and loss values are stored at the end of the training

from __future__ import print_function, division
import sys
sys.path.append('..')

import numpy as np
import time
import torch
import pickle
from torch.utils.data import DataLoader, random_split
from src.data_classes import RNADataset, Model
from src.optim_functions import train, test, loss_fn
# from sklearn.metrics import r2_score


# ignore warnings
import warnings
warnings.filterwarnings("ignore")


if __name__ == '__main__':

    # Parameters
    params = {'batch_size': 4,
          'shuffle': True,
          'num_workers': 0,
          'drop_last': True,
          'pin_memory': False}
    epochs = 200


    # Datasets and Dataloaders
    seq_data = RNADataset()

    tot_length = len(seq_data)
    test_length = int(0.2*tot_length)
    train_set, test_set = random_split(seq_data, [tot_length - test_length, test_length], generator=torch.Generator().manual_seed(42))
    print(f'Training set: {len(train_set)} elements')
    print(f'Test set: {len(test_set)} elements')

    train_dataloader = DataLoader(train_set,**params)
    test_dataloader = DataLoader(test_set,**params)


    # Model and optimizer
    top_pars = pickle.load(open('../data/CSV_minimized/top_par.p', 'rb'))
    dat_pars = pickle.load(open('../data/CSV_minimized/dat_par.p', 'rb'))
    model = Model(top_pars=top_pars, dat_pars=dat_pars)
    torch.save(model.state_dict(), '../results/initial_values.pth')
    model.load_state_dict(torch.load("../results/initial_values.pth"))

    # sd = model.state_dict()
    # sd['couplings.bonds'][:] = 10
    # sd['couplings.angles'][:] = 10
    # sd['couplings.torsions'][:] = 10
    # model.load_state_dict(sd)
    # torch.save(model.state_dict(), '../results/all_10.pth')

    weight_matrix = list(model.parameters())[:4]
    optimizer = torch.optim.Adam(weight_matrix, lr=1e-3)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', factor=0.5, patience=50, cooldown=1000, threshold=1e-7, verbose=True)
    my_loss = loss_fn  # _with_grad
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
    loss_function = np.array([train_loss, test_loss])
    
    
    # Save state_dict (parameters of the model) and loss values
    np.save('../results/b4_iv_e3/loss_200ep.npy', loss_function)
    torch.save(model.state_dict(), '../results/b4_iv_e3/200ep.pth')

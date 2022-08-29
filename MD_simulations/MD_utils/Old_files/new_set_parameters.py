### Writes new parameters in scale_RNA.dat and parametres.top files
### This version uses parameters from data/NewResults and LocalEnergy class

import numpy as np
import torch
import torch.nn as nn
from torch.nn import Parameter
import shutil
import sys

# Ignore warnings
import warnings
warnings.filterwarnings("ignore")


class Interactions(nn.Module):
    def __init__(self):
        super(Interactions, self).__init__()
        self.bonds = 0
        self.angles = 0
        self.torsions = 0


class LocalEnergy_dummy(nn.Module):

    def __init__(self, top_pars=None, dat_pars=None, device='cpu'):
        super(LocalEnergy_dummy, self).__init__()
        if top_pars is None:
            top_pars = {'bond_type': torch.zeros(15, 2), 'angle_type': torch.zeros(13, 2), 'torsion_type': torch.zeros(25, 3)}
            top_pars['torsion_type'][:, 1] = torch.tensor([*([1] * 19), *[3, 1, 1, 2, 3, 2]])
        if dat_pars is None:
            dat_pars = torch.zeros(47)
        self.device = device
        self.couplings = Interactions()
        self.equil = Interactions()
        self.opt_pars = Parameter(torch.tensor(dat_pars, dtype=torch.float, device=self.device, requires_grad=True))
        self.couplings.bonds = Parameter(torch.tensor(top_pars['bond_type'][:, 0], dtype=torch.float, device=self.device, requires_grad=True))
        self.equil.bonds = Parameter(torch.tensor(top_pars['bond_type'][:, 1], dtype=torch.float, device=self.device, requires_grad=True))
        self.couplings.angles = Parameter(torch.tensor(top_pars['angle_type'][:, 0], dtype=torch.float, device=self.device, requires_grad=True))
        self.equil.angles = Parameter(torch.tensor(top_pars['angle_type'][:, 1], dtype=torch.float, device=self.device, requires_grad=True))
        self.couplings.torsions = Parameter(torch.tensor(top_pars['torsion_type'][:, 0], dtype=torch.float, device=self.device, requires_grad=True))
        self.equil.torsions = Parameter(torch.tensor(top_pars['torsion_type'][:, 2], dtype=torch.float, device=self.device, requires_grad=True))
        self.multiplicity = torch.tensor(top_pars['torsion_type'][:, 1], dtype=torch.int64, device=self.device, requires_grad=False)




def make_dat(par_filepath,overwrite=True):   # model parameters saved as .pth file

    model = LocalEnergy_dummy()
    model.load_state_dict((torch.load(par_filepath)))
    # varied_par_indexes = [*range(0,10), *range(28,38), *[46]]
    dat_pars = model.opt_pars.data
    # new_values = np.zeros(47)
    # new_values[varied_par_indexes] = dat_pars[varied_par_indexes]

    with open('scale_RNA.dat','r') as f:
        data = np.asarray([line.strip().split(maxsplit=2) for line in f])
    data[:,1] = [f'{v:.3f}' for v in dat_pars]  # new_values

    if overwrite:
        out_file = 'scale_RNA.dat'
        shutil.copy('scale_RNA.dat', 'scale_RNA_old.dat')
    else:
        out_file = 'scale_RNA_new.dat'

    col_format = "{:>4}" + "{:>9}" + "         " + "{:<60}" + "\n"
    with open(out_file, 'w') as of:
        for row in data:
            of.write(col_format.format(*row))

    return


def parameters_line(list_of_lines):
    line_number = 0
    for idx, line in enumerate(list_of_lines):
        if 'i' in line:
            line_number = idx
            break
    for idx, line in enumerate(list_of_lines[line_number:]):
        if len(line.split()[0]) > 10:
            line_number += idx
            break
    return line_number


def pars_to_text(par_filepath):

    model = LocalEnergy_dummy()
    model.load_state_dict((torch.load(par_filepath)))
    model_parameters = [p.detach().numpy() for p in model.parameters()]
    mul = np.array(model.multiplicity.squeeze())  # .reshape(25,1)
    model_parameters.append(mul)

    # Permute the list to the right order:
    perm_idx = (0,1,4,2,5,3,7,6)
    model_parameters = [model_parameters[i] for i in perm_idx]
    text = []

    for i in range(1,8):
        j_max = len(model_parameters[i])
        str = ""
        for j in range(j_max):
            str += "{:>16}".format(f"{model_parameters[i][j]:.8E}")
            if (j+1) % 5 == 0 or j == j_max-1:
                text.append(str+"\n")
                str = ""

    return text


def make_top(par_filepath,overwrite=True):   # model parameters saved as .pth file

    in_file = 'parametres_RNA.top'
    text_pars = pars_to_text(par_filepath)

    with open(in_file,'r') as fin:
        list_of_lines = fin.readlines()
        starting_line = parameters_line(list_of_lines)
    list_of_lines[starting_line:(starting_line+27)] = text_pars

    if overwrite:
        out_file = 'parametres_RNA.top'
        shutil.copy(in_file, 'parametres_RNA_old.top')
    else:
        out_file = 'parametres_RNA_new.top'

    with open(out_file,'w') as fo:
        for line in list_of_lines:
            fo.write(line)

    return


if __name__ == "__main__":
    par_filepath = sys.argv[1]   # 200_b1_e3e3e4_psgd_sdeq.pth
    make_dat(par_filepath)
    make_top(par_filepath)

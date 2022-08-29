# Script to write the table of parameters that are then provided to the HiRE model in molecular
# dynamics simulations. It reads the dictionary state from PTH file (obtained with the model with less dof) 
# and prints the resulting parameters in parameters_table.txt file 

from __future__ import print_function, division
# import sys
import torch
import pickle
import numpy as np
from cmath import phase, rect
from data_classes import Model


def print_parameters(path_to_file):

    model = Model()
    model.load_state_dict(torch.load(path_to_file))
    globalss = model.globals.tolist()
    dat_pars = pickle.load(open('dat_par.p', 'rb'))
    print(dat_pars)
    dat_pars[0] = globalss[0]
    dat_pars[1:9] = [1]*8
    dat_pars[9] = globalss[2]
    dat_pars[28:38] = [1]*10
    dat_pars[-1] = globalss[1]
    bonds = np.array([np.abs(model.couplings.bonds.tolist()), model.equil.bonds.tolist()]).transpose()
    # angles_values = [v*180/torch.pi for v in model.equil.angles.tolist()]
    angles = np.array([np.abs(model.couplings.angles.tolist()), model.equil.angles.tolist()]).transpose()
    phases = [phase(rect(1,d)) for d in model.equil.torsions] # *180/torch.pi
    torsions = np.array([np.abs(model.couplings.torsions.tolist()), model.multiplicity.tolist(), phases]).transpose()

    with open('parameters_table.txt', 'w') as of:

        of.write('GLOBAL PARAMETERS \n')
        for i,x in enumerate(dat_pars):
            of.write('{:>9}'.format(f'{x:.4f}'))
            if (i+1) % 5 == 0:
                of.write('\n')
        of.write('\n')

        of.write('\nBONDS PARAMETERS \n')
        for row in bonds:
            of.write(('{:>9}'* 2).format(f'{row[0]:.4f}', f'{row[1]:.4f}') + '\n')

        of.write('\nANGLES PARAMETERS \n')
        for row in angles:
            of.write(('{:>9}' * 2).format(f'{row[0]:.4f}', f'{row[1]:.4f}') + '\n')

        of.write('\nTORSIONS PARAMETERS \n')
        for row in torsions:
            of.write(('{:>9}' * 3).format(f'{row[0]:.4f}', f'{row[1]:.2f}', f'{row[2]:.4f}') + '\n')

    print("File parameters_table.txt created ")
    return


if __name__ == '__main__':
    filename = '300ep.pth'   # sys.argv[1]   
    path_to_file = filename
    print_parameters(path_to_file)

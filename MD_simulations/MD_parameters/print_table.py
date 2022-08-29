from __future__ import print_function, division
import sys
import torch
import numpy as np
from cmath import phase, rect
from data_classes import LocalEnergy


def print_parameters(path_to_file):

    model = LocalEnergy()
    model.load_state_dict(torch.load(path_to_file))
    dat_pars = model.opt_pars.tolist()
    bonds = np.array([model.couplings.bonds.tolist(), model.equil.bonds.tolist()]).transpose()
    angles = np.array([model.couplings.angles.tolist(), model.equil.angles.tolist()]).transpose()
    phases = [phase(rect(1,d)) for d in model.equil.torsions]
    torsions = np.array([model.couplings.torsions.tolist(), model.multiplicity.tolist(), phases]).transpose()

    with open('parameters_table.txt', 'w') as of:
	
        of.write(path_to_file+'\n')
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
    filename = sys.argv[1]   # 200_b1_e3e3e4_psgd_sdeq.pth
    path_to_file = filename
    print_parameters(path_to_file)

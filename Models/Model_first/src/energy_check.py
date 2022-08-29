# Check if the implementation of the model with the original parameters goives the same results as HiRE
# Deprecated: old versions of functions are used 

import numpy as np
import sys
import torch
import LocalEnergy as le
import load_data as ld

if len(sys.argv) < 2:
    print('No sequence checked')
    exit()

filename = sys.argv[1]
# rootpath = '/home/flechenault/Documents/Gianluca/Structure_DB_Amber_HiRE'

types, pars = ld.load_pars()
bond_type = types['bond_type']
angle_type = types['angle_type']
torsion_type = types['torsion_type']

# dataset = pd.read_pickle('/home/flechenault/Documents/Gianluca/dataset.p')
dataset = ld.create_dataset()
for i in range(len(dataset)):
    ds = dataset[i]
    if ds['name'] == filename:
        atoms = ds['atoms']
        coords = torch.tensor(np.array(ds['coords']))
        bonds = ds['bonds']
        angles = ds['angles']
        tors = ds['torsions']

        E_bonds = le.bonds_energy(coords, bonds, bond_type, pars).item()
        E_angles = le.angles_energy(atoms, coords, angles, angle_type, pars).item()
        E_tors = le.torsions_energy(atoms, coords, tors, torsion_type, pars).item()
        print('Energies for sequence '+ds['name'])
        print(E_bonds, E_angles, E_tors)
        break

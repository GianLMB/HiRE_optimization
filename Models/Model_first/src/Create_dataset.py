# OLD VERSION -> DEPRECATED
# Create dataset, saved as dictionary, from files generated during coarse graining

import pandas as pd
import numpy as np
from Bio.PDB import *
import os


def extract_PDB(filename):
    parser = PDBParser()
    model = parser.get_structure('Y', filename)[0]
    coords = []
    for a in model.get_atoms():
        coords.append(a.get_coord())
    coords = np.array(coords)
    return coords


def extract_top(filename):
    with open(filename, 'r') as f:
        reader = f.read()

        text = reader.split("SECTION PARTICLE_MASSES")[1].split("SECTION PARTICLE_TYPE")[0].strip()
        mass = [float(i) for i in text.split()]
        # print(mass)

        text = reader.split("SECTION PARTICLE_TYPE")[1].split("SECTION CHARGES")[0].strip()
        atom_type = [int(i) for i in text.split()]
        # print(atom_type)

        text = reader.split("SECTION CHARGES")[1].split("SECTION BOND_FORCE_CONSTANT")[0].strip()
        charge = [float(i) for i in text.split()]
        # print(charge)

        text = reader.split("SECTION BONDS")[1].split("SECTION ANGLES")[0].strip()
        bonds = [int(i) for i in text.split()]
        # print(bonds)

        text = reader.split("SECTION ANGLES")[1].split("SECTION DIHEDRALS")[0].strip()
        angles = [int(i) for i in text.split()]
        # print(angles)

        text = reader.split("SECTION DIHEDRALS")[1].strip()
        tors = [int(i) for i in text.split()]

        return mass, atom_type,charge,bonds,angles,tors


def is_float(n):
    try:
        float(n)
        return True
    except:
        return False

def extract_energy(filename):
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            energy = [float(n) for n in f.read().split() if is_float(n)]
    else:
        energy = 0
    return energy


# Loop over all files in the Input directory

filelist = []
path = '/home/flechenault/Documents/Gianluca/Structure_DB_Amber_HiRE/Prep_pureHire/Input'
for f in os.listdir(path):
    if f.endswith(".pdb"):
        filelist.append(f.split('.')[0])

rootpath = '/home/flechenault/Documents/Gianluca/Structure_DB_Amber_HiRE'

for file in filelist:
    coords = extract_PDB(rootpath+'/Prep_pureHire/Output/'+file+'/'+file+'_CG.pdb')
    mass, atom_type, charge, bonds, angles, tors = extract_top(rootpath+'/Prep_pureHire/Output/'+file+'/parameters.top')
    energy = extract_energy(rootpath+'/FA_PDB/Amber_energies/'+file+'_energy')
    data = {'atom_type': atom_type,
            'mass': mass,
            'charge': charge,
            'coords': coords,
           # 'x': coords[:, 0],
           # 'y': coords[:, 1],
           # 'z': coords[:, 2],
            'bonds': bonds,
            'angles': angles,
            'torsions': tors,
            'energy': energy}
    ds = pd.Series(data)
    ds.to_pickle('/home/flechenault/Documents/Gianluca/pkl_database/'+file+'.pkl')

# ds2 = pd.read_pickle('/home/flechenault/Documents/Gianluca/pkl_database/2g1w.pkl')
# print(ds2[0])



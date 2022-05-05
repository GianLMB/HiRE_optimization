import numpy as np
import Bio.PDB as bpdb
import os
import re
import pickle
import torch


def extract_PDB(filename):
    parser = bpdb.PDBParser()
    model = parser.get_structure('Y', filename)[0]
    coords = []
    for a in model.get_atoms():
        coords.append(a.get_coord())
    coords = np.array(coords, dtype='float64')
    return coords


def extract_top(filename):
    with open(filename, 'r') as f:
        reader = f.read()

        text = reader.split("SECTION RESIDUE_LABEL")[1].split("SECTION RESIDUE_POINTER")[0].strip()
        res_label = re.split(r'[ i]+', text)

        text = reader.split("SECTION RESIDUE_POINTER")[1].split("SECTION CHAIN_POINTER")[0].strip()
        res_pointer = np.array([int(i) for i in text.split()]).reshape(-1, 2).transpose()

        text = reader.split("SECTION PARTICLE_MASSES")[1].split("SECTION PARTICLE_TYPE")[0].strip()
        mass = [float(i) for i in text.split()]
        # print(mass)

        text = reader.split("SECTION PARTICLE_TYPE")[1].split("SECTION CHARGES")[0].strip()
        atoms = [int(i) for i in text.split()]
        # print(atoms)

        text = reader.split("SECTION CHARGES")[1].split("SECTION BOND_FORCE_CONSTANT")[0].strip()
        charge = [float(i) for i in text.split()]
        # print(charge)

        text = reader.split("SECTION BONDS")[1].split("SECTION ANGLES")[0].strip()
        bonds = np.array([int(i) for i in text.split()]).reshape(-1, 3).transpose()
        # print(bonds)

        text = reader.split("SECTION ANGLES")[1].split("SECTION DIHEDRALS")[0].strip()
        angles = np.array([int(i) for i in text.split()]).reshape(-1, 4).transpose()
        # print(angles)

        text = reader.split("SECTION DIHEDRALS")[1].strip()
        tors = np.array([int(i) for i in text.split()]).reshape(-1, 5).transpose()

        top_data = {
            'atoms': atoms,
            'res_label': res_label,
            'res_pointer': res_pointer,
            'mass': mass,
            'charge': charge,
            'bonds': bonds,
            'angles': angles,
            'torsions': tors
        }

        return top_data


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


def load_pars():
    rootpath = '/home/flechenault/Documents/Gianluca/Structure_DB_Amber_HiRE'

    with open(rootpath + '/Prep_pureHire/Output/2g1w/' + 'parameters.top', 'r') as f:
        reader = f.read()

        text = reader.split("SECTION BOND_FORCE_CONSTANT")[1].split("SECTION ANGLE_FORCE_CONSTANT")[0].strip()
        bond_type = np.array([float(i) for i in text.split() if is_float(i)], dtype='float64').reshape(2, -1).transpose()

        text = reader.split("SECTION ANGLE_FORCE_CONSTANT")[1].split("SECTION DIHEDRAL_FORCE_CONSTANT")[0].strip()
        angle_type = np.array([float(i) for i in text.split() if is_float(i)], dtype='float64').reshape(2, -1).transpose()

        text = reader.split("SECTION DIHEDRAL_FORCE_CONSTANT")[1].split("SECTION BONDS")[0].strip()
        torsion_type = np.array([float(i) for i in text.split() if is_float(i)], dtype='float64').reshape(3, -1).transpose()

    fixed_pars = {
        'bond_type': bond_type,
        'angle_type': angle_type,
        'torsion_type': torsion_type
    }

    with open(rootpath + '/Prep_pureHire/scale_RNA.dat', 'r') as f:
        pars = []
        for line in f.readlines():
            pars.append(float(line.split()[1]))
        pars = np.array(pars, dtype='float64')

    return fixed_pars, pars


def create_batches(dataset,batch_size):
    n_batch = int(np.floor(len(dataset) / batch_size))
    n_seq = n_batch * batch_size
    np.random.seed(42)
    np.random.shuffle(dataset)
    for i, data in enumerate(dataset):
        data['batch_ind'] = int(np.floor(i / batch_size))
    return dataset[0:n_seq]


def create_dataset(batch=False,batch_size=32):

    try:
        dataset = pickle.load(open('/home/flechenault/Documents/Gianluca/RNAEnergy/data/RNAseq/dataset.p', 'rb'))
        print('dataset loaded')

    except:
        filelist = []
        path = '/home/flechenault/Documents/Gianluca/Structure_DB_Amber_HiRE/Prep_pureHire/Cropped_Output'

        for f in os.listdir(path):
            filelist.append(f)
        #    if f.endswith(".pdb"):
        #        filelist.append(f.split('.')[0])

        rootpath = '/home/flechenault/Documents/Gianluca/Structure_DB_Amber_HiRE'
        dataset = []

        for file in filelist:
            print(file)
            coords = extract_PDB(rootpath+'/Prep_pureHire/Cropped_Output/'+file+'/'+file+'_CG.pdb')
            top_data = extract_top(rootpath+'/Prep_pureHire/Cropped_Output/'+file+'/parameters.top')
            energy = extract_energy(rootpath+'/FA_PDB/Cropped_energies/'+file+'_energy')
            RNA_seq = {
                    'name': file,
                    'atoms': top_data['atoms'],
                    'res_label': top_data['res_label'],
                    'res_pointer': top_data['res_pointer'],
                    'mass': top_data['mass'],
                    'charge': top_data['charge'],
                    'coords': coords,
                    'bonds': top_data['bonds'],
                    'angles': top_data['angles'],
                    'torsions': top_data['torsions'],
                    'energy': energy}
            pickle.dump(RNA_seq,open('/home/flechenault/Documents/Gianluca/RNAEnergy/data/RNAseq/'+file+'.p', 'wb'))
            dataset.append(RNA_seq)

        pickle.dump(dataset, open('/home/flechenault/Documents/Gianluca/RNAEnergy/data/RNAseq/dataset.p', 'wb'))
        print(f'Created dataset of {len(dataset)} structures')
        print('dataset allocated')

    if batch:
        dataset = create_batches(dataset,batch_size)

    return dataset


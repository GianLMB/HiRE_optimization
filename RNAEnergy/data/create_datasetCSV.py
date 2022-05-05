import numpy as np
import Bio.PDB as bpdb
import os
import re
import pandas as pd
import pickle


def extract_PDB(filename):
    parser = bpdb.PDBParser()
    model = parser.get_structure('Y', filename)[0]
    coords = []
    for a in model.get_atoms():
        coords.append(a.get_coord())
    coords = np.array(coords, dtype='float64').reshape(-1)
    return coords


# converts str res_label to int
def res_to_int(res_label):
    res_dict = {
        'GUA': 1,
        'ADE': 2,
        'CYT': 3,
        'URA': 4
    }
    for i in range(len(res_label)):
        res_label[i] = res_dict[res_label[i]]
    return np.array(res_label)


def extract_top(filename):
    with open(filename, 'r') as f:
        reader = f.read()

        text = reader.split("SECTION RESIDUE_LABELS")[1].split("SECTION RESIDUE_POINTER")[0].strip()
        res_label = re.split(r'[ i]+', text)
        res_label = res_to_int(res_label)

        text = reader.split("SECTION RESIDUE_POINTER")[1].split("SECTION CHAIN_POINTER")[0].strip()
        res_pointer = np.array([int(i) for i in text.split()])

        text = reader.split("SECTION PARTICLE_MASSES")[1].split("SECTION PARTICLE_TYPE")[0].strip()
        mass = np.array([float(i) for i in text.split()])
        # print(mass)

        text = reader.split("SECTION PARTICLE_TYPE")[1].split("SECTION CHARGES")[0].strip()
        atoms = np.array([int(i) for i in text.split()])
        # print(atoms)

        text = reader.split("SECTION CHARGES")[1].split("SECTION BOND_FORCE_CONSTANT")[0].strip()
        charge = np.array([float(i) for i in text.split()])
        # print(charge)

        text = reader.split("SECTION BONDS")[1].split("SECTION ANGLES")[0].strip()
        bonds = np.array([int(i) for i in text.split()])
        # print(bonds)

        text = reader.split("SECTION ANGLES")[1].split("SECTION DIHEDRALS")[0].strip()
        angles = np.array([int(i) for i in text.split()])
        # print(angles)

        text = reader.split("SECTION DIHEDRALS")[1].strip()
        tors = np.array([int(i) for i in text.split()])

        top_data = {
            'atoms': atoms,
            'res_label': res_label,
            'res_pointer': res_pointer,
            'mass': mass,
            'charge': charge,
            'bonds': bonds,
            'angles': angles,
            'torsions': tors,
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
            energy = np.array([float(n) for n in f.read().split() if is_float(n)])
    else:
        energy = 0
    return energy


def create_dataset():

    filelist = []
    filepath = '../Structure_DB_Amber_HiRE/Prep_pureHire/Cropped_Output'

    for f in os.listdir(filepath):
        filelist.append(f)

    rootpath = '../Structure_DB_Amber_HiRE'
    dataset = []

    # create list of dictionaries for all valid sequences
    for file in filelist:
        if file.split('_')[0] == '2ho7':  # sequence for which HiRE does not work
            continue
        top_data = extract_top(rootpath+'/Prep_pureHire/Cropped_Output/'+file+'/parameters.top')
        if len(top_data['res_label'])!= 7:  # check that the length is 7 residues
            continue
        energy = extract_energy(rootpath + '/FA_PDB/Cropped_energies/' + file + '_energy')
        # print(file, len(energy), energy)
        if len(energy) == 0:  # check that energy was computed
            print(file + ': no energy was computed')
            continue
        coords = extract_PDB(rootpath+'/Prep_pureHire/Cropped_Output/'+file+'/'+file+'_CG.pdb')
        RNA_seq = {
                'atoms': pd.Series(top_data['atoms'], index=range(len(top_data['atoms']))),
                'res_label': pd.Series(top_data['res_label'], index=range(len(top_data['res_label']))),
                'res_pointer': pd.Series(top_data['res_pointer'], index=range(len(top_data['res_pointer']))),
                'mass': pd.Series(top_data['mass'], index=range(len(top_data['mass']))),
                'charge': pd.Series(top_data['charge'], index=range(len(top_data['charge']))),
                'coords': pd.Series(coords, index=range(len(coords))),
                'bonds': pd.Series(top_data['bonds'], index=range(len(top_data['bonds']))),
                'angles': pd.Series(top_data['angles'], index=range(len(top_data['angles']))),
                'torsions': pd.Series(top_data['torsions'], index=range(len(top_data['torsions']))),
                'energy': pd.Series(energy, index=range(len(energy)))
        }
        dataset.append([file, RNA_seq])

    print(f'dictionary created of {len(dataset)} sequences')
    return dataset


def save_dataset(dataset):

    os.makedirs('data/SeqCSV', exist_ok=True)
    max_dim = max(len(sample[1]['torsions']) for sample in dataset)  # torsions array is always the longer one
    seq_frame = []

    for sample in dataset:

        seq_frame_line = [sample[0]]
        RNA_seq = sample[1]
        df = pd.DataFrame(data=RNA_seq, index=range(max_dim))
        df.to_csv('data/SeqCSV/'+sample[0]+'.csv', index=False)

        for key in RNA_seq:
            seq_frame_line.append(len(RNA_seq[key]))
        seq_frame.append(seq_frame_line)

    seq_frame = np.array(seq_frame)
    df = pd.DataFrame(data=seq_frame)
    df.to_csv('data/SeqCSV/seq_frame.csv',header=False,index=False)
    print('seq_frame created')

    return


def load_pars():
    rootpath = '../Structure_DB_Amber_HiRE'

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

    pickle.dump(fixed_pars, open('data/SeqCSV/fixed_pars.p', 'wb'))
    pickle.dump(pars, open('data/SeqCSV/pars.p', 'wb'))
    print('parameters saved')

    return


if __name__ == '__main__':
    dataset = create_dataset()
    save_dataset(dataset)
    load_pars()

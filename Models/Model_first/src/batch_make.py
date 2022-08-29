# script to parse the original PDB structures in subsequences of 7 nucleotides each,
# which are then minimized to obtain the dataset for the optimization model.

import numpy as np
import Bio.PDB as bpdb
import os

path_in = '../Structure_DB_Amber_HiRE/FA_PDB/Source_PDB/'       # original structures
path_out = '../Structure_DB_Amber_HiRE/FA_PDB/Cropped_PDB/'     # output folder
seed = 42

filelist = []
for file in os.listdir(path_in):
    if file.endswith(".pdb"):
        filelist.append(file)

for file in filelist:

    # get the pdb structure and initialize variables
    s = bpdb.PDBParser().get_structure('tmp', path_in + file)
    io = bpdb.PDBIO()
    io.set_structure(s)

    chain = list(s[0].get_chains())[0].id
    first_res = list(s[0].get_residues())[0].id[1]
    n_res = (len(s[0][chain]))
    n_cuts = int(n_res/3)           
    # new_resnums = range(1, 8)

    # check if n_res >= 7, else just copy the sequence in the new folder
    if n_res < 7:
        newfile = path_out + file
        io.save(newfile)
        continue

    # loop to generate n_res/3 cropped sequences of length 7 residues for each file,
    for i in range(n_cuts):

        start_res = first_res + np.random.randint(0, n_res-6)
        end_res = start_res + 6

        # classical way to handle parsing in Bio.PDB: overwrite ResSelect class
        class ResSelect(bpdb.Select):  

            # select only the residues between start_res and end_res
            def accept_residue(self, res):
                if res.id[1] < start_res or res.id[1] > end_res:
                    return 0
                else:
                    return 1

            # eliminate P, P1, P2 from beginning of the new sequences -> necessary for Amber
            def accept_atom(self, atom):
                if atom.get_parent().id[1] == start_res and (atom.get_name() == 'P' or atom.get_name() == 'OP1' or atom.get_name() == 'OP2'):
                    return 0
                else:
                    return 1

        newfile = path_out + file.split('.')[0] + '_c' + str(i+1) + '.pdb'
        io.save(newfile, ResSelect())

    print(f'file {file} cropped {n_cuts} times')

#        replace residues id with 1:end -> deprecated, already done by Amber when performing minimization
#        s_new = bpdb.PDBParser().get_structure('tmp2', newfile)
#
#        for model in s_new:
#            for chain in model:
#                for j, residue in enumerate(chain.get_residues()):
#                    res_id = list(residue.id)
#                    res_id[1] = new_resnums[j]
#                    residue.id = tuple(res_id)
#
#        io2 = bpdb.PDBIO()
#        io2.set_structure(s_new)
#        io2.save(newfile)

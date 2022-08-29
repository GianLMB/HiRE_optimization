import numpy as np
import torch
import torch.nn as nn
from torch.nn import Parameter
import shutil

def make_dat(par_filepath,overwrite=True):   

    # varied_par_indexes = [*range(0,10), *range(28,38), *[46]]
    # new_values = np.zeros(47)
    # new_values[varied_par_indexes] = dat_pars[varied_par_indexes]
    with open(par_filepath,'r') as f:
        text = f.read().split("GLOBAL PARAMETERS")[1].split("BONDS PARAMETERS")[0].strip()
        dat_pars = [float(i) for i in text.split()]

    with open('scale_RNA.dat','r') as f:
        data = np.asarray([line.strip().split(maxsplit=2) for line in f])
    data[:,1] = [f'{v:.3f}' for v in dat_pars]

    if overwrite:
        shutil.copy('scale_RNA.dat', 'scale_RNA_old.dat')
        out_file = 'scale_RNA.dat'
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

    with open(par_filepath,'r') as f:
        reader = f.read()

        text = reader.split("BONDS PARAMETERS")[1].split("ANGLES PARAMETERS")[0].strip()
        bonds = [float(i) for i in text.split()]
        bonds_c, bonds_e = bonds[::2], bonds[1::2]

        text = reader.split("ANGLES PARAMETERS")[1].split("TORSIONS PARAMETERS")[0].strip()
        angles = [float(i) for i in text.split()]
        angles_c, angles_e = angles[::2], angles[1::2]

        text = reader.split("TORSIONS PARAMETERS")[1].strip()
        torsions = [float(i) for i in text.split()]
        torsions_c, torsions_m, torsions_e = torsions[::3],torsions[1::3], torsions[2::3]

    model_parameters = [bonds_c,bonds_e,angles_c,angles_e,torsions_c,torsions_m,torsions_e]
    text = []

    for i in range(len(model_parameters)):
        j_max = len(model_parameters[i])
        str = ""
        for j in range(j_max):
            str += "{:>16}".format(f"{model_parameters[i][j]:.8E}")
            if (j+1) % 5 == 0 or j == j_max-1:
                text.append(str+"\n")
                str = ""

    return text


def make_top(par_filepath,overwrite=True):   

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
    par_filepath = "parameters_table.txt"
    make_dat(par_filepath)
    make_top(par_filepath)


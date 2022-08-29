# All relaxation files contained in /Relaxation directory, and saved as relaxation_ + seq_name.pdb

import mdtraj as md
import numpy as np
import sys
import matplotlib.pyplot as plt
# import time
from matplotlib.ticker import NullFormatter
import matplotlib.backends.backend_pdf


# computes bonds
def compute_bonds_lengths(P,O,C,R4,R1,C1,U1,A1,G1,bdict,last_frame):

    bdict['P-O'].extend(md.compute_distances(last_frame,np.array([P,O[1:]]).transpose()).reshape(-1,).tolist())
    bdict['O-C'].extend(md.compute_distances(last_frame,np.array([O,C]).transpose()).reshape(-1,).tolist())
    bdict['C-R4'].extend(md.compute_distances(last_frame,np.array([C,R4]).transpose()).reshape(-1,).tolist())
    bdict['R4-P'].extend(md.compute_distances(last_frame,np.array([R4[:-1],P]).transpose()).reshape(-1,).tolist())
    bdict['R4-R1'].extend(md.compute_distances(last_frame,np.array([R4,R1]).transpose()).reshape(-1,).tolist())

    if C1.size > 0:
        bdict['R1-C1'].extend(md.compute_distances(last_frame,np.array([C1-1,C1]).transpose()).reshape(-1,).tolist())
    if U1.size > 0:
        bdict['R1-U1'].extend(md.compute_distances(last_frame,np.array([U1-1,U1]).transpose()).reshape(-1,).tolist())
    if A1.size > 0:
        bdict['R1-A1'].extend(md.compute_distances(last_frame,np.array([A1-1,A1]).transpose()).reshape(-1,).tolist())
        bdict['A1-A2'].extend(md.compute_distances(last_frame,np.array([A1,A1+1]).transpose()).reshape(-1,).tolist())
    if G1.size > 0:
        bdict['R1-G1'].extend(md.compute_distances(last_frame,np.array([G1-1,G1]).transpose()).reshape(-1,).tolist())
        bdict['G1-G2'].extend(md.compute_distances(last_frame,np.array([G1,G1+1]).transpose()).reshape(-1,).tolist())

    return


# computes angles
def compute_bonds_angles(P,O,C,R4,R1,C1,U1,A1,G1,adict,last_frame):

    adict['P-O-C'].extend(md.compute_angles(last_frame,np.array([P,O[1:],C[1:]]).transpose()).reshape(-1,).tolist())
    adict['O-C-R4'].extend(md.compute_angles(last_frame,np.array([O,C,R4]).transpose()).reshape(-1,).tolist())
    adict['C-R4-P'].extend(md.compute_angles(last_frame,np.array([C[:-1],R4[:-1],P]).transpose()).reshape(-1,).tolist())
    adict['R4-P-O'].extend(md.compute_angles(last_frame,np.array([R4[:-1],P,O[1:]]).transpose()).reshape(-1,).tolist())
    adict['C-R4-R1'].extend(md.compute_angles(last_frame,np.array([C,R4,R1]).transpose()).reshape(-1,).tolist())
    adict['R1-R4-P'].extend(md.compute_angles(last_frame,np.array([R1[:-1],R4[:-1],P]).transpose()).reshape(-1,).tolist())

    if A1.size > 0:
        adict['R4-R1-A1'].extend(md.compute_angles(last_frame,np.array([A1-2,A1-1,A1]).transpose()).reshape(-1,).tolist())
        adict['R1-A1-A2'].extend(md.compute_angles(last_frame,np.array([A1-1,A1,A1+1]).transpose()).reshape(-1,).tolist())
    if U1.size > 0:
        adict['R4-R1-U1'].extend(md.compute_angles(last_frame,np.array([U1-2,U1-1,U1]).transpose()).reshape(-1,).tolist())
    if G1.size > 0:
        adict['R4-R1-G1'].extend(md.compute_angles(last_frame,np.array([G1-2,G1-1,G1]).transpose()).reshape(-1,).tolist())
        adict['R1-G1-G2'].extend(md.compute_angles(last_frame,np.array([G1-1,G1,G1+1]).transpose()).reshape(-1,).tolist())
    if C1.size > 0:
        adict['R4-R1-C1'].extend(md.compute_angles(last_frame,np.array([C1-2,C1-1,C1]).transpose()).reshape(-1,).tolist())

    return


# computes dihedrals
def compute_dihedrals_angles(P,O,C,R4,R1,C1,U1,A1,G1,ddict,last_frame,top):

    ddict['C-R4-P-O'].extend(md.compute_dihedrals(last_frame,np.array([C[:-1],R4[:-1],P,O[1:]]).transpose()).reshape(-1,).tolist())
    ddict['R1-R4-P-O'].extend(md.compute_dihedrals(last_frame,np.array([R1[:-1],R4[:-1],P,O[1:]]).transpose()).reshape(-1,).tolist())
    ddict['O-C-R4-P'].extend(md.compute_dihedrals(last_frame,np.array([O[:-1],C[:-1],R4[:-1],P]).transpose()).reshape(-1,).tolist())
    ddict['O-C-R4-R1'].extend(md.compute_dihedrals(last_frame,np.array([O,C,R4,R1]).transpose()).reshape(-1,).tolist())
    ddict['P-O-C-R4'].extend(md.compute_dihedrals(last_frame,np.array([P,O[1:],C[1:],R4[1:]]).transpose()).reshape(-1,).tolist())
    ddict['R4-P-O-C'].extend(md.compute_dihedrals(last_frame,np.array([R4[:-1],P,O[1:],C[1:]]).transpose()).reshape(-1,).tolist())

    if A1.size > 0:
        ddict['R4-R1-A1-A2'].extend(md.compute_dihedrals(last_frame,np.array([A1-2,A1-1,A1,A1+1]).transpose()).reshape(-1,).tolist())
        ddict['R4-A1-A2-R1'].extend(md.compute_dihedrals(last_frame,np.array([A1-2,A1,A1+1,A1-1]).transpose()).reshape(-1,).tolist())
        ddict['C-R4-R1-A1'].extend(md.compute_dihedrals(last_frame,np.array([A1-3,A1-2,A1-1,A1]).transpose()).reshape(-1,).tolist())
        if A1[-1]+2 > top.n_atoms:
            ddict['P-R4-R1-A1'].extend(md.compute_dihedrals(last_frame,np.array([A1+2,A1-2,A1-1,A1]).transpose()).reshape(-1,).tolist())  # check it is not the last residue
        else:
            if A1.size > 1:
                ddict['P-R4-R1-A1'].extend(md.compute_dihedrals(last_frame,np.array([A1[:-1]+2,A1[:-1]-2,A1[:-1]-1,A1[:-1]]).transpose()).reshape(-1,).tolist())
    if U1.size > 0:
        ddict['C-R4-R1-U1'].extend(md.compute_dihedrals(last_frame,np.array([U1-3,U1-2,U1-1,U1]).transpose()).reshape(-1,).tolist())
        if U1[-1]+1 > top.n_atoms:
            ddict['P-R4-R1-U1'].extend(md.compute_dihedrals(last_frame,np.array([U1+2,U1-2,U1-1,U1]).transpose()).reshape(-1,).tolist())  # check it is not the last residue
        else:
            if U1.size > 1:
                ddict['P-R4-R1-U1'].extend(md.compute_dihedrals(last_frame,np.array([U1[:-1]+2,U1[:-1]-2,U1[:-1]-1,U1[:-1]]).transpose()).reshape(-1,).tolist())
    if G1.size > 0:
        ddict['R4-R1-G1-G2'].extend(md.compute_dihedrals(last_frame,np.array([G1-2,G1-1,G1,G1+1]).transpose()).reshape(-1,).tolist())
        ddict['R4-G1-G2-R1'].extend(md.compute_dihedrals(last_frame,np.array([G1-2,G1,G1+1,G1-1]).transpose()).reshape(-1,).tolist())
        ddict['C-R4-R1-G1'].extend(md.compute_dihedrals(last_frame,np.array([G1-3,G1-2,G1-1,G1]).transpose()).reshape(-1,).tolist())
        if G1[-1]+2 > top.n_atoms:
            ddict['P-R4-R1-G1'].extend(md.compute_dihedrals(last_frame,np.array([G1+2,G1-2,G1-1,G1]).transpose()).reshape(-1,).tolist())  # check it is not the last residue
        else:
            if G1.size > 1:
                ddict['P-R4-R1-G1'].extend(md.compute_dihedrals(last_frame,np.array([G1[:-1]+2,G1[:-1]-2,G1[:-1]-1,G1[:-1]]).transpose()).reshape(-1,).tolist())
    if C1.size > 0:
        ddict['C-R4-R1-C1'].extend(md.compute_dihedrals(last_frame,np.array([C1-3,C1-2,C1-1,C1]).transpose()).reshape(-1,).tolist())
        if C1[-1]+1 > top.n_atoms:
            ddict['P-R4-R1-C1'].extend(md.compute_dihedrals(last_frame,np.array([C1+2,C1-2,C1-1,C1]).transpose()).reshape(-1,).tolist())  # check it is not the last residue
        else:
            if C1.size > 1:
                ddict['P-R4-R1-C1'].extend(md.compute_dihedrals(last_frame,np.array([C1[:-1]+2,C1[:-1]-2,C1[:-1]-1,C1[:-1]]).transpose()).reshape(-1,).tolist())

    return


# print results in pdf
def print_pdf(bdict,adict,ddict):

    # distances
    print('\n')
    for key, list in bdict.items():
        print(key+' pair: %i samples' %len(list))
        plt.figure()
        plt.hist(list,bins=int(np.ceil(np.log2(len(list))) + 1))    # Sturge's rule
        plt.title("Bond "+key+" length distribution")

    pdf = matplotlib.backends.backend_pdf.PdfPages("bonds.pdf")
    for fig in range(1, plt.gcf().number + 1): ## will open an empty extra figure :(
        pdf.savefig(fig)
    pdf.close()
    plt.close('all')

    print('\n')
    for key, list in adict.items():
        print(key+' triplet: %i samples' %len(list))
        plt.figure()
        plt.hist(list,bins=int(np.ceil(np.log2(len(list))) + 1))
        plt.title("Angle "+key+" distribution")

    pdf = matplotlib.backends.backend_pdf.PdfPages("angles.pdf")
    for fig in range(1, plt.gcf().number + 1):
        pdf.savefig(fig)
    pdf.close()
    plt.close('all')

    print('\n')
    for key, list in ddict.items():
        print(key+' quadruplet: %i samples' %len(list))
        plt.figure()
        plt.hist(list,bins=int(np.ceil(np.log2(len(list))) + 1))
        plt.title("Dihedral "+key+" distribution")

    pdf = matplotlib.backends.backend_pdf.PdfPages("dihedrals.pdf")
    for fig in range(1, plt.gcf().number + 1):
        pdf.savefig(fig)
    pdf.close()
    plt.close('all')

    return


def main():

    # Dictionaries definitions
    bond_pairs = ['R4-P','R4-R1','R1-G1','R1-A1','R1-U1','R1-C1','G1-G2','A1-A2','C-R4','P-O','O-C']
    bdict = {key: [] for key in bond_pairs}

    angle_triplets = ['R4-R1-A1','R4-R1-U1','R4-R1-G1','R4-R1-C1','R1-A1-A2','R1-G1-G2','P-O-C',
    'O-C-R4','C-R4-P','R4-P-O','C-R4-R1','R1-R4-P']
    adict = {key: [] for key in angle_triplets}

    dihedrals_quadruplets = ['R4-R1-G1-G2','R4-R1-A1-A2','R4-A1-A2-R1','R4-G1-G2-R1','C-R4-R1-A1',
    'C-R4-R1-G1','C-R4-R1-C1','C-R4-R1-U1','P-R4-R1-A1','P-R4-R1-G1','P-R4-R1-C1','P-R4-R1-U1',
    'C-R4-P-O','R1-R4-P-O','O-C-R4-P','O-C-R4-R1','P-O-C-R4','R4-P-O-C']
    ddict = {key: [] for key in dihedrals_quadruplets}

    # Execution
    with open('relaxation_names.txt','r') as f:
        file_list = f.read().splitlines()

    for i,file in enumerate(file_list):

        print(file+",    %i/%i" %(i+1,len(file_list)))
        last_frame = md.load_frame('NewRelaxation/'+file+'.pdb',-1)
        top = last_frame.topology

        P = top.select("name == P")
        if len(P) < 6:
            continue

        O = top.select("symbol == O")  # all sequences start from oxygen, doesn't read name == O5*
        C = O+1                        # doesn't read name == C5* nor symbol == C
        R4 = top.select("name == CA")
        R1 = top.select("name == CY")
        C1 = top.select("name == C1")
        U1 = top.select("name == U1")
        A1 = top.select("name == A1")
        G1 = top.select("name == G1")

        compute_bonds_lengths(P,O,C,R4,R1,C1,U1,A1,G1,bdict,last_frame)
        compute_bonds_angles(P,O,C,R4,R1,C1,U1,A1,G1,adict,last_frame)
        compute_dihedrals_angles(P,O,C,R4,R1,C1,U1,A1,G1,ddict,last_frame,top)

    for key, list in bdict.items():
        bdict[key] = [el*10 for el in list]
    print_pdf(bdict,adict,ddict)

    return


if __name__ == '__main__':
    main()

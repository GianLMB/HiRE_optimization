### Computes energy for local interactions ###

from __future__ import print_function, division
import numpy as np
import torch
import math
import warnings
from numba import cuda
warnings.filterwarnings("ignore")


# TODO: optimize vectorial form of equations to avoid for loops
#### Bonds energy ####  checked with HiRE: correct

def bonds_energy(coords, bonds, bond_type, pars):
    """
    :param coords: 3*Natoms 2-array with xyz coordinates
    :param bonds: 3-array with 3*Nbonds elements: atom1, atom2, bond_type index
    :param bond_type: 2-array with bond_strength, bond_eq for Nbonds
    :param pars: array of parameters to optimize (only first one eneters calculation)
    :return: energy associated with covalent bonds
    """
    a1 = bonds[:, 0] / 3  # gives index of the 1st atom for the coordinates
    a2 = bonds[:, 1] / 3  # gives index of the 2nd atom for the coordinates
    bt = bonds[:, 2] - 1    # gives index of bond type
    n_bonds = bonds.shape[0]
    e_bonds = torch.tensor(0.0)

    for n in range(n_bonds):
        i = int(a1[n])
        j = int(a2[n])
        ind = int(bt[n])
        d = coords[i, :] - coords[j, :]
        d = (torch.dot(d, d))**0.5   # distance between atoms
        e_bonds += bond_type[ind, 0] * (d-bond_type[ind, 1])**2
    e_bonds *= pars[0]
    return e_bonds



#### Angles energy ####  checked with HiRE: correct

def assign_thetatype(atoms, i1, i2, i3):
    a1 = atoms[i1]
    a2 = atoms[i2]
    a3 = atoms[i3]
    k_type = 0
    if a1 == 1:
        if a3 == 3:
            k_type = 4
        elif a3 == 5:
            k_type = 6
    elif a1 == 2:
        if a2 == 1:
            k_type = 3
    elif a1 == 3:
        if a2 == 2:
            k_type = 2
    elif a1 == 4:
        if a2 == 5:
            k_type = 0
        elif a2 == 3:
            k_type = 5
    elif a1 == 5:
        if a2 == 6:
            k_type = 1
        elif a2 == 8:
            k_type = 1
        elif a2 == 4:
            k_type = 7
    else:
        print("ERROR: Unknown angle type")
    return k_type


def angles_energy(atoms, coords, angles, angle_type, pars):
    """
    :param atoms: Natoms 1-array with atom name
    :param coords: 3*Natoms 2-array with xyz coordinates
    :param angles: 4-array with 4*Nangles elements: atom1, atom2, atom 3, angle_type index
    :param angle_type: 2-array with angle_strength, angle_eq for Nangles
    :param pars: array of parameters to optimize (only second one eneters calculation)
    :return: energy associated with angles
    """
    a1 = angles[:, 0] / 3  # gives index of the 1st atom for the coordinates
    a2 = angles[:, 1] / 3  # gives index of the 2nd atom for the coordinates
    a3 = angles[:, 2] / 3  # gives index of the 2nd atom for the coordinates
    at = angles[:, 3] - 1  # gives index of angle type
    n_angles = angles.shape[0]
    e_angles = torch.tensor(0.0)

    for n in range(n_angles):
        i = int(a1[n])
        j = int(a2[n])
        k = int(a3[n])
        ind = int(at[n])
        rij = coords[i, :] - coords[j, :]
        rkj = coords[k, :] - coords[j, :]
        c_an = torch.dot(rij, rkj)/(torch.norm(rij)*torch.norm(rkj))  # cos_theta

        # regularization
        # eps = torch.tensor(0.001)
        # c_an = min(max(c_an, -1+eps), 1-eps)

        an = torch.arccos(c_an)

        # define angle_type multiplier index depending on atoms involved
        thetatype = assign_thetatype(atoms, i, j, k)
        e_angles += pars[46] * pars[thetatype+1] * angle_type[ind, 0] * (an - angle_type[ind, 1])**2
    return e_angles



#### Torsion energies #### checked with HiRE: correct

def assign_phitype(atoms,i1,i2,i3,i4):
    a1 = atoms[i1]
    a2 = atoms[i2]
    a3 = atoms[i3]
    a4 = atoms[i4]
    ktype = 0
    if a1 == 1:
        if a3 == 5:
            ktype = 2
        elif a3 == 3:
            ktype = 4
    elif a1 == 2:
        if a4 == 3:
            ktype = 6
        elif a4 == 5:
            ktype = 7
    elif a1 == 3:
        if a2 == 4:
            ktype = 3
        elif a2 == 2:
            ktype = 8
    elif a1 == 4:
        if a2 == 5:
            ktype = 0
        elif a2 == 6 or a2 == 8:
            ktype = 1
        elif a2 == 3:
            ktype = 9
    elif a1 == 5:
        if a2 == 4:
            ktype = 5
    else:
        print("ERROR: Unknown dihedral type")
    return ktype


def torsions_energy(atoms, coords, tors, tors_type, pars):
    """
    :param atoms: Natoms 1-array with atom name
    :param coords: 3*Natoms 2-array with xyz coordinates
    :param tors: 5-array with 5*Ntors elements: atom1, atom2, atom 3, atom 4, tors_type index
    :param tors_type: 3-array with tors_strength, tors_eq, tors_phase for Ntors
    :param pars: array of parameters to optimize
    :return: energy associated with torsions
    """
    a1 = tors[:, 0] / 3  # gives index of the 1st atom for the coordinates
    a2 = tors[:, 1] / 3  # gives index of the 2nd atom for the coordinates
    a3 = tors[:, 2] / 3  # gives index of the 2nd atom for the coordinates
    a4 = tors[:, 3] / 3  # gives index of the 2nd atom for the coordinates
    tt = tors[:, 4] - 1  # gives index of tors type
    n_tors = tors.shape[0]
    e_tors = torch.tensor(0.0)

    for n in range(n_tors):
        i = int(a1[n])
        j = int(a2[n])
        k = int(a3[n])
        l = int(a4[n])
        ind = int(tt[n])

        rij = coords[i, :] - coords[j, :]
        rkj = coords[k, :] - coords[j, :]
        rkl = coords[k, :] - coords[l, :]
        nj = torch.cross(rij, rkj)
        nj = nj / torch.norm(nj)
        nk = torch.cross(rkl, rkj)
        nk = nk / torch.norm(nk)
        phi = torch.acos(torch.dot(nj, nk))
        phi = torch.tensor(math.pi) - phi*torch.sign(torch.dot(rkj, torch.cross(nk, nj)))

        phitype = assign_phitype(atoms, i, j, k, l)
        e_tors += pars[9] * pars[28+phitype] * tors_type[ind, 0] * (1 + torch.cos(tors_type[ind, 1]*phi - tors_type[ind, 2]))

    return e_tors


def E_local(atoms, coords, bonds, bond_type, angles, angle_type, tors, tors_type, pars):
    E_bonds = bonds_energy(coords, bonds, bond_type, pars)
    E_angles = angles_energy(atoms, coords, angles, angle_type, pars)
    E_tors = torsions_energy(atoms, coords, tors, tors_type, pars)
    return E_bonds + E_angles + E_tors











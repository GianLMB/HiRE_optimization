### Computes energy for local interactions ###

from __future__ import print_function, division
import numpy as np
import torch
import math
import warnings
from numba import cuda
warnings.filterwarnings("ignore")


### New version: for loop avoided, execution time more than halved

#### Bonds energy ####  checked with HiRE: correct

def bonds_energy(coords, bonds, bond_type, pars):

    a1 = (bonds[:,0]/3).long()
    a2 = (bonds[:,1]/3).long()
    idx = bonds[:,2]-1

    d = coords[a1,:]-coords[a2,:]
    d = torch.linalg.norm(d, dim=1)
    energy = (pars[0]*bond_type[idx,0]*(d-bond_type[idx,1])**2).sum()

    return energy



#### Angles energy ####  checked with HiRE: correct

def assign_thetatype(atoms,at1,at2,at3):

    device = 'cpu' if at1.get_device() == -1 else 'cuda'
    thty = torch.zeros(len(at1),device=device).long()
    a1 = atoms[at1]
    a2 = atoms[at2]
    a3 = atoms[at3]

    thty[(a1==5) & ((a2==6) | (a2==8))] = 1
    thty[(a1==3) & (a2==2)] = 2
    thty[(a1==2) & (a2==1)] = 3
    thty[(a1==1) & (a3==3)] = 4
    thty[(a1==4) & (a2==3)] = 5
    thty[(a1==1) & (a3==5)] = 6
    thty[(a1==5) & (a2==4)] = 7

    return thty


def angles_energy(atoms, coords, angles, angle_type, pars):

    a1 = (angles[:,0]/3).long()
    a2 = (angles[:,1]/3).long()
    a3 = (angles[:,2]/3).long()
    idx = angles[:,3]-1

    rij = coords[a1,:] - coords[a2,:]
    rkj = coords[a3, :] - coords[a2, :]
    an = torch.acos((rij*rkj).sum(dim=1) / (torch.linalg.norm(rij, dim=1) * torch.linalg.norm(rkj, dim=1)))
    thetatype = assign_thetatype(atoms,a1,a2,a3)
    energy = (pars[46] * pars[thetatype+1] * angle_type[idx, 0] * (an - angle_type[idx, 1])**2).sum()

    return energy



#### Torsion energies #### checked with HiRE: correct

def assign_phitype(atoms,at1,at2,at3,at4):

    device = 'cpu' if at1.get_device() == -1 else 'cuda'
    phty = torch.zeros(len(at1),device=device).long()
    a1 = atoms[at1]
    a2 = atoms[at2]
    a3 = atoms[at3]
    a4 = atoms[at4]

    phty[(a1==4) & ((a2==6) | (a2==8))] = 1
    phty[(a1==1) & (a3==5)] = 2
    phty[(a1==3) & (a2==4)] = 3
    phty[(a1==1) & (a3==3)] = 4
    phty[(a1==5) & (a2==4)] = 5
    phty[(a1==2) & (a4==3)] = 6
    phty[(a1==2) & (a4==5)] = 7
    phty[(a1==3) & (a2==2)] = 8
    phty[(a1==4) & (a2==3)] = 9

    return phty


def torsions_energy(atoms, coords, tors, tors_type, mult, pars):

    a1 = (tors[:,0]/3).long()
    a2 = (tors[:,1]/3).long()
    a3 = (tors[:,2]/3).long()
    a4 = (tors[:,3]/3).long()
    idx = tors[:,4]-1

    rij = coords[a1,:] - coords[a2,:]
    rkj = coords[a3, :] - coords[a2, :]
    rkl = coords[a3, :] - coords[a4, :]
    nj = torch.linalg.cross(rij, rkj, dim=1)
    nj = nj / torch.linalg.norm(nj, dim=1).unsqueeze(1).expand(-1,3)
    nk = torch.linalg.cross(rkl, rkj, dim=1)
    nk = nk / torch.linalg.norm(nk,dim=1).unsqueeze(1).expand(-1,3)
    phi = torch.acos((nj*nk).sum(dim=1))
    phi = torch.tensor(math.pi) - phi*torch.sign((rkj*torch.linalg.cross(nk, nj, dim=1)).sum(dim=1))
    phitype = assign_phitype(atoms,a1,a2,a3,a4)
    energy = (pars[9] * pars[28+phitype] * tors_type[idx, 0] * (1+torch.cos(mult[idx]*phi - tors_type[idx, 1]))).sum()

    return energy



### Total bonded energy

def E_local(atoms, coords, bonds, bond_type, angles, angle_type, tors, tors_type, pars):
    E_bonds = bonds_energy(coords, bonds, bond_type, pars)
    E_angles = angles_energy(atoms, coords, angles, angle_type, pars)
    E_tors = torsions_energy(atoms, coords, tors, tors_type, pars)
    return E_bonds + E_angles + E_tors

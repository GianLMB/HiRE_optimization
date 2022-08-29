# Computes energy for local interactions: bonds, angles, torsions.
# Parameters are divided in globals, couplings and equilibrium values, and global particle dependent
# parameters are set to 1.
# All functions are written to exploit torch.linalg methods and are checked with respect to the 
# original HiRE model

from __future__ import print_function, division
import torch
import math
import warnings
warnings.filterwarnings("ignore")



### Bonds energy ###

def bonds_energy(coords, bonds, couplings, equil, globalss):

    a1 = (bonds[:,0]/3).long()
    a2 = (bonds[:,1]/3).long()
    idx = bonds[:,2]-1

    d = coords[a1,:]-coords[a2,:]
    d = torch.linalg.norm(d, dim=1)
    energy = (globalss[0] * torch.abs(couplings.bonds[idx])*(d-equil.bonds[idx])**2).sum()

    return energy



### Angles energy ### 

def angles_energy(coords, angles, couplings, equil, globalss):

    a1 = (angles[:,0]/3).long()
    a2 = (angles[:,1]/3).long()
    a3 = (angles[:,2]/3).long()
    idx = angles[:,3]-1

    rij = coords[a1,:] - coords[a2,:]
    rkj = coords[a3, :] - coords[a2, :]
    an = torch.acos((rij*rkj).sum(dim=1) / (torch.linalg.norm(rij, dim=1) * torch.linalg.norm(rkj, dim=1)))
    energy = (globalss[1] * torch.abs(couplings.angles[idx]) * (an - equil.angles[idx])**2).sum()

    return energy



### Torsion energies ###

def torsions_energy(coords, tors, couplings, equil, mult, globalss):

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
    phi = math.pi - phi*torch.sign((rkj*torch.linalg.cross(nk, nj, dim=1)).sum(dim=1))
    energy = (globalss[2] * torch.abs(couplings.torsions[idx]) * (1+torch.cos(mult[idx]*phi - equil.torsions[idx]))).sum()

    return energy

from __future__ import print_function, division
import numpy as np
import torch
import math
import warnings
warnings.filterwarnings("ignore")


# All nonbonded intercations loop over residues


def excluded_volume(res1_idx,res2_idx,res_pointer,coords,pars,Csup,Cinf):  # check Cutoff values Csup, Cinf
    '''
    Computes the excluded volume energy between two residues
    '''
    energy = torch.tensor(0.0).to(device)
    for i in range(res_pointer[res1_idx,0]-1,res_pointer[res1_idx,1]):
        for j in range(res_pointer[res2_idx,0]-1,res_pointer[res2_idx,1]):
            d = coords[i, :] - coords[j, :]
            d = (torch.dot(d, d))**0.5
            # pars[10]: EV steepness, pars[43]: EV barrier, pars[44]: ct2 ratio
            expexcl = torch.exp(pars[10]*(d-pars[44]*ct2))  # TODO: ask what ct2 is
            eexcl = pars[43]*(1 - 1/(1+expexcl))
            if d**2 >= Cinf:
                eexcl *= (Csup+2*d**2-3*Cinf)*(Csup-d**2)**2/(Csup-Cinf)**3
            energy += eexcl
    return energy



def Debye_Huckel(res1_idx,res2_idx,res_pointer,coords,charge,pars):
    '''
    Computes Debye-Huckel energy between two residues
    '''
    i = res_pointer[res1_idx,0] - 1
    j = res_pointer[res2_idx,0] - 1
    d = coords[i, :] - coords[j, :]
    d = (torch.dot(d, d))**0.5
    energy = torch.tensor(pars[11]*4.1507*torch.exp(-d/pars[12])).to(device)
    # TODO: check ionic strength for each structure and write the conversion
    return energy



def stack_params(res1_idx,res2_idx,res_label,pars):
    Itype = res_label[res1_idx]
    Jtype = res_label[res2_idx]
    if Itype < 3 and Jtype < 3:  # pur - pur
        eq, wid, sk = pars[17], pars[20], pars[14]
    elif Itype > 2 and Jtype > 2:  # pyr - pyr
        eq, wid, sk = pars[18], pars[21], pars[15]
    else:  # pur - pyr
        eq, wid, sk = pars[16], pars[19], pars[13]
    return eq, wid, sk


def stacking(res1_idx,res2_idx,res_pointer):
    # write function for stacking
    return



def planarity(res1_idx,res2_idx,res_label,res_pointer,coords,pars):
    '''
    Defines single planarity energy between 2 residues
    '''
    energy = torch.tensor(0.0).to(device)
    Itype = res_label[res1_idx]
    Jtype = res_label[res2_idx]
    d0 = plan_eq_dist(Itype,Jtype)  # function to write, gives eq distance depending on residue type
    I1 = res_pointer[res1_idx,1] - 3
    I2 = res_pointer[res1_idx,1] - 2
    I3 = res_pointer[res1_idx,1] - 1
    nI = torch.cross(coords[I1,:] - coords[I2,:], coords[I3] - coords[I2])
    nI = nI / torch.norm(nI)
    for j in range(res_pointer[res2_idx,1]-3,res_pointer[res2_idx,1]):  # probably can be vectorized
        dIj = torch.abs(torch.dot(nI, coords[I2,:] - coords[j,:])) - d0  # TODO: optimize also on d0 (see definition in NAparams)
        # pars[27]: interaction scale, pars[26]: gaussian width
        energy += pars[27]*torch.exp(-(dIj/pars[26])**2)
    return energy



def hydrogen_bond(res1_idx,res2_idx,res_label,res_pointer,coords,pars):
    # write function for hydrogen bonding and planarity
    
    return



def E_NonBonded(coords,charge,res_pointer,pars):
    energy = torch.tensor(0.0).to(device)
    for I in range(len(res_pointer)):
        for J in range(I+1,len(res_pointer)):
            if I != 0:
                energy += Debye_Huckel(I,J,res_pointer,coords,charge,pars)
            energy += excluded_volume(I,J,res_pointer,coords,pars)


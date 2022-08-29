from Bio.PDB import *
import sys
import matplotlib.pyplot as plt
import numpy as np
import sys
from operator import itemgetter, attrgetter, methodcaller
import sys, os
import string
import pylab
import math
import random 

pdbin = str(sys.argv[1])
ions = int(sys.argv[2])
box = int(sys.argv[3])

def cm(residue):
  p=[0,0,0,0]; wtot=0
  for atom in residue:
    w=1
    wtot+=w
    for i in range(3): p[i]+=atom.get_coord()[i]*w
  for i in range(3): p[i]/=wtot
  p[3]=wtot
  return p

cm

p = PDBParser()
conf = p.get_structure('X', pdbin)

io=PDBIO()
io.set_structure(conf)
io.save("out.pdb")

atnum=0
resnum=0
p=[0,0,0,0]; wtot=0
for model in conf:
    for chain in model:
        for residue in chain:
            resnum=resnum+1
            for atom in residue:
                    w=1
                    wtot+=w
                    atnum=atnum+1
                    for i in range(3): p[i]+=atom.get_coord()[i]*w
for i in range(3): p[i]/=wtot
p[3]=wtot
                    
size=[]
for model in conf:
    for chain in model:
        for residue in chain:
            for atom in residue:
                d=0
                for i in range(3): 
                    d = d + (p[i]-atom.get_coord()[i]*w)*(p[i]-atom.get_coord()[i]*w)
                size.append(math.sqrt(d))

M=max(size)
print max(size)



for i in range(ions):
    dist=0
    while dist < 8:
        random.seed()
        r = float(random.uniform(0,box))
        th = float(random.uniform(0, 3.14))
        phi = float(random.uniform(0, 6.28))
#        print r, th, phi
        x = p[0]+r*math.sin(th)*math.cos(phi)
        y = p[1]+r*math.sin(th)*math.sin(phi)
        z = p[2]+r*math.cos(th)
#        print x,y,z
        distL=[]
        for model in conf:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        for k in range(3):
                            h = np.sqrt((atom.get_coord()[0]-x)*(atom.get_coord()[0]-x) + (atom.get_coord()[1]-y)*(atom.get_coord()[1]-y) + (atom.get_coord()[2]-z)*(atom.get_coord()[2]-z))
                            distL.append(h)
        dist = min(distL)
#        print i, dist
#    print i
    with open(pdbin, "a") as myfile:
        myfile.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format('ATOM',atnum+i+1,'MG', '', ' MG', 'X', resnum+i+1, '', x,y,z,0.00, 0.00))










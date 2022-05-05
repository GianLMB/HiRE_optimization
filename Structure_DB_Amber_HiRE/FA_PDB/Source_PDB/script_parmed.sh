#!/bin/csh

parmed << EOF
parm 1bz2.prmtop 
loadCoordinates 1bz2_min.pdb 
energy 
EOF

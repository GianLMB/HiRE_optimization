#!/bin/csh

parmed << EOF
parm 1cq5.prmtop 
loadCoordinates 1cq5_min.pdb 
energy 
EOF

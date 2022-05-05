#!/bin/csh

parmed << EOF
parm xxx.prmtop 
loadCoordinates xxx_min.pdb 
energy 
EOF

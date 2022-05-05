#!/bin/csh

foreach w (`cat namestr`)
    echo $w
    sed "s/xxx/$w/g" leap.in > leap_mol.in
    sed -i "s/xxx/$w/g" script_parmed.sh
    tleap -f leap_mol.in
    echo "minimization start"
    $AMBERHOME/bin/sander -O -i min.in -o $w\_min.out -c $w.rst7 -p $w.prmtop -r $w\_min.ncrst
    echo "minimization done"
    $AMBERHOME/bin/ambpdb -p $w.prmtop -c $w\_min.ncrst > $w\_min.pdb
    ./script_parmed.sh > energy_dump
    less energy_dump |grep -A 4 "Bond" > $w\_energy 
    rm *inpcrd *rst7 *.ncrst *.out energy_dump *prmtop
    sed -i "s/$w/xxx/g" script_parmed.sh
    
    echo "energy computed"
    
    sed "s/OP1/O1P/g" $w\_min.pdb > $w.star.pdb
    sed -i "s/OP2/O2P/g" $w.star.pdb
    sed -i "s/ H21/1H2 /g" $w.star.pdb
    sed -i "s/ H22/2H2 /g" $w.star.pdb
    sed -i "s/ H41/1H4 /g" $w.star.pdb
    sed -i "s/ H42/2H4 /g" $w.star.pdb
    sed -i "s/ H61/1H6 /g" $w.star.pdb
    sed -i "s/ H62/2H6 /g" $w.star.pdb
    sed -i "s/HO5'/ H5T/g" $w.star.pdb
    sed -i "s/ H5' /1H5\* /g" $w.star.pdb
    sed -i "s/H5'' /2H5\* /g" $w.star.pdb
    sed -i "s/ H2'/1H2\*/g" $w.star.pdb
    sed -i "s/HO2'/2HO\*/g" $w.star.pdb
    sed -i "s/'/\*/g" $w.star.pdb 
    mv $w.star.pdb ../FA4CG
    mv $w\_min.pdb ../FA_minimized
    mv $w\_energy ../Amber_energies
    echo "wrote new * file for CG"
    
end

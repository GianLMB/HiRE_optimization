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
    
    sed  "s/G5/ G/g" $w\_min.pdb > $w.star.pdb
    sed -i "s/C3 / C /g" $w.star.pdb

    mv $w\_min.pdb ../FA_minimized
    mv $w.star.pdb ../FA4CG
    mv $w\_energy ../Amber_energies
end

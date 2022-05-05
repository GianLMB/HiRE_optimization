
#!/bin/basH
#Ask how many ions
    read -p 'How many MG: ' MGvar
    read -p 'How many Na: ' Navar
    read -p 'How many Cl: ' Clvar
    echo
    echo You selected $MGvar MG ions, $Navar Na ions, $Clvar Cl ions for a total of "$(($MGvar+$Navar+$Clvar))" ions.
    echo
#    read -p 'what is the box size? ' boxvar
    echo
 
#Read the number of CG particles
    atomvar=$(grep '1 ' ichain.dat|cut -c 3-5)
    echo the sytem has $atomvar particles
    
    beg=$(($atomvar+$MGvar+1))
    fin=$(($atomvar+$MGvar+$Navar))
    echo $beg $fin
    sed -i -e "${beg},${fin}s/MG/Na/g" test.pdb
    
    beg2=$(($fin+1))
    fin2=$(($fin+$Clvar))
    echo $beg2 $fin2
    sed -i -e "${beg2},${fin2}s/MG/Cl/g" test.pdb
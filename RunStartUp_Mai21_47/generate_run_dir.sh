# generate directory to run RNA-protein simulations.
# generates all necessary files for the RNA part.
# imports exectutable from code source.
# ./generate_run_dir.sh molecule_name (1F9L) source_code_dir (HiRE-OPEP_June2016)

#Ask how many ions
read -p 'How many MG: ' MGvar
read -p 'How many Na: ' Navar
read -p 'How many Cl: ' Clvar
echo
ions=$(($MGvar+$Navar+$Clvar))
echo You selected $MGvar MG ions, $Navar Na ions, $Clvar Cl ions for a total of $ions ions.
echo
read -p 'what is the box size? ' boxvar 
echo

rm -rf Output/$1\_$ions

mkdir $1
cp Input/$1.pdb $1
cp HiRE_parm $1
#cp simulator-RNA-proteins.sh $1
cp Parameter_files/* $1
cp Titfiles/* $1
cp -r SAXSfiles/* $1
cd $1
    
python place_ions.py $1.pdb $ions $boxvar 

./HiRE_parm $1.pdb

atomvar=$(grep '1 ' ichain_RNA.dat|cut -c 3-5)
echo the sytem has $atomvar particles
    
beg=$(($atomvar+$MGvar+1))
fin=$(($atomvar+$MGvar+$Navar))
echo $beg $fin
sed -i -e "${beg},${fin}s/MG/Na/g" conf_initiale_RNA.pdb
    
beg2=$(($fin+1))
fin2=$(($fin+$Clvar))
echo $beg2 $fin2
sed -i -e "${beg2},${fin2}s/MG/Cl/g" conf_initiale_RNA.pdb

python2 charges_differentions.py $MGvar $Navar $Clvar

cp ichain_RNA.dat ichain.dat           # the program does not like ichain_RNA.dat alone, but asks also for ichain.dat. Check.
#mv parametres.top parametres_RNA.top


sed -i  "s/ MG / Na /$(($MGvar+1))g" parametres_RNA.top
echo $(($MGvar+1))

sed -i  "s/ Na / Cl /$(($Navar+1))g" parametres_RNA.top
echo $(($MGvar+$Navar+1))

 
cd ../$2
cp simulator* ../$1
pwd > ../$1/code_version.txt
date >> ../$1/code_version.txt
cd ../

mv $1 Output/$1\_$ions

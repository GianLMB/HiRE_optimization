# generate directory to run RNA-protein simulations.
# generates all necessary files for the RNA part.
# imports exectutable from code source.
# ./generate_run_dir.sh molecule_name (1F9L) source_code_dir (HiRE-OPEP_June2016)

rm -rf Output/$1

mkdir $1
cp Input/$1.pdb $1
cp scale_RNA.dat $1

cd $1
./../SCRIPTS/HiRE_parm $1.pdb
touch protonated.dat
python ../SCRIPTS/charges_differentions.py 0 0 0
python ../SCRIPTS/HiREtopology/parse_oldhirefile.py
rm $1.top baselist.dat  chargeatm_RNA.dat conf_initiale_RNA.pdb ichain_RNA.dat parametres.csh parametres_RNA.top

mkdir HiRE-test
cp parameters.top scale_RNA.dat HiRE-test
mv bblist.dat protonated.dat HiRE-test
cp ../SCRIPTS/HIRE HiRE-test
awk '{print $7 "\t" $8 "\t" $9}' $1_CG.pdb > HiRE-test/start 

cd ../

mv $1 Output/$1

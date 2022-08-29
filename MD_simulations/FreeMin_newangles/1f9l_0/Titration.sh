
   here=$(pwd)

   number_bases=$(sed "1d" on-the-fly.pdb | wc -l)
#  Run titration
   pH=$(grep "setenv pH" simulator-RNA-proteins.sh | sed "s/setenv pH//" | sed "s/  * //")
#   echo "pH is ${pH}"
   sed "s/:pH_titra/:${pH}/" parametres_titration.dat > cyl.json 
   rm -fr state
   ./faunus.exe > titration.eq
   ./faunus.exe > titration.out
#  Write out results for Samuela's code
   grep -A ${number_bases} "Site            ⟨z⟩       Acceptance" titration.out | tail -${number_bases} | sed "s/ * / /g" | cut -d" " -f 4 > charges_from_titration_RNA.dat 
    

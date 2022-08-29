### Execute generate_run_dir.sh on all files in Input directory

import subprocess
import os

filelist = []
path = '/home/flechenault/Documents/Gianluca/Structure_DB_Amber_HiRE/Prep_pureHire/Input'
for f in os.listdir(path):
    if f.endswith(".pdb"):
        if os.path.getsize(path+'/'+f) > 0:  # check file is not empty
            filelist.append(f.split('.')[0])
print(filelist)

for file in filelist:
    subprocess.check_call(['./generate_run_dir.sh '+file], shell=True)
    print(f'Subprocess {file} executed')

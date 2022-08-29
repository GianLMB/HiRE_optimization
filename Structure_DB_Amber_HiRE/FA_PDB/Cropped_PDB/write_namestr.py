import os

filelist = []
for f in os.listdir():
    if f.endswith(".pdb"):
        filelist.append(f.split('.')[0])

with open('namestr2','w') as f:
    for file in filelist:
        f.write(file)
        f.write('\n')

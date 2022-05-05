import os
import sys
import numpy
import math

input_pdb = "{}.pdb".format(sys.argv[1]) #atomistic pdb
saxs_calc_type = sys.argv[2] #"vacuum", "solution" or "exvol"

if saxs_calc_type == "vacuum":
	cmd = "crysol {} -sm 0.995 -ns 200 -dro 0 -dns 0 -eh".format(input_pdb)
elif saxs_calc_type == "solution":
	cmd = "crysol {} -sm 0.995 -ns 200 -dro 0.03 -dns 0.334 -eh".format(input_pdb)
else:
	cmd = "crysol {} -sm 0.995 -ns 200 -dro 0 -dns 0.334 -eh".format(input_pdb)

print(cmd)

os.system(cmd)

q=numpy.arange(0,1,0.005)

n=0

with open ("saxs_target.dat","w") as filout:
	with open ("{}00.int".format(sys.argv[1])) as filin:
		lignes=filin.readlines()
		for ligne in lignes[1::]:		
			filout.write("   {:5.3f}  {:9.7f}E+00\n".format(q[n],math.log(float(ligne.split()[1]),10)))
			n=n+1

os.system("rm {}00.log {}00.alm".format(sys.argv[1],sys.argv[1]))

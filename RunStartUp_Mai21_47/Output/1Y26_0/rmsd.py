import mdtraj as md
import numpy as np
import sys
import subprocess
import math
import sys
import pylab
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

trjname = str(sys.argv[1])
pdbname1 = str(sys.argv[2])

strref1 = md.load(pdbname1)
traj = md.load(trjname, top=pdbname1)

rmsd1=md.rmsd(traj, strref1, 0) 
rmsdA=rmsd1*10.0

#pylab.plot(rmsd1, '-', label=pdbname1) 
#pylab.plot(rmsd2, '-', label=pdbname2) 
#plt.title(trjname, fontsize = 24,x=0.5,y=1.02)
#pylab.savefig(trjname+'_rmsd.pdf')
#pylab.show()

#dd0=md.compute_angles(traj,[idxs[0]],periodic=False)

print(np.argmin(rmsdA), np.min(rmsdA))
print(np.argmax(rmsdA), np.max(rmsdA))


plt.plot(rmsdA)
plt.ylim((0,np.max(rmsdA)+2))
plt.show()

pylab.savefig(trjname+'_rmsd.pdf')

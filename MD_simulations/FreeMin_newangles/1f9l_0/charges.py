from string import split
from math import *
import sys

g = open('baselist.dat','r')
p = g.readlines()
ch=[]
for i in p:
    ch.append(split(i))

h = open('protonated.dat','r')
q = h.readlines()
pr=[]
for i in q:
    if(q != 0):
        pr.append(split(i))
    

L = len(ch)
for i in range(L):
    ch[i].append('0')
    for j in range(len(pr)):
        if i == int(pr[j][0])-1:
          ch[i][2]='1'

f = open('bblist.dat','w')   
s=' '
for i in ch:
    f.write(s.join(i))
    f.write('\n')

L = len(ch)
q=[]
for i in range(int(ch[L-1][0])):
    q.append(0)

for i in range(L-1):
#    print int(ch[i][0])
    q[int(ch[i][0])]=-1
#    print ch[i][1],ch[i][2]
    if (ch[i][2]=='1'):
	if (ch[i][1]=='1' or ch[i][1]=='4'):
             q[int(ch[i][0])-1]=-1	
        if (ch[i][1]=='2' or ch[i][1]=='3'):
             q[int(ch[i][0])-1]=1
    if (ch[i][1]=='0'):
        q[int(ch[i][0])]=2

#for i in q:
    #print i

s=' '
y = open('chargeatm_RNA.dat', 'w')
for i in range(int(ch[L-1][0])):
    ll = [str(i+1), str(q[i])]
#    print s.join(ll)
    y.write(s.join(ll))
    y.write('\n')

#b = open(sys.argv[1],'r')
#db = b.readlines()
#pdb=[]
#for i in db:
    #pdb.append(split(i))


## write input file for widom calculations
#for i in range(int(ch[L-1][0])):
    #print "{:4s}{:7.3f}{:9.3f}{:9.3f}{:9.3f}".format('ATOM',q[i],float(pdb[i][6]),float(pdb[i][7]),float(pdb[i][8]))


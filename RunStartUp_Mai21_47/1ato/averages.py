import re, string
import sys, os

f = open('titration_MD.out','r')
tit = f.readlines()

L = len(string.split(tit[0]))
N = len(tit)
start = int(sys.argv[1])


v=[]
for i in range(L):
    v.append(0)

for j in range(start, N):
    a=(string.split(tit[j]))
    for i in range(L):
        v[i]=v[i]+int(a[i])
        
        
for i in range(L):
    print i+1, v[i]/float(N-start)
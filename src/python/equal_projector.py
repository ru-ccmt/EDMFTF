#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from scipy import *
import sys
from pylab import *

if len(sys.argv)<2:
    print 'Give the name of the projector file'
    sys.exit(0)

pfile = sys.argv[1]
fi=open(pfile,'r')
line1 = fi.readline()
Np, Nr = map(int,line1[1:].split())

proj=[]
lines=[]
for p in range(Np):
    line2 = fi.readline()
    lines.append(line2)
    Nr = int(line2[1:].split()[0])
    dat=[]
    for i in range(Nr):
        line = fi.readline()
        dat.append( map(float,line.split()) )
    proj.append(dat)
proj = array(proj)
fi.close()

projf = zeros((shape(proj)[1],shape(proj)[2]))
for i in range(Np):
    projf[:,:] += proj[i,:,:]
projf *= 1.0/Np


fo = open(pfile+'_new', 'w')
print >> fo, line1,

for p in range(Np):
    print >> fo, lines[p],
    for i in range(Nr):
        print >> fo, projf[i,0], projf[i,1], projf[i,2]
fo.close()
#plot(projf[:,0],projf[:,1])
#plot(projf[:,0],projf[:,2])
#show()






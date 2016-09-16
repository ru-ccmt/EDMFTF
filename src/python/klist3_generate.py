#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
""" Generates k-list for w2k for plotting. Needs to set "legend" and "BS"
"""
from scipy import *
import sys
from fractions import gcd


legend=['Gamma','H','P','N','Gamma','P']
BS=[[0,0,0], [0,1,0], [1/2.,1/2.,1/2.],[1/2.,1/2.,0],[0,0,0],[1/2.,1/2.,1/2.]]


#####
BS = array(BS,dtype=float)
Nt = 500
dl=zeros(len(BS)-1)
for i in range(len(BS)-1):
    dr=BS[i+1]-BS[i]
    dl[i] = sqrt(dot(dr,dr))
dN = dl/sum(dl)*Nt
Ni = [int(round(dn)) for dn in dN]

Mi = zeros(len(Ni),dtype=int)
for i in range(len(Ni)):
    fracts = 1/array(filter(lambda x:x>1e-6, hstack((BS[i+1],BS[i]))))
    fact = int( reduce(lambda x,y: x*y/gcd(x,y), fracts) )
    if Ni[i] % fact==0:
        Mi[i]=Ni[i]
    else:
        Mi[i]=Ni[i]*fact

l=0
for p in range(len(Ni)):
    NAME=legend[p]
    for i in range(Ni[p]):
        l+=1
        kint = Mi[p]*BS[p] + (BS[p+1]-BS[p])*i*Mi[p]/Ni[p]
        if i>0: NAME = '   '#+str(l)
        print "%-10s%5d%5d%5d%5d%5.1f" % (NAME, kint[0], kint[1], kint[2], Mi[p], 1.0)
NAME = legend[-1]
kint=BS[-1]*Mi[-1]
print "%-10s%5d%5d%5d%5d%5.1f" % (NAME, kint[0], kint[1], kint[2], Mi[-1], 1.0)
print 'END'



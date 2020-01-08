#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
""" Generates k-list for w2k for plotting. Needs to set "legend", "BS" and "BR2"
"""
from scipy import *
import sys
from fractions import gcd
from scipy import linalg

BR2 = None
######################################
##    Edit this part  only          ##
######################################
Nt = 500                                                  # number of all points along the path
legend=['Z','Gamma','X/2','Y/2','Gamma']                  # name of  the points
#BR2=[[0.94919,  -0.25835,   0.00000],  # reciprocal vectors from case.rotlm
#     [0.94919,   0.25835,   0.00000],  # absolute value is not important, only ratios are important
#     [0.00000,   0.00000,   0.21209]]  # Needed to find the length of each interval, and redistribute points along the path
BR1=[[0.94919,   0.00000,   0.00000],
     [0.00000,   0.25835,   0.00000],
     [0.00000,   0.00000,   0.21209]]
BS = [[0.025, 0, 0], [0.025, 1, 0], [0.1,1,0],[0,0,0]]   # CAREFUL: Has to be specified in conventional brillouine zone, not primitive
########################################
if BR2 is not None:
    Prim2Conv = dot(linalg.inv(array(BR1).T), array(BR2).T)
BS = array(BS,dtype=float)

lBS = dot(BS,array(BR1))
#print  'lBS=', lBS
dl=zeros(len(lBS)-1)
for i in range(len(lBS)-1):
    dr=lBS[i+1]-lBS[i]
    dl[i] = sqrt(dot(dr,dr))
    #print 'r1=', lBS[i],'r2=', lBS[i+1], 'dr=', dl[i]
dN = dl/sum(dl)*Nt
Ni = [int(round(dn)) for dn in dN]
#print 'Ni=', Ni

Mi = zeros(len(Ni),dtype=int)
for i in range(len(Ni)):
    fracts = 1/array(filter(lambda x:x>1e-6, hstack((BS[i+1],BS[i]))))
    fact = int( reduce(lambda x,y: x*y/gcd(int(x),int(y)), fracts) )
    Mi[i]=Ni[i]*fact


Klist = []
Name = []
l=0
for p in range(len(Ni)):
    for i in range(Ni[p]):
        l+=1
        kint = Mi[p]*BS[p] + (BS[p+1]-BS[p])*i*Mi[p]/Ni[p]
        kn = (int(kint[0]), int(kint[1]), int(kint[2]), Mi[p])
        if i==Ni[p]-1 and p < len(Ni)-1:
            # we want to avoid the situation where the points are identical, the first in the new interval, and the last in the previous interval
            knext = array(BS[p+1]*Mi[p+1],dtype=int)/float(Mi[p+1]) # first point in the new interval
            kthis = array(kn[:3])/float(Mi[p])                      # last in the current interval
            diff = sum(abs(knext-kthis))                            # how much difference
            #print 'knext=', knext, 'kthis=', kthis, 'diff=', diff
            if diff < 1e-6:                                         # the two k-points are basically identical
                continue                                            # hence skip this last point, as the begining of next interval will be identical
        if i==0:                                                    # first point of each interval is always used to print its name
            Klist.append( kn )
            Name.append(legend[p])                                  # remember name
        elif kn != Klist[-1]:                                       # store only if this point is different from previous
            Klist.append( kn )
            Name.append('   ')                                      # part of the interval, not high-symmetry point with name
kint = Mi[-1]*BS[-1]                                                # very last point
kn = (int(kint[0]), int(kint[1]), int(kint[2]), Mi[-1])             # very last point
if Klist[-1]==kn:                                                   # previous point was already identical
    Name[-1] = legend[-1]                                           # than we just rename the previous point, and not adding anything
else:
    Klist.append( kn )                                              # adding last point
    Name.append(legend[-1])                                         # with its proper name

for ik in range(len(Klist)):                                        # finally just print all obtained points
    print "%-10s%5d%5d%5d%5d%5.1f" % (Name[ik], Klist[ik][0], Klist[ik][1], Klist[ik][2], Klist[ik][3], 1.0)
print 'END'



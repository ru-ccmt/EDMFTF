import sys, re, os
from scipy import *
#from pylab import *
import string
from scipy.lib.blas import fblas
import rdVec, utils, struct1

if len(sys.argv)<2:
    print 'Please give the name of vector file[s]'
    sys.exit(0)
else:
    fnames = sys.argv[1:]

# Vector-file handle
tapes=array(range(len(fnames)))+8


w2k = utils.W2kEnvironment()
struct = struct1.Struct(w2k.case)

maxkpoints = 10000
# opens vector file

heads=[]
Eks=[]
for fname,tape in zip(fnames,tapes):
    
    (Elinear1,Elinear2) = rdVec.feopen(tape, fname, struct.nat)
    #Elinear1 = [(string.join(strng,'')).rstrip() for strng in Elinear1]
    #Elinear2 = [(string.join(strng,'')).rstrip() for strng in Elinear2]

    print 'linearization energy 1 =', Elinear1
    print 'linearization energy 2 =', Elinear2

    for ik in range(maxkpoints):
        # Reads vector file
        head = rdVec.feread1(tape)
        (k, kname, wgh, ios, n0, nb) = head
        if ios!=0: break # vector file is finished, no more k-points
        # Eigenvalues
        Ek = rdVec.feread2(tape, nb)

        heads.append(head)
        Eks.append(Ek)
        print 'k=', k
        
    rdVec.feclose(tape)
    

fh_final=8
fout = w2k.case+'.energy_out'
rdVec.fweopen(fh_final, fout, Elinear1, Elinear2)    

for head,Ek in zip(heads,Eks):
    (k, kname, wgh, ios, n0, nb) = head
    ios = rdVec.fewrite1(fh_final,k,kname,wgh,n0,nb)
    rdVec.fewrite2(fh_final,Ek)
rdVec.feclose(fh_final)

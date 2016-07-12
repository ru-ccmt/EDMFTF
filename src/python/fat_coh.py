#!/usr/bin/env python

from scipy import *
from scipy import linalg
import struct, sys, re
import rdU, utils
import findEF

if len(sys.argv)<=1:
    print 'Give dictionary of coh_orb'
    sys.exit(0)

#coh_orb = {0:2,1:2}
coh_orb = eval(sys.argv[1])
print 'coh_orb=', coh_orb




fhp = 233
cpuID = 0

dmfe = utils.DmftEnvironment()  # DMFT paths
w2k = utils.W2kEnvironment()    # W2k filenames and paths
case = w2k.case
(EF,NOE) = findEF.FindChemicalPotential(case,'')

# Read Udmft
filename = 'Udmft.'+str(cpuID)
(nkp, nsymop, norbitals) = rdU.fopen(fhp, filename)
nindo = rdU.read1(fhp, norbitals)

print 'nsymop=', nsymop

fc = open('cohfactors.dat', 'w')

print >> fc, '# Coherence factors for KS states'

for ikp in range(nkp):

    (iikp, nbands, tmaxdim2, tnorbitals, nemin) = rdU.read2(fhp)
    
    print 'ikp=', ikp,  'nbands=', nbands, 'maxdim2=', tmaxdim2, 'norbitals=', tnorbitals, 'nemin=', nemin
    
    for isym in range(nsymop):

        iisym = rdU.read3(fhp)
        if iisym!=isym+1:
            print 'ERROR: isym!=iisym', isym, iisym
        
        gs = zeros((nbands,nbands),dtype=complex)

        for iorb in range(norbitals):
            # Reading DMFT transformation UDMFT
            U = []
            for ind in range(nindo[iorb]):
                U.append( rdU.read4(fhp, nbands) )
            U = array(U)
            # print ikp, isym, iorb, shape(U) 
            
            if (isym==0):
                # Prepares coherence factors gs(i,j)=U(orb,i)*U(orb,j).conj()
                if coh_orb.has_key(iorb):
                    xorb = coh_orb[iorb]
                    for iband in range(nbands):
                        gs[iband,:] += U[xorb,iband]*conj(U[xorb,:])
                        
        if isym==0:
            # Printing the coherence factors
            print >> fc, iikp, nemin, nbands
            for i in range(len(gs)):
                for j in range(len(gs)):
                    print >> fc, "%15.10f " % abs(gs[j,i]),
                print >> fc




rdU.fclose(fhp)
            

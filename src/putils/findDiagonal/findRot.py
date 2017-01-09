#!/usr/bin/env python
import os, sys
import shutil
from scipy import *
import utils
import indmffile
import find3dRotation as f3dr
import find5dRotation as f5dr

if __name__ == '__main__':

    hybrid_infty = False  # whether to diagonalize hybridization in infinity or at zero
    
    w2k = utils.W2kEnvironment()     # W2k filenames and paths


    outdmft1 = w2k.case+'.outputdmf1'
    if not os.path.isfile(outdmft1):
        print 'Can not find file ', outdmft1, 'therefore can not correct local orbitals. Please run "dmft1" first, preferably on the real axis.'
        sys.exit(1)

    # Reads case.indmfl file
    inl = indmffile.Indmfl(w2k.case) # case.indmfl file
    inl.read()                       # case.indmfl read

    f1 = open(outdmft1, 'r')
    Found=False
    for line in f1:
        if line.strip()=='Full matrix of impurity levels follows':
            #print 'Found:', line
            Found=True
            break
    if not Found:
        print 'Could not find impurity levels in', outdmft1, '. You should consider to rerun "dmft1", preferably on the real axis.'
        sys.exit(1)

    Hc={}
    for icix in inl.cps.keys():
        line = f1.next().split()
        if int(line[1])!=icix:
            print 'It seems that '+w2k.case+'.indmfl and '+outdmft1+' are incompatible...'
            sys.exit(1)
        Sigind = inl.siginds[icix]
        cols = filter(lambda x: x!=0,Sigind.diagonal())
        dim = len(cols)
        
        Hc[icix] = zeros( (dim,dim), dtype=complex )
        for idim in range(dim):
            line = f1.next()
            if hybrid_infty:
                dat = array(map(float,line.split()))
                Hc[icix][idim,:] = dat[::2] + dat[1::2]*1j
            #print 'x: ', line,
        line = f1.next()
        line = f1.next()
        for idim in range(dim):
            line = f1.next()
            if not hybrid_infty:
                dat = array(map(float,line.split()))
                Hc[icix][idim,:] = dat[::2] + dat[1::2]*1j
            #print 'y: ', line,
        line = f1.next()
        
        #print Hc[icix]
        #print inl.cftrans[icix]
        l = inl.cps[icix][0][1]
        
        if len(Sigind) == (2*l+1):
            so=False
        elif len(Sigind)==(2*l+1)*2:
            so=True
        else:
            print ':ERROR : Sigind has dimenion', len(Sigind), 'while l=', l
            sys.exit(1)
        
        hc = 0.5*(Hc[icix]+transpose(conjugate(Hc[icix]))) # Should be Hermitian, so we enforce Hermicity
        T2C = inl.cftrans[icix][:len(hc),:]
        T2Crest = inl.cftrans[icix][len(hc):,:]

        if so:
            T2Cnew = f5dr.GiveNewT2C(hc, T2C)
        else:
            T2Cnew = f3dr.GiveNewT2C(hc, T2C)
        
        inl.cftrans[icix][:len(hc),:] = T2Cnew[:,:]
    shutil.copy(w2k.case+'.indmfl', w2k.case+'.indmfl_findRot')
    inl.write(w2k.case+'.indmfl', only_write_stored_data=True)

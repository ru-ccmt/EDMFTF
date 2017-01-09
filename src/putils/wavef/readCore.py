#!/usr/bin/env python

from scipy import *
import utils
import sys, os
import struct1
import re
from w2k_atpar import readpotential, readlinearizatione, atpar, rint13, rint13g
from pylab import *

def read_core(case, struct):
    # reading case.inc
    fc = open(case+'.inc', 'r')
    n_kappa_occup =[]
    iprint=[]
    for jatom in range(struct.nat):
        line = fc.next()
        data = line.split()
        nc = int(data[0])
        iprint.append( int(data[2]) )
        #nc = int(line[:3])
        nko = zeros((nc,3),dtype=int)
        for il in range(nc):
            line = fc.next()
            n, kappa, occup = (int(line[:1]), int(line[2:4]), int(line[5:6]))
            nko[il,:] = array([n,kappa,occup])
        n_kappa_occup.append(nko)
        
    print 'iprint=', iprint
    
    # reading case.corewf
    fi = open(case+'.corewf','r')
    wavec=[]
    cType=[]
    cEne=[]
    for jatom in range(struct.nat):
        line = fi.next().split()
        nc0 = int(line[0])
        nc = len(n_kappa_occup[jatom])
        noccup0 = n_kappa_occup[jatom][0,2]
        #print noccup0, nc, nc0
        waveci=[]
        cTypei=[]
        cEnei=[]
        if (noccup0 and iprint[jatom]):
            for ic in range(nc):
                line = fi.next()
                m1 = re.search('CORE STATES = (\d\S)', line)
                m2 = re.search('CORE ENERGY=(.*) Ry', line)
                if m1 is not None and m2 is not None:
                    ctype = m1.group(1)
                    cenergy = float(m2.group(1))

                ucore = zeros((2,struct.npt[jatom]))
                for i in range(2):
                    n=0
                    while (n<struct.npt[jatom]):
                        line = fi.next()
                        dat = [line[3:3+19], line[3+19:3+2*19], line[3+2*19:3+3*19], line[3+3*19:3+4*19]]
                        dat = filter(lambda s: len(s)>1, dat) # throw away some empty strings
                        dat = map(float,dat)
                        ucore[i,n:n+len(dat)]= dat[:]
                        n += len(dat)
                waveci.append(ucore)
                print ctype, cenergy, shape(ucore)
                cTypei.append( ctype )
                cEnei.append(cenergy)
                
        wavec.append(waveci)
        cType.append(cTypei)
        cEne.append(cEnei)

    return (wavec, cType, cEne, iprint, n_kappa_occup)

if __name__ == '__main__':
    w2k = utils.W2kEnvironment()
    struct = struct1.Struct(w2k.case)
    (wavec, cType, cEne, iprint, n_kappa_occup) = read_core(w2k.case, struct)

    for jatom in range(struct.nat):
        if iprint[jatom]:
            for ic in range(len(cType[jatom])):
                #jatom,ic = 0,4
                print cType[jatom][ic]
                print shape(wavec[jatom][ic])
                
                Nr = struct.npt[jatom]  # Number of radial points from struct file
                r0 = struct.r0[jatom]
                Rmt = struct.rmt[jatom]
                dx = log(Rmt/r0)/(Nr-1)
                Rx = r0 * exp( arange(Nr)*dx ) # radial points
                Ag = wavec[jatom][ic][0]
                Bg = wavec[jatom][ic][1]
                norm = rint13g(1.0,1.0,Ag,Bg,Ag,Bg,dx,r0)
                
                plot(Rx, Ag, label=cType[jatom][ic])
                savetxt('corewf_'+str(jatom)+'_'+str(ic)+'.dat', transpose(vstack((Rx,Ag,Bg))) )
                print cType[jatom][ic], 'norm=', norm
    title('Core wave fnunctions')
    legend(loc='best')
    show()

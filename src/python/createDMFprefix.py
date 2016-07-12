#!/usr/bin/env python
from scipy import *
#from pylab import *
import re
import sys

def FindNCpu(Nk,Ncpu_max):
    for Ncpu in range(Ncpu_max,1,-1):
        if (Nk-int(Nk/Ncpu)*Ncpu < Nk/Ncpu):
            Nrest = Nk-(Nk/Ncpu)*Ncpu
            break
    if Nrest>0 and Ncpu==Ncpu_max:
        for Ncpu in range(Ncpu_max-1,1,-1):
            if (Nk-int(Nk/Ncpu)*Ncpu < Nk/Ncpu):
                Nrest = Nk-(Nk/Ncpu)*Ncpu
                break
        return (Ncpu, Nrest)
    else:
        return (Ncpu, Nrest)

def GetMachines(fmachine):
   fm = open(fmachine,'r')
   lines=fm.readlines()
   fm.close()
   machns=[]
   for line in lines:
       if line.strip():
           machns.append(line.strip())
           #print machns[-1]
   fm.close()
   return machns

def GetNumberOfKpoints(fklist):
    fk = open(fklist,'r')
    Nk=0
    for line in fk:
        if line[:3]=='END': break
        Nk += 1
    fk.close()
    return Nk

def FindBestDistributionCores(machns, Nk):
    mold = machns[:]
    mnew = []

    
    OMP = int(round(len(machns)/float(Nk)))

    #print 'cpus=', len(machns), 'Nk=', Nk, 'OMP=', OMP
    
    if OMP>=16:  # 1 jobs per cpu (16 cores)
        OMP=16
    elif OMP>=8: # 2 jobs per cpu
        OMP=8
    elif OMP>5:  # 3 jobs per cpu
        OMP=5
    elif OMP>=4: # 4 jobs per cpu
        OMP=4
    elif OMP>3:  # 5 jobs per cpu
        OMP=3
    elif OMP>=2:  # 8 jobs per cpu
        OMP=2
    else:
        OMP=1    # 7-8 jobs per cpu
    
    while  len(mnew)<Nk :
        Of = len(mold)/float(Nk-len(mnew))
        Om = int(round(Of))
        #print 'Of=', Of, 'Om=', Om
        if Om>0:
            add = mold[::Om]
        else:
            add = mold[:]
        mnew = mnew + add
        for i in range(len(add)): mold.remove(add[i])
        #print 'mold=', mold
        #print 'mnew=', mnew
        if not mold: break  # mold empty

    #print 'OMP=', OMP
    #print 'mnew[:Nk]=', mnew
    return (mnew[:Nk], OMP)
    

if __name__ == '__main__':
    
    if (len(sys.argv)<2):
        print 'Give two arguments: case.klist, mpi_prefix.dat!'
        sys.exit(0)

    machine_out='wmachines'
    
    fklist=sys.argv[1]
    fprefix=sys.argv[2]
    mpi_prefix = open(fprefix).next().strip()

    # Finds the name of the macginefile from mpi_run.dat directive
    m = re.search('-machinefile\s+([\w|.|/]*)',mpi_prefix)
    if m is not None:
        fmachine = m.group(1)
    

    Nk = GetNumberOfKpoints(fklist)
    machns = GetMachines(fmachine)
    Ncpu_max = len(machns)

    newmach, OMP = FindBestDistributionCores(machns, Nk)

    if OMP<=1:
        OMP = 1
        mpin = mpi_prefix 
    else:
        # creates a new machinefile (wmachines) with fraction of cores in the list
        fm = open(machine_out,'w')
        for i in range(len(newmach)):
            print >> fm, newmach[i]
        fm.close()
        # replaces -machinefile directive in mpi_prefix.dat with the name of the new machinefile
        mpin = re.sub('-machinefile\s+([\w|.|/]*)', '-machinefile '+machine_out, mpi_prefix)
        # reduced the number of cores for MPI run (-np directive) in mpi_prefix.dat
        m = re.search('-np\s*(\d+)',mpin)
        if m is not None:
            np = int(m.group(1))
            mpin = re.sub('-np\s*(\d+)', '-np '+str(len(newmach)), mpin)
        else:
            mpin += '-np '+str(len(newmach))+' '
    
    m = re.search('-env\s*OMP_NUM_THREADS\s*(\d+)',mpin)
    if m is not None:
        mpin = re.sub('-env\s*OMP_NUM_THREADS\s*\d+', '-env OMP_NUM_THREADS '+str(OMP), mpin)
    else:
        mpin += '-env OMP_NUM_THREADS '+str(OMP)

    print mpin
    

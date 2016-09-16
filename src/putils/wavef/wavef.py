#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule

import os
from scipy import *
import utils
from w2k_atpar import readpotential, readlinearizatione, atpar, rint13
import struct1
from pylab import *
import optparse
import re

def SolveForContinuousFunction(A,Ae,Aee,Rx,Nr0,Nr):
    """
    Routine extends solution beyond Rmt by making value and derivative continuous and make it vanish ar R2
    """
    dA = (A[Nr0]-A[Nr0-1])/(Rx[Nr0]-Rx[Nr0-1])    # Radial derivative of ul(r)
    dAe= (Ae[Nr0]-Ae[Nr0-1])/(Rx[Nr0]-Rx[Nr0-1])  # Radial derivative of dot{ul(r)}
    dAee=(Aee[Nr0]-Aee[Nr0-1])/(Rx[Nr0]-Rx[Nr0-1])# Radial derivative of dot{dot{ul(r)}}

    am=array([[A[Nr0-1], Ae[Nr0-1], Aee[Nr0-1]],  # 3x3 equation for continuous solution
              [A[Nr-1],  Ae[Nr-1],  Aee[Nr-1]],   # for r<Rmt, we take solution of Dirac equation
              [dA, dAe, dAee]])                   # outside we use cobination of ul(r),dot{ul(r)},dot{dot{ul(r)}}
    bm=array([A[Nr0-1], 0., dA])                  # 3 Eqs.: value at Rmt, value at R2, derivative at Rmt
    cm=linalg.solve(am,bm)

    return cm

def SolveForLocalizedFunction(A,Ae,Aee,Rx,Nr):
    dA = (A[Nr-1]-A[Nr-2])/(Rx[Nr-1]-Rx[Nr-2])    # Radial derivative of ul(r)
    dAe= (Ae[Nr-1]-Ae[Nr-2])/(Rx[Nr-1]-Rx[Nr-2])  # Radial derivative of dot{ul(r)}
    dAee=(Aee[Nr-1]-Aee[Nr-2])/(Rx[Nr-1]-Rx[Nr-2])# Radial derivative of dot{dot{ul(r)}}

    am=array([[ Ae[Nr-1], Aee[Nr-1]],
              [ dAe,     dAee     ]])
    bm=array([-A[Nr-1],-dA])
    cm=linalg.solve(am,bm)
    return cm

def SolveForLocalizedFunction3(A,Ae,Nr):
    return -A[Nr-1]/Ae[Nr-1]

def ReadIndmfl(case):
    def divmodulo(x,n):
        "We want to take modulo and divide in fortran way, so that it is compatible with fortran code"
        return ( sign(x)* (abs(x)/n) , sign(x)*mod(abs(x),n))

    findmfl = open(case+'.indmfl', 'r')
    lines = [line.split('#')[0].strip() for line in findmfl.readlines()] # strip comments
    findmfl.close()
    lines = (line for line in lines if line)  # strip blank lines & create generator expression
    
    dat = lines.next().split()[:4]
    hybr_emin, hybr_emax, Qrenorm, projector = (float(dat[0]),float(dat[1]),int(dat[2]),int(dat[3]))
    dat = lines.next().split()[:6]
    matsubara, broadc, broadnc, om_npts, om_emin, om_emax = (int(dat[0]),float(dat[1]),float(dat[2]),int(dat[3]),float(dat[4]),float(dat[5]))
    natom = int(lines.next())

    Rmt2 = zeros(natom,dtype=float)
    (atms, Lsa, qsplita, icpsa) = ([],[],[],[])
    for ia in range(natom):
        dat = lines.next().split()
        iatom, nL, locrot_shift = map(int, dat[:3])
        atms.append(iatom)
        if len(dat)>3: Rmt2[ia] = float(dat[3])
        
        (shift,locrot) = divmodulo(locrot_shift,3)
        if locrot<0: 
            if locrot==-2: 
                locrot=3*nL
            else:
                locrot=3

        
        (Ls, qsplit, icps) = (zeros(nL,dtype=int), zeros(nL,dtype=int), zeros(nL,dtype=int))
        for il in range(nL):
            (Ls[il], qsplit[il], icps[il]) = map(int, lines.next().split()[:3])

        Lsa.append( Ls )
        qsplita.append( qsplit )
        icpsa.append( icps )
        
        new_zx = [[float(x) for x in lines.next().split()[:3]] for loro in range(abs(locrot))]
        vec_shift = [float(x) for x in lines.next().split()] if shift else None
        
    print 'atms=', atms
    print 'Rmt2=', Rmt2
    print 'Lsa=', Lsa
    print 'qsplita=', qsplita
    print 'icpsa=', icpsa
    return (atms, Lsa, icpsa, Rmt2)
        
def FindChemicalPotential(case, updn):
    # Looking for the LDA chemical potential                                                                                                                                                                
    Ry2eV = 13.60569193
    EF_found = False

    fname = case+".scf2"
    if os.path.isfile(fname):
        fscf = open(fname, 'r')
        print "opening", fname
    elif os.path.isfile(fname+updn):
        fscf = open(fname+updn, 'r')
        print "opening", fname+updn
    else:
        fscf = None

    if fscf is not None:
        lines = fscf.readlines()
        for line in lines:
            if re.match(r':FER', line) is not None:
                EF = float(line[38:])
                print 'EF=', EF*Ry2eV, 'eV = ', EF, 'Ry'
                EF_found = True

    # The previous DMFT chemical potential                                                                                                                                                                  
    if  os.path.isfile('EF.dat'):
        fmu = open('EF.dat','r')
        mu = float(fmu.next())
        EF = mu/Ry2eV
        print 'Found DMFT-EF=', EF*Ry2eV, 'eV'
        EF_found = True

    if not EF_found:
        fname = case+".scf"
        if os.path.isfile(fname):
            print "opening", fname
            f = open(fname, 'r')
            lines = f.readlines()
            if not EF_found:
                for line in lines:
                    if re.match(r':FER', line) is not None:
                        EF = float(line[38:])
                        print 'EF=', EF*Ry2eV
                        EF_found = True

    if not EF_found:
        raise Exception("Failed to determine chemical potential.")

    return EF

def main(case, atms, Lsa, icpsa, Rm2, localize=1, Emu='EF'):

    so=''
    if os.path.isfile(case+".inso") and os.path.getsize(case+".inso")>0 :
        print 'Found '+case+'.inso file, hence assuming so-coupling exists. Switching -so switch!'
        so='so'
        
    potential_file = case+'.vsp'
    vector_file = case+'.vector'+so

    struct = struct1.Struct(case)
    rel=True
    if struct.mode[:4]=='NREL': rel=False
    
    
    atm_l_case={}
    latom=0
    for jatom in range(struct.nat):
        ll_case=[]
        if (latom+1) in atms:           # This atom listed in indmfl file
            ia = atms.index(latom+1)    # It appears as the ia consequitive atom in indmfl file
            for il in range(len(Lsa[ia])):
                if icpsa[ia][il]>0: # icix>0, hence it is correlated
                    ll_case.append(Lsa[ia][il])
        latom += struct.mult[jatom]
        if ll_case:
            atm_l_case[jatom]= (ll_case, Rm2[ia])
        
    print 'atm_l_case=', atm_l_case

    Nrmax=0
    for jatom in atm_l_case.keys():
        Rmt2 = atm_l_case[jatom][1]
        for lc in atm_l_case[jatom][0]:
            Nr0 = struct.npt[jatom]  # Number of radial points from struct file
            Rmt = struct.rmt[jatom]
            if Rmt2>Rmt:
                # creating larger mesh to Rmt2, which can be larger then Rmt
                r0 = struct.r0[jatom]
                if Rmt2<Rmt: Rmt2 = Rmt 
                Nr = 1 + int( (Nr0-1)*log(Rmt2/r0)/log(Rmt/r0) + 0.99)  # Number of radial points after extending mesh to Rmt2
            else:
                Nr=Nr0
            if Nr>Nrmax: Nrmax=Nr

    lmax2 = 3

    Vr_all = readpotential(potential_file, max(struct.npt), struct.npt, struct.nat)               # Reading actual potential
    
    if Emu!='EF':
        E_all = readlinearizatione(vector_file, struct.nat, lmax2)      # Linearization energies from vector file
    else:
        updn=''
        EF = FindChemicalPotential(case, updn)
    
    fout = open('projectorw.dat', 'w')
    print >> fout, '#', sum([len(atm_l_case[jatom][0]) for jatom in atm_l_case.keys()]), Nrmax
    

    #####################################
    projw=[]
    for jatom in atm_l_case.keys():
        Rmt2 = atm_l_case[jatom][1]
        for lc in atm_l_case[jatom][0]:
            Nr0 = struct.npt[jatom]  # Number of radial points from struct file
            # creating larger mesh to Rmt2, which can be larger then Rmt
            r0 = struct.r0[jatom]
            Rmt = struct.rmt[jatom]
            
            if Rmt2<Rmt: Rmt2 = Rmt 
            
            
            Nr = 1 + int( (Nr0-1)*log(Rmt2/r0)/log(Rmt/r0) + 0.99)  # Number of radial points after extending mesh to Rmt2
            
            dx = log(Rmt/r0)/(Nr0-1)
            Rx = r0 * exp( arange(Nr)*dx )                        # First radial point
            
            Vr = hstack( (Vr_all[:,jatom], ones(Nr-Nr0)*Vr_all[-1,jatom] ) )             # Extending it in the interstitial with constant

            if Emu=='EF':
                Elc=EF
            else:
                El = E_all[:,jatom]
                print 'El0=', El
                lapw = ones(len(El), dtype=bool)
                for l in range(len(El)):
                    if El[l]>150:
                        El[l]-=200
                        lapw[l]=False
                print 'linearization E=', El
                Elc = El[lc]
                
            print 'dx=', dx
            print 'Rx[:]=', Rx[0], Rx[Nr0-1], Rx[Nr-1]
            print 'using E=', Elc
                
            A,B,Ae,Be,Aee,Bee,Pei = atpar(rel,lc,Vr,Rx,Elc,struct.r0[jatom],dx,struct.Znuc[jatom])
            
            
            if localize==1:

                if Nr>Nr0:
                    # Below we extend solution beyond Rmt by making value and derivative continuous and make it vanish ar Rmt2
                    cm = SolveForContinuousFunction(A,Ae,Aee,Rx,Nr0,Nr)
                    Ag = hstack(( A[:Nr0] ,  cm[0]*A[Nr0:]+cm[1]*Ae[Nr0:]+cm[2]*Aee[Nr0:] ))
                    cm = SolveForContinuousFunction(B,Be,Bee,Rx,Nr0,Nr)
                    Bg = hstack(( B[:Nr0] ,  cm[0]*B[Nr0:]+cm[1]*Be[Nr0:]+cm[2]*Bee[Nr0:] ))
                
                    # We then normalize the vawe function such that it is normalized within Rmt
                    #overlap = rint13(rel,Ag[:Nr0],Bg[:Nr0],Ag[:Nr0],Bg[:Nr0],dx,struct.r0[jatom])
                    overlap = rint13(rel,Ag[:Nr],Bg[:Nr],Ag[:Nr],Bg[:Nr],dx,struct.r0[jatom])
                    Ag *= 1/sqrt(overlap)  # normalization inside Rmt, and not Rmt2
                    Bg *= 1/sqrt(overlap)
                    A *= 1/sqrt(overlap)
                    B *= 1/sqrt(overlap)
                    print 'overlap=', overlap
                    
                else:
                    Ag = A
                    Bg = B
            elif localize==2:

                cm = SolveForLocalizedFunction(A,Ae,Aee,Rx,Nr)
                Ag = A + cm[0]*Ae + cm[1]*Aee
                cm = SolveForLocalizedFunction(B,Be,Bee,Rx,Nr)
                Bg = B + cm[0]*Be + cm[1]*Bee
                overlap = rint13(rel,Ag[:Nr],Bg[:Nr],Ag[:Nr],Bg[:Nr],dx,struct.r0[jatom])
                Ag *= 1/sqrt(overlap)  # normalization inside Rmt, and not Rmt2
                Bg *= 1/sqrt(overlap)
                print 'cm=', cm
                print 'overlap=', overlap
                
            elif localize==3:

                cm = SolveForLocalizedFunction3(A,Ae,Nr)
                Ag = A + cm*Ae
                cm = SolveForLocalizedFunction3(B,Be,Nr)
                Bg = B + cm*Be
                overlap = rint13(rel,Ag[:Nr],Bg[:Nr],Ag[:Nr],Bg[:Nr],dx,struct.r0[jatom])
                Ag *= 1/sqrt(overlap)  # normalization inside Rmt, and not Rmt2
                Bg *= 1/sqrt(overlap)
                print 'cm=', cm
                print 'overlap=', overlap
            else:
                print 'Dot yet implemented!'
                sys.exit(1)
            ############################################################################

            print >> fout, '#', Nr, Nr0, jatom+1, lc
            for ir in range(Nr):
                print >> fout, Rx[ir], Ag[ir], Bg[ir]
                
            projw.append( (Rx,Ag,Bg,A,B) )
            
    fout.close()
    return projw


if __name__ == '__main__':
    usage = """
    usage: %prog [ options ]

    Saved the projector wave function into projectorw.dat file.
    There are two Rmt used in this code. The first Rmt1 is muffin-thin sphere
    from structure file. The second is Rmt2 and is written in case.indmfl file.
    """
    parser = optparse.OptionParser(usage)
    parser.add_option("-l", "--localize",  dest="localize",    type="int", default=1, help="How localized should be the function")
    # Next, parse the arguments
    (options, args) = parser.parse_args()

    w2k = utils.W2kEnvironment()

    atms, Lsa, icpsa, Rm2 = ReadIndmfl(w2k.case)

    choice = raw_input("""Should linearization energy for the projector be\n    [1]: EF  or\n    [2]: DFT linearization energy.\n  Enter [1|2]: """)
    if choice.strip()=='1':
        Emu='EF'
    else:
        Emu='Emu'
        
    projw = main(w2k.case, atms, Lsa, icpsa, Rm2, options.localize, Emu)

    for Rx,Ag,Bg,A,B in projw:
        subplot(2,1,1)
        plot(Rx, A, 'r-')
        plot(Rx, Ag, 'y-')
        subplot(2,1,2)
        plot(Rx, B, 'r-')
        plot(Rx, Bg, 'y-')
        show()
        
        ### self.atoms.keys()  -> atms
        ### self.Lsa  -> Lsa
        ### self.icpsa -> icpsa
        ### [self.atoms[iatom][3] for iatom in self.atoms.keys()]  -> Rmt2
        

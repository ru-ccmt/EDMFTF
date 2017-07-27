#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from scipy import *
from scipy import interpolate
from scipy import integrate
from scipy import special
from scipy import optimize
#from pylab import *
import optparse
from CmpSlaterInt import CmpSlaterInt
from CmpSlaterInt import GetWaveF
from CmpSlaterInt import  SlaterF2J

Ry2eV = 13.60569193

def GetCoulombU(lmbda,Rx,Ag,Bg,l,Nk,U0):
    Fk = CmpSlaterInt(lmbda,Rx,Ag,Bg,l,Nk,False)
    UJs=SlaterF2J(Fk,l)
    return UJs[0]-U0

def GetCoulombUJ(lmbda_eps,Rx,Ag,Bg,l,Nk,U0,J0):
    (lmbda,eps) = lmbda_eps
    Fk = CmpSlaterInt(lmbda,Rx,Ag,Bg,l,Nk,False)
    UJs=SlaterF2J(array(Fk)/eps,l)
    UH = UJs[0]
    JH = sum(UJs[1:])/len(UJs[1:])
    return array([UH-U0 , JH-J0])

def GetLambdas(UJ_icase, projector_file='projectorw.dat', options_Nk=8):
    Nk=2**options_Nk+1

    #print 'UJ_icase=', UJ_icase
    
    (Rx_, Ag_, Bg_, ls_) = GetWaveF(projector_file)
    UlamJ={}
    for icase in UJ_icase.keys():
        Rx = Rx_[icase]
        Ag = Ag_[icase]
        Bg = Bg_[icase]
        l = ls_[icase]
        #lmbda = optimize.brentq(GetCoulombU, 0., 3., args=(Rx,Ag,Bg,l,Nk,U_icase[icase]))
        lmbda = 0.000001
        epsilon=1.
        if UJ_icase[icase][0]>0:
            Uc = UJ_icase[icase][0]
            lmbda = optimize.brentq(GetCoulombU, 0., 3., args=(Rx,Ag,Bg,l,Nk,Uc))
            epsilon=1.
            if UJ_icase[icase][1]>0:
                Jc = UJ_icase[icase][1]
                sol = optimize.root(GetCoulombUJ, [lmbda,epsilon], args=(Rx,Ag,Bg,l,Nk,Uc,Jc))
                if sol.success:
                    (lmbda,epsilon) = sol.x
                else:
                    print 'ERROR:', sol.message
            #print 'lmbda=', lmbda, 'epsilon=', epsilon, 'for Uc=', Uc, 'and Jc=', Jc
        if len(UJ_icase[icase])<=3 :
            #print 'lambda=', lmbda, 'epsilon=', epsilon
            Fk = CmpSlaterInt(lmbda,Rx,Ag,Bg,l,Nk,True)
            Fk = array(Fk)/epsilon
            UJs=SlaterF2J(Fk,l)
            print 'Fk=', Fk.tolist()
            print 'U,J=', UJs
            UlamJ[icase]=(UJ_icase[icase][0],lmbda,epsilon,UJs[1:])
        else:  # this is very uncommon case, where we renormalize lambda by input parameter
            nfraction = UJ_icase[icase][3]
            lmbda *= nfraction
            print 'renormalizing lmbda by', UJ_icase[icase][3],' to', lmbda
            Fk = CmpSlaterInt(lmbda,Rx,Ag,Bg,l,Nk,True)
            epsilon = Fk[0]/Uc
            print 'lambda=', lmbda, 'epsilon=', epsilon
            Fk = array(Fk)/epsilon
            print 'Fk=', Fk.tolist()
            if UJ_icase[icase][1]>0: # J existed before, hence use it
                Js = ones(l)*UJ_icase[icase][1] 
            print 'J=', Js
            UlamJ[icase]=(UJ_icase[icase][0],lmbda,epsilon,Js.tolist())
        
        #print 'lambda=', lmbda, 'epsilon=', epsilon
        #Fk = CmpSlaterInt(lmbda,Rx,Ag,Bg,l,Nk,True)
        #Fk = array(Fk)/epsilon
        #UJs=SlaterF2J(Fk,l)
        #print 'Fk=', Fk.tolist()
        #print 'U,J=', UJs
        #UlamJ[icase]=(UJ_icase[icase][0],lmbda,epsilon,UJs[1:])
    return UlamJ

def GetDielectricFunctions(UJ_icase, projector_file='projectorw.dat', options_Nk=8):
    Nk=2**options_Nk+1

    (Rx_, Ag_, Bg_, ls_) = GetWaveF(projector_file)
    UlamJ={}
    for icase in UJ_icase.keys():
        Rx = Rx_[icase]
        Ag = Ag_[icase]
        Bg = Bg_[icase]
        l = ls_[icase]
        epsilon=1.
        if UJ_icase[icase][0]>0:
            lmbda=1e-6
            if len(UJ_icase[icase])>2: lmbda=UJ_icase[icase][2]
            (Uunscr,Junscr) = GetCoulombUJ([lmbda,epsilon],Rx,Ag,Bg,l,Nk,0.0,0.0)  # First compute the unscreened interaction
            epsilon = Uunscr/UJ_icase[icase][0]                          # dielectric constant is the ration between the unscreened and screened interaction
        
        print 'lambda=', lmbda, 'epsilon=', epsilon
        Fk = CmpSlaterInt(lmbda,Rx,Ag,Bg,l,Nk,True)
        Fk = array(Fk)/epsilon
        UJs=SlaterF2J(Fk,l)
        print 'Fk=', Fk.tolist()
        print 'U,J=', UJs
        UlamJ[icase]=(UJ_icase[icase][0],lmbda,epsilon,UJs[1:])
    return UlamJ

if __name__ == '__main__':
    usage = """usage: %prog [ options ]
    Help for the command
    """
    import sys
    parser = optparse.OptionParser(usage)
    parser.add_option("-n", "--npts",  dest="Nk",    type="int", default=8, help="Number of points in the radial integration will be 2**NK+1.")
    parser.add_option("-i", "--inw", dest="finput",  default='projectorw.dat', help="filename of the input file containing projector", metavar="FILE")
    parser.add_option("-l", "--lambda", dest="lmbda",  type="float", default=0.0,  help="The screening parameter lambda (in units of inverse bohr radius)")
    parser.add_option("-U", "--U", dest="U0", type="float", default=0.0, help="Set the Hubbard U, if you want to know lambda, which gives you such U")
    parser.add_option("-J", "--J", dest="J0", type="float", default=0.0, help="Set Hund's J, if you want to know (epsilon,lambda) from V_C=exp(-lambda*r)/(r*epsilon), which gives you (U,J)")
    # Next, parse the arguments
    (options, args) = parser.parse_args()
    
    #UlamJ=GetLambdas({0:[options.U0,options.J0,options.DC_lmbda]}, projector_file=options.finput)
    #print 'UlamJ=', UlamJ
    #sys.exit(0)
    
    #UlamJ = GetDielectricFunctions({0:[options.U0,options.J0]}, projector_file=options.finput)
    #print 'UlamJ=', UlamJ
    #sys.exit(0)

    
    Nk=2**options.Nk+1

    (Rx_, Ag_, Bg_, ls_) = GetWaveF(options.finput)
    for icase in range(len(Rx_)):
        Rx = Rx_[icase]
        Ag = Ag_[icase]
        Bg = Bg_[icase]
        l = ls_[icase]

        if options.U0>0:
            lmbda = optimize.brentq(GetCoulombU, 0., 3., args=(Rx,Ag,Bg,l,Nk,options.U0))
            epsilon=1.
            if options.J0>0:
                sol = optimize.root(GetCoulombUJ, [lmbda,epsilon], args=(Rx,Ag,Bg,l,Nk,options.U0,options.J0))
                if sol.success:
                    (lmbda,epsilon) = sol.x
                else:
                    print 'ERROR:', sol.message
                 
        else:
            lmbda = options.lmbda
            epsilon=1.
            
        print 'lambda=', lmbda, 'epsilon=', epsilon
        Fk = CmpSlaterInt(lmbda,Rx,Ag,Bg,l,Nk,True)
        UJs=SlaterF2J(Fk,l)
        print 'Fk=', (array(Fk)/epsilon).tolist()
        print 'U,J=', (array(UJs)/epsilon).tolist()
        
        

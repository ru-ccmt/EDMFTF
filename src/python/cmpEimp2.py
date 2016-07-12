#!/usr/bin/env python
from scipy import *
import sys,os,re
#from pylab import *
from scipy import interpolate, special

def ReadPARAMS2(fparams, fout):
    fp = open(fparams, 'r')
    cix='imp.cix'
    Delta='Delta.inp'
    Gf='Gf.out'
    Sig='Sig.out'
    for line in fp:
        n = line.find('#')
        if n>=0:
            line = line[:n]
        dat = line.split()
        if len(dat)>=2:
            m = re.search('cix|Delta|Gf|Sig', dat[0])
            if m is not None :
                exec(dat[0]+'="'+dat[1]+'"')
                print >> fout, 'executed:', dat[0]+'="'+dat[1]+'"'
            m = re.search('beta|mu',dat[0])
            if m is not None :
                exec(dat[0]+'='+dat[1])
    return (cix,Delta,Gf,Sig,beta,mu)

def ReadCix(cix):
    fcx = open(cix,'r')
    fcx.next();fcx.next()
    (Nclust,Nm,Ns,maxm) = map(int,fcx.next().split())
    fcx.next()
    Deg={}
    for i in range(Ns):
        dat = map(int,fcx.next().split())
        dim = dat[1]
        if (dim>1): 
            print 'WARNING: The impurity functional cmpEimp2.py works for now only when baths are diagonal. You should not trust total free energy. Needs some development...!'
        for j in range(dim):
            for k in range(dim):
                if k==j:
                    ii = dat[2+k+j*dim]
                    if Deg.has_key(ii):
                        Deg[ii] += 1
                    else:
                        Deg[ii]=1
    print 'Deg=', Deg
    fcx.next()
    ss = fcx.next()
    if ss[:11]=='FL_FROM_IFL':
        for i in range(Ns):
            fcx.next()
        fcx.next()
        ss = fcx.next()
    Eks = map(float,ss.split())
    print 'Eks=', Eks
    return (Deg, array(Eks))


def ReadGandSigandDelta(Gf, Sig, Delta, fout):
    Delta_data = loadtxt(Delta).transpose()
    omD = Delta_data[0]
    Dlta = Delta_data[1::2]+Delta_data[2::2]*1j
    
    Gf_data = loadtxt(Gf).transpose()
    om = Gf_data[0]
    Gfc = Gf_data[1::2]+Gf_data[2::2]*1j

    Sig_data = loadtxt(Sig).transpose()
    om3 = Sig_data[0]
    Sigc = Sig_data[1::2]+Sig_data[2::2]*1j
    
    print 'shapes=', shape(om), shape(omD), shape(om3)
    #print 'Diff=', sum(abs(om-omD))
    #if sum(abs(om-omD))>1. or sum(abs(om-om3))>1.:
    #    print >> fout, 'frequency mesh of ',Gf,' and ',Sig, 'and', Delta,' are not equal', om-om2, om-om3
    
    fg = open(Sig, 'r')
    first_line = fg.next()
    fg.close()
    dat = first_line.split()
    for d in dat:
        m = re.search('nf|mu|T|TrSigmaG|Ekin|Epot|mom', d)
        if m is not None: exec(d)
            
    print >> fout, 'nf=', nf, 'mu=', mu, 'Ekin=', Ekin, 'Epot=', Epot
    return (om, omD, Gfc, Sigc, Dlta, nf, mu, Ekin, Epot, TrSigmaG, mom) 


def ferm(x):
    if x>300:
        return 0
    elif x<-300:
        return 1
    else:
        return 1/(exp(x)+1)

def GetHighFrequency(CC,om):
    " Approximates CC ~  A/(i*om-C) "
    A = 1./imag(1/(CC[-1]*om[-1]))
    C = -A*real(1./CC[-1])
    return (A, C)

def ReturnHighFrequency(A,C,beta,E):
    " Returns the value of Tr(A/(iom-C) 1/(iom-E) )"
    return A*(ferm(E*beta)-ferm(C*beta))/(E-C)

def LogGimp(omega,Gfc,EimpS,beta,Deg,fout):
    " Tr(log(-G_imp))"
    def NonIntF0(beta,Ene):
        if beta*Ene>200: return 0.0
        if beta*Ene<-200: return Ene
        return -log(1+exp(-Ene*beta))/beta
    
    lnGimp=0.
    for i in range(len(Deg)): 
        print >> fout, 'EimpS=', EimpS[i], ', real(iom-1/G)', real(-1/Gfc[i][-1])
        if EimpS[i] < 1000 :
            Gd = Gfc[i,:]
            eimps = EimpS[i]
            
            CC = omega*1j-eimps-1/Gd
            A,C = GetHighFrequency(CC,omega)

            ff = log(-Gd)-log(-1./(omega*1j-eimps)) - 1./(omega*1j-eimps)*A/(omega*1j-C)
            lnGimp_i = 2*sum(ff.real)/beta+NonIntF0(beta,eimps)+ReturnHighFrequency(A,C,beta,eimps)
            lnGimp_i *= Deg[i]
            
        else:
            sum_imp = 0.0
            Fnint_imp = 0.0
            lnGimp_i = 0.0
        lnGimp += lnGimp_i
        print >> fout, 'lnGimp: Deg[i]=', Deg[i], 'Tr(log(-Gimp))=', lnGimp_i/Deg[i], 'Tr(log(-Gimp))*deg=', lnGimp_i
    return lnGimp


def LogGimpSig(omega, Gfc, Sigc, Vdc, EimpS, beta, Deg, fout):
    " Tr(log(1+(Sig-Vdc)G))"
    def NonIntF0(beta,Ene):
        if beta*Ene>200: return 0.0
        if beta*Ene<-200: return Ene
        return -log(1+exp(-Ene*beta))/beta
    
    lnGSimp=0.
    for i in range(len(Deg)):#len(Gfc)):
        Gd = Gfc[i,:]
        Sg = Sigc[i,:]
        eimps = EimpS[i]
        s_oo = Sigc[i,-1].real
        vdc = Vdc[i]
        
        #print 'shape(Gd)=', shape(Gd)
        ff0 = log((Sg-vdc)*Gd+1.)
        ff1 = (s_oo-vdc)/(omega*1j-eimps)
        
        lnGimp_i = 2*sum((ff0-ff1).real)/beta + (s_oo-vdc)*ferm(eimps*beta)
        lnGimp_i *= Deg[i]

        lnGSimp += lnGimp_i
        print >> fout, 'lnGSimp: Deg[i]=', Deg[i], 'Tr(log(1+Sig*Gimp))=', lnGimp_i/Deg[i], 'Tr(log(1+Sig*Gimp))*deg=', lnGimp_i

        #plot(omega, real(ff0-ff1), label='diff-real')
        #plot(omega, imag(ff0-ff1), label='diff-imag')
        #legend(loc='best')
        #show()
    return lnGSimp
    
def fTrDeltaG(om, omD, Dlta, Gfc, Deg, EimpS, beta, fout):
    def GetInterpolation(omD, Delta):
        "Returns derivative  d Delta/d om_i"
        Dr = interpolate.UnivariateSpline(omD, Delta.real, s=0,k=3)
        Di = interpolate.UnivariateSpline(omD, Delta.imag, s=0,k=3)
        return (Dr, Di)
    lngh = min(len(Dlta),len(Gfc),len(Deg))
    GDf=0.0
    for i in range(lngh):
        if len(omD)!=len(om):
            # interpolating Delta
            Dr, Di = GetInterpolation(omD, Dlta[i])
            Delta = Dr(om)+Di(om)*1j
        else:
            Delta = Dlta[i]
        (A,C) = GetHighFrequency(Delta,om)
        print >> fout, 'A_d=', A, 'C_d=', C
        eimps = EimpS[i]
        GD1 = Gfc[i]*Delta - 1./(om*1j-eimps)*A/(om*1j-C)
        dGD = Deg[i]*(2*sum(GD1.real)/beta + ReturnHighFrequency(A,C,beta,eimps))
        GDf += dGD
        #print >> fout, 'dTr(Delta*G)=', dGD/Deg[i]
    return GDf

def fTrSigmaG(om, Gfc, Sigc, Deg, EimpS, beta, fout):
    SDf=0.0
    for i in range(len(Deg)):
        if EimpS[i] < 1000 :
            Gd = Gfc[i,:]
            eimps = EimpS[i]
            Sg = Sigc[i,:]
            s_oo = Sigc[i,-1].real
            Sg0 = Sg[:] - s_oo   # Sigma-s_oo
            (A,C) = GetHighFrequency(Sg0,om)
            #print '::: A,C, for Sigma=', A, C
            C=1.
            ff = Gd[:]*Sg0[:] - 1./(om*1j-eimps)*A/(om*1j-C)
            SDf_i = 2*sum(ff.real)/beta + ReturnHighFrequency(A,C,beta,eimps)
            SDf_i *= Deg[i]
        else:
            SDf_i = 0.0
        SDf += SDf_i
        #print >> fout, 'Tr((Sigma-s_oo)*G): Deg[i]=', Deg[i], 'Tr((Sigma-s_oo)*G)=', SDf_i/Deg[i], 'Tr((Sigma-s_oo)*G)*deg=', SDf_i
    return SDf
            
def CmpImpurityCorrection2(om, omD, Dlt, Gfc, Deg, EimpS, beta, fout):
    def GetInterpolation(omD, Delta):
        "Returns derivative  d Delta/d om_i"
        Dr = interpolate.UnivariateSpline(omD, Delta.real, s=0,k=3)
        Di = interpolate.UnivariateSpline(omD, Delta.imag, s=0,k=3)
        return (Dr, Di)

    Delta_inf=[]
    for b in range(len(Dlt)):
        Delta_inf.append( [Dlt[b,-1].real*omD[-1]**2, Dlt[b,-1].imag*omD[-1] ] )
    
    cxx = 0.0
    lngh = min(len(Dlt),len(Gfc),len(Deg))
    
    oms = om
    Gfs = Gfc
    if oms[0]<omD[0] and abs(oms[0]-omD[0])<1e-6: oms[0]=omD[0]
    
    print >> fout, 'oms[-1]=', oms[-1], 'omD[-1]=', omD[-1]
    
    for i in range(lngh):
        Dr, Di = GetInterpolation(omD, Dlt[i])
        dDlt=zeros(len(oms),dtype=complex)
        for ix,x in enumerate(oms):
            if x<=omD[-1]:
                dDlt[ix] = Dr.derivatives(x)[1]+Di.derivatives(x)[1]*1j
            else:
                dDlt[ix] = -2*Delta_inf[i][0]/x**3 - Delta_inf[i][1]/x**2 * 1j
        
        CC = oms*dDlt
        #savetxt( 'brisx', vstack((oms, real(CC), imag(CC))).transpose() )
        #print 'shape(dDlt)=', shape(dDlt), 'shape(CC)=', shape(CC), 'shape(Gfs)=', shape(Gfs)
        (A,C) = GetHighFrequency(CC,oms)
        print >> fout, 'A=', A, 'C=', C
        eimps = EimpS[i]
        GD1 = Gfs[i]*CC - 1./(oms*1j-eimps)*A/(oms*1j-C)
        dcxx = Deg[i]*(2*sum(GD1.real)/beta + ReturnHighFrequency(A,C,beta,eimps))
        cxx += dcxx
        print >> fout, 'dcxx=', dcxx/Deg[i]
        
    print >> fout,'cxx=', cxx
    return cxx

def GiveImpFunctional(dire, fparams, Vdc, fout=sys.stdout):

    if not (os.path.exists(dire+'/'+fparams) and os.path.getsize(dire+'/'+fparams)>0):
        print >> fout, 'Give input PARAMS file'
        sys.exit(0)
    
    (cix,Delta,Gf,Sig,beta,mu_qmc) = ReadPARAMS2(dire+'/'+fparams, fout)
    (Deg, Eks) = ReadCix(dire+'/'+cix)
    
    
    (om, omD, Gfc, Sigc, Dlta, nf, mu, Ekin, Epot, TrSigmaG, mom) = ReadGandSigandDelta(dire+'/'+Gf, dire+'/'+Sig, dire+'/'+Delta, fout)
    s_oo = Sigc[:,-1].real
    print 's_oo=', s_oo
    
    lengh = min(len(Eks),len(s_oo))
    EimpS = Eks[:lengh]-mu_qmc+s_oo[:lengh]
    # Tr(log(-Gimp))
    lnGimp = LogGimp(om,Gfc,EimpS,beta,Deg,fout) 

    #lnGS = LogGimpSig(om, Gfc, Sigc, Vdc, EimpS, beta, Deg, fout)

    trsigg =  fTrSigmaG(om, Gfc, Sigc, Deg, EimpS, beta, fout) # numericall Tr((Sigma-s_oo)*G)
    trsigg_oo = sum(mom[:]*s_oo[:])

    # Tr(Delta*G)
    TrDeltaG = fTrDeltaG(om, omD, Dlta, Gfc, Deg, EimpS, beta, fout)
    print >> fout, 'TrDeltaG=', TrDeltaG, '<k>/T=', Ekin
    # Tr(om*dDelta/dom)
    cxx = CmpImpurityCorrection2(om, omD, Dlta, Gfc, Deg, EimpS, beta, fout)

    # Eimp = Tr(Delta*G)+1/2*Tr(Sigma*G)+eimp*nimp-Tr(om*dDelta/dom)
    #Eimp = Ekin+Epot-cxx-mu*nf                # Impurity Energy
    Eimp = TrDeltaG+Epot-cxx-mu*nf                # Impurity Energy

    # Potthof2[G] == Phi[G]-1/2*Tr(Sigma*G)+T*S_imp = Eimp-Tr(log(-G_imp))+1/2*Tr(Sigma*G)
    #Potthof2 = Eimp-lnGimp+TrSigmaG            # Phi[G]-1/2*Tr(Sigma*G)+T*S_imp

    #nd_numeric=zeros(len(Deg))
    #for i in range(len(Deg)):
    #    Gd = Gfc[i,:]
    #    eimps = EimpS[i]
    #    ff = Gd[:]-1./(om*1j-eimps)
    #    nd_numeric[i] = (2*sum(ff.real)/beta + ferm(eimps*beta))*Deg[i]
    #print 'nd_numeric=', nd_numeric
    
    #Eimp1 = lnGS - lnGimp + Eimp  # Tr(log(1+Sigma*G))-Tr(log(-Gloc))+E_imp

    logZatom=0
    filename=dire+'/ctqmc.log'
    if os.path.isfile(filename):
        fw=open(filename,'r')
        lines=fw.readlines()
        for line in lines:
            m = re.search('log\(Zatom\)=(\d+\.\d+)',line)
            if m is not None:
                logZatom = float(m.groups()[0])
        fw.close()
        #print 'logZatom=', logZatom
    P0=0
    filename=dire+'/histogram.dat'
    if os.path.isfile(filename):
        fw=open(filename,'r')
        fw.next()# comment
        P0 = float(fw.next().split()[1])
        #print 'P0=', P0

    Eimp2 = TrSigmaG   # 1/2*Tr(Sigma*G)
    print >> fout, 'lnGimp=', lnGimp
    print >> fout, "%-12s "*6 % ('Eimp', 'Ekin', 'Epot', 'cxx', 'mu*nf', 'lnGimp')
    print >> fout, "%-12.6f "*6 % (Eimp, TrDeltaG, Epot, cxx, mu*nf, lnGimp)
    print >> fout, '#',"%-12.6f "*6 % (Eimp+Ekin-TrDeltaG, Ekin, Epot, cxx, mu*nf, lnGimp)
    print >> fout, 'Eimp=', Eimp
    if logZatom and P0>1e-7:
        logZ = logZatom-log(P0)
        Fimp = -logZ/beta
        Entropy = (Eimp-Fimp)*beta
        WARN=''
        if (P0<1e-5): WARN='   ** WARNING P0 is very small, hence Entropy is not accurate'
        print >> fout, 'Entropy=', Entropy, 'Fimp=', Fimp, 'logZatom=', logZatom, 'P0=', P0, WARN
    #print >> fout, 'P[Sigma]=', Potthof2

    Vdc_nd = sum(array(mom[:])*array(Vdc[:]))
    print >> fout, 'Vdc*nf=', Vdc_nd
    Phi_DMFT = -lnGimp+Eimp+(trsigg+trsigg_oo)
    print >> fout, 'Phi_DMFT=', Phi_DMFT
    print >> fout, 'Phi_DMFT - Vdc*nf=', Phi_DMFT-Vdc_nd
    Ry2eV = 13.60569193
    print >> fout, 'Phi_DMFT - Vdc*nf=', (Phi_DMFT-Vdc_nd)/Ry2eV, 'Ry'
    print >> fout, 'Epot - Vdc*nf=', (TrSigmaG-Vdc_nd)/Ry2eV, 'Ry'
    return (Phi_DMFT, lnGimp, Eimp, Vdc_nd)


if __name__ == '__main__':
    import re
    
    if len(sys.argv)<2:
        dire = '.'
    else:
        dire = sys.argv[1]
    Ry2eV = 13.60569193

    fi = open(dire+'/Eimp.inp','r')
    line = fi.next()  # Eimp
    line = fi.next()  # Overlap
    line = fi.next()  # Vdc
    Vdc = map(float,re.sub(r'#.*',' ',line).split())

    (Phi_DMFT, lnGimp, Fimp, Vdc_nd) = GiveImpFunctional(dire, 'PARAMS', Vdc)
    print 'ev: Phi_DMFT-Vdc*nd=', Phi_DMFT-Vdc_nd, 'Fimp+TS=', Fimp, 'logGimp=', lnGimp
    print 'Ry: Phi_DMFT-Vdc*nd=', (Phi_DMFT-Vdc_nd)/Ry2eV, 'FImp+TS=', Fimp/Ry2eV, 'logGimp=', lnGimp/Ry2eV
    


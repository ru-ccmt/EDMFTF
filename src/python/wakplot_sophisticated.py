#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule

from scipy import *
from pylab import *
import glob, os
import numpy as np
import scipy

from distutils.version import StrictVersion
if StrictVersion(scipy.__version__) > StrictVersion('0.19.0'):
    import weave
else:
    import scipy.weave as weave

code="""
     #line 11 "wakplot.py"
     using namespace std;

     double w=omega;
     double Ax=0;
     for (int ib=0; ib<nbands; ib++){
        complex<double> ew=ekw(ib);
        if (ew.imag() > -small) ew=complex<double>(ew.real(),-small);
        Ax += -(1./(w-ew)).imag()/M_PI;
     }
     return_val = Ax;
     """

code2="""
     #line 25 "wakplot.py"
     using namespace std;

     double w=omega;
     int norb=Am.size();

     for (int ib=0; ib<nbands; ib++){
        complex<double> ew=ekw(ib);
        if (ew.imag() > -small) ew=complex<double>(ew.real(),-small);
        for (int iorb=0; iorb<norb; iorb++){
           double coh = cohi(iorb,ib);
           double Ac = -abs(coh)*(1./(w-ew)).imag()/M_PI;
           Am(iorb) += Ac;
        }
     }
     """


def ReadKlist():
    wkpointi=[]
    wkpoints=[]
    wg = glob.glob('*.klist_band')
    if len(wg)>0:
        fg = open(wg[0], 'r')
        for il,line in enumerate(fg):
            if line[:3]=='END': break
            com = line[:10].split()
            if com:
                wkpoints.append(com[0])
                wkpointi.append(il)
        print wkpointi
        print wkpoints

    nkp = il #wkpointi[-1]+1
    print 'nkp=', nkp
    return (wkpoints,wkpointi,nkp,wg[0])

if __name__ == '__main__':
    aspect_scale=0.3
    small = 1e-5
    DY = 0.0 #
    intensity = [3.,3.,3.,1.]
    COHERENCE=True
    orb_plot=[1,1,2,2,1,1,0,0,0,0,0,0,0,0] #2,2,2,2,2,2]  # should have integers 0,1,2 for red,green,blue
    #orb_plot=[]  # plotting unfolded band structure, but all orbits with hot map

    formt=None
    if len(sys.argv)<2:
        print 'Give filename for eigenvalues'
        sys.exit(0)
    fname = sys.argv[1]

    if COHERENCE and os.path.isfile('UL.dat_') and os.path.getsize('UL.dat_')>0 :
        COHERENCE=True
        print 'Using coherence factors'
    else:
        COHERENCE=False
        print 'Not using coherence factors'
    
    fEF = open('EF.dat', 'r')
    mu = float(fEF.next().split()[0])
    print 'mu=', mu
    
    (wkpoints,wkpointi,nkp,fklist) = ReadKlist()

    if COHERENCE:
        if orb_plot: norb_plot=max(orb_plot)+1
        else: norb_plot=0
        
        fL = open('UL.dat_', 'r')
        fR = open('UR.dat_', 'r')
        dat = fL.next().split()
        fR.next()
        nkp2,nsymop,nom,norbitals = map(int,dat[1:5])  # dimensions
        dims = map(int,dat[5:5+norbitals])
        fL.next()
        fR.next()
        
        if nkp!=nkp2:
            print 'ERROR: UL.dat_ does not have the same number of k-points as '+fname
            sys.exit(1)
        if nsymop>1:
            print 'ERROR: More than one group operation. We do not want to symmetrize over group operations when plotting spectra!'
            sys.exit(1)
    else:
        norb_plot=0
        
    fdat = open(fname, 'r')
    strn = fdat.next()
    if strn[0]=='#':
        formt='dmft1'
    elif strn[0]=='!':
        formt='dmftgk'
    else:
        print 'Can not recognize format of', fname
        sys.exit(0)
    data = strn[1:].split()
    nomega = int(data[4])

    ikp=0
    if not COHERENCE or norb_plot==0:
        Akom=zeros((nomega,nkp),dtype=float)
    else:
        Akom=zeros((nomega,nkp,4),dtype=float)
    try:
        while True:
            if ikp>0: strn = fdat.next()
            data = strn[1:].split()
            if formt=='dmft1':
                (ikp, isym, nbands, nemin, nomega) = map(int, data[:5])
            else:
                (ikp, isym, tot_nbands, tot_nemin, nomega) = map(int, data[:5])
            wgh = float(data[5])
            
            om=[]
            for iom in range(nomega):
                dat = fdat.next().split()
                omega = float(dat[0])
                om.append( omega )
                if formt=='dmft1':
                    data = array(map(float, dat[1:]))
                else:
                    nbands = int(dat[1])
                    nemin = int(dat[2])
                    data = array(map(float, dat[3:]))
                ekw = data[::2]+data[1::2]*1j - mu

                if len(ekw)!=nbands: print 'ERROR nbands=', nbands, 'len(ekw)=', len(ekw), 'ekw=', ekw.tolist()

                if not COHERENCE:
                    cohi=ones((norb_plot+1,nbands),dtype=float)
                    Akom[iom,ikp-1] = weave.inline(code, ['nbands', 'omega', 'ekw', 'small'],
                                                   type_converters=weave.converters.blitz, compiler = 'gcc')
                else:
                    line1 = fL.next()
                    fR.next()
                    dat = line1.split()
                    omega2 = float(dat[0])
                    if abs(omega2-omega)>1e-5:
                        print 'WARNING: Frequency in UL.dat_ '+str(omega2)+' and in '+fdat+' '+str(omega)+' are not equal!'
                    nbands2 = int(dat[1])
                    if nbands!=nbands2:
                        print 'ERROR: Number of bands in UL.dat_ and '+fname+' is not consistent'
                        
                    cohi = zeros((norb_plot+1,nbands2),dtype=float)
                    for ibnd in range(nbands2):
                        datL = fL.next().split()
                        datR = fR.next().split()
                        ii = int(datL[0])
                        
                        gl = array(map(float,datL[1:]))
                        gr = array(map(float,datR[1:]))
                        UL = gl[::2]+gl[1::2]*1j
                        UR = gr[::2]+gr[1::2]*1j
                        coh = abs(UL*UR)

                        for iorb in range(len(coh)):
                            if iorb<len(orb_plot) and orb_plot[iorb]>=0:
                                cohi[orb_plot[iorb],ibnd] += coh[iorb]
                        #cohi[norb_plot,ibnd] = sum(coh)
                        cohi[norb_plot,:] = ones(nbands2)

                    #cohi = zeros((4,nbands2),dtype=float)
                    #cohi[2,:] = abs(ones(nbands)-cohi_[2,:])
                    #cohi[:2,:] = cohi_[:2,:]
                    #cohi[3,:] = ones(nbands2)
                    #print 'cohi=', cohi
                    #Am = zeros(4,dtype=float)
                    Am = zeros(norb_plot+1,dtype=float)
                    weave.inline(code2, ['Am', 'nbands', 'omega', 'ekw', 'small','cohi'],
                                 type_converters=weave.converters.blitz, compiler = 'gcc')

                    if norb_plot>0:
                        np=min(norb_plot,3)
                        Akom[iom,ikp-1,:np] = Am[:np]
                        Akom[iom,ikp-1,3] = Am[norb_plot]
                        #Akom[iom,ikp-1,:] = Am[:]
                    else:
                        Akom[iom,ikp-1]=Am[0]
                
            print 'k-point', ikp, 'finished'
            
    except StopIteration:
        pass
    if nkp!=ikp:
        print 'Number of k-points from '+fklist+' and '+fname+' are different'
        
    om=array(om)
    print 'shape(Akom)=', shape(Akom)

    #np.save('brisi', Akom)
    #Akom = np.load('brisi.npy')
    
    
    (ymin,ymax) = (om[0]+DY,om[-1]+DY)
    (xmin,xmax) = (0, shape(Akom)[1]-1)
    
    print 'xmin,xmax,ymin,ymax=', xmin, xmax, ymin, ymax

    if not COHERENCE or norb_plot==0:
        vmm = [0,max(map(max,Akom))*intensity[0]]
        imshow(Akom, interpolation='bilinear', cmap=cm.hot, origin='lower', vmin=vmm[0], vmax=vmm[1], extent=[xmin,xmax,ymin,ymax], aspect=aspect_scale*(xmax-xmin)*0.8/(ymax-ymin) )
        icolr='w'
    else:
        Akw = zeros(shape(Akom))
        for iom in range(len(Akom)):
            for ikp in range(len(Akom[iom])):
                Akw[iom,ikp,:] = 1-exp(-Akom[iom,ikp,:]*array(intensity))
        
        imshow(Akw, interpolation='bilinear', origin='lower', extent=[xmin,xmax,ymin,ymax], aspect=aspect_scale*(xmax-xmin)*0.8/(ymax-ymin) )
        icolr='k'
        
        
    for i in range(len(wkpointi)):
        print 'wp=', wkpointi[i]
        plot([wkpointi[i],wkpointi[i]], [ymin,ymax], icolr+'-')
        
    plot([xmin,xmax],[0,0], icolr+':')

    dytck=0.01*max(abs(ymin),abs(ymax))
    Ntck=5
    for j in range(len(wkpointi)-1):
        for ix in range(1,Ntck):
            x = wkpointi[j]+(wkpointi[j+1]-wkpointi[j])*ix/float(Ntck)
            plot([x,x],[-dytck,dytck], icolr+'-')
        
    axis([xmin,xmax,ymin,ymax])
    xticks( wkpointi, wkpoints, fontsize='x-large' )
    #colorbar()
    show()
    

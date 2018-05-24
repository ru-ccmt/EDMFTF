#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule

from scipy import *
from pylab import *
import glob, os, sys
import scipy
if int((scipy.__version__).split('.')[1]) > 14:
    import weave
else:
    import scipy.weave as weave

code="""
     #line 10 "wakplot.py"
     using namespace std;

     double Ry2eV = 13.6056920311654;

     double Ax = 0;
     for (int ib=0; ib<nbands; ib++){
        complex<double> ekw=ekom(ib);
        if (ekw.imag() > -small) ekw=complex<double>(ekw.real(),-small);
        complex<double> gc = abs(cohd(ib))/(omega+mu-ekw);
	//complex<double> gc = 1./(omega+mu-ekw);
        Ax += -gc.imag()/M_PI;
     }
     return_val = Ax;
     
     """
if __name__ == '__main__':

    if len(sys.argv)<2:
        intensity = 0.2
    else:
        intensity = float(sys.argv[1])
        
    small = 1e-5 # 0.01 # 1e-5
    #itensity = 0.2
    DY = 0 # 0.01318

    fEF = open('EF.dat', 'r')
    mu = float(fEF.next().split()[0])

    print 'mu=', mu


    wg = glob.glob('*.klist_band')
    if len(wg)>0:
        fg = open(wg[0], 'r')
        wkpointi=[]
        wkpoints=[]
        for il,line in enumerate(fg):
            if line[:3]=='END': break
            com = line[:10].split()
            if com:
                legnd=line.split()[0]
                wkpoints.append(legnd)
                wkpointi.append(il)
        print wkpointi
        print wkpoints

    nkp = wkpointi[-1]+1
    print 'nkp=', nkp
    fdat = open('eigvals.dat', 'r')
    if os.path.isfile('cohfactorsd.dat'):
        fcoh = open('cohfactorsd.dat', 'r')
    else:
        fcoh = None
    
    ikp=0
    Akom=[]
    try:
        while True:
            data = fdat.next().split()
            if fcoh is not None: dach = fcoh.next().split()
            
            (ikp, isym, nbands, nemin, nomega) = map(int, data[1:6])
            
            ekom = zeros(nbands, dtype=complex)
            dach=ones((nomega,nbands), dtype=complex)
            index=range(nomega)
            omw=zeros(nomega,dtype=float)
            if fcoh is not None:
                for iom in range(nomega):
                    datc = array(map(float,fcoh.next().split()))
                    omw[iom] = datc[0]
                    dach[iom,:] = datc[1::2]+datc[2::2]*1j
                    #print 'shape=', shape(dach), 'nbands=', nbands
                # need to sort frequency because open-mp mixes them up
                index=sorted(index, key=lambda i: omw[i])
                #for i in range(len(index)):
                #    print omw[index[i]],
                #print

            Aom=zeros(nomega,dtype=float)
            om=zeros(nomega,dtype=float)
            for iom in range(nomega):
                data = array(map(float, fdat.next().split()))
                omega = float(data[0])
                ekom = data[1::2]+data[2::2]*1j
                om[iom] = omega
                cohd = dach[index[iom]]
                #print 'om=', omega, omw[index[iom]]
                Aom[iom] = weave.inline(code, ['nbands', 'omega', 'mu', 'ekom', 'small', 'ikp', 'cohd'],
                                        type_converters=weave.converters.blitz, compiler = 'gcc')
            Akom.append( Aom )
            
    except StopIteration:
        pass
    
    Akom = array(Akom).transpose()
    print 'shape(Akom)=', shape(Akom)

    
    vmm = [0,max(map(max,Akom))*intensity]    
    (ymin,ymax) = (om[0]+DY,om[-1]+DY)
    (xmin,xmax) = (0, shape(Akom)[1]-1)
    #(xmin,xmax) = (0, nkp-1)
    
    print 'xmin,xmax,ymin,ymax=', xmin, xmax, ymin, ymax

    imshow(Akom, interpolation='bilinear', cmap=cm.hot, origin='lower', vmin=vmm[0], vmax=vmm[1], extent=[xmin,xmax,ymin,ymax], aspect=(xmax-xmin)*0.8/(ymax-ymin) )

    for i in range(len(wkpointi)):
        print 'wp=', wkpointi[i]
        plot([wkpointi[i],wkpointi[i]], [ymin,ymax], 'w-')
        
    plot([xmin,xmax],[0,0], 'w:')

    dytck=0.005
    Ntck=5
    for j in range(len(wkpointi)-1):
        for ix in range(1,Ntck):
            x = wkpointi[j]+(wkpointi[j+1]-wkpointi[j])*ix/float(Ntck)
            plot([x,x],[-dytck,dytck],'w-')
        
    axis([xmin,xmax,ymin,ymax])
    xticks( wkpointi, wkpoints, fontsize='x-large' )
    #colorbar()
    show()
    

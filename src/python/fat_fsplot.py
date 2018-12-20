#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from scipy import *
from pylab import *
import os, sys
import scipy

from distutils.version import StrictVersion
if StrictVersion(scipy.__version__) > StrictVersion('0.19.0'):
    import weave
else:
    import scipy.weave as weave

code="""
     using namespace std;

     double Ry2eV = 13.6056920311654;

     double Ax = 0;
     for (int ib=0; ib<nbands; ib++){
        //complex<double> ekom( data(1+2*ib),  min(data(2+2*ib), -small) );
        complex<double> ekom( data(1+2*ib),  -small );
        complex<double> cohf( dach(1+2*ib),  dach(2+2*ib) );
        complex<double> gc = cohf/(omega+mu-ekom);
        Ax += -gc.imag()/M_PI;
     }
     return_val = Ax;
     
     """

code0="""
     using namespace std;

     double Ry2eV = 13.6056920311654;

     double Ax = 0;
     for (int ib=0; ib<nbands; ib++){
        //complex<double> ekom( data(1+2*ib),  min(data(2+2*ib), -small) );
        complex<double> ekom( data(1+2*ib),  -small );
        complex<double> gc = 1./(omega+mu-ekom);
        Ax += -gc.imag()/M_PI;
     }
     return_val = Ax;
     
     """

code1="""
     using namespace std;
     bool found = false;
     int ib=0;
     for (ib=0; ib<nbands-1; ib++){
        double ekom_i = data(1+2*ib);
        double ekom_ip = data(1+2*ib+2);
        if ( (omega+mu-ekom_i)*(omega+mu-ekom_ip)<0.0 ){
           found = true;
           break;
        }
     }
     double cohf=0;
     if (found){
        double ekom_i = data(1+2*ib);
        double ekom_ip = data(1+2*ib+2);
        double x_i = omega+mu-ekom_i;
        double x_ip = omega+mu-ekom_ip;
        double coh_i  = dach(1+2*ib)*dach(1+2*ib)+dach(2+2*ib)*dach(2+2*ib);
        double coh_ip = dach(1+2*ib+2)*dach(1+2*ib+2)+dach(2+2*ib+2)*dach(2+2*ib+2);
        cohf = coh_i + (coh_ip-coh_i)*(0-x_i)/(x_ip-x_i);
        //cout<<x_i<<" "<<x_ip<<" "<<cohf<<endl;
     }
     return_val = cohf;
     
     """

def Cols(x0,x1,x2,bright):
    col = array([x0, x1, x2])
    col *= 1./sqrt(dot(col,col))
    col *= bright
    col = array([1-sqrt(col[1]**2+col[2]**2),1-sqrt(col[0]**2+col[2]**2),1-sqrt(col[0]**2+col[1]**2)])
    return col

if __name__ == '__main__':

    small = 1.e-2
    
    fEF = open('EF.dat', 'r')
    mu = float(fEF.next().split()[0])
    print 'mu=', mu
    
    fdat = open('eigvals.dat', 'r')
    #if os.path.isfile('cohfactorsd.dat'):
    #    fcoh = open('cohfactorsd.dat', 'r')
    #else:
    #    fcoh = None

    print sys.argv

    if len(sys.argv)<4:
        print 'Please give 3 files of coherence factors'
        sys.exit(0)
        
    print sys.argv[1]
    print sys.argv[2]
    print sys.argv[3]
    
    fcoh1 = open(sys.argv[1], 'r')
    fcoh2 = open(sys.argv[2], 'r')
    fcoh3 = open(sys.argv[3], 'r')
    
    wdata = fdat.readlines()

    ikp=0
    Akom=[]
    Ch1kw=[]
    Ch2kw=[]
    Ch3kw=[]
    ii=0
    while (ii<len(wdata)):
        
        data = wdata[ii].split()
        #if fcoh is not None: dach = fcoh.next().split()
        
        dach1 = fcoh1.next().split()
        dach2 = fcoh2.next().split()
        dach3 = fcoh3.next().split()
        
        ii += 1
        #print data
        
        (ikp, isym, nbands, nemin, nomega) = map(int, data[1:6])

        ekom = zeros(nbands, dtype=complex)
        #dach = array([0.]+[1.,0.]*nbands)
        Aom=[]
        Ch1w=[]
        Ch2w=[]
        Ch3w=[]
        om=[]
        for iom in range(nomega):
            data = array(map(float, wdata[ii].split()))
            #if fcoh is not None: dach = array(map(float, fcoh.next().split()))
            
            dach1 = array(map(float, fcoh1.next().split()))
            dach2 = array(map(float, fcoh2.next().split()))
            dach3 = array(map(float, fcoh3.next().split()))
            
            ii += 1
            omega = float(data[0])
            
            #Ax = weave.inline(code, ['nbands', 'omega', 'mu', 'data', 'small', 'ikp', 'dach'],
            #                  type_converters=weave.converters.blitz, compiler = 'gcc')
            Ax = weave.inline(code0, ['nbands', 'omega', 'mu', 'data', 'small', 'ikp'],
                              type_converters=weave.converters.blitz, compiler = 'gcc')
            dach = dach1
            Ch1 = weave.inline(code1, ['nbands', 'omega', 'mu', 'data', 'small', 'ikp', 'dach'],
                              type_converters=weave.converters.blitz, compiler = 'gcc')
            dach = dach2
            Ch2 = weave.inline(code1, ['nbands', 'omega', 'mu', 'data', 'small', 'ikp', 'dach'],
                              type_converters=weave.converters.blitz, compiler = 'gcc')
            dach = dach3
            Ch3 = weave.inline(code1, ['nbands', 'omega', 'mu', 'data', 'small', 'ikp', 'dach'],
                              type_converters=weave.converters.blitz, compiler = 'gcc')
            
            om.append(omega)
            Aom.append(Ax)
            Ch1w.append(Ch1)
            Ch2w.append(Ch2)
            Ch3w.append(Ch3)
            
        Akom.append( Aom )
        Ch1kw.append( Ch1w )
        Ch2kw.append( Ch2w )
        Ch3kw.append( Ch3w )

        ikp += 1

    Akom = array(Akom).transpose()
    Ch1kw = array(Ch1kw).transpose()
    Ch2kw = array(Ch2kw).transpose()
    Ch3kw = array(Ch3kw).transpose()
    Akomax=max(Akom[0,:])
    
    nk = shape(Akom)[1]
    nk1 = int((-1+sqrt(1+8*nk))/2)



    FS = zeros((nk1,nk1,3), dtype=float)
    iw=0
    for i in range(nk1):
        for j in range(i,nk1):
            
            col = Cols(Ch1kw[0,iw], Ch2kw[0,iw], Ch3kw[0,iw], Akom[0,iw]/Akomax)
            
            FS[i,j,:] = col
            FS[j,i,:] = col
            iw+=1
    
    
    #vmm = [0,20]
    #xmm = [0, shape(Akom)[1]]
    #ymm = [om[0], om[-1]]
    
    #imshow(Akom, interpolation='bilinear', cmap=cm.hot, origin='lower', vmin=vmm[0], vmax=vmm[1], extent=[xmm[0],xmm[1],ymm[0],ymm[1]], aspect=(xmm[1]-xmm[0])*0.8/(ymm[1]-ymm[0]) )
    #imshow(FS, interpolation='bilinear', cmap=cm.hot, origin='lower' , vmin=0.0, vmax=80.)
    print shape(FS)
    
    imshow(FS, origin='lower',  interpolation='bilinear', extent=[0,1,0,1], aspect=1. )
    
    show()
    
    
    

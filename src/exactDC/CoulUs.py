# @Copyright 2007 Kristjan Haule
from scipy import *
import gaunt
import scipy

from distutils.version import StrictVersion
if StrictVersion(scipy.__version__) > StrictVersion('0.14.0'):
    import weave
else:
    import scipy.weave as weave


def CoulUsC2(l, T2C):
    # Precomputes gaunt coefficients for speed
    # shape(gck)=(4,7,7,4). It contains gck(l, m4, m1, k)
    gck = gaunt.cmp_all_gaunt()

    #print 'Gaunt coefficients precomputed - shape(gck)', shape(gck)
    mw = 2*l+1
    if len(T2C) == mw:
        nw = mw
        ns = 1
    elif len(T2C) == 2*mw:
        nw = 2*(2*l+1)
        ns = 2
    else:
        print "ERROR in atom_d.py: T2C has wrong shape"
        sys.exit(0)
    
    T2Cp = conj(T2C.transpose())
    
    UC = zeros((l+1, nw, nw, nw, nw), dtype=complex)
    shft = 3-l
    shft2 = 2*l+1
    Sum1 = zeros((nw,nw,shft2*2),dtype=complex)
    Sum2 = zeros((nw,nw,shft2*2),dtype=complex)
    
    source="""
    using namespace std;

    for (int k=0; k<l+1; k++){
        Sum1=0;
        for (int i4=0; i4<nw; i4++){
            for (int i1=0; i1<nw; i1++){
                for (int m4=0; m4<mw; m4++){
                    for (int m1=0; m1<mw; m1++){
                        for (int s=0; s<ns; s++) Sum1(i4,i1,m1-m4+shft2) += T2Cp(i4,m4+s*mw)*gck(l,shft+m4,shft+m1,k)*T2C(m1+s*mw,i1);
                    }
                }
            }
        }
        Sum2=0;
        for (int i3=0; i3<nw; i3++){
            for (int i2=0; i2<nw; i2++){
                for (int m3=0; m3<mw; m3++){
                    for (int m2=0; m2<mw; m2++){
                        for (int s=0; s<ns; s++) Sum2(i3,i2,m3-m2+shft2) += T2Cp(i3,m3+s*mw)*gck(l,shft+m2,shft+m3,k)*T2C(m2+s*mw,i2);
                    }
                }
            }
        }
        for (int i4=0; i4<nw; i4++){
            for (int i3=0; i3<nw; i3++){
                for (int i2=0; i2<nw; i2++){
                    for (int i1=0; i1<nw; i1++){
                        complex<double> csum=0.0;
                        for (int dm=0; dm<shft2*2; dm++) csum += Sum1(i4,i1,dm)*Sum2(i3,i2,dm);
                        UC(k,i4,i3,i2,i1) = csum;
                    }
                }
            }
        }
    }
    """

    weave.inline(source, ['UC', 'gck', 'l', 'T2C', 'T2Cp', 'shft', 'nw', 'mw', 'ns', 'shft2', 'Sum1', 'Sum2'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')

    return UC




def CoulUsC2_diagonal(l, T2C):
    # Precomputes gaunt coefficients for speed
    # shape(gck)=(4,7,7,4). It contains gck(l, m4, m1, k)
    gck = gaunt.cmp_all_gaunt()

    #print 'Gaunt coefficients precomputed - shape(gck)', shape(gck)
    mw = 2*l+1
    if len(T2C) == mw:
        nw = mw
        ns = 1
    elif len(T2C) == 2*mw:
        nw = 2*(2*l+1)
        ns = 2
    else:
        print "ERROR in atom_d.py: T2C has wrong shape"
        sys.exit(0)
    
    T2Cp = conj(T2C.transpose())
    
    UC = zeros((l+1, nw, nw), dtype=complex)
    shft = 3-l
    shft2 = 2*l+1
    Sum1 = zeros((nw,shft2*2),dtype=complex)
    Sum2 = zeros((nw,shft2*2),dtype=complex)
    
    source="""
    using namespace std;

    for (int k=0; k<l+1; k++){
        Sum1=0;
        for (int i1=0; i1<nw; i1++){
            for (int m4=0; m4<mw; m4++){
                for (int m1=0; m1<mw; m1++){
                    for (int s=0; s<ns; s++) Sum1(i1,m1-m4+shft2) += T2Cp(i1,m4+s*mw)*gck(l,shft+m4,shft+m1,k)*T2C(m1+s*mw,i1);
                }
            }
        }
        Sum2=0;
        for (int i2=0; i2<nw; i2++){
            for (int m3=0; m3<mw; m3++){
                for (int m2=0; m2<mw; m2++){
                    for (int s=0; s<ns; s++) Sum2(i2,m3-m2+shft2) += T2Cp(i2,m3+s*mw)*gck(l,shft+m2,shft+m3,k)*T2C(m2+s*mw,i2);
                }
            }
        }
        for (int i1=0; i1<nw; i1++){
            for (int i2=0; i2<nw; i2++){
                complex<double> csum=0.0;
                for (int dm=0; dm<shft2*2; dm++) csum += Sum1(i1,dm)*Sum2(i2,dm);
                UC(k,i1,i2) = csum;
            }
        }
    }
    """

    weave.inline(source, ['UC', 'gck', 'l', 'T2C', 'T2Cp', 'shft', 'nw', 'mw', 'ns', 'shft2', 'Sum1', 'Sum2'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')

    return UC



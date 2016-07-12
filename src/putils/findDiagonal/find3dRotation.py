from scipy import *
from scipy import linalg
import copy

def mprint(Us):
    for i in range(shape(Us)[0]):
        for j in range(shape(Us)[1]):
            print "%11.8f %11.8f  " % (real(Us[i,j]), imag(Us[i,j])),
        print

def MakeOrthogonal(a, b, ii):
    a -= (a[ii]/b[ii])*b
    a *= 1/sqrt(dot(a,a.conj()))
    b -= dot(b,a.conj())*a
    b *= 1/sqrt(dot(b,b.conj()))
    return (a,b)

def StringToMatrix(cfstr):
    mm=[]
    for line in cfstr.split('\n'):
        line = line.strip()
        if line:
            data = array(map(float,line.split()))
            mm.append( data[0::2]+data[1::2]*1j )
    mm=matrix(mm)
    return mm


def RealPhase(vec):
    for j in range(len(vec)):
        v = vec[j]
        imax = 0
        vmax = abs(v[imax])
        for i in range(len(v)):
            if abs(v[i])>vmax:
                vmax=abs(v[i])
                imax = i
        vec[j,:] = array(v)*abs(v[imax])/v[imax]
    return vec

def to_normalize(a):
    return 1./sqrt(abs(dot(conj(a), a)))

def swap(a,b):
    an = copy.deepcopy(a)
    bn = copy.deepcopy(b)
    return (bn,an)


def GiveTrans(strHc, strT2C, strT2Crest, dim=6):
    Hc=StringToMatrix(strHc)[:dim,:dim]
    T2C=StringToMatrix(strT2C)
    T2Crest = StringToMatrix(strT2Crest)
    
    ee = linalg.eigh(Hc)
    Es = ee[0]
    Us = matrix(ee[1])
    
    #print 'In Eigensystem:'
    #mprint(Us.H * Hc * Us)
    # Us.H * Hc * Us  === diagonal
    
    print 'Eigenvalues=', Es.tolist()
    
    for i0 in range(0,dim,2): 
        if i0+1<dim and abs(Es[i0]-Es[i0+1])<1e-4: # if two states are degenerate. I need to orthogonalize them

            print 'Degenerate eigenvalues: ', Es[i0], Es[i0+1]
            
            i2=i0+2
            vects = Us[:,i0:i2]
            O = transpose(conj(vects))[0:2,i0:i2]
            (u_,s_,v_) = linalg.svd(O)
            print 'S=', s_.tolist()
            
            m = min(shape(u_)[1],shape(v_)[0])
            R = dot(u_[:,:m],v_[:m,:])
            vectn = dot(vects,R)
            Us[:,i0:i2] = vectn[:,:]
        
    
    #Us = u_ * s_ * v_
    #print 'Eigenvalues'
    #print "%10.5f "*len(Es) % tuple(Es)
    
    print 'Transformation in crystal harmonics='
    mprint(Us)
    print
    
    final = Us.T*T2C
    
    final = array(final)
    final2 = RealPhase(final)
    final=copy.deepcopy(final2)
    return (final, T2C, Hc, T2Crest)

def Check(final, T2C, Hc):
    # the modified final transofrmation is rotated back to t2g-eg base to see how weell diagonal remains
    Us_new = transpose(matrix(final)*T2C.H)
    print 'Check-diagonal:'
    mprint(Us_new.H * Hc * Us_new)
    
    print 'Check unitary:'
    mprint( matrix(final) * matrix(final).H )
    print

def CheckDet(final, T2Crest):
    totalfinal = vstack((final,T2Crest))
    Det = linalg.det(totalfinal)
    print 'Determinant=', Det
    if abs(Det+1)<1e-3:
        print 'Determinant is -1, you need to change an eigenvector, to make the rotation proper!'
    return Det




if __name__ == '__main__':

    strHc1="""
    0.57872164    -0.00000000       0.28797188     0.00000000      -0.00000000     0.28799122   
    0.28797188    -0.00000000       0.57871482    -0.00000000       0.00000000     0.28798664   
   -0.00000000    -0.28799122       0.00000000    -0.28798664       0.57873288     0.00000000   
    """

    strT2C="""
     0.00000000 0.00000000   0.70710679 0.00000000   0.00000000 0.00000000  -0.70710679 0.00000000   0.00000000 0.00000000 
     0.00000000 0.00000000   0.00000000 0.70710679   0.00000000 0.00000000   0.00000000 0.70710679   0.00000000 0.00000000 
    -0.70710679 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.70710679 0.00000000 
    """
    strT2Crest="""
     0.00000000 0.00000000   0.00000000 0.00000000   1.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000 
     0.70710679 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.70710679 0.00000000 
    """

    Hc1=StringToMatrix(strHc1)
    T2C=StringToMatrix(strT2C)
    T2Crest = StringToMatrix(strT2Crest)
    
    print len(Hc1)
    
    (final1, T2C, Hc1, T2Crest) = GiveTrans(strHc1, strT2C, strT2Crest, len(Hc1))
    
    print 'Rotation to input : '
    mprint( final1 )
    print
    print 'rest='
    mprint( T2Crest )
    

# @Copyright 2007 Kristjan Haule
from scipy import *
from scipy import linalg
import sys
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
        #print 'checkin j'
        imax = 0
        vmax = abs(v[imax])
        for i in range(len(v)):
            if abs(v[i])>vmax:
                vmax=abs(v[i])
                imax = i
        #print 'imax', imax, v[imax]
        vec[j,:] = array(v)*abs(v[imax])/v[imax]
    return vec

def to_normalize(a):
    return 1./sqrt(abs(dot(conj(a), a)))

def swap(a,b):
    an = copy.deepcopy(a)
    bn = copy.deepcopy(b)
    return (bn,an)

strHc1="""
   -2.00018980     0.00000000       0.00000000     0.00000000      -0.00004096    -0.24670234       0.00000000     0.00000000       0.00000000     0.00000000       0.17125395     0.18220260   
    0.00000000     0.00000000      -2.00018979     0.00000000       0.00000000     0.00000000      -0.00004096     0.24670234       0.17125395    -0.18220260       0.00000000     0.00000000   
   -0.00004096     0.24670234       0.00000000     0.00000000      -1.99983239     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000      -0.18221011     0.17125654   
    0.00000000     0.00000000      -0.00004096    -0.24670234       0.00000000     0.00000000      -1.99983239     0.00000000      -0.18221011    -0.17125654       0.00000000     0.00000000   
    0.00000000     0.00000000       0.17125395     0.18220260       0.00000000     0.00000000      -0.18221011     0.17125654      -1.65479986     0.00000000       0.00000000     0.00000000   
    0.17125395    -0.18220260       0.00000000     0.00000000      -0.18221011    -0.17125654       0.00000000     0.00000000       0.00000000     0.00000000      -1.65479986     0.00000000   
"""
strHc2="""
   -2.00018980     0.00000000       0.00000000     0.00000000      -0.00004096    -0.24670234       0.00000000     0.00000000       0.00000000     0.00000000       0.17125395     0.18220260   
    0.00000000     0.00000000      -2.00018979     0.00000000       0.00000000     0.00000000      -0.00004096     0.24670234       0.17125395    -0.18220260       0.00000000     0.00000000   
   -0.00004096     0.24670234       0.00000000     0.00000000      -1.99983239     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000      -0.18221011     0.17125654   
    0.00000000     0.00000000      -0.00004096    -0.24670234       0.00000000     0.00000000      -1.99983239     0.00000000      -0.18221011    -0.17125654       0.00000000     0.00000000   
    0.00000000     0.00000000       0.17125395     0.18220260       0.00000000     0.00000000      -0.18221011     0.17125654      -1.65479986     0.00000000       0.00000000     0.00000000   
    0.17125395    -0.18220260       0.00000000     0.00000000      -0.18221011    -0.17125654       0.00000000     0.00000000       0.00000000     0.00000000      -1.65479986     0.00000000   
"""
strHc3="""
   -2.00016629     0.00000000       0.00000000     0.00000000      -0.00008153    -0.24670254       0.00000000     0.00000000       0.00000000     0.00000000       0.17125296     0.18219979   
    0.00000000     0.00000000      -2.00016991     0.00000000       0.00000000     0.00000000      -0.00008153     0.24670261       0.17125291    -0.18219987       0.00000000     0.00000000   
   -0.00008153     0.24670254       0.00000000     0.00000000      -1.99961597     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000      -0.18220613     0.17125778   
    0.00000000     0.00000000      -0.00008153    -0.24670261       0.00000000     0.00000000      -1.99961959     0.00000000      -0.18220620    -0.17125773       0.00000000     0.00000000   
    0.00000000     0.00000000       0.17125291     0.18219987       0.00000000     0.00000000      -0.18220620     0.17125773      -1.65468685     0.00000000       0.00000000     0.00000000   
    0.17125296    -0.18219979       0.00000000     0.00000000      -0.18220613    -0.17125778       0.00000000     0.00000000       0.00000000     0.00000000      -1.65469210     0.00000000   
"""
strHc4="""
   -2.00016629     0.00000000       0.00000000     0.00000000      -0.00008153    -0.24670254       0.00000000     0.00000000       0.00000000     0.00000000       0.17125296     0.18219979   
    0.00000000     0.00000000      -2.00016991     0.00000000       0.00000000     0.00000000      -0.00008153     0.24670261       0.17125291    -0.18219987       0.00000000     0.00000000   
   -0.00008153     0.24670254       0.00000000     0.00000000      -1.99961597     0.00000000       0.00000000     0.00000000       0.00000000     0.00000000      -0.18220613     0.17125778   
    0.00000000     0.00000000      -0.00008153    -0.24670261       0.00000000     0.00000000      -1.99961959     0.00000000      -0.18220620    -0.17125773       0.00000000     0.00000000   
    0.00000000     0.00000000       0.17125291     0.18219987       0.00000000     0.00000000      -0.18220620     0.17125773      -1.65468685     0.00000000       0.00000000     0.00000000   
    0.17125296    -0.18219979       0.00000000     0.00000000      -0.18220613    -0.17125778       0.00000000     0.00000000       0.00000000     0.00000000      -1.65469210     0.00000000   
"""


strT2C="""
 0.00000000 0.00000000   0.70710679 0.00000000   0.00000000 0.00000000  -0.70710679 0.00000000   0.00000000 0.00000000    0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000
 0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000    0.00000000 0.00000000   0.70710679 0.00000000   0.00000000 0.00000000  -0.70710679 0.00000000   0.00000000 0.00000000  
 0.00000000 0.00000000   0.00000000 0.70710679   0.00000000 0.00000000   0.00000000 0.70710679   0.00000000 0.00000000    0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000
 0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000    0.00000000 0.00000000   0.00000000 0.70710679   0.00000000 0.00000000   0.00000000 0.70710679   0.00000000 0.00000000  
-0.70710679 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.70710679 0.00000000    0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000
 0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   -0.70710679 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.70710679 0.00000000
"""
strT2Crest="""
 0.00000000 0.00000000   0.00000000 0.00000000   1.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000    0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000
 0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000    0.00000000 0.00000000   0.00000000 0.00000000   1.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000 
 0.70710679 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.70710679 0.00000000    0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000
 0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000    0.70710679 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.70710679 0.00000000
"""

# -2-1 0 1 2  -2-1 0 1 2
# [0,1,2,3,4,  5,6,7,8,9]
# [9,8,7,6,5,  4,3,2,1,0]
# time-reversal


def GiveTrans(strHc, strT2C, strT2Crest):
    Hc=StringToMatrix(strHc)[:6,:6]
    T2C=StringToMatrix(strT2C)
    T2Crest = StringToMatrix(strT2Crest)
    
    #print 'T2Crest='
    #mprint(T2Crest)
    #print
    
    ee = linalg.eigh(Hc)
    Es = ee[0]
    Us = matrix(ee[1])
    
    #print 'In Eigensystem:'
    #mprint(Us.H * Hc * Us)
    # Us.H * Hc * Us  === diagonal
    
    print 'Eigenvalues=', Es.tolist()
    
    for i0 in range(0,6,2):
        i2=i0+2
        vects = Us[:,i0:i2]
        O = transpose(conj(vects))[0:2,i0:i2]
        (u_,s_,v_) = linalg.svd(O)
        print 'S=', s_.tolist()
        
        m = min(shape(u_)[1],shape(v_)[0])
        R = dot(u_[:,:m],v_[:m,:])
        vectn = dot(vects,R)
        Us[:,i0:i2] = vectn[:,:]
        #print
        #mprint( vectn )
    
    #Us = u_ * s_ * v_
    
    
    #c_ = zeros((shape(u_)[0],shape(v_)[1]),dtype=complex)
    #for i in range(shape(u_)[0]):
    #    for l in range(shape(v_)[1]):
    #        for j in range(shape(u_)[1]):
    #            c_[i,l] += u_[i,j]*s_[j]*v_[j,l]
    
    #print 'Eigenvalues'
    #print "%10.5f "*len(Es) % tuple(Es)
    
    print 'Transformation in crystal harmonics='
    mprint(Us)
    print
    
    final = Us.T*T2C
    
    final = array(final)
    
    
    
    final2 = RealPhase(final)

    # final2[0,:] -= final2[1,:] * final2[0,5]/final2[1,5]
    # final2[1,:] -= final2[0,:] * final2[1,0]/final2[0,0]
    # final2[0,:] *= to_normalize(final2[0,:])
    # final2[1,:] *= to_normalize(final2[1,:])
    # 
    # final2[2,:] -= final2[3,:] * final2[2,5]/final2[3,5]
    # final2[3,:] -= final2[2,:] * final2[3,0]/final2[2,0]
    # final2[2,:] *= to_normalize(final2[2,:])
    # final2[3,:] *= to_normalize(final2[3,:])
    # 
    # final2[4,:] -= final2[5,:] * final2[4,5]/final2[5,5]
    # final2[5,:] -= final2[4,:] * final2[5,0]/final2[4,0]
    # final2[4,:] *= to_normalize(final2[4,:])
    # final2[5,:] *= to_normalize(final2[5,:])
    # final2[0:2,:] = swap(final2[0,:],final2[1,:])
    
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


(final1, T2C, Hc1, T2Crest) = GiveTrans(strHc1, strT2C, strT2Crest)
(final2, T2C, Hc2, T2Crest) = GiveTrans(strHc2, strT2C, strT2Crest)
(final3, T2C, Hc3, T2Crest) = GiveTrans(strHc3, strT2C, strT2Crest)
(final4, T2C, Hc4, T2Crest) = GiveTrans(strHc4, strT2C, strT2Crest)
    



print 'Rotation to input : '
mprint( final1 )
print
mprint( final2 )
print
mprint( final3 )
print
mprint( final4 )
print

print 'rest='
mprint( T2Crest )






# correcting the transformation to have higher symmetry
# there is some freedom due to degenerate eigenvalues
# first the phase factors
#  final[0,:] *= 1j
#  final[1,:] *= abs(final[1,8])/final[1,8]
#  final[2,:] *= -1j
#  cc= final[3,6]/abs(final[3,6])
#  final[3,:] *= 1/cc
#  final[4,:] *= -1j
#  cc = final[5,6]/abs(final[5,6])
#  final[5,:] *= 1/cc
#  # second linear combinations
#  final[2,:],final[3,:] = MakeOrthogonal(final[2,:], final[3,:], 0)
#  final[4,:],final[5,:] = MakeOrthogonal(final[4,:], final[5,:], 0)


# final[0,:] *= abs(final[0,8])/final[0,8]
# final[1,:] *= 1j
# final[3,:] *= -1j
# cc= final[2,6]/abs(final[2,6])
# final[2,:] *= 1/cc
# cc= final[4,6]/abs(final[4,6])
# final[4,:] *= 1/cc
# final[5,:] *= 1j


    
#a = final2[0]
#b = final2[1]
#print 'a*b=', dot(conj(a),b)
#print 'a*a=', dot(conj(a),a)
#print 'b*b=', dot(conj(b),b)

#print a.tolist()
#print b.tolist()
#print a[::-1].tolist()

#c1 = (b[::-1]-a)
#c2 = (b+a[::-1])
#for i in range(len(a)):
#    print "%11.8f %11.8f  " % (c1[i].real,c1[i].imag),
#print
#for i in range(len(a)):
#    print "%11.8f %11.8f  " % (c2[i].real,c2[i].imag),
#print
#print c2/c1



# y=sqrt(1-x^2)
# a_new =  x*a+y*b
# b_new = -y.c*a+x.c*b
# 
# T(x*a+y*b) = (-y.c*a+x.c*b)
# T(-y.c*a+x.c*b) = x*a+y*b
#
# x.c*Ta + y.c*Tb = -y.c*a+x.c*b
# -y*Ta+x*Tb = x*a+y*b
#
# x.c*(b-Ta) = y.c*(Tb+a)
# x*(Tb-a) = y*(b+Ta)
#



# y*(Tb+a) = x*(b-Ta)
# x*(Tb-a) = y*(b+Ta)
# 
# x/y = (Tb+a)/(b-Ta)
# (x/y)**2 + 1 = 1/y**2
# 
# y = sqrt(1/(1+(x/y)**2))

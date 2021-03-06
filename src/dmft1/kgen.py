from scipy import *
import sys

dM = 120
dN = 8

#kpath=[[0,0,0],[0.5,0.5,0],[0.75,0.5,0.25],[1,0.5,0.5],[0.75,0.375,0.375],[0,0,0],[0.5,0.5,0.5]]
kpath=[[0,0,0],[1.,0,0],[0.75,0.5,0.25],[1,0.5,0.5],[0.75,0.375,0.375],[0,0,0],[0.5,0.5,0.5]]
spath=['GAMMA', 'X', 'W', 'X', 'K', 'GAMMA', 'L']

Kb = array([[-1,1,1],[1,-1,1],[1,1,-1]])

#for k in kpath:
#    w=zeros(3)
#    for i in range(3):
#        w += k[i]*Kb[i]
#    print w
#sys.exit(0)




for k in kpath:
    for i in range(3):
        if abs(k[i]*dN - int(k[i]*dN)) > 1e-3:
            print 'dN should be changed: ', dN, k[i], dN*k[i]
            sys.exit(0)


kw0 = zeros(3,dtype=int)
for i in range(3):
    kw0[i] = int(kpath[0][i]*dN+0.5)*dM

for ik in range(1,len(kpath)):
    k = kpath[ik]
    
    kw = zeros(3,dtype=int)
    for i in range(3):
        kw[i] = int(k[i]*dN+0.5)*dM
    dkw = kw - kw0

    for i in range(dM):
        kk = [kw0[0] + dkw[0]*i/dM, kw0[1] + dkw[1]*i/dM, kw0[2] + dkw[2]*i/dM ]
        #print array(kk), dM*dN
        kname=' '
        if i==0: kname = spath[ik-1]
        if (ik==1 and i==0):
            print '%-10s%5d%5d%5d%5d%5.1f%5.2f%5.2f'%(kname,kk[0],kk[1],kk[2],dM*dN,2.0,-8.0,8.0),
            print '   k-list generated by python'
        else:
            print '%-10s%5d%5d%5d%5d%5.1f'%(kname,kk[0],kk[1],kk[2],dM*dN,2.0)
            
    kw0 = kw

kname = spath[-1]
print '%-10s%5d%5d%5d%5d%5.1f'%(kname,kw[0],kw[1],kw[2],dM*dN,2.0)
print 'END'

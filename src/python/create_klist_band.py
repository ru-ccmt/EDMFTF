from scipy import *

P = [ ('$Y$', [0,1,0]),
      ('$\Gamma$', [0,0,0]),
      ('$X$', [1,0,0])
     ]

Ni = [50,50]

SCALE = max(Ni)
l=0
for i in range(len(P)-1):
    p0 = array(P[i][1])
    p1 = array(P[i+1][1])
    for j in range(Ni[i]):
        l+=1
        kk = (p0 + (p1-p0)*j/(Ni[i]+0.0))*SCALE
        if j==0: NAME = P[i][0]
        else: NAME = ''
        print "%-10s%5d%5d%5d%5d%5.1f" % (NAME, kk[0], kk[1], kk[2], SCALE, 1.0)

NAME = P[-1][0]
kk = array(P[-1][1])*SCALE
print "%-10s%5d%5d%5d%5d%5.1f" % (NAME, kk[0], kk[1], kk[2], SCALE, 1.0)

print 'END'


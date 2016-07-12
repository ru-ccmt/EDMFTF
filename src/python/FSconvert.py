from scipy import *
import sys, re
from utils import W2kEnvironment, Ry_in_eV

# Finds what is case
w2k = W2kEnvironment()

band0=0
band1=50

f = open(w2k.case+'.outputkgen', 'r')

head0=27 #23
head1=8 # 4


for line in f:
    if re.search('G1\s*G2\s*G3',line) is not None: break

bs=[]
for i in range(3):
    bs.append( map(float,f.next().split()) )
bs=array(bs).transpose()
b0 = bs[0]
b1 = bs[1]
b2 = bs[2]
print b0, b1, b2

for line in f:
    if re.search('NO.\s*OF\s*MESH\s*POINTS\s*IN\s*THE\s*BRILLOUIN\s*ZONE', line) is not None: break


nk = map(int, f.next().split()[6:10])
f.next()

print 'nk=', nk

index = zeros((nk[0]+1,nk[1]+1,nk[2]+1),dtype=int)
for k0 in range(nk[0]+1):
    for k1 in range(nk[1]+1):
        for k2 in range(nk[2]+1):
            data = f.next().split()
            #print data, k0, k1, k2
            (ii, i0, i1, i2, tindex) = map(int, data[:5])
            if (k0!=i0): print 'ERROR0', k0, i0
            if (k1!=i1): print 'ERROR1', k1, i1
            if (k2!=i2): print 'ERROR2', k2, i2
            index[k0,k1,k2] = tindex
            #print k0, k1, k2, tindex

f.next()
Nkmax = max(index.flatten())

wind={}
for ik in range(Nkmax):
    line = f.next()
    if line[2:5]=='sum': break
    (r0,r1) = map(int,line.split()[:2])
    if r0-1!=ik: print 'ERROR3', r0, ik
    wind[r1]=r0-1
    #print r0, r1

for k0 in range(nk[0]+1):
    for k1 in range(nk[1]+1):
        for k2 in range(nk[2]+1):
            index[k0,k1,k2] = wind[index[k0,k1,k2]]
            #print "%3d %3d %3d %4d" % (k0, k1, k2, index[k0,k1,k2])

Nkmax = max(index.flatten())+1
print Nkmax


g = open('eigvals.dat')
n_emin=0
n_emax=1000
for ik in range(Nkmax):
    (ikp, isym, nbands, nemin, nomega) = map(int, g.next().split()[1:6])
    if ikp-1!=ik: print 'ERROR4', ikp-1, ik
    if isym!=1: print 'ERROR isym!=1'
    for im in range(nomega): g.next()
    n_emin = max(n_emin, nemin)
    n_emax = min(n_emax, nemin+nbands)

print 'n_emin(eigvals)=', n_emin
print 'n_emax(eigvals)=', n_emax

g.close()
g = open('eigvals.dat')
eigenvals = zeros((Nkmax,n_emax-n_emin),dtype=float)
for ik in range(Nkmax):
    (ikp, isym, nbands, nemin, nomega) = map(int, g.next().split()[1:6])
    if ikp-1!=ik: print 'ERROR4', ikp-1, ik
    if isym!=1: print 'ERROR isym!=1'
    dat = map(float,g.next().split())
    for im in range(1,nomega): g.next()
    evals = [dat[1+2*i] for i in range(nbands)]
    eigenvals[ik,:] = evals[n_emin-nemin:n_emax-nemin]


ef = open('EF.dat')
EF = float(ef.next())
ef.close()

#print 'eigvals=', eigenvals
print 'EF=', EF

band0 = max(band0,0)
band1 = min(band1,nbands)


h = open(w2k.case+'.bxsf', 'w')

nbands=n_emax-n_emin
band1 = min(band1,n_emax-n_emin-1)
print >> h, ' BEGIN_INFO'
print >> h, '   #'
print >> h, '   # Launch as: xcrysden --bxsf '+w2k.case+'.bxsf'
print >> h, '   #'
print >> h, '   Fermi Energy:', '%7.5f' % (EF/Ry_in_eV)
print >> h, ' END_INFO'
print >> h
print >> h, ' BEGIN_BLOCK_BANDGRID_3D'
print >> h, '   bands_energies'
print >> h, '   BEGIN_BANDGRID_3D_BANDS'
print >> h, '    ', band1-band0+1
print >> h, '    ', nk[0]+1, nk[1]+1, nk[2]+1
print >> h, '  ', '%10.6f '*3 % (0.0, 0.0, 0.0)
print >> h, '  ', '%10.6f '*3 % tuple(b0)
print >> h, '  ', '%10.6f '*3 % tuple(b1)
print >> h, '  ', '%10.6f '*3 % tuple(b2)

#print 'band0=', band0, 'band1=', band1
print 'band0=', band0, 'band1=', band1

for iband in range(band0,band1+1):
    print >> h, '    BAND:', iband+1
    for k0 in range(nk[0]+1):
        for k1 in range(nk[1]+1):
            print >> h, '    ',
            for k2 in range(nk[2]+1):
                #print shape(eigenvals), index[k0,k1,k2], iband
                de=0.0
                print >> h, "%10.6f " % (eigenvals[index[k0,k1,k2],iband]/Ry_in_eV+de),
            print >> h
        if k0<nk[0]: print >> h
        
print >> h, '   END_BANDGRID_3D'
print >> h, ' END_BLOCK_BANDGRID_3D'

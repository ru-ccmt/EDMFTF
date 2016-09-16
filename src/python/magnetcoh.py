# @Copyright 2007 Kristjan Haule
import sys, re, os
from scipy import *
from pylab import *
from scipy.lib.blas import fblas
import rdU, utils, findEF, struct1

########################################
# input that we need to have
delta = array([0.0, 0.5,0.5])

# input for LDA plot if needed
plotLDA = True
gamma = 0.05
omega = linspace(-1,1,100)

#########################################

# Vector-file handle
vtape=9

w2k = utils.W2kEnvironment()
print 'case=', w2k.case
EF = findEF.FindChemicalPotential(w2k.case, '')[0]
struct = struct1.Struct(w2k.case)

vectortype=float
vecread = rdU.fvread3
so=''
if os.path.isfile(w2k.SCRATCH+"/"+w2k.case+".vectorso") and os.path.getsize(w2k.SCRATCH+"/"+w2k.case+".vectorso")>0 :
    print 'Found '+w2k.case+'.vectorso file, hence assuming so-coupling exists. Switching -so switch!'
    so = 'so'
    vectortype=complex
    vecread = rdU.fvread3c

fc = open('cohfactors.dat', 'w')
print >> fc, '# Coherence factors for KS states'

maxkpoints = 10000
# opens vector file
rdU.fvopen(vtape, w2k.case+'.vector'+so, struct.nat) 
gw=[]
for ik in range(maxkpoints):
    # Reads vector file
    (k, kname, wgh, ios, n0, nb) = rdU.fvread1(vtape)
    if ios!=0: break # vector file is finished, no more k-points
    print 'k=', k

    # Reciprocal vectors
    Gs = rdU.fvread2(vtape, n0)

    # Reading KS eigensystem
    As=zeros((nb,n0), dtype=vectortype)
    Ek=zeros(nb, dtype=float)
    for i in range(nb):
        #(num, ek, A) = rdU.fvread3(vtape, n0)
        (num, ek, A) = vecread(vtape, n0)
        As[i,:] = A       # KS eigenvector
        Ek[i] = utils.Ry2eV(ek) # KS eigenvalue

    # Building overlap
    # A = u * s * v
    # A * O * A^{-1} = 1
    # Op = v^+ * 1/s * u^+ * u * 1/s * v
    (u,s,v) = linalg.svd( array(As) )
    uu = dot(transpose(conjugate(u)), u)
    for i in range(len(uu)):
        for j in range(len(uu)):
            uu[i,j] = 1/s[i] * uu[i,j] * 1/s[j]

    Op = dot(dot(transpose(conjugate(v))[:,:len(uu)], uu), v[:len(uu),:])


    # Building phase factors for the two sublatices
    exph=zeros(n0, dtype=float)
    for ig in range(n0):
        ph = dot(Gs[:,ig], delta)
        exph[ig] = exp(-1j*ph*2*pi)
        #print "%3d %4.1f %8.5f " % (ig, ph, exph[ig])
    # We want to use matrix multiplications only. Need the eigenvector with the phase factor As1.
    As0 = array(As, dtype=complex)
    As1 = zeros(shape(As), dtype=complex)
    for ig in range(n0):
        As1[:,ig] = As[:,ig]*exph[ig]

    gs = identity(nb) + dot(dot(As1, Op), transpose(conjugate(As1))) + dot(dot(As1, Op), transpose(conjugate(As0))) + dot(dot(As0, Op), transpose(conjugate(As1)))

    # Printing the coherence factors
    print >> fc, ik+1, 1, len(gs)
    for i in range(len(gs)):
        for j in range(len(gs)):
            print >> fc, "%15.10f " % gs[j,i].real,
        print >> fc
    
    if plotLDA:
        gt = zeros(len(omega), dtype=float)
        for iom,om in enumerate(omega):
            dsum = 0.0
            for i in range(nb):
                dsum -= gs[i,i]*(1./(om+gamma*1j+EF-Ek[i])).imag
            gt[iom] = dsum
        gw.append( gt )

nkp=ik
print 'nkp=', nkp

if plotLDA:
    nkp=ik
    gw = transpose(array(gw))
    
    xmm = [0, nkp]
    ymm = [omega[0], omega[-1]]
    
    imshow(gw, interpolation='bilinear', cmap=cm.hot, origin='lower', extent=[xmm[0],xmm[1],ymm[0],ymm[1]], aspect=(xmm[1]-xmm[0])*0.8/(ymm[1]-ymm[0]) )
    show()

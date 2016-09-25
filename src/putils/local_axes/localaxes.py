#!/usr/bin/env python
from numpy import * 
from scipy import linalg
from scipy.linalg import norm
from utils import W2kEnvironment
import sys
from numpy.linalg import svd
import w2k_nn

def ReadNNFile(case):
    f = open(case+'.outputnn_', 'r')
    lines = f.readlines()
    f.close()
    lines = [line.strip() for line in lines if (not line.startswith(' RMT(')) and (not line.startswith(' SUMS TO'))]
    return lines

def GetUserInput(lines):
    print 'ATOMIC INFORMATION FROM STRUCT FILE'
    print

    headers = [(n,line) for n,line in enumerate(lines) if ('ATOM:' in line) and ('EQUIV.' in line)]
    for n,header in enumerate(headers):
        linenum, text = header
        print '%3d %s' % (n+1, text)
    
    atoms = raw_input('Enter the atom for which you want to determine local axes: ')
    return int(atoms)
    
def get_bravais_matrix2(case):
    """Returns a matrix S2C, which converts vector in latice units to cartesian coordinates"""
    extn = '.rotlm_'
    fi = open(case+extn, 'r')
    fi.next()
    BR1=[]
    for i in range(3):
        BR1.append( map(float,fi.next().split()) )
    BR1=array(BR1)
    fi.close()
    S2C = linalg.inv(BR1)*2*pi
    return S2C

def FindCageBasis(iatom,S2C,lines):
    """Find the best directions for the local coordinate system, which is non-orthogonal.
    """
    # supported cages which are currently detected
    # (type, N, angles)
    cases = [('cube',8,[arccos(1/3.)]*3+[arccos(-1/3.)]*3+[pi]), ('octahedra',6,[pi/2]*4+[pi]), ('tetrahedra',4,[arccos(-1/3.)]*3)]

    
    headers = [(n,line) for n,line in enumerate(lines) if ('ATOM:' in line) and ('EQUIV.' in line)]

    startline, text = headers[iatom-1]
    x,y,z = [float(coord) for coord in text.split()[-3:]]  # coordinates of central atom

    # Reads all atoms around the central atom
    neighs = [] # these are in lattice units
    cneighs=[]  # these are in cartesian
    for line in lines[startline+1:]:
        #print line
        if not line: break # empty line stops the environment of this atom
        neigh = [float(line[23:37])-x, float(line[37:51])-y, float(line[51:65])-z]
        cneigh =  dot(S2C,neigh) # convert coordinates into orthonormal basis
        neighs.append(neigh)
        cneighs.append( cneigh )

    
    lis0=[] # precompute all distances
    negs0=[] # precompute all unit vectors
    for neigh in cneighs:
        d = linalg.norm(neigh)
        lis0.append( d )
        negs0.append( neigh/d )
        
    # Now go over all atoms in the cage, and guess what type of cage it is
    criteria=[]
    lmin = lis0[0]
    for (ctype,N,angles) in cases: # over possible environments: cube, octahedra, tetrahedra
        lis = lis0[:N]
        negs = negs0[:N]
        # First, compute effective coordination number (see VESTA manual)
        zw=0.
        zi=0.
        for li in lis:
            w = exp(1-(li/lmin)**6)
            zi += li*w
            zw += w
        lav = zi/zw
        ws = 0
        for li in lis:
            ws += exp(1-(li/lav)**6)
        # Second, compute all angles between all vectors
        angs=[]
        for i in range(N):
            an=[]
            for j in range(N):
                if i==j: continue
                cosp = dot(negs[i],negs[j])
                if cosp>1: cosp=1
                if cosp<-1: cosp=-1
                phi = arccos(cosp)
                an.append(phi)
                #print i, j, phi*180/pi, cosp
            angs.append(an)
        angs = [array(sorted(ang)) for ang in angs] # These are all angles, sorted
        #
        angles = array(angles)
        chi2 = 0 # Now we can compute bond angle variance, i.e., distance between expected angles and actual angles
        for i in range(len(angs)):
            chi2 += sum((angs[i]-angles)**2)
        chi2 *= 1./(len(angs)*len(angles))

        criteria.append( (ctype, abs(ws-N), chi2) ) # now remember how good is this guess in terms of coordination number and angle variance
        print 'trying ', ctype, ': <w>-N=', ws-N, 'chi2=', chi2, '<l>=', lav
        #print 'phi='
        #for i in range(len(angs)):
        #    print angs[i]*180/pi
    # On the basis of these criteria, we should be able to tell which environment we have
    criteria = sorted(criteria, key=lambda ct: ct[2]) # sort according to the angle variance
    # First take the cage, which has most similar angles.
    if criteria[0][1]<0.1: # If the bond variance is not too bad, i.e., |<w>-N|< 0.1, we found the type of cage.
        ctype = criteria[0]
    elif criteria[1][1]<0.1 and criteria[1][2]<0.01: # If the bond variance was very bad for the first case, we go down the list
        ctype = criteria[1]
    else: # If the second is not OK, we boil out at the moment
        print 'Could not detect the type of environment, Boiling out.'
        sys.exit(0)
    print 'Found the environment is', ctype

    if ctype[0]=='octahedra':
        # Now that we know it is octahedra, we take all vectors to atoms, and construct the best coordinate system that
        # goes through these atoms
        nis = negs0[:6] # all unit vectors
        R0=[]
        for i in range(3):
            nu = nis.pop(0) # first unit vector
            cosp = [dot(nu,ni) for ni in nis] # which one is opposite?
            which = argmin(cosp) # this is index of the opposite
            nv = nis.pop(which)  # nv is opposite to nu
            n1 = nu-nv           # now take the average of the two, which are opposite
            n1 *= 1./linalg.norm(n1) # and normalize
            R0.append(n1)        # we have a good unit vector
            #print n1, nu, nv, arccos(cosp[which])*180/pi
        return R0
    elif ctype[0]=='tetrahedra':
        # Now that we know it is tetrahedra, we take vectors which go through the middle of the two atoms, and construct
        # the best coordinate system with such vectors.
        n0 = negs0[:4] # all unit vectors
        nis=[]
        for i in range(4):
            for j in range(i+1,4):
                ni = n0[i]+n0[j]           # new unit vectors are between each pair of vertices
                ni *=  1./linalg.norm(ni)  # and normalize
                nis.append(ni)
        R0=[]
        for i in range(3):
            nu = nis.pop(0) # first unit vector
            cosp = [dot(nu,ni) for ni in nis] # which one is opposite?
            which = argmin(cosp) # this is index of the opposite
            nv = nis.pop(which)  # nv is opposite to nu
            n1 = nu-nv           # now take the average of the two, which are opposite
            n1 *= 1./linalg.norm(n1) # and normalize
            R0.append(n1)       # we have a good unit vector
            #print n1, nu, nv, arccos(cosp[which])*180/pi
        return R0
    elif ctype=='cube':
        # For each vector to atom, we should find 3 other atoms which have smallest angle between them
        # Now we have four vectors with small angle, which define the phase of a cube.
        # We then sum this four vectors, and get unit vector along x,y, or z
        print 'Noy yet implemented'
        sys.exit(0)
    else:
        print 'Noy yet implemented'
        sys.exit(0)

def ResortToDiagonal(R):
    # We now resort, so that the rotation is close to identity
    # This is not necessary, but is convenient
    permutations=[(0,1,2),(0,2,1),(1,0,2),(1,2,0),(2,1,0),(2,0,1)]
    ii, wi = 0, 1000
    for ip,p in enumerate(permutations):
        Rt = array( [R[p[0],:], R[p[1],:], R[p[2],:]] )
        wj = sum([abs(abs(Rt[i,i])-1.) for i in range(3)])
        if wj<wi:
            ii=ip
            wi=wj
    p=permutations[ii]
    Rt = array( [R[p[0],:], R[p[1],:], R[p[2],:]] )
    Rt *= linalg.det(Rt)
    return Rt
                
if __name__ == '__main__':
    case=W2kEnvironment().case

    # executes nn in python module to obtain case.outputnn_ and case.rotlm_
    w2k_nn.w2knn(case)
    
    # reads information from case.outputnn_
    lines = ReadNNFile(case)
    # get conversion from lattice to cartesian coordinates
    S2C = get_bravais_matrix2(case)
    # gets users input
    iatom  = GetUserInput(lines)

    # Main part of the algorithm
    R0 = FindCageBasis(iatom,S2C,lines)
    # Now orthogonalizing the set of vectors
    U,S,V = linalg.svd(R0)
    R = dot(U,V)
    
    R = ResortToDiagonal(R)
    
    print
    print 'Rotation to input into case.indmfl by locrot=-1 : '
    print
    print "%12.8f "*3 % tuple(R[0,:])
    print "%12.8f "*3 % tuple(R[1,:])
    print "%12.8f "*3 % tuple(R[2,:])
    print

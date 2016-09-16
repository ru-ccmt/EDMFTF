#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule

from numpy import * 
from scipy import linalg, optimize
from scipy.linalg import norm
import operator, optparse
import debug, wienfile
from utils import W2kEnvironment
import sys
from numpy.linalg import svd
import w2k_nn

'''Rotation Convention

We use active rotations in this file.  Example: if a is the vector coordinate
of an atom, then R.a should give something close to (1,0,0) or similar.  Here
a is a column vector, and we multiply R from the left.
'''

def readcages2(case, nnn, atoms=''):
    '''nnn - number of nearest neighbor atoms used to form cage'''
    extn = '.outputnn_'
    f = open(case+extn, 'r')
    lines = f.readlines()
    f.close()

    lines = [line.strip() for line in lines if (not line.startswith(' RMT(')) and (not line.startswith(' SUMS TO'))]

    print 'ATOMIC INFORMATION FROM STRUCT FILE'
    print

    headers = [(n,line) for n,line in enumerate(lines) if ('ATOM:' in line) and ('EQUIV.' in line)]
    for n,header in enumerate(headers):
        linenum, text = header
        print '%3d %s' % (n+1, text)

    if not atoms:
        atoms = raw_input('Enter indexes of atoms for which you want to determine local axes: ')
    chosen_iatoms = expand_intlist(atoms)

    # extract cage coordinates for each atom relative to the central atom
    # these are given in units of lattice vectors
    cages = []
    for iatom in chosen_iatoms:
        startline, text = headers[iatom-1]
        x,y,z = [float(coord) for coord in text.split()[-3:]]  # coordinates of central atom
        offset = 1
        coordlines = lines[startline+offset:startline+offset+nnn]
        cage = [[float(line[23:37])-x, float(line[37:51])-y, float(line[51:65])-z] for line in coordlines]
        cages.append((text, cage))

    # convert cage coordinates into orthonormal basis
    S2C = get_bravais_matrix2(case)
    bmat = S2C.transpose()
    print 'Bravais Matrix:'
    print '[',bmat[0].tolist(),',',bmat[1].tolist(),',',bmat[2].tolist(),']'
    print 'S2C = {'+','.join(['{'+','.join([str(e) for e in v])+'}' for v in S2C])+'};'
    print
    #if print_mathematica:
    #    print 'bmat = {'+','.join(['{'+','.join([str(e) for e in v])+'}' for v in bmat])+'};'
    

    #for i,itaom in enumerate(chosen_iatoms):
    #    text, cage = cages[i]
    #    for atom in cage:
    #        print atom, ";", dot(atom,S2C)

    return [(text, [dot(S2C,atom) for atom in cage]) for text, cage in cages]


def expand_intlist(input):
    '''Expand out any ranges in user input of integer lists.
    Example: input  = "1,2,4-6"
             output = [1, 2, 4, 5, 6]'''
    def parse1(x):
        y = x.split('-')
        return [int(x)] if len(y) == 1 else range(*[int(y[0]), int(y[1])+1])

    return reduce(operator.add, [parse1(x) for x in input.split(',')])


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

def sort_octahedral_cage(cage, coscut = 0.707):
    '''Given an octahedral cage of atoms, sorts them so they are ordered
    in the same manner as the vertices are defined (in the list `octahedron`).

    Tries to place closest atoms in x-y plane.'''
    norms = [linalg.norm(v) for v in cage]
    unitvecs = [asarray(v)/linalg.norm(v) for v in cage]

    # choose closest atom along x-axis
    nx = norms.index(min(norms))

    # find atom along minus-x axis
    xangles = dot(unitvecs, unitvecs[nx])
    [nminusx] = [n for n,angle in enumerate(xangles) if angle < -coscut]

    # the other four must lie in the yz plane
    nyz = [n for n,angle in enumerate(xangles) if (angle < coscut) and (angle > -coscut)]
    _,ny = min([(norms[n], n) for n in nyz])  # choose closest atom in yz plane as y-axis

    # find z, -y, -z vectors by taking cross products with y-vector and dotting with x-vector
    xcrossangles = [dot(unitvecs[nx], cross(unitvecs[ny], unitvecs[n])) for n in nyz]
    [nz] = [n for n,angle in zip(nyz, xcrossangles) if angle > coscut]  # y-axis cross z-axis is along x-axis
    [nminusz] = [n for n,angle in zip(nyz, xcrossangles) if angle < -coscut]  # y-axis cross -z-axis is along -x-axis
    [nminusy] = [n for n in range(len(cage)) if n not in [nx, nminusx, ny, nz, nminusz]]

    return [cage[nx], cage[ny], cage[nz], cage[nminusx], cage[nminusy], cage[nminusz]]

def sort_square_cage(cage, coscut = 0.707):
    '''Places z axis perpendicular to square arrangement of atoms.
    x-axis passes through corner of square.'''
    norms = [linalg.norm(v) for v in cage]
    unitvecs = [asarray(v)/linalg.norm(v) for v in cage]

    # choose closest atom along x-axis
    nx = norms.index(min(norms))

    # find atom along minus-x axis
    relative_angles = dot(unitvecs, unitvecs[nx])
    [nminusx] = [n for n,angle in enumerate(relative_angles) if angle < -coscut]

    # find remaining two atoms (which must lie on y-axis)
    nys = [n for n,angle in enumerate(relative_angles) if (angle < coscut) and (angle > -coscut)]
    _,ny = min([(norms[n], n) for n in nys]) # place next closest atom on positive y-axis
    [nminusy] = set(range(4)) - set([nx, nminusx, ny])

    return [cage[nx], cage[ny], cage[nminusx], cage[nminusy]]

def sort_cubic_cage(cage, coscut = 0.707):
    pass

def octahedral_group():
    '''Generates 3x3 representation of octahedral group in {x,y,z} basis.'''
    angles = [0, pi/2, pi, 3*pi/2]

    # brute force: generate all 64 combinations of rotations
    # and filter out redundant matrices to get the 24 elements
    rep = []
    for a in angles:
        for b in angles:
            for c in angles:
                R0 = mat(rotation(0,0,a)) * mat(rotation(b,0,0)) * mat(rotation(0,0,c)) # z-x-z convention
                R = [[int(round(e)) for e in row] for row in asarray(R0)]  # coerce rotation matrix to integer entries
                if R not in rep:  # filter out redundant matrices
                    rep.append(R)
    return rep

def square_group():
    '''Generates 3x3 representations of symmetries of square.'''
    C4angles = [0, pi/2, pi, 3*pi/2]
    C2angles = [0, pi]
    rep = []
    for a in C2angles:
        for c in C4angles:
            R0 = mat(rotation(0,0,c)) * mat(rotation(a,0,0))
            R = [[int(round(e)) for e in row] for row in asarray(R0)]  # coerce rotation matrix to integer entries
            if R not in rep:  # filter out redundant matrices
                rep.append(R)
    return rep


def print_cage(sortedcage, R):
    print len(sortedcage)
    print 'cage'
    format = '  O ' + '%10f '*3
    for v in dot(sortedcage, R.T):
        print format % tuple(v)


if __name__ == '__main__':
    case=W2kEnvironment().case

    # executes nn in python module to obtain case.outputnn_ and case.rotlm_
    w2k_nn.w2knn(case)
    
    sortfunc = sort_octahedral_cage

    # reads case.outputnn_ and case.rotlm_
    cages=readcages2(case,6)

    for text, cage in cages:
        #print text, 'cage=', cage

        # cage must be sorted so atoms are listed in same order
        # as vertices of polyhedron
        sortedcage = sortfunc(cage)
        
        nv=[]
        for v in sortedcage:
            nv.append( v/sqrt(dot(v,v)) )
            
        R0 = array([nv[0],nv[1],nv[2]])

        U,S,V = linalg.svd(R0)
        #print 'S=', S
        R = dot(U,V)
        #print
        print  text
        print 'Rotation to input into case.indmfl by locrot=-1 : '
        print "%12.8f "*3 % tuple(R[0,:])
        print "%12.8f "*3 % tuple(R[1,:])
        print "%12.8f "*3 % tuple(R[2,:])
        print

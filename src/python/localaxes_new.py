#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule

from numpy import asarray,mat,sum,dot,cross
from scipy import linalg, optimize
from scipy.linalg import norm
from math import pi,sin,cos,sqrt
import operator, optparse
import debug, wienfile
from utils import W2kEnvironment


'''Rotation Convention

We use active rotations in this file.  Example: if a is the vector coordinate
of an atom, then R.a should give something close to (1,0,0) or similar.  Here
a is a column vector, and we multiply R from the left.
'''

# software design note:
# in retrospect, this thing should have been a class,
# since i am passing too many extraneous arguments (like `case` or
# `print_mathematica` into all these functions.  


# polyhedra

octahedron = [[ 1,  0,  0],
              [ 0,  1,  0],
              [ 0,  0,  1],
              [-1,  0,  0],
              [ 0, -1,  0],
              [ 0,  0, -1]]

cube = [[ 1,  1,  1],  # not supported yet (need sort_cubic_cage)
        [-1,  1,  1],
        [ 1, -1,  1],
        [ 1,  1, -1],
        [-1, -1,  1],
        [ 1, -1, -1],
        [-1,  1, -1],
        [-1, -1, -1]]

square = [[ 1, 0, 0],
          [ 0, 1, 0],
          [-1, 0, 0],
          [ 0,-1, 0]]

def rotation(a, b, c):
    '''Construct a SO(3) matrix from three parameters by
    exponentiating the following antisymmetric matrix:
    [ 0 -c  b]
    [ c  0 -a]
    [-b  a  0]
    '''
    return linalg.expm([[0, -c, b], [c, 0, -a], [-b, a, 0]])

def utility(x, poly, atoms, radialweight=lambda(r): 1):
    '''Computes utility function for a given atomic cage compared against
    a given polyhedron.
    x     - 3 component vector of parameters defining rotation matrix
    poly  - [(x,y,z), (x,y,z), ...] list of tuples specifying polyhedron vertices
    atoms - [(x,y,z), (x,y,z), ...] list of tuples specifying atoms forming cage
    radialweight - function taking one scalar argument `R` determining how
            strongly an atom at radius `R` should be weighted in utility

    utility = \sum_i f(|a_i|) \frac{x_i \cdot R(a,b,c) \cdot a_i}{|x_i||a_i|}

    f   - radial weighting function
    a_i - vector coordinate of the ith atom forming cage
    x_i - vector coordinate of the ith vertex of polyhedron
    '''
    a,b,c = x
    R = rotation(a,b,c)
    dotprods = sum(asarray(poly).T * dot(R, asarray(atoms).T), axis=0)
    polynorms = asarray([linalg.norm(v) for v in poly])
    atomnorms = asarray([linalg.norm(v) for v in atoms])
    radialweights = asarray([radialweight(r) for r in atomnorms])
    return -sum(radialweights*dotprods/(polynorms*atomnorms))

def expand_intlist(input):
    '''Expand out any ranges in user input of integer lists.
    Example: input  = "1,2,4-6"
             output = [1, 2, 4, 5, 6]'''
    def parse1(x):
        y = x.split('-')
        return [int(x)] if len(y) == 1 else range(*[int(y[0]), int(y[1])+1])

    return reduce(operator.add, [parse1(x) for x in input.split(',')])

def readcages(case, nnn, atoms='', print_mathematica=False):
    '''nnn - number of nearest neighbor atoms used to form cage'''
    extn = '.outputnn'
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
        cage = [[float(line[23:31])-x, float(line[31:39])-y, float(line[39:47])-z] for line in coordlines]
        cages.append((text, cage))

    # convert cage coordinates into orthonormal basis
    bmat = get_bravais_matrix(case)
    print 'Bravais Matrix:'
    print '[',bmat[0].tolist(),',',bmat[1].tolist(),',',bmat[2].tolist(),']'
    S2C = bmat.transpose()
    print 'S2C = {'+','.join(['{'+','.join([str(e) for e in v])+'}' for v in S2C])+'};'
    print
    if print_mathematica:
        print 'bmat = {'+','.join(['{'+','.join([str(e) for e in v])+'}' for v in bmat])+'};'
    return [(text, [dot(atom, bmat) for atom in cage]) for text, cage in cages]

def get_bravais_matrix(case):
    '''Returns a matrix B_{ij} = <i|j> where i indexes the cartesian basis
    vectors and j indexes the crystal lattice vectors.

    For example, in a rhombohedral system, have BR1_DIR =
        2.97706  -5.15641   8.28091
        2.97706   5.15641   8.28091
       -5.95411   0.00000   8.28091

    The crystal lattice vectors are the rows:
        a = ( 2.97706,  -5.15641,   8.28091)  = (1 0 0) . B
        b = ( 2.97706,   5.15641,   8.28091)  = (0 1 0) . B
        c = (-5.95411,   0.00000,   8.28091)  = (0 0 1) . B
    The z-axis is the 3-fold symmetry.  Looking down the z-axis, the lattice
    vectors appear roughly as:

           / b         ^ y     
      c   /            |       
    -----O             |       
          \            *----> x
           \ a

    BR1_DIR converts a vector (v) written in units of the lattice units, into the
    global cartensian coordinates (x) when dotted from the LEFT:
        x = v.B
    To do the opposite, we compute the (simple matrix-inverse) Binv = B^{-1}
        v = x.Binv
    '''
    extn = '.outputd'
    f = open(case+extn, 'r')
    lines = f.readlines()
    f.close()

    [i] = [n for n,line in enumerate(lines) if 'BR1_DIR' in line]
    return asarray([[float(x) for x in line.split()] for line in lines[i+1:i+4]])

def sort_octahedral_cage1(cage, coscut = 0.707):
    '''Given an octahedral cage of atoms, sorts them so they are ordered
    in the same manner as the vertices are defined (in the list `octahedron`).

    Tries to place closest atoms along z-axis.'''
    norms = [linalg.norm(v) for v in cage]
    unitvecs = [asarray(v)/linalg.norm(v) for v in cage]

    # choose closest atom as "top" of cage
    ntop = norms.index(min(norms))

    # find "bottom" of cage
    topangles = dot(unitvecs, unitvecs[ntop])
    [nbottom] = [n for n,angle in enumerate(topangles) if angle < -coscut]

    # the other four must lie in the xy plane
    nxy = [n for n,angle in enumerate(topangles) if (angle < coscut) and (angle > -coscut)]
    nx = norms.index(min([norms[n] for n in nxy]))  # choose closest atom in xy plane as x-axis

    # find y, -x, -y vectors by taking cross products and dotting with "top" vector
    topcrossangles = [dot(unitvecs[ntop], cross(unitvecs[nx], unitvecs[n])) for n in nxy]
    [ny] = [n for n,angle in zip(nxy, topcrossangles) if angle > coscut]
    [nminusy] = [n for n,angle in zip(nxy, topcrossangles) if angle < -coscut]
    [nminusx] = [n for n in range(len(cage)) if n not in [ntop, nbottom, nx, ny, nminusy]]

    return [cage[nx], cage[ny], cage[ntop], cage[nminusx], cage[nminusy], cage[nbottom]]
    
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

def optimize_utility(polyhedron, cage, radialweight, print_mathematica=False):
    x0 = [0, 0, 0]  # initial guess
    args = (polyhedron, cage, radialweight)
    x, f, iter, funcalls, warnflag = optimize.fmin(utility, x0, args=args, full_output=True, disp=False)

    if warnflag != 0:
        print 'WARNING: optimization failed.  Continuing.'

    if print_mathematica:
        a,b,c = x
        print 'R' + ' = {'+','.join(['{'+','.join([str(e) for e in v])+'}' for v in rotation(a,b,c)])+'};'

    return x

def newxz_from_rotation(R, case):

    #print R
    # First, extract new x and z axes in global (orthogonal) basis
    # NOTE: WIEN2K uses passive rotations, hence the transpose (R.T)
    [newx0, newy0, newz0] = R.T
    bmatinv = linalg.inv(get_bravais_matrix(case))
    # write new axis in units of lattice vectors
    newx = dot(newx0, bmatinv)
    newz = dot(newz0, bmatinv)
    return newx/norm(newx), newz/norm(newz)

def parse_cmdline_args():
    usage = """usage: %prog [options]
    Generates all possible local axes for a given atom."""

    parser = optparse.OptionParser(usage)

    parser.add_option("-a", "--atoms", default='', help="atoms in unit cell for which local rotation should be found")
    parser.add_option("-c", "--cage", default='', help="arrangement of atoms forming cage: " + ', '.join(polyhedra.keys()))
    parser.add_option("-m", "--mathematica", action='store_true', default=False, help='print output for pasting into Mathematica notebook')
    parser.add_option("-x", "--xyz", action='store_true', default=False, help='print output for pasting into .xyz file (for use in XCrySDen)')
    parser.add_option("-r", "--rotations", action='store_true', default=False, help='compute all possible rotations of local axes generated by point group')
    parser.add_option("-t", "--transposes", action='store_true', default=False, help='print transposed newx/z as well')

    # Next, parse the arguments
    (options, args) = parser.parse_args()

    if len(args) > 1:
        parser.error("Too many arguments specified (zero or one expected).")
    else:
        case = args[0] if args else W2kEnvironment().case

    if options.cage not in polyhedra.keys() and options.cage is not '':
        parser.error("Unsupported polyhedron specified for cage.")

    return case, options

def print_cage(sortedcage, R):
    print len(sortedcage)
    print 'cage'
    format = '  O ' + '%10f '*3
    for v in dot(sortedcage, R.T):
        print format % tuple(v)

polyhedra = {   #  Vertex list  Sort function         Symmetry group
    'octahedron': (octahedron,  sort_octahedral_cage, octahedral_group),
    'cube'      : (cube,        sort_cubic_cage,      octahedral_group),
    'square'    : (square,      sort_square_cage,     square_group),
    }

def main():
    case, opts = parse_cmdline_args()

    if not opts.cage:
        opts.cage = raw_input('Enter cage type (' + ', '.join(polyhedra.keys()) + '): ')

    polyhedron, sortfunc, groupfunc = polyhedra[opts.cage]

    cages = readcages(case, len(polyhedron), opts.atoms, opts.mathematica)

    if opts.xyz:
        for text, cage in cages:
            print text
            # print cage in .xyz format for xcrysden
            print len(cage)
            print 'original cage'
            format = '  O ' + '%10f '*3
            for v in cage:
                print format % tuple(v)
            print

    print
    print 'LOCAL ROTATIONS FOR SELECTED ATOMS'
    print

    for text, cage in cages:
        #print text, 'cage=', cage

        # cage must be sorted so atoms are listed in same order
        # as vertices of polyhedron
        sortedcage = sortfunc(cage)

        # optimize utility function - atoms inversely weighted by distance from central atom
        a,b,c = optimize_utility(polyhedron, sortedcage, lambda(r): 1/r)

        # compute new x and z
        R = rotation(a,b,c)
        print 'Rotation to input into case.indmfl by locrot=-1 : '
        print "%12.8f "*3 % tuple(R[0,:])
        print "%12.8f "*3 % tuple(R[1,:])
        print "%12.8f "*3 % tuple(R[2,:])
        print
        print
        newx, newz = newxz_from_rotation(R, case)
        if opts.transposes:
            newxT, newzT = newxz_from_rotation(R.T, case)
        
        if opts.mathematica:
            print
            print 'R = {'+','.join(['{'+','.join([str(e) for e in v])+'}' for v in R])+'};'
            print 'Cage = {'+','.join(['{'+','.join([str(e) for e in v])+'}' for v in cage])+'};'
            print 'SortedCage = {'+','.join(['{'+','.join([str(e) for e in v])+'}' for v in sortedcage])+'};'
            print 'newx = {'+','.join([str(e) for e in newx])+'};'
            print 'newz = {'+','.join([str(e) for e in newz])+'};'

        if opts.rotations:
            for irot,rot in enumerate(groupfunc()):
                print irot+1
                Rnew = dot(rot,R)
                # print cage in .xyz format for xcrysden
                if opts.xyz:
                    print_cage(sortedcage, Rnew)
                rot_newx, rot_newz = newxz_from_rotation(Rnew, case)
                print ('%10f '*3) % tuple(rot_newz), '  # new z-axis'
                print ('%10f '*3) % tuple(rot_newx), '  # new x-axis'
                if opts.transposes:
                    rot_newxT, rot_newzT = newxz_from_rotation(Rnew.T, case)
                    print ('%10f '*3) % tuple(rot_newzT), '  # transposed new z-axis'
                    print ('%10f '*3) % tuple(rot_newxT), '  # transposed new x-axis'
        else:
            # print cage in .xyz format for xcrysden
            if opts.xyz:
                print_cage(sortedcage, R)

            print ('%10f '*3) % tuple(newz), '  # new z-axis'
            print ('%10f '*3) % tuple(newx), '  # new x-axis'

            if opts.transposes:
                print ('%10f '*3) % tuple(newzT), '  # transposed new z-axis'
                print ('%10f '*3) % tuple(newxT), '  # transposed new x-axis'




    print


if __name__ == '__main__':
    main()

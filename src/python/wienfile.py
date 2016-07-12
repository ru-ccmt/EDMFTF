#!/usr/bin/env python


'''
Classes to handle reading of WIEN2K files.
'''


import sys, re, operator
from numpy import zeros, arange


def vectorso_exists(case):
    filename = env.SCRATCH+"/"+case+".vectorso"
    return os.path.isfile(filename) and os.path.getsize(filename) > 0


def in1c_exists(case):
    filename = case+".in1c"
    return os.path.isfile(filename) and os.path.getsize(filename) > 0


class Struct:
    def __init__(self, case, fh_info=sys.stdout):
        self.case = case
        self.extn = 'struct'
        self.parse()

    def parse(self):
        '''The file is structured for reading by fortran code,
        so data is positioned by line and by column.'''
        f = open(self.case + '.' + self.extn, 'r')
        self.title = f.next().strip()

        line = f.next()
        self.lattice = line[0:4].strip()
        self.nat = int(line[27:30])
        self.latticename = line[30:].strip()

        self.mode = f.next()[13:17]

        line = f.next()
        self.a, self.b, self.c, self.alpha, self.beta, self.gamma = [float(line[i:i+10]) for i in range(0,60,10)]

        self.iatom  = []
        self.mult   = []
        self.isplit = []
        self.aname  = []
        self.npt    = []
        self.r0     = []
        self.rmt    = []
        self.Znuc   = []
        self.pos    = []
        self.rotloc = []

        for iat in range(self.nat):
            line = f.next()
            self.iatom.append( int(line[4:8]) )
            pos = [[float(line[col:col+10]) for col in 12+13*arange(3)]]
            
            line = f.next()
            mult = int(line[15:17])
            self.mult.append(mult)

            self.isplit.append( int(line[34:37]) )

            for mu in range(mult-1):
                line = f.next()
                pos.append( [float(line[col:col+10]) for col in 12+13*arange(3)] )

            self.pos.append(pos)
                
            line = f.next()
            self.aname.append( line[0:10].strip() )
            self.npt.append( int(line[15:20]) )
            self.r0.append( float(line[25:35]) )
            self.rmt.append( float(line[40:50]) )
            self.Znuc.append( int(float(line[55:65])) )

            rt=[]
            for i in range(3):
                line = f.next()
                rt.append( [float(line[col:col+10]) for col in 20+10*arange(3)] )
            self.rotloc.append(rt)

        self.aZnuc=[]
        for iat in range(self.nat):
            for mu in range(self.mult[iat]):
                self.aZnuc.append(self.Znuc[iat])

        f.close()

    def debugprint(self, f=sys.stdout):
        print >> f, '***** Structure File Contents *****'
        print >> f, 'title =', self.title
        print >> f, 'Number of sorts, nat =', self.nat
        print >> f, 'unit cell parameters (a,b,c) =', self.a, self.b, self.c
        print >> f, 'angles (alpha, beta, gamma) =', self.alpha, self.beta, self.gamma


        printlist = [
            (self.aname,  'Atom name, aname ='),
            (self.iatom,  'Atom index, iatom ='),
            (self.mult,   'Number of atoms of this type, mult ='),
            (self.isplit, 'Symmetry of the atom, isplit ='),
            (self.npt,    'Number of radial points, npt ='),
            (self.r0,     'First point in radial mesh, r0 ='),
            (self.rmt,    'Muffin tin radius, rmt ='),
            (self.Znuc,   'Atomic number, Znuc ='),
            ]

        for i in range(self.nat):
            print >> f, '---------- atom type', i, '------------'

            for var,docstr in printlist:
                print >> f, docstr, var[i]

            print >> f, 'Position(s) inside the unit cell:'
            for m in range(self.mult[i]):
                print >> f, '    pos =', self.pos[i][m]

            print >> f, 'Local rotation matrix:'
            for row in self.rotloc[i]:
                print >> f, row

        print >> f, '***********************************'

    def flat(self, notflat):
        '''Return a flat view of given data as a list.
        Example: if w.mult = [2,4] and w.aname = ['V', 'O']
        w.flatten(w.aname) -> ['V', 'V', 'O', 'O', 'O', 'O']'''
        if notflat is self.pos:
            listoflists = self.pos
        else:
            listoflists = [[elem]*mult for elem,mult in zip(notflat, self.mult)]
        return reduce(operator.add, listoflists)


class SCF2:
    def __init__(self, case, spin=''):
        self.case = case
        self.extn = 'scf2'+spin  # spin can be '', 'up' or 'dn'
        self.parse()

    def filename(self):
        return self.case + '.' + self.extn

    def parse(self):
        f = open(self.filename(), 'r')
        lines = f.readlines()
        f.close()

        # use [-1] to take values from last iteration
        self.NOE = [float(line[38:]) for line in lines if re.match(r':NOE', line)][-1]
        self.EF  = [float(line[38:]) for line in lines if re.match(r':FER', line)][-1]

    def debugprint(self, f=sys.stdout):
        print >> f, 'From file `%s`', self.filename()
        print >> f, 'NOE =', self.NOE
        print >> f, 'EF  =', self.EF, 'Ry'

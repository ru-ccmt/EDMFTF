#!/usr/bin/env python


'''
Classes to handle reading/writing of case.indmf* files.
'''



import operator, os, re
from copy import deepcopy
from scipy import *
from numpy import array, log
import wienfile
from utils import L2str, L2num


qsplit_doc = '''    Qsplit  Description
------  ------------------------------------------------------------
     0  average GF, non-correlated
     1  |j,mj> basis, no symmetry, except time reversal (-jz=jz)
    -1  |j,mj> basis, no symmetry, not even time reversal (-jz=jz)
     2  real harmonics basis, no symmetry, except spin (up=dn)
    -2  real harmonics basis, no symmetry, not even spin (up=dn)
     3  t2g orbitals 
    -3  eg orbitals
     4  |j,mj>, only l-1/2 and l+1/2
     5  axial symmetry in real harmonics
     6  hexagonal symmetry in real harmonics
     7  cubic symmetry in real harmonics
     8  axial symmetry in real harmonics, up different than down
     9  hexagonal symmetry in real harmonics, up different than down
    10  cubic symmetry in real harmonics, up different then down
    11  |j,mj> basis, non-zero off diagonal elements
    12  real harmonics, non-zero off diagonal elements
    13  J_eff=1/2 basis for 5d ions, non-magnetic with symmetry
    14  J_eff=1/2 basis for 5d ions, no symmetry
------  ------------------------------------------------------------'''
projector_doc="""  Projector  Drscription
------  ------------------------------------------------------------
     1  projection to the solution of Dirac equation (to the head)
     2  projection to the Dirac solution, its energy derivative, 
          LO orbital, as described by P2 in PRB 81, 195107 (2010)
     4  similar to projector-2, but takes fixed number of bands in
          some energy range, even when chemical potential and 
          MT-zero moves (folows band with certain index)
     5  fixed projector, which is written to projectorw.dat. You can
        generate projectorw.dat with the tool wavef.py
------  ------------------------------------------------------------
"""

def expand_intlist(input):
    '''Expand out any ranges in user input of integer lists.
    Example: input  = "1,2,4-6"
             output = [1, 2, 4, 5, 6]'''
    def parse1(x):
        y = x.split('-')
        return [int(x)] if len(y) == 1 else range(*[int(y[0]), int(y[1])+1])
            
    return reduce(operator.add, [parse1(x) for x in input.split(',')])

def divmodulo(x,n):
    "We want to take modulo and divide in fortran way, so that it is compatible with fortran code"
    return ( sign(x)* (abs(x)/n) , sign(x)*mod(abs(x),n))

class IndmfBase:
    '''Conventions used in naming data structures stored in this class:
    i  = index
    u  = unique
    cp = correlated problem (either single-site or cluster)

    The data structures are dictionaries:
    self.atoms[iatom] = (locrot, locrot_veclist, shift_vec)
    self.cps[icp] = [(iatom_1, L_1, qsplit_1), (iatom_2, L_2, qsplit_2), ...]
    self.ucps[iucp] = [icp_1, icp_2, ...]

    Derived classes are responsible for filling in self.ucps
    '''
    def __init__(self, case):
        self.case = case
        self.extn = 'indmf'      # derived classes should override this
        self.initvars()
        self.__create_inverses()

    def member_vars(self):
        # list of tuples (varname, default value)
        # these are deepcopied when we copy_constryct()
        return [
            ('hybr_emin', -10.0 ),  # (eV) range of hybridization to pass on to impurity problem
            ('hybr_emax',  10.0  ), # (eV)
            ('Qrenorm',    1    ),  # whether or not to renormalize (what are we renormalizing?)
            ('projector',  2    ),  # type of projection onto correlated space (0,1,2,3,4)
            ('matsubara',  0    ),  # 0 = real axis, 1 = imaginary axis
            ('broadc',     0.025),  # (eV) broadening for correlated (applied to what -- self-energy?)
            ('broadnc',    0.025),  # (eV) broadening for noncorrelated
            ('om_npts',    200  ),  # number and range of default omega mesh (if no sig.inp file given)
            ('om_emin',    -3.0 ),  # (eV)
            ('om_emax',    1.0  ),  # (eV)
            ('broken_sym', 0    ),  # FM, AFM or ferrimagnetic run
            ('atoms',      {}   ),
            ('cps',        {}   ), 
            ('ucps',       {}   ),
            ('symclasses', {}   ),  # group cps forming each ucp into symmetry classes (e.g. spin-up vs. spin-down)
            ('Lsa',        []   ), 
            ('icpsa',      []   ),
            ]

    def initvars(self):
        for attr,val in self.member_vars():
            setattr(self, attr, val)

    def copy_construct(self, c):
        myattr = [attr for attr,val in self.member_vars()]
        for attr in dir(c):
            if attr in myattr:
                setattr(self, attr, deepcopy(getattr(c, attr)))

    def __create_inverses(self):
        class Iucps:
            def __getitem__(s,icp):
                return [iucp for iucp,icps in self.ucps.iteritems() if icp in icps][0]
        self.iucps = Iucps()

        class Icps:
            def __getitem__(s,iatom):
                return [icp for icp,cp in self.cps.iteritems() if iatom in [iat for iat,L,qsplit in cp]]
        self.icps = Icps()

    def filename(self):
        return self.case + '.' + self.extn

    def file_exists(self):
        return os.path.isfile(self.filename())

    def readlines(self, filename = None):
        fname = filename if filename else self.filename()
        findmf = open(fname, 'r')
        lines = [line.split('#')[0].strip() for line in findmf.readlines()] # strip comments
        findmf.close()
        return (line for line in lines if line)  # strip blank lines & create generator expression

    def parse_head(self, lines):
        self.hybr_emin, self.hybr_emax, self.Qrenorm, self.projector = [float(x) for x in lines.next().split()]
        self.matsubara, self.broadc, self.broadnc, self.om_npts, self.om_emin, self.om_emax = [float(e) for e in lines.next().split()]
        self.matsubara = int(self.matsubara)  # recast these to integers
        self.om_npts   = int(self.om_npts) 

    def parse_atomlist(self, lines):
        self.Lsa=[]
        self.icpsa=[]
        natom = int(lines.next())
        for i in range(natom):
            
            dat=lines.next().split()
            iatom, nL, locrot_shift = [int(x) for x in dat[:3]]
            
            Rmt2=0
            if len(dat)>3:
                Rmt2 = float(dat[3])
                
            #locrot = locrot_shift % 3
            #shift  = locrot_shift / 3
            #if locrot_shift<0: 
            #    locrot=3
            #    shift= (locrot_shift+3)/3

            (shift,locrot) = divmodulo(locrot_shift,3)
            if locrot<0: locrot=3
            #print 'shift=', shift, 'locrot=', locrot
            

            # Ls, qsplits, icps = array([[int(x) for x in lines.next().split()] for i in range(nL)]).T
            
            (Ls, qsplits, icps) = (zeros(nL,dtype=int), zeros(nL,dtype=int), zeros(nL,dtype=int))
            for il in range(nL):
                (Ls[il], qsplits[il], icps[il]) = map(int, lines.next().split()[:3])
                
            self.Lsa.append( Ls )
            self.icpsa.append( icps )
            
            new_zx = [[float(x) for x in lines.next().split()] for loro in range(abs(locrot))]
            vec_shift = [float(x) for x in lines.next().split()] if shift else None

            #print 'new_zx=', new_zx
            #print 'vec_shift=', vec_shift

            self.atoms[iatom] = (locrot, new_zx, vec_shift, Rmt2)
            for icp, L, qsplit in zip(icps, Ls, qsplits):
                if self.cps.has_key(icp):
                    self.cps[icp] += [(iatom, L, qsplit)]
                else:
                    self.cps[icp] = [(iatom, L, qsplit)]

    def write_head(self, lines):
        lines += [
            ("%f %f %d %d" % (self.hybr_emin, self.hybr_emax, self.Qrenorm, self.projector), "hybridization Emin and Emax, measured from FS, renormalize for interstitials, projection type"),
            ("%1d %g %g %d %f %f" % (self.matsubara, self.broadc, self.broadnc, self.om_npts, self.om_emin, self.om_emax),
             "matsubara, broadening-corr, broadening-noncorr, nomega, omega_min, omega_max (in eV)")
            ]

    def write_atomlist(self, lines):
        # create flat list of correlated orbitals (tricky because may have cluster problems)
        corbs = [[(icp,iatom,L,qsplit) for iatom,L,qsplit in v] for icp,v in self.cps.iteritems()]
        corbs = reduce(operator.add, corbs)

        # list of atom-indices of correlated atoms
        icatoms = list(set(iatom for icp,iatom,L,qsplit in corbs))
        icatoms.sort()

        lines.append((str(len(icatoms)), "number of correlated atoms"))

        for iatom in icatoms:
            locrot, locrot_veclist, shift_vec, Rmt2 = self.atoms[iatom]
            orbs = [(icp,L,qsplit) for icp,iat,L,qsplit in corbs if iat==iatom]

            locrot_shift = 3+locrot if shift_vec else locrot
            if Rmt2>0:
                atom_header = ("%-3d %3d %3d  %f" % (iatom, len(orbs), locrot_shift, Rmt2), "iatom, nL, locrot Rmt2")
            else:
                atom_header = ("%-3d %3d %3d" % (iatom, len(orbs), locrot_shift), "iatom, nL, locrot")
            lines.append(atom_header)

            for icp,L,qsplit in orbs:
                orbstring = ("%3d %3d %3d" % (L, qsplit, icp), "L, qsplit, cix")
                lines.append(orbstring)

            locrot_labels = ["new z-axis", "new x-axis"][:len(locrot_veclist)]

            for vec,label in zip(locrot_veclist, locrot_labels):
                lines.append( ("%10.7f %10.7f %10.7f" % tuple(vec), label) )

            if shift_vec:
                lines.append( ("%4s %4s %4s" % tuple(shift_vec), "real-space shift in atom position") )

    def format(self, lines):
        # merge comments with values
        comment_column = max([len(entry) for entry,comment in lines])
        format = '%-' + str(comment_column) + 's  # %s\n'
        return [format % line for line in lines]

    def writelines(self, text, filename = None):
        fname = filename if filename else self.filename()
        f = open(fname, 'w')
        f.writelines(text)
        f.close()
        return fname


class Indmf(IndmfBase):
    def __init__(self, case):
        IndmfBase.__init__(self, case)
        self.extn = 'indmf'

    def read(self):
        lines = self.readlines()
        self.parse_head(lines)

        # read ucp = {cp1, cp2, ...} arrays
        nucp = int(lines.next())
        for i in range(nucp):
            line = [int(x) for x in lines.next().split()]
            iucp = line[0]
            self.ucps[iucp] = line[1:]

        self.parse_atomlist(lines)

    def write(self):
        lines = []
        self.write_head(lines)

        # write ucp = {cp1, cp2, ...} arrays
        lines.append((str(len(self.ucps)), "number of nonequivalent correlated problems"))
        for iucp,icps in self.ucps.iteritems():
            entry = ("%3d   " % (iucp,) + ' '.join([str(icp) for icp in icps]), "iucp, cix's")
            lines.append(entry)

        self.write_atomlist(lines)

        text = self.format(lines)
        self.writelines(text)

    def user_continue(self, prompt = "Do you want to continue; or edit again? (c/e): "):
        while True:
            userin = raw_input(prompt).strip().lower()
            if userin in ['c', '', 'e']:
                break
            else:
                print 'Invalid input.'
        return userin == 'c' or userin == ''

    def orb_strings(self, cp, anames):
        '''Given list of orbitals, creates string with atomname, iatom and L for each orbital.'''
        orbstrings = []
        for iatom,L,qsplit in cp:
            orbstrings.append("%s%d %s" % (anames[iatom], iatom, L2str(L)))
        return orbstrings

    def user_input(self):
        '''Conventions used in this function:
        n = nonequivalent        cp  = correlated problem
        c = correlated           orb = orbital
        i = index (into list)

        The intermediate (temporary) data structures:
        catoms[icatom] = iatom
        corbs[icorb] = (iatom, L)
        qsplits[icorb] = qsplit

        Internal indices run from 0, user input indices run from 1.
        '''
        self.initvars()  # clear old data (if any)

        w = wienfile.Struct(self.case)     # parse WIEN2k struct file
        anames = [None] + w.flat(w.aname)  # pad to index from 1; flat list of atom names

        print "There are %d atoms in the unit cell:" % sum(w.mult)
        for i,name in enumerate(anames[1:]):
            print "%3d %s" % (i+1, name)

        while True:
            userin = raw_input("Specify correlated atoms (ex: 1-4,7,8): ")
            catoms = expand_intlist(userin)

            print "You have chosen the following atoms to be correlated:"
            for iatom in catoms:
                print "%3d %s" % (iatom, anames[iatom])

            if self.user_continue():
                break

        # currently there's no user interface to input local rotations
        for iatom in catoms:
            locrot = 0
            locrot_veclist = []
            shift_vec = []
            Rmt2=0
            self.atoms[iatom] = (locrot, locrot_veclist, shift_vec, Rmt2)

        print
        while True:
            print "For each atom, specify correlated orbital(s) (ex: d,f):"
            corbs = []
            for iatom in catoms:
                prompt = "%3d %s: " % (iatom, anames[iatom])
                userin = raw_input(prompt)
                for orb in userin.split(','):
                    entry = (iatom, L2num(orb.strip()))
                    corbs.append(entry)

            print "You have chosen to apply correlations to the following orbitals:"
            for icorb, (iatom, L) in enumerate(corbs):
                print "%3d  %s-%d %s" % (icorb+1, anames[iatom], iatom, L2str(L))

            if self.user_continue():
                break

        print
        while True:
            print "Specify qsplit for each correlated orbital (default = 0):"
            print qsplit_doc
            qsplits = []
            for icorb, (iatom, L) in enumerate(corbs):
                prompt = "%3d  %s-%d %s: " % (icorb+1, anames[iatom], iatom, L2str(L))
                userin = raw_input(prompt).strip()
                qsplit = 0 if userin == '' else int(userin)
                qsplits.append(qsplit)

            print "You have chosen the following qsplits:"
            for icorb, (iatom, L) in enumerate(corbs):
                print "%3d  %s-%d %s: %d" % (icorb+1, anames[iatom], iatom, L2str(L), qsplits[icorb])
            if self.user_continue():
                break
        
        print
        while True:
            print "Specify projector type (default = 2):"
            print projector_doc
            userin = raw_input("> ").strip()
            self.projector = 2 if userin == '' else int(userin)
            if self.projector > 4:
                import glob
                strfile = self.case+'.struct'
                enefiles = glob.glob(self.case+'.energyso')+glob.glob(self.case+'.energyso_'+'*')+glob.glob(self.case+'.energy') + glob.glob(self.case+'.energy_'+'*')
                enefiles = filter(lambda fil: os.path.getsize(fil)>0, enefiles) # Remove empty files
                if len(enefiles)==0:
                    print 'WARNING: Energy files are not present in this directory. Please generate/copy case.energy files here when using projector 5.'
                    print
            if self.user_continue():
                break
        
        print
        userin = raw_input("Do you want to group any of these orbitals into cluster-DMFT problems? (y/n): ").strip().lower()
        if userin == 'n':
            for icorb,(iatom,L) in enumerate(corbs):
                icp = icorb+1
                self.cps[icp] = [(iatom, L, qsplits[icorb])]
        else:
            print
            while True:
                print "Enter the orbitals forming each cluster-DMFT problem, separated by spaces"
                userin = raw_input("(ex: 1,2 3,4 5-8): ")
                expanded = [expand_intlist(group) for group in userin.split()]
                expandedflat = reduce(operator.add, expanded)

                # add orbitals not in CDMFT problems
                icp = 1
                for icorb,(iatom,L) in enumerate(corbs):
                    if icorb+1 not in expandedflat:
                        self.cps[icp] = [(iatom, L, qsplits[icorb])]
                        icp += 1

                # then add orbitals that are part of CDMFT problems
                for group in expanded:
                    self.cps[icp] = []
                    for icorb in group:
                        iatom, L = corbs[icorb-1]
                        self.cps[icp] += [(iatom, L, qsplits[icorb-1])]
                    icp += 1

                print "Your choices give the following correlated problems:"
                for icp,cp in self.cps.iteritems():
                    orbstrings = self.orb_strings(cp, anames)
                    print "%2d  (%s)" % (icp, ', '.join(orbstrings))

                if self.user_continue():
                    break

        print
        while True:
            print "Enter the correlated problems forming each unique correlated"
            userin = raw_input("problem, separated by spaces (ex: 1,3 2,4 5-8): ")
            for i,group in enumerate(userin.split()):
                self.ucps[i+1] = expand_intlist(group)
            print
            print "Each set of equivalent correlated problems are listed below:"
            for iucp,ucp in self.ucps.iteritems():
                cpstrings = ['(%s)' % ', '.join(self.orb_strings(self.cps[icp], anames)) for icp in ucp]
                print "%3d   %s are equivalent." % (iucp, ' '.join(cpstrings))
            if self.user_continue():
                break
            self.ucps = {}  # reset

        print
        #userin = raw_input("Broken symmetry run? (y/n): ").strip().lower()
        #if userin == 'y':
        #    self.broken_sym = True
        #    userin = int(raw_input("What type of broken symmetry (1 = FM, 2 = AFM, 3 = spiral, ferrimagnetic, etc.)?: ").strip())
        #    if userin == 1:
        #        # FM run
        #        pass
        #    elif userin == 2:
        #        # AFM
        #        print "Not FM, so must be AFM, spiral or ferrimagnetic run."
        #        while True:
        #            for iucp,ucp in self.ucps.iteritems():
        #                print "For unique correlated problem %d containing:" % iucp
        #                cpstrings = ['%s' % ', '.join(self.orb_strings(self.cps[icp], anames)) for icp in ucp]
        #                # print out enumerated list of orbitals forming ucp
        #                print "   ", '\n   '.join(cpstrings)
        #                userin = raw_input("Group correlated orbitals into symmetry classes, separated by spaces (ex: 1,3 2,4 5-8): ").strip().lower()
        #                self.symclasses[iucp] = expand_intlist(userin)
        #    elif userin == 3:
        #        # spiral or ferrimagnetic run
        #        pass
        #    else:
        #        # bad user input
        #        pass

        print
        print "Range (in eV) of hybridization taken into account in impurity"
        userin = raw_input("problems; default %.1f, %.1f: " % (self.hybr_emin, self.hybr_emax))
        if userin.strip():
            self.hybr_emin, self.hybr_emax = [float(e) for e in userin.split(',')]

        print
        userin = raw_input("Perform calculation on real; or imaginary axis? (r/i): ").strip().lower()
        self.matsubara = 1 if userin == 'i' else 0



class Indmfl(IndmfBase):
    '''Class for case.indmfl file.

    Additional member variables/data structures:
    self.siginds[icp] = sigind
    self.cftrans[icp] = cftrans
    self.legends[icp] = legends
    EF = fermi level in eV
    '''
    def __init__(self, case, extn='indmfl'):
        IndmfBase.__init__(self, case)
        self.extn = extn #'indmfl'

        # Finding the chemical potential
        EF_exists = os.path.isfile('EF.dat')
        scf2_exists = os.path.isfile(case+".scf2")
        scf2up_exists = os.path.isfile(case+".scf2up")
        scf_exists = os.path.isfile(case+".scf")
        self.EF = None
        
        if EF_exists:
            # The previous DMFT chemical potential
            self.EF = float( open('EF.dat','r').readline() )
            
        if self.EF is None and (scf2_exists or scf2up_exists):
            fname = case+".scf2" if scf2_exists else case+".scf2up"
            fscf = open(fname, 'r')
            lines = fscf.readlines()
            for line in lines:
                if re.match(r':FER', line) is not None:
                    Ry2eV = 13.60569193
                    self.EF = float(line[38:])*Ry2eV
                    break
                
        if self.EF is None and scf_exists:
            fname = case+".scf"
            fscf = open(fname, 'r')
            lines = fscf.readlines()
            for line in lines:
                if re.match(r':FER', line) is not None:
                    Ry2eV = 13.60569193
                    self.EF = float(line[38:])*Ry2eV
                    
        if self.EF is None: self.EF = 0
        
        
    def member_vars(self):
        myvars = [
            ('siginds', {} ),
            ('cftrans', {} ),
            ('legends', {} ),
            ]
        return IndmfBase.member_vars(self) + myvars

    def read(self, filename = None):
        lines = self.readlines(filename)
        self.parse_head(lines)
        self.parse_atomlist(lines)

        # read the big block of siginds and cftrans
        ncp, maxdim, maxsize = [int(e) for e in lines.next().split()]
        for i in range(ncp):
            icp, dim, size = [int(e) for e in lines.next().split()]
            self.legends[icp] = lines.next().split("'")[1::2]
            self.siginds[icp] = array([[int(e) for e in lines.next().split()] for row in range(dim)])
            raw_cftrans = array([[float(e) for e in lines.next().split()] for row in range(dim)])
            self.cftrans[icp] = raw_cftrans[:,0::2] + raw_cftrans[:,1::2]*1j

    def write_head(self, lines):

        if abs(self.projector)<4: # This is the old scheme, where hybridization is cut-off by energy
            emin = self.hybr_emin+self.EF
            emax = self.hybr_emax+self.EF
            styp="%f "
            sdoc = "hybridization Emin and Emax, measured from FS, renormalize for interstitials, projection type"
        else:    # In the new scheme, we cut-off at certain band index
            import findNbands
            import glob
            import sys
            strfile = self.case+'.struct'
            enefiles = glob.glob(self.case+'.energyso')+glob.glob(self.case+'.energyso_'+'*')
            if not enefiles:  # Not spin-orbit run
                enefiles = glob.glob(self.case+'.energy') + glob.glob(self.case+'.energy_'+'*')
            enefiles = filter(lambda fil: os.path.getsize(fil)>0, enefiles) # Remove empty files
            
            print 'all enefiles=', enefiles
            
            #for fil in enefiles:
            #    if re.match(self.case+'.energyso', fil): # Spin-orbit on, remove non-spin-orbit files
            #        enefiles = filter(lambda fil: re.match(self.case+'.energyso', fil) is not None, enefiles) # Remove empty files
            #        break
            
            if len(enefiles)==0:
                print "ERROR : The case.energy* files should be present in this directory when using projector 5. Exiting...."
                sys.exit(1)
            
            (nemin,nemax) = findNbands.findNbands(self.hybr_emin+self.EF,self.hybr_emax+self.EF,enefiles,strfile)
            emin,emax = nemin,nemax
            styp="%d "
            sdoc = "hybridization band index nemin and nemax, renormalize for interstitials, projection type"
        
        if abs(self.projector)==5:
            import wavef
            Rm2=[self.atoms[iatom][3] for iatom in self.atoms.keys()] 

            print 'Rm2=', Rm2, 'atms=', self.atoms.keys(), 'Lsa=', self.Lsa
            wavef.main(self.case, self.atoms.keys(), self.Lsa, self.icpsa, Rm2)
            
        lines += [
            ( (styp+styp+"%d %d") % (emin, emax, self.Qrenorm, self.projector), sdoc),
            ("%1d %g %g %d %f %f" % (self.matsubara, self.broadc, self.broadnc, self.om_npts, self.om_emin, self.om_emax),
             "matsubara, broadening-corr, broadening-noncorr, nomega, omega_min, omega_max (in eV)")
            ]

            

    def write(self, filename = None):
        # generate text in two chunks, stored in text and text2
        #   text contains basic information about correlated problems
        #   text2 contains all the siginds, legends and crystal-field transformation matrices
        lines = []
        self.write_head(lines)
        self.write_atomlist(lines)
        text = self.format(lines)

        maxdim = max(len(s) for s in self.siginds.values()) # dimension of largest sigind matrix
        sizes = {}
        for icp,sigind in self.siginds.iteritems():
            sizes[icp] = len([x for x in set(sigind.flat) if x != 0])
        maxsize = max(sizes.values())  # number of columns in largest

        text2 = [
            '#================ # Siginds and crystal-field transformations for correlated orbitals ================',
            '%-3d %3d %3d       # Number of independent kcix blocks, max dimension, max num-independent-components' % (len(self.cps), maxdim, maxsize)
            ]

        for icp,sigind in self.siginds.iteritems():
            legend = self.legends[icp]
            cftrans = self.cftrans[icp]

            text2 += [
                "%-3d %3d %3d       # %s" % (icp, len(sigind), sizes[icp], 'cix-num, dimension, num-independent-components'),
                '#---------------- # Independent components are --------------',
                "'%s' "*len(legend) % tuple(legend),
                ]

            text2.append('#---------------- # Sigind follows --------------------------')
            max_sigfig = 1 + int(log(max(sigind.flat))/log(10))
            format = '%' + str(max_sigfig) + 'd'
            for row in sigind:
                text2.append(' '.join([format % elem for elem in row]))

            # print local transform matrix (real & imag)
            text2.append('#---------------- # Transformation matrix follows -----------')
            for row in self.cftrans[icp]:
                text2.append(' '.join(["%11.8f %-11.8f" % (elem.real, elem.imag) for elem in row]))

        # join with first half; add \n to each line in text2
        text += [line+'\n' for line in text2]
        self.writelines(text, filename)
        

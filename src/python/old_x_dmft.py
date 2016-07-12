#!/usr/bin/env python
import sys, os, glob, re, numpy, shutil, subprocess
from scipy import *
from scipy import optimize, interpolate
import copy,optparse
# custom imports
import utils,indmffile
import chemicalP as c_p

def FindChemicalPotential(case, updn):
    # Looking for the LDA chemical potential
    Ry2eV = 13.60569193
    EF_found = False

    if os.path.isfile(case+'.in2'):
        fname = case+".in2"
    elif os.path.isfile(case+'.in2c'):
        fname = case+".in2c"
    else:
        raise Exception("Failed to determine number of electrons. You should have case.in2 or case.in2c file!")
        
    lines = open(fname, 'r').readlines()
    NOE = float(lines[1].split()[1])
    print 'NOE=', NOE
        
    fname = case+".scf2"
    if os.path.isfile(fname):
        fscf = open(fname, 'r')
        print "opening", fname
    elif os.path.isfile(fname+updn):
        fscf = open(fname+updn, 'r')
        print "opening", fname+updn
    else:
        fscf = None

    if fscf is not None:
        lines = fscf.readlines()
        for line in lines:
            if re.match(r':FER', line) is not None:
                EF = float(line[38:])*Ry2eV
                print 'Ef=', EF
                EF_found = True

    # The previous DMFT chemical potential
    if  os.path.isfile('EF.dat'):
        fmu = open('EF.dat','r')
        mu = float(fmu.next())
        print 'Found DMFT-mu=', mu
        EF = mu
        EF_found = True

    if not EF_found:
        fname = case+".scf"
        if os.path.isfile(fname):
            print "opening", fname
            f = open(fname, 'r')
            lines = f.readlines()
            if not EF_found:
                for line in lines:
                    if re.match(r':FER', line) is not None:
                        EF = float(line[38:])*Ry2eV
                        print 'Ef=', EF
                        EF_found = True

    if not EF_found:
        raise Exception("Failed to determine chemical potential.")

    return (EF, NOE)

def CheckMatsubara(filename):
    finq =open(filename, 'r')
    finq.next()
    (matsubara, gamma, gammac, nom_def, aom_def, bom_def) = finq.next().split()[:6]
    return int(matsubara)

def PrepareDefinitionFile(idmf, mode, case, cixs, updn, dnup, so, scratch, cmplx, _band='', m_ext=''):

    fdef = open('dmft'+idmf+'.def', 'w')
    
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 2, "'"+case+".indmf"+idmf+"'",         "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 4, "'"+case+".inso"+"'",               "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 5, "'"+case+".indmfl"+m_ext+"'",             "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 6, "'"+case+".outputdmf"+idmf+updn+"'","'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 7, "'"+case+".in1c"+"'",               "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 8, "'"+case+".scf2"+updn+"'",          "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 9, "'"+scratch+"/"+case+".vector"+so+updn+"'", "'unknown'","'unformatted'",9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (10, "'"+scratch+"/"+case+".vector"+so+dnup+"'", "'unknown'","'unformatted'",9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (14, "'"+case+".klist"+_band+"'",        "'old'",    "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (18, "'"+case+".vsp"+updn+"'",           "'old'",    "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (19, "'"+case+".vsp"+dnup+"'",           "'unknown'","'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (20, "'"+case+".struct"+"'",             "'old'",    "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (22, "'"+case+".rotlm"+"'",              "'unknown'","'formatted'",0)

    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (59, "'"+case+".energy"+so+dnup+"'", "'unknown'","'formatted'",0) #
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (60, "'"+case+".energy"+so+updn+"'", "'unknown'","'formatted'",0) # 

    if mode=='g':
        print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (100,"'"+case+".cdos"+m_ext+"'",  "'unknown'","'formatted'",0) # total DOS written here

    # Appending to the def file for correlated atoms
    # for each correlated (atom,l)-block one file
    for c in cixs:
        if mode=='g':
            print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (120+c,"'"+case+'.gc'+str(c)+m_ext+"'",  "'unknown'","'formatted'",0) # only correlated green's functions
            print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (140+c,"'"+case+'.dlt'+str(c)+m_ext+"'", "'unknown'","'formatted'",0) # delta's for correlated
        
        print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (80+c, "'"+'sig.inp'+str(c)+m_ext+_band+"'",        "'old'",    "'formatted'",0) # self-energy for correlated

        if mode=='g':
            print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (180+c,"'"+case+'.Eimp'+str(c)+m_ext+"'", "'unknown'","'formatted'",0) # Impurity levels for correlated
            
    if mode in ['e', 'p']:
        print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (180,"'eigvals"+m_ext+".dat'", "'unknown'","'formatted'",0) # eigenvalues will be stored

    if mode=='u':
        print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (399,"'Udmft"+m_ext+updn+".0'", "'unknown'","'unformatted'",0) # eigenvalues will be stored

    fdef.close()


def PrepareDefinitionFile2(idmf, mode, case, cixs, updn, dnup, so, scratch, cmplx='', m_ext=''):

    fdef = open('dmft'+idmf+'.def', 'w')
    if so=='so':
        sodum = 'dum'
    else:
        sodum=dnup

    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 2, "'"+case+".indmf"+idmf+"'",    "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 3, "'"+case+".in1c"+"'",          "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 4, "'"+case+".inso"+"'",          "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 5, "'"+case+".in2"+cmplx+"'",     "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 6, "'"+case+".outputdmf"+idmf+updn+"'","'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 7, "'"+case+'.indmfl'+m_ext+"'",        "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 8, "'"+case+".clmval"+updn+m_ext+"'",   "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 9, "'"+scratch+"/"+case+".vector"+so+updn+"'", "'unknown'","'unformatted'",9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (10, "'"+scratch+"/"+case+".vector"+so+dnup+"'", "'unknown'","'unformatted'",9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (13, "'"+case+".recprlist"+"'",     "'unknown'", "'formatted'", 9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (14, "'"+case+".kgen"+"'",          "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (18, "'"+case+".vsp"+updn+"'",      "'old'",     "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (20, "'"+case+".struct"+"'",        "'old'",     "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (21, "'"+case+".scf2"+updn+m_ext+"'",     "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (22, "'"+case+".rotlm"+"'",         "'unknown'", "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (29, "'"+case+".energy"+sodum+"'",  "'unknown'", "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (30, "'"+case+".energy"+so+updn+"'","'unknown'", "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (31, "'"+case+".energy"+so+dnup+"'","'unknown'", "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (12, "'"+case+".norm"+so+"'",       "'unknown'", "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (100, "'"+case+".cdos3"+m_ext+"'",           "'unknown'", "'formatted'",0)
    
    # Appending to the def file for correlated atoms
    # for each correlated (atom,l)-block one file
    for c in cixs:
        print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (80+c, "'"+'sig.inp'+str(c)+m_ext+"'","'old'","'formatted'",0) # self-energy for correlated
    fdef.close()

def PrepareInFile(idmf, mode, case, updn, nsymop=0, recomputeEF=None, mixEF=1.0, WL=None):

    fdin = open(case+'.indmf'+idmf, 'w')
    ## Looking for the chemical potential
    (EF, NOE) = FindChemicalPotential(case, updn)

    print >> fdin, '#-------------- mode and the current chemical potential --------------------'
    if recomputeEF is not None:
        print >> fdin, mode, EF, recomputeEF, mixEF, '        #  mode, EF, recomputeEF, mixing-EF'
        if (recomputeEF>1):
            print >> fdin, WL, '       # frequency range where to look for a pole'
    else:
        print >> fdin, mode, EF, nsymop, '        #  mode, EF, nsymop'
    
    fdin.close()

def checkSigmas(cixs, sig='sig.inp'):
    """ Checks that all self-energies exist and have enough columns
    """
    for c in cixs:
        sigfile = sig+str(c)
        fs = open(sigfile, 'r')
        data = fs.next().split()
        m = re.match(r'\s*#',data[0])
        if m is not None:            # We found comment in first line
            data = fs.next().split() # Need to read next line
        if len(data)<3: #1+2*SigSize[c]:
            print 'ERROR: Self-energy file', sigfile, 'does not have enough columns!'
            print 'ERROR: Required number of independent components is', SigSize[c], 'which requires', 1+2*SigSize[c], 'columns. Found', len(data)
            print 'ERROR: You can continue to run the code, but the ouput is likely wrong!'



def runExternal(idmf, ROOT, MPI, nom, ntail, dmft_exe='dmft', m_ext=''):
    
    cmd = ROOT+'/ssplit.py -n '+str(nom)+' -t '+str(ntail)
    if len(m_ext)>0: cmd += ' -l '+m_ext

    print '..... running:', cmd
    stdin, stdout, stderr = os.popen3( cmd )
    print stderr.read().strip(), stdout.read().strip()

    deffile = 'dmft'+idmf+m_ext+'.def'
    
    if MPI: # In parallel mode needs to copy the executable to the current directory
        shutil.copy(ROOT+'/'+dmft_exe, os.getcwd() )
        exe = MPI + ' ./'+dmft_exe
    else:
        exe = dmfe.ROOT+'/'+dmft_exe

    cmd = exe + ' ' + deffile

    print '..... running:', cmd
    subprocess.call(cmd, bufsize = 1, shell=True)  # line-buffered

    # OLD VERSION: annoyingly fully buffered, so no status on dmft fortran executable
    # until dmft exits.
    #
    # sout, serr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()
    # print sout, 
    # if serr is not None:
    #     print serr,
    
    #stdin, stdout, stderr = os.popen3( cmd )
    #print stdout.read().strip(), stderr.read().strip()
    

if __name__=='__main__':
    
    gammamu = 1e-14     # cmall broadening for chemical potential
    sdmu=0.1            # steps when looking for the chemical potential
    mix_mu = 1.0        # mixing of the chemical potential
    com = 0             # when computing density on imaginary axis, we can subtract Sigma_oo (com=0) or Sigma(0) (com=1)
    
    usage = """usage: %prog [ options ] mode

    Executes one of the dmft steps (dmft0|dmft1|dmft2|mu|dmftp|dmftu):
      dmft0  -- computes all complex frequency dependent eigenvalues. Stores them in 'eigenvalues.dat'
      dmft1  -- computes local green's function (case.gc?), hybridization (case.dlt?) and impurity levels (case.Eimp?)
      dmft2  -- computes the valence LDA+DMFT charge
      mu     -- using eigenvalues computed by dmft0 recomptes the chemical potential
      dmftp  -- prepares the data for the plot of A(k,omega)
      dmftu  -- prints out the DMFT transformation needed by other programs, such as optics.
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("--so", dest="so", action='store_true', default=False, help="switches on spin-orbit coupling")
    parser.add_option("--up", dest="up", action='store_true', default=False, help="magnetic LDA calculation with vector-up first")
    parser.add_option("--dn", dest="dn", action='store_true', default=False, help="magnetic LDA calculation with vector-dn first")
    parser.add_option("-c",  dest="c",  action='store_true', default=False, help="complex version")
    parser.add_option("-m",  dest="recomputeEF",  type="int", default=1, help="Recompute EF in dmft2 step [0|1|2]")
    parser.add_option("-d", dest="def_only", action='store_true', default=False, help="prepared only def and input files, but do not run")
    parser.add_option("-n", "--nom",  dest="nom",    type="int", default=100, help="Number of matsubara points in log mesh before the tail")
    parser.add_option("-t", "--tail", dest="ntail",  type="int", default=30,  help="Number of matsubara points in the tail of the log mesh")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="print verbose messages to stdout")
    parser.add_option("-l", "--lext", dest="m_ext",  default='',  help="For magnetic calculation, it can be 'dn'.")
    parser.add_option("-x", "--mixEF", dest="mixEF",  default=1.0,  help="Chemical potential can be mixed with mixEF<1.0")
    parser.add_option("-w", "--WL", dest="WL",  default=1.0,  help="If recomputeEF=2, pole is searched in the interval [-WL,WL]")
    parser.add_option("-b", "--band", action="store_true", dest="_band", default=False, help="use case.klist_band instead of case.klist")
    parser.add_option("-g", "--no_group", action="store_true", dest="_no_group", default=False, help="use case.klist_band instead of case.klist")
    # Next, parse the arguments
    (options, args) = parser.parse_args()

    if len(args)!=1:
        print 'Need exactly one argument. One of the dmft0|dmft1|dmft2|mu|dmftp|dmftu'
        sys.exit(1)
    
    dmfe = utils.DmftEnvironment()  # DMFT paths
    w2k = utils.W2kEnvironment()    # W2k filenames and paths
    
    if not options.so and os.path.isfile(w2k.case+".inso") and os.path.getsize(w2k.case+".inso")>0 :
        print 'Found '+w2k.case+'.vectorso file, hence assuming so-coupling exists. Switching -so switch!'
        options.so = True

    if not options.c and os.path.isfile(w2k.case+".in1c") and os.path.getsize(w2k.case+".in1c")>0 :
        options.c = True

    # processing arguments
    updn=''
    dnup='dn'
    so=''
    cmplx=''
    
    if options.so:
        so = 'so'
        sodum = 'dum'
        cmplx = 'c'
    if options.up:
        spin=2
        updn = 'up'
        dnup = 'dn'
    if options.dn:
        spin=2
        updn = 'dn'
        dnup = 'up'
    if options.c:
        cmplx = 'c'
        
    # Processing 'case.indmfl' file
    inl = indmffile.Indmfl(w2k.case) # case.indmfl file
    inl.read()                       # case.indmfl read
    cixs = inl.siginds.keys()        # all columns for cix

    #print 'name=', args[0]
    print 'case=', w2k.case
    #print 'c=', cmplx
    #print 'scratch=', w2k.SCRATCH
    print 'root=', dmfe.ROOT
    print 'cixs=', cixs
    #print 'def_only=', options.def_only
    
    if args[0] == 'dmft0':
        idmf = '0'
        mode = 'e'               # mode for computing eigenvalues
        
        PrepareDefinitionFile(idmf, mode, w2k.case, cixs, updn, dnup, so, w2k.SCRATCH, cmplx)
        PrepareInFile(idmf, mode, w2k.case, updn, 1)

        if not options.def_only:
            #checkSigmas(cixs)
            runExternal(idmf, dmfe.ROOT, dmfe.MPI, options.nom, options.ntail, 'dmft')
            
    elif args[0] == 'dmft1':
        idmf = '1'
        mode = 'g'               # mode for computing eigenvalues
        
        PrepareDefinitionFile(idmf+options.m_ext, mode, w2k.case, cixs, updn, dnup, so, w2k.SCRATCH, cmplx, '', options.m_ext)
        PrepareInFile(idmf+options.m_ext, mode, w2k.case, updn)
        
        if not options.def_only:
            #checkSigmas(cixs)
            runExternal(idmf, dmfe.ROOT, dmfe.MPI, options.nom, options.ntail, 'dmft', options.m_ext)
            
    elif args[0] == 'dmft2':
        idmf = '2'
        mode = 'c'

        PrepareDefinitionFile2(idmf+options.m_ext, mode, w2k.case, cixs, updn, dnup, so, w2k.SCRATCH, cmplx, options.m_ext)
        PrepareInFile(idmf+options.m_ext, mode, w2k.case, updn, 0, options.recomputeEF, options.mixEF, options.WL)

        if not options.def_only:
            #checkSigmas(cixs)
            runExternal(idmf, dmfe.ROOT, dmfe.MPI, options.nom, options.ntail, 'dmft2',options.m_ext)
        
    elif args[0] == 'dmftu':
        idmf = 'u'
        mode = 'u'

	k_band = ''
        nsymop=0
	if (options._band): 
            k_band = '_band'
            nsymop = 1
	if (options._no_group):
            nsymop = 1
            
        PrepareDefinitionFile(idmf, mode, w2k.case, cixs, updn, dnup, so, w2k.SCRATCH, cmplx, k_band)
        PrepareInFile(idmf, mode, w2k.case, updn, nsymop)
        if not options.def_only:
            runExternal(idmf, dmfe.ROOT, dmfe.MPI, options.nom, options.ntail, 'dmft')
        
    elif args[0] == 'mu':
        # See if working on imaginary or real axis
        matsubara = CheckMatsubara(w2k.case+'.indmfl')
        
        # The previous DMFT chemical potential
        (EF, NOE) = FindChemicalPotential(w2k.case, updn)

        if (not options.so) and (not options.up) and (not options.dn): NOE/=2. # spin-up and spin-down equivalent. Only 1/2 of all bands is stored and computed
        
        wupdn=['']
        if updn: wupdn = ['-up', '-dn']
        
        valC=[]
        for iupd,upd in enumerate(wupdn):
            # reading eigenvectors
            (Ek, omega, wkp, nbands, nemin) = c_p.ReadEigenvals('eigvals.dat'+upd, matsubara, gammamu)
            
            # find roughly where the chemical potential is
            valC.append( c_p.ValCharge(matsubara, wkp, omega, Ek, nemin, nbands, com) )
            
        (mu0, mu1) = c_p.FindSignChange(c_p.ChargeDiff, EF, sdmu, sys.stdout, args=(valC, NOE))
        
        # find precise chemical potential
        mu_new = optimize.brentq(c_p.ChargeDiff, mu0, mu1, args=(valC, NOE))
        
        # difference in mu
        dmu = mu_new - EF
        # mixing for the chemical potential
        mu = EF*(1-mix_mu) + mix_mu*mu_new
        print 'New chemical potential found at ', dmu, ' -> Chemical potential becomes ', mu
        fc = open('EF.dat','w')
        print >> fc, mu_new
        fc.close()
        
    elif args[0] == 'dmftp':
        idmf = 'p'
        mode = 'p'               # mode for computing eigenvalues

        PrepareInFile(idmf, mode, w2k.case, updn,1)
        PrepareDefinitionFile(idmf, mode, w2k.case, cixs, updn, dnup, so, w2k.SCRATCH, cmplx, '_band')

        nom = inl.om_npts 
        aom = inl.om_emin  
        bom = inl.om_emax  
 
        eq_om = linspace(inl.om_emin, inl.om_emax, inl.om_npts)  # Equidistant mesh

        # Interpolates sigmas on equidistant mesh!
        for c in cixs:
            fsig = open('sig.inp'+str(c), 'r')
            data = fsig.readlines()
            fsig.close()
            om0 =[]
            sig0=[]
            for line in data:
                dd = map(float,line.split())
                om0.append( dd[0] )
                sig0.append( dd[1:] )
            sig0=array(sig0).transpose()

            sig1=[]
            for i in range(len(sig0)):
                tck = interpolate.splrep(om0, sig0[i])
                sig1.append( interpolate.splev(eq_om, tck) )
            sig1 = array(sig1).transpose()
            
            osig = open('sig.inp'+str(c)+'_band', 'w')
            for i in range(len(eq_om)):
                print >> osig, '%15.8f  ' % eq_om[i], ('%15.8f '*len(sig1[i])) % tuple(sig1[i])
            osig.close()

        if not options.def_only:
            runExternal(idmf, dmfe.ROOT, dmfe.MPI, options.nom, options.ntail, 'dmft')
    else:
        print 'not yet implemented!'

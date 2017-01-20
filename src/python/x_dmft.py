#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import sys, os, glob, re, numpy, shutil, subprocess
from scipy import *
from scipy import optimize, interpolate
import copy,optparse
# custom imports
import utils,indmffile
import chemicalP as c_p

def getTemperature():
    if os.path.exists('params.dat'):
        execfile('params.dat')
        # get locals variables from params.dat
        vglobl={}
        vlocal={}
        execfile('params.dat',vglobl, vlocal)
        iparams0 = vlocal['iparams0']
        if iparams0.has_key('T'):
            Temperature = iparams0['T'][0]
        elif iparams0.has_key('beta'):
            Temperature = 1./iparams0['beta'][0]
        else:
            Temperature = 1e-10
    else:
        fi = open('sig.inp')
        fi.next() # s_oo
        fi.next() # Edc
        om0 = float(fi.next().split()[0])
        Temperature=om0/pi
        print 'Temperature set to ', Temperature
    return Temperature
    
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

def PrepareDefinitionFile_dmft1(idmf, mode, case, cixs, updn, dnup, so, para, scratch, cmplx, _band='', m_ext=''):

    IND_CDOS = 500
    IND_UDMF = 501
    IND_S = 80
    IND_G = 180
    IND_D = 280   
    IND_IMP = 380

    fdef = open('dmft'+idmf+'.def', 'w')
    idmf_updn = idmf+updn
    if idmf[-2:]==updn:
        idmf_updn = idmf

    print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % ( 2, "'"+case+".indmf"+idmf+"'",         "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % ( 4, "'"+case+".inso"+"'",               "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % ( 5, "'"+case+".indmfl"+m_ext+"'",             "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % ( 6, "'"+case+".outputdmf"+idmf_updn+"'","'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % ( 7, "'"+case+".in1c"+"'",               "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % ( 8, "'"+case+".scf2"+updn+"'",          "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % ( 9, "'"+scratch+"/"+case+".vector"+so+updn+para+"'", "'old'","'unformatted'",9000)
    print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % (10, "'"+scratch+"/"+case+".vector"+so+dnup+para+"'", "'unknown'","'unformatted'",9000)
    print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % (14, "'"+case+".klist"+_band+"'",        "'old'",    "'formatted'",0)
    print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % (18, "'"+case+".vsp"+updn+"'",           "'old'",    "'formatted'",0)
    print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % (19, "'"+case+".vsp"+dnup+"'",           "'unknown'","'formatted'",0)
    print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % (20, "'"+case+".struct"+"'",             "'old'",    "'formatted'",0)
    print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % (21, "'"+case+".scfdmft1"+m_ext+"'",           "'unknown'","'formatted'",0)
    print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % (22, "'"+case+".rotlm"+"'",              "'unknown'","'formatted'",0)

    if mode=='g':
        print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % (IND_CDOS,"'"+case+".cdos"+m_ext+"'",  "'unknown'","'formatted'",0) # total DOS written here

    print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % (59, "'"+case+".energy"+so+dnup+para+"'", "'unknown'","'formatted'",0) #
    print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % (60, "'"+case+".energy"+so+updn+para+"'", "'old'","'formatted'",0) # 

    # Appending to the def file for correlated atoms
    # for each correlated (atom,l)-block one file
    for c in cixs:
        print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % (IND_S+c, "'"+'sig.inp'+str(c)+m_ext+_band+"'",        "'old'",    "'formatted'",0) # self-energy for correlated

        if mode=='g':
            print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % (IND_G+c,"'"+case+'.gc'+str(c)+m_ext+"'",  "'unknown'","'formatted'",0) # only correlated green's functions
            print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % (IND_D+c,"'"+case+'.dlt'+str(c)+m_ext+"'", "'unknown'","'formatted'",0) # delta's for correlated
            print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % (IND_IMP+c,"'"+case+'.Eimp'+str(c)+m_ext+"'", "'unknown'","'formatted'",0) # Impurity levels for correlated
            
    if mode in ['e', 'p']:
        print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % (IND_IMP,"'eigvals"+m_ext+".dat'", "'unknown'","'formatted'",0) # eigenvalues will be stored

    if mode=='u':
        print >> fdef, "%3d, %-20s, %-10s, %-13s, %-4d" % (IND_UDMF,"'Udmft"+m_ext+updn+".0'", "'unknown'","'unformatted'",0) # eigenvalues will be stored

    fdef.close()


def PrepareDefinitionFile_dmft2(idmf, mode, case, cixs, updn, dnup, so, para, scratch, cmplx='', m_ext=''):

    IND_CDOS = 500
    IND_S = 80

    fdef = open('dmft'+idmf+'.def', 'w')
    if so=='so':
        sodum = 'dum'
    else:
        sodum=dnup

    updn_m_ext = updn+m_ext
    if m_ext and updn and m_ext==updn:  # if both are set to ither up or dn
        updn_m_ext = m_ext              # we should not repeat
    idmf_updn = idmf+updn
    if idmf[-2:]==updn:
        idmf_updn = idmf
        
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 2, "'"+case+".indmf"+idmf+"'",    "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 3, "'"+case+".in1c"+"'",          "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 4, "'"+case+".inso"+"'",          "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 5, "'"+case+".in2"+cmplx+"'",     "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 6, "'"+case+".outputdmf"+idmf_updn+"'","'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 7, "'"+case+'.indmfl'+m_ext+"'",        "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 8, "'"+case+".clmval"+updn_m_ext+"'",   "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 9, "'"+scratch+"/"+case+".vector"+so+updn+para+"'", "'old'","'unformatted'",9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (10, "'"+scratch+"/"+case+".vector"+so+dnup+para+"'", "'unknown'","'unformatted'",9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (12, "'"+case+".norm"+so+"'",       "'unknown'", "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (13, "'"+case+".recprlist"+"'",     "'unknown'", "'formatted'", 9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (14, "'"+case+".kgen"+"'",          "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (15, "'"+scratch+"/"+case+".bnds"+m_ext+"'",    "'unknown'","'unformatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (18, "'"+case+".vsp"+updn+"'",      "'old'",     "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (20, "'"+case+".struct"+"'",        "'old'",     "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (21, "'"+case+".scf2"+updn_m_ext+"'",     "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (22, "'"+case+".rotlm"+"'",         "'unknown'", "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (IND_CDOS, "'"+case+".cdos3"+m_ext+"'",           "'unknown'", "'formatted'",0)    
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (29, "'"+case+".energy"+sodum+"'",  "'unknown'", "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (30, "'"+case+".energy"+so+updn+para+"'","'old'", "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (31, "'"+case+".energy"+so+dnup+para+"'","'unknown'", "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (9902, "'"+case+".nsh"+updn+"'","'unknown'", "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (9919, "'"+case+".vns"+updn+"'","'unknown'", "'formatted'",0)
    
    # Appending to the def file for correlated atoms
    # for each correlated (atom,l)-block one file
    for c in cixs:
        print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (IND_S+c, "'"+'sig.inp'+str(c)+m_ext+"'","'old'","'formatted'",0) # self-energy for correlated
    fdef.close()

def PrepareDefinitionFile_lapw1(widmf, case, updn, para, scratch, band, cmplx):
    fdef = open(widmf+updn+'.def', 'w')
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 4, "'"+case+".klist"+band+"'",              "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 5, "'"+case+".in1"+cmplx+"'",               "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 6, "'"+case+".output1"+updn+"'",            "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (10, "'"+scratch+"/"+case+".vector"+updn+para+"'","'unknown'", "'unformatted'",9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (11, "'"+case+".energy"+updn+para+"'",        "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (18, "'"+case+".vsp"+updn+"'",                "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (19, "'"+case+".vns"+updn+"'",                "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (20, "'"+case+".struct"+"'",                  "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (21, "'"+case+".scf1"+updn+"'",               "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (55, "'"+case+".vec"+"'",                     "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (71, "'"+case+".nsh"+updn+"'",                "'unknown'", "'formatted'", 0)
    #print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (71, "'"+case+".nsh"+updn+"'",                "'unknown'", "'formatted'", 0)
    #print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (200,"'"+scratch+"/"+case+".storeHinv"+updn+"'",  "'replace'", "'unformatted'", 9000)
    writehs=False
    nmr=False
    if writehs:
        print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (12, "'"+scratch+"/"+case+".ham"+"'",         "'unknown'", "'unformatted'", 0)
    if nmr:
        print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (3,  "'"+case+".qvec"+"'",                "'unknown'", "'formatted'", 0)
    fdef.close()
    
def PrepareDefinitionFile_lapwso(widmf, case, para, scratch, cmplx, updn='',orb=False):
    fdef = open(widmf+updn+'.def', 'w')
    dnup=''
    if (updn=='up'): dnup='dn'
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 4, "'"+case+".in1"+cmplx+"'",                   "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 5, "'"+case+".inso'",                           "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 6, "'"+case+".outputso'",                       "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 8, "'"+case+".scfso'",                          "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 9, "'"+scratch+"/"+case+".vector"+dnup+para+"'","'unknown'", "'unformatted'",9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (10, "'"+scratch+"/"+case+".vectorup"+para+"'",   "'unknown'", "'unformatted'",9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (18, "'"+case+".vsp"+dnup+"'",                    "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (19, "'"+case+".vspup'",                          "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (20, "'"+case+".struct"+"'",                      "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (22, "'"+case+".vns"+dnup+"'",                    "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (23, "'"+case+".vnsup'",                          "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (41, "'"+scratch+"/"+case+".vectorsodn"+para+"'",   "'unknown'", "'unformatted'",9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (42, "'"+scratch+"/"+case+".vectorso"+updn+para+"'","'unknown'", "'unformatted'",9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (45, "'"+case+".normsodn'",                       "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (46, "'"+case+".normsoup'",                       "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (51, "'"+case+".energysodn"+para+"'",             "'unknown'", "'formatted'", 9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (52, "'"+case+".energyso"+updn+para+"'",          "'unknown'", "'formatted'", 9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (53, "'"+case+".energydum'",                      "'unknown'", "'formatted'", 9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (54, "'"+case+".energy"+dnup+para+"'",            "'unknown'", "'formatted'", 9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (55, "'"+case+".energyup"+para+"'",               "'unknown'", "'formatted'", 9000)
    if (orb):
        print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (11, "'"+case+".vorbdn'",                     "'unknown'", "'formatted'", 0)
        print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (12, "'"+case+".vorbup'",                     "'unknown'", "'formatted'", 0)
        print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (13, "'"+case+".vorbdu'",                     "'unknown'", "'formatted'", 0)


def PrepareInFile(idmf, mode, case, updn, nsymop=0, recomputeEF=None, mixEF=1.0, WL=None,Temperature=None, Q_ETOT=False):
    if Q_ETOT:
        Q_ETOT_ = '.T.'
    else:
        Q_ETOT_ = '.F.'
    fdin = open(case+'.indmf'+idmf, 'w')
    ## Looking for the chemical potential
    (EF, NOE) = FindChemicalPotential(case, updn)

    print >> fdin, '#-------------- mode and the current chemical potential --------------------'
    if recomputeEF is not None:
        print >> fdin, mode, EF, recomputeEF, mixEF, Q_ETOT_, Temperature, '        #  mode, EF, recomputeEF, mixing-EF, Temperature'
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



def runExternal(widmf, ROOT, MPI, nom, ntail, dmft_exe='dmft', m_ext='',performSplit=True):
    "Runs external executable"
    if performSplit:
        cmd = ROOT+'/ssplit.py -n '+str(nom)+' -t '+str(ntail)
        if len(m_ext)>0: cmd += ' -l '+m_ext
        print '..... running:', cmd
        subprocess.call(cmd, bufsize = 1,  shell=True)

    deffile = widmf+m_ext+'.def'
    
    if MPI: # In parallel mode needs to copy the executable to the current directory
        shutil.copy(ROOT+'/'+dmft_exe, os.getcwd() )
        exe = MPI + ' ./'+dmft_exe
    else:
        exe = dmfe.ROOT+'/'+dmft_exe

    cmd = exe + ' ' + deffile

    print '..... running:', cmd
    subprocess.call(cmd, bufsize = 1, shell=True)  # line-buffered

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
      --p option for paralllel execution
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("--so", dest="so", action='store_true', default=False, help="switches on spin-orbit coupling")
    parser.add_option("--up", dest="up", action='store_true', default=False, help="magnetic LDA calculation with vector-up first")
    parser.add_option("--dn", dest="dn", action='store_true', default=False, help="magnetic LDA calculation with vector-dn first")
    parser.add_option("-p", dest="para", action='store_true', default=False, help="turns on parallel mode")
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
    parser.add_option("-q", "--q_etot",  dest="Q_ETOT",  type="int", default=1, help="Compute total energy of free energy in dmft2 [0|1]")
    parser.add_option("-u", "--mode", dest="mode",  default='c',  help="Mode of calculation for dmft2. Can be 'c' for charge and 'm' for calculating EF only")


    # Next, parse the arguments
    (options, args) = parser.parse_args()


    print 'sys.argv=', sys.argv
    
    if len(args)!=1:
        print 'Need exactly one argument. One of the dmft0|dmft1|dmft2|mu|dmftp|dmftu'
        sys.exit(1)
    
    dmfe = utils.DmftEnvironment()  # DMFT paths
    w2k = utils.W2kEnvironment()    # W2k filenames and paths
    
    if not options.c and os.path.isfile(w2k.case+".in1c") and os.path.getsize(w2k.case+".in1c")>0 :
        print 'Found '+w2k.case+'.in1c file, hence assuming complex version needed. Switching -c switch!'
        options.c = True

    # processing arguments
    updn=''
    dnup='dn'
    so=''
    cmplx=''
    para=''
    if options.para:
        para = '_x'
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

    if args[0] == 'lapw1':
        exe = 'lapw1'+cmplx
        k_band = ''
        if (options._band): k_band = '_band'
        if para: print 'Parallel vector/energy files written.'
        PrepareDefinitionFile_lapw1('lapw1', w2k.case, updn, para, w2k.SCRATCH, k_band, cmplx)
        if not options.def_only:
            __updn = '' if options.m_ext else updn  # if m_ext adds up/dn then updn should not. Otherwise it should.... maybe you can find a better solution.
            runExternal('lapw1'+__updn, dmfe.ROOT, dmfe.MPI2, None, None, exe, options.m_ext, False)
        sys.exit(0)

    if args[0]=='lapwso':
        exe = 'lapwso'
        if para: print 'Parallel vector/energy files written.'
        PrepareDefinitionFile_lapwso('lapwso', w2k.case, para, w2k.SCRATCH, cmplx)
        if not options.def_only:
            runExternal(exe, dmfe.ROOT, dmfe.MPI2, None, None, exe, options.m_ext, False)
        sys.exit(0)
        
    # This is needed for dmft but not for lapw1
    if not options.so and os.path.isfile(w2k.case+".inso") and os.path.getsize(w2k.case+".inso")>0 :
        print 'Found '+w2k.case+'.vectorso file, hence assuming so-coupling exists. Switching -so switch!'
        options.so = True

    if options.so:
        so = 'so'
        sodum = 'dum'
        cmplx = 'c'
        
    # Processing 'case.indmfl' file
    inl = indmffile.Indmfl(w2k.case) # case.indmfl file
    inl.read()                       # case.indmfl read
    cixs = inl.siginds.keys()        # all columns for cix

    print 'case=', w2k.case
    #print 'c=', cmplx
    #print 'scratch=', w2k.SCRATCH
    print 'root=', dmfe.ROOT
    print 'cixs=', cixs
    #print 'def_only=', options.def_only
    
    if args[0] == 'dmft0':
        idmf = '0'
        mode = 'e'               # mode for computing eigenvalues
        
        PrepareDefinitionFile_dmft1(idmf, mode, w2k.case, cixs, updn, dnup, so, para, w2k.SCRATCH, cmplx)
        PrepareInFile(idmf, mode, w2k.case, updn, 1)

        if not options.def_only:
            #checkSigmas(cixs)
            runExternal('dmft'+idmf, dmfe.ROOT, dmfe.MPI2, options.nom, options.ntail, 'dmft')
            
    elif args[0] == 'dmft1':
        idmf = '1'
        mode = 'g'               # mode for computing eigenvalues
        #__updn = '' if options.m_ext else updn  # if m_ext adds up/dn then updn should not. Otherwise it should.... maybe you can find a better solution.
        #dn = options.m_ext+__updn
        PrepareDefinitionFile_dmft1(idmf+options.m_ext, mode, w2k.case, cixs, updn, dnup, so, para, w2k.SCRATCH, cmplx, '', options.m_ext)
        PrepareInFile(idmf+options.m_ext, mode, w2k.case, updn)
        
        if not options.def_only:
            #checkSigmas(cixs)
            runExternal('dmft'+idmf, dmfe.ROOT, dmfe.MPI2, options.nom, options.ntail, 'dmft', options.m_ext)
            
    elif args[0] == 'dmft2':
        idmf = '2'
        #mode = 'c'
        mode = options.mode
        Temperature = getTemperature()
        
        PrepareDefinitionFile_dmft2(idmf+options.m_ext, mode, w2k.case, cixs, updn, dnup, so, para, w2k.SCRATCH, cmplx, options.m_ext)
        PrepareInFile(idmf+options.m_ext, mode, w2k.case, updn, 0, options.recomputeEF, options.mixEF, options.WL,Temperature,options.Q_ETOT)

        if not options.def_only:
            #checkSigmas(cixs)
            runExternal('dmft'+idmf, dmfe.ROOT, dmfe.MPI2, options.nom, options.ntail, 'dmft2',options.m_ext)
        
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
            
        PrepareDefinitionFile_dmft1(idmf, mode, w2k.case, cixs, updn, dnup, so, para, w2k.SCRATCH, cmplx, k_band)
        PrepareInFile(idmf, mode, w2k.case, updn, nsymop)
        if not options.def_only:
            runExternal('dmft'+idmf, dmfe.ROOT, dmfe.MPI2, options.nom, options.ntail, 'dmft')
        
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

        PrepareInFile(idmf+options.m_ext, mode, w2k.case, updn,1)
        PrepareDefinitionFile_dmft1(idmf+options.m_ext, mode, w2k.case, cixs, updn, dnup, so, para, w2k.SCRATCH, cmplx, '_band', options.m_ext)

        nom = inl.om_npts 
        aom = inl.om_emin  
        bom = inl.om_emax  

        cmd = dmfe.ROOT+'/ssplit.py '
        if len(options.m_ext)>0: cmd += ' -l '+options.m_ext

        print '..... running:', cmd
        stdin, stdout, stderr = os.popen3( cmd )
        print stderr.read().strip(), stdout.read().strip()

        eq_om = linspace(inl.om_emin, inl.om_emax, inl.om_npts)  # Equidistant mesh
        #print 'eq_om=', eq_om
        
        # Interpolates sigmas on equidistant mesh!
        for c in cixs:
            fsig = open('sig.inp'+str(c)+options.m_ext, 'r')
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
                fs = interpolate.UnivariateSpline(om0, sig0[i],s=0)
                sig1.append( fs(eq_om) )
                #print 'sig1=', sig1
                #tck = interpolate.splrep(om0, sig0[i])
                #sig1.append( interpolate.splev(eq_om, tck) )
            sig1 = array(sig1).transpose()
            
            osig = open('sig.inp'+str(c)+options.m_ext+'_band', 'w')
            for i in range(len(eq_om)):
                print >> osig, '%15.8f  ' % eq_om[i], ('%15.8f '*len(sig1[i])) % tuple(sig1[i])
            osig.close()

        if not options.def_only:
            runExternal('dmft'+idmf, dmfe.ROOT, dmfe.MPI2, options.nom, options.ntail, 'dmft',options.m_ext,False)


        ####
        
    else:
        print 'not yet implemented!'

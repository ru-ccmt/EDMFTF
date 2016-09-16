#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import optparse
import os, shutil
import subprocess
from numpy import zeros, log
import indmffile
from utils import DmftEnvironment, W2kEnvironment
import sys



def print_prompt(text):
    print "----->", text


dmfenv = DmftEnvironment()
w2kenv = W2kEnvironment()

indmf = indmffile.Indmf(w2kenv.case)

# First, check if case.indmf exists
create_indmf = True
if indmf.file_exists():
    prompt = "File `%s` already exists.  Do you want to create new inputs? (y/n): " % indmf.filename()
    userin = raw_input(prompt).lower().strip()
    if userin.lower() == 'n':
        create_indmf = False

# create new case.indmf (if desired)
if create_indmf:
    indmf.user_input()
    indmf.write()
    fname = indmf.filename()
    fnametemp = fname + '_st'

    finished_editing = False
    while not finished_editing:
        shutil.copy(fname, fnametemp)  # copy to temp file for user to edit
        print_prompt("Edit %s file." % (fnametemp,))
        subprocess.call(w2kenv.EDITOR.split() + [fnametemp])

        # read user-edited file in just to check for any syntax errors
        indmf_st = indmffile.Indmf(w2kenv.case)
        indmf_st.extn = 'indmf_st'
        try:
            indmf_st.read()
            finished_editing = True
            shutil.copy(fnametemp, fname)  # move it back            
        except:
            print_prompt("Your edits have introduced syntax errors in %s." % (fnametemp,))
            prompt = "Do you want to edit the file again? (y/n): "
            userin = raw_input(prompt).lower().strip()
            if userin.lower() == 'n':
                finished_editing = True

# ask user if this is a SO run
so_run = False
prompt = "Is this a spin-orbit run? (y/n): "
userin = raw_input(prompt).lower().strip()
if userin.lower() == 'y':
    so_run = True

# run sigen to calculate siginds, cftrans and write output to case.indmfl file
sigen = os.path.join( dmfenv.ROOT, 'sigen.py' )
#sigen = './sigen.py'
args = ["--so"] if so_run else []
ireturn = subprocess.call([sigen] + args, shell=True,stdout=sys.stdout,stderr=sys.stdout)
if ireturn:
    print "Error executing sigen.py."
    exit(1)

# copy case.indmfl --> case.indmfl_st to allow user to make changes
findmfl = w2kenv.case + '.indmfl'
findmfltemp = findmfl + '_st'

finished_editing = False
while not finished_editing:
    shutil.copy(findmfl, findmfltemp)
    print_prompt("Edit %s file." % (findmfltemp,))
    subprocess.call(w2kenv.EDITOR.split() + [findmfltemp])

    # read user-edited file in just to check for any syntax errors
    inl = indmffile.Indmfl(w2kenv.case)
    inl.extn = 'indmfl_st'
    try:
        inl.read()
        finished_editing = True
        shutil.copy(findmfltemp, findmfl)  # move it back
    except:
        print_prompt("Your edits have introduced syntax errors in %s." % (findmfltemp,))
        prompt = "Do you want to edit the file again? (y/n): "
        userin = raw_input(prompt).lower().strip()
        if userin.lower() == 'n':
            finished_editing = True



# Generate case.indmfi file -----------------------------------

indmf.read()
inl = indmffile.Indmfl(w2kenv.case)
inl.read()

f = open(w2kenv.case + '.indmfi', 'w')
print >> f, len(indmf.ucps), "  # number of sigind blocks"
for iucp,icps in indmf.ucps.iteritems():
    icp = icps[0]  # just take sigind of first correlated problem  (TODO: modify appropriately for broken symmetry runs)
    sigind = inl.siginds[icp]
    print >> f, sigind.shape[0], "  # dimension of this sigind block"
    max_sigfig = 1 + int(log(max(sigind.flat))/log(10))
    format = '%' + str(max_sigfig) + 'd'
    for row in sigind:
        print >> f, ' '.join([format % elem for elem in row])
f.close()


# Generate blank sig.inp --------------------------------------

#prompt = "Generate blank sig.inp? (y/n): "
#userin = raw_input(prompt).lower().strip()
#if userin.lower() == 'y':
#
#    cmd = os.path.join( dmfenv.ROOT, 'szero.py' )
#    stdoutdata, stderrdata = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()
    

    #
    #
    ## use siginds to determine number of columns needed in sig.inp
    #inl = indmffile.Indmfl(w2kenv.case)
    #inl.read()
    #nsigs = max([max(sigind.flat) for sigind in inl.siginds.values()])  # assumes siginds start with 1
    #ncols = 2*nsigs   # real/imag for each sigma column
    #
    ## use gaumesh.py to generate standard mesh
    #gaumesh = os.path.join( dmfenv.ROOT, 'gaumesh.py' )
    ## gaumesh = './gaumesh.py'
    #args = [
    #    'x0=[0,0]',
    #    'dx0=[0.25,0.01]',
    #    'fwhm=[30, 3]',
    #    'xmin=-20',
    #    'xmax=20'
    #    ]
    #stdoutdata, stderrdata = subprocess.Popen([gaumesh] + args, stdout=subprocess.PIPE).communicate()
    #omegas = [x for x in stdoutdata.split('\n') if x.strip() and not x.strip().startswith('#')]
    #
    ## if sig.inp already exists, make backup
    #if os.path.isfile('sig.inp'):
    #    print "File `sig.inp` already exists, backing up to 'sig.inp.bak'"
    #    shutil.move('sig.inp', 'sig.inp.bak')
    #
    ## now generate zeros and write to file
    #f = open('sig.inp', 'w')
    #nsig_zeros = ', '.join(['0.0']*nsigs)       # first generate s_oo and Edc
    #print >> f, '# s_oo = [' + nsig_zeros + ']'
    #print >> f, '# Edc  = [' + nsig_zeros + ']'
    #ncol_zeros = ' '.join(['0.0']*ncols)       # now write out main file
    #for omega in omegas:
    #    print >> f, omega, ncol_zeros
    #f.close()

#    print 'Blank self-energy written to `sig.inp`.'

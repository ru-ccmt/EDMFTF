#!/usr/bin/env python
import sys
import re, os, glob, shutil
from os.path import getsize

if len(sys.argv)<2:
    print """The following arguments are needed:
             > dmft_copy.py dir_name -[l|d|a]

             dir_name        -- the name of the source directory with a wien2k/dmft run
             shitch can be d,l, a :  l -- the files for LDA execution will be copied
                                     d -- the files for DMFT will be compied
                                     a -- files for both LDA and DMFT will be copied
    """      
    sys.exit(0)

cpdr = sys.argv[1]

if not os.path.exists(cpdr):
    print 'The directory '+cprd+' does not exists!'
    sys.exit(1)

if len(sys.argv)==2:
    opt = '-a'
else:
    opt = sys.argv[2]
    if opt not in ['-l','-d','-a']:
        print 'Switch can be either  -l,  -d or  -a'

#if len(sys.argv)<=3:
    
strcase = glob.glob(cpdr+'/*.struct')[0]
if not strcase:
    print 'Can not find structure file!'
    sys.exit(1)

(dirName, fileName) = os.path.split(strcase); (fileBaseName, fileExtension)=os.path.splitext(fileName)
case = fileBaseName

#else:
#    case = sys.argv[3]
#    if not os.path.exists(cpdr+'/'+case+'.struct'):
#        print 'Can not find structure file '+cpdr+'/'+case+'.struct'
#        sys.exit(1)
        

print 'case=', case

w0 = [case+'.struct',case+'.in0',case+'.clmsum', case+'.inm']
w1 = [case+'.in1', case+'.in1c', case+'.klist', case+'.inso', case+'.in2c']
w2 = [case+'.in2', case+'.kgen']
wc = [case+'.inc', case+'.scf2']

if opt in ['-l', '-a']:
    for f in w0+w1+w2+wc:
        fle = glob.glob(cpdr+'/'+f)
        if fle:
            print 'Copying ... ', fle[0]
            shutil.copy2(fle[0], '.')
            
d0 = ['params.dat', 'EF.dat', 'sig.inp', case+'.indmfl', case+'.indmfldn', case+'.indmfi', 'Sigma.000', 'projectorw.dat', 'Eorb_imp.dat'] # case+'.vsp', case+'.energy']
dn = ['sfx.*']
if opt in ['-d', '-a']:
    for f in d0:
        fle = glob.glob(cpdr+'/'+f)
        if fle:
            print 'Copying ... ', fle[0]
            shutil.copy2(fle[0], '.')

    for f in dn:
        fle = glob.glob(cpdr+'/'+f)
        for fl in fle:
            print 'Copying ... ', fl
            shutil.copy2(fl, '.')

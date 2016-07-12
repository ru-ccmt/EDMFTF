#!/usr/bin/env python

import re,glob,os,sys,copy,shutil

fkeep0=['params.dat', 'info.iterate', 'wnohup.dat','dmft_info.out', 'ksum_info.out','sig.inp', 'dclean.py', 'optics.dat']
fkeep=['sig.inp.*', '*.tgz']
ckeep=['.dlt?', '.gc?', '.cdos', '.outputdmf?.0', '.indmf', '.indmf[0-2]','.struct']

fremove0=['dmft']
fremove= ['dmft?.error.*']
cremove=['.outputdmf?.[1-9]*', '.vector', '.vectordn', '.vectorup', '.vector.gz', '.vectordn.gz', '.vectorup.gz']

impfiles = ['ctqmc', 'Probability\.dat\..*', 'status\..*', 'check.dat.*', 'g_hb\d\..*', 'new.cix', 'nohup_imp.out.*', 'nohup.out', 'Delta.tau\..*', 'Sigma\..*', 'oca_log\..*', 'Ac\.inp\..*']


if len(sys.argv)<2:
    print 'Give directory name'
    sys.exit(1)
else:
    drn = sys.argv[1]
    print 'Directory is', drn


case = os.path.splitext(os.path.split(glob.glob(drn+"/*.struct")[0])[1])[0]
print 'case=', case

#case = os.path.splitext(glob.glob(drn+"/*.struct")[0])[0] # finds file xxx.struct and sets case to xxx
#print 'case=', case


# First, find all files on thid directory tree
allf=[]
imps=[]
impzip={}
for root, dirs, files in os.walk(drn):
    for name in files:
        if root==drn:         # files on current directory
            allf.append(name)
        elif re.match(drn+'/?imp.', root) is not None:
            #print root, name
            for imp in impfiles:
                if re.match(imp, name):  # should be zipped
                    
                    if impzip.has_key(root):
                        impzip[root].append(name)
                    else:
                        impzip[root]=[name]
                        
                    #print 'matches!'
            
    for name in dirs:
        pass
        #print 'directory=', os.path.join(root, name)

#print 'impzip=', impzip
for root in impzip.keys():
    cmd = 'cd '+root+'; tar czvf ofiles.tgz '+(' '.join(impzip[root]))
    stdin, stdout, stderr = os.popen3(cmd)
    print stdout.read(), stderr.read()
    cmd = 'cd '+root+'; rm '+(' '.join(impzip[root]))
    stdin, stdout, stderr = os.popen3(cmd)
    print stdout.read(), stderr.read()
    
# Files to keep
keep = []
for fk in fkeep+fkeep0:
    gg = glob.glob(drn+'/'+fk)
    print 'fk=', drn+'/'+fk, 'gg=', gg
    keep.extend( gg )

for fk in ckeep:
    gg = glob.glob(drn+'/'+case+fk)
    print 'fk=', drn+'/'+case+fk, 'gg=', gg
    keep.extend( gg )
    
# Files to remove
remove = [] #copy.deepcopy( fremove0 )
for dk in fremove+fremove0:
    gg = glob.glob(drn+'/'+dk)
    print 'dk=', drn+'/'+dk, 'gg=', gg
    remove.extend( gg )
    
for dk in cremove:
    gg = glob.glob(drn+'/'+case+dk)
    print 'dk=', drn+'/'+case+dk, 'gg=', gg
    remove.extend( gg )


print 'keep=', keep


for k in keep+remove:
    name = os.path.split(k)[1]
    if name in allf:
        allf.remove(name)


print 'Zipping ', len(allf), 'files'
if len(allf)>0:
   root=drn
   cmd = 'cd '+root+'; tar czvf ofiles.tgz '+(' '.join(allf))
   stdin, stdout, stderr = os.popen3(cmd)
   print stdout.read(), stderr.read()

print 'Removing ', len(remove), 'files'
if len(remove)>0:
   cmd = 'rm '+(' '.join(remove))
   stdin, stdout, stderr = os.popen3(cmd)
   print stdout.read(), stderr.read()

print 'Cleaning ', len(allf), 'files'
if len(allf)>0:
   cmd = 'cd '+drn+'; rm '+(' '.join(allf))
   stdin, stdout, stderr = os.popen3(cmd)
   print stdout.read(), stderr.read()

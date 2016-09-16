#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule

from scipy import *
from scipy import linalg
import utils
import struct1
import indmffile
import latgen2
from MagneticMoment import *
from Wigner import *
import transformations as tr
import sys,optparse

def toUnitCell(r):
    tol=1e-6
    for i in range(3):
        if r[i]<-tol: r[i]+=1.
        if r[i]>1-tol: r[i]-=1.
    return r

def GetRotij(atmu,st, fil=sys.stdout):
    iat,imu = atmu
    at_first = st.pos[iat][0]
    at_curr  = st.pos[iat][imu]
    lattic = st.lattice
    
    #print 'lattice type=', st.lattice
    print >> fil, 'first atom position=', at_first, 'current atom position', at_curr
    rotij = zeros((3,3))
    Found=False
    for i in range(st.Ng):
        Det=linalg.det(st.iz[i])
        x = toUnitCell(dot(st.iz[i].T,at_first) + st.tau[i])
        x1 = toUnitCell(abs(x-at_curr))
        #print "Operation %-2d gives vector " % (i+1), '%4.2f'*3 % tuple(x), 'and difference', ' %4.2f'*3 % tuple(x1)
        if sum(abs(x1))<1e-4:
            rotij = st.iz[i]
            Found = True
            break
        if lattic[0] == 'B':
            x2 = toUnitCell(x1 + array([0.5,0.5,0.5]))
            #print 'x2=', ' %4.2f'*3 % tuple(x2)
            if sum(abs(x2))<1e-4:
                rotij = st.iz[i]
                Found = True
                break
        if lattic[0] == 'F' or lattic[0:3] == 'CXY':
            x2 = toUnitCell(x1 + array([0.5,0.5,0.0]))
            #print 'x2=', ' %4.2f'*3 % tuple(x2)
            if sum(abs(x2))<1e-4:
                rotij = st.iz[i]
                Found = True
                break
        if lattic[0] == 'F' or lattic[0:3] == 'CXZ':
            x2 = toUnitCell(x1 + array([0.5,0.0,0.5]))
            #print 'x2=', ' %4.2f'*3 % tuple(x2)
            if sum(abs(x2))<1e-4:
                rotij = st.iz[i]
                Found = True
                break
        if lattic[0] == 'F' or lattic[0:3] == 'CYZ':
            x2 = toUnitCell(x1 + array([0.0,0.5,0.5]))
            #print 'x2=', ' %4.2f'*3 % tuple(x2)
            if sum(abs(x2))<1e-4:
                rotij = st.iz[i]
                Found = True
                break
    if Found:
        print >> fil, 'Found Operation %-2d transforms first atom to atom %2d:%-2d  det=%4.1f' % (i+1,iat,imu,Det)
    else:
        print >> fil, 'ERROR: Can not find group operation which would give the atom in the unit cell'
        print 'ERROR: Can not find group operation which would give the atom in the unit cell'
        exit(1)
    return rotij

def Angles(n):
  nm=sqrt(n[0]**2+n[1]**2+n[2]**2)
  theta = arccos(n[2]/nm)
  xx = sqrt(n[0]**2+n[1]**2)
  if xx<1e-5:
    phi = 0.0
  else:
    phi = arccos(n[0]/xx)
    if abs(n[1])>1e-5:
        phi = phi * n[1]/abs(n[1])
  return (theta,phi)

def MatrixToMix_2x2(n, mc_x,mc_y,mc_z):
    # first compute 3x3 matrix to transform to a basis with vector n pointing in new z
    (theta,phi) = Angles(n)
    cf,sf,ct,st = cos(phi),sin(phi),cos(theta),sin(theta)
    # rot transforms to the basis in which n points in new z-direction
    rot=[[ct*cf,  -sf,  st*cf],
        [ct*sf,   cf,  st*sf],
        [  -st,    0,  ct   ]]
    Ri = transpose(rot)
    # Magnetization in the rotated basis where n is pointing in new z
    mr_x = Ri[0,0]*mc_x + Ri[0,1]*mc_y + Ri[0,2]*mc_z
    mr_y = Ri[1,0]*mc_x + Ri[1,1]*mc_y + Ri[1,2]*mc_z
    mr_z = Ri[2,0]*mc_x + Ri[2,1]*mc_y + Ri[2,2]*mc_z
    #print 'Mx='
    #fprint(mr_x)
    #print 'My='
    #fprint(mr_y)
    #print 'Mz='
    #fprint(mr_z)
    # Now we solve a system of non-linear equations to find 2x2 Unitary matrix of the form
    # Ut = [[cos(thta),sin(thta)*exp(phii*1j)],[-sin(thta)*exp(-phii*1j),cos(thta)]]
    # which makes mr_x[0,0]=0 and mr_y[0,0]=0
    # We are looking for two angles thta and phii, and we have two conditions
    a11 = mr_x[0,0].real
    a12 = mr_x[0,1]
    b11 = mr_y[0,0].real
    b12 = mr_y[0,1]
    #print 'a11,b11,a12,b12=', a11, b11, a12, b12
    phii = arctan2((a11*b12.real - b11*a12.real), (b11*a12.imag - a11*b12.imag))
    thta = arctan2( a11, (cos(phii)*a12.real + sin(phii)*a12.imag) )/2.
    Ut = matrix([[cos(thta),sin(thta)*exp(phii*1j)],[-sin(thta)*exp(-phii*1j),cos(thta)]])
    return Ut

def PointMomentInDirection(direction, M, cfx, orb_ud):
    """ Takes time-reversal pairs of vectors in cfx and creates linear combination of such two vectors 
    (members of the timer reversal pairs)
    such that the resulting magnetic moment points in direction specified by direction vector.
    """
    T0 = matrix(conj(cfx))
    M_x = T0 * M.Mx * T0.H
    M_y = T0 * M.My * T0.H
    M_z = T0 * M.Mz * T0.H
    
    CF = copy(cfx)  # origonal transformation
    for l,m in orb_ud:
        m_x = array([[M_x[l,l],M_x[l,m]],[M_x[m,l],M_x[m,m]]])
        m_y = array([[M_y[l,l],M_y[l,m]],[M_y[m,l],M_y[m,m]]])
        m_z = array([[M_z[l,l],M_z[l,m]],[M_z[m,l],M_z[m,m]]])
        c = MatrixToMix_2x2(direction, m_x,m_y,m_z)
        ct = transpose(c)
        CF[l] = ct[0,0]*cfx[l,:]+ct[0,1]*cfx[m,:]
        CF[m] = ct[1,0]*cfx[l,:]+ct[1,1]*cfx[m,:]
    return matrix(CF)

def Direction(x,y,z):
    n = array([x,-y,z])
    nz = n/sqrt(sum(n*n))
    if abs(nz[2]-1.)<1e-4: # pointing in z -- no rotation
        return array([[1,0],[0,1]])
    ez = array([0,0,1.])
    nx = ez - dot(ez,nz)*nz
    nx *= 1./sqrt(sum(nx**2))
    ny = cross(nz,nx)
    R = array([nx,ny,nz])
    (al,be,ga) = tr.euler_from_matrix(R, 'rzyz')
    Us = WignerMatrix(0.0,be,ga, 0.5)
    return Us

def Direction_old(x,y,z):
    nm = sqrt(z**2+y**2+x**2)
    sz = abs(z)/nm
    sx = x/nm
    sy = y/nm
    alp = sqrt(0.5*(1+sz))
    bet = sqrt(0.5*(1-sz))*exp(1j*arctan2(sx,sy))
    c=1.
    if (sz!=sx): c = sz/(sz-sx)
    res = [[alp,bet],[alp+c*(bet-alp),-bet+c*(alp+bet)]]
    if  (x<0 and y<0) or z<0 :
        return array([res[1],res[0]])
    return array(res)

def Direction2(z1,z2,z3):
   def TZ(newz_):
       "Given new z-axis, gives orthogonal 3x3 matrix of coordinate system"
       newz = array(newz_)/sqrt(dot(newz_,newz_))
       rn=newz.tolist().index(min(newz)) # the smallest component
       newx = zeros(3); newx[rn]=1.
       newx -= dot(newx,newz)*newz # orthogonalize to z
       newx = newx/sqrt(dot(newx,newx)) # normalize
       newy=cross(newz,newx)
       newy=newy/sqrt(dot(newy,newy))
       #print 'TZ=', array([newx,newy,newz])
       return array([newx,newy,newz])

   return SpinRotation(TZ([z1,z2,z3]))

def Mix2Eigenvectorsn(TT,c,orb_ud):
    """
    T[0,:] = c[0,0]*TT[0,:]+c[0,1]*TT[1,:]
    T[1,:] = c[1,0]*TT[0,:]+c[1,1]*TT[1,:]
    T[2,:] = c[0,0]*TT[2,:]+c[0,1]*TT[3,:]
    T[3,:] = c[1,0]*TT[2,:]+c[1,1]*TT[3,:]
    T[4,:] = c[0,0]*TT[4,:]+c[0,1]*TT[5,:]
    T[5,:] = c[1,0]*TT[4,:]+c[1,1]*TT[5,:]
    T[6,:] = c[0,0]*TT[6,:]+c[0,1]*TT[7,:]
    T[7,:] = c[1,0]*TT[6,:]+c[1,1]*TT[7,:]
    T[8,:] = c[0,0]*TT[8,:]+c[0,1]*TT[9,:]
    T[9,:] = c[1,0]*TT[8,:]+c[1,1]*TT[9,:]
    """
    T = zeros(shape(TT),dtype=complex)
    for n in range(len(TT)/2):
        for i in range(2):
            #T[2*n+i,:] = c[i,0]*TT[2*n+0,:]+c[i,1]*TT[2*n+1,:]
            T[orb_ud[n][i],:] = c[i,0]*TT[orb_ud[n][0],:]+c[i,1]*TT[orb_ud[n][1],:]
            #print n,i, '1: ', 2*n+i,orb_ud[n][i], '2:', 2*n, orb_ud[n][0], '3:', 2*n+1,orb_ud[n][1]
    T = matrix(T)
    return T


def Find_orb_ud(Tu):
    """Finds time reversal pairs of orbitals
    """
    def  AreVectorsTimeReversal(Tl,Tm):
         N = len(Tl)
         im = index_max(abs(Tl[:]))
         phase = (-1)**(im)*conjugate(Tm[N-1-im])/Tl[im]
         if abs(abs(phase)-1)>0.5:
             return 100.
         phase = phase/abs(phase)
         dsum=0.0
         for i in range(N):
             c1, c2 = Tl[i]*phase, (-1)**(i)*conjugate(Tm[N-1-i])
             dsum += abs(c1-c2)
             #print ('%3d '+('%10.6f '*6)) % (i, (c1-c2).real, (c1-c2).imag, c1.real, c1.imag, c2.real, c2.imag)
         return dsum
    def index_max(values):
        return max(xrange(len(values)),key=values.__getitem__)
    def index_min(values):
        return min(xrange(len(values)),key=values.__getitem__)
    
    Tu=array(Tu)
    N = len(Tu)
    orb_ud=[]
    still_left=range(N)
    while(still_left):
        l = still_left[0]
        cm=ones(N)*100
        for m in still_left[1:]:
            cm[m]=AreVectorsTimeReversal(array(Tu)[l],array(Tu)[m])
        partner=index_min(cm)
        orb_ud.append([l,partner])
        #print l,partner, still_left
        still_left.remove(l)
        still_left.remove(partner)
    return orb_ud

def fprint(Us, fil=sys.stdout):
    Qcomplex =  type(Us[0,0])==complex or type(Us[0,0])==complex128
    for i in range(shape(Us)[0]):
        for j in range(shape(Us)[1]):
            if Qcomplex:
                print >> fil, "%11.8f %11.8f  " % (real(Us[i,j]), imag(Us[i,j])),
            else:
                print >> fil, "%11.8f  " % Us[i,j],
        print  >> fil



if __name__ == '__main__':


    s2 = sqrt(2.)
    cy=array([[ (1-1j)/2.,1/sqrt(2.)], [-(1-1j)/2.,1/sqrt(2.)]])
    cx=array([[ (1+1j)/2.,1/sqrt(2.)], [-(1+1j)/2.,1/sqrt(2.)]])
    cz=array([[1,0],[0,1]])
    mcy=array([[-(1-1j)/2.,1/sqrt(2.)], [(1-1j)/2.,1/sqrt(2.)]])
    #cc=array([[1./sqrt(2.),1/sqrt(2.)],[1./sqrt(2.),-1./sqrt(2.)]])  #  local-x
    
    filename = 'moments.dat'

    usage = """usage: %prog [ options ] mode
    Given DMFT input files, returns the matrix of the magnetic moment matrix elements in global or local (attached to each atom) coordinate system.
    Can also be used to point magnetic moment in arbitrary direction.
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("-a", dest="all", action='store_true', default=False, help="prints all moments, global, local, spin, orbital")
    #parser.add_option("-p",  dest="prnt",  action='store_true', default=False, help="prints to file magnetic moment matrix elements")
    # Next, parse the arguments
    (options, args) = parser.parse_args()


    w2k = utils.W2kEnvironment()

    inl = indmffile.Indmfl(w2k.case)
    inl.read()

    st = struct1.Struct(w2k.case)

    BR1,BR2,vol,ortho = latgen2.latgen([st.alpha,st.beta,st.gamma], st.a, st.b, st.c, st.lattice)
    BR1 = matrix(BR1)
    #print 'BR1=', BR1
    
    M = Momentd()
    
    atmu=[]
    for iat in range(len(st.pos)):
        for imu in range(len(st.pos[iat])):
            atmu.append( (iat,imu) )
            #print iat, imu, st.pos[iat][imu]

    print '*** You have the following atoms treated by DMFT: ***'
    for ii,icp in enumerate(inl.cps.keys()):
        iat,l,qsplit = inl.cps[icp][0]
        iat,imu = atmu[iat-1]
        position  = st.pos[iat][imu]
        print st.aname[iat], "[%3d]" % ii, "%7.4f "*3 % tuple(position)
    Ans1 = raw_input('Do you want to point magnetic moments in a new specified direction? [y/n]:')
    if Ans1 in ['y','Y']:
        while True:
            direction=[]
            print '*** Please give your answer in python notaions, such as the example: [1,1,1]'
            print '    Valid answer is also None, which does not change the current direction'
            for ii,icp in enumerate(inl.cps.keys()):
                iat,l,qsplit = inl.cps[icp][0]
                iat,imu = atmu[iat-1]
                position  = st.pos[iat][imu]
                statom = st.aname[iat] + ("[%3d]" % ii) + ("%7.4f "*3 % tuple(position))
                Ans2 = raw_input('Direction for atom '+statom+'?  ')
                direction.append( eval(Ans2) )
            sdirection = str(direction)
            Ans3 = raw_input('Your input was '+sdirection+'. Is that OK [y/n]?  ')
            if Ans3 in ['y','Y']: break
    else:
        direction = [None]*len(inl.cps.keys())
    
    Ans1 = raw_input('Do you want short or long output [s/l]?  ')
    if Ans1 == 'l':
        options.all = True
    print 'Resuts written to output file '+filename+' .'
    fo = open(filename, 'w')
    
    for ii,icp in enumerate(inl.cps.keys()):
        iat,l,qsplit = inl.cps[icp][0]

        locrot,crotloc,shft = inl.atoms[iat][:3]
        #locrot,crotloc,shft = inl.atoms[icp]
        if locrot==0 or crotloc==[]: 
            crotloc=identity(3)
        
        print >> fo, '**** Correlated atom', iat, ' *****'
        #print 'crotloc=', crotloc
        
        Sigind = inl.siginds[icp]
        cfx = inl.cftrans[icp]
        rotij = GetRotij(atmu[iat-1], st, fo)
        

        orb_ud = Find_orb_ud(cfx)
        print >> fo, 'time reversal pairs of states orb_ud=', orb_ud
        print >> fo, 'crotloc=', crotloc

        rotij_ = BR1*rotij*BR1.I
        global2local = matrix(crotloc) * rotij_

        if (direction[ii] is not None):
            CF = PointMomentInDirection(direction[ii], M, cfx, orb_ud)
        else:
            CF = matrix(cfx)
        
        T2 = conj(CF)
        
        Mx_local = T2 * M.Mx * T2.H
        My_local = T2 * M.My * T2.H
        Mz_local = T2 * M.Mz * T2.H
        
        
        Ds = OrbitSpinRotation(global2local,2)
        Dt = T2*Ds
        
        Mx_global = Dt * M.Mx * Dt.H
        My_global = Dt * M.My * Dt.H
        Mz_global = Dt * M.Mz * Dt.H


        Lx_local = T2 * M.Lx *  T2.H
        Ly_local = T2 * M.Ly *  T2.H
        Lz_local = T2 * M.Lz *  T2.H
        Sx_local = T2 * M.Sgx * T2.H
        Sy_local = T2 * M.Sgy * T2.H
        Sz_local = T2 * M.Sgz * T2.H

        Lx_global = Dt * M.Lx * Dt.H
        Ly_global = Dt * M.Ly * Dt.H
        Lz_global = Dt * M.Lz * Dt.H
        Sx_global = Dt * M.Sgx * Dt.H
        Sy_global = Dt * M.Sgy * Dt.H
        Sz_global = Dt * M.Sgz * Dt.H
        
        print >> fo, 'rotij='
        fprint(rotij, fo)
        print >> fo, 'rotij_='
        fprint(rotij_, fo)
        print >> fo, 'global2local='
        fprint(global2local, fo)
        print >> fo, 'Wigner rotation='
        fprint( Ds , fo)
        print >> fo, 'New CF Transformation='
        fprint(CF, fo)
        
        #print 'cfx transformation matrix:'
        #fprint(CF)
        
        print >> fo, 'Mx_global='
        fprint(Mx_global, fo)
        print >> fo, 'My_global='
        fprint(My_global, fo)
        print >> fo, 'Mz_global='
        fprint(Mz_global, fo)

        #if options.prnt:
        #    fo = open('Moment.x.'+str(icp),'w')
        #    fprint(Mx_global,fo)
        #    fo.close()
        #    fo = open('Moment.y.'+str(icp),'w')
        #    fprint(My_global,fo)
        #    fo.close()
        #    fo = open('Moment.z.'+str(icp),'w')
        #    fprint(Mz_global,fo)
        #    fo.close()
            

        if options.all:
            print >> fo, 'Lx_global='
            fprint(Lx_global, fo)
            print >> fo, 'Ly_global='
            fprint(Ly_global, fo)
            print >> fo, 'Lz_global='
            fprint(Lz_global, fo)
            
            print >> fo, '2Sx_global='
            fprint(Sx_global, fo)
            print >> fo, '2Sy_global='
            fprint(Sy_global, fo)
            print >> fo, '2Sz_global='
            fprint(Sz_global, fo)
            
            print >> fo, 'Mx_local='
            fprint(Mx_local, fo)
            print >> fo, 'My_local='
            fprint(My_local, fo)
            print >> fo, 'Mz_local='
            fprint(Mz_local, fo)
            print >> fo

            print >> fo, 'Lx_local='
            fprint(Lx_local, fo)
            print >> fo, 'Ly_local='
            fprint(Ly_local, fo)
            print >> fo, 'Lz_local='
            fprint(Lz_local, fo)
            
            print >> fo, '2Sx_local='
            fprint(Sx_local, fo)
            print >> fo, '2Sy_local='
            fprint(Sy_local, fo)
            print >> fo, '2Sz_local='
            fprint(Sz_local, fo)

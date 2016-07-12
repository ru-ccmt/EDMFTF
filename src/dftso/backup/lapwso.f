PROGRAM lapwso
  USE abcd
  USE loabc; USE loabcr; USE lolog; USE rlolog; USE rpars; USE orb
  USE struct
  USE rotmat
  USE vns
  USE ams
  USE hmsout
  USE couples
  IMPLICIT REAL*8(A-H,O-Z)
  
  COMMON/GENER/  BR1(3,3),BR2(3,3)

  character*80    deffn,errfn,fname
  CHARACTER*11    STATUS,FORM
  character*10    bname
  logical         fl,flcom

  !n matrix elements <L,ML,MS|l.s|L,ML',MS'> calculated in COUPLE
  !n MS=1 spin down, MS=2 spin up

  dimension       nv(2),ne(2),emm(3)
  complex*16,ALLOCATABLE :: meigve(:,:,:)
  REAL*8,ALLOCATABLE     :: vru(:,:,:),e(:,:,:)
  REAL*8,ALLOCATABLE     :: p(:,:,:),dp(:,:,:),pe(:,:,:),dpe(:,:,:)
  real*8,ALLOCATABLE     :: ri_mat(:,:,:,:,:,:),ee(:,:)
  real*8,ALLOCATABLE     :: ri_orb(:,:,:,:,:,:)
  INTEGER,ALLOCATABLE    :: kv(:,:,:)
  dimension       ss(3),cp(5)
  flcom=.false.
  call gtfnam(deffn,errfn)
  !     call errflg(errfn,'Error in SPINORB')
  OPEN(1,FILE=DEFFN,STATUS='OLD',ERR=8000)
8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
  !
  !        check for in1c file for complex case
  !
  DO I = LEN(fname),1,-1
     IF (fname(I:I) .EQ. '.') THEN
        IF (fname(i+1:i+4) .eq. 'in1c') flcom=.true.
        EXIT
     ENDIF
  ENDDO
     !
11 OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=8002)
  GOTO 8003
8000 WRITE(*,*) ' ERROR IN OPENING LAPWSO.DEF !!!!'
  STOP 'LAPWSO.DEF'
8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
  WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, '  FORM:',FORM
  STOP 'OPEN FAILED'
8001 CONTINUE
  CLOSE (1)
     
  CALL init_struct
  CALL init_ams
  jspin=1
  DO I=1,NAT 
     READ(54,'(f9.5)') EMIST
     READ(54,'(f9.5)') EMIST
  ENDDO
  DO I=1,NAT                                                  
     READ(55,'(f9.5)',IOSTAT=ios) EMIST
     READ(55,'(f9.5)',IOSTAT=ios) EMIST
  ENDDO
  if (ios.eq.0) jspin=2
  DO
     READ(54,'(3e19.12,a10,2i6)',IOSTAT=ios) s,t,z,bname,n,nen
     IF (ios /= 0) EXIT
     nmat=MAX(n,nmat)
     nume=MAX(nen,nume)
     DO ii=1,nen
        READ(54,*)
     ENDDO
     IF(jspin.EQ.2) THEN
        READ(55,'(3e19.12,a10,2i6)',IOSTAT=ios) S,T,Z,bNAME,N,NEn
        nmat=MAX(n,nmat)
        nume=MAX(nen,nume)
        DO ii=1,nen
           READ(55,*) NUM,E1
        ENDDO
     ENDIF
  ENDDO
  CALL init_rotmat(ndif)
  CALL init_loabc(lomax,nloat,nato)
  CALL init_loabcr(lomax,nloat,nato)
  CALL init_lolog(nato,lomax,lmx)
  CALL init_rlolog(nato,lomax)
  CALL init_rpars(nato,lomax)
  CALL init_orb(nato,ndif)
  CALL init_vns(nato)
  ALLOCATE(vru(nrad,nato,2),e(lmx,nato,2))
  read(18,2030,end=3031) 
  read(19,2030,end=3031) 
  jspin=2
  goto 3032
3031 jspin=1
3032 continue
  
2030 format(///)
  write(6,*)'*********** spin-orbit calculation **************'
  if(jspin.eq.1)write(6,*)' non-spin polarized case'
  if(jspin.eq.2)write(6,*)' spin polarized case'
  rewind(18)
  rewind(19)
  
  do ispin=1,2
     do jatom=1,nato
        do l=0,lomax
           elor(l,jatom,ispin)=1.e+6
        enddo
     enddo
  enddo
  
  call cputim(dt0)
  call init(deffn,errfn,e,vru,emm,fl,jspin,kpot,ipr,theta,phi,irlotot)
  ALLOCATE(p(labc+1,nato,2),dp(labc+1,nato,2),pe(labc+1,nato,2),dpe(labc+1,nato,2))
  ALLOCATE(ri_mat(4,4,0:labc,nato,2,2))
  ALLOCATE(ri_orb(4,4,0:labc,nato,2,2))
  write(6,*)nmat,nnrlo
  nmat=nmat+nnrlo
  nume=nume+nnrlo
  nume2=  2*nume
  num2= nume2*(nume2+1)/2+(nume2/hblock+1)*(hblock*(hblock-1)/2)   
  
  CALL init_abcd(labc2,ndif,nume)
  CALL init_hmsout(nmat,nume2)
  ALLOCATE(meigve(nmat,nume,2),kv(3,nmat,2),ee(nume,2))
  call garadme(e,vru,p,dp,pe,dpe,ri_mat,jspin,kpot,ipr)
  
  if(iorbpot.ne.0)then
     call garadorb(e,vru,p,dp,pe,dpe,ri_orb,jspin,kpot,ipr)
  endif
  call gaunt2
  call cputim(dtime1)
  ktap=40
  do isi=1,jspin
     ktap=ktap+1
     do ity=1,nat
        do l=1,lmx
           if(.not.lapw(l-1,ity)) e(l,ity,isi)=e(l,ity,isi)+200.d0
        enddo
        write(ktap)(e(l,ity,isi),l=1,lmx)
        write(ktap)((elo(l,nn,ity,isi),l=0,lomax),nn=1,nloat)
        write(ktap+10,'(100(f9.5))')(e(l,ity,isi),l=1,lmx)
        write(ktap+10,'(100(f9.5))')((elo(l,nn,ity,isi),l=0,lomax),nn=1,nloat)   
        if (jspin.eq.1) then
           write(ktap+1)(e(l,ity,isi),l=1,lmx)
           write(ktap+1)((elo(l,nn,ity,isi),l=0,lomax),nn=1,nloat)
           write(ktap+11,'(100(f9.5))')(e(l,ity,isi),l=1,lmx)
           write(ktap+11,'(100(f9.5))')((elo(l,nn,ity,isi),l=0,lomax),nn=1,nloat)
        endif
        do l=1,lmx
           if(.not.lapw(l-1,ity)) e(l,ity,isi)=e(l,ity,isi)-200.d0 
        enddo
     end do
  end do
  
     
  DO kkk=1,1000000
     do isi=1,jspin
        call cputim(dtime0)
        call kptin(flcom,meigve,ss,ne(isi),nv(isi),ee,kv,bname,weight,kkk,isi)  ! meigve contains eigenvector.
        call cputim(dtime1)
        cp(1)=cp(1)+dtime1-dtime0
        if (kkk.lt.0) EXIT
        call cputim(dtime0)
        call abclm(meigve,ss,ne(isi),nv(isi),kv,p,dp,pe,dpe,isi)
        call cputim(dtime1)
        cp(2)=cp(2)+dtime1-dtime0
     enddo
     if(jspin.eq.1)then
        nv(2)=nv(1)
        ne(2)=ne(1)
     endif
     
     nban2=ne(1)+ne(2)+2*nnrlo
     if(.not.fl) then
        nv(1)=0
        nv(2)=0
     endif
     call cputim(dtime0)	
     call hmsec(fl,emm,ne,nv,ee,meigve,ri_mat,ri_orb,kkk,jspin,ipr)
     call cputim(dtime1)
     cp(3)=cp(3)+dtime1-dtime0
     call cputim(dtime0)
     call kptout(ss,bname,weight,kkk,kv,jspin,nv,ne)
     call cputim(dtime1)
     cp(4)=cp(4)+dtime1-dtime0
  ENDDO
  
  call cputim(dt4)
  write(6,*) 'TOTAL NUMBER OF K-POINTS:',abs(kkk)-1
  write(6,*) 'TOTAL TIME:',dt4-dt0
  write(6,*) 'read:',cp(1)
  write(6,*) 'alm:',cp(2)
  write(6,*) 'hmsec:',cp(3)
  write(6,*) 'write:',cp(4)
  
  CALL ERRCLR(ERRFN)
  STOP 'LAPWSO END'
END PROGRAM lapwso

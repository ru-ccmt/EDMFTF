MODULE param
  INTEGER,PARAMETER :: LMAX=   12
  INTEGER           :: nato=    0
  INTEGER           :: ndif=    0
  INTEGER           :: nume=    0
  INTEGER           :: nmat=    0
  INTEGER,PARAMETER :: NRAD=  881
  ! LOMAX must be = LOMAX in LAPW1 otherwise conflict in INIT
  INTEGER,PARAMETER :: LOMAX=   3
  INTEGER,PARAMETER :: FLMAX=   3 
  INTEGER,PARAMETER :: LMAX2=  14
  ! LMX (not LMAX!) must be = LMAX in LAPW1 otherwise conflict in INIT
  INTEGER,PARAMETER :: LMX= LMAX+1
  INTEGER           :: LABC =    0
  INTEGER           :: LABC2=    0
  INTEGER,PARAMETER :: LMX2=  LMX*LMX
  INTEGER,PARAMETER :: MMAX= 2*LMAX+1
  INTEGER           :: nume2=    0
  INTEGER           :: num2=     0
  REAL*8,PARAMETER  :: CLIGHT=137.0359895D0
  INTEGER,PARAMETER :: nloat=3
  INTEGER           :: nsym=0
!rschmid
!  Extension of relativistic local orbitals.
!rschmid
  INTEGER,PARAMETER :: HBLOCK = 32
  INTEGER,PARAMETER :: NSLMAX=  5
END MODULE param

MODULE abcd
  COMPLEX*16,ALLOCATABLE  :: abcdlm(:,:,:,:,:)

  CONTAINS 
  SUBROUTINE init_abcd(labc2,ndif,nume)
    ALLOCATE(abcdlm(4,labc2,ndif,nume,2))
  END SUBROUTINE init_abcd
END MODULE abcd

MODULE ams
  REAL*8,ALLOCATABLE          :: atom_mass(:)

 CONTAINS
  SUBROUTINE init_ams
  REAL*8          :: a_m(103)
  ALLOCATE(atom_mass(103))
     DATA a_m /1.,4.,6.9,9.,10.8,12.,14.,16.,19.,20.2, &
          23.,24.3,27.,28.1,31.,32.,35.4,40.,39.1,40.,45., &
          47.9,50.9,52.,54.9,55.8,58.9,58.7,63.5,65.4,69.7, &
          72.6,74.9,79.,79.9,83.8,85.5,87.6,88.9,91.2,92.9, &
          95.9,98.,101.1,102.9,106.4,107.9,112.4,114.8, &
          118.7,121.8,127.6,126.9,131.3,132.9,137.3,138.9,140.1, &
          140.9,144.2,145.,150.4,152.,157.3,158.9,162.5, &
          164.9,167.3,168.9,173.,175.,178.5,180.9,183.8,186.2, &
          190.2,192.2,195.1,197.,200.6,204.4,207.2,209.,209., &
          210.,222.,223.,226.,227.,232.,231.,238.,237.,244.,243., &
          247.,247.,251.,252.,257.,258.,259.,262./     
   atom_mass(1:103) = a_m(1:103) 
   END SUBROUTINE init_ams
 END MODULE ams

MODULE couples
  COMPLEX*16,ALLOCATABLE  ::couplo(:,:,:,:,:,:)
 CONTAINS
 SUBROUTINE init_couplo(ndif,labc)
 ALLOCATE(couplo(ndif,labc,-labc:labc,-labc:labc,2,2))
 END SUBROUTINE init_couplo
END MODULE couples
	

MODULE hmsout
  INTEGER                :: neig
  REAL*8,ALLOCATABLE     :: en(:),vnorm(:,:)
  COMPLEX*16,ALLOCATABLE :: vect(:,:,:)

CONTAINS
  SUBROUTINE init_hmsout(nmat,nume2)
    ALLOCATE(en(nume2),vnorm(nume2,2))
    ALLOCATE(vect(nmat,nume2,2))
    en=0.0d0; vnorm=0.0d0
    vect=(0.0d0,0.0d0)
  END SUBROUTINE init_hmsout
END MODULE hmsout

MODULE orb
  INTEGER                :: nmod,nsp,natorb,nonso,iorbpot
  INTEGER,ALLOCATABLE    :: iat(:),iiat(:),ll(:)
  REAL*8                 :: Bext
  COMPLEX*16,ALLOCATABLE :: vv(:,:,:,:)

CONTAINS
  SUBROUTINE init_orb(nato,ndif)
    ALLOCATE(iat(3*ndif),iiat(3*ndif),ll(3*ndif))
    ALLOCATE(vv(3*ndif,-3:3,-3:3,3))
    iatom=0; nlorb=0; lorb=0; indorb=0
    vorb=(0.0d0,0.0d0)
  END SUBROUTINE init_orb
END MODULE orb

MODULE loabc
  REAL*8,ALLOCATABLE  :: alo(:,:,:,:),blo(:,:,:,:),clo(:,:,:,:)
  REAL*8,ALLOCATABLE  :: dplo(:,:,:),dpelo(:,:,:),elo(:,:,:,:)
  REAL*8,ALLOCATABLE  :: peilo(:,:,:),pelo(:,:,:)
  REAL*8,ALLOCATABLE  :: plo(:,:,:),pi12lo(:,:,:),pe12lo(:,:,:)
  
CONTAINS
  SUBROUTINE init_loabc(lomax,nloat,nato)
    ALLOCATE(alo(0:lomax,nloat,nato,2),blo(0:lomax,nloat,nato,2),clo(0:lomax,nloat,nato,2))
    ALLOCATE(dplo(0:lomax,nato,2),dpelo(0:lomax,nato,2),elo(0:lomax,nloat,nato,2))
    ALLOCATE(peilo(0:lomax,nato,2),pelo(0:lomax,nato,2))
    ALLOCATE(plo(0:lomax,nato,2),pi12lo(0:lomax,nato,2),pe12lo(0:lomax,nato,2))
  END SUBROUTINE init_loabc
END MODULE loabc

MODULE loabcr
  real*8,ALLOCATABLE  :: alor(:,:,:), blor(:,:,:)
  real*8,ALLOCATABLE  :: clor(:,:,:), dplor(:,:,:)
  real*8,ALLOCATABLE  :: dpelor(:,:,:), elor(:,:,:)
  real*8,ALLOCATABLE  :: peilor(:,:,:), pelor(:,:,:)
  real*8,ALLOCATABLE  :: plor(:,:,:)
  real*8,ALLOCATABLE  :: pi2lor(:,:,:),pe2lor(:,:,:)
  real*8,ALLOCATABLE  :: elor2(:,:,:)

CONTAINS
  SUBROUTINE init_loabcr(lomax,nloat,nato)
    ALLOCATE(alor(0:lomax,nato,2), blor(0:lomax,nato,2))
    ALLOCATE(clor(0:lomax,nato,2), dplor(0:lomax,nato,2))
    ALLOCATE(dpelor(0:lomax,nato,2), elor(0:lomax,nato,2))
    ALLOCATE(peilor(0:lomax,nato,2), pelor(0:lomax,nato,2))
    ALLOCATE(plor(0:lomax,nato,2))
    ALLOCATE(pi2lor(0:lomax,nato,2),pe2lor(0:lomax,nato,2))
    ALLOCATE(elor2(0:lomax,nato,2))
  END SUBROUTINE init_loabcr
END MODULE loabcr

MODULE lolog
  INTEGER,ALLOCATABLE :: nlo(:),nlov(:),nlon(:),ilo(:,:),mrf(:,:)
  LOGICAL,ALLOCATABLE :: loor(:,:),lapw(:,:),lso(:)

CONTAINS
  SUBROUTINE init_lolog(nato,lomax,labc)
    ALLOCATE(nlo(nato),nlov(nato),nlon(nato),ilo(0:lomax,nato),mrf(0:labc,nato))
    ALLOCATE(loor(0:lomax,nato),lapw(0:labc,nato),lso(nato))
    lso=.true.
  END SUBROUTINE init_lolog
END MODULE lolog

MODULE radovlp
  real*8,ALLOCATABLE :: rup  (:,:,:,:)
  real*8,ALLOCATABLE :: rudp (:,:,:,:)
  real*8,ALLOCATABLE :: ru2ud(:,:,:,:)
  real*8,ALLOCATABLE :: ru2p (:,:,:,:)
  real*8,ALLOCATABLE :: ru2u (:,:,:,:)
  real*8,ALLOCATABLE :: rrr  (:,:,:,:)
  real*8,ALLOCATABLE :: rupd (:,:,:,:)
  real*8,ALLOCATABLE :: rudpd(:,:,:,:)
  real*8,ALLOCATABLE :: ru2pd(:,:,:,:)
  real*8,ALLOCATABLE :: rpdpd(:,:,:,:)
  real*8,ALLOCATABLE :: rppd (:,:,:,:) 

CONTAINS
  SUBROUTINE init_radovlp(lmax,nato)
    ALLOCATE(rup  (0:lmax,nato,2,2))
    ALLOCATE(rudp (0:lmax,nato,2,2))
    ALLOCATE(ru2ud(0:lmax,nato,2,2))
    ALLOCATE(ru2p (0:lmax,nato,2,2))
    ALLOCATE(ru2u (0:lmax,nato,2,2))
    ALLOCATE(rrr  (0:lmax,nato,2,2))
    ALLOCATE(rupd (0:lmax,nato,2,2))
    ALLOCATE(rudpd(0:lmax,nato,2,2))
    ALLOCATE(ru2pd(0:lmax,nato,2,2))
    ALLOCATE(rpdpd(0:lmax,nato,2,2))
    ALLOCATE(rppd (0:lmax,nato,2,2))
  END SUBROUTINE init_radovlp
END MODULE radovlp

MODULE rlolog
  INTEGER              :: nnrlo
  integer,ALLOCATABLE  :: nrlo(:),nrlov(:),nrlon(:)
  logical,ALLOCATABLE  :: loorext(:,:)

CONTAINS
  SUBROUTINE init_rlolog(nato,lomax)
    ALLOCATE(nrlo(nato),nrlov(nato),nrlon(nato))
    ALLOCATE(loorext(0:lomax,nato))
  END SUBROUTINE init_rlolog
END MODULE rlolog

MODULE rpars
  INTEGER,allocatable :: extl(:,:),nlr(:)
  REAL*8,ALLOCATABLE  :: extei(:,:),extde(:,:)
  CHARACTER*4,allocatable  :: extem(:,:)

CONTAINS
  SUBROUTINE init_rpars(nato,lomax)
    ALLOCATE(extl(nato, lomax),nlr(nato))
    ALLOCATE(extei(nato, lomax),extde(nato, lomax))
    ALLOCATE(extem(nato, lomax))
  END SUBROUTINE init_rpars
END MODULE rpars

MODULE rotmat
  REAL*8,ALLOCATABLE       :: ROTIJ(:,:,:),TAUIJ(:,:),det(:),phase(:)

  CONTAINS
    SUBROUTINE init_rotmat(ndif)
      ALLOCATE(rotij(3,3,ndif),tauij(3,ndif),det(ndif),phase(ndif))
    END SUBROUTINE init_rotmat
END MODULE rotmat

MODULE struct
  USE param
  REAL*8                   :: AA,BB,CC,VOL,pia(3),alpha(3)
  REAL*8,ALLOCATABLE       :: R0(:),DX(:),RMT(:),zz(:),rotloc(:,:,:)
  REAL*8,ALLOCATABLE       :: tau(:,:)
  REAL*8,POINTER           :: pos(:,:)
  CHARACTER*4              :: lattic,irel,cform
  CHARACTER*80             :: title
  CHARACTER*10,ALLOCATABLE :: aname(:)
  INTEGER                  :: nat,iord
  INTEGER,ALLOCATABLE      :: mult(:),jrj(:),iatnr(:),isplit(:)
  INTEGER,ALLOCATABLE      :: iz(:,:,:),inum(:)

CONTAINS
  SUBROUTINE init_struct
    USE reallocate
    IMPLICIT NONE

    INTEGER                :: ios
    REAL*8                 :: test,ninety
!loop indexs
    INTEGER                :: index,i,j,j1,j2,m,jatom

    test=1.D-5
    ninety=90.0D0

    read (20,1000) title
    read (20,1010)  lattic,nat,cform,irel
    nato=nat
    ALLOCATE(aname(nato),mult(0:nato),jrj(nato),r0(nato),dx(nato),rmt(nato),zz(nato),rotloc(3,3,nato),iatnr(nato),isplit(nato))
    ALLOCATE (pos(3,48*nat))
    mult(0)=0
    read (20,1020) aa,bb,cc,alpha(1),alpha(2),alpha(3)
    IF(ABS(ALPHA(1)).LT.test) ALPHA(1)=ninety
    IF(ABS(ALPHA(2)).LT.test) ALPHA(2)=ninety
    IF(ABS(ALPHA(3)).LT.test) ALPHA(3)=ninety
    INDEX=0
    DO jatom=1,NAT
       INDEX=INDEX+1
       READ(20,1030,iostat=ios) iatnr(jatom),( pos(j,index),j=1,3 ), &
            mult(jatom),isplit(jatom) 
       IF(ios /= 0) THEN
          WRITE(6,*) iatnr(jatom),( pos(j,index),j=1,3 ), &
               mult(jatom),isplit(jatom) 
          WRITE(6,*) 'ERROR IN STRUCT FILE READ'
          STOP
       ENDIF
       IF (mult(jatom) .EQ. 0) THEN
          WRITE (6,6000) jatom, index, mult(jatom)
          STOP
       ENDIF
       DO m=1,mult(jatom)-1                                     
          index=index+1                                            
          READ(20,1031) iatnr(jatom),( pos(j,index),j=1,3)         
       ENDDO
       READ(20,1050) aname(jatom),jrj(jatom),r0(jatom),rmt(jatom), &
            zz(jatom)
       dx(jatom)=LOG(rmt(jatom)/r0(jatom)) / (jrj(jatom)-1)           
       rmt(jatom)=r0(jatom)*EXP( dx(jatom)*(jrj(jatom)-1) )           
       READ(20,1051) ((rotloc(i,j,jatom),i=1,3),j=1,3)                
    ENDDO
    ndif=index
    CALL doreallocate(pos, 3, ndif)
    READ(20,1151) iord
    nsym=iord
    ALLOCATE(iz(3,3,nsym),tau(3,nsym),inum(nsym))
    DO j=1,iord
       READ(20,1101) ( (iz(j1,j2,j),j1=1,3),tau(j2,j),j2=1,3 ),inum(j)
    ENDDO
  
1000 FORMAT(A80)                                                       
1010 FORMAT(A4,23X,I3,1x,a4,/,13X,A4,18X,A4)                                 
1020 FORMAT(6F10.7,10X,F10.7)                                          
1030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
1031 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
1051 FORMAT(20X,3F10.8)                                             
1101 FORMAT(3(3I2,F11.8/),I8)
1151 FORMAT(I4)
6000 FORMAT(///,3X,'ERROR IN LAPW0 : MULT(JATOM)=0 ...', &
          /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
  END SUBROUTINE init_struct

END MODULE struct

MODULE vns
  LOGICAL,ALLOCATABLE   :: lvns(:)
  INTEGER,ALLOCATABLE   :: nvns(:),mvns(:,:)
  REAL*8,ALLOCATABLE    :: vaa(:,:,:),vab(:,:,:),vad(:,:,:), &
                           vbb(:,:,:),vbd(:,:,:),vdd(:,:,:)

CONTAINS
  SUBROUTINE init_vns(nato)
    ALLOCATE(lvns(nato),nvns(nato),mvns(5,nato))
    ALLOCATE(vaa(nato,5,2),vab(nato,5,2),vad(nato,5,2), &
             vbb(nato,5,2),vbd(nato,5,2),vdd(nato,5,2))
    lvns=.FALSE.
  END SUBROUTINE init_vns
END MODULE vns

MODULE peic
  REAL*8,ALLOCATABLE    :: pei(:,:,:)

CONTAINS
  SUBROUTINE init_peic(lmx,nato)
    ALLOCATE(pei(lmx,nato,2))
  END SUBROUTINE init_peic
END MODULE peic

MODULE hexpt
  REAL*8,ALLOCATABLE    :: hexl(:,:,:,:)

CONTAINS
  SUBROUTINE init_hexpt(lomax,nato)
    ALLOCATE(hexl(0:lomax,nato,2,9))
    hexl = 0.d0
  END SUBROUTINE init_hexpt
END MODULE hexpt

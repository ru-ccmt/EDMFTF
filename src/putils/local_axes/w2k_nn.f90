!program nn
!  REAL*8 :: dlimit, dstmax, dfac
!  CHARACTER*80 :: case
!  
!  dlimit=0.d0 
!  dstmax=-1   
!  dfac=2      
!  case='Sr2IrO4'
!  call w2k_nn(case,dfac,dlimit,dstmax)
!  
!end program nn


module variable_fields
  real*8, allocatable :: pos(:,:)
end module variable_fields

SUBROUTINE w2knn(case,dfac_,dlimit_,dstmax_)
  use variable_fields
  IMPLICIT REAL*8 (A-H,O-Z)
  !     Constant parameter definition
  PARAMETER (NNN=100000)                                              
  PARAMETER (mshell= 100)
  !
  CHARACTER*80, intent(in) :: case
  REAL*8, intent(in)       :: dfac_, dlimit_,dstmax_
  !
  !f2py real*8  optional,intent(in)        :: dfac_ = 2.0
  !f2py real*8  optional,intent(in)        :: dlimit_ = 0.0
  !f2py real*8  optional,intent(in)        :: dstmax_ = -1
  !
  !     NATO    NUMBER OF  NONEQUIVALENT ATOMS            
  !     NNN     NUMBER OF NEIGHBOURS                                      
  !     NDIF    NUMBER OF ALL ATOMS IN UNITCELL                           
  !     MSHELL  Max. Nr. of SHELLS
  !
  ! from struk
  real*8            BR2(3,3)
  real*8            BR1(3,3)
  CHARACTER*4       LATTIC                                            
  CHARACTER*10, allocatable :: name(:)
  ! 
  LOGICAL           ORTHO,error
  LOGICAL           SHOWALL
  CHARACTER*10      KNAME
  CHARACTER*11      STATUS,FORM                                     
  CHARACTER*79      TITLE                                           
  CHARACTER*100     FNAME,fname1                                           
  !-----------------------------------------------------------------------
  real*8  A(3),PIA(3),VI   
  DIMENSION DIF(3),XX(3),PP(3),P(3),DISTS(NNN),NR(NNN),PNN(3,NNN)   
  DIMENSION NNAT(NNN),help(3)                                            
  
  LOGICAL, allocatable :: ovlap(:),used(:)
  real*8,  allocatable :: rmt(:),zz(:),shdist(:,:)
  real*8,  allocatable :: zzorig(:),shellz(:,:,:)
  real*8,  allocatable :: V(:),POSN(:,:)
  real*8,  allocatable :: R0(:),R0N(:),DH(:)  
  real*8,  allocatable :: rmtt(:),rmttn(:),zzn(:),zzo(:) 
  real*8,  allocatable :: ROTLOC(:,:,:),ROTIJ(:,:,:),TAUIJ(:,:)  
  integer,  allocatable :: JRJ(:),JRJN(:)
  integer,  allocatable :: IATNR(:),MULT(:),icnt(:,:,:) 
  integer,  allocatable :: iz(:,:),ishell(:),ityp(:),imult(:)
  CHARACTER*10, allocatable :: namen(:)
  CHARACTER*100 :: f_20, f_66, f_67, f_21, f_68, f_06
  !-----------------------------------------------------------------------
  logical there
  dfac=dfac_
  dlimit=dlimit_
  dstmax=dstmax_
  
  ! 20,'Sr2IrO4.struct',   'old',    'formatted',0
  ! 66,'Sr2IrO4.outputnn', 'unknown','formatted',0
  ! 67,'Sr2IrO4.nnshells', 'unknown','formatted',0
  ! 68,'Sr2IrO4.rotlm_',   'unknown','formatted',0
  ! 21,'Sr2IrO4.struct_nn','unknown','formatted',0
  f_20 = TRIM(ADJUSTL(case))//'.struct'
  f_66 = TRIM(ADJUSTL(case))//'.outputnn_'
  f_67 = TRIM(ADJUSTL(case))//'.nnshells_'
  f_21 = TRIM(ADJUSTL(case))//'.struct_nn'
  f_68 = TRIM(ADJUSTL(case))//'.rotlm_'
  f_06 = TRIM(ADJUSTL(case))//'.info_'

  WRITE(6,*) 'structure file=', f_90
  
  open(20, FILE=f_20, STATUS='old', FORM='formatted')
  open(66, FILE=f_66, STATUS='unknown', FORM='formatted')
  open(67, FILE=f_67, STATUS='unknown', FORM='formatted')

  open(6, FILE=f_06, STATUS='unknown', FORM='formatted')
  
  fname1 = f_20
  
  !     This should be how much larger (in 1D) the DFT distances are relative to true
  DFTERR = 1.01D0        !This is about right for PBE -- could look in case.in0?
  !
  !     See if the user has calibrated the lattice parameter scaling difference
  !inquire(file='.latcalib',exist=there)
  !if(there)then
  !   open(unit=99,file='.latcalib')
  !   !       This should be how much larger (in 1D) the DFT distances are relative to true
  !   read(99,*)DFTERR
  !   close(unit=99)
  !endif
  SHOWALL=.TRUE.
  
  if(dfac .lt. 0.)then
     SHOWALL=.FALSE.
     dfac=abs(dfac)
  endif
  if(dlimit.lt.1.d-7) dlimit=1.d-5
  ishellmax=99
  
  RTB2=SQRT(3.)/2.                                                  
  !                                                                       
  !                                                                       
  !.....START READING FILE STRUCT, writing new struct (21)                   
  READ(20,1510) TITLE                                               
  write(66,1510) TITLE
  READ(20,1511) LATTIC,NAT,title                                          
  write(66,1511) LATTIC,NAT,title                                          
  !    allocate nato-arrays
  allocate (  name(nat) )
  allocate ( ovlap(nat) )
  allocate ( rmt(nat),zz(nat) )
  allocate ( zzorig(nat),V(nat),R0(nat),DH(nat),rmtt(nat) )
  allocate ( ROTLOC(3,3,nat),JRJ(nat),IATNR(nat),MULT(nat) ) 
  allocate ( POS(3,48*nat*2) )
  !     READ IN LATTICE CONSTANTS                                         
  read(20,'(6F10.6)') A(1),A(2),A(3),alpha,beta,gamma
  if(alpha.eq.0.d0) alpha=90.d0                                       
  if(beta .eq.0.d0) beta =90.d0                                       
  if(gamma.eq.0.d0) gamma=90.d0                                       
  write(66,'(6F10.6)') A(1),A(2),A(3),alpha,beta,gamma
  !17 FORMAT(6F10.6)                                                    
  if(dstmax .lt. 0)dstmax=max(20.d0,0.61d0*max(a(1),a(2),a(3)))  
  if(nat .eq. 1) then
     dstmax=max(20.d0,a(1),a(2),a(3))*1.1
  else
     dstmax=min(40.d0,dstmax)
  endif

  WRITE(6,*) 'A[1]=', A(1), 'A[2]=', A(2), 'A[3]=', A(3)

  WRITE(6,*) 'DSTMAX:',dstmax                  
  iix=5
  iiy=5
  iiz=5
  if(lattic(1:1).eq.'P'.or.lattic(1:1).eq.'H') then
     iix=max(1,nint(dstmax/a(1)))+1
     iiy=max(1,nint(dstmax/a(2)))+1
     iiz=max(1,nint(dstmax/a(3)))+1
  endif
  
  DO i=1,1000
     if(nat .gt. 1)then
        if( (iix-1)*a(1) .gt. dstmax*2.d0)then
           iix=max(1,iix-1)
        else if( (iiy-1)*a(2) .gt. dstmax*2.d0)then
           iiy=max(1,iiy-1)
        else if( (iiz-1)*a(3) .gt. dstmax*2.d0)then
           iiz=max(1,iiz-1)
           EXIT
        endif
     endif
  ENDDO
  
  if(lattic(1:3).eq.'CXY') then     !fix for orthorombic CXY with very different a,b
     iix=max(iix,iiy)
     iiy=iix
  endif
  WRITE(6,*) 'iix,iiy,iiz',iix,iiy,iiz,a(1)*iix,a(2)*iiy,a(3)*iiz

  !     INDEX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE    
  !     NONEQUIVALENT ATOMS                                               
  INDEX=0                                                           
  DO JATOM = 1,NAT        ! 50
     INDEX=INDEX+1                                                  
     READ(20,1012) IATNR(JATOM),( POS(J,INDEX),J=1,3 ),MULT(JATOM)  
     write(66,1012) IATNR(JATOM),( POS(J,INDEX),J=1,3 ),MULT(JATOM)  
     IF (MULT(JATOM).EQ.0) THEN                                     
        write(66,1020) JATOM,INDEX,MULT(JATOM)                       
        STOP ' NNN: MULT EQ 0'                                      
     ENDIF
     DO M=1,MULT(JATOM)-1 
        INDEX=INDEX+1                                            
        READ(20,1011) IATNR(JATOM),( POS(J,INDEX),J=1,3)         
        write(66,1011) IATNR(JATOM),( POS(J,INDEX),J=1,3)         
     ENDDO

     READ(20,1050) NAME(JATOM),JRJ(JATOM),R0(JATOM),RMTT(jatom),  zz(jatom)  
     zzorig(jatom)=zz(jatom)      
     if(name(jatom)(3:3).ne.' ') then
        write(66,*) 'NAMED ATOM: ',name(jatom), 'Z changed to ', 'IATNR+999 to determine equivalency'
        write(*,*) 'NAMED ATOM: ',name(jatom), 'Z changed to ', 'IATNR+999 to determine equivalency'
        zz(jatom)=999+jatom
     endif
     write(66,1049) NAME(JATOM),JRJ(JATOM),R0(JATOM),RMTT(jatom), zzorig(jatom)        
     if((jrj(jatom)/2)*2.eq.jrj(jatom)) then
        write(*,*) 'WARNING: JRJ of atom',jatom,' is even:',jrj(jatom)
        write(*,*) 'CHANGE it to ODD number !!!!'
        write(66,*) 'WARNING: JRJ of atom',jatom,' is even:', jrj(jatom)
        write(66,*) 'CHANGE it to ODD number !!!!'
     endif
     DH(JATOM)=LOG(RMTT(jatom)/R0(JATOM)) / (JRJ(JATOM)-1)                 
     RMT(JATOM)=R0(JATOM)*EXP( DH(JATOM)*(JRJ(JATOM)-1) )           
     READ(20,1051) ((ROTLOC(I1,I2,JATOM),I1=1,3),I2=1,3)               
     write(66,1051) ((ROTLOC(I1,I2,JATOM),I1=1,3),I2=1,3)               
  ENDDO !50   CONTINUE                                                          
  
  call reduce_alloc_r8_2d(3,48*nat,3,index)
  allocate ( shdist(index,mshell),shellz(index,mshell,48) )
  allocate ( namen(index),POSN(3,index),R0N(index),rmttn(index) )
  allocate ( zzn(index),zzo(index),ROTIJ(3,3,index),TAUIJ(3,index) )
  allocate ( JRJN(index),icnt(index,mshell,48),iz(index,mshell) )   
  allocate ( used(nat*index),ishell(nat*index),ityp(nat*index) )
  allocate ( imult(nat*index) )
  !                                                                       
  !                                                                       
  !.....SET UP LATTICE, AND IF REQUIRED ROTATION MATRICES                 
  CALL DIRLAT (BR1,BR2,ortho,lattic,nat,alpha,beta,gamma,A)
  !                                           
  write(66,*) ' '
  write(66,*) 'Bond-Valence Sums are calculated for current lattice parameters'
  write(66,*) 'and rescaled ones by 1 %. (You can put scaling into  .latcalib)' 
  pi=4.d0*atan(1.d0)
  cosgam=cos(gamma/180.d0*pi)
  singam=sin(gamma/180.d0*pi)           
  INDEX=0                                                           
  DO JATOM=1,NAT       ! 200
     do i=1,nat
        ovlap(i)=.true.
     enddo
     DO M=1,MULT(JATOM)  ! 190
        INDEX=INDEX+1                                                     
        DO J=1,3                                                      
           XX(J)=POS(J,INDEX)
        ENDDO
if(SHOWALL.or.(M.EQ.1)) write(66,'(/,A,I3,2X,A,I3,2X,A10,A,3F10.5)') ' ATOM:',JATOM,'EQUIV.',M,NAME(JATOM),' AT',XX(1),XX(2),XX(3)
        NC=0             
        DO I1=-iix,iix              ! 180
           DO I2=-iiy,iiy           ! 180
              DO I3=-iiz,iiz        ! 180
                 IF(ortho) THEN                                   
                    P(1)=I1*BR2(1,1)+I2*BR2(2,1)+I3*BR2(3,1)                    
                    P(2)=I1*BR2(1,2)+I2*BR2(2,2)+I3*BR2(3,2)                    
                    P(3)=I1*BR2(1,3)+I2*BR2(2,3)+I3*BR2(3,3)                    
                 ELSE                                                        
                    P(1)=I1                                                     
                    P(2)=I2                                                     
                    P(3)=I3                                                     
                    IF(LATTIC(1:3).eq.'CXZ') THEN
                       P(1)=I1*0.5d0+i3*0.5d0
                       P(2)=I2
                       P(3)=-I1*0.5d0+i3*0.5d0
                    END IF
                 ENDIF
                 K=0                                                           
                 DO JAT=1,NAT           ! 120
                    DO MM=1,MULT(JAT)   ! 110
                       K=K+1                                                             
                       DIST=0.                                                           
                       DO L=1,3                                                      
                          PP(L)=POS(L,K)+P(L)
                          DIF(L)=XX(L)-PP(L)
                       ENDDO
                       
                       IF (.not.ortho) THEN                                       
                          help(1)=dif(1)  
                          help(2)=dif(2)  
                          help(3)=dif(3)  
                          if(lattic(1:1).eq.'R') then
                             dif(1)=help(1)*BR2(1,1)+help(2)*BR2(2,1)+help(3)*BR2(3,1)             
                             dif(2)=help(1)*BR2(1,2)+help(2)*BR2(2,2)+help(3)*BR2(3,2)             
                             dif(3)=help(1)*BR2(1,3)+help(2)*BR2(2,3)+help(3)*BR2(3,3)           
                          elseif(lattic(1:3).eq.'CXZ') then
                             dif(1)=help(1)*singam            
                             dif(2)=(help(1)*cosgam*a(1)+help(2)*a(2))/a(2)             
                             dif(3)=help(3)           
                          else
                             dif(1)=(help(1)*BR2(1,1)*a(1)+help(2)*BR2(2,1)*a(2)+help(3)*BR2(3,1)*a(3))/a(1)
                             dif(2)=(help(1)*BR2(1,2)*a(1)+help(2)*BR2(2,2)*a(2)+help(3)*BR2(3,2)*a(3))/a(2)
                             dif(3)=(help(1)*BR2(1,3)*a(1)+help(2)*BR2(2,3)*a(2)+help(3)*BR2(3,3)*a(3))/a(3)
                          endif
                       ENDIF
                       DO L=1,3                                                      
                          DIST=DIST+DIF(L)*DIF(L)*A(L)*A(L)
                       ENDDO
                       DIST=SQRT(DIST)
                       
                       IF(DIST.GT.dstmax) CYCLE
                       IF(DIST.LT..001) CYCLE
                       NC=NC+1         
                       if(nc.gt.nnn) stop ' nnn too small'                                     
                       DISTS(NC)=DIST                                                    
                       NNAT(NC)=JAT                                                      
                       DO L=1,3                                                      
                          PNN(L,NC)=PP(L)
                       ENDDO
                    ENDDO !110 CONTINUE
                 ENDDO !120 CONTINUE                                                          
              ENDDO !  180 CONTINUE
           ENDDO    ! 180
        ENDDO       ! 180
        CALL ORD2(DISTS,NR,NC)                                            
        N1=1                                                              
        N2=NR(N1)                                                         
        N3=NNAT(N2)             
        SUMRAD=RMT(JATOM)+RMT(N3)                                         
        IF(M.EQ.1) THEN                                                   
           IF(SUMRAD.LE.DISTS(N1))THEN                                       
              WRITE(*,'(/,3X,A,I3,2X,A10,A,I3,2X,A10)') ' ATOM',JATOM,NAME(JATOM),' ATOM',N3,NAME(N3)
    WRITE(*, '(A,I3,A,F7.5,A,I3,A,F7.5,/,A,F8.5,2X,A,F8.5)') ' RMT(',JATOM,')=',RMT(JATOM),' AND RMT(',N3,')=',RMT(N3),' &
      & SUMS TO',SUMRAD,'LT.  NN-DIST=',DISTS(1)
              WRITE(66,'(A,I3,A,F7.5,A,I3,A,F7.5,/,A,F8.5,2X,A,F8.5)') ' RMT(',JATOM,')=',RMT(JATOM),' & 
         &AND RMT(',N3,')=',RMT(N3),' SUMS TO',SUMRAD,'LT.  NN-DIST=',DISTS(1)
           ELSE                                                              
              ovlap(n3)=.false.                        
              write(66,4) JATOM,RMT(JATOM),N3,RMT(N3), SUMRAD, DISTS(N1)
4             FORMAT(/,'   ERROR !!!!!!!!!!!!!!!', /,' RMT(',I3,')=',F7.5,' AND RMT(',I3,')=',F7.5,/,' &
              &SUMS TO',F8.5,' GT NNN-DIST=',F8.5)
              WRITE(*,4) JATOM,RMT(JATOM),N3,RMT(N3), SUMRAD, DISTS(N1)         
           ENDIF
        ENDIF
        !.....determination of "equal" atoms
        !
        olddis=0.d0
        ishell(index)=0
        !
        bva=0.0
        bvaerr=0.0

        !!!! XX(:) = POS(:,INDEX)
        !!!! N1=1,...NC
        !!!! DISTS(N1) -- distance
        !!!! NNAT(NC)  -- isort
        !!!! PNN(:3, N2=NR(N1) ) -- neigbor
        
        DO N1=1,NC
           N2=NR(N1)                                                         
           N3=NNAT(N2)                                                       
           !
           !     Include BVA
           i1=nint(zzorig(jatom))
           i2=nint(zzorig(N3))

           !! You can restore this setting this statement to .true.
           if (.false.) then
              call bvan(i1,i2,scale,dbva)
              if(scale .gt. 0)then
                 val=exp( (dbva - dists(n1)*0.529177) /scale)
                 bva=bva+val
                 valerr=exp( (dbva - dists(n1)*0.529177/DFTERR) /scale)
                 bvaerr=bvaerr+valerr
              endif
           endif
              
           !
           SUMRAD=RMT(JATOM)+RMT(N3)                                        
           if(dists(n1).lt.dfac*dists(1)) then
              if(SHOWALL.or.(M.EQ.1)) write(66,3) N3,NAME(N3),(PNN(L,N2),L=1,3),DISTS(N1),DISTS(N1)*0.529177  
              IF(ovlap(n3).and.SUMRAD.GE.DISTS(N1)) THEN
                 ovlap(n3)=.false.                        
                 write(66,4) JATOM,RMT(JATOM),N3,RMT(N3), SUMRAD, DISTS(N1)         
                 WRITE(*,4) JATOM,RMT(JATOM),N3,RMT(N3), SUMRAD, DISTS(N1)         
              end if
           end if
           !
           if((dists(n1)-olddis).gt.dlimit) then
              !.....new shell
              ishell(index)=ishell(index)+1
              iz(index,ishell(index))=1
              olddis=dists(n1)  
              if(ishell(index).ge.mshell) EXIT !goto 187
              icnt(index,ishell(index),iz(index,ishell(index)))=1
              shellz(index,ishell(index),iz(index,ishell(index)))=zz(n3)
              shdist(index,ishell(index))=olddis
           else
              !.....old shell
              do i=1,iz(index,ishell(index))
                 if(zz(n3).eq.shellz(index,ishell(index),i)) then
                    icnt(index,ishell(index),i)=icnt(index,ishell(index),i)+1
                    EXIT !goto 186
                 endif
              enddo
              iz(index,ishell(index))=iz(index,ishell(index))+1
              shellz(index,ishell(index),iz(index,ishell(index)))=zz(n3)
              icnt(index,ishell(index),iz(index,ishell(index)))=1          
           endif
           !186         continue
           !
        ENDDO !185 CONTINUE
        !187      CONTINUE
        if( (SHOWALL.or.(M.EQ.1)).and.(bva .gt. 0.001D0) )then
           write(66,'(A,i3,A,i2,A,A,A,2F8.2)') 'Atom ',JATOM,' equiv ',M,' ',NAME(JATOM),' Bond-Valence Sum ',bva,bvaerr
        endif
        !
        !.....limit shells to some maximum
        write(67,*) '---------------- ATOM',index
        do i=1,ishell(index)
           if(shdist(index,i).gt.dstmax-dstmax/3.d0) then
              ishell(index)=i
              !goto 190
              EXIT
           endif
           write(67,*) ' SHELL:',i,shdist(index,i)
           do j=1,iz(index,i)
              write(67,*) icnt(index,i,j),shellz(index,i,j)
           enddo
        enddo
        !
     ENDDO !  190 CONTINUE
  ENDDO !200 CONTINUE                                                          
  !                      
  write(66,552)
  INDEX=0                                                           
  DO  JATOM=1,NAT                
     DO  M=1,MULT(JATOM)                                            
        INDEX=INDEX+1
        zzo(index)=zz(jatom)                                                  
        write(66,553)index,((shdist(index,i),icnt(index,i,j),shellz(index,i,j),j=1,iz(index,i)),i=1,4)
     enddo
  enddo
  write(66,*)
  !
  inat=0
  do i=1,index
     ityp(i)=i
     imult(i)=0
     used(i)=.false.  
  enddo
  write(66,554)
  ityp1=0    
  ij=0          
  i=0                           
  do i0=1,nat            ! 500
     do i00=1,mult(i0)   ! 500
        i=i+1
        if(used(i)) CYCLE !goto 500
        write(66,*)
        do j=i,index     ! 501
           if(zz(i0).ne.zzo(j)) CYCLE !goto 501
           !     compare atom with index i with shells of other atoms (Fuhr)
           if (ishell(i).ne.ishell(j)) then
              write(66,559) i,j,ishell(i),ishell(j)
              CYCLE !goto 501
           endif
           !     compare atom with index i with all other atoms
           do i1=1,ishell(i)-1   ! 510
              if(abs(shdist(i,i1)-shdist(j,i1)).gt.dlimit.and.shdist(i,i1).lt.2.d000*max(a(1),a(2),a(3))) then
                 write(66,550) i,j,i1,shdist(i,i1),shdist(j,i1)
                 goto 501
              endif
              do i2=1,iz(i,i1)   ! 511
                 do j2=1,iz(j,i1)
                    if((icnt(i,i1,i2).eq.icnt(j,i1,j2)).and.(shellz(i,i1,i2).eq.shellz(j,i1,j2))) goto 511
                    if(shdist(i,i1).gt.2.d0000*max(a(1),a(2),a(3))) goto 511
                 enddo
                 write(66,551) i,j,i1,icnt(i,i1,i2),icnt(j,i1,j2-1),shellz(i,i1,i2),shellz(j,i1,j2-1),shdist(i,i1)
                 goto 501
511              CONTINUE
              enddo !511              continue
           enddo       !510           continue
           write(66,555) i,j
           if(i.eq.j) then
              ityp1=ityp1+1
              namen(ityp1)=name(i0)
              JRJN(ityp1)=jrj(i0)
              R0N(ityp1)=r0(i0)
              RMTTN(ityp1)=rmtt(i0)
              zzn(ityp1)=zzorig(i0)
           endif
           ityp(j)=ityp1
           if(inat.lt.ityp(j)) inat=ityp(j)
           imult(ityp1)=imult(ityp1)+1
           used(j)=.true.
           ij=ij+1
           posn(1,ij)=pos(1,j)
           posn(2,ij)=pos(2,j)
           posn(3,ij)=pos(3,j)
501        continue
        ENDDO
     ENDDO ! 500
  ENDDO ! 500

  write(66,*)      
  error=.false.
  INDEX=0                                                           
  DO  JATOM=1,NAT                
     write(66,556) jatom,mult(jatom),imult(jatom)  
     if(mult(jatom).ne.imult(jatom)) then
        error=.true.
        write(66,*)'WARNING: MULT not equal. The new multiplicity is',' different from the old one'  
        write(6,*)'WARNING: Mult not equal. PLEASE CHECK outputnn-file'
     end if
     DO  M=1,MULT(JATOM)                                            
        INDEX=INDEX+1                                                     
        if(jatom.ne.ityp(index)) then
           error=.true.
           write(66,557) index,jatom,ityp(index)
           write(66,*) 'WARNING: ITYP not equal. The new type is',' different from the old one'  
           write(6,*)'WARNING: ityp not equal. PLEASE CHECK outputnn-file'
        endif
     enddo
  enddo
  !

  open(68, FILE=f_68, STATUS='unknown', FORM='formatted')
  WRITE(68,*) 'BR1'
  DO JR=1,3
     WRITE(68, '(3F15.10)') BR1(1,JR), BR1(2,JR), BR1(3,JR) !--- Writting BR1 ----!
  ENDDO
  !WRITE(68,*) 'BR2'
  !DO JR=1,3
  !   WRITE(68, '(3F15.10)') BR2(1,JR), BR2(2,JR), BR2(3,JR) !--- Writting BR1 ----!
  !ENDDO
  
  close(20)
  close(66)
  close(67)
  close(68)
  if(.not.error) return !stop 'NN ENDS'
  
  write(66,*)
  write(66,*) 'NEW LIST of EQUIVALENT POSITIONS written to',' case.struct_nn'   
  !.....write new struct file
  !

  open(21, FILE=f_21, STATUS='unknown', FORM='formatted')
  write(21,1510) TITLE
  write(21,1511) LATTIC,inat,title                                  
  write(21,'(6F10.6)') A(1),A(2),A(3),alpha,beta,gamma
  index=0
  do jatom=1,inat
     INDEX=INDEX+1                                                     
     write(21,1012) -JATOM,( POSN(J,INDEX),J=1,3 ),iMULT(JATOM),8
     write(66,*)
     write(66,1011) -JATOM,( POSN(J,INDEX),J=1,3 ),NAMEN(JATOM),zzn(jatom) 
     DO M=1,iMULT(JATOM)-1                                     
        INDEX=INDEX+1                                            
        write(21,1011) -JATOM,( POSN(J,INDEX),J=1,3)         
        write(66,1011) -JATOM,( POSN(J,INDEX),J=1,3)         
     ENDDO
     if(jatom.lt.10) then
        write(21,1149) NAMEN(JATOM)(1:2),jatom,JRJN(JATOM),R0N(JATOM),RMTTN(jatom),zzn(jatom)       
     else if(jatom.lt.100) then
        write(21,1148) NAMEN(JATOM)(1:2),jatom,JRJN(JATOM),R0N(JATOM),RMTTN(jatom),zzn(jatom)       
     else
        write(21,1147) NAMEN(JATOM)(1:2),jatom,JRJN(JATOM),R0N(JATOM),RMTTN(jatom),zzn(jatom)       
     endif
     write(21,1512) 1.0,0.0,0.0
     write(21,1512) 0.0,1.0,0.0
     write(21,1512) 0.0,0.0,1.0
  enddo
  write(21,*) 0
  !                 
  WRITE(6,*) ' '
  WRITE(6,*) " NN created a new ",trim(fname1),"_nn file"   
                                                     
  close(21)

!3 FORMAT(' ATOM:',I3,2X,A10,'AT',3F8.4,' IS',F9.5,' A.U.',f10.5,' ANG')
3 FORMAT(' ATOM:',I3,2X,A10,'AT',3F14.10,' IS',F9.5,' A.U.',f10.5,' ANG')

550 format(' atom:',i4,' and ATOM:',i4,' differ in shell',i3,2f10.5)
551 format(' atom:',i4,' and ATOM:',i4,' differ in shell',i3,2i4,2f7.1,f10.5)
552 format(//,'SHELL STRUCTURE for all ATOMS:',/,'ATOM  | DISTANCE   #of NEIGHBORS   Z |')
553 format(i3,1x,8('|',f6.3,i3,f6.1))
554 format(/,'LISTING of INDEX of all EQUIVALENT ATOMS:')
555 format(' ATOM:',i4,' and ATOM:',i4,' are equivalent')
556 format(/,' ATOM KIND:',i4,'  OLD and NEW MULTIPLICITY:  ',2i4)
557 format(5x,'ATOM INDEX:',i4,'  OLD and NEW ATOM KIND:',2i4)
559 format(' atom:',i4,' and ATOM:',i4,' differ in number of shells',i4,'ne',i4)

1011 FORMAT(4X,I4,4X,F10.8,3X,F10.8,3X,F10.8,5x,a10,f5.1)         
1012 FORMAT(4X,I4,4X,F10.8,3X,F10.8,3X,F10.8,/,15X,I2,17X,I2)     
1020 FORMAT(///,3X,'ERROR IN EBND  : MULT(JATOM)=0 ...',/, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
1049 FORMAT(A10,' NPT=',I5,'  R0=',F10.8,' RMT=',F10.4,'   Z:',f5.1)      
1149 FORMAT(A2,i1,7x,' NPT=',I5,'  R0=',F10.8,' RMT=',F10.4,'   Z:',f5.1)  
1148 FORMAT(A2,i2,6x,' NPT=',I5,'  R0=',F10.8,' RMT=',F10.4,'   Z:',f5.1)  
1147 FORMAT(A2,i3,5x,' NPT=',I5,'  R0=',F10.8,' RMT=',F10.4,'   Z:',f5.1)  
1050 FORMAT(A10,5X,I5,5X,F10.8,5X,F10.4,5x,f5.1)                           
1051 FORMAT(20X,3F10.7)                                                
1510 FORMAT(A79)                                                       
1511 FORMAT(A4,23X,I3,/,13X,A4)                                        
1512 FORMAT(20x,3f10.7)


  
END SUBROUTINE w2knn



SUBROUTINE REDUCE_ALLOC_R8_2D(N1OLD,N2OLD,N1NEW,N2NEW)
  !..REDUCES UP TO 4 DIMENSIONAL REAL*8 ARRAYS FROM OLD TO NEW ALLOCATIONS
  !
  use variable_fields,  a => pos 
  real*8,  allocatable ::  HELP(:,:)
  ALLOCATE ( HELP(N1NEW,n2new) )
  DO I=1,N1NEW
     DO j=1,N2NEW
        HELP(I,j)=A(I,j)
     ENDDO
  ENDDO
  DEALLOCATE ( A )
  ALLOCATE ( A(N1NEW,n2new) )
  DO I=1,N1NEW
     DO j=1,N2NEW
        A(I,j)=HELP(I,j)
     ENDDO
  ENDDO
  DEALLOCATE ( HELP )      
END SUBROUTINE REDUCE_ALLOC_R8_2D


SUBROUTINE DIRLAT (BR1,BR2,ortho,lattic,nat,alpha,beta,gamma,A)
  !                                                                       
  !     LATGEN GENERATES THE BRAVAIS MATRIX BR2(3,3), WHICH TRANSFORMS    
  !     A RECIPROCAL LATTICE VECTOR OF A SPECIAL COORDINATE SYSTEM (IN    
  !     UNITS OF 2 PI/A(I)) INTO CARTESIAN SYSTEM                         
  !     Convention:  R_i = (i,*)
  !                                                                       
  !use struk, only: br2,lattic
  IMPLICIT REAL*8 (A-H,O-Z)
  REAL*8, intent(out)    :: BR1(3,3), BR2(3,3)
  LOGICAL, intent(out)   :: ORTHO
  CHARACTER*4, intent(in):: lattic
  REAL*8, intent(in)     :: alpha,beta,gamma,A(3)
  REAL*8 :: PIA(3)

  WRITE(6,*) 'A=', A
  
  pi=4.d0*atan(1.d0)
  alpha0=alpha
  beta0=beta
  gamma1=gamma*pi/180.d0
  beta1=beta*pi/180.d0
  alpha1=alpha*pi/180.d0
  cosg1=(cos(gamma1)-cos(alpha1)*cos(beta1))/sin(alpha1)/sin(beta1)
  gamma0=acos(cosg1)
  SQRT3=SQRT(3.D0)

  PIA(1)=2.D0*PI/A(1)
  PIA(2)=2.D0*PI/A(2)
  PIA(3)=2.D0*PI/A(3)

  WRITE(6,*) 'PIA=', PIA, 'LATTIC=', LATTIC(1:1)

  IF(LATTIC(1:1).EQ.'H') GOTO 10                                    
  IF(LATTIC(1:1).EQ.'S') GOTO 20                                    
  IF(LATTIC(1:1).EQ.'P') GOTO 20                                    
  IF(LATTIC(1:1).EQ.'B') GOTO 30                                    
  IF(LATTIC(1:1).EQ.'F') GOTO 40                                    
  IF(LATTIC(1:1).EQ.'C') GOTO 50                                    
  IF(LATTIC(1:1).EQ.'R') GOTO 80                                    
  STOP 'LATTIC WRONG'                                               
  !                                                                       
  !.....HEXAGONAL CASE                                                    
10 CONTINUE
  BR2(1,1)=SQRT(3.d0)/2.d0                                              
  BR2(1,2)=-.5d0                                                      
  BR2(1,3)=0.0d0                                                      
  BR2(2,1)=0.0d0                                                     
  BR2(2,2)=1.0d0                                                      
  BR2(2,3)=0.0d0                                                      
  BR2(3,1)=0.0d0                                                      
  BR2(3,2)=0.0d0                                                      
  BR2(3,3)=1.d0                                                       
  ORTHO=.FALSE.
  !
  BR1(1,1)=2.D0/SQRT3*PIA(1)                                        
  BR1(1,2)=1.D0/SQRT3*PIA(1)                                        
  BR1(1,3)=0.0D0                                                    
  BR1(2,1)=0.0D0                                                    
  BR1(2,2)=PIA(2)                                                   
  BR1(2,3)=0.0D0                                                    
  BR1(3,1)=0.0D0                                                    
  BR1(3,2)=0.0D0                                                    
  BR1(3,3)=PIA(3)                                                   

  WRITE(6,*) 'BR1=', BR1

  GOTO 100                                                          
  !                                                                       
  !.....PRIMITIVE LATTICE CASE 
20 continue
  !
  BR2(1,1)=1.0d0*sin(gamma0)*sin(beta1)                
  BR2(1,2)=1.0d0*cos(gamma0)*sin(beta1)                 
  BR2(1,3)=1.0d0*cos(beta1)                                   
  BR2(2,1)=0.0d0                                                      
  BR2(2,2)=1.0d0*sin(alpha1)
  BR2(2,3)=1.0d0*cos(alpha1)
  BR2(3,1)=0.0d0                                                      
  BR2(3,2)=0.0d0                                                      
  BR2(3,3)=1.0d0                                                      
  ORTHO=.TRUE. 
  if(gamma.ne.90.d0) ortho=.false.                                      
  if(beta.ne.90.d0) ortho=.false.                                      
  if(alpha.ne.90.d0) ortho=.false.                                      
  !        write(*,*) alpha0,beta0,gamma0,ortho,br2
  !
  SINBC=SIN(alpha1)
  COSAB=COS(gamma1)
  COSAC=COS(beta1)
  COSBC=COS(alpha1)
  WURZEL=SQRT(SINBC**2-COSAC**2-COSAB**2+2*COSBC*COSAC*COSAB)
  BR1(1,1)= SINBC/WURZEL*PIA(1)
  BR1(1,2)= (-COSAB+COSBC*COSAC)/(SINBC*WURZEL)*PIA(2)
  BR1(1,3)= (COSBC*COSAB-COSAC)/(SINBC*WURZEL)*PIA(3)
  BR1(2,1)= 0.0
  BR1(2,2)= PIA(2)/SINBC
  BR1(2,3)= -PIA(3)*COSBC/SINBC
  BR1(3,1)= 0.0
  BR1(3,2)= 0.0
  BR1(3,3)= PIA(3)
  GOTO 100     
  !                                                                       
  !.....BC CASE (DIRECT LATTICE)                                          
30 CONTINUE                                                          
  BR2(1,1)=-0.5d0                                                     
  BR2(1,2)=0.5d0                                                      
  BR2(1,3)=0.5d0                                                      
  BR2(2,1)=0.5d0                                                     
  BR2(2,2)=-0.5d0                                                     
  BR2(2,3)=0.5d0                                                      
  BR2(3,1)=0.5d0                                                      
  BR2(3,2)=0.5d0                                                      
  BR2(3,3)=-0.5d0                                                     
  ORTHO=.TRUE.
  !
  BR1(1,1)=PIA(1)                                                   
  BR1(1,2)=0.0D0                                                    
  BR1(1,3)=0.0D0                                                    
  BR1(2,1)=0.0D0                                                    
  BR1(2,2)=PIA(2)                                                   
  BR1(2,3)=0.0D0                                                    
  BR1(3,1)=0.0D0                                                    
  BR1(3,2)=0.0D0                                                    
  BR1(3,3)=PIA(3)                                                   
  GOTO 100                                                          
  !                                                                       
  !.....FC CASE (DIRECT LATTICE)                                          
40 CONTINUE                                                          
  BR2(1,1)=0.0d0                                                      
  BR2(1,2)=0.5d0                                                      
  BR2(1,3)=0.5d0                                                      
  BR2(2,1)=0.5d0                                                      
  BR2(2,2)=0.0d0                                                      
  BR2(2,3)=0.5d0                                                      
  BR2(3,1)=0.5d0                                                      
  BR2(3,2)=0.5d0                                                      
  BR2(3,3)=0.0d0                                                      
  ORTHO=.TRUE.
  !
  BR1(1,1)=PIA(1)                                                   
  BR1(1,2)=0.0D0                                                    
  BR1(1,3)=0.0D0                                                    
  BR1(2,1)=0.0D0                                                    
  BR1(2,2)=PIA(2)                                                   
  BR1(2,3)=0.0D0                                                    
  BR1(3,2)=0.0D0                                                    
  BR1(3,1)=0.0D0                                                    
  BR1(3,3)=PIA(3)                                                   
  GOTO 100                                                          
  !                                                                       
  !.....CXY  CASE (DIRECT LATTICE)                                          
50 CONTINUE                                                          
  IF(LATTIC(2:3).EQ.'XZ') GOTO 60                                    
  IF(LATTIC(2:3).EQ.'YZ') GOTO 70                                    
  BR2(1,1)=0.5d0                                                      
  BR2(1,2)=-0.5d0                                                     
  BR2(1,3)=0.0d0                                                      
  BR2(2,1)=0.5d0                                                      
  BR2(2,2)=0.5d0                                                      
  BR2(2,3)=0.0d0                                                      
  BR2(3,1)=0.0d0                                                      
  BR2(3,2)=0.0d0                                                      
  BR2(3,3)=1.0d0                                                      
  ORTHO=.TRUE.
  !
  BR1(1,1)=PIA(1)                                                   
  BR1(1,2)=0.0D0                                                    
  BR1(1,3)=0.0D0                                                    
  BR1(2,1)=0.0D0                                                    
  BR1(2,2)=PIA(2)                                                   
  BR1(2,3)=0.0D0                                                    
  BR1(3,1)=0.0D0                                                    
  BR1(3,2)=0.0D0                                                    
  BR1(3,3)=PIA(3)                                                   
  GOTO 100                                                          
  !                                                                       
  !.....CXZ  CASE (DIRECT LATTICE)                                          
60 CONTINUE 
  !.....CXZ ORTHOROMBIC CASE
  if(gamma.eq.90.d0) then
     BR2(1,1)=0.5d0                                                      
     BR2(1,2)=0.0d0                                                     
     BR2(1,3)=-0.5d0                                                      
     BR2(2,1)=0.0d0                                                      
     BR2(2,2)=1.0d0                                                      
     BR2(2,3)=0.0d0                                                      
     BR2(3,1)=0.5d0                                                     
     BR2(3,2)=0.0d0                                                      
     BR2(3,3)=0.5d0                                                      
     ORTHO=.TRUE.
     !
     BR1(1,1)=PIA(1)                                                   
     BR1(1,2)=0.0D0                                                    
     BR1(1,3)=0.0D0                                                    
     BR1(2,1)=0.0D0                                                    
     BR1(2,2)=PIA(2)                                                   
     BR1(2,3)=0.0D0                                                    
     BR1(3,1)=0.0D0                                                    
     BR1(3,2)=0.0D0                                                    
     BR1(3,3)=PIA(3)                                                   
     GOTO 100                                                          
  ELSE
     !.....CXZ MONOCLINIC CASE
     write(*,*) 'gamma not equal 90'
     SINAB=SIN(gamma1)
     COSAB=COS(gamma1)
     !
     BR2(1,1)=0.5d0*sinab                                                
     BR2(1,2)=0.5d0*cosab                                               
     BR2(1,3)=-0.5d0                                                      
     BR2(2,1)=0.0d0                                                      
     BR2(2,2)=1.0d0                                                      
     BR2(2,3)=0.0d0                                                      
     BR2(3,1)=0.5d0*sinab                                               
     BR2(3,2)=0.5d0*cosab                                                
     BR2(3,3)=0.5d0                                                      
     ORTHO=.FALSE.
     !
     BR1(1,1)= PIA(1)/SINAB 
     BR1(1,2)= -PIA(2)*COSAB/SINAB
     BR1(1,3)= 0.0                                                   
     BR1(2,1)= 0.0                                                      
     BR1(2,2)= PIA(2)                                                     
     BR1(2,3)= 0.0                                                     
     BR1(3,1)= 0.0                                                     
     BR1(3,2)= 0.0                                                     
     BR1(3,3)= PIA(3)                                                     
     GOTO 100
  ENDIF

  !                                                                       
  !.....CYZ  CASE (DIRECT LATTICE)                                          
70 CONTINUE                                                          
  BR2(1,1)=1.0d0                                                      
  BR2(1,2)=0.0d0                                                     
  BR2(1,3)=0.0d0                                                      
  BR2(2,1)=0.0d0                                                      
  BR2(2,2)=0.5d0                                                      
  BR2(2,3)=0.5d0                                                      
  BR2(3,1)=0.0d0                                                      
  BR2(3,2)=-0.5d0                                                     
  BR2(3,3)=0.5d0                                                      
  ORTHO=.TRUE.
  !!
  BR1(1,1)=PIA(1)                                                   
  BR1(1,2)=0.0D0                                                    
  BR1(1,3)=0.0D0                                                    
  BR1(2,1)=0.0D0                                                    
  BR1(2,2)=PIA(2)                                                   
  BR1(2,3)=0.0D0                                                    
  BR1(3,1)=0.0D0                                                    
  BR1(3,2)=0.0D0                                                    
  BR1(3,3)=PIA(3)                                                   
  GOTO 100
  
  !.....RHOMBOHEDRAL CASE
80 CONTINUE
  BR2(1,1)=1/2.d0/sqrt(3.d0)
  BR2(1,2)=-1/2.d0                                                     
  BR2(1,3)=1/3.d0                                                      
  BR2(2,1)=1/2.d0/SQRT(3.d0)                                          
  BR2(2,2)=1*0.5d0                                                
  BR2(2,3)=1/3.d0                                                      
  BR2(3,1)=-1/SQRT(3.d0)                                         
  BR2(3,2)=0.d0                                                
  BR2(3,3)=1/3.d0                                                      
  ORTHO=.FALSE.
  !
  BR1(1,1)=1.D0/SQRT(3.D0)*PIA(1)                                          
  BR1(1,2)=1.D0/SQRT(3.D0)*PIA(1)                                          
  BR1(1,3)=-2.d0/sqrt(3.d0)*PIA(1)                                         
  BR1(2,1)=-1.0d0*PIA(2)                                                  
  BR1(2,2)=1.0d0*PIA(2)                                                    
  BR1(2,3)=0.0d0*PIA(2)                                                    
  BR1(3,1)=1.0d0*PIA(3)                                                    
  BR1(3,2)=1.0d0*PIA(3)                                                    
  BR1(3,3)=1.0d0*PIA(3)                                                    
  GOTO 100                                                         
  !                                                                       
100 CONTINUE                                                          
  write(66,*) 'Bravais Matrix:'
  write(66,999) br2
999 format(3f15.5)
  !                                                                       
  RETURN                                                            
END SUBROUTINE DIRLAT


SUBROUTINE ORD2(A,NR,IMAX)                                        
  !     ORDERS ELEMENTS IN ARRAY A INCREASING IN SIZE                     
  !       REORDERS CORRESPONDING INDICES (NR)                             
  IMPLICIT REAL*8 (A-H,O-Z)
  LOGICAL CONT                                                      
  DIMENSION A(*),NR(*)                                              
  DO I=1,IMAX  ! 50
     NR(I)=I
  ENDDO
  if(imax.eq.1) return
  
  CONT=.TRUE.
  do while(CONT)
     CONT=.FALSE.
     DO I=2,IMAX
        IF(A(I).LT.A(I-1)) THEN
           !       INTERCHANGE I AND (I-1) ELEMENT IN ARRAYS A AND NR              
           HLP=A(I)                                                        
           A(I)=A(I-1)                                                     
           A(I-1)=HLP                                                      
           NHLP=NR(I)                                                      
           NR(I)=NR(I-1)                                                   
           NR(I-1)=NHLP                                                    
           CONT=.TRUE.                                                     
        ENDIF
     ENDDO
  enddo
  RETURN                                                
END SUBROUTINE

!subroutine bvan(i1,i2,scale,dist) 
!  !     Code for bond-valence parameters
!  !     Form is exp( (dist-R)/scale ), R=distancd in Angers 
!  !     from http://www.ccp14.ac.uk/ccp/web-mirrors/i_d_brown/
!  !     L. D. Marks & J. A. Enterkin, Sept 2009
!  ! 
!  implicit real*8 (a-h,o-z) 
!  dimension BVA(2,98,98) 
!  save BVA 
!      data ( ( BVA(i1, 1,i2),i1=1,2),i2=1,98) /  &  ! H 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   1.110,   0.370,   1.140, & 
!         0.370,   1.100,  -9.000,  -9.000,   0.280,   0.907,  -9.000,  -9.000,  -9.000,  -9.000, & 
!         0.370,   1.680,   0.370,   1.530,   0.370,   1.450,   0.370,   1.470,   0.370,   1.410, & 
!         0.370,   1.380,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.100,   0.370,   1.830, & 
!         0.370,   1.680,   0.370,   1.610,   0.370,   1.580,   0.370,   1.520,   0.370,   1.550, & 
!         0.370,   1.530,   0.370,   1.440,   0.370,   1.400,   0.370,   1.210,   0.370,   1.420, & 
!         0.370,   1.510,   0.370,   1.550,  -9.000,  -9.000,   0.370,   1.540,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.260,   0.370,   2.010,   0.370,   1.860,   0.370,   1.790, & 
!         0.370,   1.750,   0.370,   1.730,  -9.000,  -9.000,   0.370,   1.610,   0.370,   1.550, & 
!         0.370,   1.470,   0.370,   1.500,   0.370,   1.660,   0.370,   1.720,   0.370,   1.850, & 
!         0.370,   2.770,   0.370,   1.830,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.440, & 
!         0.370,   2.220,   0.370,   2.060,   0.370,   2.040,   0.370,   2.020,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.960,   0.370,   1.950,   0.370,   1.930,   0.370,   1.910, & 
!         0.370,   1.890,   0.370,   1.880,   0.370,   1.860,   0.370,   1.850,   0.370,   1.820, & 
!         0.370,   1.820,   0.370,   1.780,   0.370,   1.760,   0.370,   1.760,   0.370,   1.750, & 
!        -9.000,  -9.000,   0.370,   1.760,   0.370,   1.400,   0.370,   1.370,   0.370,   1.710, & 
!         0.370,   2.050,   0.370,   1.970,   0.370,   1.970,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.070, & 
!        -9.000,  -9.000,   0.370,   1.970,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1, 2,i2),i1=1,2),i2=1,98) /  &  ! He
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1, 3,i2),i1=1,2),i2=1,98) /  &  ! Li
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.610,   0.370,   1.466,   0.370,   1.360,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!         0.370,   1.940,   0.370,   1.910,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.090,   0.370,   2.020, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.300,   0.370,   2.220,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1, 4,i2),i1=1,2),i2=1,98) /  &  ! Be
!         0.370,   1.110,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.500,   0.370,   1.381,   0.370,   1.281,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   1.950, & 
!         0.370,   1.830,   0.370,   1.760,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.000,   0.370,   1.970,   0.370,   1.900, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.210,   0.370,   2.100,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1, 5,i2),i1=1,2),i2=1,98) /  &  ! B 
!         0.370,   1.140,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   1.402, & 
!        -9.000,  -9.000,   0.370,   1.470,   0.370,   1.371,   0.370,   1.310,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   1.880, & 
!         0.370,   1.820,   0.370,   1.740,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   1.970,   0.370,   1.950,   0.370,   1.880, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.200,   0.370,   2.100,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1, 6,i2),i1=1,2),i2=1,98) /  &  ! C 
!         0.370,   1.100,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!         0.370,   1.540,   0.370,   1.470,   0.370,   1.390,   0.370,   1.320,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   1.883,   0.370,   1.890, & 
!         0.370,   1.820,   0.370,   1.760,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!         0.370,   1.689,   0.370,   1.634,  -9.000,  -9.000,   0.370,   1.720,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   1.930,   0.370,   1.970,   0.370,   1.900, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!         0.370,   1.730,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.210,   0.370,   2.120,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   1.760,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1, 7,i2),i1=1,2),i2=1,98) /  &  ! N 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   1.610,   0.370,   1.500,   0.370,   1.470, & 
!         0.370,   1.470,   0.350,   1.440,   0.370,   1.432,   0.370,   1.360,  -9.000,  -9.000, & 
!         0.370,   1.930,   0.370,   1.850,   0.370,   1.790,   0.370,   1.770,   0.370,   1.730, & 
!         0.370,   1.740,   0.370,   1.800,  -9.000,  -9.000,   0.370,   2.260,   0.370,   2.140, & 
!         0.370,   1.980,   0.370,   1.930,   0.370,   1.860,   0.370,   1.850,   0.370,   1.870, & 
!         0.370,   1.860,   0.370,   1.840,   0.370,   1.750,   0.370,   1.610,   0.370,   1.720, & 
!         0.370,   1.840,   0.370,   1.880,  -9.000,  -9.000,   0.350,   1.900,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.370,   0.370,   2.230,   0.370,   2.170,   0.370,   2.110, & 
!         0.370,   2.060,   0.370,   2.040,  -9.000,  -9.000,   0.370,   1.880,   0.370,   1.880, & 
!         0.370,   1.810,   0.370,   1.850,   0.370,   1.960,   0.370,   2.030,   0.370,   2.140, & 
!         0.370,   2.120,   0.370,   2.120,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.530, & 
!         0.370,   2.470,   0.370,   2.340,   0.370,   2.340,   0.370,   2.300,   0.370,   2.300, & 
!        -9.000,  -9.000,   0.370,   2.240,   0.370,   2.240,   0.370,   2.220,   0.370,   2.200, & 
!         0.370,   2.180,   0.370,   2.180,   0.370,   2.160,   0.370,   2.140,   0.370,   2.120, & 
!         0.370,   2.110,   0.370,   2.090,   0.370,   2.010,   0.370,   2.060,   0.370,   2.060, & 
!        -9.000,  -9.000,   0.370,   2.060,   0.370,   1.770,   0.370,   1.720,   0.370,   2.020, & 
!         0.370,   2.290,   0.370,   2.220,   0.370,   2.240,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.340, & 
!        -9.000,  -9.000,   0.370,   2.240,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1, 8,i2),i1=1,2),i2=1,98) /  &  ! O 
!         0.280,   0.907,  -9.000,  -9.000,   0.370,   1.466,   0.370,   1.381,   0.370,   1.371, & 
!         0.370,   1.390,   0.370,   1.432,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!         0.370,   1.803,   0.370,   1.693,   0.370,   1.620,   0.370,   1.624,   0.370,   1.604, & 
!         0.370,   1.644,   0.370,   1.632,  -9.000,  -9.000,   0.370,   2.132,   0.370,   1.967, & 
!         0.370,   1.849,   0.370,   1.790,   0.320,   1.788,   0.370,   1.724,   0.370,   1.754, & 
!         0.370,   1.745,   0.370,   1.661,   0.370,   1.654,   0.370,   1.679,   0.370,   1.704, & 
!         0.370,   1.730,   0.370,   1.748,   0.370,   1.767,   0.370,   1.788,   0.370,   1.810, & 
!        -9.000,  -9.000,   0.370,   2.263,   0.370,   2.118,   0.370,   2.014,   0.370,   1.937, & 
!         0.370,   1.911,   0.370,   1.907,   0.370,   1.900,   0.370,   1.834,   0.370,   1.793, & 
!         0.370,   1.792,   0.370,   1.805,   0.370,   1.904,   0.370,   1.902,   0.370,   1.905, & 
!         0.370,   1.942,   0.370,   1.917,   0.370,   1.930,   0.370,   2.000,   0.370,   2.417, & 
!         0.370,   2.285,   0.370,   2.172,   0.370,   2.151,   0.370,   2.138,   0.370,   2.117, & 
!        -9.000,  -9.000,   0.370,   2.088,   0.370,   2.074,   0.370,   2.065,   0.370,   2.049, & 
!         0.370,   2.001,   0.370,   2.025,   0.370,   1.988,   0.370,   2.000,   0.370,   1.985, & 
!         0.370,   1.971,   0.370,   1.923,   0.370,   1.920,   0.370,   1.921,   0.370,   1.970, & 
!         0.370,   1.811,   0.370,   1.916,   0.370,   1.879,   0.370,   1.833,   0.370,   1.930, & 
!         0.370,   2.003,   0.370,   2.042,   0.370,   2.060,   0.370,   2.190,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.240,   0.370,   2.167, & 
!         0.350,   2.090,   0.370,   2.112,   0.370,   2.180,   0.370,   2.110,   0.370,   2.110, & 
!         0.370,   2.230,   0.370,   2.080,   0.370,   2.070/
!      data ( ( BVA(i1, 9,i2),i1=1,2),i2=1,98) /  &  ! F 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   1.360,   0.370,   1.281,   0.370,   1.310, & 
!         0.370,   1.320,   0.370,   1.360,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!         0.370,   1.677,   0.370,   1.578,   0.370,   1.545,   0.370,   1.580,   0.370,   1.540, & 
!         0.370,   1.600,   0.370,   1.550,  -9.000,  -9.000,   0.370,   1.992,   0.370,   1.842, & 
!         0.370,   1.760,   0.370,   1.760,   0.370,   1.700,   0.370,   1.640,   0.370,   1.710, & 
!         0.370,   1.650,   0.370,   1.620,   0.370,   1.596,   0.370,   1.594,   0.370,   1.620, & 
!         0.370,   1.620,   0.370,   1.660,   0.370,   1.620,   0.370,   1.690,   0.370,   1.720, & 
!         0.370,   1.880,   0.370,   2.160,   0.370,   2.019,   0.370,   1.904,   0.370,   1.854, & 
!         0.370,   1.870,   0.370,   1.810,   0.400,   1.880,   0.370,   1.740,   0.370,   1.710, & 
!         0.370,   1.740,   0.370,   1.800,   0.370,   1.811,   0.370,   1.792,   0.370,   1.843, & 
!         0.370,   1.797,   0.370,   1.820,   0.370,   1.830,   0.370,   1.890,   0.370,   2.330, & 
!         0.370,   2.188,  -9.000,  -9.000,   0.370,   2.036,   0.370,   2.022,   0.370,   2.008, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   1.961,   0.370,   1.950,   0.370,   1.936, & 
!         0.370,   1.922,   0.370,   1.908,   0.370,   1.904,   0.370,   1.842,   0.370,   1.875, & 
!         0.370,   1.876,   0.370,   1.850,   0.370,   1.880,   0.370,   1.830,   0.370,   1.860, & 
!         0.370,   1.720,   0.370,   1.820,   0.370,   1.759,   0.370,   1.810,   0.370,   1.900, & 
!         0.370,   1.880,   0.370,   1.940,   0.370,   1.970,   0.370,   2.380,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.130,   0.370,   2.068, & 
!         0.370,   2.040,   0.370,   2.034,   0.370,   2.020,   0.370,   2.000,   0.370,   2.000, & 
!         0.370,   2.120,   0.370,   1.960,   0.370,   1.950/
!      data ( ( BVA(i1,10,i2),i1=1,2),i2=1,98) /  &  ! Ne
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,11,i2),i1=1,2),i2=1,98) /  &  ! Na
!         0.370,   1.680,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.930,   0.370,   1.803,   0.370,   1.677,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.360, & 
!         0.370,   2.280,   0.370,   2.150,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.530,   0.370,   2.410,   0.370,   2.330, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.640,   0.370,   2.560,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,12,i2),i1=1,2),i2=1,98) /  &  ! Mg
!         0.370,   1.530,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.850,   0.370,   1.693,   0.370,   1.578,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.290, & 
!         0.370,   2.180,   0.370,   2.080,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.380,   0.370,   2.320,   0.370,   2.280, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.530,   0.370,   2.460,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,13,i2),i1=1,2),i2=1,98) /  &  ! Al
!         0.370,   1.450,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.790,   0.370,   1.620,   0.370,   1.545,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.240, & 
!         0.370,   2.130,   0.370,   2.032,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.300,   0.370,   2.270,   0.370,   2.200, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.480,   0.370,   2.410,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,14,i2),i1=1,2),i2=1,98) /  &  ! Si
!         0.370,   1.470,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!         0.370,   1.883,   0.370,   1.770,   0.370,   1.624,   0.370,   1.580,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.230, & 
!         0.370,   2.126,   0.370,   2.030,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.310,   0.370,   2.260,   0.370,   2.200, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.490,   0.370,   2.410,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,15,i2),i1=1,2),i2=1,98) /  &  ! P 
!         0.370,   1.410,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   1.950,   0.370,   1.880, & 
!         0.370,   1.890,   0.370,   1.730,   0.370,   1.604,   0.370,   1.540,  -9.000,  -9.000, & 
!         0.370,   2.360,   0.370,   2.290,   0.370,   2.240,   0.370,   2.230,   0.370,   2.190, & 
!         0.370,   2.150,   0.370,   2.020,  -9.000,  -9.000,   0.370,   2.640,   0.370,   2.550, & 
!         0.370,   2.400,   0.370,   2.360,   0.370,   2.310,   0.370,   2.270,   0.370,   2.240, & 
!         0.370,   2.270,   0.370,   2.210,   0.370,   2.170,   0.370,   1.970,   0.370,   2.150, & 
!         0.370,   2.260,   0.370,   2.320,   0.370,   2.250,   0.370,   2.340,   0.370,   2.150, & 
!        -9.000,  -9.000,   0.370,   2.760,   0.370,   2.670,   0.370,   2.570,   0.370,   2.520, & 
!         0.370,   2.460,   0.370,   2.440,  -9.000,  -9.000,   0.370,   2.290,   0.370,   2.290, & 
!         0.370,   2.220,   0.370,   2.220,   0.370,   2.340,   0.370,   2.430,   0.370,   2.450, & 
!         0.370,   2.520,   0.370,   2.520,   0.370,   2.400,  -9.000,  -9.000,   0.370,   2.930, & 
!         0.370,   2.880,   0.370,   2.730,   0.370,   2.700,   0.370,   2.680,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.630,   0.370,   2.620,   0.370,   2.610,   0.370,   2.590, & 
!         0.370,   2.570,   0.370,   2.560,   0.370,   2.550,   0.370,   2.530,   0.370,   2.530, & 
!         0.370,   2.510,   0.370,   2.480,   0.370,   2.470,   0.370,   2.460,   0.370,   2.460, & 
!        -9.000,  -9.000,   0.370,   2.460,   0.370,   2.190,   0.370,   2.140,   0.370,   2.420, & 
!         0.370,   2.710,   0.370,   2.640,   0.370,   2.630,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.730, & 
!        -9.000,  -9.000,   0.370,   2.640,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,16,i2),i1=1,2),i2=1,98) /  &  ! S 
!         0.370,   1.380,  -9.000,  -9.000,   0.370,   1.940,   0.370,   1.830,   0.370,   1.820, & 
!         0.370,   1.820,   0.370,   1.740,   0.370,   1.644,   0.370,   1.600,  -9.000,  -9.000, & 
!         0.370,   2.280,   0.370,   2.180,   0.370,   2.130,   0.370,   2.126,   0.370,   2.150, & 
!         0.370,   2.070,   0.370,   2.020,  -9.000,  -9.000,   0.370,   2.590,   0.370,   2.450, & 
!         0.370,   2.321,   0.370,   2.240,   0.370,   2.230,   0.370,   2.180,   0.370,   2.200, & 
!         0.370,   2.160,   0.370,   2.060,   0.370,   2.040,   0.370,   1.860,   0.370,   2.090, & 
!         0.370,   2.170,   0.370,   2.230,   0.370,   2.272,   0.370,   2.250,   0.370,   2.170, & 
!        -9.000,  -9.000,   0.370,   2.700,   0.370,   2.590,   0.370,   2.480,   0.370,   2.410, & 
!         0.370,   2.370,   0.370,   2.350,  -9.000,  -9.000,   0.370,   2.160,   0.370,   2.150, & 
!         0.370,   2.100,   0.370,   2.119,   0.370,   2.304,   0.370,   2.360,   0.370,   2.450, & 
!         0.370,   2.450,   0.370,   2.450,   0.370,   2.360,  -9.000,  -9.000,   0.370,   2.890, & 
!         0.370,   2.769,   0.370,   2.643,   0.370,   2.620,   0.370,   2.600,   0.370,   2.590, & 
!        -9.000,  -9.000,   0.370,   2.550,   0.370,   2.530,   0.370,   2.530,   0.370,   2.510, & 
!         0.370,   2.470,   0.370,   2.480,   0.370,   2.460,   0.370,   2.450,   0.370,   2.430, & 
!         0.370,   2.430,   0.370,   2.390,   0.370,   2.390,   0.370,   2.390,   0.370,   2.370, & 
!         0.370,   2.210,   0.370,   2.380,   0.370,   2.080,   0.370,   2.030,   0.370,   2.320, & 
!         0.370,   2.630,   0.370,   2.550,   0.370,   2.550,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.640, & 
!        -9.000,  -9.000,   0.370,   2.560,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,17,i2),i1=1,2),i2=1,98) /  &  ! Cl
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   1.910,   0.370,   1.760,   0.370,   1.740, & 
!         0.370,   1.760,   0.370,   1.800,   0.370,   1.632,   0.370,   1.550,  -9.000,  -9.000, & 
!         0.370,   2.150,   0.370,   2.080,   0.370,   2.032,   0.370,   2.030,   0.370,   2.020, & 
!         0.370,   2.020,   0.370,   2.000,  -9.000,  -9.000,   0.370,   2.519,   0.370,   2.370, & 
!         0.370,   2.230,   0.370,   2.184,   0.370,   2.160,   0.370,   2.080,   0.370,   2.130, & 
!         0.370,   2.060,   0.370,   2.050,   0.370,   2.020,   0.370,   2.000,   0.370,   2.010, & 
!         0.370,   2.070,   0.370,   2.140,   0.370,   2.140,   0.370,   2.160,   0.370,   2.190, & 
!        -9.000,  -9.000,   0.370,   2.652,   0.370,   2.510,   0.370,   2.400,   0.370,   2.330, & 
!         0.370,   2.270,   0.370,   2.280,   0.370,   2.210,   0.370,   2.210,   0.370,   2.170, & 
!         0.370,   2.050,   0.370,   2.090,   0.370,   2.230,   0.370,   2.280,   0.370,   2.276, & 
!         0.370,   2.300,   0.370,   2.300,   0.370,   2.310,  -9.000,  -9.000,   0.370,   2.791, & 
!         0.370,   2.690,   0.370,   2.545,   0.370,   2.410,   0.370,   2.500,   0.370,   2.492, & 
!        -9.000,  -9.000,   0.370,   1.977,   0.370,   2.530,   0.370,   2.445,   0.370,   2.427, & 
!         0.370,   2.410,   0.370,   2.401,   0.370,   2.390,   0.370,   2.380,   0.370,   2.371, & 
!         0.370,   2.361,   0.370,   2.300,   0.370,   2.300,   0.370,   2.270,   0.370,   2.230, & 
!         0.370,   2.190,   0.370,   2.300,   0.370,   2.170,   0.370,   2.170,   0.370,   2.250, & 
!         0.370,   2.320,   0.370,   2.430,   0.370,   2.440,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.630,   0.370,   2.550, & 
!        -9.000,  -9.000,   0.370,   2.460,   0.400,   2.460,   0.370,   2.480,   0.370,   2.480, & 
!         0.370,   2.620,   0.370,   2.460,   0.370,   2.450/
!      data ( ( BVA(i1,18,i2),i1=1,2),i2=1,98) /  &  ! Ar
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,19,i2),i1=1,2),i2=1,98) /  &  ! K 
!         0.370,   2.100,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.260,   0.370,   2.132,   0.370,   1.992,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.640, & 
!         0.370,   2.590,   0.370,   2.519,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.830,   0.370,   2.720,   0.370,   2.660, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.930,   0.370,   2.880,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,20,i2),i1=1,2),i2=1,98) /  &  ! Ca
!         0.370,   1.830,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.140,   0.370,   1.967,   0.370,   1.842,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.550, & 
!         0.370,   2.450,   0.370,   2.370,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.620,   0.370,   2.560,   0.370,   2.490, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.760,   0.370,   2.720,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,21,i2),i1=1,2),i2=1,98) /  &  ! Sc
!         0.370,   1.680,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.980,   0.370,   1.849,   0.370,   1.760,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.400, & 
!         0.370,   2.321,   0.370,   2.230,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.480,   0.370,   2.440,   0.370,   2.380, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.640,   0.370,   2.590,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,22,i2),i1=1,2),i2=1,98) /  &  ! Ti
!         0.370,   1.610,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.930,   0.370,   1.790,   0.370,   1.760,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.360, & 
!         0.370,   2.240,   0.370,   2.184,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.420,   0.370,   2.380,   0.370,   2.320, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.600,   0.370,   2.540,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,23,i2),i1=1,2),i2=1,98) /  &  ! V 
!         0.370,   1.580,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.860,   0.320,   1.788,   0.370,   1.700,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.310, & 
!         0.370,   2.230,   0.370,   2.160,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.390,   0.370,   2.330,   0.370,   2.300, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.570,   0.370,   2.510,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,24,i2),i1=1,2),i2=1,98) /  &  ! Cr
!         0.370,   1.520,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.850,   0.370,   1.724,   0.370,   1.640,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.270, & 
!         0.370,   2.180,   0.370,   2.080,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.340,   0.370,   2.290,   0.370,   2.260, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.520,   0.370,   2.450,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,25,i2),i1=1,2),i2=1,98) /  &  ! Mn
!         0.370,   1.550,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.870,   0.370,   1.754,   0.370,   1.710,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.240, & 
!         0.370,   2.200,   0.370,   2.130,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.360,   0.370,   2.320,   0.370,   2.260, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.350,   2.604,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.550,   0.370,   2.490,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,26,i2),i1=1,2),i2=1,98) /  &  ! Fe
!         0.370,   1.530,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!         0.370,   1.689,   0.370,   1.860,   0.370,   1.745,   0.370,   1.650,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.270, & 
!         0.370,   2.160,   0.370,   2.060,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.350,   0.370,   2.280,   0.370,   2.260, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.530,   0.370,   2.470,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,27,i2),i1=1,2),i2=1,98) /  &  ! Co
!         0.370,   1.440,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!         0.370,   1.634,   0.370,   1.840,   0.370,   1.661,   0.370,   1.620,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.210, & 
!         0.370,   2.060,   0.370,   2.050,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.280,   0.370,   2.240,   0.370,   2.180, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.350,   2.593,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.460,   0.350,   2.370,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,28,i2),i1=1,2),i2=1,98) /  &  ! Ni
!         0.370,   1.400,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.750,   0.370,   1.654,   0.370,   1.596,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.170, & 
!         0.370,   2.040,   0.370,   2.020,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.240,   0.370,   2.140,   0.370,   2.160, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.430,   0.370,   2.340,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,29,i2),i1=1,2),i2=1,98) /  &  ! Cu
!         0.370,   1.210,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!         0.370,   1.720,   0.370,   1.610,   0.370,   1.679,   0.370,   1.594,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   1.970, & 
!         0.370,   1.860,   0.370,   2.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.080,   0.370,   2.020,   0.370,   1.990, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.270,   0.370,   2.160,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,30,i2),i1=1,2),i2=1,98) /  &  ! Zn
!         0.370,   1.420,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.720,   0.370,   1.704,   0.370,   1.620,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.150, & 
!         0.370,   2.090,   0.370,   2.010,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.240,   0.370,   2.220,   0.370,   2.150, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.450,   0.370,   2.360,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,31,i2),i1=1,2),i2=1,98) /  &  ! Ga
!         0.370,   1.510,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.840,   0.370,   1.730,   0.370,   1.620,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.260, & 
!         0.370,   2.170,   0.370,   2.070,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.340,   0.370,   2.300,   0.370,   2.240, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.540,   0.370,   2.450,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,32,i2),i1=1,2),i2=1,98) /  &  ! Ge
!         0.370,   1.550,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.880,   0.370,   1.748,   0.370,   1.660,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.320, & 
!         0.370,   2.230,   0.370,   2.140,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.430,   0.370,   2.350,   0.370,   2.300, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.560,   0.370,   2.500,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,33,i2),i1=1,2),i2=1,98) /  &  ! As
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.000,   0.370,   1.970, & 
!         0.370,   1.930,  -9.000,  -9.000,   0.370,   1.767,   0.370,   1.620,  -9.000,  -9.000, & 
!         0.370,   2.530,   0.370,   2.380,   0.370,   2.300,   0.370,   2.310,   0.370,   2.250, & 
!         0.370,   2.272,   0.370,   2.140,  -9.000,  -9.000,   0.370,   2.830,   0.370,   2.620, & 
!         0.370,   2.480,   0.370,   2.420,   0.370,   2.390,   0.370,   2.340,   0.370,   2.360, & 
!         0.370,   2.350,   0.370,   2.280,   0.370,   2.240,   0.370,   2.080,   0.370,   2.240, & 
!         0.370,   2.340,   0.370,   2.430,  -9.000,  -9.000,   0.370,   2.420,   0.370,   2.350, & 
!        -9.000,  -9.000,   0.370,   2.870,   0.370,   2.760,   0.370,   2.640,   0.370,   2.570, & 
!         0.370,   2.540,   0.370,   2.520,  -9.000,  -9.000,   0.370,   2.360,   0.370,   2.370, & 
!         0.370,   2.300,   0.370,   2.300,   0.370,   2.430,   0.370,   2.510,   0.370,   2.620, & 
!         0.370,   2.600,   0.370,   2.600,   0.370,   2.580,  -9.000,  -9.000,   0.370,   3.040, & 
!         0.370,   2.960,   0.370,   2.800,   0.370,   2.780,   0.370,   2.750,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.700,   0.370,   2.700,   0.370,   2.680,   0.370,   2.660, & 
!         0.370,   2.640,   0.370,   2.640,   0.370,   2.630,   0.370,   2.620,   0.370,   2.590, & 
!         0.370,   2.590,   0.370,   2.560,   0.370,   2.550,   0.370,   2.540,   0.370,   2.540, & 
!        -9.000,  -9.000,   0.370,   2.540,   0.370,   2.260,   0.370,   2.220,   0.370,   2.500, & 
!         0.370,   2.790,   0.370,   2.720,   0.370,   2.720,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.800, & 
!        -9.000,  -9.000,   0.370,   2.720,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,34,i2),i1=1,2),i2=1,98) /  &  ! Se
!         0.370,   1.540,  -9.000,  -9.000,   0.370,   2.090,   0.370,   1.970,   0.370,   1.950, & 
!         0.370,   1.970,   0.350,   1.900,   0.370,   1.788,   0.370,   1.690,  -9.000,  -9.000, & 
!         0.370,   2.410,   0.370,   2.320,   0.370,   2.270,   0.370,   2.260,   0.370,   2.340, & 
!         0.370,   2.250,   0.370,   2.160,  -9.000,  -9.000,   0.370,   2.720,   0.370,   2.560, & 
!         0.370,   2.440,   0.370,   2.380,   0.370,   2.330,   0.370,   2.290,   0.370,   2.320, & 
!         0.370,   2.280,   0.370,   2.240,   0.370,   2.140,   0.370,   2.020,   0.370,   2.220, & 
!         0.370,   2.300,   0.370,   2.350,   0.370,   2.420,   0.370,   2.360,   0.370,   2.330, & 
!        -9.000,  -9.000,   0.370,   2.810,   0.370,   2.720,   0.370,   2.610,   0.370,   2.530, & 
!         0.370,   2.510,   0.370,   2.490,  -9.000,  -9.000,   0.370,   2.330,   0.370,   2.330, & 
!         0.370,   2.220,   0.370,   2.260,   0.370,   2.400,   0.370,   2.470,   0.370,   2.590, & 
!         0.370,   2.570,   0.370,   2.530,   0.370,   2.540,  -9.000,  -9.000,   0.370,   2.980, & 
!         0.370,   2.880,   0.370,   2.740,   0.370,   2.740,   0.370,   2.720,   0.370,   2.710, & 
!        -9.000,  -9.000,   0.370,   2.670,   0.370,   2.660,   0.370,   2.650,   0.370,   2.630, & 
!         0.370,   2.610,   0.370,   2.610,   0.370,   2.590,   0.370,   2.580,   0.370,   2.560, & 
!         0.370,   2.560,   0.370,   2.520,   0.370,   2.510,   0.370,   2.510,   0.370,   2.500, & 
!        -9.000,  -9.000,   0.370,   2.510,   0.370,   2.190,   0.370,   2.180,   0.370,   2.470, & 
!         0.370,   2.700,   0.370,   2.670,   0.370,   2.720,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.760, & 
!        -9.000,  -9.000,   0.370,   2.700,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,35,i2),i1=1,2),i2=1,98) /  &  ! Br
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.020,   0.370,   1.900,   0.370,   1.880, & 
!         0.370,   1.900,  -9.000,  -9.000,   0.370,   1.810,   0.370,   1.720,  -9.000,  -9.000, & 
!         0.370,   2.330,   0.370,   2.280,   0.370,   2.200,   0.370,   2.200,   0.370,   2.150, & 
!         0.370,   2.170,   0.370,   2.190,  -9.000,  -9.000,   0.370,   2.660,   0.370,   2.490, & 
!         0.370,   2.380,   0.370,   2.320,   0.370,   2.300,   0.370,   2.260,   0.370,   2.260, & 
!         0.370,   2.260,   0.370,   2.180,   0.370,   2.160,   0.370,   1.990,   0.370,   2.150, & 
!         0.370,   2.240,   0.370,   2.300,   0.370,   2.350,   0.370,   2.330,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.780,   0.370,   2.680,   0.370,   2.550,   0.370,   2.480, & 
!         0.370,   2.450,   0.370,   2.430,  -9.000,  -9.000,   0.370,   2.260,   0.370,   2.250, & 
!         0.370,   2.190,   0.370,   2.220,   0.370,   2.350,   0.370,   2.410,   0.370,   2.550, & 
!         0.370,   2.500,   0.370,   2.530,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.950, & 
!         0.370,   2.880,   0.370,   2.720,   0.370,   2.690,   0.370,   2.670,   0.370,   2.660, & 
!        -9.000,  -9.000,   0.370,   2.660,   0.370,   2.610,   0.370,   2.600,   0.370,   2.580, & 
!         0.370,   2.560,   0.370,   2.550,   0.370,   2.540,   0.370,   2.530,   0.370,   2.451, & 
!         0.370,   2.500,   0.370,   2.470,   0.370,   2.450,   0.370,   2.450,   0.370,   2.450, & 
!         0.370,   2.370,   0.370,   2.450,   0.370,   2.180,   0.370,   2.120,   0.370,   2.400, & 
!         0.370,   2.700,   0.370,   2.640,   0.370,   2.620,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.710, & 
!        -9.000,  -9.000,   0.370,   2.630,   0.400,   2.620,   0.400,   2.600,   0.400,   2.590, & 
!        -9.000,  -9.000,   0.400,   2.560,   0.400,   2.550/
!      data ( ( BVA(i1,36,i2),i1=1,2),i2=1,98) /  &  ! Kr
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   1.880,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,37,i2),i1=1,2),i2=1,98) /  &  ! Rb
!         0.370,   2.260,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.370,   0.370,   2.263,   0.370,   2.160,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.760, & 
!         0.370,   2.700,   0.370,   2.652,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.870,   0.370,   2.810,   0.370,   2.780, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   3.000,   0.370,   3.010,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,38,i2),i1=1,2),i2=1,98) /  &  ! Sr
!         0.370,   2.010,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.230,   0.370,   2.118,   0.370,   2.019,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.670, & 
!         0.370,   2.590,   0.370,   2.510,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.760,   0.370,   2.720,   0.370,   2.680, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.870,   0.370,   2.880,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,39,i2),i1=1,2),i2=1,98) /  &  ! Y 
!         0.370,   1.860,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.170,   0.370,   2.014,   0.370,   1.904,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.570, & 
!         0.370,   2.480,   0.370,   2.400,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.640,   0.370,   2.610,   0.370,   2.550, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.800,   0.370,   2.770,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,40,i2),i1=1,2),i2=1,98) /  &  ! Zr
!         0.370,   1.790,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.110,   0.370,   1.937,   0.370,   1.854,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.520, & 
!         0.370,   2.410,   0.370,   2.330,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.570,   0.370,   2.530,   0.370,   2.480, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.670,   0.370,   2.690,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,41,i2),i1=1,2),i2=1,98) /  &  ! Nb
!         0.370,   1.750,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.060,   0.370,   1.911,   0.370,   1.870,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.460, & 
!         0.370,   2.370,   0.370,   2.270,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.540,   0.370,   2.510,   0.370,   2.450, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.700,   0.370,   2.680,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,42,i2),i1=1,2),i2=1,98) /  &  ! Mo
!         0.370,   1.730,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.040,   0.370,   1.907,   0.370,   1.810,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.440, & 
!         0.370,   2.350,   0.370,   2.280,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.520,   0.370,   2.490,   0.370,   2.430, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.690,   0.370,   2.640,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,43,i2),i1=1,2),i2=1,98) /  &  ! Tc
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   1.900,   0.400,   1.880,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.210,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,44,i2),i1=1,2),i2=1,98) /  &  ! Ru
!         0.370,   1.610,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.880,   0.370,   1.834,   0.370,   1.740,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.290, & 
!         0.370,   2.160,   0.370,   2.210,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.360,   0.370,   2.330,   0.370,   2.260, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.540,   0.370,   2.480,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,45,i2),i1=1,2),i2=1,98) /  &  ! Rh
!         0.370,   1.550,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.880,   0.370,   1.793,   0.370,   1.710,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.290, & 
!         0.370,   2.150,   0.370,   2.170,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.370,   0.370,   2.330,   0.370,   2.250, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.550,   0.370,   2.480,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,46,i2),i1=1,2),i2=1,98) /  &  ! Pd
!         0.370,   1.470,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!         0.370,   1.730,   0.370,   1.810,   0.370,   1.792,   0.370,   1.740,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.220, & 
!         0.370,   2.100,   0.370,   2.050,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.300,   0.370,   2.220,   0.370,   2.190, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.480,   0.370,   2.380,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,47,i2),i1=1,2),i2=1,98) /  &  ! Ag
!         0.370,   1.500,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.850,   0.370,   1.805,   0.370,   1.800,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.220, & 
!         0.370,   2.119,   0.370,   2.090,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.300,   0.370,   2.260,   0.370,   2.220, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.510,   0.370,   2.380,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,48,i2),i1=1,2),i2=1,98) /  &  ! Cd
!         0.370,   1.660,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.960,   0.370,   1.904,   0.370,   1.811,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.340, & 
!         0.370,   2.304,   0.370,   2.230,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.430,   0.370,   2.400,   0.370,   2.350, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.590,   0.370,   2.570,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,49,i2),i1=1,2),i2=1,98) /  &  ! In
!         0.370,   1.720,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.030,   0.370,   1.902,   0.370,   1.792,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.430, & 
!         0.370,   2.360,   0.370,   2.280,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.350,   2.604, & 
!        -9.000,  -9.000,   0.350,   2.593,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.510,   0.370,   2.470,   0.370,   2.410, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.690,   0.370,   2.630,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,50,i2),i1=1,2),i2=1,98) /  &  ! Sn
!         0.370,   1.850,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.140,   0.370,   1.905,   0.370,   1.843,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.450, & 
!         0.370,   2.450,   0.370,   2.276,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.620,   0.370,   2.590,   0.370,   2.550, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.760,   0.370,   2.760,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,51,i2),i1=1,2),i2=1,98) /  &  ! Sb
!         0.370,   2.770,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.120,   0.370,   1.942,   0.370,   1.797,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.520, & 
!         0.370,   2.450,   0.370,   2.300,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.600,   0.370,   2.570,   0.370,   2.500, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.780,   0.370,   2.720,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,52,i2),i1=1,2),i2=1,98) /  &  ! Te
!         0.370,   1.830,  -9.000,  -9.000,   0.370,   2.300,   0.370,   2.210,   0.370,   2.200, & 
!         0.370,   2.210,   0.370,   2.120,   0.370,   1.917,   0.370,   1.820,  -9.000,  -9.000, & 
!         0.370,   2.640,   0.370,   2.530,   0.370,   2.480,   0.370,   2.490,   0.370,   2.520, & 
!         0.370,   2.450,   0.370,   2.300,  -9.000,  -9.000,   0.370,   2.930,   0.370,   2.760, & 
!         0.370,   2.640,   0.370,   2.600,   0.370,   2.570,   0.370,   2.520,   0.370,   2.550, & 
!         0.370,   2.530,   0.370,   2.460,   0.370,   2.430,   0.370,   2.270,   0.370,   2.450, & 
!         0.370,   2.540,   0.370,   2.560,   0.370,   2.600,   0.370,   2.530,   0.370,   2.530, & 
!        -9.000,  -9.000,   0.370,   3.000,   0.370,   2.870,   0.370,   2.800,   0.370,   2.670, & 
!         0.370,   2.700,   0.370,   2.690,  -9.000,  -9.000,   0.370,   2.540,   0.370,   2.550, & 
!         0.370,   2.480,   0.370,   2.510,   0.370,   2.590,   0.370,   2.690,   0.370,   2.760, & 
!         0.370,   2.780,   0.370,   2.760,   0.370,   2.760,  -9.000,  -9.000,   0.370,   3.160, & 
!         0.370,   3.080,   0.370,   2.940,   0.370,   2.920,   0.370,   2.900,   0.370,   2.890, & 
!        -9.000,  -9.000,   0.370,   2.860,   0.370,   2.850,   0.370,   2.840,   0.370,   2.820, & 
!         0.370,   2.800,   0.370,   2.800,   0.370,   2.780,   0.370,   2.770,   0.370,   2.760, & 
!         0.370,   2.750,   0.370,   2.720,   0.370,   2.700,   0.370,   2.710,   0.370,   2.700, & 
!        -9.000,  -9.000,   0.370,   2.710,   0.370,   2.450,   0.370,   2.410,   0.370,   2.610, & 
!         0.370,   2.930,   0.370,   2.840,   0.370,   2.870,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.940, & 
!        -9.000,  -9.000,   0.370,   2.860,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,53,i2),i1=1,2),i2=1,98) /  &  ! I 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.220,   0.370,   2.100,   0.370,   2.100, & 
!         0.370,   2.120,  -9.000,  -9.000,   0.370,   1.930,   0.370,   1.830,  -9.000,  -9.000, & 
!         0.370,   2.560,   0.370,   2.460,   0.370,   2.410,   0.370,   2.410,   0.370,   2.400, & 
!         0.370,   2.360,   0.370,   2.310,  -9.000,  -9.000,   0.370,   2.880,   0.370,   2.720, & 
!         0.370,   2.590,   0.370,   2.540,   0.370,   2.510,   0.370,   2.450,   0.370,   2.490, & 
!         0.370,   2.470,   0.350,   2.370,   0.370,   2.340,   0.370,   2.160,   0.370,   2.360, & 
!         0.370,   2.450,   0.370,   2.500,   0.370,   2.580,   0.370,   2.540,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   3.010,   0.370,   2.880,   0.370,   2.770,   0.370,   2.690, & 
!         0.370,   2.680,   0.370,   2.640,  -9.000,  -9.000,   0.370,   2.480,   0.370,   2.480, & 
!         0.370,   2.380,   0.370,   2.380,   0.370,   2.570,   0.370,   2.630,   0.370,   2.760, & 
!         0.370,   2.720,   0.370,   2.760,   0.350,   2.195,  -9.000,  -9.000,   0.370,   3.180, & 
!         0.370,   3.130,   0.370,   2.930,   0.370,   2.920,   0.370,   2.890,   0.370,   2.870, & 
!        -9.000,  -9.000,   0.370,   2.840,   0.370,   2.830,   0.370,   2.820,   0.370,   2.800, & 
!         0.370,   2.770,   0.370,   2.770,   0.370,   2.750,   0.370,   2.740,   0.370,   2.720, & 
!         0.370,   2.730,   0.370,   2.680,   0.370,   2.660,   0.370,   2.660,   0.370,   2.610, & 
!        -9.000,  -9.000,   0.370,   2.660,   0.370,   2.370,   0.370,   2.340,   0.370,   2.590, & 
!         0.370,   2.910,   0.370,   2.780,   0.370,   2.840,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.930, & 
!        -9.000,  -9.000,   0.370,   2.840,   0.400,   2.850,   0.400,   2.840,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,54,i2),i1=1,2),i2=1,98) /  &  ! Xe
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.000,   0.370,   1.890,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,55,i2),i1=1,2),i2=1,98) /  &  ! Cs
!         0.370,   2.440,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.530,   0.370,   2.417,   0.370,   2.330,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.930, & 
!         0.370,   2.890,   0.370,   2.791,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   3.040,   0.370,   2.980,   0.370,   2.950, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   3.160,   0.370,   3.180,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,56,i2),i1=1,2),i2=1,98) /  &  ! Ba
!         0.370,   2.220,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.470,   0.370,   2.285,   0.370,   2.188,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.880, & 
!         0.370,   2.769,   0.370,   2.690,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.960,   0.370,   2.880,   0.370,   2.880, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   3.080,   0.370,   3.130,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,57,i2),i1=1,2),i2=1,98) /  &  ! La
!         0.370,   2.060,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.340,   0.370,   2.172,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.730, & 
!         0.370,   2.643,   0.370,   2.545,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.800,   0.370,   2.740,   0.370,   2.720, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.940,   0.370,   2.930,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,58,i2),i1=1,2),i2=1,98) /  &  ! Ce
!         0.370,   2.040,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.340,   0.370,   2.151,   0.370,   2.036,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.700, & 
!         0.370,   2.620,   0.370,   2.410,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.780,   0.370,   2.740,   0.370,   2.690, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.920,   0.370,   2.920,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,59,i2),i1=1,2),i2=1,98) /  &  ! Pr
!         0.370,   2.020,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.300,   0.370,   2.138,   0.370,   2.022,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.680, & 
!         0.370,   2.600,   0.370,   2.500,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.750,   0.370,   2.720,   0.370,   2.670, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.900,   0.370,   2.890,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,60,i2),i1=1,2),i2=1,98) /  &  ! Nd
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.300,   0.370,   2.117,   0.370,   2.008,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!         0.370,   2.590,   0.370,   2.492,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.710,   0.370,   2.660, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.890,   0.370,   2.870,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,61,i2),i1=1,2),i2=1,98) /  &  ! Pm
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,62,i2),i1=1,2),i2=1,98) /  &  ! Sm
!         0.370,   1.960,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.240,   0.370,   2.088,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.630, & 
!         0.370,   2.550,   0.370,   1.977,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.700,   0.370,   2.670,   0.370,   2.660, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.860,   0.370,   2.840,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,63,i2),i1=1,2),i2=1,98) /  &  ! Eu
!         0.370,   1.950,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.240,   0.370,   2.074,   0.370,   1.961,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.620, & 
!         0.370,   2.530,   0.370,   2.530,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.700,   0.370,   2.660,   0.370,   2.610, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.850,   0.370,   2.830,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,64,i2),i1=1,2),i2=1,98) /  &  ! Gd
!         0.370,   1.930,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.220,   0.370,   2.065,   0.370,   1.950,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.610, & 
!         0.370,   2.530,   0.370,   2.445,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.680,   0.370,   2.650,   0.370,   2.600, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.840,   0.370,   2.820,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,65,i2),i1=1,2),i2=1,98) /  &  ! Tb
!         0.370,   1.910,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.200,   0.370,   2.049,   0.370,   1.936,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.590, & 
!         0.370,   2.510,   0.370,   2.427,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.660,   0.370,   2.630,   0.370,   2.580, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.820,   0.370,   2.800,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,66,i2),i1=1,2),i2=1,98) /  &  ! Dy
!         0.370,   1.890,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.180,   0.370,   2.001,   0.370,   1.922,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.570, & 
!         0.370,   2.470,   0.370,   2.410,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.640,   0.370,   2.610,   0.370,   2.560, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.800,   0.370,   2.770,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,67,i2),i1=1,2),i2=1,98) /  &  ! Ho
!         0.370,   1.880,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.180,   0.370,   2.025,   0.370,   1.908,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.560, & 
!         0.370,   2.480,   0.370,   2.401,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.640,   0.370,   2.610,   0.370,   2.550, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.800,   0.370,   2.770,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,68,i2),i1=1,2),i2=1,98) /  &  ! Er
!         0.370,   1.860,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.160,   0.370,   1.988,   0.370,   1.904,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.550, & 
!         0.370,   2.460,   0.370,   2.390,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.630,   0.370,   2.590,   0.370,   2.540, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.780,   0.370,   2.750,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,69,i2),i1=1,2),i2=1,98) /  &  ! Tm
!         0.370,   1.850,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.140,   0.370,   2.000,   0.370,   1.842,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.530, & 
!         0.370,   2.450,   0.370,   2.380,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.620,   0.370,   2.580,   0.370,   2.530, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.770,   0.370,   2.740,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,70,i2),i1=1,2),i2=1,98) /  &  ! Yb
!         0.370,   1.820,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.120,   0.370,   1.985,   0.370,   1.875,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.530, & 
!         0.370,   2.430,   0.370,   2.371,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.590,   0.370,   2.560,   0.370,   2.451, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.760,   0.370,   2.720,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,71,i2),i1=1,2),i2=1,98) /  &  ! Lu
!         0.370,   1.820,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.110,   0.370,   1.971,   0.370,   1.876,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.510, & 
!         0.370,   2.430,   0.370,   2.361,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.590,   0.370,   2.560,   0.370,   2.500, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.750,   0.370,   2.730,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,72,i2),i1=1,2),i2=1,98) /  &  ! Hf
!         0.370,   1.780,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.090,   0.370,   1.923,   0.370,   1.850,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.480, & 
!         0.370,   2.390,   0.370,   2.300,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.560,   0.370,   2.520,   0.370,   2.470, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.720,   0.370,   2.680,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,73,i2),i1=1,2),i2=1,98) /  &  ! Ta
!         0.370,   1.760,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.010,   0.370,   1.920,   0.370,   1.880,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.470, & 
!         0.370,   2.390,   0.370,   2.300,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.550,   0.370,   2.510,   0.370,   2.450, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.700,   0.370,   2.660,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,74,i2),i1=1,2),i2=1,98) /  &  ! W 
!         0.370,   1.760,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.060,   0.370,   1.921,   0.370,   1.830,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.460, & 
!         0.370,   2.390,   0.370,   2.270,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.540,   0.370,   2.510,   0.370,   2.450, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.710,   0.370,   2.660,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,75,i2),i1=1,2),i2=1,98) /  &  ! Re
!         0.370,   1.750,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.060,   0.370,   1.970,   0.370,   1.860,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.460, & 
!         0.370,   2.370,   0.370,   2.230,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.540,   0.370,   2.500,   0.370,   2.450, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.700,   0.370,   2.610,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,76,i2),i1=1,2),i2=1,98) /  &  ! Os
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   1.811,   0.370,   1.720,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!         0.370,   2.210,   0.370,   2.190,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.370, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,77,i2),i1=1,2),i2=1,98) /  &  ! Ir
!         0.370,   1.760,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.060,   0.370,   1.916,   0.370,   1.820,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.460, & 
!         0.370,   2.380,   0.370,   2.300,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.540,   0.370,   2.510,   0.370,   2.450, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.710,   0.370,   2.660,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,78,i2),i1=1,2),i2=1,98) /  &  ! Pt
!         0.370,   1.400,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!         0.370,   1.760,   0.370,   1.770,   0.370,   1.879,   0.370,   1.759,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.190, & 
!         0.370,   2.080,   0.370,   2.170,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.260,   0.370,   2.190,   0.370,   2.180, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.450,   0.370,   2.370,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,79,i2),i1=1,2),i2=1,98) /  &  ! Au
!         0.370,   1.370,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   1.720,   0.370,   1.833,   0.370,   1.810,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.140, & 
!         0.370,   2.030,   0.370,   2.170,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.220,   0.370,   2.180,   0.370,   2.120, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.410,   0.370,   2.340,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,80,i2),i1=1,2),i2=1,98) /  &  ! Hg
!         0.370,   1.710,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.020,   0.370,   1.930,   0.370,   1.900,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.420, & 
!         0.370,   2.320,   0.370,   2.250,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.500,   0.370,   2.470,   0.370,   2.400, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.610,   0.370,   2.590,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.350,   2.510, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,81,i2),i1=1,2),i2=1,98) /  &  ! Tl
!         0.370,   2.050,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.290,   0.370,   2.003,   0.370,   1.880,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.710, & 
!         0.370,   2.630,   0.370,   2.320,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.790,   0.370,   2.700,   0.370,   2.700, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.930,   0.370,   2.910,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,82,i2),i1=1,2),i2=1,98) /  &  ! Pb
!         0.370,   1.970,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.220,   0.370,   2.042,   0.370,   1.940,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.640, & 
!         0.370,   2.550,   0.370,   2.430,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.720,   0.370,   2.670,   0.370,   2.640, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.840,   0.370,   2.780,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,83,i2),i1=1,2),i2=1,98) /  &  ! Bi
!         0.370,   1.970,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.240,   0.370,   2.060,   0.370,   1.970,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.630, & 
!         0.370,   2.550,   0.370,   2.440,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.720,   0.370,   2.720,   0.370,   2.620, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.870,   0.370,   2.840,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,84,i2),i1=1,2),i2=1,98) /  &  ! Po
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.190,   0.370,   2.380,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,85,i2),i1=1,2),i2=1,98) /  &  ! At
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,86,i2),i1=1,2),i2=1,98) /  &  ! Rn
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,87,i2),i1=1,2),i2=1,98) /  &  ! Fr
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,88,i2),i1=1,2),i2=1,98) /  &  ! Ra
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,89,i2),i1=1,2),i2=1,98) /  &  ! Ac
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.240,   0.370,   2.130,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.630,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,90,i2),i1=1,2),i2=1,98) /  &  ! Th
!         0.370,   2.070,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.340,   0.370,   2.167,   0.370,   2.068,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.730, & 
!         0.370,   2.640,   0.370,   2.550,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.800,   0.370,   2.760,   0.370,   2.710, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.940,   0.370,   2.930,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,91,i2),i1=1,2),i2=1,98) /  &  ! Pa
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.350,   2.090,   0.370,   2.040,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,92,i2),i1=1,2),i2=1,98) /  &  ! U 
!         0.370,   1.970,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.240,   0.370,   2.112,   0.370,   2.034,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.640, & 
!         0.370,   2.560,   0.370,   2.460,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.720,   0.370,   2.700,   0.370,   2.630, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.860,   0.370,   2.840,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,93,i2),i1=1,2),i2=1,98) /  &  ! Np
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.180,   0.370,   2.020,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.400,   2.460,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.400,   2.620, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.400,   2.850,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,94,i2),i1=1,2),i2=1,98) /  &  ! Pu
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.110,   0.370,   2.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.480,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.400,   2.600, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.400,   2.840,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,95,i2),i1=1,2),i2=1,98) /  &  ! Am
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.110,   0.370,   2.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.480,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.400,   2.590, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,96,i2),i1=1,2),i2=1,98) /  &  ! Cm
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.230,   0.370,   2.120,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.620,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,97,i2),i1=1,2),i2=1,98) /  &  ! Bk
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.080,   0.370,   1.960,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.460,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.400,   2.560, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!      data ( ( BVA(i1,98,i2),i1=1,2),i2=1,98) /  &  ! Cf
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,   0.370,   2.070,   0.370,   1.950,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,   0.370,   2.450,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,   0.400,   2.550, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000, & 
!        -9.000,  -9.000,  -9.000,  -9.000,  -9.000,  -9.000/
!  if( (i1.gt.98).or.(i2.gt.98))then 
!     scale=-1
!     dist =-1
!     return
!  endif
!  scale=BVA(1,i1,i2) 
!  dist =BVA(2,i1,i2) 
!  if(scale .lt. 0)then 
!     scale=BVA(1,i2,i1)
!     dist =BVA(2,i2,i1)
!  endif
!  return
!end subroutine bvan

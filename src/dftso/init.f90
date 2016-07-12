! Next job
subroutine init(deffn,errfn,e,vru,emm,fl,jspin,kpot,ipr,th,ph,irlotot)

  USE param
  USE loabc; USE loabcr; USE lolog; USE rlolog; USE rpars; USE orb
  USE struct
  USE rotmat
  USE vns
  USE couples
  implicit real*8(a-h,o-z)
  !
  character*80    deffn,errfn,fname
  character*11    status,form
  character*5     vect
  common /gener/  br1(3,3),br2(3,3)
  logical         fl
  dimension       vru(nrad,nato,2),lm(2,9),iatoff(nato)
  dimension       emm(2),e(lmx,nato,2),xms(3)
  character*25 charop(3)
  character*5  charsp (-1:2)
  complex*16   vorb(-3:3,-3:3)
  
  parameter( pi     = 3.141592653589793238462643d0 )
  
  !rschmid
  !   Read parameters from unit 5,i.e. case.inso
  !     ipr     parameter for debug output.
  !     kpot    determines how the potential for spin-polarized
  !             calculation is determined.
  !rschmid
  
  read(5,549) vect
  read(5,*)labc,ipr,kpot
  labc2=(labc+1)**2
  
  if (vect.eq.'WFFIL') fl=.true.
  if (jspin.eq.2) then
     if (kpot.eq.1) then
        write(8,*)' Averaged potential used when calculating dV/dr'
        write(6,*)' Averaged potential used when calculating dV/dr'
     else
        write(8,*)' Potential not averaged when calculating dV/dr'
        write(6,*)' Potential not averaged when calculating dV/dr'
     endif
  endif
  
  
  if (fl) then
     write(6,*)'LMAX=',labc,' S-O eigenvectors and eigenvalues'
  else
     write(6,*)'LMAX=',labc,' S-O eigenvalues only'
  endif
  
  !rschmid
  !  Read energy ranges for which H_so is calculated and
  !  orientation for polarization axis for s-p-systems.
  !rschmid
  
  read(5,*) emm(1),emm(2)
  read(5,*) xms(1),xms(2),xms(3)
  
111 format(2f6.1,' angle (M,z), angle (M,x) deg')
  write(6,110) (emm(i),i=1,2)
110 format(' Emin=',f9.4,' Emax=',f9.4)
  
  if (emm(2) .lt. emm(1)) stop 'EMAX < EMIN'
  
  !rschmid
  !  read in parameters for basis extension by relativistic basis function.
  !rschmid
  do ity=1,nat
     nlr(ity)=0
  end do
  read(5,*,end=644)nrelcase
  do i=1,nrelcase 
     read(5,*) ity,edum1,edum2
     extei(ity,1)=edum1
     extde(ity,1)=edum2
     nlr(ity)=1 
  enddo
644 read(5,*,end=345)noff,(iatoff(i),i=1,noff)
  do i=1,noff
     lso(iatoff(i))=.false.
     write(6,*)'SOC on atom ',iatoff(i),'SWITCHED OFF !!!'
  enddo
345 continue
5000 format (1x,i1,2f10.5,a4)
  
  !rschmid
  ! sort entries for rlo extension by angular momentum.
  !rschmid
  do ity = 1, nat
     do l = 0, lomax
        loorext(l,ity) = .false.
        do j = 1, nlr(ity)
           loorext(1,ity) = .true.
           irlotot=(2*1+1)*mult(nat)
        enddo
     enddo
  enddo
  !
  !.....Generate symmop
  call init_couplo(ndif,labc)
  call latgen2
  call angle(xms,th,ph)
  nn=0
  do ity=1,nat
     nn=nn+mult(ity)
  end do
  call symop(th,ph)
  index=0
  do i=1,nat
     do j=1,mult(i)
	index=index+1
 !       if(det(index).gt.0.)then
 !write(6,557)i,j
 !        else
 !write(6,558)i,j
 !       endif
557 format(' for atom type',i3,' eq. atom',i3,' sym. oper. conserves spin')
558     format(' for atom type',i3,' eq. atom',i3,' sym. oper.   inverts spin')
     end do
  end do
  theta=th*180/pi
  phi=ph*180/pi
  write(8,111)theta,phi
  write(6,111)theta,phi 
  !
  mtap=17
  ntap=8
  do isi=1,jspin   ! 999
     mtap=mtap+1
     ntap=ntap+1
     read(mtap,532)
     nnlo=0
     nnrlo=0
     do ity=1,nat
        read(mtap,531)
        read(mtap,533)(vru(i,ity,isi),i=1,jrj(ity))
        read(mtap,531)
        read(ntap)(e(l,ity,isi),l=1,lmx)
        read(ntap)((elo(l,nn,ity,isi),l=0,lomax),nn=1,nloat)
        ! write also on dummy vector file 
        if(isi.eq.1)then
           write(53,'(100(f9.5))')(e(l,ity,isi),l=1,lmx) 
           write(53,'(100(f9.5))')((elo(l,nn,ity,isi),l=0,lomax),nn=1,nloat)  
        endif
        !rschmid
        !     Counting of rlo's included.
        !        nnlo    total number of local orbitals
        !        nlov    index-1 of the first local orbital belonging
        !                to atom type ity
        !        nnrlo   is nnlo for rlo's
        !        nrlov   is nlov for rlo's
        !rschmid
        DO l=1,lmax
           lapw(l-1,ity)=.true.
           IF(e(l,ity,isi).gt.150) THEN
              e(l,ity,isi)=e(l,ity,isi)-200.0D0
              lapw(l-1,ity)=.FALSE.
           ENDIF
        ENDDO
        
        DO l=0,lomax
           ilo(l,ity)=0
        ENDDO
        nlov(ity)=nnlo
        nrlov(ity) = nnrlo
        
        nlo(ity)  = 0
        nrlo(ity) = 0
        do l = 0,lomax
           loor(l,ity)=.false.
           DO k=1,nloat 
              if (elo(l,k,ity,isi).lt.(995.D+0)) then
                 !              loor(l,ity) = .true.
                 ilo(l,ity) = k
                 IF((lapw(l,ity).AND.ilo(l,ity).EQ.1).OR.ilo(l,ity).eq.2) loor(l,ity)=.true.
                 nnlo        =  nnlo       +((2*l+1))*mult(ity)
                 nlo(ity)    =  nlo(ity)   +((2*l+1))*mult(ity)
                 !           else
                 !              loor(l,ity) = .false.
              endif
           ENDDO
           if (loorext(l, ity)) then
              nnrlo          = nnrlo    +((2*l+1))*mult(ity)
              nrlo(ity)      = nrlo(ity)+((2*l+1))*mult(ity)
           endif
        enddo
     enddo
     
     !rschmid
     !  Determines the position of the first k-vector in klist associated
     !  with the local orbitals belonging to atom type ity.
     !rschmid
     do ity=1,nat
        nlon(ity)=nnlo-nlo(ity)-nlov(ity)
     enddo
  enddo !999  continue
  !  orbital potential part: start
  ! find out whether orb potential should be added (iorbpot=1)
  iorbpot=0
  CALL GTFNAM(DEFFN,ERRFN)
  OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=920)
10 CONTINUE
  READ (1,*,END=920,ERR=920) IUNIT
  if(iunit.eq.12)then    
     iorbpot=1
  endif
  goto 10
920 continue

  natorb=0
  natorb0=0
  nonso=0
  nmod=0
  if(iorbpot.eq.1)then
     ! read the orbital potential
     charop(1)=' LDA+U potential'
     charop(2)=' Orbital polarization'
     charop(3)=' Interaction with Bext'
     charsp(-1)=' down'
     charsp(0)= ' para'
     charsp(1)= ' up  '
     charsp(2)= ' dnup'
     do isi=1,3
        INDEX=0
        read(10+isi,*,end=555)nmod,nsp,natorb,Bext
        DO JAT = 1, NATORB
           read(10+isi,*)iatom,nlorb
           do nl=1,nlorb
              read(10+isi,*)l
              write(6,556)charop(nmod),iatom,l,charsp(nsp)
              write(8,556)charop(nmod),iatom,l,charsp(nsp)
556           format(a22,' added for atom type',i3,' L=',i3,' spin block',a5)
              do i=-l,l
                 do j=-l,l
                    read(10+isi,*) rval,cval
                    vorb(i,j)=cmplx(rval,cval)
103                 format(2f18.11)
                 enddo
              enddo
              call vorblo(index,iatom,l,isi,vorb)
              if(isi.eq.1)NATORB0=INDEX
102           format(' M=',i3,7f11.5)
              write(6,*)
              write(6,*)' Orbital potential real part'
              do m=-l,l
                 write(6,102)m,(dble(vorb(m,m1)),m1=-l,l)
              enddo
              write(6,*)
              write(6,*)' Orbital potential imaginary part'
              do m=-l,l
                 write(6,102)m,(dimag(vorb(m,m1)),m1=-l,l)
              enddo
              write(6,*)
           enddo
           ! end of orbital potential input
        endDO
        if(nonso.gt.0)then
           write(6,*)
           write(6,*) ' **** spin-orbit is excluded from the calculation ****'
           write(6,*)
           write(8,*)
           write(8,*) ' **** spin-orbit is excluded from the calculation ****'
           write(8,*)
        endif
555     continue
     enddo
     !  orbital potential part: end
  endif
  natorb=natorb0
  return
  
501 FORMAT(A10)     
531 FORMAT(/////)
532 FORMAT(//)
533 FORMAT((3X,4E19.12))
546 FORMAT(4I3)
547 FORMAT(3F10.4)
542 FORMAT(3f10.7)
541 FORMAT(a4,24x,i2)
543 FORMAT(12x,f10.7,3x,f10.7,3x,f10.7)
544 FORMAT(15x,i2)
545 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
549 FORMAT(A5)
553 FORMAT(2x,13f8.3)
1020 FORMAT(6F10.7,10X,F10.7)
5010 FORMAT(A80)
  !
END subroutine init

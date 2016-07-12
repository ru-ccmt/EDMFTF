SUBROUTINE HMSEC(FL,emm,ne,nv,ee,meigve,ri_mat,ri_orb,kkk,jspin,ipr)
  USE param
  USE loabc; USE loabcr; USE lolog; USE rlolog; USE rpars; USE orb
  USE struct
  USE peic
  USE hexpt
  USE abcd
  USE radovlp
  USE hmsout
  USE rotmat
  USE couples
  USE vns
  implicit real*8 (a-h,o-z)
  !
  !     THIS IS THE MAIN PROGRAM FOR SPIN-ORBIT INTERACTION. A modified
  !     version of the original B.N. Harmon program.
  !
  !     The input is the wavefunction outside Mt which is read in from file
  !     #10. Other inputs include file5(main input) and file11(potential
  !     which is *ONLY* spherical part of self-consistent flapw potential) and
  !     file20(the master structure file for flapw package)
  !
  !     RADIAL PART OF MATRIX ELEMENT OF Hs.o. IS CALCULATED IN GARADME.
  !     ANGULAR PART IS DONE BY COUPPLE. WAVEFUNCTION INSIDE MT (Alm and Blm)
  !     IS CALCULATED BY ABLM.   ............Jun Ye 1991
  !**********************************************************************
  dimension       emm(2)
  logical         FL
  dimension       ee(nume,2),ne(2),nv(2)
  complex*16      meigve(nmat,nume,2)
  real*8          ri_mat(4,4,0:labc,nato,2,2)
  real*8          ri_orb(4,4,0:labc,nato,2,2)
  character*1     job2
  !      COMPLEX*16      hso(num2)
  COMPLEX*16,allocatable ::  hso(:)
  !      complex*16      h_(nume2,nume2),s_(nume2,nume2)
  complex*16,allocatable ::  h_(:,:),s_(:,:)
  complex*16      ovrlc,hsoc,vnstot
  !      real*8          rwork(7*nume2)
  !      integer         iwork(5*nume2),ifail(nume2)
  !      complex*16      dwork((2+hblock)*nume2)
  !      complex*16      ovrful(2*nnrlo,2*nnrlo)
  real*8,allocatable ::          rwork(:)
  integer,allocatable ::         iwork(:),ifail(:)
  complex*16,allocatable ::      dwork(:)
  complex*16,allocatable ::      ovrful(:,:)
  complex*16      czero,imag,dd
  !      complex*16      vec(nume2,nume2)
  complex*16,allocatable ::  vec(:,:)
  
  real*8  cp(7)
  
  data   czero/(0.0d0,0.0d0)/,imag/(0.0d0,1.0d0)/
  
  !**********************************************************************
  !rschmid
  !  If required output for diagnostic purposes is provided.
  !rschmid

  if (ipr.gt.2) then
     write(6,*)'  ************ HMSEC **************'
     write(6,*)' l mlp ml    -,-    -,+    +,-    +,+'
     do l=1,2 
        do mlp=-l,l
           do ml =-l,l
              write(6,529) l,mlp,ml,(couplo(1,l,mlp,ml,1,ms),ms=1,2),(couplo(1,l,mlp,ml,2,ms),ms=1,2)
           enddo
        enddo
     enddo
     
     write(6,*)'*********** ALM,BLM,CLM **********'
     write(6,*)'TEST:',abcdlm(4,2,1,1,1)
     do i=1,ne(1)+nnrlo
        write(6,602) i,ee(i,1)
        l=2
        do m=-l,l
           ilm=l*l+l+m+1
           write(6,601) m,abcdlm(1,ilm,1,i,1),abcdlm(2,ilm,1,i,1),abcdlm(4,ilm,1,i,1)
        enddo
     enddo
  endif

529 format(3i3,4(2f7.3,1X))
602 format(' level',i3,' energy',f20.14)
601 format(i3,3(2f12.6,2x))

  PI=3.1415926535898D0
  TOL=1.D-6
  CALL CPUTIM(DTIME1)

  nban2=ne(1)+ne(2)+nnrlo*2
  n_scr=ne(1)+ne(2)
  nmatr=nban2*(nban2+1)/2
  
  allocate      (h_(nume2,nume2))
  if (nnrlo.eq.0) then
     allocate      (s_(1,1))
     allocate      (ovrful(1,1))
  else
     allocate      (s_(nume2,nume2))
     allocate      (ovrful(2*nnrlo,2*nnrlo))
  endif

  do ibi=1,nban2    ! 300
     DO ibf=1,ibi   ! 301
        hsoc  = czero
        ovrlc = czero
        if (ibi.le.n_scr) then
           if (ibi.gt.ne(1)) then
              iei=ibi-ne(1)
              isi=2
           else
              iei=ibi
              isi=1
           end if
        else
           if ((ibi-n_scr).gt.nnrlo) then
              iei=ibi-ne(1)-nnrlo
              isi=2
           else
              iei=ibi-ne(2)
              isi=1
           end if
        end if
        
        if (ibf.le.n_scr) then
           if (ibf.gt.ne(1)) then
              ief=ibf-ne(1)
              isf=2
           else
              ief=ibf
              isf=1
           end if
        else
           if ((ibf-n_scr).gt.nnrlo) then
              ief=ibf-ne(1)-nnrlo
              isf=2
           else
              ief=ibf-ne(2)
              isf=1
           end if
        end if
        
        if (jspin.eq.1) then
           isd=1
           isu=1
        else
           isd=isi
           isu=isf
        endif

        if (isf.ne.isi) goto 700
        
        if (ibi.le.n_scr) then
           ! For no extra p_{1/2} states, it is really simple
           if (ibi.eq.ibf)  then
              hsoc  = ee(iei,isd)
              ovrlc = (1.d0, 0.d0)
           end if
        else
           ! For p_{1/2} states, more complicated
           lfirst=1
           do  jatom=1,nat
              lfirst=lfirst+mult(jatom-1)
              do  l=1,1
                 if (.not.loorext(l,jatom)) cycle
                 do mi=-l,l                 
                    ilmi = l*l +l +mi +1
                    do ina=0,mult(jatom)-1    
                       latom=lfirst+ina
                       ! This is overlap when we add extra states p_{1/2}
                       ! Psi_{k,i} = Y_{lm}(r) ( alm ul + blm dot{ul} + clm u_LO)
                       ! Therefore <Psi_{ief}|Psi_{iei}> = alm^*_{ief} alm_{iei} + blm^*_{ief} blm_{iei}<udot|udot> + clm....
                       ovrlc = ovrlc + abcdlm(1,ilmi,latom,iei,isd)*dconjg(abcdlm(1,ilmi,latom,ief,isu)) &
                            +   abcdlm(2,ilmi,latom,iei,isd)*dconjg(abcdlm(2,ilmi,latom,ief,isu))*pei(l+1,jatom,isu) &
                            +   abcdlm(1,ilmi,latom,iei,isd)*dconjg(abcdlm(3,ilmi,latom,ief,isu))*ru2u(l,jatom,isu,isd) &
                            +   abcdlm(3,ilmi,latom,iei,isd)*dconjg(abcdlm(1,ilmi,latom,ief,isu))*ru2u(l,jatom,isu,isd) &
                            +   abcdlm(2,ilmi,latom,iei,isd)*dconjg(abcdlm(3,ilmi,latom,ief,isu))*ru2ud(l,jatom,isu,isd)&
                            +   abcdlm(3,ilmi,latom,iei,isd)*dconjg(abcdlm(2,ilmi,latom,ief,isu))*ru2ud(l,jatom,isu,isd)&
                            +   abcdlm(3,ilmi,latom,iei,isd)*dconjg(abcdlm(3,ilmi,latom,ief,isu))  &
                            +   abcdlm(1,ilmi,latom,iei,isd)*dconjg(abcdlm(4,ilmi,latom,ief,isu))*rup(l,jatom,isu,isd) &
                            +   abcdlm(4,ilmi,latom,iei,isd)*dconjg(abcdlm(1,ilmi,latom,ief,isu))*rup(l,jatom,isu,isd) &
                            +   abcdlm(2,ilmi,latom,iei,isd)*dconjg(abcdlm(4,ilmi,latom,ief,isu))*rudp(l,jatom,isu,isd)&
                            +   abcdlm(4,ilmi,latom,iei,isd)*dconjg(abcdlm(2,ilmi,latom,ief,isu))*rudp(l,jatom,isu,isd)&
                            +   abcdlm(4,ilmi,latom,iei,isd)*dconjg(abcdlm(3,ilmi,latom,ief,isu))*ru2p(l,jatom,isu,isd)&
                            +   abcdlm(3,ilmi,latom,iei,isd)*dconjg(abcdlm(4,ilmi,latom,ief,isu))*ru2p(l,jatom,isu,isd)&
                            +   abcdlm(4,ilmi,latom,iei,isd)*dconjg(abcdlm(4,ilmi,latom,ief,isu))*rrr (l,jatom,isu,isd)
                    enddo
                 enddo
              enddo
           enddo
           
           if (ibf.le.n_scr) then
              hsoc = ee(ief,isu)*ovrlc
           else
              lfirst=1
              do jatom=1,nat           ! 600
                 lfirst=lfirst+mult(jatom-1)
                 do l=1,1            ! 600
                    if (.not.loorext(l,jatom)) cycle
                    do ina=0,mult(jatom)-1
                       latom=lfirst+ina
                       do m1=-l,l   ! 620
                          ilmi = l*l +l + m1 +1
                          hsoc = hsoc + abcdlm(1,ilmi,latom,iei,isd)*dconjg(abcdlm(1,ilmi,latom,ief,isu))*hexl(l,jatom,isu,1) &
                                      + abcdlm(1,ilmi,latom,iei,isd)*dconjg(abcdlm(2,ilmi,latom,ief,isu))*hexl(l,jatom,isu,2) &
                                      + abcdlm(1,ilmi,latom,iei,isd)*dconjg(abcdlm(4,ilmi,latom,ief,isu))*hexl(l,jatom,isu,3) &
                                      + abcdlm(2,ilmi,latom,iei,isd)*dconjg(abcdlm(4,ilmi,latom,ief,isu))*hexl(l,jatom,isu,6) &
                                      + abcdlm(4,ilmi,latom,iei,isd)*dconjg(abcdlm(1,ilmi,latom,ief,isu))*hexl(l,jatom,isu,7) &
                                      + abcdlm(4,ilmi,latom,iei,isd)*dconjg(abcdlm(2,ilmi,latom,ief,isu))*hexl(l,jatom,isu,8) &
                                      + abcdlm(4,ilmi,latom,iei,isd)*dconjg(abcdlm(4,ilmi,latom,ief,isu))*hexl(l,jatom,isu,9)
                       enddo
                       if (lvns(jatom)) then
                          call hns(jatom,latom,isf,iei,ief,vnstot)
                          hsoc = hsoc +vnstot
                       end if
                    enddo
                 enddo
              enddo
           end if
        end if
700     continue 
        h_(ibf,ibi)=hsoc
        if (nnrlo.eq.0) CYCLE !goto 301
        s_(ibf,ibi)=ovrlc
        s_(ibi,ibf)=dconjg(ovrlc)
     ENDDO ! 301     CONTINUE
  enddo    ! 300  CONTINUE

  CALL CPUTIM(DTIME4)
  ! This is the crucial SO term, which adds l*s term inside MT-spheres.
  ! Uses precalculated array cuplo(latom,l,mf,mi,sf,si)
  call hsocalc(h_,ri_mat,ne,nnrlo,jspin)

  CALL CPUTIM(DTIME5)
  if (natorb.gt.0) call horbcalc(h_,ri_orb,ne,jspin)
  CALL CPUTIM(DTIME6)
  
  if (nnrlo.ne.0) then
     do i=1,2*nnrlo
        do j=1,2*nnrlo
           ovrful(i,j)=s_(n_scr+i,n_scr+j)
        end do
     end do
  endif
  !rschmid
  !   The hamiltonmatrix has been set up now.
  !rschmid
  
  CALL CPUTIM(DTIME2)
  
  !rschmid
  !    Provides diagnostic ouput of the Hamiltonmatrix.
  !rschmid
  
  if (ipr.gt.3) then
     write(6,*)' ******* real. part of h_ *****'
     do i=1,nban2
        write(6,568)(real(h_(i,j)),j=1,i)
     enddo
     write(6,*)' ******* imag. part of h_ *****'
     do i=1,nban2
        write(6,568)(dimag(h_(i,j)),j=1,i)
     enddo
     if (nnrlo.ne.0) then
        write(6,*)' ******* real  part of s_ *****'
        do i=n_scr+1,nban2
           write(6,568)(real(s_(i,j)),j=1,i)
        enddo
        write(6,*)' ******* imag. part of s_ *****'
        do i=n_scr+1,nban2
           write(6,568)(dimag(s_(i,j)),j=1,i)
        enddo
     endif
  endif
  
568 format(12f10.6)

  call cputim(cp(2))
  call cputim(cp(1))
  
  if (nnrlo.eq.0) goto 555
  
  !**************************************************************
  !.....SOLVE THE SECULAR EQUATION HSO*VEC=EN*VEC
  !
  !....computes the subblock of Cholesky transformation of the
  !....overlap matrix
  M=2*nnrlo
  K=n_scr

  !>>>>>>>>> start CHOLESKY
  cp=0
  lda=nume2
  call zgemm('N','C',M,M,K,(-1.d0,0.d0),s_(n_scr+1,1),lda,s_(n_scr+1,1),lda,(1.d0,0.d0),s_(n_scr+1,n_scr+1),lda)
  if(.not.fl) then
     job2='N'
  else
     job2='V'
  endif
  il=0
  iu=0
  abstol=0.0d0
  call zpotrf('L',M,s_(n_scr+1,n_scr+1),lda,info)
  !>>>>>>>>> end CHOLESKY
  call cputim(cp(2))
  !>>>>>  H~_ba=H_ba-S_ba*H_aa
  call zgemm('N','N',M,K,K,(-1d0,0.d0),s_(n_scr+1,1),lda,h_,lda,(1.d0,0.d0),h_(n_scr+1,1),lda)
  !>>>>>
  call zgemm('C','C',M,M,K,(-1.d0,0.d0),h_(1,n_scr+1),lda,s_(n_scr+1,1),lda,(1.d0,0.d0),h_(n_scr+1,n_scr+1),lda)
  call zgemm('N','C',M,M,K,(-1.d0,0.d0),s_(n_scr+1,1),lda,h_(n_scr+1,1),lda,(1.d0,0.d0),h_(n_scr+1,n_scr+1),lda)
  call ztrtri('L','N',M,s_(n_scr+1,n_scr+1),lda,info)
  !>>>>>  H`_ba=(U^+_bb)^-1*H~_ba
  call ztrmm('L','L','N','N',M,K,(1.d0,0.d0),s_(n_scr+1,n_scr+1),lda,h_(n_scr+1,1),lda)
  !>>>>>
  
  !>>>>>  H`_bb=(U^+_bb)^-1*(H_bb-S_ba*H~^+_ba-H_ba*S_ba^+)*U^-1
  call ztrmm('L','L','N','N',M,M,(1.d0,0.d0),s_(n_scr+1,n_scr+1),lda,h_(n_scr+1,n_scr+1),lda)
  call ztrmm('R','L','C','N',M,M,(1.d0,0.d0),s_(n_scr+1,n_scr+1),lda,h_(n_scr+1,n_scr+1),lda)
  !>>>>>

555 continue
  if(.not.fl) then
     job2='N'
  else
     job2='V'
  endif
  il=0
  iu=0
  abstol=0.0d0

!!! Here we create matrix hso, which contains lower triangular part of the matrix only, but it conjugates it, hence it is
!!! actually like upper triangular part.
  allocate (hso(num2))
  index=0
  do ibi=1,nban2
     do ibf=1,ibi
        index=index+1
        hso(index)=dconjg(h_(ibi,ibf))
     end do
  end do
  deallocate (h_)
  allocate (vec(nume2,nume2))
  allocate (rwork(7*nume2),iwork(5*nume2),ifail(nume2),dwork((2+hblock)*nume2))
  call cputim(cp(3))
  
  call zhhevx(job2,'V','U',nban2,hso,emm(1),emm(2),il,iu,abstol,neig,en,vec,nume2,dwork,rwork,iwork,ifail,hblock,nume2,nbelw,info)
  deallocate   (hso)
  deallocate   (rwork,iwork,ifail,dwork)
  call cputim(cp(4))
  if (nnrlo.eq.0) goto 1555
  !rschmid
  !  backsubstitution  
  !rschmid
  
  call ztrmm('L','L','C','N',M,neig,(1.d0,0.d0),s_(n_scr+1,n_scr+1),lda,vec(n_scr+1,1),lda)
  ldc=nume2
  call zgemm('C','N',k,neig,m,(-1.d0,0.d0),s_(n_scr+1,1),lda,vec(n_scr+1,1),lda,(1.d0,0.d0),vec(1,1),lda)

1555 continue
  do is = 1, 2
     do i = 1, nv(is)+nnrlo
        do ibi = 1,neig
           vect(i,ibi,is) = (0.d0,0.d0)
        enddo
     enddo
  enddo

!!! Here we calculate the eigenvector in (K,s) space. Namely, the eigenvector for SO was calculated in
!!! band space (vec). Now we back transform to (K,s) space, so that an eigenvector has a form
!!!   A(K,s,iband) == vect
!!!   
  do is=1,2
     if (jspin.eq.2) then
        ism=is
     else
        ism=1
     end if
     ! vect(K,iband_s,1) = A(K,iband1) * vec(1:iband1,       iband_is) 
     ! vect(K,iband_s,2) = A(K,iband1) * vec(iband1:2*iband1,iband_is)
     ! here iband_s is the band index in double-hilbert space
     call zgemm('N','N',nv(is),neig,ne(is),(1.d0,0.d0),meigve(1,1,ism),nmat,vec(1+ne(1)*(is -1),1),nume2,(1.d0,0.d0),vect(1,1,is),nmat)
     do j=1,nnrlo
        do ibi=1, neig
           nf = ne(1)+ne(2)+(is-1)*nnrlo
           vect(nv(is)+j,ibi,is) = vect(nv(is)+j,ibi,is)+(vec(nf+j,ibi))
        enddo
     enddo
  enddo
  call cputim(cp(5))
  
  !  calculate the norm of the spin parts of vector
  !      if (ipr.ge.1) then
  if (nnrlo.ne.0) then
     do i=1,2*nnrlo
        do j=1,2*nnrlo
           s_(n_scr+i,n_scr+j)=ovrful(i,j)
        end do
     end do
  endif
  
  do j=1,neig
     vnorm(j,1) = 0.d0
     vnorm(j,2) = 0.d0
     
     do i=1,nban2
        if (i.gt.ne(1).and.i.le.n_scr.or.i.gt.(n_scr+nnrlo)) then
           is = 2
        else
           is = 1
        endif
        if (nnrlo.eq.0) then
           vnorm(j,is) = vnorm(j,is) + dconjg(vec(i,j))*vec(i,j)
        else
           do k=1,nban2
              !          if(j.eq.1.and.abs(s_(i,k)).gt.1.d-10) print*, i,k,s_(i,k)
              vnorm(j,is) = vnorm(j,is)+dconjg(vec(i,j))*s_(i,k)*vec(k,j)
           enddo
        endif
     enddo
  enddo
  !       end if
  
  deallocate  (s_,ovrful)
  deallocate  (vec)
  call cputim(cp(6))
  
  do j = 1,5
     cp(j) = cp(j+1)-cp(j)
  enddo
  WRITE (6,1001) 'Cholesky complete' , CP(1) 
  WRITE (6,1001) 'Transform to eigenwertproblem' , CP(2) 
  WRITE (6,1001) 'Compute eigenvalues' , CP(3) 
  WRITE (6,1001) 'Backtransform' , CP(4) 
  WRITE (6,1001) 'norm' , CP(5) 
1001 FORMAT (1X,'Seclr4(',A,') :', t50, f9.3)
  
  CALL CPUTIM(DTIME3)
  write(6,6020)dtime2-dtime1,dtime3-dtime2
  write(6,6021)dtime5-dtime4,dtime6-dtime5
  !
  !**********************************************************************
  !
  RETURN
  
  !
  
201 format(i4,4e16.8)
202 format(4e20.12)
502 FORMAT(20I4)
530 FORMAT(8(7X,5F13.7/))
531 FORMAT(8(2X,5F13.7/))
533 FORMAT((3X,4E19.12))
1818 FORMAT(12f10.5)
6010 FORMAT(I13,' EIGENVALUES BELOW THE ENERGY ',F10.5)
6020 FORMAT(7X,'CPTIME HAMILT=',F9.3,', DIAG=',F9.3/)
6021 FORMAT(7X,'CPTIME Hsoset=',F9.3,', Horb=',F9.3/)
6030 FORMAT(7X,14('****'))
  !
END SUBROUTINE HMSEC

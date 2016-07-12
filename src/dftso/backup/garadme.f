subroutine garadme(e,vru,p,dp,pe,dpe,ri_mat,jspin,kpot,ipr)
  USE param
  USE loabc; USE loabcr; USE lolog; USE rlolog; USE rpars
  USE struct
  USE hexpt
  USE peic
  USE radovlp
  USE vns
  implicit real*8(a-h,o-z)
  !rschmid
  !   The overlapp matrix and the radial matrix elements of the spin-orbit
  !   term as given in Koelling, Harmon (J. Phys. C (1976) ????) are calculated.
  !   A basis extension for an improved representation of p_{1/2} has been included.  
  !
  !   Input:
  !        e
  !        vru          potential * r in Rydberg.  
  !        p
  !        dp
  !        dpe
  !        zz(nato)     nuclear charge of the different atoms
  !        jspin        = 1 for unspinpolarized calculations
  !                     = 2 for spin-polarized calculations
  !         (Note that the present version is not intended for spin-polarized calculations,
  !          as some of the programming for this case might be wrong. )
  !        kpot           
  !        ipr          = debugging option
  !
  !   Output
  !        ri_mat      matrix containing the radial overlap integrals.
  !
  !rschmid   
  real*8 ri_mat(4,4,0:labc,nato,2,2)
  real*8 vrr(nrad,2),vder(nrad,2),rx(nrad)
  real*8 ameaam(NRAD,2), amebbm(NRAD,2), ameabm(NRAD,2)
  real*8 ameadm(NRAD,2), amebdm(NRAD,2), amecdm(NRAD,2)
  real*8 amedam(NRAD,2), amedbm(NRAD,2), amedcm(NRAD,2)
  real*8 ameddm(NRAD,2), ameccm(NRAD,2), ameacm(NRAD,2)
  real*8 amebcm(NRAD,2), amebam(NRAD,2), amecam(NRAD,2)
  real*8 amecbm(NRAD,2)

  real*8 bup   (nrad,2)
  real*8 budp  (nrad,2)
  real*8 bu2ud (nrad,2)
  real*8 bu2p  (nrad,2)
  real*8 bu2u  (nrad,2)
  real*8 brr   (nrad,2)
  real*8 bupd  (nrad,2)
  real*8 budpd (nrad,2)
  real*8 bu2pd (nrad,2)
  real*8 bpdpd (nrad,2)
  real*8 bppd  (nrad,2)
  
  real*8 tota1(3), tota2(3), tota3(3)
  real*8 tota4(3), tota5(3), tota6(3)
  real*8 tota7(3), tota8(3), tota9(3)
  real*8 tota10(3),tota11(3)
  
  real*8 vru,e
  DIMENSION       VRU(NRAD,NATO,2),E(LMX,NATO,2)
  
  integer jspin,kpot,ipr
  
  real*8 p,dp,pe,dpe
  DIMENSION       P(Labc+1,NATO,2),DP(Labc+1,NATO,2),PE(Labc+1,NATO,2),DPE(Labc+1,NATO,2)
  
  real*8 radf  (nrad,labc+1,2,2)
  real*8 radfde(nrad,labc+1,2,2)
  real*8 radflo(nrad,lomax+1,2,2)
  real*8 radfrl(nrad,lomax+1,2,2)
  
  
  integer i,ic
  integer isj,isi
  integer jri,ity
  integer ir,lk,l,lpot
  integer lspin
  
  real*8 coc,coe,cin
  real*8 rmt2,dx2
  real*8 dnomi,dnomj,dnom,dnom1,dnom2
  real*8 vde
  
  CALL init_peic(labc+1,nato)
  CALL init_hexpt(lomax,nato)
  CALL init_radovlp(labc,nato)
  !rschmid
  !    coc  : velocity of light in Rydberg units
  !    cin  : inverse velocity of light squared. (Hartree units)
  !rschmid 
  COC=2.d0*CLIGHT
  COE=2.D0
  CIN=1.d0/CLIGHT**2

  ri_mat(:,:,:,:,:,:)=0.d0
  !do 1 il=0,labc
  !   do 1 ia=1,nato
  !      do 1 i1=1,2
  !         do 1 i2=1,2
  !            do 335 j1=1,4
  !               do 335 j2=1,4
  !                  ri_mat(j1,j2,il,ia,i1,i2) = 0.d0
  !335   continue
  !1     continue

  !rschmid 
  !  ity is a loop over all non-equivalent atoms.
  !rschmid
  DO ITY=1,nat        ! 109
     DX2  = DX(ITY)
     JRI = JRJ(ITY)
     RMT2= RMT(ITY)
     !rschmid
     !  The logarithmic mesh is setup in RX for the respective atom.
     !rschmid
     do i = 1, jri
        rx(i) = rmt2 * dexp(dfloat(i-jri)*dx2)
     enddo
     ! if lpot=3, averaged potential is used in spinpolarized calculations
     ! when calculating dV/dr (kpot=1, jspin=2)
     lpot=jspin+kpot
     if (lpot.eq.3) then
        lspin=1
     else
        lspin=jspin
     endif
     
     !rschmid
     !  Setup potential ( V*r -> V ), do averaging and calculate
     !  derivative of potential
     !rschmid
     call vderiv(lpot,ity,lspin,rx,vru,vrr,vder,jri)
     do isi=1,jspin
        call atpar(e,vru,labc+1,ity,radf,radfde,radflo,radfrl,p,dp,pe,dpe,pei,isi)
        !honza contribution of V(2,m) to <p1/2|H|p1/2>
        if (nrlo(ity).gt.0) call vnsrint(ity,isi,radf,radfde,radfrl)
     enddo
     write(6,*)
     !rschmid
     !  Calculate some radial integrals needed for the construction
     !  of the overlap matrix.
     !rschmid
     do isi=1,jspin
        do isj = 1,jspin
           do l=0,labc
              lk = l+1
              do ic=1,2
                 do ir=1,jri
                    if (l .le. lomax) then
                       if (loorext(l,ity)) then
                          bup(ir,ic)  = radf  (ir,lk,ic,isi)*radfrl(ir,lk,ic,isj)
                          budp(ir,ic) = radfde(ir,lk,ic,isi)*radfrl(ir,lk,ic,isj)
                          brr(ir,ic)  = radfrl(ir,lk,ic,isi)*radfrl(ir,lk,ic,isj)  
                       else
                          bup(ir,ic)   = 0.d0
                          budp(ir,ic)  = 0.d0
                          brr(ir,ic)   = 0.d0
                          bupd(ir,ic)  = 0.d0
                          budpd(ir,ic) = 0.d0
                          bpdpd(ir,ic) = 0.d0
                          bppd(ir,ic)  = 0.d0
                       endif
                    else
                       bup(ir,ic)   = 0.d0 
                       budp(ir,ic)  = 0.d0
                       brr(ir,ic)   = 0.d0
                       bupd(ir,ic)  = 0.d0
                       budpd(ir,ic) = 0.d0
                       bpdpd(ir,ic) = 0.d0
                       bppd(ir,ic)  = 0.d0
                    endif
                    if (l .le. lomax) then
                       if (loor(l,ity)) then  
                          bu2ud(ir,ic)= radflo(ir,lk,ic,isi)*radfde(ir,lk,ic,isj)
                          bu2u(ir,ic) = radflo(ir,lk,ic,isi)*radf(ir,lk,ic,isj)
                       else
                          bu2ud(ir,ic) = 0.d0
                          bu2u(ir,ic) = 0.d0
                       endif
                    else
                       bu2ud(ir,ic) = 0.d0
                       bu2u(ir,ic) = 0.d0 
                    endif
                    if (l .le. lomax) then
                       if (loor(l,ity) .and. loorext(l,ity)) then
                          bu2p (ir,ic) = radflo(ir,lk,ic,isi)*radfrl(ir,lk,ic,isj)
                       else
                          bu2p (ir,ic) = 0.d0
                          bu2pd(ir,ic) = 0.d0
                       endif
                    else
                       bu2p(ir,ic)  = 0.d0
                       bu2pd(ir,ic) = 0.d0
                    endif
                 enddo
                 call cali(  bup(1,ic),tota1(ic) ,rx(1),dx2,jri)
                 call cali( budp(1,ic),tota2(ic) ,rx(1),dx2,jri)
                 call cali(bu2ud(1,ic),tota3(ic) ,rx(1),dx2,jri)
                 call cali( bu2p(1,ic),tota4(ic) ,rx(1),dx2,jri)
                 call cali( bu2u(1,ic),tota5(ic) ,rx(1),dx2,jri)
                 call cali( brr (1,ic),tota6(ic) ,rx(1),dx2,jri)
                 call cali( bupd(1,ic),tota7(ic) ,rx(1),dx2,jri)  
                 call cali(budpd(1,ic),tota8(ic) ,rx(1),dx2,jri)  
                 call cali(bu2pd(1,ic),tota9(ic) ,rx(1),dx2,jri)  
                 call cali(bpdpd(1,ic),tota10(ic),rx(1),dx2,jri)  
                 call cali(bppd(1,ic) ,tota11(ic),rx(1),dx2,jri)
              enddo
              !rschmid
              !   Note that ic = 2 is not the whole expression for
              !   the small component. If eq. 10 of Koelling, Harmon
              !   is used for the overlapp matrix and S_l is discarded
              !   an additional angular momentum dependent contribution
              !   has to be added.
              !   This contribution amouts to approx. 0.05 mRyd for p states
              !   in Au. (It has not been taken into account in other parts
              !   of the WIEN package.) For excited states there are contributions
              !   up to 0.1 mRyd.  
              !rschmid
              !            ic = 1
              !            do ir=1,jri
              !             fac = l * (l+1) / rx(ir) / rx(ir)
              !               bup(ir,ic)  =
              !     c          fac*radf  (ir,lk,ic,isi)*radfrl(ir,lk,ic,isj)
              !     c         / ( coc + (e   (lk,ity,isi) - vrr(ir,isi))/coc)
              !     c         / ( coc + (elor(l ,ity,isj) - vrr(ir,isj))/coc)      
              !
              !               budp(ir,ic) =
              !     c          fac*radfde  (ir,lk,ic,isi)*radfrl(ir,lk,ic,isj)
              !     c        / ( coc + (e   (lk,ity,isi) - vrr(ir,isi))/coc)
              !     c        / ( coc + (elor(l ,ity,isj) - vrr(ir,isj))/coc) 
              !
              !                bu2ud(ir,ic) =
              !     c         fac*radflo(ir,lk,ic,isi)*radfde  (ir,lk,ic,isj)
              !     c        / ( coc + (elo (l ,ity,isi) - vrr(ir,isi))/coc)
              !     c        / ( coc + (e   (lk,ity,isj) - vrr(ir,isj))/coc) 
              !                bu2u (ir,ic) =
              !     c         fac*radflo(ir,lk,ic,isi)*radf  (ir,lk,ic,isj)
              !     c        / ( coc + (elo (l ,ity,isi) - vrr(ir,isi))/coc)
              !     c        / ( coc + (e   (lk,ity,isj) - vrr(ir,isj))/coc) 
              !
              !             bu2p(ir,ic) =
              !     c            fac*radflo(ir,lk,ic,isi)*radfrl(ir,lk,ic,isj)
              !     c        / ( coc + (elo (l ,ity,isi) - vrr(ir,isi))/coc)
              !     c        / ( coc + (elor(l ,ity,isj) - vrr(ir,isj))/coc) 
              !
              !           enddo
              !           call cali(  bup(1,ic),tota1(3),rx(1),dx,jri)
              !           call cali( budp(1,ic),tota2(3),rx(1),dx,jri)
              !           call cali(bu2ud(1,ic),tota3(3),rx(1),dx,jri)
              !           call cali( bu2p(1,ic),tota4(3),rx(1),dx,jri)
              !           call cali( bu2u(1,ic),tota5(3),rx(1),dx,jri)
              !          else
              tota1 (3) = 0.d0
              tota2 (3) = 0.d0
              tota3 (3) = 0.d0
              tota4 (3) = 0.d0
              tota5 (3) = 0.d0     
              tota6 (3) = 0.d0
              tota7 (3) = 0.d0
              tota8 (3) = 0.d0
              tota9 (3) = 0.d0
              tota10(3) = 0.d0
              tota11(3) = 0.d0
              !          endif
              
              
              rup(  l,ity,isi,isj)= tota1(1) + sqrt(cin)*tota1(2) + tota1(3)
              rudp( l,ity,isi,isj)= tota2(1) + sqrt(cin)*tota2(2) + tota2(3)
              ru2ud(l,ity,isi,isj)= tota3(1) +      cin *tota3(2) + tota3(3)
              ru2p( l,ity,isi,isj)= tota4(1) + sqrt(cin)*tota4(2) + tota4(3)
              ru2u( l,ity,isi,isj)= tota5(1) +      cin *tota5(2) + tota5(3)
              rrr ( l,ity,isi,isj)= tota6(1) + tota6(2) + tota6(3)
              rupd ( l,ity,isi,isj) = tota7(1)+ sqrt(cin)*tota7(2)+tota7(3)
              rudpd( l,ity,isi,isj) = tota8(1)+ sqrt(cin)*tota8(2)+tota8(3)
              ru2pd( l,ity,isi,isj) = tota9(1)+ sqrt(cin)*tota9(2)+tota9(3)
              rpdpd( l,ity,isi,isj) = tota10(1) + tota10(2) + tota10(3)
              rppd(  l,ity,isi,isj) = tota11(1) + tota11(2) + tota11(3)
           enddo
           !         //  l  //
        enddo
        !       //  isj //
     enddo
     !     // isi //
     !
     do isi=1,jspin     ! 91
        do isj=1,jspin  ! 91
           !.........................................................................
           !.....CALCULATE THE MATRIX ELEMENT
           !
           DO L=0,LABC          ! 100
              LK=L+1
              DO IC=1,1         ! 104
                 DO IR=1,JRI    ! 160
                    !....      DNOM=(COC+(E(LK,ITY)-VRR(IR)/COC))**2
                    DNOMi=(COC+(E(LK,ITY,isi)-VRR(IR,isi))/COC)
                    DNOMj=(COC+(E(LK,ITY,isj)-VRR(IR,isj))/COC)
                    dnom=dnomi*dnomj
                    
                    if (isi.eq.isj) then
                       vde=vder(ir,isi)
                    else
                       vde=(vder(ir,isi)+vder(ir,isj))/2.
                    endif
                    
                    AMEAAM(IR,IC)=RADF(IR,LK,IC,isi)*RADF(IR,LK,IC,isj)*vde/RX(IR)/DNOM
                    AMEBBM(IR,IC)=RADfde(IR,LK,IC,isi)*RADFDE(IR,LK,IC,isj)*vde/RX(IR)/DNOM
                    AMEABM(IR,IC)=RADF(IR,LK,IC,isi)*RADFDE(IR,LK,IC,isj)*vde/RX(IR)/DNOM
                    AMEBAM(IR,IC)=RADFDE(IR,LK,IC,isi)*RADF(IR,LK,IC,isj)*vde/RX(IR)/DNOM
                 ENDDO !106 CONTINUE
                 !
                 CALL CALI(AMEAAM(1,IC),TOTAA,RX(1),DX2,JRI)
                 CALL CALI(AMEBBM(1,IC),TOTBB,RX(1),DX2,JRI)
                 CALL CALI(AMEABM(1,IC),TOTAB,RX(1),DX2,JRI)
                 CALL CALI(AMEBAM(1,IC),TOTBA,RX(1),DX2,JRI)
              ENDDO ! 104              CONTINUE
              !
              ri_mat(1,1,L,ITY,isi,isj)=COE*(TOTAA)
              ri_mat(2,2,L,ITY,isi,isj)=COE*(TOTBB)
              ri_mat(1,2,L,ITY,isi,isj)=COE*(TOTAB)
              ri_mat(2,1,L,ITY,isi,isj)=COE*(TOTBA)
           ENDDO  !100           CONTINUE
           !
           !.........................................................................
           !.....CALCULATE THE MATRIX ELEMENT of LO
           !
           DO L=0,lomax            ! 200
              if (loor(l,ity)) then
                 LK=L+1
                 DO IC=1,1          ! 204
                    DO IR=1,JRI     ! 206
                       if(isi.eq.isj)then
                          vde=vder(ir,isi)
                       else
                          vde=(vder(ir,isi)+vder(ir,isj))/2.
                       endif
                       if(ir.eq.1.and.l.eq.lomax) write(882,*) ir   ! stupid fix for ifort12 bug
                       DNOM= (COC+(ELO(L,ilo(l,ity),ity,isj)-VRR(IR,isj))/COC)*(COC+(E (LK,ity,isi)-VRR(IR,isi))/COC)
                       DNOM2= (COC+(ELO(L,ilo(l,ity),ity,isi)-VRR(IR,isi))/COC)*(COC+(E (LK,ity,isj)-VRR(IR,isj))/COC)
                       DNOM1=(COC+(ELO(L,ilo(l,ity),ity,isi)-VRR(IR,isi))/COC)*(COC+(ELO(L,ilo(l,ity),ity,isj)-VRR(IR,isj))/COC)
                       ameccm(ir,ic)= radflo(ir,lk,ic,isi)*radflo(ir,lk,ic,isj)*vde/rx(ir)/dnom1
                       ameacm(ir,ic)= radf  (ir,lk,ic,isi)*radflo(ir,lk,ic,isj)*vde/rx(ir)/dnom
                       amebcm(ir,ic)=radfde(ir,lk,ic,isi)*radflo(ir,lk,ic,isj)*vde/rx(ir)/dnom
                       amecam(ir,ic)= radflo  (ir,lk,ic,isi)*radf(ir,lk,ic,isj)*vde/rx(ir)/dnom2
                       amecbm(ir,ic)=radflo(ir,lk,ic,isi)*radfde(ir,lk,ic,isj)*vde/rx(ir)/dnom2
                    ENDDO  !206                    CONTINUE
                    CALL CALI(AMECCM(1,IC),TOTCC,RX(1),DX2,JRI)
                    CALL CALI(AMEACM(1,IC),TOTAC,RX(1),DX2,JRI)
                    CALL CALI(AMEBCM(1,IC),TOTBC,RX(1),DX2,JRI)
                    CALL CALI(AMECAM(1,IC),TOTCA,RX(1),DX2,JRI)
                    CALL CALI(AMECBM(1,IC),TOTCB,RX(1),DX2,JRI)
                 ENDDO             ! 204   CONTINUE
                 !
                 ri_mat(3,3,L,ITY,isi,isj)=COE*TOTCC
                 ri_mat(1,3,L,ITY,isi,isj)=COE*TOTAC
                 ri_mat(3,1,L,ITY,isi,isj)=COE*TOTCA
                 ri_mat(2,3,L,ITY,isi,isj)=COE*TOTBC
                 ri_mat(3,2,L,ITY,isi,isj)=COE*TOTCB
              endif
           ENDDO ! 200 CONTINUE
           !rschmid
           !   Compute matrix elements for basis extension.
           !rschmid
           do l=0,lomax
              if (loorext(l,ity)) then
                 lk = l+1
                 do ic=1,1
                    do ir = 1,jri
                       if (isi.eq.isj) then
                          vde = vder(ir,isi)
                       else
                          vde = (vder(ir,isi)+vder(ir,isj))/2.
                       endif
                       if(ir.eq.1.and.l.eq.lomax) write(882,*) ir   ! stupid fix for ifort12 bug
                       dnom = (coc+(elor2(l,ity,isj)-vrr(ir,isj))/coc)*(coc+(e(lk,ity,isi)-vrr(ir,isi))/coc)
                       dnom21 = (coc+(elor2(l,ity,isi)-vrr(ir,isi))/coc)*(coc+(e(lk,ity,isj)-vrr(ir,isj))/coc)
                       dnom1= (coc+(elor2(l,ity,isj)-vrr(ir,isj))/coc)*(coc+(elor2(l,ity,isi)-vrr(ir,isi))/coc)
                       dnom2 =(coc+(elo(l,ilo(l,ity),ity,isj)-vrr(ir,isj))/coc)*(coc+(elor2(l,ity,isi)-vrr(ir,isi))/coc)
                       dnom12 =(coc+(elo(l,ilo(l,ity),ity,isi)-vrr(ir,isi))/coc)*(coc+(elor2(l,ity,isj)-vrr(ir,isj))/coc)
                       ameadm(ir,ic) = radf  (ir,lk,ic,isi)*radfrl(ir,lk,ic,isj)*vde/rx(ir)/dnom
                       amebdm(ir,ic) = radfde(ir,lk,ic,isi)*radfrl(ir,lk,ic,isj)*vde/rx(ir)/dnom
                       amecdm(ir,ic) = radflo(ir,lk,ic,isi)*radfrl(ir,lk,ic,isj)*vde/rx(ir)/dnom2
                       amedam(ir,ic) = radfrl  (ir,lk,ic,isi)*radf(ir,lk,ic,isj)*vde/rx(ir)/dnom21
                       amedbm(ir,ic) = radfrl(ir,lk,ic,isi)*radfde(ir,lk,ic,isj)*vde/rx(ir)/dnom21
                       amedcm(ir,ic) = radfrl(ir,lk,ic,isi)*radflo(ir,lk,ic,isj)*vde/rx(ir)/dnom12
                       ameddm(ir,ic) = radfrl(ir,lk,ic,isi)*radfrl(ir,lk,ic,isj)*vde/rx(ir)/dnom1
                    enddo
                    call cali(ameddm(1,ic),totdd,rx(1),dx2,jri)
                    call cali(ameadm(1,ic),totad,rx(1),dx2,jri)
                    call cali(amebdm(1,ic),totbd,rx(1),dx2,jri)
                    call cali(amecdm(1,ic),totcd,rx(1),dx2,jri)
                    call cali(amedam(1,ic),totda,rx(1),dx2,jri)
                    call cali(amedbm(1,ic),totdb,rx(1),dx2,jri)
                    call cali(amedcm(1,ic),totdc,rx(1),dx2,jri)
                 enddo
                 ri_mat(4,4,l,ity,isi,isj) = coe*totdd
                 ri_mat(1,4,l,ity,isi,isj) = coe*totad
                 ri_mat(2,4,l,ity,isi,isj) = coe*totbd
                 ri_mat(3,4,l,ity,isi,isj) = coe*totcd
                 ri_mat(4,1,l,ity,isi,isj) = coe*totda 
                 ri_mat(4,2,l,ity,isi,isj) = coe*totdb 
                 ri_mat(4,3,l,ity,isi,isj) = coe*totdc 
              endif
           enddo
           !rschmid
           !   Calculate some integrals required for computation of 
           !   < \phi_sc | H_dirac | \phi_sc >.
           !rschmid
           do l=0,lomax
              if (isi.eq.isj.and.loorext(l,ity)) then
                 lk = l+1
                 !honza
                 hexl(l,ity,isi,1)=e(lk,ity,isi)
                 hexl(l,ity,isi,2)=1.d0
                 hexl(l,ity,isi,3)=pi2lor(l,ity,isi)*elor2(l,ity,isi)+hscalc(radf(1,lk,1,isi),radfrl(1,lk,1,isi),vrr(1,isi),elor2(l,ity,isi),l,jri,dx2,rx(1))
                 hexl(l,ity,isi,4)=0.d0
                 hexl(l,ity,isi,5)=0.d0
                 hexl(l,ity,isi,6)=pe2lor(l,ity,isi)*elor2(l,ity,isi)+hscalc(radfde(1,lk,1,isi),radfrl(1,lk,1,isi),vrr(1,isi),elor2(l,ity,isi),l,jri,dx2,rx(1))
                 hexl(l,ity,isi,7)=e(lk,ity,isi)*pi2lor(l,ity,isi)
                 hexl(l,ity,isi,8)=pi2lor(l,ity,isi)
                 hexl(l,ity,isi,9)=elor2(l,ity,isi)+hscalc(radfrl(1,lk,1,isi),radfrl(1,lk,1,isi),vrr(1,isi),elor2(l,ity,isi),l,jri,dx2,rx(1))
              end if
           enddo
           !rschmid
!!! Only the upper part of the spin-matrix has yet been calculated.
           !rschmid
        enddo  ! 91
     enddo     ! 91
  ENDDO    !109 CONTINUE
  RETURN
END subroutine garadme

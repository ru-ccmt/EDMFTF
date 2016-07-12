      subroutine garadorb(e,vru,p,dp,pe,dpe,ri_orb,   &
                         jspin,kpot,ipr)

        USE param
      USE loabc; USE loabcr; USE lolog; USE rlolog; USE rpars; USE orb
      USE struct
      USE peic
      implicit real*8(a-h,o-z)

! subroutine based on garadme, calculate radial elements of orb potential
! in this version vorb(r)=1 is assumed
! local orbitals included, but p1/2 probably is not correct
! P.Novak October 2001
      real*8 ri_orb(4,4,0:labc,nato,2,2)
      real*8 vrr(nrad,2),rx(nrad)

      real*8 totaa(2),totbb(2),totcc(2),totab(2),totba(2),       &
             totac(2),totca(2),totbc(2),totcb(2)
      real*8 ameaam(NRAD,2), amebbm(NRAD,2), ameabm(NRAD,2)
      real*8 ameadm(NRAD,2), amebdm(NRAD,2), amecdm(NRAD,2)
      real*8 amedam(NRAD,2), amedbm(NRAD,2), amedcm(NRAD,2)
      real*8 ameddm(NRAD,2), ameccm(NRAD,2), ameacm(NRAD,2)
      real*8 amebcm(NRAD,2), amebam(NRAD,2), amecam(NRAD,2)
      real*8 amecbm(NRAD,2)

      real*8 vru,e
      DIMENSION       VRU(NRAD,NATO,2),E(LMX,NATO,2)
     
      integer jspin,kpot,ipr

      real*8 p,dp,pe,dpe
      DIMENSION       P(Labc+1,NATO,2),DP(Labc+1,NATO,2),PE(Labc+1,NATO,2), &
                      DPE(Labc+1,NATO,2)

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

!rschmid
!    coc  : velocity of light in Rydberg units
!    cin  : inverse velocity of light squared. (Hartree units)
!rschmid 
      COE=1.D0
      CIN=1.d0/CLIGHT**2
      do 1 il=0,labc
       do 1 ia=1,nato
        do 1 i1=1,2
         do 1 i2=1,2
           do 335 j1=1,4
            do 335 j2=1,4
      ri_orb(j1,j2,il,ia,i1,i2) = 0.d0
 335   continue
 1     continue

!rschmid 
!  ity runs over all non-equivalent atoms.
!rschmid
      DO 109 index=1,natorb 
      ity=iat(index)
        DX2 = DX(ITY)
        JRI = JRJ(ITY)
        RMT2= RMT(ITY)

!rschmid
!  The logarithmic mesh is setup in RX for the respective atom.
!rschmid
      do i = 1, jri
        rx(i) = rmt2 * dexp(dfloat(i-jri)*dx2)
      enddo

      do isi=1,jspin
        call atpar(e,vru,labc+1,ity,radf,radfde,    &
                   radflo,radfrl,                &
                   p,dp,pe,dpe,pei,isi)
      enddo

      do 91 isi=1,jspin
      do 91 isj=1,jspin
!.........................................................................
!.....CALCULATE THE MATRIX ELEMENT
!
      L=ll(index)
      LK=L+1
      DO 104 IC=1,2
      DO 106 IR=1,JRI
      AMEAAM(IR,IC)=RADF(IR,LK,IC,isi)*RADF(IR,LK,IC,isj)
      AMEBBM(IR,IC)=RADfde(IR,LK,IC,isi)*RADFDE(IR,LK,IC,isj)
      AMEABM(IR,IC)=RADF(IR,LK,IC,isi)*RADFDE(IR,LK,IC,isj)
      AMEBAM(IR,IC)=RADFDE(IR,LK,IC,isi)*RADF(IR,LK,IC,isj)
  106 CONTINUE
!
      CALL CALI(AMEAAM(1,IC),TOTAA(IC),RX(1),DX,JRI)
      CALL CALI(AMEBBM(1,IC),TOTBB(IC),RX(1),DX,JRI)
      CALL CALI(AMEABM(1,IC),TOTAB(IC),RX(1),DX,JRI)
      CALL CALI(AMEBAM(1,IC),TOTBA(IC),RX(1),DX,JRI)
  104 CONTINUE
!

      ri_orb(1,1,L,ITY,isi,isj)=COE*(TOTAA(1)+CIN*TOTAA(2))
      ri_orb(2,2,L,ITY,isi,isj)=COE*(TOTBB(1)+CIN*TOTBB(2))
      ri_orb(1,2,L,ITY,isi,isj)=COE*(TOTAB(1)+CIN*TOTAB(2))
      ri_orb(2,1,L,ITY,isi,isj)=COE*(TOTBA(1)+CIN*TOTBA(2))
!
!.........................................................................
!.....CALCULATE THE MATRIX ELEMENT of LO
!
      if (loor(l,ity)) then
        DO 204 IC=1,2
        DO 206 IR=1,JRI
          ameccm(ir,ic)= radflo(ir,lk,ic,isi)*    &
                          radflo(ir,lk,ic,isj)
          ameacm(ir,ic)= radf  (ir,lk,ic,isi)*    &
                         radflo(ir,lk,ic,isj)
          amebcm(ir,ic)=radfde(ir,lk,ic,isi)*     &
                        radflo(ir,lk,ic,isj)
          amecam(ir,ic)= radflo  (ir,lk,ic,isi)*  &
                         radf(ir,lk,ic,isj)
          amecbm(ir,ic)=radflo(ir,lk,ic,isi)*     &
                        radfde(ir,lk,ic,isj)
  206   CONTINUE
        CALL CALI(AMECCM(1,IC),TOTCC(IC),RX(1),DX,JRI)
        CALL CALI(AMEACM(1,IC),TOTAC(IC),RX(1),DX,JRI)
        CALL CALI(AMEBCM(1,IC),TOTBC(IC),RX(1),DX,JRI)
        CALL CALI(AMECAM(1,IC),TOTCA(IC),RX(1),DX,JRI)
        CALL CALI(AMECBM(1,IC),TOTCB(IC),RX(1),DX,JRI)
  204   CONTINUE
!
        ri_orb(3,3,L,ITY,isi,isj)=COE*(TOTCC(1)+CIN*TOTCC(2))
        ri_orb(1,3,L,ITY,isi,isj)=COE*(TOTAC(1)+CIN*TOTAC(2))
        ri_orb(3,1,L,ITY,isi,isj)=COE*(TOTCA(1)+CIN*TOTCA(1))
        ri_orb(2,3,L,ITY,isi,isj)=COE*(TOTBC(1)+CIN*TOTBC(1))
        ri_orb(3,2,L,ITY,isi,isj)=COE*(TOTCB(1)+CIN*TOTCB(1))
      endif

91    continue                              
  109 CONTINUE

      RETURN
      END

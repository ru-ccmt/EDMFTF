      subroutine vderiv(lpot,ity,lspin,rx,vru,vrr,vder,jri)

        USE param
      implicit none

      real*8 vru(nrad,nato,2), vrr(nrad,2)
      real*8 vder(nrad,2)
      real*8 rx(nrad)
      real*8 vaver,der
      real*8 x(5),y(5) 

      integer lpot
      integer lspin
      integer isi,i,imt
      integer jri,ity

      do 30 isi=1,lspin
!.....use 5-order-langrange-interpretaion to calculate dV/dr.
        do i=1,jri
          if (lpot.eq.3) then
            vaver = ( vru(i,ity,1) + vru(i,ity,2) )/2.
            vrr(i,1) = vaver / rx(i)
            vrr(i,2) = vrr(i,1)
          else 
            vrr(i,isi) = vru(i,ity,isi) / rx(i)
          endif
        enddo

!.....VRR is the real potential in units of Ryd.
      IMT=JRI
      DO 30 I=1,IMT
      IF(I.EQ.1) THEN
        X(1)=RX(1)
        X(2)=RX(2)
        X(3)=RX(3)
        X(4)=RX(4)
        X(5)=RX(5)
        Y(1)=VRR(1,isi)
        Y(2)=VRR(2,isi)
        Y(3)=VRR(3,isi)
        Y(4)=VRR(4,isi)
        Y(5)=VRR(5,isi)
      ELSE IF(I.EQ.2) THEN
          X(1)=RX(2)
          X(2)=RX(1)
          X(3)=RX(3)
          X(4)=RX(4)
          X(5)=RX(5)
          Y(1)=VRR(2,isi)
          Y(2)=VRR(1,isi)
          Y(3)=VRR(3,isi)
          Y(4)=VRR(4,isi)
          Y(5)=VRR(5,isi)
        ELSE IF(I.EQ.IMT-1) THEN
          X(1)=RX(IMT-1)
          X(2)=RX(IMT)  
          X(3)=RX(IMT-2)
          X(4)=RX(IMT-3)
          X(5)=RX(IMT-4)
          Y(1)=VRR(IMT-1,isi)
          Y(2)=VRR(IMT,isi)  
          Y(3)=VRR(IMT-2,isi)
          Y(4)=VRR(IMT-3,isi)
          Y(5)=VRR(IMT-4,isi)
        ELSE IF(I.EQ.IMT) THEN
          X(1)=RX(IMT)
          X(2)=RX(IMT-1)
          X(3)=RX(IMT-2)
          X(4)=RX(IMT-3)
          X(5)=RX(IMT-4)
          Y(1)=VRR(IMT,isi)
          Y(2)=VRR(IMT-1,isi)
          Y(3)=VRR(IMT-2,isi)
          Y(4)=VRR(IMT-3,isi)
          Y(5)=VRR(IMT-4,isi)
        ELSE
          X(1)=RX(I)
          X(2)=RX(I-1)
          X(3)=RX(I-2)
          X(4)=RX(I+1)
          X(5)=RX(I+2)
          Y(1)=VRR(I,isi)
          Y(2)=VRR(I-1,isi)
          Y(3)=VRR(I-2,isi)
          Y(4)=VRR(I+1,isi)
          Y(5)=VRR(I+2,isi)
      END IF
      CALL LAGDER(5,X,Y,DER)
!  vder is dV/dr in units of Ryd/a.u.

      if(lpot.eq.3)then
        vder(i,1)=der  
        vder(i,2)=der  
      else
        VDER(I,isi)=DER
      endif
   30 CONTINUE

      return 
      end

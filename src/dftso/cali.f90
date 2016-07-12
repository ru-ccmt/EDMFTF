      SUBROUTINE CALI(AT,SUM,RX,DX,JWI)
        USE param
      IMPLICIT REAL*8(A-H,O-Z)
!---------------------------------------------------
!
      DIMENSION AT(NRAD),RX(NRAD)
!*****************************************************
      SUM=0.0D0
      SUM1=0.0D0
      JS=JWI-2*(JWI/2)+1
      JF=JWI-1
      DO 1 I=JS,JF,2
      SUM1=SUM1+(2.D0*RX(I)*AT(I)+RX(I+1)*AT(I+1))
    1 CONTINUE
      IF(JS.EQ.1) THEN
      SUM=(SUM1+SUM1-RX(JWI)*AT(JWI))*DX/3.D0
      ELSE IF(JS.EQ.2) THEN
      SUM=(SUM1+SUM1-RX(JWI)*AT(JWI)+RX(1)*AT(1))*DX/3.D0
      END IF
      RETURN
      END

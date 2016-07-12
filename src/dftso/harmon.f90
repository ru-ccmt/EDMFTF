      SUBROUTINE HARMON(YR,YI,N,X,Y,Z)
        USE param
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER J,JJ,L,M
      COMPLEX*16 YY(LMX2)
      DIMENSION YR(LMX*(LMX+1)/2,N),YI(LMX*(LMX+1)/2,N)
      DIMENSION X(1),Y(1),Z(1)
      DIMENSION A(3)
      DO 1 I=1,N
      A(1)=X(I)
      A(2)=Y(I)
      A(3)=Z(I)
      CALL YLM(A,LMAX,YY)
!
      JJ=0
      J=0
      DO 2 L=0,LMAX
        DO 2 M=-L,L
          J=J+1
          IF(M.GE.0) THEN
            JJ=JJ+1
             YR(JJ,I)=Dble(YY(J))
             YI(JJ,I)=aIMAG(YY(J))
          ENDIF
    2 CONTINUE
    1 CONTINUE
      RETURN
      END

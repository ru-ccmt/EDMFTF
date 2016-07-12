      SUBROUTINE LAGDER(N,X,Y,DER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),Y(N)
      XU=X(1)
      DE=0.D0
      DO 10 K=2,N
      DE=DE+1.D0/(XU-X(K))
   10 CONTINUE
      DE=DE*Y(1)
      DO 20 I=2,N
      ANU=1.D0
      DO 30 J=2,N
      IF(J.EQ.I) GO TO 30
      ANU=ANU*(XU-X(J))
   30 CONTINUE
      ADE=1.D0
      DO 40 J=1,N
      IF(J.EQ.I) GO TO 40
      ADE=ADE*(X(I)-X(J))
   40 CONTINUE
      DE=DE+Y(I)*ANU/ADE
   20 CONTINUE
      DER=DE
      RETURN
      END

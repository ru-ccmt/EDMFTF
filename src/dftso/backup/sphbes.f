      SUBROUTINE SPHBES(N,X,FJ)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FJ(*)
      DATA XLIM/1.D0/,HF/.5D0/,ZERO/0.D0/,ONE/1.D0/,TNHF/1.05D1/, &
      FFT/1.5D1/,T25/1.D18/,TN25/1.D-18/,TN50/1.D-36/,two/2.D0/
      IF(N.EQ.0) THEN
      IF(X.EQ.0.D0) THEN
      FJ(1)=1.D0
      RETURN
      END IF
      FJ(1)=DSIN(X)/X
      RETURN
      ELSE IF(N.EQ.1) THEN
      IF(X.EQ.0.D0) THEN
      FJ(1)=1.D0
      FJ(2)=0.D0
      RETURN
      END IF
      FJ(1)=DSIN(X)/X
      FJ(2)=DSIN(X)/X**2-DCOS(X)/X
      RETURN
      END IF
      IF (N.GE. 0) GO TO 7
    1 WRITE(6,2)
    2 FORMAT (33H1 ERROR, N SHOULD NOT BE NEGATIVE  )
      GO TO 99
    7 IF (X.GE.ZERO)  GO TO 10
    8 WRITE(6,9)
    9 FORMAT (33H1 ERROR, X SHOULD NOT BE NEGATIVE  )
      GO TO 99
   10 IF (X.GT.XLIM) GO TO 25
!
!     THE FOLLOWING SECTION IS ADDED TO HANDLE THE X<1.0 CASES
!     USING THE POLYNOMIAL EXPANSION OF SPHERICAL BESSEL FUNCTION.
!
!                                                  JUN YE    05/31/90
      DO 100 I=1,N+1
      BRU=SPHBRU(I-1,X)
      COB=1.D0
      DO 110 J=1,I-1
      COB=COB*X/(2.D0*J+1)
  110 CONTINUE
      FJ(I)=BRU*COB
  100 CONTINUE
      RETURN
!
   25 CUFAC=4.2
      IF (X.LT.(N-2)) CUFAC=TNHF/(N+HF-X)
      NS=N+5+X*CUFAC
!
!     ADD ADDITIONAL FACTOR
!
      NS=NS + (FFT/(ONE+DSQRT(X)))
!
      IF ((NS+2).LE.100) GO TO 30
      WRITE(6,28) X,NS
  28  FORMAT (1X,'* WARNING. FOR X=',G13.6,'   BESSH WANTS TO START AT N=',I5)
      NS=98
      IF (X.GT.FFT) GO TO 99
  30  CONTINUE
 313  FFO=ZERO
      FFN=TN25
      M=NS-1
      XI=ONE/X
      FM=(M+M)+ONE
      SDR=FM*TN50
 314  FFP=FM*XI*FFN-FFO
      IF (DABS(FFP).LT.T25) GO TO 315
      SDR=SDR*TN50
      FFP=FFP*TN25
      FFN=FFN*TN25
 315  SDR=SDR + (FM-TWO)*FFP*FFP
      FFO=FFN
      FFN=FFP
      IF (M.LE.N) GO TO 316
      M=M-1
      FM=FM-TWO
      GO TO 314
 316  FJ(M)=FFN
      FJ(M+1)=FFO
      GO TO 33
  32  FJ(M)=FM*XI*FJ(M+1)-FJ(M+2)
      IF(DABS(FJ(M)).GE.T25) GO TO 56
      SDR=SDR + (FM-TWO)*FJ(M)*FJ(M)
      IF (M.LE.1) GO TO 34
   33 M = M-1
      FM=FM-TWO
      GO TO 32
   34 SER=ONE/DSQRT(SDR)
   39 MM = N+1
      DO 40 M=1,MM
      FJ(M)=FJ(M)*SER
   40 CONTINUE
      GO TO 98
   56 JJ= M+1
      NS=N+1
      DO 57 J = JJ,NS
      FJ(J)=FJ(J)*TN25
   57 CONTINUE
      SDR=SDR*TN50
      GO TO 32
   99 CALL EXIT(1)
   98 continue
      return
      END

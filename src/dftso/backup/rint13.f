      SUBROUTINE RINT13(C1,C2,A,B,X,Y,S,JATOM,RNOT,DX,JRI)
!
        USE param
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
!        Arguments   
!
      INTEGER            JATOM
      DOUBLE PRECISION   C1, C2, S
      DOUBLE PRECISION   A(*), B(*), X(*), Y(*)
!
!     ..................................................................
!
!        perform radial integrals required  (D.D.Koelling)
!
!     ..................................................................
!       
!        Common blocks   
!
!
!        DH(i)   - step size of the logarithmic radial mesh for atom i
!        JRJ(i)  - number of radial (logarithmic) mesh points for atom i
!        VR(j,i) - spherical part (l=0, m=0) of the total potential r*V
!                  at mesh point j for atom i 
!
      INTEGER            JRJ(NATO)
      DOUBLE PRECISION   DH(NATO), VR(NRAD,NATO)
!       
!        Local Scalars   
!
      INTEGER            J, J1, JRI
      DOUBLE PRECISION   D, DX, GTEM1, GTEM2, P1, P2, R, R1, RNOT, Z2
      DOUBLE PRECISION   Z4
!       
!        Intrinsic Functions   
!
      INTRINSIC          EXP, MOD
!       
!      RNOT = RO(JATOM)
!      DX = DH(JATOM)
!      JRI = JRJ(JATOM)
      D = EXP(DX)
!
      J = 3 - MOD(JRI,2)
      J1 = J - 1
      R = RNOT*(D**(J-1))
      R1 = R/D
      Z4 = 0.0D+0
      Z2 = 0.0D+0
!
!        LOOP ... IF (J .GE. JRI) EXIT ... END LOOP
!
!   10 Z4=Z4+R*(C1*A(J)*X(J)+C2*B(J)*Y(J))
   10 CONTINUE
         GTEM1 = C1*A(J)*X(J)
         GTEM2 = C2*B(J)*Y(J)
         Z4 = Z4 + R*(GTEM1+GTEM2)
         IF (Z4 .EQ. 0.0D+0) Z4 = 2.0D-55
         R = R*D
         J = J + 1
         IF (J .GE. JRI) GOTO 20
         Z2 = Z2 + R*(C1*A(J)*X(J) + C2*B(J)*Y(J))
         R = R*D
         J = J + 1
!
!        END LOOP
!
      GOTO 10
   20 CONTINUE
      P1 = RNOT*(C1*A(1)*X(1) + C2*B(1)*Y(1))
      P2 = R1*(C1*A(J1)*X(J1) + C2*B(J1)*Y(J1))
      S = 2*Z2 + 4*Z4 + R*(C1*A(J)*X(J) + C2*B(J)*Y(J)) + P2
      S = (DX*S+P1)/3.0D+0
      IF (J1 .GT. 1) S = S + 0.5D+0*DX*(P1+P2)
!
      RETURN
!
!        End of 'RINT13'
!
      END

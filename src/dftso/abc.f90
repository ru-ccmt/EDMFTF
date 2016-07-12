SUBROUTINE ABC (JATOM,J,p,dp,pe,dpe,pei,rmt,isi,jlo,lapw)
  !
  USE param
  USE loabc
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  !
  !        Arguments
  !
  INTEGER            J,  JATOM, jlo
  logical lapw
  !n  isi=1,2 is the spin index down, up
  DOUBLE PRECISION   DP(Labc+1,NATO,2), DPE(Labc+1,NATO,2)
  DOUBLE PRECISION   P(Labc+1,NATO,2), PE(Labc+1,NATO,2), PEI(Labc+1,NATO,2)
  !        Local Parameters
  !
  DOUBLE PRECISION   CUTOFF
  PARAMETER          (CUTOFF = 200.0D+0)
  !
  !        Local Scalars
  !
  DOUBLE PRECISION   XAC, XBC
  !
  !        Intrinsic Functions
  !
  INTRINSIC          SQRT
  !
  if(lapw) then
     XAC = PLO(J,JATOM,isi)*DPE(J+1,JATOM,isi) - DPLO(J,JATOM,isi)*PE(J+1,JATOM,isi)
     XAC = XAC*RMT*RMT
     XBC = PLO(J,JATOM,isi)*DP(J+1,JATOM,isi) - DPLO(J,JATOM,isi)*P(J+1,JATOM,isi)
     XBC = -XBC*RMT*RMT
     CLO(J,jlo,JATOM,isi) = XAC*(XAC + 2.0D+0*PI12LO(J,JATOM,isi)) + XBC*(XBC*PEI(J+1,JATOM,isi) + 2.0D+0*PE12LO(J,JATOM,isi)) + 1.0D+0
     CLO(J,jlo,JATOM,isi) = 1.0D0/SQRT(CLO(J,jlo,JATOM,isi))
     CLO(J,jlo,JATOM,isi) = MIN(CLO(J,jlo,JATOM,isi),CUTOFF)
     ALO(J,jlo,JATOM,isi) = CLO(J,jlo,JATOM,isi)*XAC
     BLO(J,jlo,JATOM,isi) = CLO(J,jlo,JATOM,isi)*XBC
     !     WRITE (6,6000) J,ALO(J,JATOM,isi), BLO(J,JATOM,isi),
     !    &               CLO(J,JATOM,isi)
  else
     !.....APW definitions
     if(jlo.eq.1) then
        alonorm=sqrt(1.d0+(P(J+1,JATOM,isi)/PE(J+1,JATOM,isi))**2*PEI(J+1,JATOM,isi))
        ALO(J,jlo,JATOM,isi) = 1.d0 /alonorm 
        BLO(J,jlo,JATOM,isi) = -P(J+1,JATOM,isi)/PE(J+1,JATOM,isi)/alonorm
        CLO(J,jlo,JATOM,isi) = 0.d0
     else 
        xbc=-P(J+1,JATOM,isi)/PLO(J,JATOM,isi)
        xac=sqrt(1+xbc**2+2*xbc*PI12LO(J,JATOM,isi))
        ALO(J,jlo,JATOM,isi) = 1.d0/xac 
        BLO(J,jlo,JATOM,isi) = 0.d0
        CLO(J,jlo,JATOM,isi) = xbc/xac
     endif
  endif
  write (6,10) j,alo(j,jlo,jatom,isi),blo(j,jlo,jatom,isi),clo(j,jlo,jatom,isi)
10 FORMAT ('LO COEFFICIENT: l,A,B,C  ',i2,5X,3F12.5)
  
  
  RETURN
  !
  !
  !        End of 'ABC'
  !
END SUBROUTINE ABC

       subroutine abc_r (jatom,j,p,dp,pe,dpe,pei,rmt,isi,lapw)
         USE loabcr
         USE param
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      logical lapw
      DOUBLE PRECISION   DP(Labc+1,NATO,2),DPE(Labc+1,NATO,2)
      DOUBLE PRECISION   P(Labc+1,NATO,2),PE(Labc+1,NATO,2),PEI(Labc+1,NATO,2)
      DOUBLE PRECISION   CUTOFF
      PARAMETER          (CUTOFF = 200.0D+0)
      DOUBLE PRECISION   XAC, XBC
      INTRINSIC          SQRT

      if (lapw) then
      XAC =  PLOR(J,JATOM,isi)*DPE(J+1,JATOM,isi) - &
            DPLOR(J,JATOM,isi)* PE(J+1,JATOM,isi)
      XAC = XAC*RMT*RMT
      XBC =  PLOR(J,JATOM,isi)* DP(J+1,JATOM,isi) - &
            DPLOR(J,JATOM,isi)*  P(J+1,JATOM,isi)
      XBC = -XBC*RMT*RMT
      CLOR(J,JATOM,isi) = XAC*(XAC + 2.0D0*PI2LOR(J,JATOM,isi)) &
                        + XBC*(XBC*PEI(J+1,JATOM,isi) &
                        + 2.0D0*PE2LOR(J,JATOM,isi)) &
                        + 1.0D0
      CLOR(J,JATOM,isi) = 1.0D0/SQRT(CLOR(J,JATOM,isi))
      CLOR(J,JATOM,isi) = MIN(CLOR(J,JATOM,isi),CUTOFF)
      ALOR(J,JATOM,isi) = CLOR(J,JATOM,isi)*XAC
      BLOR(J,JATOM,isi) = CLOR(J,JATOM,isi)*XBC
      else
      xbc=-P(J+1,JATOM,isi)/PLOR(J,JATOM,isi)
      xac=sqrt(1+xbc**2+2*xbc*PI2LOR(J,JATOM,isi))
      ALOR(J,JATOM,isi) = 1.d0/xac
      BLOR(J,JATOM,isi) = 0.d0
      CLOR(J,JATOM,isi) = xbc/xac
      end if
      write (6,10)j,alor(j,jatom,isi),blor(j,jatom,isi), &
                  clor(j,jatom,isi)
      RETURN
!
 10   FORMAT ('RLO COEFFICIENT: l,A,B,C  ',i2,5X,3F12.5)
!
!        End of 'ABC_r'
!
      END

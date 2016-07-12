      subroutine atpar(e,vru,nt,jatom,radf,radfde, &
                         radflo,radfrl, &
                         p,dp,pe,dpe,pei,isi)
!
      USE param
      USE loabc; USE loabcr; USE lolog; USE rlolog; USE rpars
      USE struct
      IMPLICIT  REAL*8 (A-H,O-Z)
!      implicit none
!
!rschmid
!  Note that vru is  v*r in Rydberg unit on input.
!  It is transformed to v*r in Hartree units.
!rschmid
      real*8 vru
!
!        Arguments
!
      INTEGER            NT
      LOGICAL            REL
!
!     ..................................................................
!
!        ATPAR calculates the solutions u(l) of the radial Schroedinger
!        Equation, the energy derivatives ue(l), the radial non muffin-
!        tin integrals uvu, uevu, and the Gaunt-coefficients.
!        Spin-orbitless relativistic equations are used.
!
!     ..................................................................
!
!        Common blocks
!
!
!                            
!
!n     DIMENSION       RADFLO(NRAD,LOMAX,2,2)    error in original code???

      real*8  radf  (nrad,labc+1,2,2)
      real*8  radfde(nrad,labc+1,2,2)
      real*8  radflo(nrad,lomax+1,2,2)
      real*8  radfrl(nrad,lomax+1,2,2)


      real*8  a(nrad),b(nrad),ap(nrad),bp(nrad),ae(nrad),be(nrad)
      common  /work/     a,b,ap,bp,ae,be
      save    /work/

      dimension      vru(nrad,nato,2),vr(nrad),e(lmx,nato,2)
      dimension      p(labc+1,nato,2),dp(labc+1,nato,2),pe(labc+1,nato,2), &
                     dpe(labc+1,nato,2),pei(labc+1,nato,2)

      character*4 emain

!
!        Data statements
!
      DATA DELE /2.0D-3/
      DATA IMAG /(0.0D+0,1.0D+0)/
!
!        initialize constants
!        CFEIN1 and CFEIN2 are used in RINT13 to calculate the
!        integrals of radialfunctions and are derivates of the
!        Feinstruktur-constant (in Hartrees)
!
      REL=(.true.)
      IF (REL) THEN
         CFEIN1 = 1.0D+0
         CFEIN2 = 1.0D0/(CLIGHT)**2
      ELSE
         CFEIN1 = 1.0D+0
         CFEIN2 = 1.0D-22
      ENDIF
!
      DX2 = DX(JATOM)
      JRI = JRJ(JATOM)
      RMT2=RMT(JATOM)
      RNOT = RMT2*DEXP(DX2*(1.D0-JRI))

!rschmid
!  The potential*r in Hartree units is calculated.  
!rschmid
      WRITE(6,7) JATOM,isi
    7 FORMAT(/10X,'ATOMIC PARAMETERS FOR ',I2,'  SPIN=',I1/)
      WRITE(6,5)(E(ll,jatom,isi),ll=1,7)
      write(6,*)'for lo:'
      do nl=1,nloat
      write(6,6)(elo(ll,nl,jatom,isi),ll=0,3)
      end do
    6 FORMAT(10X,7F7.2)
    5 FORMAT(10X,' ENERGY PARAMETERS ARE',7F7.2)
      WRITE(6,14)
      do j=1,jri
        vr(j)=vru(j,jatom,isi)/2
      enddo
      delei = .25d0/dele
      DO 150 J = 1, NT
         L = J - 1
      mrf(l,jatom)=2
         FL = L
         EI = E(J,JATOM,isi)/2.D0
!        calculate energy-derivative by finite difference
!        DELE is the up and downward energy-shift in Hartrees
!
         E1 = EI - DELE
         CALL OUTWIN(REL,VR(1),RNOT,DX2,JRI,E1,FL,UVB,DUVB, &
                     NODEL,ZZ(JATOM))
         CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM,RNOT,DX2,JRI)
         TRX = 1.0D+0/SQRT(OVLP)
         DO 100 M = 1, JRI
            AE(M) = TRX*A(M)
            BE(M) = TRX*B(M)
  100    CONTINUE
         UVB = TRX*UVB
         DUVB = TRX*DUVB
         E1 = EI + DELE
         CALL OUTWIN(REL,VR(1),RNOT,DX2,JRI,E1,FL,UVE,DUVE, &
                     NODEU,ZZ(JATOM))
         CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM,RNOT,DX2,JRI)
         TRX = 1.0D+0/SQRT(OVLP)
         UVE = DELEI*(TRX*UVE - UVB)
         DUVE = DELEI*(TRX*DUVE - DUVB)
         DO 110 M = 1, JRI
            AE(M) = DELEI*(TRX*A(M) - AE(M))
            BE(M) = DELEI*(TRX*B(M) - BE(M))
  110    CONTINUE

!        calculate function at EI
!
         CALL OUTWIN(REL,VR(1),RNOT,DX2,JRI,EI,FL,UV,DUV,NODES &
                     ,ZZ(JATOM))
         CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM,RNOT,DX2,JRI)
         TRX = 1.0D+0/SQRT(OVLP)
         P(J,JATOM,isi) = TRX*UV
         DP(J,JATOM,isi) = TRX*DUV
         DO 120 M = 1, JRI
            A(M) = TRX*A(M)
            B(M) = TRX*B(M)
  120    CONTINUE
!
!        insure orthogonalization
!
         CALL RINT13(CFEIN1,CFEIN2,A,B,AE,BE,CROSS,JATOM,RNOT,DX2,JRI)
         TRY = -CROSS
!        IF (TRY .LT. -0.05D+0) WRITE (6,6070) L, TRY, OVLP
         DO 130 M = 1, JRI
            AE(M) = AE(M) + TRY*A(M)
            BE(M) = BE(M) + TRY*B(M)
  130    CONTINUE
         PE(J,JATOM,isi) = UVE + TRY*P(J,JATOM,isi)
         DPE(J,JATOM,isi) = DUVE + TRY*DP(J,JATOM,isi)
         CALL RINT13(CFEIN1,CFEIN2,AE,BE,AE,BE, &
                     PEI(J,JATOM,isi),JATOM,RNOT,DX2,JRI)
	 if (l.le.6) then
         write(6,8) l,p(j,jatom,isi),dp(j,jatom,isi),pe(j,jatom,isi), &
                 dpe(j,jatom,isi),pei(j,jatom,isi)
	 end if
         DO 140 M = 1, JRI
            RADF  (M,J,1,isi) = A(M)
            RADF  (M,J,2,isi) = B(M)
            RADFDE(M,J,1,isi) = AE(M)
            RADFDE(M,J,2,isi) = BE(M)
  140    CONTINUE
  150 CONTINUE
!         write(6,*)'Lo:'
!
!        and now for local orbitals
!
      DO 220 L = 0, LOMAX
         IF (LOOR(L,JATOM)) THEN
            J = L + 1
         mrf(l,jatom)=3
            FL = L
            EI = ELO(L,ilo(l,jatom),JATOM,isi)/2.D0

            CALL OUTWIN(REL,VR(1),RNOT,DX2,JRI,EI,FL,UV,DUV, &
                        NODES,ZZ(JATOM))
            CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM,RNOT,DX2,JRI)
            TRX = 1.0D+0/SQRT(OVLP)
            PLO(L,JATOM,isi) = TRX*UV
            DPLO(L,JATOM,isi) = TRX*DUV
            DO 180 M = 1, JRI
               A(M) = TRX*A(M)
               B(M) = TRX*B(M)
  180       CONTINUE
            CALL RINT13(CFEIN1,CFEIN2,RADF(1,J,1,isi),RADF(1,J,2,isi), &
                        A,B,PI12LO(L,JATOM,isi),JATOM,RNOT,DX2,JRI)
            CALL RINT13(CFEIN1,CFEIN2,RADFDE(1,J,1,isi), &
                         RADFDE(1,J,2,isi), &
                        A,B,PE12LO(L,JATOM,isi),JATOM,RNOT,DX2,JRI)
            DO 200 M = 1, JRI
               RADFLO(M,J,1,isi) = A(M)
               RADFLO(M,J,2,isi) = B(M)
  200       CONTINUE
         ENDIF
  220 CONTINUE
         DO l=0,lomax
            DO jlo=1,ilo(l,jatom)
              call abc (JATOM,L,p,dp,pe,dpe,pei,rmt(jatom),isi, &
                             jlo,lapw(l,jatom))
            ENDDO
         ENDDO
!       write(6,*)'for rlo:'
!rschmid
!    Now for the additional local orbital used in second
!    diagonalisation.
!rschmid
      DO L = 0, LOMAX
         IF (LOOREXT(L,JATOM)) THEN
            mrf(l,jatom)=4
            J = L + 1
            FL = L   
            kappa = l

!        calculate energy-derivative by finite difference
!        DELE is the up and downward energy-shift in hartrees
!

            de = extde(jatom,l)
            emain = 'CONT'
!rschmid
! Test for RLO energy in  valence band
!rschmid

            ei = extei (jatom, l)
            write(6,*) ei
            CALL SELECT(L,ei,DE,EMAIN,REL,VR(1),RNOT,DX2,JRI &
                           ,ZZ(JATOM))
            
            elor(l,jatom,isi) = ei
	    elo(l,nloat,jatom,isi)=ei
	    write(6,*)'relativistic lo:'
	    write(6,*)l,elor(l,jatom,isi)
            write(8,*)'RELATIVISTIC LOs:'
            write(8,*)'on atom ',jatom,'e(',l,')=',elor(l,jatom,isi)
            write(8,*)
            eparam = elor(l,jatom,isi)

            ei = ei/2.d0

            CALL diracout(REL,VR,RNOT,DX2,JRI,EI,L,kappa,UV,DUV, &
                        NODES,ZZ(JATOM))

!            do m = 1, jri
!                r_m = rnot * exp(dx2*(m-1))
!                write(70,*)r_m,A(m),B(m)
!            end do

            call dergl(a,b,rnot,dx2,jri)

            do m = 1, jri
                r_m = rnot * exp(dx2*(m-1))
                b(m) =  &
                     b(m)*r_m/ &
             ( 2.d0*clight+ &
               (eparam-2.d0*vr(m)/r_m)/(2.d0*clight)) 
!            write(71,*)r_m,b(m)                                        
!            b(m)=0.0d0
            enddo


            CALL RINT13(CFEIN1,1.d0,A,B,A,B,OVLP,JATOM,RNOT,DX2,JRI)

            TRX = 1.0D+0/SQRT(OVLP)
             PLOR(L,JATOM,isi) = TRX*UV
            DPLOR(L,JATOM,isi) = TRX*DUV

            DO M = 1, JRI
               A(M) = TRX*A(M)
               B(M) = TRX*B(M)
            enddo
           
            CALL RINT13(CFEIN1,1.d0/clight,RADF(1,J,1,isi), &
                        RADF(1,J,2,isi),a,b, &
                        PI2LOR(L,JATOM,isi),JATOM,RNOT,DX2,JRI)
            CALL RINT13(CFEIN1,1.d0/clight,RADFDE(1,J,1,isi), &
                        RADFDE(1,J,2,isi),a,b, &
                        PE2LOR(L,JATOM,isi),JATOM,RNOT,DX2,JRI)
            call rint13(cfein1,1.d0/clight,radflo(1,j,1,isi), &
                        radflo(1,j,2,isi),a,b,r22,jatom,rnot,dx2,jri)


            do m = 1, jri
               radfrl(m,j,1,isi) =a(m)
               radfrl(m,j,2,isi) =b(m)
 444       format(i4,2e15.6)
            enddo

 40      continue

         elor2(l,jatom,isi) = elor(l,jatom,isi)
         write(6,*)p(l+1,jatom,isi)
         call abc_r(jatom,L,p,dp,pe,dpe,pei,rmt(jatom),isi,lapw)
         else
           elor(l,jatom,isi) = 1.d+6
         ENDIF
      enddo

!
!        end of loop over atoms in unitcell
!
!
      RETURN
!
!
!
   14 FORMAT(/11X,1HL,5X,4HU(R),10X, &
             5HU'(R),9X,5HDU/DE,8X,6HDU'/DE,6X,7HNORM-U')
 5000 FORMAT (1X,I1,2F10.5,A4)
 5010 FORMAT (3X)
 5020 FORMAT (16X,I2,/,/)
 5030 FORMAT (16X,I2,5X,I2,/)
 5040 FORMAT (3X,4E19.12)
 5050 FORMAT (/,/,/)
 5060 FORMAT (/)
 5070 FORMAT (/,/)
 6000 FORMAT (/,10X,'ATOMIC SPHERE DEPENDENT PARAMETERS FOR ATOM  ',A10, &
              /,/,10X,'OVERALL ENERGY PARAMETER IS',F10.4)
 6010 FORMAT (10X,'E(',I2,2H)=,F10.4)
 6020 FORMAT (/,10X,'POTENTIAL PARAMETERS ',/,11X,'L',5X,'U(R)',10X, &
              'U''(R)',9X,'DU/DE',8X,'DU''/DE',6X,'NORM-U''')
 6030 FORMAT (10X,I2,5E14.6,5X,3I2)
 6040 FORMAT (/,/,/,' AT.NR.:',I2,5X,'POSITION: ',3F8.3,5X, &
              'MULTIPLICITY:',I5)
 6050 FORMAT (3X,6I3,5F15.8,I3)
 6060 FORMAT ('lo ',6I3,5F15.8,I3)
 6070 FORMAT (10X,'FOR L=',I2,' CORRECTION=',E14.5,' OVERLAP=',E14.5)
 6080 FORMAT (/,/,3X,'NOT EQUIV ATOM ',A10,'  LOCAL ROTATION MATRIX')
 6090 FORMAT (30X,3F10.5)
 6100 FORMAT (13X,'EQUIV ATOM ',I3,3X,'POSITION: ',3F8.3)
 6110 FORMAT (30X,3F10.5)
!
!        End of 'ATPAR'
!
    8 FORMAT(10X,I2,5E21.10,5X,3I2)

      END

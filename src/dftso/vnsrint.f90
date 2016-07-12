	subroutine vnsrint(jatom,isi,radf,radfde,radfrl)

          USE param
          USE struct
          USE vns
      implicit real*8(a-h,o-z)

      dimension vlm(nrad,5),a(nrad),b(nrad)
      real*8 radf  (nrad,labc+1,2,2)
      real*8 radfde(nrad,labc+1,2,2)
      real*8 radfrl(nrad,lomax+1,2,2)

      itape=21+isi
      REWIND itape
      READ (itape,5070)
	do 340 ity=1,jatom
       READ (itape,5010)
            READ (itape,5020)LMMAX
            index = 0
            DO 230 LM1=1,LMMAX
               READ (itape,5030)LL,MM
               if (iabs(ll).eq.2.and.ity.eq.jatom) then
               lvns(ity)=.true.
               index=index+1
               mvns(index,ity)=isign(mm,ll)
               READ(itape,5040)(VLM(J,index),J=1,JRJ(JATOM))
	       else
               READ(itape,5040)(dummy,J=1,JRJ(JATOM))
               endif
               nvns(jatom)=index
               READ (itape,5060)
  230       CONTINUE
            READ (itape,5050)
  340	    continue

 5070 FORMAT (/,/)
 5010 FORMAT (3X)
 5020 FORMAT (15X,I3,/,/)
 5030 FORMAT (15X,I3,5X,I2,/)
 5040 FORMAT (3X,4E19.12)
 5050 FORMAT (/,/,/)
 5060 FORMAT (/)
	
      DX2 = DX(JATOM)
      JRI = JRJ(JATOM)
      RMT2=RMT(JATOM)
      RNOT = RMT2*DEXP(DX2*(1.D0-JRI))
	DO index=1,nvns(jatom)
	do 200 m=1,jri
	a(m)=radf(m,2,1,isi)*vlm(m,index)
        b(m)=radf(m,2,2,isi)*vlm(m,index)
 200    continue
        CALL RINT13(1.d0,1.0/(clight**2),A,B,radf(1,2,1,isi), &
        radf(1,2,2,isi),vaa(jatom,index,isi),JATOM,RNOT,DX2,JRI)
        CALL RINT13(1.d0,1.0/(clight**2),A,B,radfde(1,2,1,isi), &
        radfde(1,2,2,isi),vab(jatom,index,isi),JATOM,RNOT,DX2,JRI)
        CALL RINT13(1.d0,1.0/clight,A,B,radfrl(1,2,1,isi), &
        radfrl(1,2,2,isi),vad(jatom,index,isi),JATOM,RNOT,DX2,JRI)

        do 210 m=1,jri
        a(m)=radfde(m,2,1,isi)*vlm(m,index)
        b(m)=radfde(m,2,2,isi)*vlm(m,index)
 210    continue
        CALL RINT13(1.d0,1.0/(clight**2),A,B,radfde(1,2,1,isi), &
        radfde(1,2,2,isi),vbb(jatom,index,isi),JATOM,RNOT,DX2,JRI)
        CALL RINT13(1.d0,1.0/clight,A,B,radfrl(1,2,1,isi), &
        radfrl(1,2,2,isi),vbd(jatom,index,isi),JATOM,RNOT,DX2,JRI)

        do 220 m=1,jri
        a(m)=radfrl(m,2,1,isi)*vlm(m,index)
        b(m)=radfrl(m,2,2,isi)*vlm(m,index)
 220    continue
        CALL RINT13(1.d0,1.d0,A,B,radfrl(1,2,1,isi), &
        radfrl(1,2,2,isi),vdd(jatom,index,isi),JATOM,RNOT,DX2,JRI)
        ENDDO
!write(70,*)jatom,nvns(jatom)
!write(70,*)(mvns(i,jatom),i=1,nvns(jatom))
!       do index=1,nvns(jatom)
!       write(70,*)vaa(jatom,index,isi),vab(jatom,index,isi),&
!        vad(jatom,index,isi),vbb(jatom,index,isi),vbd(jatom,index,isi), &
!        vdd(jatom,index,isi)
!end do
	END

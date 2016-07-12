	subroutine hns(ity,iat,isi,iei,ief,tot)
	USE param
        USE abcd
	USE vns
        implicit real*8 (a-h,o-z)

        COMPLEX*16 rsum,tot,imag,imag1



	tot=(0.d0,0.d0)
        imag=(0.d0,1.d0)
      do 200 index=1,nvns(ity)
	m=mvns(index,ity)
        mm=iabs(mvns(index,ity))
 100    continue
	do 201 indi=2,4
         indf=indi+mm
          if ((indf).ge.2.and.(indf).le.4) then
	   rsum=(vaa(ity,index,isi)*conjg(abcdlm(1,indf,iat,ief,isi))* &
            abcdlm(1,indi,iat,iei,isi)+  &
            vbb(ity,index,isi)*conjg(abcdlm(2,indf,iat,ief,isi))*  &
            abcdlm(2,indi,iat,iei,isi)+  &
            vdd(ity,index,isi)*conjg(abcdlm(4,indf,iat,ief,isi))*  &
            abcdlm(4,indi,iat,iei,isi)+  &
            vab(ity,index,isi)*conjg(abcdlm(1,indf,iat,ief,isi))*  &
            abcdlm(2,indi,iat,iei,isi)+  &
            vab(ity,index,isi)*conjg(abcdlm(2,indf,iat,ief,isi))*  &
            abcdlm(1,indi,iat,iei,isi)+  &
            vad(ity,index,isi)*conjg(abcdlm(1,indf,iat,ief,isi))*  &
            abcdlm(4,indi,iat,iei,isi)+  &
            vad(ity,index,isi)*conjg(abcdlm(4,indf,iat,ief,isi))*  &
            abcdlm(1,indi,iat,iei,isi)+  &
            vbd(ity,index,isi)*conjg(abcdlm(2,indf,iat,ief,isi))*  &
            abcdlm(4,indi,iat,iei,isi)+  &
            vbd(ity,index,isi)*conjg(abcdlm(4,indf,iat,ief,isi))*  &
            abcdlm(2,indi,iat,iei,isi))  &
            *gaunt1(1,2,1,indf-3,mm,indi-3)
      
             IF (MM .NE. 0) THEN
                MINU = 1
                IMAG1 = (1.0D+0,0.0D+0)
                  IF (M.LT.0) THEN
                         IMAG1=(-1.0D+0,0.0D+0)
                         MINU=-1
                  ENDIF
                  IF (MOD(MM,2).EQ.1) THEN
                         IMAG1=-IMAG1
                         MINU=-MINU
                  ENDIF
                IF (MM.GT.0) MINU = 1
                rsum=rsum*imag1*dble(minu)/sqrt(2.d0)
             END IF
           tot=tot+rsum
          end if
 201    continue
        if (mm.gt.0) then
	    mm=-mm
	    goto 100
        end if
 200   continue
	end

	

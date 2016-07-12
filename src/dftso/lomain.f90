      subroutine lomain(a,bkx,bky,bkz,fac,             &
                        ne,lfirst,indj,jatom,n,isi)
        USE param
        USE abcd
      USE loabc; USE loabcr; USE lolog; USE rlolog; USE rpars
      USE struct
      USE rotmat
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16     A(NMAT,NUME,2),YL(LABC2),ALM,BLM,CLM
      COMPLEX*16     PHS,PHSHEL,CFAC,IMAG,CZERO,PT
      COMMON/GENER/  BR1(3,3),BR2(3,3)

      DIMENSION      PHS(NUME)
      DIMENSION      BK(3),BKROT(3),BKRLOC(3)
      DIMENSION      BKX(NMAT),BKY(NMAT),BKZ(NMAT)
      DATA           CZERO/(0.0D0,0.0D0)/,IMAG/(0.0D0,1.0D0)/
      DATA           ZERO/0.0D0/
      DATA           TOL1/1.0D-6/
!                                                                       
!------------------------------------------------------------------     
!                                                                       
!.initiales a,b,c of lo                                      
!	
      PI=ACOS(-1.0D0)                                                   
      TWOPI=2.D0*PI
      i=n-(nlo(jatom)+nlon(jatom)) 
      DO 10 L=0,LOMAX
        DO 15 jlo=1,ilo(l,jatom)
        CFAC=FAC*(IMAG**L)
        do 20 jneq=1,mult(jatom)
          DO 25 M1=-l,+l                                                    
          i=i+1                                   
          BK(1)=BKX(I)                                                      
          BK(2)=BKY(I)                                                      
          BK(3)=BKZ(I)     
        CALL ROTATE (BK,ROTIJ(1,1,indj),BKROT)                                                 
          BK(1)=BKROT(1)*BR1(1,1)+BKROT(2)*BR1(1,2)+BKROT(3)*BR1(1,3)   
          BK(2)=BKROT(1)*BR1(2,1)+BKROT(2)*BR1(2,2)+BKROT(3)*BR1(2,3)   
          BK(3)=BKROT(1)*BR1(3,1)+BKROT(2)*BR1(3,2)+BKROT(3)*BR1(3,3) 
        CALL ROTATE (BK,ROTLOC(1,1,jatom),BKRLOC)
        CALL YLM(BKRLOC,LABC,YL)
!
        ARG1=BKROT(1)*POS(1,LFIRST)*TWOPI
        ARG2=BKROT(2)*POS(2,LFIRST)*TWOPI
        ARG3=BKROT(3)*POS(3,LFIRST)*TWOPI
        ARGT=(BKX(I)*TAUIJ(1,indj)+BKY(I)*TAUIJ(2,indj)+  &
              BKZ(I)*TAUIJ(3,indj))*TWOPI

        PHSHEL=EXP(IMAG*(ARG1+ARG2+ARG3+ARGT))
          DO 50 NUM=1,NE
          PHS(NUM)=ZERO
          PT=PHSHEL*A(I,NUM,isi)
          if(dabs(dreal(pt)).gt.tol1)phs(num)=phs(num)+dreal(pt)
          if(dabs(dimag(pt)).gt.tol1)phs(num)=phs(num)+dimag(pt)*imag
  50      continue                                                        
          do 30 m=-l,+l                                                    
            index=l*(l+1)+m+1 
            do 40 num=1,ne
      abcdlm(1,index,indj,num,isi)=abcdlm(1,index,indj,num,isi)  &
         +alo(l,jlo,jatom,isi)*conjg(yl(index))*phs(num)*cfac
      abcdlm(2,index,indj,num,isi)=abcdlm(2,index,indj,num,isi)  &
         +blo(l,jlo,jatom,isi)*conjg(yl(index))*phs(num)*cfac
      abcdlm(3,index,indj,num,isi)=abcdlm(3,INDEX,indj,num,isi)  &
         +clo(l,jlo,jatom,isi)*conjg(yl(index))*phs(num)*cfac
  40        CONTINUE    
  30      CONTINUE
  25    CONTINUE                                                          
  20    CONTINUE                                                          
  15    CONTINUE
  10  CONTINUE     
      return                   
      END        

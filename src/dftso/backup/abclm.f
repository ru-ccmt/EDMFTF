SUBROUTINE ABCLM(meigve,SS,NE,NV,KV,P,DP,PE,DPE,isi)   
  USE param
  USE loabc; USE loabcr; USE lolog; USE rlolog; USE rpars
  USE abcd
  USE struct
  USE vns
  USE rotmat
  IMPLICIT REAL*8(A-H,O-Z)
  !
  !-----X CALCULATE RADIAL ENERGY DERIVATIVE ANALYTICALLY.
  !     THIS PROGRAM CALCULATE THE WAVEFUNCTION INSIDE MT SPHERE
  !     (TRANSFERED FROM OUTSIDE MT). THE WAVEFUNCTION INSIDE MT IS
  !     EXPRESSED IN TERMS OF Alm Blm Clm FOR EACH ATOM.
  !
  !**********************************************************************
  complex*16     meigve(nmat,nume,2)
  complex*16,allocatable ::  yl(:),h_yl(:,:),h_alyl(:,:),h_blyl(:,:),alm(:,:),blm(:,:)
  dimension      al(hblock),bl(hblock)
  complex*16     phs,phshel,cfac,imag,czero,pt
  
  COMMON/GENER/  BR1(3,3),BR2(3,3)
  
  DIMENSION      PHS(NUME),FJ(0:LMAX,NMAT),DFJ(0:LMAX,NMAT)
  DIMENSION      BK(3),BKROT(3),bkrloc(3)
  DIMENSION      FCT(100),BKX(NMAT),BKY(NMAT),BKZ(NMAT)
  DIMENSION      P(Labc+1,NATO,2),DP(Labc+1,NATO,2),PE(Labc+1,NATO,2),DPE(Labc+1,NATO,2)
  DIMENSION      SS(3),KV(3,NMAT,2)
  
  DATA           CZERO/(0.0D0,0.0D0)/,IMAG/(0.0D0,1.0D0)/
  DATA           ZERO/0.0D0/
  DATA           TOL1/1.0D-6/
 
  allocate     ( yl(labc2),h_yl(labc2,hblock),h_alyl(labc2,hblock),h_blyl(labc2,hblock),alm(labc2,nume),blm(labc2,nume) )
  PI=ACOS(-1.0D0)                                                   
  TWOPI=2.D0*PI
  
  !!    FCT(2*j+1) = j!
  FCT(:)=0.d0
  Y=1.0D0                                                           
  DO i=1,50
     j=2*i-1                                                           
     FCT(j)=y
     y=y*i                                                             
  ENDDO
  
  DO I=1,NV
     BKX(I)=(SS(1)+KV(1,I,isi))                                               
     BKY(I)=(SS(2)+KV(2,I,isi))                                               
     BKZ(I)=(SS(3)+KV(3,I,isi))         
  enddo
  
  INDJ=0
  call cputim(dt0)
  DO JA=1,NAT                                                   ! 777
     lfirst=indj+1
     FAC=4.0D0*PI*RMT(JA)**2/SQRT(VOL)
     CALL HARMON2(NV,BKX,BKY,BKZ,LMAX,FJ,DFJ,RMT(JA))
     DO MU=1,MULT(JA)                                           ! 777
        INDJ=INDJ+1
     
        do ix=1,labc2                                           
           do num=1,ne+nnrlo 
              alm(ix,num)=czero
              blm(ix,num)=czero
              abcdlm(1,ix,indj,num,isi)=czero
              abcdlm(2,ix,indj,num,isi)=czero
              abcdlm(3,ix,indj,num,isi)=czero
              abcdlm(4,ix,indj,num,isi)=czero
           enddo
        enddo

        DO II=1,NV-(nlo(1)+nlon(1)+nlov(1)),hblock                ! 120
           i3=0
           DO I=ii,min(ii+hblock-1,NV-(nlo(1)+nlon(1)+nlov(1)))   ! 121
              i3=i3+1
              BK(1)=BKX(I)                                                      
              BK(2)=BKY(I)                                                      
              BK(3)=BKZ(I) 
              CALL ROTATE (BK,ROTIJ(1,1,indj),BKROT)
              BK(1)=BKROT(1)*BR1(1,1)+BKROT(2)*BR1(1,2)+BKROT(3)*BR1(1,3)   
              BK(2)=BKROT(1)*BR1(2,1)+BKROT(2)*BR1(2,2)+BKROT(3)*BR1(2,3)   
              BK(3)=BKROT(1)*BR1(3,1)+BKROT(2)*BR1(3,2)+BKROT(3)*BR1(3,3)  
              CALL ROTATE (BK,ROTLOC(1,1,JA),BKRLOC)
              CALL YLM(BKRLOC,LABC,YL)
              
              ARG1=BKROT(1)*POS(1,LFIRST)*TWOPI
              ARG2=BKROT(2)*POS(2,LFIRST)*TWOPI
              ARG3=BKROT(3)*POS(3,LFIRST)*TWOPI
              ARGT=(BKX(I)*TAUIJ(1,indj)+BKY(I)*TAUIJ(2,indj)+BKZ(I)*TAUIJ(3,indj))*TWOPI
              PHSHEL=EXP(IMAG*(ARG1+ARG2+ARG3+ARGT))
              do index=1,LABC2
                 h_yl(index,i3)=conjg(yl(index))*phshel
              end do
              !121           CONTINUE
           ENDDO
           index=0
           do L=0,LABC
              i3=0
              do i=ii,min(ii+hblock-1,NV-(nlo(1)+nlon(1)+nlov(1)))
                 i3=i3+1
                 IF(lapw(l,ja)) THEN
                    al(i3)=dfj(l,i)*pe(l+1,ja,isi)- fj(l,i)*dpe(l+1,ja,isi)
                    bl(i3)= fj(l,i)*dp(l+1,ja,isi)-dfj(l,i)*  p(l+1,ja,isi)
                 ELSE
                    al(i3)= fj(l,i)/p(l+1,ja,isi)/rmt(ja)**2
                    bl(i3)= 0.d0
                 ENDIF
              end do

              do m=1,2*l+1
                 index=index+1
                 i3=0
                 do i=ii,min(ii+hblock-1,NV-(nlo(1)+nlon(1)+nlov(1)))
                    i3=i3+1
                    h_alyl(index,i3)=AL(i3)*h_YL(INDEX,i3)
                    h_blyl(index,i3)=BL(i3)*h_YL(INDEX,i3)
                 enddo
              enddo
           enddo

           ibb=min(hblock,NV-(nlo(1)+nlon(1)+nlov(1))-ii+1)
           lda=labc2
           ldb=nmat
           
           call zgemm('N','N',index,ne,ibb,(1.d0,0.d0),h_alyl,lda,meigve(ii,1,isi),ldb,(1.d0,0.d0),alm,lda)
           call zgemm('N','N',index,ne,ibb,(1.d0,0.d0),h_blyl,lda,meigve(ii,1,isi),ldb,(1.d0,0.d0),blm,lda)
        ENDDO !120        CONTINUE
        index=0
        do l=0,LABC
           cfac=fac*(imag**l)
           do m=1,2*l+1
              index=index+1
              do num=1,ne
                 abcdlm(1,index,indj,num,isi)=alm(index,num)*cfac
                 abcdlm(2,index,indj,num,isi)=blm(index,num)*cfac
              enddo
           enddo
        enddo
        if (nlo(ja).ne.0) call lomain(meigve,bkx,bky,bkz,fac,ne,lfirst,indj,JA,nv,isi)
        if (nrlo(ja).ne.0) call rlomain(bkx,bky,bkz,fac,ne,lfirst,indj,JA,nv,isi,kv)
     ENDDO
  ENDDO
  !777 CONTINUE
  deallocate   ( yl,h_yl, h_alyl,h_blyl, alm,blm)
  RETURN
END SUBROUTINE ABCLM

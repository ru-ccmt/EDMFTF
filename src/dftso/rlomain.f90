subroutine rlomain(bkx,bky,bkz,fac,ne,lfirst,indj,jatom,n,isi,kv)
  USE param
  USE abcd
  USE loabcr; USE lolog; USE rlolog
  USE struct
  USE abcd
  USE rotmat
  implicit none
  complex*16     yl(labc2)
  complex*16     phs,phshel,cfac,imag,czero,pt
  integer ii,kv(3,nmat,2),iat,lfirst
  real*8  br1,br2,argt
  common/gener/  br1(3,3),br2(3,3)
  dimension      phs(nume)

  real*8 bk,bkrot,bkx,bky,bkz,bkrloc
  dimension      bk(3),bkrot(3),bkrloc(3)
  dimension      bkx(nmat),bky(nmat),bkz(nmat)
  integer jatom
  real*8 zero,tol1
  data           czero/(0.0d0,0.0d0)/,imag/(0.0d0,1.0d0)/
  data           zero/0.0d0/
  data           tol1/1.0d-6/
  !                                                                       
  real*8 arg1,arg2,arg3
  integer i,l,m,m1,jneq,num,n,ne
  real*8 pi,twopi,fac
  integer indj,isi
  integer index
  !                                                                       
  !                                                                       
  !.initiales a,b,c of rlo                                      
  !
  pi=acos(-1.0d0)                                                   
  twopi=2.d0*pi
  ii=0
  do iat=1,jatom-1
     do l=0,lomax
        if (loorext(l,iat)) then
           ii=ii+(2*l+1)*mult(iat)
        end if
     end do
  end do
  
  i=n-(nlo(jatom)+nlon(jatom))
  do l=0,lomax                                                       ! 10
     if (.not.loorext(l,jatom)) then
        i=i+(2*l+1)*mult(jatom)*ilo(l,jatom)
        CYCLE !goto 10
     endif
     cfac=fac*(imag**l)
     do jneq=1,mult(jatom)                                             ! 20
        do m1=-l,+l                                                    ! 25
           i=i+1 
           ii=ii+1 
           bk(1)=bkx(i)
           bk(2)=bky(i)
           bk(3)=bkz(i)
	   kv(1,n+ii,isi)=kv(1,i,isi)
           kv(2,n+ii,isi)=kv(2,i,isi)
           kv(3,n+ii,isi)=kv(3,i,isi)
           CALL ROTATE (BK,ROTIJ(1,1,indj),BKROT)
           bk(1)=bkrot(1)*br1(1,1)+bkrot(2)*br1(1,2)+bkrot(3)*br1(1,3)   
           bk(2)=bkrot(1)*br1(2,1)+bkrot(2)*br1(2,2)+bkrot(3)*br1(2,3)   
           bk(3)=bkrot(1)*br1(3,1)+bkrot(2)*br1(3,2)+bkrot(3)*br1(3,3)  
           BK(1)=BKROT(1)*BR1(1,1)+BKROT(2)*BR1(1,2)+BKROT(3)*BR1(1,3)
           BK(2)=BKROT(1)*BR1(2,1)+BKROT(2)*BR1(2,2)+BKROT(3)*BR1(2,3)
           BK(3)=BKROT(1)*BR1(3,1)+BKROT(2)*BR1(3,2)+BKROT(3)*BR1(3,3)
           CALL ROTATE (BK,ROTLOC(1,1,jatom),BKRLOC)
           CALL YLM(BKRLOC,LABC,YL)
           !
           ARG1=BKROT(1)*POS(1,LFIRST)*TWOPI
           ARG2=BKROT(2)*POS(2,LFIRST)*TWOPI
           ARG3=BKROT(3)*POS(3,LFIRST)*TWOPI
           ARGT=(BKX(I)*TAUIJ(1,indj)+BKY(I)*TAUIJ(2,indj)+BKZ(I)*TAUIJ(3,indj))*TWOPI
           
           PHSHEL=EXP(IMAG*(ARG1+ARG2+ARG3+ARGT))
           num = ne + nrlov(jatom) + (m1 + l + 1) + (jneq-1) * (2*l+1)
           phs(num) = zero
           pt = phshel*1.d0
           
           if (dabs(dreal(pt)).gt.tol1)  phs(num) = phs(num) + dreal(pt)
            if (dabs(dimag(pt)).gt.tol1) phs(num) = phs(num) + dimag(pt)*imag
            do m = -l,l
               index = l*(l+1)+m+1 
               abcdlm(1,index,indj,num,isi) = abcdlm(1,index,indj,num,isi) + alor(l,jatom,isi)*conjg(yl(index))*phs(num)*cfac
               abcdlm(2,index,indj,num,isi) = abcdlm(2,index,indj,num,isi) + blor(l,jatom,isi)*conjg(yl(index))*phs(num)*cfac
               abcdlm(4,index,indj,num,isi) = abcdlm(4,index,indj,num,isi) + clor(l,jatom,isi)*conjg(yl(index))*phs(num)*cfac
            enddo
         enddo ! 25          continue
      enddo !20       continue
   enddo !10          continue     
   return                   
 END subroutine rlomain

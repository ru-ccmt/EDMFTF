subroutine kptout(ss,bname,weight,kkk,kv,jspin,nv,ne)
  USE param
  USE lolog; USE rlolog
  USE hmsout
  USE struct
  implicit real*8 (a-h,o-z)
  character*10    bname
  real*8 weight
  real*8 ss(3)
  integer kv(3,nmat,2)
  integer kt(3,nmat)
  integer nv(2),ne(2)
  complex*16 vecdummy
  real*8 totnorm,ene
  complex*16 vt(nmat)

  WRITE(6,534) SS(1),SS(2),SS(3),BNAME,ne(1)+ne(2)+2*nnrlo,WEIGHT
  if(kkk.eq.1) then
     write(8,*)
     write(8,*) '       SPIN-ORBIT EIGENVALUES:'
     WRITE(8,534) SS(1),SS(2),SS(3),BNAME,ne(1)+ne(2)+2*nnrlo,WEIGHT
  end if
  
  do isi=1,2
     mfile=40+isi
     is=isi
     if (jspin.eq.1) is = 1
     !...rearrange
     do n=1,nv(is)-(nlo(1)+nlon(1)+nlov(1))
	kt(1,n)=kv(1,n,is)
        kt(2,n)=kv(2,n,is)
        kt(3,n)=kv(3,n,is)
     end do
     n=n-1
     nn=n
     nl=0
     do ity=1,nat
	do l=0,lomax
           do jlo=1,ilo(l,ity)
              do index=1,(2*l+1)*mult(ity)
                 n=n+1
                 nn=nn+1
                 kt(1,n)=kv(1,nn,is)
                 kt(2,n)=kv(2,nn,is)
                 kt(3,n)=kv(3,nn,is)
              enddo
           enddo
           if (loorext(l,ity)) then
              do index=1,(2*l+1)*mult(ity)
                 n=n+1
                 nl=nl+1
                 kt(1,n)=kv(1,nv(is)+nl,is)
                 kt(2,n)=kv(2,nv(is)+nl,is)
                 kt(3,n)=kv(3,nv(is)+nl,is)
              enddo
           end if
	end do
     end do
     !....rearrange
     write(mfile)ss(1),ss(2),ss(3),bname,nv(is)+nnrlo,neig,weight
     write(10+mfile,'(3e19.12,a10,2i6,f5.1)')ss(1),ss(2),ss(3),bname,nv(is)+nnrlo,neig,weight
     write(mfile)(kt(1,i),kt(2,i),kt(3,i),i=1,nv(is)+nnrlo)
     !       do i=1,nv(is)+nnrlo
     !       write(10+mfile,'(3i5)')kt(1,i),kt(2,i),kt(3,i)
     !       end do
  enddo
      
  ! write also on dummy vector file 43 to be used in LAPW2 (see comment).
  nvd=1
  WRITE(53,'(3e19.12,a10,2i6,f5.1)') SS(1),SS(2),SS(3),BNAME,nvd,neig,WEIGHT
      
  !rschmid
  !  Write output to case.outputso
  !  and case.scfso
  !rschmid
  write(6,531) (en(ie),ie=1,neig)
  write(6,6010) 0
  write(6,6030)  
  if(kkk.eq.1) then
     write(8,530) (en(ie),ie=1,neig)
     write(8,6030)
  endif
  
  !rschmid
  !  The expansion coefficients are written to the vector files
  !  coresponding to fort.41 and fort.42.
  !rschmid
  do isi=1,2
     ms = 40+isi
     do j=1,neig 
        write(ms)j,en(j)
        write(ms+10,*)j,en(j)
        !...rearrange
        do n=1,nv(is)-(nlo(1)+nlon(1)+nlov(1))
           vt(n)=vect(n,j,isi)
        end do
	n=n-1
        nn=n
        nl=0
        do ity=1,nat
           do l=0,lomax
              do jlo=1,ilo(l,ity)
                 do index=1,(2*l+1)*mult(ity)
                    n=n+1
                    nn=nn+1
                    vt(n)=vect(nn,j,isi)
                 enddo
              enddo
              if (loorext(l,ity)) then
                 do index=1,(2*l+1)*mult(ity)
                    n=n+1
                    nl=nl+1
                    vt(n)=vect(nv(is)+nl,j,isi)
                 enddo
              end if

           end do
        end do
        !...rearrange
        write(ms)(vt(i),i=1,nv(isi)+nnrlo)
        !             do ii=1,nv(isi)+nnrlo
        !               write(ms+10,5)ii,vt(ii)
        !             end do
     enddo
  enddo
5 format(I4,2x,2f14.8)
  
  !rschmid
  ! The norms respective norms of the eigenstates are written
  !  to unit 45 and unit 46.
  !rschmid
  ! write norms of spin down part of vectors on file 45
  ! write norms of spin up part of vectors on file 46  
  write(45,202) (vnorm(j,1),j=1,neig)
  write(46,202) (vnorm(j,2),j=1,neig)
  
  ene = 0.5 + en(neig)  
  do j=1,neig
     !        totnorm = vnorm(j,1) + vnorm(j,2)
     !        write(47,201) j,en(j),totnorm,(vnorm(j,is),is=1,2)
     write(53,*) j,ene
     ene=ene+0.001
  enddo
  RETURN
  
  
201 format(i4,4e16.8)
202 format(4e20.12)
502 FORMAT(20I4)
530 FORMAT(8(7X,5F13.7/))
531 FORMAT(8(2X,5F13.7/))
533 FORMAT((3X,4E19.12))
534 FORMAT(5x,'K=',3f10.5,3x,a10,/,6x,'MATRIX SIZE=',i5,3x,'WEIGHT=',f5.2/,5x,'EIGENVALUES ARE:')
1818 FORMAT(12f10.5)
6010 FORMAT(I13,' EIGENVALUES BELOW THE ENERGY ',F10.5)
6030 FORMAT(7X,14('****'))
end subroutine kptout

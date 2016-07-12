subroutine hsocalc(h,ri_mat,ne,nnrlo,jspin)

  USE abcd
  USE couples
  USE param
  USE lolog
  USE struct
  implicit real*8 (a-h,o-z)

  complex*16      h(nume2,nume2),hso
  complex*16      czero,imag,dd
  real*8          ri_mat(4,4,0:labc,nato,2,2)   
!  complex*16      t(ndif*(labc+1)**2*4,2,nume2)   
  complex*16,allocatable ::  t(:,:,:)   
  dimension ne(2) 

  allocate(t(ndif*(labc+1)**2*4,2,nume2))   
  czero=(0.0d0,0.0d0)
  nban2=ne(1)+ne(2)+nnrlo*2
  n_scr=ne(1)+ne(2)
  nmatr=nban2*(nban2+1)/2

  do ibi = 1, nban2
     if (ibi.le.n_scr) then
        if (ibi.gt.ne(1)) then
           iei=ibi-ne(1)
           isi=2
        else
           iei=ibi
           isi=1
        end if
     else
        if ((ibi-n_scr).gt.nnrlo) then
           iei=ibi-ne(1)-nnrlo
           isi=2
        else
           iei=ibi-ne(2)
           isi=1
        end if
     end if

     if (jspin.eq.1) then
        isd=1
     else
        isd=isi
     endif
     
     do isf=1,2
        if (jspin.eq.1) then
           isu=1
        else
           isu=isf
        endif
        lfirst=1
        index=0
        DO ITY = 1,nat
           lfirst=lfirst+mult(ity-1)
           if (.not.lso(ity)) cycle
           DO INA = 1,mult(ITY)
              INDJ=lfirst+ina-1
              DO L = 1,LABC
                 DO MF = -L,L
                    ILMF = L*L +L +MF +1
                    DO l_mat2=1,mrf(l,ity)
                       index=index+1
                       t(index,isf,ibi)=czero
                       
                       do mi=MF-1,MF+1
                          ILMI = L*L +L +MI +1
                          IF(ABS(mi).le.L) THEN
                             dd=czero
                             do l_mat1=1,mrf(l,ity)
                                dd=dd+abcdlm(l_mat1,ilmi,indj,iei,isd)*   &
                                     ri_mat(l_mat1,l_mat2,l,ity,isd,isu)
                             end do
                             t(index,isf,ibi)=t(index,isf,ibi)+dd* &
                                  couplo(indj,l,mf,mi,isf,isi)
                          END IF
                       end do
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  do ibi = 1, nban2
     do ibf = 1, ibi
        hsoc  = czero
        ovrlc = czero
        if (ibi.le.n_scr) then
           if (ibi.gt.ne(1)) then
              iei=ibi-ne(1)
              isi=2
           else
              iei=ibi
              isi=1
           end if
        else
           if ((ibi-n_scr).gt.nnrlo) then
              iei=ibi-ne(1)-nnrlo
              isi=2
           else
              iei=ibi-ne(2)
              isi=1
           end if
        end if

        if (ibf.le.n_scr) then
           if (ibf.gt.ne(1)) then
              ief=ibf-ne(1)
              isf=2
           else
              ief=ibf
              isf=1
           end if
        else
           if ((ibf-n_scr).gt.nnrlo) then
              ief=ibf-ne(1)-nnrlo
              isf=2
           else
              ief=ibf-ne(2)
              isf=1
           end if
        end if
        
        if (jspin.eq.1) then
           isd=1
           isu=1
        else
           isd=isi
           isu=isf
        endif
        
        lfirst=1
        index=0
        hso=czero
        DO ITY = 1,nat
           lfirst=lfirst+mult(ity-1)
           if (.not.lso(ity)) cycle
           DO INA = 1,mult(ITY)
              INDJ=lfirst+ina-1
              DO L = 1,LABC
                 DO MF = -L,L
                    ILMF = L*L +L +MF +1
                    DO l_mat2=1,mrf(l,ity)
                       index=index+1
                       hso=hso+dconjg(abcdlm(l_mat2,ilmf,indj,ief,isu))* &
                            t(index,isf,ibi)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        h(ibf,ibi)=h(ibf,ibi)+hso
        h(ibi,ibf)=dconjg(h(ibf,ibi))
     ENDDO
  ENDDO
  deallocate (t) 
END subroutine hsocalc

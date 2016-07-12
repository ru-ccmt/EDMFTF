subroutine hsocalc(h,ri_mat,ne,nnrlo,jspin)
  USE abcd
  USE couples
  USE param
  USE lolog
  USE struct
  IMPLICIT NONE
  !implicit real*8 (a-h,o-z)
  complex*16, intent(inout) :: h(nume2,nume2)
  real*8,     intent(in)    :: ri_mat(4,4,0:labc,nato,2,2)
  integer,    intent(in)    :: ne(2), nnrlo, jspin
  ! locals
  !                         t(ndif*(labc+1)**2*4,2,nume2)   
  complex*16,allocatable ::  t(:,:,:)   
  complex*16 :: dd
  complex*16 :: hso
  integer    :: ibi, ibf, iei, ief, ilmf, ilmi, mu, latom, index, isd, isi, isu, isf, jatom, l, l_mat1, l_mat2, lfirst, mf, mi, n_scr, nban2, nmatr
  
  allocate(t(ndif*(labc+1)**2*4,2,nume2))
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
        DO jatom = 1,nat
           lfirst=lfirst+mult(jatom-1)
           if (.not.lso(jatom)) cycle
           DO mu = 1,mult(jatom)
              latom=lfirst+mu-1
              DO l = 1,LABC
                 DO MF = -L,L
                    ILMF = L*L+L+MF+1
                    DO l_mat2=1,mrf(l,jatom)  ! this is (1,2,3,4) for (alm,blm,clm_1,clm_2)
                       index=index+1
                       t(index,isf,ibi)=0.d0 
                       do mi=MF-1,MF+1
                          ILMI = L*L+L+MI+1
                          IF(ABS(mi).le.L) THEN
                             dd=0.d0
                             do l_mat1=1,mrf(l,jatom)
                                dd=dd+abcdlm(l_mat1,ilmi,latom,iei,isd)*ri_mat(l_mat1,l_mat2,l,jatom,isd,isu)
                             end do
                             ! t(jatom,mu,l,mf,kappa2,ibi,sf) = sum_{mi,kappa1} cuplo(latom,l,mf,mi,sf,si) * alm^{kappa1}_{l,mi,iei,si} <u^{kappa1}_{jatom,l,si}|u^{kappa2}_{jatom,l,si}>
                             t(index,isf,ibi)=t(index,isf,ibi)+dd*couplo(latom,l,mf,mi,isf,isi)
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
        hso=0.d0
        DO jatom= 1,nat
           lfirst=lfirst+mult(jatom-1)
           if (.not.lso(jatom)) cycle
           DO mu = 1,mult(jatom)
              LATOM=lfirst+mu-1
              DO L = 1,LABC
                 DO MF = -L,L
                    ILMF = L*L +L +MF +1
                    DO l_mat2=1,mrf(l,jatom)
                       index=index+1
                       ! h(ief,sf,iei,si) = \sum_{jatom,mu,l,mf,kappa2} alm^{* kappa2}_{l,mf,ibi,sf} * t(jatom,mu,l,mf,kappa2,ibi,sf)
                       ! or
                       ! h(ief,sf,iei,si) = sum_{latom,l,mf,kappa2,mi,kappa1} alm^{* kappa2}_{latom,l,mf,ibi,sf}*cuplo(latom,l,mf,mi,sf,si)*alm^{kappa1}_{latom,l,mi,iei,si} <u^{kappa1}_{jatom,l,si}|u^{kappa2}_{jatom,l,si}>
                       ! or
                       ! h(ief,sf,iei,si) = <psi_{ief,sf}|Y_{lmf}chi_{sf}cuplot(l,mf,mi,sf,si)Y^*_{lmi}chi_{si}|psi_{iei,si}>
                       hso=hso+dconjg(abcdlm(l_mat2,ilmf,latom,ief,isu))*t(index,isf,ibi)
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

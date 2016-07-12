SUBROUTINE GetSpinRotation(Rspin,rotij,crotloc,norbitals,natom,natm,iso,lmaxp,iatom,nl,cix,ll,iorbital)
  USE com_mpi, ONLY: myrank, master, Qprint
  USE structure, ONLY: BR1, rot_spin_quantization
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: Rspin(2,2,norbitals)
  REAL*8, intent(in)  :: rotij(3,3,natm), crotloc(3,3,natom)
  INTEGER, intent(in) :: norbitals, natom, natm, iso, lmaxp
  INTEGER, intent(in) :: iatom(natom), nl(natom), cix(natom,4), ll(natom,4), iorbital(natom,lmaxp+1)
  ! locals
  INTEGER :: icase, latom, lcase, icix, l, nind, iorb, i, j
  REAL*8  :: Det, phi1, the1, psi1, phi2, the2, psi2, rotij_proper(3,3), crotloc_proper(3,3)
  COMPLEX*16 :: R2_n(2,2), R2_o(2,2), Rtmp(2,2)
  REAL*8 :: R3tmp(3,3), tmp3(3,3), Trans3(3,3), BR1inv(3,3), rotij_cartesian(3,3), Id(3,3)!, x
  ! external
  REAL*8 :: detx
  INTEGER            STDERR
  PARAMETER          (STDERR = 99)
  !
  call INVERSSYMDEF(BR1,BR1inv)
  Id(:,:)=0.d0
  Id(1,1)=1.d0
  Id(2,2)=1.d0
  Id(3,3)=1.d0
  DO icase=1,natom
     latom = iatom(icase)   ! The succesive number of atom (all atoms counted)
     do lcase=1,nl(icase)
        icix = cix(icase,lcase)
        if ( icix.EQ.0 ) CYCLE
        l = ll(icase,lcase)
        nind = (2*l+1)*iso
        iorb = iorbital(icase,lcase)
        phi2=0.d0
        the2=0.d0
        psi2=0.d0
        
        if (.True.) then
           !!  local_axis_defined_by_locrot  <- local_axis_of_equivalent_atom <- from_spin_quantization_to_global_cartesian
           !!* Trans3 = crotloc(:,:,icase) * rotij_cartesian * rot_spin_quantization
           tmp3 = matmul(BR1,rotij(:,:,latom))
           rotij_cartesian = matmul(tmp3,BR1inv)
           tmp3 = matmul(crotloc(:,:,icase),rotij_cartesian)
           Trans3 = matmul(tmp3, rot_spin_quantization)
           !
           !WRITE(6,*) 'rotij_cartesian'
           !do i=1,3
           !   do j=1,3
           !      write(6,'(g20.10)',advance='no') rotij_cartesian(i,j)
           !   enddo
           !   write(6,*)
           !enddo
           !WRITE(6,*) 'rot_spin_quantization'
           !do i=1,3
           !   do j=1,3
           !      write(6,'(g20.10)',advance='no') rot_spin_quantization(i,j)
           !   enddo
           !   write(6,*)
           !enddo
           !WRITE(6,*) 'crotloc'
           !do i=1,3
           !   do j=1,3
           !      write(6,'(g20.10)',advance='no') crotloc(i,j,icase)
           !   enddo
           !   write(6,*)
           !enddo
           !WRITE(6,*) 'Trans3'
           !do i=1,3
           !   do j=1,3
           !      write(6,'(g20.10)',advance='no') Trans3(i,j)
           !   enddo
           !   write(6,*)
           !enddo
           !
           Det = detx(Trans3)
           Trans3 = transpose(Trans3*Det)
           if (sum(abs(Trans3-Id)).LT.1e-10) Trans3(:,:)=Id(:,:)
           CALL Angles(phi1,the1,psi1, Trans3 )
           CALL Spin_Rotation(Rtmp,phi1,the1,psi1)
        else  ! This is the old version, which is wrong, as rotij was not transformed to cartesian system, moreover, case.inso file was not taken into account
           Det = detx(rotij(:,:,latom))
           rotij_proper = transpose( rotij(:,:,latom)*Det )
           CALL Angles(phi1,the1,psi1, rotij_proper )
           CALL Spin_Rotation(R2_n,phi1,the1,psi1)
           CALL Space_Rotation(R3tmp, phi1,the1,psi1 )
           R3tmp = R3tmp-rotij_proper
           if (abs(R3tmp(1,1))+abs(R3tmp(1,2))+abs(R3tmp(1,3))+abs(R3tmp(2,1))+abs(R3tmp(2,2))+abs(R3tmp(2,3))+abs(R3tmp(3,1))+abs(R3tmp(3,2))+abs(R3tmp(3,3)) > 1e-6) then
              WRITE(STDERR,*) 'ERROR: difference does not vanish', R3tmp
           endif
           Det = detx(crotloc(:,:,icase))
           crotloc_proper = transpose( crotloc(:,:,icase)*Det )
           CALL Angles(phi2,the2,psi2, crotloc_proper )
           CALL Spin_Rotation(R2_o,phi2,the2,psi2)
           Rtmp = matmul(R2_n,R2_o)
        endif
        
        Rspin(:,:,iorb) = Rtmp
        
        if (Qprint) then
           !WRITE(6,*) 'crotloc_proper='
           !WRITE(6,*) crotloc_proper(:,:)
           WRITE(6,'(A,I3)') 'Spin-rotations for iorb=', iorb
           WRITE(6,*) 'angles=', phi1, the1, psi1!, phi2, the2, psi2
           WRITE(6,*) 'Rspin='
           do i=1,2
              write(6,'(2f12.6,2x,2f12.6)') (Rspin(i,j,iorb),j=1,2)
           enddo
        endif
     enddo
  ENDDO
END SUBROUTINE GetSpinRotation

REAL*8 FUNCTION acos_(x)
  IMPLICIT NONE
  REAL*8, intent(in) :: x
  if (abs(x).lt.1.d0) then
     acos_ = acos(x)
     return
  endif
  acos_ = 0.d0
  if (x.ge.1.d0) return
  if (x.le.-1.d0) acos_ = acos(-1.d0)
  return
END FUNCTION acos_


SUBROUTINE Angles(phi,the,psi,U)
  ! Routine gets rotation matrix U -- (only proper rotations are allowed)
  ! It returns three Euler angles (phi,the,psi), which generate rotation U by
  ! U = Rz(psi)*Rx(the)*Rz(phi)
  IMPLICIT NONE
  REAL*8, intent(out):: phi,the,psi
  REAL*8, intent(in) :: U(3,3)
  !
  REAL*8 :: acos_
  REAL*8 :: PI
  PARAMETER (PI = 3.141592653589793D0)

  the = acos_(U(3,3))
  
  if (abs(abs(U(3,3))-1)>1e-4) then
     phi = atan2(U(3,1),U(3,2))
     psi = PI-atan2(U(1,3),U(2,3))
  else
     psi=0.
     if (U(2,1)==0 .and. U(1,2)==0) then
        phi = acos_(U(1,1))-cos(the)*psi
     else
        phi = atan2( U(2,1)*cos(the), U(2,2)*cos(the) )-psi*cos(the)
     endif
  endif

  !WRITE(6,*) 'Routine Angles called for rotation='
  !do i=1,3
  !   write(6,'(2f15.9,2x,2f15.9,2x,2f15.9)') (U(i,j),j=1,3)
  !enddo
  !WRITE(6,*) 'angles are', phi,the,psi
  
END SUBROUTINE Angles

!!!!!! Rotations of real space vectors !!!!!!!!!!!!!!!
SUBROUTINE Rz_3(R, psi)
  IMPLICIT NONE
  REAL*8, intent(out) :: R(3,3)
  REAL*8, intent(in)  :: psi
  !
  R(1,1)=cos(psi); R(1,2)=-sin(psi); R(1,3)=0
  R(2,1)=sin(psi); R(2,2)=cos(psi); R(2,3)=0
  R(3,1)=0; R(3,2)=0; R(3,3)=1
END SUBROUTINE Rz_3

SUBROUTINE Ry_3(R,the)
  IMPLICIT NONE
  REAL*8, intent(out) :: R(3,3)
  REAL*8, intent(in)  :: the
  !
  R(1,1)=cos(the); R(1,2)=0.0; R(1,3)=sin(the)
  R(2,1)=0.0; R(2,2)=1.0; R(2,3)=0.0
  R(3,1)=-sin(the); R(3,2)=0.0; R(3,3)=cos(the)
END SUBROUTINE Ry_3

SUBROUTINE Rx_3(R, phi)
  IMPLICIT NONE
  REAL*8, intent(out) :: R(3,3)
  REAL*8, intent(in)  :: phi
  !
  R(1,1)=1.0; R(1,2)=0.0; R(1,3)=0.0
  R(2,1)=0.0; R(2,2)=cos(phi); R(2,3)=-sin(phi)
  R(3,1)=0.0; R(3,2)=sin(phi); R(3,3)=cos(phi)
END SUBROUTINE Rx_3

!!!!!! Rotations of spinors !!!!!!!!!!!!!!!
SUBROUTINE Rz_2(R, psi)
  IMPLICIT NONE
  COMPLEX*16, intent(out):: R(2,2)
  REAL*8, intent(in)     :: psi
  ! locals
  COMPLEX*16 :: i
  i = (0,1)
  R(1,1)=exp(psi/2.*i); R(1,2)=0.0;
  R(2,1)=0.0;  R(2,2)=exp(-psi/2.*i)
END SUBROUTINE Rz_2

SUBROUTINE Ry_2(R, the)
  IMPLICIT NONE
  COMPLEX*16, intent(out):: R(2,2)
  REAL*8, intent(in)     :: the
  ! 
  R(1,1)=cos(the/2.); R(1,2)=sin(the/2.)
  R(2,1)=-sin(the/2.); R(2,2)=cos(the/2.)
END SUBROUTINE Ry_2
  
SUBROUTINE Rx_2(R, phi)
  IMPLICIT NONE
  COMPLEX*16, intent(out):: R(2,2)
  REAL*8, intent(in)     :: phi
  ! locals
  COMPLEX*16 :: i
  i = (0,1)
  R(1,1)=cos(phi/2.);  R(1,2)=sin(phi/2.)*i
  R(2,1)=sin(phi/2.)*i; R(2,2)=cos(phi/2.)
END SUBROUTINE Rx_2

SUBROUTINE Space_Rotation( R, phi,the,psi )
  IMPLICIT NONE
  REAL*8, intent(in)  :: phi, the, psi
  REAL*8, intent(out) :: R(3,3)
  ! locals
  REAL*8 :: R1(3,3), R2(3,3), R3(3,3)
  CALL Rz_3(R1, psi)
  CALL Rx_3(R2, the)
  CALL Rz_3(R3, phi)
  !R = matmul(matmul(R1,R2),R3)
  R = matmul(R1,R2)
  R = matmul(R,R3)
END SUBROUTINE Space_Rotation

SUBROUTINE Spin_Rotation( R, phi,the,psi )
  IMPLICIT NONE
  REAL*8, intent(in)     :: phi, the, psi
  COMPLEX*16, intent(out):: R(2,2)
  ! locals
  COMPLEX*16 :: R1(2,2), R2(2,2), R3(2,2)
  CALL Rz_2(R1, psi)
  CALL Rx_2(R2, the)
  CALL Rz_2(R3, phi)
  R = matmul(matmul(R1,R2),R3)
END SUBROUTINE Spin_Rotation

REAL*8 FUNCTION detx(a)
  IMPLICIT NONE
  REAL*8, intent(in) :: a(3,3)
  !
  REAL*8 :: cc
  cc = a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(3,1)*a(2,2)*a(1,3)-a(1,1)*a(3,2)*a(2,3)-a(2,1)*a(1,2)*a(3,3)
  detx = cc
  return
END FUNCTION detx

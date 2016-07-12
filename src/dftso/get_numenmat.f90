SUBROUTINE get_numenmat
  CHARACTER*10  :: kname
  DO
     READ(54,'(3e19.12,a10,2i6)',IOSTAT=ios) S,T,Z,KNAME,N,NE
     IF (ios /= 0) EXIT
     nmat=MAX(n,nmat)
     nume=MAX(nen,nume)
     DO ii=1,ne
        READ(30,*)
     ENDDO
  ENDDO

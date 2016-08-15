SUBROUTINE ERRFLG(FNAME,MSG)
  use mpi, ONLY: stop_MPI, myrank, master
  IMPLICIT NONE
  CHARACTER*(*)      FNAME, MSG
  if (myrank.eq.master) then
     OPEN (99,FILE=FNAME,ERR=900)
     WRITE (99,9000) MSG
     CLOSE (99)
     OPEN (99,FILE=FNAME,ERR=900)
  else
     OPEN (99,FILE="."//trim(FNAME),ERR=900)
  endif
  RETURN
  !        Errors
900 write(*,*)'Cannot open error-file'
  call stop_MPI
  STOP 'ERRFLG - couldn''t open errorflag-file.'
  !
9000 FORMAT (A)
END SUBROUTINE ERRFLG

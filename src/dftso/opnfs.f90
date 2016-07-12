      SUBROUTINE OPNFS(DEFFN,ERRFN,FL)
        USE param
      IMPLICIT REAL*8(A-H,O-Z)
!
      CHARACTER*80    DEFFN,ERRFN,FNAME
      CHARACTER*11    STATUS,FORM
      LOGICAL         FL
!**********************************************************************
!
      FL=.FALSE.
      CALL GTFNAM(DEFFN,ERRFN)
      OPEN(1,FILE=DEFFN,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
!
!        check for in1c file for complex case
!
!
 11   OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=8002)
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING SPINORB.DEF !!!!'
      STOP 'SPINOB.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
      '  FORM:',FORM
      STOP 'OPEN FAILED'
 8001 CONTINUE
      CLOSE (1)
!
      END  

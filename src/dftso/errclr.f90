SUBROUTINE ERRCLR(FNAME)
  CHARACTER*(*)      FNAME
  close(99, status='delete')
  RETURN
END SUBROUTINE ERRCLR

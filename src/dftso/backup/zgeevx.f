      SUBROUTINE ZGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL, &
                         LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, &
                         RCONDV, WORK, LWORK, RWORK, INFO )
!
!  -- LAPACK driver routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          BALANC, JOBVL, JOBVR, SENSE
      INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
      DOUBLE PRECISION   ABNRM
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   RCONDE( * ), RCONDV( * ), RWORK( * ), &
                         SCALE( * )
      COMPLEX*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), &
                         W( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZGEEVX computes for an N-by-N complex nonsymmetric matrix A, the
!  eigenvalues and, optionally, the left and/or right eigenvectors.
!
!  Optionally also, it computes a balancing transformation to improve
!  the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
!  SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues
!  (RCONDE), and reciprocal condition numbers for the right
!  eigenvectors (RCONDV).
!
!  The right eigenvector v(j) of A satisfies
!                   A * v(j) = lambda(j) * v(j)
!  where lambda(j) is its eigenvalue.
!  The left eigenvector u(j) of A satisfies
!                u(j)**H * A = lambda(j) * u(j)**H
!  where u(j)**H denotes the conjugate transpose of u(j).
!
!  The computed eigenvectors are normalized to have Euclidean norm
!  equal to 1 and largest component real.
!
!  Balancing a matrix means permuting the rows and columns to make it
!  more nearly upper triangular, and applying a diagonal similarity
!  transformation D * A * D**(-1), where D is a diagonal matrix, to
!  make its rows and columns closer in norm and the condition numbers
!  of its eigenvalues and eigenvectors smaller.  The computed
!  reciprocal condition numbers correspond to the balanced matrix.
!  Permuting rows and columns will not change the condition numbers
!  (in exact arithmetic) but diagonal scaling will.  For further
!  explanation of balancing, see section 4.10.2 of the LAPACK
!  Users' Guide.
!
!  Arguments
!  =========
!
!  BALANC  (input) CHARACTER*1
!          Indicates how the input matrix should be diagonally scaled
!          and/or permuted to improve the conditioning of its
!          eigenvalues.
!          = 'N': Do not diagonally scale or permute;
!          = 'P': Perform permutations to make the matrix more nearly
!                 upper triangular. Do not diagonally scale;
!          = 'S': Diagonally scale the matrix, ie. replace A by
!                 D*A*D**(-1), where D is a diagonal matrix chosen
!                 to make the rows and columns of A more equal in
!                 norm. Do not permute;
!          = 'B': Both diagonally scale and permute A.
!
!          Computed reciprocal condition numbers will be for the matrix
!          after balancing and/or permuting. Permuting does not change
!          condition numbers (in exact arithmetic), but balancing does.
!
!  JOBVL   (input) CHARACTER*1
!          = 'N': left eigenvectors of A are not computed;
!          = 'V': left eigenvectors of A are computed.
!          If SENSE = 'E' or 'B', JOBVL must = 'V'.
!
!  JOBVR   (input) CHARACTER*1
!          = 'N': right eigenvectors of A are not computed;
!          = 'V': right eigenvectors of A are computed.
!          If SENSE = 'E' or 'B', JOBVR must = 'V'.
!
!  SENSE   (input) CHARACTER*1
!          Determines which reciprocal condition numbers are computed.
!          = 'N': None are computed;
!          = 'E': Computed for eigenvalues only;
!          = 'V': Computed for right eigenvectors only;
!          = 'B': Computed for eigenvalues and right eigenvectors.
!
!          If SENSE = 'E' or 'B', both left and right eigenvectors
!          must also be computed (JOBVL = 'V' and JOBVR = 'V').
!
!  N       (input) INTEGER
!          The order of the matrix A. N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the N-by-N matrix A.
!          On exit, A has been overwritten.  If JOBVL = 'V' or
!          JOBVR = 'V', A contains the Schur form of the balanced
!          version of the matrix A.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  W       (output) COMPLEX*16 array, dimension (N)
!          W contains the computed eigenvalues.
!
!  VL      (output) COMPLEX*16 array, dimension (LDVL,N)
!          If JOBVL = 'V', the left eigenvectors u(j) are stored one
!          after another in the columns of VL, in the same order
!          as their eigenvalues.
!          If JOBVL = 'N', VL is not referenced.
!          u(j) = VL(:,j), the j-th column of VL.
!
!  LDVL    (input) INTEGER
!          The leading dimension of the array VL.  LDVL >= 1; if
!          JOBVL = 'V', LDVL >= N.
!
!  VR      (output) COMPLEX*16 array, dimension (LDVR,N)
!          If JOBVR = 'V', the right eigenvectors v(j) are stored one
!          after another in the columns of VR, in the same order
!          as their eigenvalues.
!          If JOBVR = 'N', VR is not referenced.
!          v(j) = VR(:,j), the j-th column of VR.
!
!  LDVR    (input) INTEGER
!          The leading dimension of the array VR.  LDVR >= 1; if
!          JOBVR = 'V', LDVR >= N.
!
!  ILO,IHI (output) INTEGER
!          ILO and IHI are integer values determined when A was
!          balanced.  The balanced A(i,j) = 0 if I > J and
!          J = 1,...,ILO-1 or I = IHI+1,...,N.
!
!  SCALE   (output) DOUBLE PRECISION array, dimension (N)
!          Details of the permutations and scaling factors applied
!          when balancing A.  If P(j) is the index of the row and column
!          interchanged with row and column j, and D(j) is the scaling
!          factor applied to row and column j, then
!          SCALE(J) = P(J),    for J = 1,...,ILO-1
!                   = D(J),    for J = ILO,...,IHI
!                   = P(J)     for J = IHI+1,...,N.
!          The order in which the interchanges are made is N to IHI+1,
!          then 1 to ILO-1.
!
!  ABNRM   (output) DOUBLE PRECISION
!          The one-norm of the balanced matrix (the maximum
!          of the sum of absolute values of elements of any column).
!
!  RCONDE  (output) DOUBLE PRECISION array, dimension (N)
!          RCONDE(j) is the reciprocal condition number of the j-th
!          eigenvalue.
!
!  RCONDV  (output) DOUBLE PRECISION array, dimension (N)
!          RCONDV(j) is the reciprocal condition number of the j-th
!          right eigenvector.
!
!  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  If SENSE = 'N' or 'E',
!          LWORK >= max(1,2*N), and if SENSE = 'V' or 'B',
!          LWORK >= N*N+2*N.
!          For good performance, LWORK must generally be larger.
!
!  RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = i, the QR algorithm failed to compute all the
!                eigenvalues, and no eigenvectors or condition numbers
!                have been computed; elements 1:ILO-1 and i+1:N of W
!                contain eigenvalues which have converged.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            SCALEA, WANTVL, WANTVR, WNTSNB, WNTSNE, WNTSNN, &
                         WNTSNV
      CHARACTER          JOB, SIDE
      INTEGER            HSWORK, I, ICOND, IERR, ITAU, IWRK, K, MAXB, &
                         MAXWRK, MINWRK, NOUT
      DOUBLE PRECISION   ANRM, BIGNUM, CSCALE, EPS, SCL, SMLNUM
      COMPLEX*16         TMP
!     ..
!     .. Local Arrays ..
      LOGICAL            SELECT( 1 )
      DOUBLE PRECISION   DUM( 1 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLABAD, DLASCL, XERBLA, ZDSCAL, ZGEBAK, ZGEBAL, &
                         ZGEHRD, ZHSEQR, ZLACPY, ZLASCL, ZSCAL, ZTREVC, &
                         ZTRSNA, ZUNGHR
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX, ILAENV
      DOUBLE PRECISION   DLAMCH, DZNRM2, ZLANGE
      EXTERNAL           LSAME, IDAMAX, ILAENV, DLAMCH, DZNRM2, ZLANGE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, DCONJG, DIMAG, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      WANTVL = LSAME( JOBVL, 'V' )
      WANTVR = LSAME( JOBVR, 'V' )
      WNTSNN = LSAME( SENSE, 'N' )
      WNTSNE = LSAME( SENSE, 'E' )
      WNTSNV = LSAME( SENSE, 'V' )
      WNTSNB = LSAME( SENSE, 'B' )
      IF( .NOT.( LSAME( BALANC, 'N' ) .OR. LSAME( BALANC, &
          'S' ) .OR. LSAME( BALANC, 'P' ) .OR. LSAME( BALANC, 'B' ) ) ) &
           THEN
         INFO = -1
      ELSE IF( ( .NOT.WANTVL ) .AND. ( .NOT.LSAME( JOBVL, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( ( .NOT.WANTVR ) .AND. ( .NOT.LSAME( JOBVR, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( .NOT.( WNTSNN .OR. WNTSNE .OR. WNTSNB .OR. WNTSNV ) .OR. &
               ( ( WNTSNE .OR. WNTSNB ) .AND. .NOT.( WANTVL .AND. &
               WANTVR ) ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDVL.LT.1 .OR. ( WANTVL .AND. LDVL.LT.N ) ) THEN
         INFO = -10
      ELSE IF( LDVR.LT.1 .OR. ( WANTVR .AND. LDVR.LT.N ) ) THEN
         INFO = -12
      END IF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       CWorkspace refers to complex workspace, and RWorkspace to real
!       workspace. NB refers to the optimal block size for the
!       immediately following subroutine, as returned by ILAENV.
!       HSWORK refers to the workspace preferred by ZHSEQR, as
!       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
!       the worst case.)
!
      MINWRK = 1
      IF( INFO.EQ.0 .AND. LWORK.GE.1 ) THEN
         MAXWRK = N + N*ILAENV( 1, 'ZGEHRD', ' ', N, 1, N, 0 )
         IF( ( .NOT.WANTVL ) .AND. ( .NOT.WANTVR ) ) THEN
            MINWRK = MAX( 1, 2*N )
            IF( .NOT.( WNTSNN .OR. WNTSNE ) ) &
               MINWRK = MAX( MINWRK, N*N+2*N )
            MAXB = MAX( ILAENV( 8, 'ZHSEQR', 'SN', N, 1, N, -1 ), 2 )
            IF( WNTSNN ) THEN
               K = MIN( MAXB, N, MAX( 2, ILAENV( 4, 'ZHSEQR', 'EN', N, &
                   1, N, -1 ) ) )
            ELSE
               K = MIN( MAXB, N, MAX( 2, ILAENV( 4, 'ZHSEQR', 'SN', N, &
                   1, N, -1 ) ) )
            END IF
            HSWORK = MAX( K*( K+2 ), 2*N )
            MAXWRK = MAX( MAXWRK, 1, HSWORK )
            IF( .NOT.( WNTSNN .OR. WNTSNE ) ) &
               MAXWRK = MAX( MAXWRK, N*N+2*N )
         ELSE
            MINWRK = MAX( 1, 2*N )
            IF( .NOT.( WNTSNN .OR. WNTSNE ) ) &
               MINWRK = MAX( MINWRK, N*N+2*N )
            MAXB = MAX( ILAENV( 8, 'ZHSEQR', 'SN', N, 1, N, -1 ), 2 )
            K = MIN( MAXB, N, MAX( 2, ILAENV( 4, 'ZHSEQR', 'EN', N, 1, &
                N, -1 ) ) )
            HSWORK = MAX( K*( K+2 ), 2*N )
            MAXWRK = MAX( MAXWRK, 1, HSWORK )
            MAXWRK = MAX( MAXWRK, N+( N-1 )* &
                     ILAENV( 1, 'ZUNGHR', ' ', N, 1, N, -1 ) )
            IF( .NOT.( WNTSNN .OR. WNTSNE ) ) &
               MAXWRK = MAX( MAXWRK, N*N+2*N )
            MAXWRK = MAX( MAXWRK, 2*N, 1 )
         END IF
         WORK( 1 ) = MAXWRK
      END IF
      IF( LWORK.LT.MINWRK ) THEN
         INFO = -20
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEEVX', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Get machine constants
!
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      ICOND = 0
      ANRM = ZLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEA ) &
         CALL ZLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
!
!     Balance the matrix and compute ABNRM
!
      CALL ZGEBAL( BALANC, N, A, LDA, ILO, IHI, SCALE, IERR )
      ABNRM = ZLANGE( '1', N, N, A, LDA, DUM )
      IF( SCALEA ) THEN
         DUM( 1 ) = ABNRM
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, 1, 1, DUM, 1, IERR )
         ABNRM = DUM( 1 )
      END IF
!
!     Reduce to upper Hessenberg form
!     (CWorkspace: need 2*N, prefer N+N*NB)
!     (RWorkspace: none)
!
      ITAU = 1
      IWRK = ITAU + N
      CALL ZGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), &
                   LWORK-IWRK+1, IERR )
!
      IF( WANTVL ) THEN
!
!        Want left eigenvectors
!        Copy Householder vectors to VL
!
         SIDE = 'L'
         CALL ZLACPY( 'L', N, N, A, LDA, VL, LDVL )
!
!        Generate unitary matrix in VL
!        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
!        (RWorkspace: none)
!
         CALL ZUNGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ), &
                      LWORK-IWRK+1, IERR )
!
!        Perform QR iteration, accumulating Schur vectors in VL
!        (CWorkspace: need 1, prefer HSWORK (see comments) )
!        (RWorkspace: none)
!
         IWRK = ITAU
         CALL ZHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, W, VL, LDVL, &
                      WORK( IWRK ), LWORK-IWRK+1, INFO )
!
         IF( WANTVR ) THEN
!
!           Want left and right eigenvectors
!           Copy Schur vectors to VR
!
            SIDE = 'B'
            CALL ZLACPY( 'F', N, N, VL, LDVL, VR, LDVR )
         END IF
!
      ELSE IF( WANTVR ) THEN
!
!        Want right eigenvectors
!        Copy Householder vectors to VR
!
         SIDE = 'R'
         CALL ZLACPY( 'L', N, N, A, LDA, VR, LDVR )
!
!        Generate unitary matrix in VR
!        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
!        (RWorkspace: none)
!
         CALL ZUNGHR( N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ), &
                      LWORK-IWRK+1, IERR )
!
!        Perform QR iteration, accumulating Schur vectors in VR
!        (CWorkspace: need 1, prefer HSWORK (see comments) )
!        (RWorkspace: none)
!
         IWRK = ITAU
         CALL ZHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, W, VR, LDVR, &
                      WORK( IWRK ), LWORK-IWRK+1, INFO )
!
      ELSE
!
!        Compute eigenvalues only
!        If condition numbers desired, compute Schur form
!
         IF( WNTSNN ) THEN
            JOB = 'E'
         ELSE
            JOB = 'S'
         END IF
!
!        (CWorkspace: need 1, prefer HSWORK (see comments) )
!        (RWorkspace: none)
!
         IWRK = ITAU
         CALL ZHSEQR( JOB, 'N', N, ILO, IHI, A, LDA, W, VR, LDVR, &
                      WORK( IWRK ), LWORK-IWRK+1, INFO )
      END IF
!
!     If INFO > 0 from ZHSEQR, then quit
!
      IF( INFO.GT.0 ) &
         GO TO 50
!
      IF( WANTVL .OR. WANTVR ) THEN
!
!        Compute left and/or right eigenvectors
!        (CWorkspace: need 2*N)
!        (RWorkspace: need N)
!
         CALL ZTREVC( SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, &
                      N, NOUT, WORK( IWRK ), RWORK, IERR )
      END IF
!
!     Compute condition numbers if desired
!     (CWorkspace: need N*N+2*N unless SENSE = 'E')
!     (RWorkspace: need 2*N unless SENSE = 'E')
!
      IF( .NOT.WNTSNN ) THEN
         CALL ZTRSNA( SENSE, 'A', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, &
                      RCONDE, RCONDV, N, NOUT, WORK( IWRK ), N, RWORK, &
                      ICOND )
      END IF
!
      IF( WANTVL ) THEN
!
!        Undo balancing of left eigenvectors
!
         CALL ZGEBAK( BALANC, 'L', N, ILO, IHI, SCALE, N, VL, LDVL, &
                      IERR )
!
!        Normalize left eigenvectors and make largest component real
!
         DO 20 I = 1, N
            SCL = ONE / DZNRM2( N, VL( 1, I ), 1 )
            CALL ZDSCAL( N, SCL, VL( 1, I ), 1 )
            DO 10 K = 1, N
               RWORK( K ) = DBLE( VL( K, I ) )**2 + &
                            DIMAG( VL( K, I ) )**2
   10       CONTINUE
            K = IDAMAX( N, RWORK, 1 )
            TMP = DCONJG( VL( K, I ) ) / SQRT( RWORK( K ) )
            CALL ZSCAL( N, TMP, VL( 1, I ), 1 )
            VL( K, I ) = DCMPLX( DBLE( VL( K, I ) ), ZERO )
   20    CONTINUE
      END IF
!
      IF( WANTVR ) THEN
!
!        Undo balancing of right eigenvectors
!
         CALL ZGEBAK( BALANC, 'R', N, ILO, IHI, SCALE, N, VR, LDVR, &
                      IERR )
!
!        Normalize right eigenvectors and make largest component real
!
         DO 40 I = 1, N
            SCL = ONE / DZNRM2( N, VR( 1, I ), 1 )
            CALL ZDSCAL( N, SCL, VR( 1, I ), 1 )
            DO 30 K = 1, N
               RWORK( K ) = DBLE( VR( K, I ) )**2 + &
                            DIMAG( VR( K, I ) )**2
   30       CONTINUE
            K = IDAMAX( N, RWORK, 1 )
            TMP = DCONJG( VR( K, I ) ) / SQRT( RWORK( K ) )
            CALL ZSCAL( N, TMP, VR( 1, I ), 1 )
            VR( K, I ) = DCMPLX( DBLE( VR( K, I ) ), ZERO )
   40    CONTINUE
      END IF
!
!     Undo scaling if necessary
!
   50 CONTINUE
      IF( SCALEA ) THEN
         CALL ZLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, W( INFO+1 ), &
                      MAX( N-INFO, 1 ), IERR )
         IF( INFO.EQ.0 ) THEN
            IF( ( WNTSNV .OR. WNTSNB ) .AND. ICOND.EQ.0 ) &
               CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, RCONDV, N, &
                            IERR )
         ELSE
            CALL ZLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, W, N, IERR )
         END IF
      END IF
!
      WORK( 1 ) = MAXWRK
      RETURN
!
!     End of ZGEEVX
!
      END
      SUBROUTINE DLABAD( SMALL, LARGE )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   LARGE, SMALL
!     ..
!
!  Purpose
!  =======
!
!  DLABAD takes as input the values computed by SLAMCH for underflow and
!  overflow, and returns the square root of each of these values if the
!  log of LARGE is sufficiently large.  This subroutine is intended to
!  identify machines with a large exponent range, such as the Crays, and
!  redefine the underflow and overflow limits to be the square roots of
!  the values computed by DLAMCH.  This subroutine is needed because
!  DLAMCH does not compensate for poor arithmetic in the upper half of
!  the exponent range, as is found on a Cray.
!
!  Arguments
!  =========
!
!  SMALL   (input/output) DOUBLE PRECISION
!          On entry, the underflow threshold as computed by DLAMCH.
!          On exit, if LOG10(LARGE) is sufficiently large, the square
!          root of SMALL, otherwise unchanged.
!
!  LARGE   (input/output) DOUBLE PRECISION
!          On entry, the overflow threshold as computed by DLAMCH.
!          On exit, if LOG10(LARGE) is sufficiently large, the square
!          root of LARGE, otherwise unchanged.
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          LOG10, SQRT
!     ..
!     .. Executable Statements ..
!
!     If it looks like we're on a Cray, take the square root of
!     SMALL and LARGE to avoid overflow and underflow problems.
!
      IF( LOG10( LARGE ).GT.2000.D0 ) THEN
         SMALL = SQRT( SMALL )
         LARGE = SQRT( LARGE )
      END IF
!
      RETURN
!
!     End of DLABAD
!
      END
      SUBROUTINE ZGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            IHI, ILO, INFO, LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   SCALE( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGEBAL balances a general complex matrix A.  This involves, first,
!  permuting A by a similarity transformation to isolate eigenvalues
!  in the first 1 to ILO-1 and last IHI+1 to N elements on the
!  diagonal; and second, applying a diagonal similarity transformation
!  to rows and columns ILO to IHI to make the rows and columns as
!  close in norm as possible.  Both steps are optional.
!
!  Balancing may reduce the 1-norm of the matrix, and improve the
!  accuracy of the computed eigenvalues and/or eigenvectors.
!
!  Arguments
!  =========
!
!  JOB     (input) CHARACTER*1
!          Specifies the operations to be performed on A:
!          = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0
!                  for i = 1,...,N;
!          = 'P':  permute only;
!          = 'S':  scale only;
!          = 'B':  both permute and scale.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the input matrix A.
!          On exit,  A is overwritten by the balanced matrix.
!          If JOB = 'N', A is not referenced.
!          See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  ILO     (output) INTEGER
!  IHI     (output) INTEGER
!          ILO and IHI are set to integers such that on exit
!          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N.
!          If JOB = 'N' or 'S', ILO = 1 and IHI = N.
!
!  SCALE   (output) DOUBLE PRECISION array, dimension (N)
!          Details of the permutations and scaling factors applied to
!          A.  If P(j) is the index of the row and column interchanged
!          with row and column j and D(j) is the scaling factor
!          applied to row and column j, then
!          SCALE(j) = P(j)    for j = 1,...,ILO-1
!                   = D(j)    for j = ILO,...,IHI
!                   = P(j)    for j = IHI+1,...,N.
!          The order in which the interchanges are made is N to IHI+1,
!          then 1 to ILO-1.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The permutations consist of row and column interchanges which put
!  the matrix in the form
!
!             ( T1   X   Y  )
!     P A P = (  0   B   Z  )
!             (  0   0   T2 )
!
!  where T1 and T2 are upper triangular matrices whose eigenvalues lie
!  along the diagonal.  The column indices ILO and IHI mark the starting
!  and ending columns of the submatrix B. Balancing consists of applying
!  a diagonal similarity transformation inv(D) * B * D to make the
!  1-norms of each row of B and its corresponding column nearly equal.
!  The output matrix is
!
!     ( T1     X*D          Y    )
!     (  0  inv(D)*B*D  inv(D)*Z ).
!     (  0      0           T2   )
!
!  Information about the permutations P and the diagonal matrix D is
!  returned in the vector SCALE.
!
!  This subroutine is based on the EISPACK routine CBAL.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   SCLFAC
      PARAMETER          ( SCLFAC = 1.0D+1 )
      DOUBLE PRECISION   FACTOR
      PARAMETER          ( FACTOR = 0.95D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOCONV
      INTEGER            I, ICA, IEXC, IRA, J, K, L, M
      DOUBLE PRECISION   C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1, &
                         SFMIN2
      COMPLEX*16         CDUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IZAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, IZAMAX, DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZDSCAL, ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX, MIN
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. &
          .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEBAL', -INFO )
         RETURN
      END IF
!
      K = 1
      L = N
!
      IF( N.EQ.0 ) &
         GO TO 210
!
      IF( LSAME( JOB, 'N' ) ) THEN
         DO 10 I = 1, N
            SCALE( I ) = ONE
   10    CONTINUE
         GO TO 210
      END IF
!
      IF( LSAME( JOB, 'S' ) ) &
         GO TO 120
!
!     Permutation to isolate eigenvalues if possible
!
      GO TO 50
!
!     Row and column exchange.
!
   20 CONTINUE
      SCALE( M ) = J
      IF( J.EQ.M ) &
         GO TO 30
!
      CALL ZSWAP( L, A( 1, J ), 1, A( 1, M ), 1 )
      CALL ZSWAP( N-K+1, A( J, K ), LDA, A( M, K ), LDA )
!
   30 CONTINUE
      GO TO ( 40, 80 )IEXC
!
!     Search for rows isolating an eigenvalue and push them down.
!
   40 CONTINUE
      IF( L.EQ.1 ) &
         GO TO 210
      L = L - 1
!
   50 CONTINUE
      DO 70 J = L, 1, -1
!
         DO 60 I = 1, L
            IF( I.EQ.J ) &
               GO TO 60
            IF( DBLE( A( J, I ) ).NE.ZERO .OR. DIMAG( A( J, I ) ).NE. &
                ZERO )GO TO 70
   60    CONTINUE
!
         M = L
         IEXC = 1
         GO TO 20
   70 CONTINUE
!
      GO TO 90
!
!     Search for columns isolating an eigenvalue and push them left.
!
   80 CONTINUE
      K = K + 1
!
   90 CONTINUE
      DO 110 J = K, L
!
         DO 100 I = K, L
            IF( I.EQ.J ) &
               GO TO 100
            IF( DBLE( A( I, J ) ).NE.ZERO .OR. DIMAG( A( I, J ) ).NE. &
                ZERO )GO TO 110
  100    CONTINUE
!
         M = K
         IEXC = 2
         GO TO 20
  110 CONTINUE
!
  120 CONTINUE
      DO 130 I = K, L
         SCALE( I ) = ONE
  130 CONTINUE
!
      IF( LSAME( JOB, 'P' ) ) &
         GO TO 210
!
!     Balance the submatrix in rows K to L.
!
!     Iterative loop for norm reduction
!
      SFMIN1 = DLAMCH( 'S' ) / DLAMCH( 'P' )
      SFMAX1 = ONE / SFMIN1
      SFMIN2 = SFMIN1*SCLFAC
      SFMAX2 = ONE / SFMIN2
  140 CONTINUE
      NOCONV = .FALSE.
!
      DO 200 I = K, L
         C = ZERO
         R = ZERO
!
         DO 150 J = K, L
            IF( J.EQ.I ) &
               GO TO 150
            C = C + CABS1( A( J, I ) )
            R = R + CABS1( A( I, J ) )
  150    CONTINUE
         ICA = IZAMAX( L, A( 1, I ), 1 )
         CA = ABS( A( ICA, I ) )
         IRA = IZAMAX( N-K+1, A( I, K ), LDA )
         RA = ABS( A( I, IRA+K-1 ) )
!
!        Guard against zero C or R due to underflow.
!
         IF( C.EQ.ZERO .OR. R.EQ.ZERO ) &
            GO TO 200
         G = R / SCLFAC
         F = ONE
         S = C + R
  160    CONTINUE
         IF( C.GE.G .OR. MAX( F, C, CA ).GE.SFMAX2 .OR. &
             MIN( R, G, RA ).LE.SFMIN2 )GO TO 170
         F = F*SCLFAC
         C = C*SCLFAC
         CA = CA*SCLFAC
         R = R / SCLFAC
         G = G / SCLFAC
         RA = RA / SCLFAC
         GO TO 160
!
  170    CONTINUE
         G = C / SCLFAC
  180    CONTINUE
         IF( G.LT.R .OR. MAX( R, RA ).GE.SFMAX2 .OR. &
             MIN( F, C, G, CA ).LE.SFMIN2 )GO TO 190
         F = F / SCLFAC
         C = C / SCLFAC
         G = G / SCLFAC
         CA = CA / SCLFAC
         R = R*SCLFAC
         RA = RA*SCLFAC
         GO TO 180
!
!        Now balance.
!
  190    CONTINUE
         IF( ( C+R ).GE.FACTOR*S ) &
            GO TO 200
         IF( F.LT.ONE .AND. SCALE( I ).LT.ONE ) THEN
            IF( F*SCALE( I ).LE.SFMIN1 ) &
               GO TO 200
         END IF
         IF( F.GT.ONE .AND. SCALE( I ).GT.ONE ) THEN
            IF( SCALE( I ).GE.SFMAX1 / F ) &
               GO TO 200
         END IF
         G = ONE / F
         SCALE( I ) = SCALE( I )*F
         NOCONV = .TRUE.
!
         CALL ZDSCAL( N-K+1, G, A( I, K ), LDA )
         CALL ZDSCAL( L, F, A( 1, I ), 1 )
!
  200 CONTINUE
!
      IF( NOCONV ) &
         GO TO 140
!
  210 CONTINUE
      ILO = K
      IHI = L
!
      RETURN
!
!     End of ZGEBAL
!
      END
      SUBROUTINE ZGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( LWORK )
!     ..
!
!  Purpose
!  =======
!
!  ZGEHRD reduces a complex general matrix A to upper Hessenberg form H
!  by a unitary similarity transformation:  Q' * A * Q = H .
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          It is assumed that A is already upper triangular in rows
!          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!          set by a previous call to ZGEBAL; otherwise they should be
!          set to 1 and N respectively. See Further Details.
!          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the N-by-N general matrix to be reduced.
!          On exit, the upper triangle and the first subdiagonal of A
!          are overwritten with the upper Hessenberg matrix H, and the
!          elements below the first subdiagonal, with the array TAU,
!          represent the unitary matrix Q as a product of elementary
!          reflectors. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  TAU     (output) COMPLEX*16 array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to
!          zero.
!
!  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The length of the array WORK.  LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of (ihi-ilo) elementary
!  reflectors
!
!     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a complex scalar, and v is a complex vector with
!  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
!  exit in A(i+2:ihi,i), and tau in TAU(i).
!
!  The contents of A are illustrated by the following example, with
!  n = 7, ilo = 2 and ihi = 6:
!
!  on entry,                        on exit,
!
!  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
!  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
!  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
!  (                         a )    (                          a )
!
!  where a denotes an element of the original matrix A, h denotes a
!  modified element of the upper Hessenberg matrix H, and vi denotes an
!  element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), &
                         ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IB, IINFO, IWS, LDWORK, NB, NBMIN, NH, NX
      COMPLEX*16         EI
!     ..
!     .. Local Arrays ..
      COMPLEX*16         T( LDT, NBMAX )
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEHD2, ZGEMM, ZLAHRD, ZLARFB
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEHRD', -INFO )
         RETURN
      END IF
!
!     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
!
      DO 10 I = 1, ILO - 1
         TAU( I ) = ZERO
   10 CONTINUE
      DO 20 I = MAX( 1, IHI ), N - 1
         TAU( I ) = ZERO
   20 CONTINUE
!
!     Quick return if possible
!
      NH = IHI - ILO + 1
      IF( NH.LE.1 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
!     Determine the block size.
!
      NB = MIN( NBMAX, ILAENV( 1, 'ZGEHRD', ' ', N, ILO, IHI, -1 ) )
      NBMIN = 2
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.NH ) THEN
!
!        Determine when to cross over from blocked to unblocked code
!        (last block is always handled by unblocked code).
!
         NX = MAX( NB, ILAENV( 3, 'ZGEHRD', ' ', N, ILO, IHI, -1 ) )
         IF( NX.LT.NH ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            IWS = N*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  determine the
!              minimum value of NB, and reduce NB or force use of
!              unblocked code.
!
               NBMIN = MAX( 2, ILAENV( 2, 'ZGEHRD', ' ', N, ILO, IHI, &
                       -1 ) )
               IF( LWORK.GE.N*NBMIN ) THEN
                  NB = LWORK / N
               ELSE
                  NB = 1
               END IF
            END IF
         END IF
      END IF
      LDWORK = N
!
      IF( NB.LT.NBMIN .OR. NB.GE.NH ) THEN
!
!        Use unblocked code below
!
         I = ILO
!
      ELSE
!
!        Use blocked code
!
         DO 30 I = ILO, IHI - 1 - NX, NB
            IB = MIN( NB, IHI-I )
!
!           Reduce columns i:i+ib-1 to Hessenberg form, returning the
!           matrices V and T of the block reflector H = I - V*T*V'
!           which performs the reduction, and also the matrix Y = A*V*T
!
            CALL ZLAHRD( IHI, I, IB, A( 1, I ), LDA, TAU( I ), T, LDT, &
                         WORK, LDWORK )
!
!           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
!           right, computing  A := A - Y * V'. V(i+ib,ib-1) must be set
!           to 1.
!
            EI = A( I+IB, I+IB-1 )
            A( I+IB, I+IB-1 ) = ONE
            CALL ZGEMM( 'No transpose', 'Conjugate transpose', IHI, &
                        IHI-I-IB+1, IB, -ONE, WORK, LDWORK, &
                        A( I+IB, I ), LDA, ONE, A( 1, I+IB ), LDA )
            A( I+IB, I+IB-1 ) = EI
!
!           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
!           left
!
            CALL ZLARFB( 'Left', 'Conjugate transpose', 'Forward', &
                         'Columnwise', IHI-I, N-I-IB+1, IB, A( I+1, I ), &
                         LDA, T, LDT, A( I+1, I+IB ), LDA, WORK, &
                         LDWORK )
   30    CONTINUE
      END IF
!
!     Use unblocked code to reduce the rest of the matrix
!
      CALL ZGEHD2( N, I, IHI, A, LDA, TAU, WORK, IINFO )
      WORK( 1 ) = IWS
!
      RETURN
!
!     End of ZGEHRD
!
      END
      SUBROUTINE ZUNGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( LWORK )
!     ..
!
!  Purpose
!  =======
!
!  ZUNGHR generates a complex unitary matrix Q which is defined as the
!  product of IHI-ILO elementary reflectors of order N, as returned by
!  ZGEHRD:
!
!  Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix Q. N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          ILO and IHI must have the same values as in the previous call
!          of ZGEHRD. Q is equal to the unit matrix except in the
!          submatrix Q(ilo+1:ihi,ilo+1:ihi).
!          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the vectors which define the elementary reflectors,
!          as returned by ZGEHRD.
!          On exit, the N-by-N unitary matrix Q.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,N).
!
!  TAU     (input) COMPLEX*16 array, dimension (N-1)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by ZGEHRD.
!
!  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= IHI-ILO.
!          For optimum performance LWORK >= (IHI-ILO)*NB, where NB is
!          the optimal blocksize.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), &
                         ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IINFO, J, NH
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZUNGQR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, IHI-ILO ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGHR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
!     Shift the vectors which define the elementary reflectors one
!     column to the right, and set the first ilo and the last n-ihi
!     rows and columns to those of the unit matrix
!
      DO 40 J = IHI, ILO + 1, -1
         DO 10 I = 1, J - 1
            A( I, J ) = ZERO
   10    CONTINUE
         DO 20 I = J + 1, IHI
            A( I, J ) = A( I, J-1 )
   20    CONTINUE
         DO 30 I = IHI + 1, N
            A( I, J ) = ZERO
   30    CONTINUE
   40 CONTINUE
      DO 60 J = 1, ILO
         DO 50 I = 1, N
            A( I, J ) = ZERO
   50    CONTINUE
         A( J, J ) = ONE
   60 CONTINUE
      DO 80 J = IHI + 1, N
         DO 70 I = 1, N
            A( I, J ) = ZERO
   70    CONTINUE
         A( J, J ) = ONE
   80 CONTINUE
!
      NH = IHI - ILO
      IF( NH.GT.0 ) THEN
!
!        Generate Q(ilo+1:ihi,ilo+1:ihi)
!
         CALL ZUNGQR( NH, NH, NH, A( ILO+1, ILO+1 ), LDA, TAU( ILO ), &
                      WORK, LWORK, IINFO )
      END IF
      RETURN
!
!     End of ZUNGHR
!
      END
      SUBROUTINE ZHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ, &
                         WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          COMPZ, JOB
      INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  ZHSEQR computes the eigenvalues of a complex upper Hessenberg
!  matrix H, and, optionally, the matrices T and Z from the Schur
!  decomposition H = Z T Z**H, where T is an upper triangular matrix
!  (the Schur form), and Z is the unitary matrix of Schur vectors.
!
!  Optionally Z may be postmultiplied into an input unitary matrix Q,
!  so that this routine can give the Schur factorization of a matrix A
!  which has been reduced to the Hessenberg form H by the unitary
!  matrix Q:  A = Q*H*Q**H = (QZ)*T*(QZ)**H.
!
!  Arguments
!  =========
!
!  JOB     (input) CHARACTER*1
!          = 'E': compute eigenvalues only;
!          = 'S': compute eigenvalues and the Schur form T.
!
!  COMPZ   (input) CHARACTER*1
!          = 'N': no Schur vectors are computed;
!          = 'I': Z is initialized to the unit matrix and the matrix Z
!                 of Schur vectors of H is returned;
!          = 'V': Z must contain an unitary matrix Q on entry, and
!                 the product Q*Z is returned.
!
!  N       (input) INTEGER
!          The order of the matrix H.  N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          It is assumed that H is already upper triangular in rows
!          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!          set by a previous call to ZGEBAL, and then passed to CGEHRD
!          when the matrix output by ZGEBAL is reduced to Hessenberg
!          form. Otherwise ILO and IHI should be set to 1 and N
!          respectively.
!          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!
!  H       (input/output) COMPLEX*16 array, dimension (LDH,N)
!          On entry, the upper Hessenberg matrix H.
!          On exit, if JOB = 'S', H contains the upper triangular matrix
!          T from the Schur decomposition (the Schur form). If
!          JOB = 'E', the contents of H are unspecified on exit.
!
!  LDH     (input) INTEGER
!          The leading dimension of the array H. LDH >= max(1,N).
!
!  W       (output) COMPLEX*16 array, dimension (N)
!          The computed eigenvalues. If JOB = 'S', the eigenvalues are
!          stored in the same order as on the diagonal of the Schur form
!          returned in H, with W(i) = H(i,i).
!
!  Z       (input/output) COMPLEX*16 array, dimension (LDZ,N)
!          If COMPZ = 'N': Z is not referenced.
!          If COMPZ = 'I': on entry, Z need not be set, and on exit, Z
!          contains the unitary matrix Z of the Schur vectors of H.
!          If COMPZ = 'V': on entry Z must contain an N-by-N matrix Q,
!          which is assumed to be equal to the unit matrix except for
!          the submatrix Z(ILO:IHI,ILO:IHI); on exit Z contains Q*Z.
!          Normally Q is the unitary matrix generated by ZUNGHR after
!          the call to ZGEHRD which formed the Hessenberg matrix H.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.
!          LDZ >= max(1,N) if COMPZ = 'I' or 'V'; LDZ >= 1 otherwise.
!
!  WORK    (workspace) COMPLEX*16 array, dimension (N)
!
!  LWORK   (input) INTEGER
!          This argument is currently redundant.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, ZHSEQR failed to compute all the
!                eigenvalues in a total of 30*(IHI-ILO+1) iterations;
!                elements 1:ilo-1 and i+1:n of W contain those
!                eigenvalues which have been successfully computed.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), &
                         ONE = ( 1.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION   RZERO, RONE, CONST
      PARAMETER          ( RZERO = 0.0D+0, RONE = 1.0D+0, &
                         CONST = 1.5D+0 )
      INTEGER            NSMAX, LDS
      PARAMETER          ( NSMAX = 15, LDS = NSMAX )
!     ..
!     .. Local Scalars ..
      LOGICAL            INITZ, WANTT, WANTZ
      INTEGER            I, I1, I2, IERR, II, ITEMP, ITN, ITS, J, K, L, &
                         MAXB, NH, NR, NS, NV
      DOUBLE PRECISION   OVFL, RTEMP, SMLNUM, TST1, ULP, UNFL
      COMPLEX*16         CDUM, TAU, TEMP
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   RWORK( 1 )
      COMPLEX*16         S( LDS, NSMAX ), V( NSMAX+1 ), VV( NSMAX+1 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV, IZAMAX
      DOUBLE PRECISION   DLAMCH, DLAPY2, ZLANHS
      EXTERNAL           LSAME, ILAENV, IZAMAX, DLAMCH, DLAPY2, ZLANHS
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLABAD, XERBLA, ZCOPY, ZDSCAL, ZGEMV, ZLACPY, &
                         ZLAHQR, ZLARFG, ZLARFX, ZLASET, ZSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, MAX, MIN
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      WANTT = LSAME( JOB, 'S' )
      INITZ = LSAME( COMPZ, 'I' )
      WANTZ = INITZ .OR. LSAME( COMPZ, 'V' )
!
      INFO = 0
      IF( .NOT.LSAME( JOB, 'E' ) .AND. .NOT.WANTT ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( COMPZ, 'N' ) .AND. .NOT.WANTZ ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -5
      ELSE IF( LDH.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDZ.LT.1 .OR. WANTZ .AND. LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHSEQR', -INFO )
         RETURN
      END IF
!
!     Initialize Z, if necessary
!
      IF( INITZ ) &
         CALL ZLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
!
!     Store the eigenvalues isolated by ZGEBAL.
!
      DO 10 I = 1, ILO - 1
         W( I ) = H( I, I )
   10 CONTINUE
      DO 20 I = IHI + 1, N
         W( I ) = H( I, I )
   20 CONTINUE
!
!     Quick return if possible.
!
      IF( N.EQ.0 ) &
         RETURN
      IF( ILO.EQ.IHI ) THEN
         W( ILO ) = H( ILO, ILO )
         RETURN
      END IF
!
!     Set rows and columns ILO to IHI to zero below the first
!     subdiagonal.
!
      DO 40 J = ILO, IHI - 2
         DO 30 I = J + 2, N
            H( I, J ) = ZERO
   30    CONTINUE
   40 CONTINUE
      NH = IHI - ILO + 1
!
!     I1 and I2 are the indices of the first row and last column of H
!     to which transformations must be applied. If eigenvalues only are
!     being computed, I1 and I2 are re-set inside the main loop.
!
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      ELSE
         I1 = ILO
         I2 = IHI
      END IF
!
!     Ensure that the subdiagonal elements are real.
!
      DO 50 I = ILO + 1, IHI
         TEMP = H( I, I-1 )
         IF( DIMAG( TEMP ).NE.RZERO ) THEN
            RTEMP = DLAPY2( DBLE( TEMP ), DIMAG( TEMP ) )
            H( I, I-1 ) = RTEMP
            TEMP = TEMP / RTEMP
            IF( I2.GT.I ) &
               CALL ZSCAL( I2-I, DCONJG( TEMP ), H( I, I+1 ), LDH )
            CALL ZSCAL( I-I1, TEMP, H( I1, I ), 1 )
            IF( I.LT.IHI ) &
               H( I+1, I ) = TEMP*H( I+1, I )
            IF( WANTZ ) &
               CALL ZSCAL( NH, TEMP, Z( ILO, I ), 1 )
         END IF
   50 CONTINUE
!
!     Determine the order of the multi-shift QR algorithm to be used.
!
      NS = ILAENV( 4, 'ZHSEQR', JOB // COMPZ, N, ILO, IHI, -1 )
      MAXB = ILAENV( 8, 'ZHSEQR', JOB // COMPZ, N, ILO, IHI, -1 )
      IF( NS.LE.1 .OR. NS.GT.NH .OR. MAXB.GE.NH ) THEN
!
!        Use the standard double-shift algorithm
!
         CALL ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, &
                      LDZ, INFO )
         RETURN
      END IF
      MAXB = MAX( 2, MAXB )
      NS = MIN( NS, MAXB, NSMAX )
!
!     Now 1 < NS <= MAXB < NH.
!
!     Set machine-dependent constants for the stopping criterion.
!     If norm(H) <= sqrt(OVFL), overflow should not occur.
!
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = RONE / UNFL
      CALL DLABAD( UNFL, OVFL )
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( NH / ULP )
!
!     ITN is the total number of multiple-shift QR iterations allowed.
!
      ITN = 30*NH
!
!     The main loop begins here. I is the loop index and decreases from
!     IHI to ILO in steps of at most MAXB. Each iteration of the loop
!     works with the active submatrix in rows and columns L to I.
!     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
!     H(L,L-1) is negligible so that the matrix splits.
!
      I = IHI
   60 CONTINUE
      IF( I.LT.ILO ) &
         GO TO 180
!
!     Perform multiple-shift QR iterations on rows and columns ILO to I
!     until a submatrix of order at most MAXB splits off at the bottom
!     because a subdiagonal element has become negligible.
!
      L = ILO
      DO 160 ITS = 0, ITN
!
!        Look for a single small subdiagonal element.
!
         DO 70 K = I, L + 1, -1
            TST1 = CABS1( H( K-1, K-1 ) ) + CABS1( H( K, K ) )
            IF( TST1.EQ.RZERO ) &
               TST1 = ZLANHS( '1', I-L+1, H( L, L ), LDH, RWORK )
            IF( ABS( DBLE( H( K, K-1 ) ) ).LE.MAX( ULP*TST1, SMLNUM ) ) &
               GO TO 80
   70    CONTINUE
   80    CONTINUE
         L = K
         IF( L.GT.ILO ) THEN
!
!           H(L,L-1) is negligible.
!
            H( L, L-1 ) = ZERO
         END IF
!
!        Exit from loop if a submatrix of order <= MAXB has split off.
!
         IF( L.GE.I-MAXB+1 ) &
            GO TO 170
!
!        Now the active submatrix is in rows and columns L to I. If
!        eigenvalues only are being computed, only the active submatrix
!        need be transformed.
!
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
!
         IF( ITS.EQ.20 .OR. ITS.EQ.30 ) THEN
!
!           Exceptional shifts.
!
            DO 90 II = I - NS + 1, I
               W( II ) = CONST*( ABS( DBLE( H( II, II-1 ) ) )+ &
                         ABS( DBLE( H( II, II ) ) ) )
   90       CONTINUE
         ELSE
!
!           Use eigenvalues of trailing submatrix of order NS as shifts.
!
            CALL ZLACPY( 'Full', NS, NS, H( I-NS+1, I-NS+1 ), LDH, S, &
                         LDS )
            CALL ZLAHQR( .FALSE., .FALSE., NS, 1, NS, S, LDS, &
                         W( I-NS+1 ), 1, NS, Z, LDZ, IERR )
            IF( IERR.GT.0 ) THEN
!
!              If ZLAHQR failed to compute all NS eigenvalues, use the
!              unconverged diagonal elements as the remaining shifts.
!
               DO 100 II = 1, IERR
                  W( I-NS+II ) = S( II, II )
  100          CONTINUE
            END IF
         END IF
!
!        Form the first column of (G-w(1)) (G-w(2)) . . . (G-w(ns))
!        where G is the Hessenberg submatrix H(L:I,L:I) and w is
!        the vector of shifts (stored in W). The result is
!        stored in the local array V.
!
         V( 1 ) = ONE
         DO 110 II = 2, NS + 1
            V( II ) = ZERO
  110    CONTINUE
         NV = 1
         DO 130 J = I - NS + 1, I
            CALL ZCOPY( NV+1, V, 1, VV, 1 )
            CALL ZGEMV( 'No transpose', NV+1, NV, ONE, H( L, L ), LDH, &
                        VV, 1, -W( J ), V, 1 )
            NV = NV + 1
!
!           Scale V(1:NV) so that max(abs(V(i))) = 1. If V is zero,
!           reset it to the unit vector.
!
            ITEMP = IZAMAX( NV, V, 1 )
            RTEMP = CABS1( V( ITEMP ) )
            IF( RTEMP.EQ.RZERO ) THEN
               V( 1 ) = ONE
               DO 120 II = 2, NV
                  V( II ) = ZERO
  120          CONTINUE
            ELSE
               RTEMP = MAX( RTEMP, SMLNUM )
               CALL ZDSCAL( NV, RONE / RTEMP, V, 1 )
            END IF
  130    CONTINUE
!
!        Multiple-shift QR step
!
         DO 150 K = L, I - 1
!
!           The first iteration of this loop determines a reflection G
!           from the vector V and applies it from left and right to H,
!           thus creating a nonzero bulge below the subdiagonal.
!
!           Each subsequent iteration determines a reflection G to
!           restore the Hessenberg form in the (K-1)th column, and thus
!           chases the bulge one step toward the bottom of the active
!           submatrix. NR is the order of G.
!
            NR = MIN( NS+1, I-K+1 )
            IF( K.GT.L ) &
               CALL ZCOPY( NR, H( K, K-1 ), 1, V, 1 )
            CALL ZLARFG( NR, V( 1 ), V( 2 ), 1, TAU )
            IF( K.GT.L ) THEN
               H( K, K-1 ) = V( 1 )
               DO 140 II = K + 1, I
                  H( II, K-1 ) = ZERO
  140          CONTINUE
            END IF
            V( 1 ) = ONE
!
!           Apply G' from the left to transform the rows of the matrix
!           in columns K to I2.
!
            CALL ZLARFX( 'Left', NR, I2-K+1, V, DCONJG( TAU ), &
                         H( K, K ), LDH, WORK )
!
!           Apply G from the right to transform the columns of the
!           matrix in rows I1 to min(K+NR,I).
!
            CALL ZLARFX( 'Right', MIN( K+NR, I )-I1+1, NR, V, TAU, &
                         H( I1, K ), LDH, WORK )
!
            IF( WANTZ ) THEN
!
!              Accumulate transformations in the matrix Z
!
               CALL ZLARFX( 'Right', NH, NR, V, TAU, Z( ILO, K ), LDZ, &
                            WORK )
            END IF
  150    CONTINUE
!
!        Ensure that H(I,I-1) is real.
!
         TEMP = H( I, I-1 )
         IF( DIMAG( TEMP ).NE.RZERO ) THEN
            RTEMP = DLAPY2( DBLE( TEMP ), DIMAG( TEMP ) )
            H( I, I-1 ) = RTEMP
            TEMP = TEMP / RTEMP
            IF( I2.GT.I ) &
               CALL ZSCAL( I2-I, DCONJG( TEMP ), H( I, I+1 ), LDH )
            CALL ZSCAL( I-I1, TEMP, H( I1, I ), 1 )
            IF( WANTZ ) THEN
               CALL ZSCAL( NH, TEMP, Z( ILO, I ), 1 )
            END IF
         END IF
!
  160 CONTINUE
!
!     Failure to converge in remaining number of iterations
!
      INFO = I
      RETURN
!
  170 CONTINUE
!
!     A submatrix of order <= MAXB in rows and columns L to I has split
!     off. Use the double-shift QR algorithm to handle it.
!
      CALL ZLAHQR( WANTT, WANTZ, N, L, I, H, LDH, W, ILO, IHI, Z, LDZ, &
                   INFO )
      IF( INFO.GT.0 ) &
         RETURN
!
!     Decrement number of remaining iterations, and return to start of
!     the main loop with a new value of I.
!
      ITN = ITN - ITS
      I = L - 1
      GO TO 60
!
  180 CONTINUE
      RETURN
!
!     End of ZHSEQR
!
      END
      SUBROUTINE ZTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, &
                         LDVR, MM, M, WORK, RWORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N
!     ..
!     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), &
                         WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZTREVC computes some or all of the right and/or left eigenvectors of
!  a complex upper triangular matrix T.
!
!  The right eigenvector x and the left eigenvector y of T corresponding
!  to an eigenvalue w are defined by:
!
!               T*x = w*x,     y'*T = w*y'
!
!  where y' denotes the conjugate transpose of the vector y.
!
!  If all eigenvectors are requested, the routine may either return the
!  matrices X and/or Y of right or left eigenvectors of T, or the
!  products Q*X and/or Q*Y, where Q is an input unitary
!  matrix. If T was obtained from the Schur factorization of an
!  original matrix A = Q*T*Q', then Q*X and Q*Y are the matrices of
!  right or left eigenvectors of A.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'R':  compute right eigenvectors only;
!          = 'L':  compute left eigenvectors only;
!          = 'B':  compute both right and left eigenvectors.
!
!  HOWMNY  (input) CHARACTER*1
!          = 'A':  compute all right and/or left eigenvectors;
!          = 'B':  compute all right and/or left eigenvectors,
!                  and backtransform them using the input matrices
!                  supplied in VR and/or VL;
!          = 'S':  compute selected right and/or left eigenvectors,
!                  specified by the logical array SELECT.
!
!  SELECT  (input) LOGICAL array, dimension (N)
!          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
!          computed.
!          If HOWMNY = 'A' or 'B', SELECT is not referenced.
!          To select the eigenvector corresponding to the j-th
!          eigenvalue, SELECT(j) must be set to .TRUE..
!
!  N       (input) INTEGER
!          The order of the matrix T. N >= 0.
!
!  T       (input/output) COMPLEX*16 array, dimension (LDT,N)
!          The upper triangular matrix T.  T is modified, but restored
!          on exit.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= max(1,N).
!
!  VL      (input/output) COMPLEX*16 array, dimension (LDVL,MM)
!          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
!          contain an N-by-N matrix Q (usually the unitary matrix Q of
!          Schur vectors returned by ZHSEQR).
!          On exit, if SIDE = 'L' or 'B', VL contains:
!          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
!          if HOWMNY = 'B', the matrix Q*Y;
!          if HOWMNY = 'S', the left eigenvectors of T specified by
!                           SELECT, stored consecutively in the columns
!                           of VL, in the same order as their
!                           eigenvalues.
!          If SIDE = 'R', VL is not referenced.
!
!  LDVL    (input) INTEGER
!          The leading dimension of the array VL.  LDVL >= max(1,N) if
!          SIDE = 'L' or 'B'; LDVL >= 1 otherwise.
!
!  VR      (input/output) COMPLEX*16 array, dimension (LDVR,MM)
!          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
!          contain an N-by-N matrix Q (usually the unitary matrix Q of
!          Schur vectors returned by ZHSEQR).
!          On exit, if SIDE = 'R' or 'B', VR contains:
!          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
!          if HOWMNY = 'B', the matrix Q*X;
!          if HOWMNY = 'S', the right eigenvectors of T specified by
!                           SELECT, stored consecutively in the columns
!                           of VR, in the same order as their
!                           eigenvalues.
!          If SIDE = 'L', VR is not referenced.
!
!  LDVR    (input) INTEGER
!          The leading dimension of the array VR.  LDVR >= max(1,N) if
!           SIDE = 'R' or 'B'; LDVR >= 1 otherwise.
!
!  MM      (input) INTEGER
!          The number of columns in the arrays VL and/or VR. MM >= M.
!
!  M       (output) INTEGER
!          The number of columns in the arrays VL and/or VR actually
!          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M
!          is set to N.  Each selected eigenvector occupies one
!          column.
!
!  WORK    (workspace) COMPLEX*16 array, dimension (2*N)
!
!  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The algorithm used in this program is basically backward (forward)
!  substitution, with scaling to make the the code robust against
!  possible overflow.
!
!  Each eigenvector is normalized so that the element of largest
!  magnitude has magnitude 1; here the magnitude of a complex number
!  (x,y) is taken to be |x| + |y|.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CMZERO, CMONE
      PARAMETER          ( CMZERO = ( 0.0D+0, 0.0D+0 ), &
                         CMONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            ALLV, BOTHV, LEFTV, OVER, RIGHTV, SOMEV
      INTEGER            I, II, IS, J, K, KI
      DOUBLE PRECISION   OVFL, REMAX, SCALE, SMIN, SMLNUM, ULP, UNFL
      COMPLEX*16         CDUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IZAMAX
      DOUBLE PRECISION   DLAMCH, DZASUM
      EXTERNAL           LSAME, IZAMAX, DLAMCH, DZASUM
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLABAD, XERBLA, ZCOPY, ZDSCAL, ZGEMV, ZLATRS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      BOTHV = LSAME( SIDE, 'B' )
      RIGHTV = LSAME( SIDE, 'R' ) .OR. BOTHV
      LEFTV = LSAME( SIDE, 'L' ) .OR. BOTHV
!
      ALLV = LSAME( HOWMNY, 'A' )
      OVER = LSAME( HOWMNY, 'B' ) .OR. LSAME( HOWMNY, 'O' )
      SOMEV = LSAME( HOWMNY, 'S' )
!
!     Set M to the number of columns required to store the selected
!     eigenvectors.
!
      IF( SOMEV ) THEN
         M = 0
         DO 10 J = 1, N
            IF( SELECT( J ) ) &
               M = M + 1
   10    CONTINUE
      ELSE
         M = N
      END IF
!
      INFO = 0
      IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -1
      ELSE IF( .NOT.ALLV .AND. .NOT.OVER .AND. .NOT.SOMEV ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDVL.LT.1 .OR. ( LEFTV .AND. LDVL.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LDVR.LT.1 .OR. ( RIGHTV .AND. LDVR.LT.N ) ) THEN
         INFO = -10
      ELSE IF( MM.LT.M ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTREVC', -INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Set the constants to control overflow.
!
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL DLABAD( UNFL, OVFL )
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )
!
!     Store the diagonal elements of T in working array WORK.
!
      DO 20 I = 1, N
         WORK( I+N ) = T( I, I )
   20 CONTINUE
!
!     Compute 1-norm of each column of strictly upper triangular
!     part of T to control overflow in triangular solver.
!
      RWORK( 1 ) = ZERO
      DO 30 J = 2, N
         RWORK( J ) = DZASUM( J-1, T( 1, J ), 1 )
   30 CONTINUE
!
      IF( RIGHTV ) THEN
!
!        Compute right eigenvectors.
!
         IS = M
         DO 80 KI = N, 1, -1
!
            IF( SOMEV ) THEN
               IF( .NOT.SELECT( KI ) ) &
                  GO TO 80
            END IF
            SMIN = MAX( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM )
!
            WORK( 1 ) = CMONE
!
!           Form right-hand side.
!
            DO 40 K = 1, KI - 1
               WORK( K ) = -T( K, KI )
   40       CONTINUE
!
!           Solve the triangular system:
!              (T(1:KI-1,1:KI-1) - T(KI,KI))*X = SCALE*WORK.
!
            DO 50 K = 1, KI - 1
               T( K, K ) = T( K, K ) - T( KI, KI )
               IF( CABS1( T( K, K ) ).LT.SMIN ) &
                  T( K, K ) = SMIN
   50       CONTINUE
!
            IF( KI.GT.1 ) THEN
               CALL ZLATRS( 'Upper', 'No transpose', 'Non-unit', 'Y', &
                            KI-1, T, LDT, WORK( 1 ), SCALE, RWORK, &
                            INFO )
               WORK( KI ) = SCALE
            END IF
!
!           Copy the vector x or Q*x to VR and normalize.
!
            IF( .NOT.OVER ) THEN
               CALL ZCOPY( KI, WORK( 1 ), 1, VR( 1, IS ), 1 )
!
               II = IZAMAX( KI, VR( 1, IS ), 1 )
               REMAX = ONE / CABS1( VR( II, IS ) )
               CALL ZDSCAL( KI, REMAX, VR( 1, IS ), 1 )
!
               DO 60 K = KI + 1, N
                  VR( K, IS ) = CMZERO
   60          CONTINUE
            ELSE
               IF( KI.GT.1 ) &
                  CALL ZGEMV( 'N', N, KI-1, CMONE, VR, LDVR, WORK( 1 ), &
                              1, DCMPLX( SCALE ), VR( 1, KI ), 1 )
!
               II = IZAMAX( N, VR( 1, KI ), 1 )
               REMAX = ONE / CABS1( VR( II, KI ) )
               CALL ZDSCAL( N, REMAX, VR( 1, KI ), 1 )
            END IF
!
!           Set back the original diagonal elements of T.
!
            DO 70 K = 1, KI - 1
               T( K, K ) = WORK( K+N )
   70       CONTINUE
!
            IS = IS - 1
   80    CONTINUE
      END IF
!
      IF( LEFTV ) THEN
!
!        Compute left eigenvectors.
!
         IS = 1
         DO 130 KI = 1, N
!
            IF( SOMEV ) THEN
               IF( .NOT.SELECT( KI ) ) &
                  GO TO 130
            END IF
            SMIN = MAX( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM )
!
            WORK( N ) = CMONE
!
!           Form right-hand side.
!
            DO 90 K = KI + 1, N
               WORK( K ) = -DCONJG( T( KI, K ) )
   90       CONTINUE
!
!           Solve the triangular system:
!              (T(KI+1:N,KI+1:N) - T(KI,KI))'*X = SCALE*WORK.
!
            DO 100 K = KI + 1, N
               T( K, K ) = T( K, K ) - T( KI, KI )
               IF( CABS1( T( K, K ) ).LT.SMIN ) &
                  T( K, K ) = SMIN
  100       CONTINUE
!
            IF( KI.LT.N ) THEN
               CALL ZLATRS( 'Upper', 'Conjugate transpose', 'Non-unit', &
                            'Y', N-KI, T( KI+1, KI+1 ), LDT, &
                            WORK( KI+1 ), SCALE, RWORK, INFO )
               WORK( KI ) = SCALE
            END IF
!
!           Copy the vector x or Q*x to VL and normalize.
!
            IF( .NOT.OVER ) THEN
               CALL ZCOPY( N-KI+1, WORK( KI ), 1, VL( KI, IS ), 1 )
!
               II = IZAMAX( N-KI+1, VL( KI, IS ), 1 ) + KI - 1
               REMAX = ONE / CABS1( VL( II, IS ) )
               CALL ZDSCAL( N-KI+1, REMAX, VL( KI, IS ), 1 )
!
               DO 110 K = 1, KI - 1
                  VL( K, IS ) = CMZERO
  110          CONTINUE
            ELSE
               IF( KI.LT.N ) &
                  CALL ZGEMV( 'N', N, N-KI, CMONE, VL( 1, KI+1 ), LDVL, &
                              WORK( KI+1 ), 1, DCMPLX( SCALE ), &
                              VL( 1, KI ), 1 )
!
               II = IZAMAX( N, VL( 1, KI ), 1 )
               REMAX = ONE / CABS1( VL( II, KI ) )
               CALL ZDSCAL( N, REMAX, VL( 1, KI ), 1 )
            END IF
!
!           Set back the original diagonal elements of T.
!
            DO 120 K = KI + 1, N
               T( K, K ) = WORK( K+N )
  120       CONTINUE
!
            IS = IS + 1
  130    CONTINUE
      END IF
!
      RETURN
!
!     End of ZTREVC
!
      END
      SUBROUTINE ZTRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, &
                         LDVR, S, SEP, MM, M, WORK, LDWORK, RWORK, &
                         INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          HOWMNY, JOB
      INTEGER            INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N
!     ..
!     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   RWORK( * ), S( * ), SEP( * )
      COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), &
                         WORK( LDWORK, * )
!     ..
!
!  Purpose
!  =======
!
!  ZTRSNA estimates reciprocal condition numbers for specified
!  eigenvalues and/or right eigenvectors of a complex upper triangular
!  matrix T (or of any matrix Q*T*Q**H with Q unitary).
!
!  Arguments
!  =========
!
!  JOB     (input) CHARACTER*1
!          Specifies whether condition numbers are required for
!          eigenvalues (S) or eigenvectors (SEP):
!          = 'E': for eigenvalues only (S);
!          = 'V': for eigenvectors only (SEP);
!          = 'B': for both eigenvalues and eigenvectors (S and SEP).
!
!  HOWMNY  (input) CHARACTER*1
!          = 'A': compute condition numbers for all eigenpairs;
!          = 'S': compute condition numbers for selected eigenpairs
!                 specified by the array SELECT.
!
!  SELECT  (input) LOGICAL array, dimension (N)
!          If HOWMNY = 'S', SELECT specifies the eigenpairs for which
!          condition numbers are required. To select condition numbers
!          for the j-th eigenpair, SELECT(j) must be set to .TRUE..
!          If HOWMNY = 'A', SELECT is not referenced.
!
!  N       (input) INTEGER
!          The order of the matrix T. N >= 0.
!
!  T       (input) COMPLEX*16 array, dimension (LDT,N)
!          The upper triangular matrix T.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= max(1,N).
!
!  VL      (input) COMPLEX*16 array, dimension (LDVL,M)
!          If JOB = 'E' or 'B', VL must contain left eigenvectors of T
!          (or of any Q*T*Q**H with Q unitary), corresponding to the
!          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
!          must be stored in consecutive columns of VL, as returned by
!          ZHSEIN or ZTREVC.
!          If JOB = 'V', VL is not referenced.
!
!  LDVL    (input) INTEGER
!          The leading dimension of the array VL.
!          LDVL >= 1; and if JOB = 'E' or 'B', LDVL >= N.
!
!  VR      (input) COMPLEX*16 array, dimension (LDVR,M)
!          If JOB = 'E' or 'B', VR must contain right eigenvectors of T
!          (or of any Q*T*Q**H with Q unitary), corresponding to the
!          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
!          must be stored in consecutive columns of VR, as returned by
!          ZHSEIN or ZTREVC.
!          If JOB = 'V', VR is not referenced.
!
!  LDVR    (input) INTEGER
!          The leading dimension of the array VR.
!          LDVR >= 1; and if JOB = 'E' or 'B', LDVR >= N.
!
!  S       (output) DOUBLE PRECISION array, dimension (MM)
!          If JOB = 'E' or 'B', the reciprocal condition numbers of the
!          selected eigenvalues, stored in consecutive elements of the
!          array. Thus S(j), SEP(j), and the j-th columns of VL and VR
!          all correspond to the same eigenpair (but not in general the
!          j-th eigenpair, unless all eigenpairs are selected).
!          If JOB = 'V', S is not referenced.
!
!  SEP     (output) DOUBLE PRECISION array, dimension (MM)
!          If JOB = 'V' or 'B', the estimated reciprocal condition
!          numbers of the selected eigenvectors, stored in consecutive
!          elements of the array.
!          If JOB = 'E', SEP is not referenced.
!
!  MM      (input) INTEGER
!          The number of elements in the arrays S (if JOB = 'E' or 'B')
!           and/or SEP (if JOB = 'V' or 'B'). MM >= M.
!
!  M       (output) INTEGER
!          The number of elements of the arrays S and/or SEP actually
!          used to store the estimated condition numbers.
!          If HOWMNY = 'A', M is set to N.
!
!  WORK    (workspace) COMPLEX*16 array, dimension (LDWORK,N+1)
!          If JOB = 'E', WORK is not referenced.
!
!  LDWORK  (input) INTEGER
!          The leading dimension of the array WORK.
!          LDWORK >= 1; and if JOB = 'V' or 'B', LDWORK >= N.
!
!  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
!          If JOB = 'E', RWORK is not referenced.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The reciprocal of the condition number of an eigenvalue lambda is
!  defined as
!
!          S(lambda) = |v'*u| / (norm(u)*norm(v))
!
!  where u and v are the right and left eigenvectors of T corresponding
!  to lambda; v' denotes the conjugate transpose of v, and norm(u)
!  denotes the Euclidean norm. These reciprocal condition numbers always
!  lie between zero (very badly conditioned) and one (very well
!  conditioned). If n = 1, S(lambda) is defined to be 1.
!
!  An approximate error bound for a computed eigenvalue W(i) is given by
!
!                      EPS * norm(T) / S(i)
!
!  where EPS is the machine precision.
!
!  The reciprocal of the condition number of the right eigenvector u
!  corresponding to lambda is defined as follows. Suppose
!
!              T = ( lambda  c  )
!                  (   0    T22 )
!
!  Then the reciprocal condition number is
!
!          SEP( lambda, T22 ) = sigma-min( T22 - lambda*I )
!
!  where sigma-min denotes the smallest singular value. We approximate
!  the smallest singular value by the reciprocal of an estimate of the
!  one-norm of the inverse of T22 - lambda*I. If n = 1, SEP(1) is
!  defined to be abs(T(1,1)).
!
!  An approximate error bound for a computed right eigenvector VR(i)
!  is given by
!
!                      EPS * norm(T) / SEP(i)
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D0+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            SOMCON, WANTBH, WANTS, WANTSP
      CHARACTER          NORMIN
      INTEGER            I, IERR, IX, J, K, KASE, KS
      DOUBLE PRECISION   BIGNUM, EPS, EST, LNRM, RNRM, SCALE, SMLNUM, &
                         XNORM
      COMPLEX*16         CDUM, PROD
!     ..
!     .. Local Arrays ..
      COMPLEX*16         DUMMY( 1 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IZAMAX
      DOUBLE PRECISION   DLAMCH, DZNRM2
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, IZAMAX, DLAMCH, DZNRM2, ZDOTC
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLABAD, XERBLA, ZDRSCL, ZLACON, ZLACPY, ZLATRS, &
                         ZTREXC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTSP = LSAME( JOB, 'V' ) .OR. WANTBH
!
      SOMCON = LSAME( HOWMNY, 'S' )
!
!     Set M to the number of eigenpairs for which condition numbers are
!     to be computed.
!
      IF( SOMCON ) THEN
         M = 0
         DO 10 J = 1, N
            IF( SELECT( J ) ) &
               M = M + 1
   10    CONTINUE
      ELSE
         M = N
      END IF
!
      INFO = 0
      IF( .NOT.WANTS .AND. .NOT.WANTSP ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( HOWMNY, 'A' ) .AND. .NOT.SOMCON ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDVL.LT.1 .OR. ( WANTS .AND. LDVL.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LDVR.LT.1 .OR. ( WANTS .AND. LDVR.LT.N ) ) THEN
         INFO = -10
      ELSE IF( MM.LT.M ) THEN
         INFO = -13
      ELSE IF( LDWORK.LT.1 .OR. ( WANTSP .AND. LDWORK.LT.N ) ) THEN
         INFO = -16
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTRSNA', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( N.EQ.1 ) THEN
         IF( SOMCON ) THEN
            IF( .NOT.SELECT( 1 ) ) &
               RETURN
         END IF
         IF( WANTS ) &
            S( 1 ) = ONE
         IF( WANTSP ) &
            SEP( 1 ) = ABS( T( 1, 1 ) )
         RETURN
      END IF
!
!     Get machine constants
!
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
!
      KS = 1
      DO 50 K = 1, N
!
         IF( SOMCON ) THEN
            IF( .NOT.SELECT( K ) ) &
               GO TO 50
         END IF
!
         IF( WANTS ) THEN
!
!           Compute the reciprocal condition number of the k-th
!           eigenvalue.
!
            PROD = ZDOTC( N, VR( 1, KS ), 1, VL( 1, KS ), 1 )
            RNRM = DZNRM2( N, VR( 1, KS ), 1 )
            LNRM = DZNRM2( N, VL( 1, KS ), 1 )
            S( KS ) = ABS( PROD ) / ( RNRM*LNRM )
!
         END IF
!
         IF( WANTSP ) THEN
!
!           Estimate the reciprocal condition number of the k-th
!           eigenvector.
!
!           Copy the matrix T to the array WORK and swap the k-th
!           diagonal element to the (1,1) position.
!
            CALL ZLACPY( 'Full', N, N, T, LDT, WORK, LDWORK )
            CALL ZTREXC( 'No Q', N, WORK, LDWORK, DUMMY, 1, K, 1, IERR )
!
!           Form  C = T22 - lambda*I in WORK(2:N,2:N).
!
            DO 20 I = 2, N
               WORK( I, I ) = WORK( I, I ) - WORK( 1, 1 )
   20       CONTINUE
!
!           Estimate a lower bound for the 1-norm of inv(C'). The 1st
!           and (N+1)th columns of WORK are used to store work vectors.
!
            SEP( KS ) = ZERO
            EST = ZERO
            KASE = 0
            NORMIN = 'N'
   30       CONTINUE
            CALL ZLACON( N-1, WORK( 1, N+1 ), WORK, EST, KASE )
!
            IF( KASE.NE.0 ) THEN
               IF( KASE.EQ.1 ) THEN
!
!                 Solve C'*x = scale*b
!
                  CALL ZLATRS( 'Upper', 'Conjugate transpose', &
                               'Nonunit', NORMIN, N-1, WORK( 2, 2 ), &
                               LDWORK, WORK, SCALE, RWORK, IERR )
               ELSE
!
!                 Solve C*x = scale*b
!
                  CALL ZLATRS( 'Upper', 'No transpose', 'Nonunit', &
                               NORMIN, N-1, WORK( 2, 2 ), LDWORK, WORK, &
                               SCALE, RWORK, IERR )
               END IF
               NORMIN = 'Y'
               IF( SCALE.NE.ONE ) THEN
!
!                 Multiply by 1/SCALE if doing so will not cause
!                 overflow.
!
                  IX = IZAMAX( N-1, WORK, 1 )
                  XNORM = CABS1( WORK( IX, 1 ) )
                  IF( SCALE.LT.XNORM*SMLNUM .OR. SCALE.EQ.ZERO ) &
                     GO TO 40
                  CALL ZDRSCL( N, SCALE, WORK, 1 )
               END IF
               GO TO 30
            END IF
!
            SEP( KS ) = ONE / MAX( EST, SMLNUM )
         END IF
!
   40    CONTINUE
         KS = KS + 1
   50 CONTINUE
      RETURN
!
!     End of ZTRSNA
!
      END
      SUBROUTINE ZGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, &
                         INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            IHI, ILO, INFO, LDV, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   SCALE( * )
      COMPLEX*16         V( LDV, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGEBAK forms the right or left eigenvectors of a complex general
!  matrix by backward transformation on the computed eigenvectors of the
!  balanced matrix output by ZGEBAL.
!
!  Arguments
!  =========
!
!  JOB     (input) CHARACTER*1
!          Specifies the type of backward transformation required:
!          = 'N', do nothing, return immediately;
!          = 'P', do backward transformation for permutation only;
!          = 'S', do backward transformation for scaling only;
!          = 'B', do backward transformations for both permutation and
!                 scaling.
!          JOB must be the same as the argument JOB supplied to ZGEBAL.
!
!  SIDE    (input) CHARACTER*1
!          = 'R':  V contains right eigenvectors;
!          = 'L':  V contains left eigenvectors.
!
!  N       (input) INTEGER
!          The number of rows of the matrix V.  N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          The integers ILO and IHI determined by ZGEBAL.
!          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!
!  SCALE   (input) DOUBLE PRECISION array, dimension (N)
!          Details of the permutation and scaling factors, as returned
!          by ZGEBAL.
!
!  M       (input) INTEGER
!          The number of columns of the matrix V.  M >= 0.
!
!  V       (input/output) COMPLEX*16 array, dimension (LDV,M)
!          On entry, the matrix of right or left eigenvectors to be
!          transformed, as returned by ZHSEIN or ZTREVC.
!          On exit, V is overwritten by the transformed eigenvectors.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V. LDV >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV
      INTEGER            I, II, K
      DOUBLE PRECISION   S
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZDSCAL, ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Decode and Test the input parameters
!
      RIGHTV = LSAME( SIDE, 'R' )
      LEFTV = LSAME( SIDE, 'L' )
!
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. &
          .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -7
      ELSE IF( LDV.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEBAK', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
      IF( M.EQ.0 ) &
         RETURN
      IF( LSAME( JOB, 'N' ) ) &
         RETURN
!
      IF( ILO.EQ.IHI ) &
         GO TO 30
!
!     Backward balance
!
      IF( LSAME( JOB, 'S' ) .OR. LSAME( JOB, 'B' ) ) THEN
!
         IF( RIGHTV ) THEN
            DO 10 I = ILO, IHI
               S = SCALE( I )
               CALL ZDSCAL( M, S, V( I, 1 ), LDV )
   10       CONTINUE
         END IF
!
         IF( LEFTV ) THEN
            DO 20 I = ILO, IHI
               S = ONE / SCALE( I )
               CALL ZDSCAL( M, S, V( I, 1 ), LDV )
   20       CONTINUE
         END IF
!
      END IF
!
!     Backward permutation
!
!     For  I = ILO-1 step -1 until 1,
!              IHI+1 step 1 until N do --
!
   30 CONTINUE
      IF( LSAME( JOB, 'P' ) .OR. LSAME( JOB, 'B' ) ) THEN
         IF( RIGHTV ) THEN
            DO 40 II = 1, N
               I = II
               IF( I.GE.ILO .AND. I.LE.IHI ) &
                  GO TO 40
               IF( I.LT.ILO ) &
                  I = ILO - II
               K = SCALE( I )
               IF( K.EQ.I ) &
                  GO TO 40
               CALL ZSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   40       CONTINUE
         END IF
!
         IF( LEFTV ) THEN
            DO 50 II = 1, N
               I = II
               IF( I.GE.ILO .AND. I.LE.IHI ) &
                  GO TO 50
               IF( I.LT.ILO ) &
                  I = ILO - II
               K = SCALE( I )
               IF( K.EQ.I ) &
                  GO TO 50
               CALL ZSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   50       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of ZGEBAK
!
      END
      SUBROUTINE ZLAHRD( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            K, LDA, LDT, LDY, N, NB
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), T( LDT, NB ), TAU( NB ), &
                         Y( LDY, NB )
!     ..
!
!  Purpose
!  =======
!
!  ZLAHRD reduces the first NB columns of a complex general n-by-(n-k+1)
!  matrix A so that elements below the k-th subdiagonal are zero. The
!  reduction is performed by a unitary similarity transformation
!  Q' * A * Q. The routine returns the matrices V and T which determine
!  Q as a block reflector I - V*T*V', and also the matrix Y = A * V * T.
!
!  This is an auxiliary routine called by ZGEHRD.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.
!
!  K       (input) INTEGER
!          The offset for the reduction. Elements below the k-th
!          subdiagonal in the first NB columns are reduced to zero.
!
!  NB      (input) INTEGER
!          The number of columns to be reduced.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N-K+1)
!          On entry, the n-by-(n-k+1) general matrix A.
!          On exit, the elements on and above the k-th subdiagonal in
!          the first NB columns are overwritten with the corresponding
!          elements of the reduced matrix; the elements below the k-th
!          subdiagonal, with the array TAU, represent the matrix Q as a
!          product of elementary reflectors. The other columns of A are
!          unchanged. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  TAU     (output) COMPLEX*16 array, dimension (NB)
!          The scalar factors of the elementary reflectors. See Further
!          Details.
!
!  T       (output) COMPLEX*16 array, dimension (NB,NB)
!          The upper triangular matrix T.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T.  LDT >= NB.
!
!  Y       (output) COMPLEX*16 array, dimension (LDY,NB)
!          The n-by-nb matrix Y.
!
!  LDY     (input) INTEGER
!          The leading dimension of the array Y. LDY >= max(1,N).
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of nb elementary reflectors
!
!     Q = H(1) H(2) . . . H(nb).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a complex scalar, and v is a complex vector with
!  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in
!  A(i+k+1:n,i), and tau in TAU(i).
!
!  The elements of the vectors v together form the (n-k+1)-by-nb matrix
!  V which is needed, with T and Y, to apply the transformation to the
!  unreduced part of the matrix, using an update of the form:
!  A := (I - V*T*V') * (A - Y*V').
!
!  The contents of A on exit are illustrated by the following example
!  with n = 7, k = 3 and nb = 2:
!
!     ( a   h   a   a   a )
!     ( a   h   a   a   a )
!     ( a   h   a   a   a )
!     ( h   h   a   a   a )
!     ( v1  h   a   a   a )
!     ( v1  v2  a   a   a )
!     ( v1  v2  a   a   a )
!
!  where a denotes an element of the original matrix A, h denotes a
!  modified element of the upper Hessenberg matrix H, and vi denotes an
!  element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), &
                         ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      COMPLEX*16         EI
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZAXPY, ZCOPY, ZGEMV, ZLACGV, ZLARFG, ZSCAL, &
                         ZTRMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.LE.1 ) &
         RETURN
!
      DO 10 I = 1, NB
         IF( I.GT.1 ) THEN
!
!           Update A(1:n,i)
!
!           Compute i-th column of A - Y * V'
!
            CALL ZLACGV( I-1, A( K+I-1, 1 ), LDA )
            CALL ZGEMV( 'No transpose', N, I-1, -ONE, Y, LDY, &
                        A( K+I-1, 1 ), LDA, ONE, A( 1, I ), 1 )
            CALL ZLACGV( I-1, A( K+I-1, 1 ), LDA )
!
!           Apply I - V * T' * V' to this column (call it b) from the
!           left, using the last column of T as workspace
!
!           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
!                    ( V2 )             ( b2 )
!
!           where V1 is unit lower triangular
!
!           w := V1' * b1
!
            CALL ZCOPY( I-1, A( K+1, I ), 1, T( 1, NB ), 1 )
            CALL ZTRMV( 'Lower', 'Conjugate transpose', 'Unit', I-1, &
                        A( K+1, 1 ), LDA, T( 1, NB ), 1 )
!
!           w := w + V2'*b2
!
            CALL ZGEMV( 'Conjugate transpose', N-K-I+1, I-1, ONE, &
                        A( K+I, 1 ), LDA, A( K+I, I ), 1, ONE, &
                        T( 1, NB ), 1 )
!
!           w := T'*w
!
            CALL ZTRMV( 'Upper', 'Conjugate transpose', 'Non-unit', I-1, &
                        T, LDT, T( 1, NB ), 1 )
!
!           b2 := b2 - V2*w
!
            CALL ZGEMV( 'No transpose', N-K-I+1, I-1, -ONE, A( K+I, 1 ), &
                        LDA, T( 1, NB ), 1, ONE, A( K+I, I ), 1 )
!
!           b1 := b1 - V1*w
!
            CALL ZTRMV( 'Lower', 'No transpose', 'Unit', I-1, &
                        A( K+1, 1 ), LDA, T( 1, NB ), 1 )
            CALL ZAXPY( I-1, -ONE, T( 1, NB ), 1, A( K+1, I ), 1 )
!
            A( K+I-1, I-1 ) = EI
         END IF
!
!        Generate the elementary reflector H(i) to annihilate
!        A(k+i+1:n,i)
!
         EI = A( K+I, I )
         CALL ZLARFG( N-K-I+1, EI, A( MIN( K+I+1, N ), I ), 1, &
                      TAU( I ) )
         A( K+I, I ) = ONE
!
!        Compute  Y(1:n,i)
!
         CALL ZGEMV( 'No transpose', N, N-K-I+1, ONE, A( 1, I+1 ), LDA, &
                     A( K+I, I ), 1, ZERO, Y( 1, I ), 1 )
         CALL ZGEMV( 'Conjugate transpose', N-K-I+1, I-1, ONE, &
                     A( K+I, 1 ), LDA, A( K+I, I ), 1, ZERO, T( 1, I ), &
                     1 )
         CALL ZGEMV( 'No transpose', N, I-1, -ONE, Y, LDY, T( 1, I ), 1, &
                     ONE, Y( 1, I ), 1 )
         CALL ZSCAL( N, TAU( I ), Y( 1, I ), 1 )
!
!        Compute T(1:i,i)
!
         CALL ZSCAL( I-1, -TAU( I ), T( 1, I ), 1 )
         CALL ZTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, LDT, &
                     T( 1, I ), 1 )
         T( I, I ) = TAU( I )
!
   10 CONTINUE
      A( K+NB, NB ) = EI
!
      RETURN
!
!     End of ZLAHRD
!
      END
      SUBROUTINE ZGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZGEHD2 reduces a complex general matrix A to upper Hessenberg form H
!  by a unitary similarity transformation:  Q' * A * Q = H .
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          It is assumed that A is already upper triangular in rows
!          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!          set by a previous call to ZGEBAL; otherwise they should be
!          set to 1 and N respectively. See Further Details.
!          1 <= ILO <= IHI <= max(1,N).
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the n by n general matrix to be reduced.
!          On exit, the upper triangle and the first subdiagonal of A
!          are overwritten with the upper Hessenberg matrix H, and the
!          elements below the first subdiagonal, with the array TAU,
!          represent the unitary matrix Q as a product of elementary
!          reflectors. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  TAU     (output) COMPLEX*16 array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace) COMPLEX*16 array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of (ihi-ilo) elementary
!  reflectors
!
!     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a complex scalar, and v is a complex vector with
!  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
!  exit in A(i+2:ihi,i), and tau in TAU(i).
!
!  The contents of A are illustrated by the following example, with
!  n = 7, ilo = 2 and ihi = 6:
!
!  on entry,                        on exit,
!
!  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
!  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
!  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
!  (                         a )    (                          a )
!
!  where a denotes an element of the original matrix A, h denotes a
!  modified element of the upper Hessenberg matrix H, and vi denotes an
!  element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      COMPLEX*16         ALPHA
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF, ZLARFG
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEHD2', -INFO )
         RETURN
      END IF
!
      DO 10 I = ILO, IHI - 1
!
!        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
!
         ALPHA = A( I+1, I )
         CALL ZLARFG( IHI-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAU( I ) )
         A( I+1, I ) = ONE
!
!        Apply H(i) to A(1:ihi,i+1:ihi) from the right
!
         CALL ZLARF( 'Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ), &
                     A( 1, I+1 ), LDA, WORK )
!
!        Apply H(i)' to A(i+1:ihi,i+1:n) from the left
!
         CALL ZLARF( 'Left', IHI-I, N-I, A( I+1, I ), 1, &
                     DCONJG( TAU( I ) ), A( I+1, I+1 ), LDA, WORK )
!
         A( I+1, I ) = ALPHA
   10 CONTINUE
!
      RETURN
!
!     End of ZGEHD2
!
      END
      SUBROUTINE ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, &
                         IHIZ, Z, LDZ, INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      LOGICAL            WANTT, WANTZ
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), W( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  ZLAHQR is an auxiliary routine called by ZHSEQR to update the
!  eigenvalues and Schur decomposition already computed by ZHSEQR, by
!  dealing with the Hessenberg submatrix in rows and columns ILO to IHI.
!
!  Arguments
!  =========
!
!  WANTT   (input) LOGICAL
!          = .TRUE. : the full Schur form T is required;
!          = .FALSE.: only eigenvalues are required.
!
!  WANTZ   (input) LOGICAL
!          = .TRUE. : the matrix of Schur vectors Z is required;
!          = .FALSE.: Schur vectors are not required.
!
!  N       (input) INTEGER
!          The order of the matrix H.  N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          It is assumed that H is already upper triangular in rows and
!          columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1).
!          ZLAHQR works primarily with the Hessenberg submatrix in rows
!          and columns ILO to IHI, but applies transformations to all of
!          H if WANTT is .TRUE..
!          1 <= ILO <= max(1,IHI); IHI <= N.
!
!  H       (input/output) COMPLEX*16 array, dimension (LDH,N)
!          On entry, the upper Hessenberg matrix H.
!          On exit, if WANTT is .TRUE., H is upper triangular in rows
!          and columns ILO:IHI, with any 2-by-2 diagonal blocks in
!          standard form. If WANTT is .FALSE., the contents of H are
!          unspecified on exit.
!
!  LDH     (input) INTEGER
!          The leading dimension of the array H. LDH >= max(1,N).
!
!  W       (output) COMPLEX*16 array, dimension (N)
!          The computed eigenvalues ILO to IHI are stored in the
!          corresponding elements of W. If WANTT is .TRUE., the
!          eigenvalues are stored in the same order as on the diagonal
!          of the Schur form returned in H, with W(i) = H(i,i).
!
!  ILOZ    (input) INTEGER
!  IHIZ    (input) INTEGER
!          Specify the rows of Z to which transformations must be
!          applied if WANTZ is .TRUE..
!          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
!
!  Z       (input/output) COMPLEX*16 array, dimension (LDZ,N)
!          If WANTZ is .TRUE., on entry Z must contain the current
!          matrix Z of transformations accumulated by ZHSEQR, and on
!          exit Z has been updated; transformations are applied only to
!          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
!          If WANTZ is .FALSE., Z is not referenced.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z. LDZ >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          > 0: if INFO = i, ZLAHQR failed to compute all the
!               eigenvalues ILO to IHI in a total of 30*(IHI-ILO+1)
!               iterations; elements i+1:ihi of W contain those
!               eigenvalues which have been successfully computed.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), &
                         ONE = ( 1.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION   RZERO, RONE, HALF
      PARAMETER          ( RZERO = 0.0D+0, RONE = 1.0D+0, &
                         HALF = 0.5D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I1, I2, ITN, ITS, J, K, L, M, NH, NZ
      DOUBLE PRECISION   H10, H21, OVFL, RTEMP, S, SMLNUM, T2, TST1, &
                         ULP, UNFL
      COMPLEX*16         CDUM, H11, H11S, H22, SUM, T, T1, TEMP, U, V2, &
                         X, Y
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   RWORK( 1 )
      COMPLEX*16         V( 2 )
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2, ZLANHS
      COMPLEX*16         ZLADIV
      EXTERNAL           DLAMCH, DLAPY2, ZLANHS, ZLADIV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLABAD, ZCOPY, ZLARFG, ZSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, MAX, MIN, SQRT
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
      IF( ILO.EQ.IHI ) THEN
         W( ILO ) = H( ILO, ILO )
         RETURN
      END IF
!
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
!
!     Set machine-dependent constants for the stopping criterion.
!     If norm(H) <= sqrt(OVFL), overflow should not occur.
!
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = RONE / UNFL
      CALL DLABAD( UNFL, OVFL )
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( NH / ULP )
!
!     I1 and I2 are the indices of the first row and last column of H
!     to which transformations must be applied. If eigenvalues only are
!     being computed, I1 and I2 are set inside the main loop.
!
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
!
!     ITN is the total number of QR iterations allowed.
!
      ITN = 30*NH
!
!     The main loop begins here. I is the loop index and decreases from
!     IHI to ILO in steps of 1. Each iteration of the loop works
!     with the active submatrix in rows and columns L to I.
!     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
!     H(L,L-1) is negligible so that the matrix splits.
!
      I = IHI
   10 CONTINUE
      IF( I.LT.ILO ) &
         GO TO 130
!
!     Perform QR iterations on rows and columns ILO to I until a
!     submatrix of order 1 splits off at the bottom because a
!     subdiagonal element has become negligible.
!
      L = ILO
      DO 110 ITS = 0, ITN
!
!        Look for a single small subdiagonal element.
!
         DO 20 K = I, L + 1, -1
            TST1 = CABS1( H( K-1, K-1 ) ) + CABS1( H( K, K ) )
            IF( TST1.EQ.RZERO ) &
               TST1 = ZLANHS( '1', I-L+1, H( L, L ), LDH, RWORK )
            IF( ABS( DBLE( H( K, K-1 ) ) ).LE.MAX( ULP*TST1, SMLNUM ) ) &
               GO TO 30
   20    CONTINUE
   30    CONTINUE
         L = K
         IF( L.GT.ILO ) THEN
!
!           H(L,L-1) is negligible
!
            H( L, L-1 ) = ZERO
         END IF
!
!        Exit from loop if a submatrix of order 1 has split off.
!
         IF( L.GE.I ) &
            GO TO 120
!
!        Now the active submatrix is in rows and columns L to I. If
!        eigenvalues only are being computed, only the active submatrix
!        need be transformed.
!
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
!
         IF( ITS.EQ.10 .OR. ITS.EQ.20 ) THEN
!
!           Exceptional shift.
!
            T = ABS( DBLE( H( I, I-1 ) ) ) + &
                ABS( DBLE( H( I-1, I-2 ) ) )
         ELSE
!
!           Wilkinson's shift.
!
            T = H( I, I )
            U = H( I-1, I )*DBLE( H( I, I-1 ) )
            IF( U.NE.ZERO ) THEN
               X = HALF*( H( I-1, I-1 )-T )
               Y = SQRT( X*X+U )
               IF( DBLE( X )*DBLE( Y )+DIMAG( X )*DIMAG( Y ).LT.RZERO ) &
                  Y = -Y
               T = T - ZLADIV( U, ( X+Y ) )
            END IF
         END IF
!
!        Look for two consecutive small subdiagonal elements.
!
         DO 40 M = I - 1, L, -1
!
!           Determine the effect of starting the single-shift QR
!           iteration at row M, and see if this would make H(M,M-1)
!           negligible.
!
            H11 = H( M, M )
            H22 = H( M+1, M+1 )
            H11S = H11 - T
            H21 = H( M+1, M )
            S = CABS1( H11S ) + ABS( H21 )
            H11S = H11S / S
            H21 = H21 / S
            V( 1 ) = H11S
            V( 2 ) = H21
            IF( M.EQ.L ) &
               GO TO 50
            H10 = H( M, M-1 )
            TST1 = CABS1( H11S )*( CABS1( H11 )+CABS1( H22 ) )
            IF( ABS( H10*H21 ).LE.ULP*TST1 ) &
               GO TO 50
   40    CONTINUE
   50    CONTINUE
!
!        Single-shift QR step
!
         DO 100 K = M, I - 1
!
!           The first iteration of this loop determines a reflection G
!           from the vector V and applies it from left and right to H,
!           thus creating a nonzero bulge below the subdiagonal.
!
!           Each subsequent iteration determines a reflection G to
!           restore the Hessenberg form in the (K-1)th column, and thus
!           chases the bulge one step toward the bottom of the active
!           submatrix.
!
!           V(2) is always real before the call to ZLARFG, and hence
!           after the call T2 ( = T1*V(2) ) is also real.
!
            IF( K.GT.M ) &
               CALL ZCOPY( 2, H( K, K-1 ), 1, V, 1 )
            CALL ZLARFG( 2, V( 1 ), V( 2 ), 1, T1 )
            IF( K.GT.M ) THEN
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
            END IF
            V2 = V( 2 )
            T2 = DBLE( T1*V2 )
!
!           Apply G from the left to transform the rows of the matrix
!           in columns K to I2.
!
            DO 60 J = K, I2
               SUM = DCONJG( T1 )*H( K, J ) + T2*H( K+1, J )
               H( K, J ) = H( K, J ) - SUM
               H( K+1, J ) = H( K+1, J ) - SUM*V2
   60       CONTINUE
!
!           Apply G from the right to transform the columns of the
!           matrix in rows I1 to min(K+2,I).
!
            DO 70 J = I1, MIN( K+2, I )
               SUM = T1*H( J, K ) + T2*H( J, K+1 )
               H( J, K ) = H( J, K ) - SUM
               H( J, K+1 ) = H( J, K+1 ) - SUM*DCONJG( V2 )
   70       CONTINUE
!
            IF( WANTZ ) THEN
!
!              Accumulate transformations in the matrix Z
!
               DO 80 J = ILOZ, IHIZ
                  SUM = T1*Z( J, K ) + T2*Z( J, K+1 )
                  Z( J, K ) = Z( J, K ) - SUM
                  Z( J, K+1 ) = Z( J, K+1 ) - SUM*DCONJG( V2 )
   80          CONTINUE
            END IF
!
            IF( K.EQ.M .AND. M.GT.L ) THEN
!
!              If the QR step was started at row M > L because two
!              consecutive small subdiagonals were found, then extra
!              scaling must be performed to ensure that H(M,M-1) remains
!              real.
!
               TEMP = ONE - T1
               TEMP = TEMP / DLAPY2( DBLE( TEMP ), DIMAG( TEMP ) )
               H( M+1, M ) = H( M+1, M )*DCONJG( TEMP )
               IF( M+2.LE.I ) &
                  H( M+2, M+1 ) = H( M+2, M+1 )*TEMP
               DO 90 J = M, I
                  IF( J.NE.M+1 ) THEN
                     IF( I2.GT.J ) &
                        CALL ZSCAL( I2-J, TEMP, H( J, J+1 ), LDH )
                     CALL ZSCAL( J-I1, DCONJG( TEMP ), H( I1, J ), 1 )
                     IF( WANTZ ) THEN
                        CALL ZSCAL( NZ, DCONJG( TEMP ), Z( ILOZ, J ), &
                                    1 )
                     END IF
                  END IF
   90          CONTINUE
            END IF
  100    CONTINUE
!
!        Ensure that H(I,I-1) is real.
!
         TEMP = H( I, I-1 )
         IF( DIMAG( TEMP ).NE.RZERO ) THEN
            RTEMP = DLAPY2( DBLE( TEMP ), DIMAG( TEMP ) )
            H( I, I-1 ) = RTEMP
            TEMP = TEMP / RTEMP
            IF( I2.GT.I ) &
               CALL ZSCAL( I2-I, DCONJG( TEMP ), H( I, I+1 ), LDH )
            CALL ZSCAL( I-I1, TEMP, H( I1, I ), 1 )
            IF( WANTZ ) THEN
               CALL ZSCAL( NZ, TEMP, Z( ILOZ, I ), 1 )
            END IF
         END IF
!
  110 CONTINUE
!
!     Failure to converge in remaining number of iterations
!
      INFO = I
      RETURN
!
  120 CONTINUE
!
!     H(I,I-1) is negligible: one eigenvalue has converged.
!
      W( I ) = H( I, I )
!
!     Decrement number of remaining iterations, and return to start of
!     the main loop with new value of I.
!
      ITN = ITN - ITS
      I = L - 1
      GO TO 10
!
  130 CONTINUE
      RETURN
!
!     End of ZLAHQR
!
      END
      DOUBLE PRECISION FUNCTION ZLANHS( NORM, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   WORK( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZLANHS  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  Hessenberg matrix A.
!
!  Description
!  ===========
!
!  ZLANHS returns the value
!
!     ZLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in ZLANHS as described
!          above.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, ZLANHS is
!          set to zero.
!
!  A       (input) COMPLEX*16 array, dimension (LDA,N)
!          The n by n upper Hessenberg matrix A; the part of A below the
!          first sub-diagonal is not referenced.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(N,1).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
!          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
!          referenced.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   SCALE, SUM, VALUE
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZLASSQ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, MIN( N, J+1 )
               VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
!
!        Find norm1(A).
!
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, MIN( N, J+1 )
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            VALUE = MAX( VALUE, SUM )
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
         DO 50 I = 1, N
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, MIN( N, J+1 )
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, N
            VALUE = MAX( VALUE, WORK( I ) )
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL ZLASSQ( MIN( N, J+1 ), A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
!
      ZLANHS = VALUE
      RETURN
!
!     End of ZLANHS
!
      END
      SUBROUTINE ZLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            LDC, M, N
      COMPLEX*16         TAU
!     ..
!     .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZLARFX applies a complex elementary reflector H to a complex m by n
!  matrix C, from either the left or the right. H is represented in the
!  form
!
!        H = I - tau * v * v'
!
!  where tau is a complex scalar and v is a complex vector.
!
!  If tau = 0, then H is taken to be the unit matrix
!
!  This version uses inline code if H has order < 11.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': form  H * C
!          = 'R': form  C * H
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  V       (input) COMPLEX*16 array, dimension (M) if SIDE = 'L'
!                                        or (N) if SIDE = 'R'
!          The vector v in the representation of H.
!
!  TAU     (input) COMPLEX*16
!          The value tau in the representation of H.
!
!  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDA >= max(1,M).
!
!  WORK    (workspace) COMPLEX*16 array, dimension (N) if SIDE = 'L'
!                                            or (M) if SIDE = 'R'
!          WORK is not referenced if H has order < 11.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), &
                         ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            J
      COMPLEX*16         SUM, T1, T10, T2, T3, T4, T5, T6, T7, T8, T9, &
                         V1, V10, V2, V3, V4, V5, V6, V7, V8, V9
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZGEMV, ZGERC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
!     ..
!     .. Executable Statements ..
!
      IF( TAU.EQ.ZERO ) &
         RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  H * C, where H has order m.
!
         GO TO ( 10, 30, 50, 70, 90, 110, 130, 150, &
                 170, 190 )M
!
!        Code for general M
!
!        w := C'*v
!
         CALL ZGEMV( 'Conjugate transpose', M, N, ONE, C, LDC, V, 1, &
                     ZERO, WORK, 1 )
!
!        C := C - tau * v * w'
!
         CALL ZGERC( M, N, -TAU, V, 1, WORK, 1, C, LDC )
         GO TO 410
   10    CONTINUE
!
!        Special code for 1 x 1 Householder
!
         T1 = ONE - TAU*V( 1 )*DCONJG( V( 1 ) )
         DO 20 J = 1, N
            C( 1, J ) = T1*C( 1, J )
   20    CONTINUE
         GO TO 410
   30    CONTINUE
!
!        Special code for 2 x 2 Householder
!
         V1 = DCONJG( V( 1 ) )
         T1 = TAU*DCONJG( V1 )
         V2 = DCONJG( V( 2 ) )
         T2 = TAU*DCONJG( V2 )
         DO 40 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
   40    CONTINUE
         GO TO 410
   50    CONTINUE
!
!        Special code for 3 x 3 Householder
!
         V1 = DCONJG( V( 1 ) )
         T1 = TAU*DCONJG( V1 )
         V2 = DCONJG( V( 2 ) )
         T2 = TAU*DCONJG( V2 )
         V3 = DCONJG( V( 3 ) )
         T3 = TAU*DCONJG( V3 )
         DO 60 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
   60    CONTINUE
         GO TO 410
   70    CONTINUE
!
!        Special code for 4 x 4 Householder
!
         V1 = DCONJG( V( 1 ) )
         T1 = TAU*DCONJG( V1 )
         V2 = DCONJG( V( 2 ) )
         T2 = TAU*DCONJG( V2 )
         V3 = DCONJG( V( 3 ) )
         T3 = TAU*DCONJG( V3 )
         V4 = DCONJG( V( 4 ) )
         T4 = TAU*DCONJG( V4 )
         DO 80 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
   80    CONTINUE
         GO TO 410
   90    CONTINUE
!
!        Special code for 5 x 5 Householder
!
         V1 = DCONJG( V( 1 ) )
         T1 = TAU*DCONJG( V1 )
         V2 = DCONJG( V( 2 ) )
         T2 = TAU*DCONJG( V2 )
         V3 = DCONJG( V( 3 ) )
         T3 = TAU*DCONJG( V3 )
         V4 = DCONJG( V( 4 ) )
         T4 = TAU*DCONJG( V4 )
         V5 = DCONJG( V( 5 ) )
         T5 = TAU*DCONJG( V5 )
         DO 100 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
  100    CONTINUE
         GO TO 410
  110    CONTINUE
!
!        Special code for 6 x 6 Householder
!
         V1 = DCONJG( V( 1 ) )
         T1 = TAU*DCONJG( V1 )
         V2 = DCONJG( V( 2 ) )
         T2 = TAU*DCONJG( V2 )
         V3 = DCONJG( V( 3 ) )
         T3 = TAU*DCONJG( V3 )
         V4 = DCONJG( V( 4 ) )
         T4 = TAU*DCONJG( V4 )
         V5 = DCONJG( V( 5 ) )
         T5 = TAU*DCONJG( V5 )
         V6 = DCONJG( V( 6 ) )
         T6 = TAU*DCONJG( V6 )
         DO 120 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
  120    CONTINUE
         GO TO 410
  130    CONTINUE
!
!        Special code for 7 x 7 Householder
!
         V1 = DCONJG( V( 1 ) )
         T1 = TAU*DCONJG( V1 )
         V2 = DCONJG( V( 2 ) )
         T2 = TAU*DCONJG( V2 )
         V3 = DCONJG( V( 3 ) )
         T3 = TAU*DCONJG( V3 )
         V4 = DCONJG( V( 4 ) )
         T4 = TAU*DCONJG( V4 )
         V5 = DCONJG( V( 5 ) )
         T5 = TAU*DCONJG( V5 )
         V6 = DCONJG( V( 6 ) )
         T6 = TAU*DCONJG( V6 )
         V7 = DCONJG( V( 7 ) )
         T7 = TAU*DCONJG( V7 )
         DO 140 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                  V7*C( 7, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
  140    CONTINUE
         GO TO 410
  150    CONTINUE
!
!        Special code for 8 x 8 Householder
!
         V1 = DCONJG( V( 1 ) )
         T1 = TAU*DCONJG( V1 )
         V2 = DCONJG( V( 2 ) )
         T2 = TAU*DCONJG( V2 )
         V3 = DCONJG( V( 3 ) )
         T3 = TAU*DCONJG( V3 )
         V4 = DCONJG( V( 4 ) )
         T4 = TAU*DCONJG( V4 )
         V5 = DCONJG( V( 5 ) )
         T5 = TAU*DCONJG( V5 )
         V6 = DCONJG( V( 6 ) )
         T6 = TAU*DCONJG( V6 )
         V7 = DCONJG( V( 7 ) )
         T7 = TAU*DCONJG( V7 )
         V8 = DCONJG( V( 8 ) )
         T8 = TAU*DCONJG( V8 )
         DO 160 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                  V7*C( 7, J ) + V8*C( 8, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
            C( 8, J ) = C( 8, J ) - SUM*T8
  160    CONTINUE
         GO TO 410
  170    CONTINUE
!
!        Special code for 9 x 9 Householder
!
         V1 = DCONJG( V( 1 ) )
         T1 = TAU*DCONJG( V1 )
         V2 = DCONJG( V( 2 ) )
         T2 = TAU*DCONJG( V2 )
         V3 = DCONJG( V( 3 ) )
         T3 = TAU*DCONJG( V3 )
         V4 = DCONJG( V( 4 ) )
         T4 = TAU*DCONJG( V4 )
         V5 = DCONJG( V( 5 ) )
         T5 = TAU*DCONJG( V5 )
         V6 = DCONJG( V( 6 ) )
         T6 = TAU*DCONJG( V6 )
         V7 = DCONJG( V( 7 ) )
         T7 = TAU*DCONJG( V7 )
         V8 = DCONJG( V( 8 ) )
         T8 = TAU*DCONJG( V8 )
         V9 = DCONJG( V( 9 ) )
         T9 = TAU*DCONJG( V9 )
         DO 180 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                  V7*C( 7, J ) + V8*C( 8, J ) + V9*C( 9, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
            C( 8, J ) = C( 8, J ) - SUM*T8
            C( 9, J ) = C( 9, J ) - SUM*T9
  180    CONTINUE
         GO TO 410
  190    CONTINUE
!
!        Special code for 10 x 10 Householder
!
         V1 = DCONJG( V( 1 ) )
         T1 = TAU*DCONJG( V1 )
         V2 = DCONJG( V( 2 ) )
         T2 = TAU*DCONJG( V2 )
         V3 = DCONJG( V( 3 ) )
         T3 = TAU*DCONJG( V3 )
         V4 = DCONJG( V( 4 ) )
         T4 = TAU*DCONJG( V4 )
         V5 = DCONJG( V( 5 ) )
         T5 = TAU*DCONJG( V5 )
         V6 = DCONJG( V( 6 ) )
         T6 = TAU*DCONJG( V6 )
         V7 = DCONJG( V( 7 ) )
         T7 = TAU*DCONJG( V7 )
         V8 = DCONJG( V( 8 ) )
         T8 = TAU*DCONJG( V8 )
         V9 = DCONJG( V( 9 ) )
         T9 = TAU*DCONJG( V9 )
         V10 = DCONJG( V( 10 ) )
         T10 = TAU*DCONJG( V10 )
         DO 200 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                  V7*C( 7, J ) + V8*C( 8, J ) + V9*C( 9, J ) + &
                  V10*C( 10, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
            C( 8, J ) = C( 8, J ) - SUM*T8
            C( 9, J ) = C( 9, J ) - SUM*T9
            C( 10, J ) = C( 10, J ) - SUM*T10
  200    CONTINUE
         GO TO 410
      ELSE
!
!        Form  C * H, where H has order n.
!
         GO TO ( 210, 230, 250, 270, 290, 310, 330, 350, &
                 370, 390 )N
!
!        Code for general N
!
!        w := C * v
!
         CALL ZGEMV( 'No transpose', M, N, ONE, C, LDC, V, 1, ZERO, &
                     WORK, 1 )
!
!        C := C - tau * w * v'
!
         CALL ZGERC( M, N, -TAU, WORK, 1, V, 1, C, LDC )
         GO TO 410
  210    CONTINUE
!
!        Special code for 1 x 1 Householder
!
         T1 = ONE - TAU*V( 1 )*DCONJG( V( 1 ) )
         DO 220 J = 1, M
            C( J, 1 ) = T1*C( J, 1 )
  220    CONTINUE
         GO TO 410
  230    CONTINUE
!
!        Special code for 2 x 2 Householder
!
         V1 = V( 1 )
         T1 = TAU*DCONJG( V1 )
         V2 = V( 2 )
         T2 = TAU*DCONJG( V2 )
         DO 240 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
  240    CONTINUE
         GO TO 410
  250    CONTINUE
!
!        Special code for 3 x 3 Householder
!
         V1 = V( 1 )
         T1 = TAU*DCONJG( V1 )
         V2 = V( 2 )
         T2 = TAU*DCONJG( V2 )
         V3 = V( 3 )
         T3 = TAU*DCONJG( V3 )
         DO 260 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
  260    CONTINUE
         GO TO 410
  270    CONTINUE
!
!        Special code for 4 x 4 Householder
!
         V1 = V( 1 )
         T1 = TAU*DCONJG( V1 )
         V2 = V( 2 )
         T2 = TAU*DCONJG( V2 )
         V3 = V( 3 )
         T3 = TAU*DCONJG( V3 )
         V4 = V( 4 )
         T4 = TAU*DCONJG( V4 )
         DO 280 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
  280    CONTINUE
         GO TO 410
  290    CONTINUE
!
!        Special code for 5 x 5 Householder
!
         V1 = V( 1 )
         T1 = TAU*DCONJG( V1 )
         V2 = V( 2 )
         T2 = TAU*DCONJG( V2 )
         V3 = V( 3 )
         T3 = TAU*DCONJG( V3 )
         V4 = V( 4 )
         T4 = TAU*DCONJG( V4 )
         V5 = V( 5 )
         T5 = TAU*DCONJG( V5 )
         DO 300 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
  300    CONTINUE
         GO TO 410
  310    CONTINUE
!
!        Special code for 6 x 6 Householder
!
         V1 = V( 1 )
         T1 = TAU*DCONJG( V1 )
         V2 = V( 2 )
         T2 = TAU*DCONJG( V2 )
         V3 = V( 3 )
         T3 = TAU*DCONJG( V3 )
         V4 = V( 4 )
         T4 = TAU*DCONJG( V4 )
         V5 = V( 5 )
         T5 = TAU*DCONJG( V5 )
         V6 = V( 6 )
         T6 = TAU*DCONJG( V6 )
         DO 320 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
  320    CONTINUE
         GO TO 410
  330    CONTINUE
!
!        Special code for 7 x 7 Householder
!
         V1 = V( 1 )
         T1 = TAU*DCONJG( V1 )
         V2 = V( 2 )
         T2 = TAU*DCONJG( V2 )
         V3 = V( 3 )
         T3 = TAU*DCONJG( V3 )
         V4 = V( 4 )
         T4 = TAU*DCONJG( V4 )
         V5 = V( 5 )
         T5 = TAU*DCONJG( V5 )
         V6 = V( 6 )
         T6 = TAU*DCONJG( V6 )
         V7 = V( 7 )
         T7 = TAU*DCONJG( V7 )
         DO 340 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                  V7*C( J, 7 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
  340    CONTINUE
         GO TO 410
  350    CONTINUE
!
!        Special code for 8 x 8 Householder
!
         V1 = V( 1 )
         T1 = TAU*DCONJG( V1 )
         V2 = V( 2 )
         T2 = TAU*DCONJG( V2 )
         V3 = V( 3 )
         T3 = TAU*DCONJG( V3 )
         V4 = V( 4 )
         T4 = TAU*DCONJG( V4 )
         V5 = V( 5 )
         T5 = TAU*DCONJG( V5 )
         V6 = V( 6 )
         T6 = TAU*DCONJG( V6 )
         V7 = V( 7 )
         T7 = TAU*DCONJG( V7 )
         V8 = V( 8 )
         T8 = TAU*DCONJG( V8 )
         DO 360 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                  V7*C( J, 7 ) + V8*C( J, 8 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
            C( J, 8 ) = C( J, 8 ) - SUM*T8
  360    CONTINUE
         GO TO 410
  370    CONTINUE
!
!        Special code for 9 x 9 Householder
!
         V1 = V( 1 )
         T1 = TAU*DCONJG( V1 )
         V2 = V( 2 )
         T2 = TAU*DCONJG( V2 )
         V3 = V( 3 )
         T3 = TAU*DCONJG( V3 )
         V4 = V( 4 )
         T4 = TAU*DCONJG( V4 )
         V5 = V( 5 )
         T5 = TAU*DCONJG( V5 )
         V6 = V( 6 )
         T6 = TAU*DCONJG( V6 )
         V7 = V( 7 )
         T7 = TAU*DCONJG( V7 )
         V8 = V( 8 )
         T8 = TAU*DCONJG( V8 )
         V9 = V( 9 )
         T9 = TAU*DCONJG( V9 )
         DO 380 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                  V7*C( J, 7 ) + V8*C( J, 8 ) + V9*C( J, 9 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
            C( J, 8 ) = C( J, 8 ) - SUM*T8
            C( J, 9 ) = C( J, 9 ) - SUM*T9
  380    CONTINUE
         GO TO 410
  390    CONTINUE
!
!        Special code for 10 x 10 Householder
!
         V1 = V( 1 )
         T1 = TAU*DCONJG( V1 )
         V2 = V( 2 )
         T2 = TAU*DCONJG( V2 )
         V3 = V( 3 )
         T3 = TAU*DCONJG( V3 )
         V4 = V( 4 )
         T4 = TAU*DCONJG( V4 )
         V5 = V( 5 )
         T5 = TAU*DCONJG( V5 )
         V6 = V( 6 )
         T6 = TAU*DCONJG( V6 )
         V7 = V( 7 )
         T7 = TAU*DCONJG( V7 )
         V8 = V( 8 )
         T8 = TAU*DCONJG( V8 )
         V9 = V( 9 )
         T9 = TAU*DCONJG( V9 )
         V10 = V( 10 )
         T10 = TAU*DCONJG( V10 )
         DO 400 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                  V7*C( J, 7 ) + V8*C( J, 8 ) + V9*C( J, 9 ) + &
                  V10*C( J, 10 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
            C( J, 8 ) = C( J, 8 ) - SUM*T8
            C( J, 9 ) = C( J, 9 ) - SUM*T9
            C( J, 10 ) = C( J, 10 ) - SUM*T10
  400    CONTINUE
         GO TO 410
      END IF
  410 CONTINUE
      RETURN
!
!     End of ZLARFX
!
      END
      SUBROUTINE ZLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE, &
                         CNORM, INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          DIAG, NORMIN, TRANS, UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   SCALE
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   CNORM( * )
      COMPLEX*16         A( LDA, * ), X( * )
!     ..
!
!  Purpose
!  =======
!
!  ZLATRS solves one of the triangular systems
!
!     A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
!
!  with scaling to prevent overflow.  Here A is an upper or lower
!  triangular matrix, A**T denotes the transpose of A, A**H denotes the
!  conjugate transpose of A, x and b are n-element vectors, and s is a
!  scaling factor, usually less than or equal to 1, chosen so that the
!  components of x will be less than the overflow threshold.  If the
!  unscaled problem will not cause overflow, the Level 2 BLAS routine
!  ZTRSV is called. If the matrix A is singular (A(j,j) = 0 for some j),
!  then s is set to 0 and a non-trivial solution to A*x = 0 is returned.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the matrix A is upper or lower triangular.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  TRANS   (input) CHARACTER*1
!          Specifies the operation applied to A.
!          = 'N':  Solve A * x = s*b     (No transpose)
!          = 'T':  Solve A**T * x = s*b  (Transpose)
!          = 'C':  Solve A**H * x = s*b  (Conjugate transpose)
!
!  DIAG    (input) CHARACTER*1
!          Specifies whether or not the matrix A is unit triangular.
!          = 'N':  Non-unit triangular
!          = 'U':  Unit triangular
!
!  NORMIN  (input) CHARACTER*1
!          Specifies whether CNORM has been set or not.
!          = 'Y':  CNORM contains the column norms on entry
!          = 'N':  CNORM is not set on entry.  On exit, the norms will
!                  be computed and stored in CNORM.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input) COMPLEX*16 array, dimension (LDA,N)
!          The triangular matrix A.  If UPLO = 'U', the leading n by n
!          upper triangular part of the array A contains the upper
!          triangular matrix, and the strictly lower triangular part of
!          A is not referenced.  If UPLO = 'L', the leading n by n lower
!          triangular part of the array A contains the lower triangular
!          matrix, and the strictly upper triangular part of A is not
!          referenced.  If DIAG = 'U', the diagonal elements of A are
!          also not referenced and are assumed to be 1.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max (1,N).
!
!  X       (input/output) COMPLEX*16 array, dimension (N)
!          On entry, the right hand side b of the triangular system.
!          On exit, X is overwritten by the solution vector x.
!
!  SCALE   (output) DOUBLE PRECISION
!          The scaling factor s for the triangular system
!             A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b.
!          If SCALE = 0, the matrix A is singular or badly scaled, and
!          the vector x is an exact or approximate solution to A*x = 0.
!
!  CNORM   (input or output) DOUBLE PRECISION array, dimension (N)
!
!          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
!          contains the norm of the off-diagonal part of the j-th column
!          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
!          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
!          must be greater than or equal to the 1-norm.
!
!          If NORMIN = 'N', CNORM is an output argument and CNORM(j)
!          returns the 1-norm of the offdiagonal part of the j-th column
!          of A.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -k, the k-th argument had an illegal value
!
!  Further Details
!  ======= =======
!
!  A rough bound on x is computed; if that is less than overflow, ZTRSV
!  is called, otherwise, specific code is used which checks for possible
!  overflow or divide-by-zero at every operation.
!
!  A columnwise scheme is used for solving A*x = b.  The basic algorithm
!  if A is lower triangular is
!
!       x[1:n] := b[1:n]
!       for j = 1, ..., n
!            x(j) := x(j) / A(j,j)
!            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]
!       end
!
!  Define bounds on the components of x after j iterations of the loop:
!     M(j) = bound on x[1:j]
!     G(j) = bound on x[j+1:n]
!  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.
!
!  Then for iteration j+1 we have
!     M(j+1) <= G(j) / | A(j+1,j+1) |
!     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |
!            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )
!
!  where CNORM(j+1) is greater than or equal to the infinity-norm of
!  column j+1 of A, not counting the diagonal.  Hence
!
!     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )
!                  1<=i<=j
!  and
!
!     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )
!                                   1<=i< j
!
!  Since |x(j)| <= M(j), we use the Level 2 BLAS routine ZTRSV if the
!  reciprocal of the largest M(j), j=1,..,n, is larger than
!  max(underflow, 1/overflow).
!
!  The bound on x(j) is also used to determine when a step in the
!  columnwise method can be performed without fear of overflow.  If
!  the computed bound is greater than a large constant, x is scaled to
!  prevent overflow, but if the bound overflows, x is set to 0, x(j) to
!  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.
!
!  Similarly, a row-wise scheme is used to solve A**T *x = b  or
!  A**H *x = b.  The basic algorithm for A upper triangular is
!
!       for j = 1, ..., n
!            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)
!       end
!
!  We simultaneously compute two bounds
!       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j
!       M(j) = bound on x(i), 1<=i<=j
!
!  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we
!  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.
!  Then the bound on x(j) is
!
!       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |
!
!            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )
!                      1<=i<=j
!
!  and we can safely call ZTRSV if 1/M(n) and 1/G(n) are both greater
!  than max(underflow, 1/overflow).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0, &
                         TWO = 2.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN, NOUNIT, UPPER
      INTEGER            I, IMAX, J, JFIRST, JINC, JLAST
      DOUBLE PRECISION   BIGNUM, GROW, REC, SMLNUM, TJJ, TMAX, TSCAL, &
                         XBND, XJ, XMAX
      COMPLEX*16         CSUMJ, TJJS, USCAL, ZDUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX, IZAMAX
      DOUBLE PRECISION   DLAMCH, DZASUM
      COMPLEX*16         ZDOTC, ZDOTU, ZLADIV
      EXTERNAL           LSAME, IDAMAX, IZAMAX, DLAMCH, DZASUM, ZDOTC, &
                         ZDOTU, ZLADIV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLABAD, DSCAL, XERBLA, ZAXPY, ZDSCAL, ZTRSV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, MIN
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   CABS1, CABS2
!     ..
!     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      CABS2( ZDUM ) = ABS( DBLE( ZDUM ) / 2.D0 ) + &
                      ABS( DIMAG( ZDUM ) / 2.D0 )
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )
!
!     Test the input parameters.
!
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
               LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.LSAME( NORMIN, 'Y' ) .AND. .NOT. &
               LSAME( NORMIN, 'N' ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLATRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Determine machine dependent parameters to control overflow.
!
      SMLNUM = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SMLNUM / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM
      SCALE = ONE
!
      IF( LSAME( NORMIN, 'N' ) ) THEN
!
!        Compute the 1-norm of each column, not including the diagonal.
!
         IF( UPPER ) THEN
!
!           A is upper triangular.
!
            DO 10 J = 1, N
               CNORM( J ) = DZASUM( J-1, A( 1, J ), 1 )
   10       CONTINUE
         ELSE
!
!           A is lower triangular.
!
            DO 20 J = 1, N - 1
               CNORM( J ) = DZASUM( N-J, A( J+1, J ), 1 )
   20       CONTINUE
            CNORM( N ) = ZERO
         END IF
      END IF
!
!     Scale the column norms by TSCAL if the maximum element in CNORM is
!     greater than BIGNUM/2.
!
      IMAX = IDAMAX( N, CNORM, 1 )
      TMAX = CNORM( IMAX )
      IF( TMAX.LE.BIGNUM*HALF ) THEN
         TSCAL = ONE
      ELSE
         TSCAL = HALF / ( SMLNUM*TMAX )
         CALL DSCAL( N, TSCAL, CNORM, 1 )
      END IF
!
!     Compute a bound on the computed solution vector to see if the
!     Level 2 BLAS routine ZTRSV can be used.
!
      XMAX = ZERO
      DO 30 J = 1, N
         XMAX = MAX( XMAX, CABS2( X( J ) ) )
   30 CONTINUE
      XBND = XMAX
!
      IF( NOTRAN ) THEN
!
!        Compute the growth in A * x = b.
!
         IF( UPPER ) THEN
            JFIRST = N
            JLAST = 1
            JINC = -1
         ELSE
            JFIRST = 1
            JLAST = N
            JINC = 1
         END IF
!
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 60
         END IF
!
         IF( NOUNIT ) THEN
!
!           A is non-unit triangular.
!
!           Compute GROW = 1/G(j) and XBND = 1/M(j).
!           Initially, G(0) = max{x(i), i=1,...,n}.
!
            GROW = HALF / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 40 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               IF( GROW.LE.SMLNUM ) &
                  GO TO 60
!
               TJJS = A( J, J )
               TJJ = CABS1( TJJS )
!
               IF( TJJ.GE.SMLNUM ) THEN
!
!                 M(j) = G(j-1) / abs(A(j,j))
!
                  XBND = MIN( XBND, MIN( ONE, TJJ )*GROW )
               ELSE
!
!                 M(j) could overflow, set XBND to 0.
!
                  XBND = ZERO
               END IF
!
               IF( TJJ+CNORM( J ).GE.SMLNUM ) THEN
!
!                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
!
                  GROW = GROW*( TJJ / ( TJJ+CNORM( J ) ) )
               ELSE
!
!                 G(j) could overflow, set GROW to 0.
!
                  GROW = ZERO
               END IF
   40       CONTINUE
            GROW = XBND
         ELSE
!
!           A is unit triangular.
!
!           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
!
            GROW = MIN( ONE, HALF / MAX( XBND, SMLNUM ) )
            DO 50 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               IF( GROW.LE.SMLNUM ) &
                  GO TO 60
!
!              G(j) = G(j-1)*( 1 + CNORM(j) )
!
               GROW = GROW*( ONE / ( ONE+CNORM( J ) ) )
   50       CONTINUE
         END IF
   60    CONTINUE
!
      ELSE
!
!        Compute the growth in A**T * x = b  or  A**H * x = b.
!
         IF( UPPER ) THEN
            JFIRST = 1
            JLAST = N
            JINC = 1
         ELSE
            JFIRST = N
            JLAST = 1
            JINC = -1
         END IF
!
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 90
         END IF
!
         IF( NOUNIT ) THEN
!
!           A is non-unit triangular.
!
!           Compute GROW = 1/G(j) and XBND = 1/M(j).
!           Initially, M(0) = max{x(i), i=1,...,n}.
!
            GROW = HALF / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 70 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               IF( GROW.LE.SMLNUM ) &
                  GO TO 90
!
!              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
!
               XJ = ONE + CNORM( J )
               GROW = MIN( GROW, XBND / XJ )
!
               TJJS = A( J, J )
               TJJ = CABS1( TJJS )
!
               IF( TJJ.GE.SMLNUM ) THEN
!
!                 M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
!
                  IF( XJ.GT.TJJ ) &
                     XBND = XBND*( TJJ / XJ )
               ELSE
!
!                 M(j) could overflow, set XBND to 0.
!
                  XBND = ZERO
               END IF
   70       CONTINUE
            GROW = MIN( GROW, XBND )
         ELSE
!
!           A is unit triangular.
!
!           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
!
            GROW = MIN( ONE, HALF / MAX( XBND, SMLNUM ) )
            DO 80 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               IF( GROW.LE.SMLNUM ) &
                  GO TO 90
!
!              G(j) = ( 1 + CNORM(j) )*G(j-1)
!
               XJ = ONE + CNORM( J )
               GROW = GROW / XJ
   80       CONTINUE
         END IF
   90    CONTINUE
      END IF
!
      IF( ( GROW*TSCAL ).GT.SMLNUM ) THEN
!
!        Use the Level 2 BLAS solve if the reciprocal of the bound on
!        elements of X is not too small.
!
         CALL ZTRSV( UPLO, TRANS, DIAG, N, A, LDA, X, 1 )
      ELSE
!
!        Use a Level 1 BLAS solve, scaling intermediate results.
!
         IF( XMAX.GT.BIGNUM*HALF ) THEN
!
!           Scale X so that its components are less than or equal to
!           BIGNUM in absolute value.
!
            SCALE = ( BIGNUM*HALF ) / XMAX
            CALL ZDSCAL( N, SCALE, X, 1 )
            XMAX = BIGNUM
         ELSE
            XMAX = XMAX*TWO
         END IF
!
         IF( NOTRAN ) THEN
!
!           Solve A * x = b
!
            DO 120 J = JFIRST, JLAST, JINC
!
!              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
!
               XJ = CABS1( X( J ) )
               IF( NOUNIT ) THEN
                  TJJS = A( J, J )*TSCAL
               ELSE
                  TJJS = TSCAL
                  IF( TSCAL.EQ.ONE ) &
                     GO TO 110
               END IF
               TJJ = CABS1( TJJS )
               IF( TJJ.GT.SMLNUM ) THEN
!
!                    abs(A(j,j)) > SMLNUM:
!
                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                          Scale x by 1/b(j).
!
                        REC = ONE / XJ
                        CALL ZDSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X( J ) = ZLADIV( X( J ), TJJS )
                  XJ = CABS1( X( J ) )
               ELSE IF( TJJ.GT.ZERO ) THEN
!
!                    0 < abs(A(j,j)) <= SMLNUM:
!
                  IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
!                       to avoid overflow when dividing by A(j,j).
!
                     REC = ( TJJ*BIGNUM ) / XJ
                     IF( CNORM( J ).GT.ONE ) THEN
!
!                          Scale by 1/CNORM(j) to avoid overflow when
!                          multiplying x(j) times column j.
!
                        REC = REC / CNORM( J )
                     END IF
                     CALL ZDSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
                  X( J ) = ZLADIV( X( J ), TJJS )
                  XJ = CABS1( X( J ) )
               ELSE
!
!                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                    scale = 0, and compute a solution to A*x = 0.
!
                  DO 100 I = 1, N
                     X( I ) = ZERO
  100             CONTINUE
                  X( J ) = ONE
                  XJ = ONE
                  SCALE = ZERO
                  XMAX = ZERO
               END IF
  110          CONTINUE
!
!              Scale x if necessary to avoid overflow when adding a
!              multiple of column j of A.
!
               IF( XJ.GT.ONE ) THEN
                  REC = ONE / XJ
                  IF( CNORM( J ).GT.( BIGNUM-XMAX )*REC ) THEN
!
!                    Scale x by 1/(2*abs(x(j))).
!
                     REC = REC*HALF
                     CALL ZDSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                  END IF
               ELSE IF( XJ*CNORM( J ).GT.( BIGNUM-XMAX ) ) THEN
!
!                 Scale x by 1/2.
!
                  CALL ZDSCAL( N, HALF, X, 1 )
                  SCALE = SCALE*HALF
               END IF
!
               IF( UPPER ) THEN
                  IF( J.GT.1 ) THEN
!
!                    Compute the update
!                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
!
                     CALL ZAXPY( J-1, -X( J )*TSCAL, A( 1, J ), 1, X, &
                                 1 )
                     I = IZAMAX( J-1, X, 1 )
                     XMAX = CABS1( X( I ) )
                  END IF
               ELSE
                  IF( J.LT.N ) THEN
!
!                    Compute the update
!                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
!
                     CALL ZAXPY( N-J, -X( J )*TSCAL, A( J+1, J ), 1, &
                                 X( J+1 ), 1 )
                     I = J + IZAMAX( N-J, X( J+1 ), 1 )
                     XMAX = CABS1( X( I ) )
                  END IF
               END IF
  120       CONTINUE
!
         ELSE IF( LSAME( TRANS, 'T' ) ) THEN
!
!           Solve A**T * x = b
!
            DO 170 J = JFIRST, JLAST, JINC
!
!              Compute x(j) = b(j) - sum A(k,j)*x(k).
!                                    k<>j
!
               XJ = CABS1( X( J ) )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN
!
!                 If x(j) could overflow, scale x by 1/(2*XMAX).
!
                  REC = REC*HALF
                  IF( NOUNIT ) THEN
                     TJJS = A( J, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                  END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.ONE ) THEN
!
!                       Divide by A(j,j) when scaling x if A(j,j) > 1.
!
                     REC = MIN( ONE, REC*TJJ )
                     USCAL = ZLADIV( USCAL, TJJS )
                  END IF
                  IF( REC.LT.ONE ) THEN
                     CALL ZDSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
               END IF
!
               CSUMJ = ZERO
               IF( USCAL.EQ.DCMPLX( ONE ) ) THEN
!
!                 If the scaling needed for A in the dot product is 1,
!                 call ZDOTU to perform the dot product.
!
                  IF( UPPER ) THEN
                     CSUMJ = ZDOTU( J-1, A( 1, J ), 1, X, 1 )
                  ELSE IF( J.LT.N ) THEN
                     CSUMJ = ZDOTU( N-J, A( J+1, J ), 1, X( J+1 ), 1 )
                  END IF
               ELSE
!
!                 Otherwise, use in-line code for the dot product.
!
                  IF( UPPER ) THEN
                     DO 130 I = 1, J - 1
                        CSUMJ = CSUMJ + ( A( I, J )*USCAL )*X( I )
  130                CONTINUE
                  ELSE IF( J.LT.N ) THEN
                     DO 140 I = J + 1, N
                        CSUMJ = CSUMJ + ( A( I, J )*USCAL )*X( I )
  140                CONTINUE
                  END IF
               END IF
!
               IF( USCAL.EQ.DCMPLX( TSCAL ) ) THEN
!
!                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
!                 was not used to scale the dotproduct.
!
                  X( J ) = X( J ) - CSUMJ
                  XJ = CABS1( X( J ) )
                  IF( NOUNIT ) THEN
                     TJJS = A( J, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                     IF( TSCAL.EQ.ONE ) &
                        GO TO 160
                  END IF
!
!                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
!
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.SMLNUM ) THEN
!
!                       abs(A(j,j)) > SMLNUM:
!
                     IF( TJJ.LT.ONE ) THEN
                        IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                             Scale X by 1/abs(x(j)).
!
                           REC = ONE / XJ
                           CALL ZDSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
                     X( J ) = ZLADIV( X( J ), TJJS )
                  ELSE IF( TJJ.GT.ZERO ) THEN
!
!                       0 < abs(A(j,j)) <= SMLNUM:
!
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
!
                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL ZDSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                     X( J ) = ZLADIV( X( J ), TJJS )
                  ELSE
!
!                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                       scale = 0 and compute a solution to A**T *x = 0.
!
                     DO 150 I = 1, N
                        X( I ) = ZERO
  150                CONTINUE
                     X( J ) = ONE
                     SCALE = ZERO
                     XMAX = ZERO
                  END IF
  160             CONTINUE
               ELSE
!
!                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
!                 product has already been divided by 1/A(j,j).
!
                  X( J ) = ZLADIV( X( J ), TJJS ) - CSUMJ
               END IF
               XMAX = MAX( XMAX, CABS1( X( J ) ) )
  170       CONTINUE
!
         ELSE
!
!           Solve A**H * x = b
!
            DO 220 J = JFIRST, JLAST, JINC
!
!              Compute x(j) = b(j) - sum A(k,j)*x(k).
!                                    k<>j
!
               XJ = CABS1( X( J ) )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN
!
!                 If x(j) could overflow, scale x by 1/(2*XMAX).
!
                  REC = REC*HALF
                  IF( NOUNIT ) THEN
                     TJJS = DCONJG( A( J, J ) )*TSCAL
                  ELSE
                     TJJS = TSCAL
                  END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.ONE ) THEN
!
!                       Divide by A(j,j) when scaling x if A(j,j) > 1.
!
                     REC = MIN( ONE, REC*TJJ )
                     USCAL = ZLADIV( USCAL, TJJS )
                  END IF
                  IF( REC.LT.ONE ) THEN
                     CALL ZDSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
               END IF
!
               CSUMJ = ZERO
               IF( USCAL.EQ.DCMPLX( ONE ) ) THEN
!
!                 If the scaling needed for A in the dot product is 1,
!                 call ZDOTC to perform the dot product.
!
                  IF( UPPER ) THEN
                     CSUMJ = ZDOTC( J-1, A( 1, J ), 1, X, 1 )
                  ELSE IF( J.LT.N ) THEN
                     CSUMJ = ZDOTC( N-J, A( J+1, J ), 1, X( J+1 ), 1 )
                  END IF
               ELSE
!
!                 Otherwise, use in-line code for the dot product.
!
                  IF( UPPER ) THEN
                     DO 180 I = 1, J - 1
                        CSUMJ = CSUMJ + ( DCONJG( A( I, J ) )*USCAL )* &
                                X( I )
  180                CONTINUE
                  ELSE IF( J.LT.N ) THEN
                     DO 190 I = J + 1, N
                        CSUMJ = CSUMJ + ( DCONJG( A( I, J ) )*USCAL )* &
                                X( I )
  190                CONTINUE
                  END IF
               END IF
!
               IF( USCAL.EQ.DCMPLX( TSCAL ) ) THEN
!
!                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
!                 was not used to scale the dotproduct.
!
                  X( J ) = X( J ) - CSUMJ
                  XJ = CABS1( X( J ) )
                  IF( NOUNIT ) THEN
                     TJJS = DCONJG( A( J, J ) )*TSCAL
                  ELSE
                     TJJS = TSCAL
                     IF( TSCAL.EQ.ONE ) &
                        GO TO 210
                  END IF
!
!                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
!
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.SMLNUM ) THEN
!
!                       abs(A(j,j)) > SMLNUM:
!
                     IF( TJJ.LT.ONE ) THEN
                        IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                             Scale X by 1/abs(x(j)).
!
                           REC = ONE / XJ
                           CALL ZDSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
                     X( J ) = ZLADIV( X( J ), TJJS )
                  ELSE IF( TJJ.GT.ZERO ) THEN
!
!                       0 < abs(A(j,j)) <= SMLNUM:
!
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
!
                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL ZDSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                     X( J ) = ZLADIV( X( J ), TJJS )
                  ELSE
!
!                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                       scale = 0 and compute a solution to A**H *x = 0.
!
                     DO 200 I = 1, N
                        X( I ) = ZERO
  200                CONTINUE
                     X( J ) = ONE
                     SCALE = ZERO
                     XMAX = ZERO
                  END IF
  210             CONTINUE
               ELSE
!
!                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
!                 product has already been divided by 1/A(j,j).
!
                  X( J ) = ZLADIV( X( J ), TJJS ) - CSUMJ
               END IF
               XMAX = MAX( XMAX, CABS1( X( J ) ) )
  220       CONTINUE
         END IF
         SCALE = SCALE / TSCAL
      END IF
!
!     Scale the column norms by 1/TSCAL for return.
!
      IF( TSCAL.NE.ONE ) THEN
         CALL DSCAL( N, ONE / TSCAL, CNORM, 1 )
      END IF
!
      RETURN
!
!     End of ZLATRS
!
      END
      SUBROUTINE ZTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          COMPQ
      INTEGER            IFST, ILST, INFO, LDQ, LDT, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         Q( LDQ, * ), T( LDT, * )
!     ..
!
!  Purpose
!  =======
!
!  ZTREXC reorders the Schur factorization of a complex matrix
!  A = Q*T*Q**H, so that the diagonal element of T with row index IFST
!  is moved to row ILST.
!
!  The Schur form T is reordered by a unitary similarity transformation
!  Z**H*T*Z, and optionally the matrix Q of Schur vectors is updated by
!  postmultplying it with Z.
!
!  Arguments
!  =========
!
!  COMPQ   (input) CHARACTER*1
!          = 'V':  update the matrix Q of Schur vectors;
!          = 'N':  do not update Q.
!
!  N       (input) INTEGER
!          The order of the matrix T. N >= 0.
!
!  T       (input/output) COMPLEX*16 array, dimension (LDT,N)
!          On entry, the upper triangular matrix T.
!          On exit, the reordered upper triangular matrix.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= max(1,N).
!
!  Q       (input/output) COMPLEX*16 array, dimension (LDQ,N)
!          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
!          On exit, if COMPQ = 'V', Q has been postmultiplied by the
!          unitary transformation matrix Z which reorders T.
!          If COMPQ = 'N', Q is not referenced.
!
!  LDQ     (input) INTEGER
!          The leading dimension of the array Q.  LDQ >= max(1,N).
!
!  IFST    (input) INTEGER
!  ILST    (input) INTEGER
!          Specify the reordering of the diagonal elements of T:
!          The element with row index IFST is moved to row ILST by a
!          sequence of transpositions between adjacent elements.
!          1 <= IFST <= N; 1 <= ILST <= N.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            WANTQ
      INTEGER            K, M1, M2, M3
      DOUBLE PRECISION   CS
      COMPLEX*16         SN, T11, T22, TEMP
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARTG, ZROT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters.
!
      INFO = 0
      WANTQ = LSAME( COMPQ, 'V' )
      IF( .NOT.LSAME( COMPQ, 'N' ) .AND. .NOT.WANTQ ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.MAX( 1, N ) ) ) THEN
         INFO = -6
      ELSE IF( IFST.LT.1 .OR. IFST.GT.N ) THEN
         INFO = -7
      ELSE IF( ILST.LT.1 .OR. ILST.GT.N ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTREXC', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.1 .OR. IFST.EQ.ILST ) &
         RETURN
!
      IF( IFST.LT.ILST ) THEN
!
!        Move the IFST-th diagonal element forward down the diagonal.
!
         M1 = 0
         M2 = -1
         M3 = 1
      ELSE
!
!        Move the IFST-th diagonal element backward up the diagonal.
!
         M1 = -1
         M2 = 0
         M3 = -1
      END IF
!
      DO 10 K = IFST + M1, ILST + M2, M3
!
!        Interchange the k-th and (k+1)-th diagonal elements.
!
         T11 = T( K, K )
         T22 = T( K+1, K+1 )
!
!        Determine the transformation to perform the interchange.
!
         CALL ZLARTG( T( K, K+1 ), T22-T11, CS, SN, TEMP )
!
!        Apply transformation to the matrix T.
!
         IF( K+2.LE.N ) &
            CALL ZROT( N-K-1, T( K, K+2 ), LDT, T( K+1, K+2 ), LDT, CS, &
                       SN )
         CALL ZROT( K-1, T( 1, K ), 1, T( 1, K+1 ), 1, CS, &
                    DCONJG( SN ) )
!
         T( K, K ) = T22
         T( K+1, K+1 ) = T11
!
         IF( WANTQ ) THEN
!
!           Accumulate transformation in the matrix Q.
!
            CALL ZROT( N, Q( 1, K ), 1, Q( 1, K+1 ), 1, CS, &
                       DCONJG( SN ) )
         END IF
!
   10 CONTINUE
!
      RETURN
!
!     End of ZTREXC
!
      END
      SUBROUTINE ZLACON( N, V, X, EST, KASE )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            KASE, N
      DOUBLE PRECISION   EST
!     ..
!     .. Array Arguments ..
      COMPLEX*16         V( N ), X( N )
!     ..
!
!  Purpose
!  =======
!
!  ZLACON estimates the 1-norm of a square, complex matrix A.
!  Reverse communication is used for evaluating matrix-vector products.
!
!  Arguments
!  =========
!
!  N      (input) INTEGER
!         The order of the matrix.  N >= 1.
!
!  V      (workspace) COMPLEX*16 array, dimension (N)
!         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
!         (W is not returned).
!
!  X      (input/output) COMPLEX*16 array, dimension (N)
!         On an intermediate return, X should be overwritten by
!               A * X,   if KASE=1,
!               A' * X,  if KASE=2,
!         where A' is the conjugate transpose of A, and ZLACON must be
!         re-called with all the other parameters unchanged.
!
!  EST    (output) DOUBLE PRECISION
!         An estimate (a lower bound) for norm(A).
!
!  KASE   (input/output) INTEGER
!         On the initial call to ZLACON, KASE should be 0.
!         On an intermediate return, KASE will be 1 or 2, indicating
!         whether X should be overwritten by A * X  or A' * X.
!         On the final return from ZLACON, KASE will again be 0.
!
!  Further Details
!  ======= =======
!
!  Contributed by Nick Higham, University of Manchester.
!  Originally named CONEST, dated March 16, 1988.
!
!  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of
!  a real or complex matrix, with applications to condition estimation",
!  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      DOUBLE PRECISION   ONE, TWO
      PARAMETER          ( ONE = 1.0D0, TWO = 2.0D0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ), &
                         CONE = ( 1.0D0, 0.0D0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, ITER, J, JLAST, JUMP
      DOUBLE PRECISION   ALTSGN, ESTOLD, SAFMIN, TEMP
!     ..
!     .. External Functions ..
      INTEGER            IZMAX1
      DOUBLE PRECISION   DLAMCH, DZSUM1
      EXTERNAL           IZMAX1, DLAMCH, DZSUM1
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZCOPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX
!     ..
!     .. Save statement ..
      SAVE
!     ..
!     .. Executable Statements ..
!
      SAFMIN = DLAMCH( 'Safe minimum' )
      IF( KASE.EQ.0 ) THEN
         DO 10 I = 1, N
            X( I ) = DCMPLX( ONE / DBLE( N ) )
   10    CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      END IF
!
      GO TO ( 20, 40, 70, 90, 120 )JUMP
!
!     ................ ENTRY   (JUMP = 1)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
!
   20 CONTINUE
      IF( N.EQ.1 ) THEN
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
!        ... QUIT
         GO TO 130
      END IF
      EST = DZSUM1( N, X, 1 )
!
      DO 30 I = 1, N
         IF( ABS( X( I ) ).GT.SAFMIN ) THEN
            X( I ) = X( I ) / DCMPLX( ABS( X( I ) ) )
         ELSE
            X( I ) = CONE
         END IF
   30 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
!
!     ................ ENTRY   (JUMP = 2)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY ZTRANS(A)*X.
!
   40 CONTINUE
      J = IZMAX1( N, X, 1 )
      ITER = 2
!
!     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
!
   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = CZERO
   60 CONTINUE
      X( J ) = CONE
      KASE = 1
      JUMP = 3
      RETURN
!
!     ................ ENTRY   (JUMP = 3)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
   70 CONTINUE
      CALL ZCOPY( N, X, 1, V, 1 )
      ESTOLD = EST
      EST = DZSUM1( N, V, 1 )
!
!     TEST FOR CYCLING.
      IF( EST.LE.ESTOLD ) &
         GO TO 100
!
      DO 80 I = 1, N
         IF( ABS( X( I ) ).GT.SAFMIN ) THEN
            X( I ) = X( I ) / DCMPLX( ABS( X( I ) ) )
         ELSE
            X( I ) = CONE
         END IF
   80 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
!
!     ................ ENTRY   (JUMP = 4)
!     X HAS BEEN OVERWRITTEN BY ZTRANS(A)*X.
!
   90 CONTINUE
      JLAST = J
      J = IZMAX1( N, X, 1 )
      IF( ( DBLE( X( JLAST ) ).NE.ABS( DBLE( X( J ) ) ) ) .AND. &
          ( ITER.LT.ITMAX ) ) THEN
         ITER = ITER + 1
         GO TO 50
      END IF
!
!     ITERATION COMPLETE.  FINAL STAGE.
!
  100 CONTINUE
      ALTSGN = ONE
      DO 110 I = 1, N
         X( I ) = DCMPLX( ALTSGN*( ONE+DBLE( I-1 ) / DBLE( N-1 ) ) )
         ALTSGN = -ALTSGN
  110 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
!
!     ................ ENTRY   (JUMP = 5)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
  120 CONTINUE
      TEMP = TWO*( DZSUM1( N, X, 1 ) / DBLE( 3*N ) )
      IF( TEMP.GT.EST ) THEN
         CALL ZCOPY( N, X, 1, V, 1 )
         EST = TEMP
      END IF
!
  130 CONTINUE
      KASE = 0
      RETURN
!
!     End of ZLACON
!
      END
      SUBROUTINE ZDRSCL( N, SA, SX, INCX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SA
!     ..
!     .. Array Arguments ..
      COMPLEX*16         SX( * )
!     ..
!
!  Purpose
!  =======
!
!  ZDRSCL multiplies an n-element complex vector x by the real scalar
!  1/a.  This is done without overflow or underflow as long as
!  the final result x/a does not overflow or underflow.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of components of the vector x.
!
!  SA      (input) DOUBLE PRECISION
!          The scalar a which is used to divide each component of x.
!          SA must be >= 0, or the subroutine will divide by zero.
!
!  SX      (input/output) COMPLEX*16 array, dimension
!                         (1+(N-1)*abs(INCX))
!          The n-element vector x.
!
!  INCX    (input) INTEGER
!          The increment between successive values of the vector SX.
!          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            DONE
      DOUBLE PRECISION   BIGNUM, CDEN, CDEN1, CNUM, CNUM1, MUL, SMLNUM
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLABAD, ZDSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.LE.0 ) &
         RETURN
!
!     Get machine parameters
!
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
!
!     Initialize the denominator to SA and the numerator to 1.
!
      CDEN = SA
      CNUM = ONE
!
   10 CONTINUE
      CDEN1 = CDEN*SMLNUM
      CNUM1 = CNUM / BIGNUM
      IF( ABS( CDEN1 ).GT.ABS( CNUM ) .AND. CNUM.NE.ZERO ) THEN
!
!        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM.
!
         MUL = SMLNUM
         DONE = .FALSE.
         CDEN = CDEN1
      ELSE IF( ABS( CNUM1 ).GT.ABS( CDEN ) ) THEN
!
!        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM.
!
         MUL = BIGNUM
         DONE = .FALSE.
         CNUM = CNUM1
      ELSE
!
!        Multiply X by CNUM / CDEN and return.
!
         MUL = CNUM / CDEN
         DONE = .TRUE.
      END IF
!
!     Scale the vector X by MUL
!
      CALL ZDSCAL( N, MUL, SX, INCX )
!
      IF( .NOT.DONE ) &
         GO TO 10
!
      RETURN
!
!     End of ZDRSCL
!
      END
      SUBROUTINE ZLARTG( F, G, CS, SN, R )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   CS
      COMPLEX*16         F, G, R, SN
!     ..
!
!  Purpose
!  =======
!
!  ZLARTG generates a plane rotation so that
!
!     [  CS  SN  ]     [ F ]     [ R ]
!     [  __      ]  .  [   ]  =  [   ]   where CS**2 + |SN|**2 = 1.
!     [ -SN  CS  ]     [ G ]     [ 0 ]
!
!  This is a faster version of the BLAS1 routine ZROTG, except for
!  the following differences:
!     F and G are unchanged on return.
!     If G=0, then CS=1 and SN=0.
!     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
!        floating point operations.
!
!  Arguments
!  =========
!
!  F       (input) COMPLEX*16
!          The first component of vector to be rotated.
!
!  G       (input) COMPLEX*16
!          The second component of vector to be rotated.
!
!  CS      (output) DOUBLE PRECISION
!          The cosine of the rotation.
!
!  SN      (output) COMPLEX*16
!          The sine of the rotation.
!
!  R       (output) COMPLEX*16
!          The nonzero component of the rotated vector.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   D, DI, F1, F2, FA, G1, G2, GA
      COMPLEX*16         FS, GS, SS, T
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, SQRT
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   ABS1, ABSSQ
!     ..
!     .. Statement Function definitions ..
      ABS1( T ) = ABS( DBLE( T ) ) + ABS( DIMAG( T ) )
      ABSSQ( T ) = DBLE( T )**2 + DIMAG( T )**2
!     ..
!     .. Executable Statements ..
!
!     [ 25 or 38 ops for main paths ]
!
      IF( G.EQ.CZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.CZERO ) THEN
         CS = ZERO
!
         SN = DCONJG( G ) / ABS( G )
         R = ABS( G )
!
!         SN = ONE
!         R = G
!
      ELSE
         F1 = ABS1( F )
         G1 = ABS1( G )
         IF( F1.GE.G1 ) THEN
            GS = G / F1
            G2 = ABSSQ( GS )
            FS = F / F1
            F2 = ABSSQ( FS )
            D = SQRT( ONE+G2 / F2 )
            CS = ONE / D
            SN = DCONJG( GS )*FS*( CS / F2 )
            R = F*D
         ELSE
            FS = F / G1
            F2 = ABSSQ( FS )
            FA = SQRT( F2 )
            GS = G / G1
            G2 = ABSSQ( GS )
            GA = SQRT( G2 )
            D = SQRT( ONE+F2 / G2 )
            DI = ONE / D
            CS = ( FA / GA )*DI
            SS = ( DCONJG( GS )*FS ) / ( FA*GA )
            SN = SS*DI
            R = G*SS*D
         END IF
      END IF
      RETURN
!
!     End of ZLARTG
!
      END
      DOUBLE PRECISION FUNCTION DZSUM1( N, CX, INCX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         CX( * )
!     ..
!
!  Purpose
!  =======
!
!  DZSUM1 takes the sum of the absolute values of a complex
!  vector and returns a double precision result.
!
!  Based on DZASUM from the Level 1 BLAS.
!  The change is to use the 'genuine' absolute value.
!
!  Contributed by Nick Higham for use with ZLACON.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements in the vector CX.
!
!  CX      (input) COMPLEX*16 array, dimension (N)
!          The vector whose elements will be summed.
!
!  INCX    (input) INTEGER
!          The spacing between successive values of CX.  INCX > 0.
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, NINCX
      DOUBLE PRECISION   STEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      DZSUM1 = 0.0D0
      STEMP = 0.0D0
      IF( N.LE.0 ) &
         RETURN
      IF( INCX.EQ.1 ) &
         GO TO 20
!
!     CODE FOR INCREMENT NOT EQUAL TO 1
!
      NINCX = N*INCX
      DO 10 I = 1, NINCX, INCX
!
!        NEXT LINE MODIFIED.
!
         STEMP = STEMP + ABS( CX( I ) )
   10 CONTINUE
      DZSUM1 = STEMP
      RETURN
!
!     CODE FOR INCREMENT EQUAL TO 1
!
   20 CONTINUE
      DO 30 I = 1, N
!
!        NEXT LINE MODIFIED.
!
         STEMP = STEMP + ABS( CX( I ) )
   30 CONTINUE
      DZSUM1 = STEMP
      RETURN
!
!     End of DZSUM1
!
      END
      INTEGER          FUNCTION IZMAX1( N, CX, INCX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         CX( * )
!     ..
!
!  Purpose
!  =======
!
!  IZMAX1 finds the index of the element whose real part has maximum
!  absolute value.
!
!  Based on IZAMAX from Level 1 BLAS.
!  The change is to use the 'genuine' absolute value.
!
!  Contributed by Nick Higham for use with ZLACON.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements in the vector CX.
!
!  CX      (input) COMPLEX*16 array, dimension (N)
!          The vector whose elements will be summed.
!
!  INCX    (input) INTEGER
!          The spacing between successive values of CX.  INCX >= 1.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IX
      DOUBLE PRECISION   SMAX
      COMPLEX*16         ZDUM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function definitions ..
!
!     NEXT LINE IS THE ONLY MODIFICATION.
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) )
!     ..
!     .. Executable Statements ..
!
      IZMAX1 = 0
      IF( N.LT.1 ) &
         RETURN
      IZMAX1 = 1
      IF( N.EQ.1 ) &
         RETURN
      IF( INCX.EQ.1 ) &
         GO TO 30
!
!     CODE FOR INCREMENT NOT EQUAL TO 1
!
      IX = 1
      SMAX = CABS1( CX( 1 ) )
      IX = IX + INCX
      DO 20 I = 2, N
         IF( CABS1( CX( IX ) ).LE.SMAX ) &
            GO TO 10
         IZMAX1 = I
         SMAX = CABS1( CX( IX ) )
   10    CONTINUE
         IX = IX + INCX
   20 CONTINUE
      RETURN
!
!     CODE FOR INCREMENT EQUAL TO 1
!
   30 CONTINUE
      SMAX = CABS1( CX( 1 ) )
      DO 40 I = 2, N
         IF( CABS1( CX( I ) ).LE.SMAX ) &
            GO TO 40
         IZMAX1 = I
         SMAX = CABS1( CX( I ) )
   40 CONTINUE
      RETURN
!
!     End of IZMAX1
!
      END

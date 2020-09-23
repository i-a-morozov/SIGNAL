! GFORTRAN -O DSVD -STD=F2018 -WALL -PEDANTIC -O3 -FFAST-MATH -MARCH=NATIVE -FRECURSIVE DSVD.F90 -LBLAS -LARPACK
! ORIGINAL CODE, HTTPS://GITHUB.COM/OPENCOLLAB/ARPACK-NG/BLOB/MASTER/EXAMPLES/SVD/DSVD.F
! P110, HTTP://LI.MIT.EDU/ARCHIVE/ACTIVITIES/ARCHIVE/COURSEWORK/JU_LI/MITCOURSES/18.335/DOC/ARPACK/LEHOUCQ97.PDF

MODULE SVD

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: IK => C_INT, RK => C_DOUBLE
  IMPLICIT NONE

  PRIVATE

  REAL(RK), PUBLIC, PARAMETER :: SVD_LEVEL = 1.0E-9_RK  ! SINGULAR VALUE THRESHOLD
  INTEGER,  PUBLIC, PARAMETER :: SVD_ROW = 2_IK**12_IK  ! MAX NUMBER OF ROWS
  INTEGER,  PUBLIC, PARAMETER :: SVD_COL = 2_IK**12_IK  ! MAX NUMBER OF ROWS
  INTEGER,  PUBLIC, PARAMETER :: SVD_BASIS = 64_IK      ! MAX NUMBER OF BASIS VECTORS IN THE IMPLICITLY RESTARTED ARNOLDI PROCESS, OPTIMAL VALUE (?)
  INTEGER,  PUBLIC, PARAMETER :: SVD_LENGH = 16_IK      ! LENGTH OF ARNOLDI FACTORIZATION
  INTEGER,  PUBLIC, PARAMETER :: SVD_LOOP = 256_IK      ! MAX NUMBER OF MAIN LOOP ITERAIONS


  PUBLIC :: SVD_TRUNCATED_
  PUBLIC :: MATRIX_

  EXTERNAL :: DSAUPD
  EXTERNAL :: DSEUPD
  EXTERNAL :: DGEMV
  EXTERNAL :: DGEMM

  CONTAINS

  SUBROUTINE SVD_TRUNCATED_(NR,NC,NS,MATRIX,LIST,RVEC,LVEC)
    INTEGER(IK), INTENT(IN) :: NR
    INTEGER(IK), INTENT(IN) :: NC
    INTEGER(IK), INTENT(IN) :: NS
    REAL(RK), DIMENSION(NR, NC), INTENT(IN) :: MATRIX
    REAL(RK), DIMENSION(NS), INTENT(OUT) :: LIST
    REAL(RK), DIMENSION(NC, NS), INTENT(OUT) :: RVEC
    REAL(RK), DIMENSION(NR, NS), INTENT(OUT) :: LVEC
    REAL(RK), DIMENSION(NS, NS) :: DIAG
    REAL(RK), DIMENSION(NC, NS) :: COPY
    INTEGER(IK) :: I
    REAL(RK), DIMENSION(SVD_COL, SVD_BASIS) :: V = 0.0_RK
    REAL(RK), DIMENSION(SVD_BASIS*(SVD_BASIS+8_IK)) :: WL
    REAL(RK), DIMENSION(3_IK*SVD_COL) :: WD
    REAL(RK), DIMENSION(SVD_BASIS,2_IK) :: S = 0.0_RK
    REAL(RK), DIMENSION(SVD_COL) :: ID
    REAL(RK), DIMENSION(SVD_ROW) :: AX
    INTEGER(IK) :: PAR(11_IK), PNT(11_IK)
    LOGICAL :: SELECT(SVD_BASIS)
    CHARACTER :: MAT*1
    CHARACTER :: MODE*2
    INTEGER(IK) :: IDO, NEV, NCV, WORK, INFO, IERR, NCONV
    REAL(RK) :: TOL, SIGMA
    ! ARPACK
    NEV    = NS                    ! # OF SINGULAR VALUES TO COMPUTE, NEV < N
    NCV    = MIN(SVD_LENGH, NC)    ! LENGTH OF ARNOLDI FACTORIZATION
    MAT    = 'I'                   ! OP = A^T.A
    MODE   = 'LM'                  ! COMPUTE LARGEST (MAGNITUDE) SINGULAR VALUES, ALSO LA
    WORK   = NCV*(NCV+8_IK)        ! WORK ARRAY SIZE
    TOL    = 0.0_RK                ! TOLERANCE
    INFO   = 0_IK                  ! INITIAL ERROR CODE
    IDO    = 0_IK                  ! REVERSE COMMUNICATION PARAMETER
    PAR(1) = 1_IK                  ! SHIFT, 0/1
    PAR(3) = SVD_LOOP              ! MAX NUMBER OF ITERAIONS
    PAR(7) = 1_IK                  ! MODE
    ! MAIN LOOP
    MAIN : DO
      CALL DSAUPD(IDO,MAT,NC,MODE,NEV,TOL,ID,NCV,V,SVD_COL,PAR,PNT,WD,WL,WORK,INFO)
      IF (.NOT.(IDO .EQ. -1_IK .OR. IDO .EQ. 1_IK)) EXIT
      CALL DGEMV('N',NR,NC,1.0_RK,MATRIX,NR,WD(PNT(1)),1_IK,0.0_RK,AX,1_IK)
      CALL DGEMV('T',NR,NC,1.0_RK,MATRIX,NR,AX,1_IK,0.0_RK,WD(PNT(2)),1_IK)
    END DO MAIN
    ! RIGHT SINGULAR VERCTORS
    CALL DSEUPD (.TRUE.,'ALL',SELECT,S,V,SVD_COL,SIGMA,MAT,NC,MODE,NEV,TOL,ID,NCV,V,SVD_COL,PAR,PNT,WD,WL,WORK,IERR)
    ! SCALE, REVERSE AND SET SINGULAR VALUES
    NCONV =  PAR(5)
    LIST = SQRT(ABS(S(NEV:1:-1,1)))
    ! TRIM
    WHERE (LIST < SVD_LEVEL) LIST = 0.0_RK
    ! REVERSE AND SET RIGHT VECTORS
    RVEC(:,1:NEV:1) = V(1:NC,NEV:1:-1)
    ! COMPUTE LEFT VECTOR, U = A.V.S^-1, DIRECT COMPUTATION
    DIAG = 0.0_RK
    DO I = 1 , NEV, 1
      IF(LIST(I) > SVD_LEVEL) DIAG(I, I) = 1.0_RK/LIST(I)
    END DO
    CALL DGEMM('N','N',NC,NS,NS,1.0_RK,RVEC,NC,DIAG,NS,0.0_RK,COPY,NC)
    CALL DGEMM('N','N',NR,NS,NC,1.0_RK,MATRIX,NR,COPY,NC,0.0_RK,LVEC,NR)
  END SUBROUTINE SVD_TRUNCATED_

  SUBROUTINE MATRIX_(LENGTH, SEQUENCE, MATRIX)
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(LENGTH), INTENT(IN) :: SEQUENCE
    REAL(RK), DIMENSION(LENGTH/2_IK+1_IK, LENGTH/2_IK), INTENT(OUT) :: MATRIX
    INTEGER(IK) :: I
    DO I = 1_IK, LENGTH/2_IK+1_IK, 1_IK
      MATRIX(I, :) = SEQUENCE(I:I-1_IK+LENGTH/2_IK)
    END DO
  END SUBROUTINE MATRIX_

END MODULE SVD

PROGRAM DSVD
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: IK => C_INT, RK => C_DOUBLE
  USE SVD
  IMPLICIT NONE
  BLOCK
    INTEGER(IK), PARAMETER :: LENGTH = 2_IK**4
    INTEGER(IK), PARAMETER :: NR = LENGTH/2_IK + 1_IK
    INTEGER(IK), PARAMETER :: NC = LENGTH/2_IK
    INTEGER(IK), PARAMETER :: NS = 3_IK

    REAL(RK) :: SEQUENCE(LENGTH)
    REAL(RK) :: MATRIX(NR,NC)
    REAL(RK) :: SVD_LIST(NS)
    REAL(RK) :: DIAG(NS,NS)
    REAL(RK) :: RVEC(NC,NS)
    REAL(RK) :: LVEC(NR,NS)

    INTEGER :: I

    SEQUENCE = REAL([(I, I = 1, LENGTH, 1)], RK)
    CALL MATRIX_(LENGTH, SEQUENCE, MATRIX)

    WRITE(*, *) "INPUT MATRIX"
    WRITE(*, *) SIZE(MATRIX,1)
    WRITE(*, *) SIZE(MATRIX,2)
    DO I=1,NR,1
      WRITE(*, *) INT(MATRIX(I,:))
    END DO
    WRITE(*, *)

    CALL SVD_TRUNCATED_(NR, NC, NS, MATRIX, SVD_LIST, RVEC, LVEC)
    DIAG = 0.0_RK
    WRITE(*, *) "SINGULAR VALUES"
    DO I = 1, NS, 1
      DIAG(I, I) = SVD_LIST(I)
      WRITE(*, *) DIAG(I, :)
    END DO
    WRITE(*, *)

    WRITE(*, *) "RIGHT VECTORS"
    DO I = 1, NC, 1
      WRITE(*, *) RVEC(I,:)
    END DO
    WRITE(*, *)

    WRITE(*, *) "LEFT VECTORS"
    DO I = 1, NR, 1
      WRITE(*, *) LVEC(I, :)
    END DO
    WRITE(*, *)

    MATRIX = MATMUL(LVEC, MATMUL(DIAG, TRANSPOSE(RVEC)))
    WRITE(*, *) "U.S.V^T"
    DO I=1, NR, 1
      WRITE(*, *) CEILING(MATRIX(I, :))
    END DO
    WRITE(*, *)

  END BLOCK
END PROGRAM DSVD



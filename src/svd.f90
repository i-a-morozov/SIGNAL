#include "signal.inc"

SUBMODULE (SIGNAL) SVD
  IMPLICIT NONE
  CONTAINS
  ! ############################################################################################################################# !
  ! SVD (DGESVD)
  ! (SUBROUTINE) SVD_(<NR>, <NC>, <MATRIX>(<NR>, <NC>), <SVD_LIST>(MIN(<NR>, <NC>)), <U_MATRIX>(<NR>, <NR>), <V_MATRIX>(<NC>, <NC>))
  ! <NR>                   -- (IN)     NUMBER OF ROWS (IK)
  ! <NC>                   -- (IN)     NUMBER OF COLS (IK)
  ! <MATRIX>               -- (IN)     INPUT MATRIX(<NR>, <NC>) (RK)
  ! <SVD_LIST>             -- (OUT)    LIST OF SINGULAR VALUES WITH SIZE MIN(<NR>, <NC>) (RK)
  ! <U_MATRIX>             -- (OUT)    L-SINGULAR VECTORS (<NR>, <NR>) (RK)
  ! <V_MATRIX>             -- (OUT)    R-SINGULAR VECTORS (<NC>, <NC>) (RK)
  MODULE SUBROUTINE SVD_(NR, NC, MATRIX, SVD_LIST, U_MATRIX, V_MATRIX)
    INTEGER(IK), INTENT(IN) :: NR
    INTEGER(IK), INTENT(IN) :: NC
    REAL(RK), DIMENSION(NR, NC), INTENT(IN) :: MATRIX
    REAL(RK), DIMENSION(MIN(NR, NC)), INTENT(OUT) :: SVD_LIST
    REAL(RK), DIMENSION(NR, NR), INTENT(OUT) :: U_MATRIX
    REAL(RK), DIMENSION(NC, NC), INTENT(OUT) :: V_MATRIX
    REAL(RK), DIMENSION(2*MAX(1, 3*MIN(INT(NR), INT(NC))+MAX(INT(NR), INT(NC)), 5*MIN(INT(NR), INT(NC)))) :: WORK
    INTEGER :: WORK_SIZE
    INTEGER :: INFO
    WORK_SIZE = SIZE(WORK)
    CALL DGESVD('A','A',NR,NC,MATRIX,NR,SVD_LIST,U_MATRIX,NR,V_MATRIX,NC,WORK,WORK_SIZE,INFO)
    V_MATRIX = TRANSPOSE(V_MATRIX)
    WHERE (SVD_LIST < SVD_LEVEL) SVD_LIST = 0.0_RK
  END SUBROUTINE SVD_
  ! ############################################################################################################################# !
  ! SVD LIST (DGESVD)
  ! (SUBROUTINE) SVD_LIST_(<NR>, <NC>, <MATRIX>(<NR>, <NC>), <SVD_LIST>(MIN(<NR>, <NC>)))
  ! <NR>                   -- (IN)     NUMBER OF ROWS (IK)
  ! <NC>                   -- (IN)     NUMBER OF COLS (IK)
  ! <MATRIX>               -- (IN)     INPUT MATRIX(<NR>, <NC>) (RK)
  ! <SVD_LIST>             -- (OUT)    LIST OF SINGULAR VALUES WITH SIZE MIN(<NR>, <NC>) (RK)
  MODULE SUBROUTINE SVD_LIST_(NR, NC, MATRIX, SVD_LIST)
    INTEGER(IK), INTENT(IN) :: NR
    INTEGER(IK), INTENT(IN) :: NC
    REAL(RK), DIMENSION(NR, NC), INTENT(IN) :: MATRIX
    REAL(RK), DIMENSION(MIN(NR, NC)), INTENT(OUT) :: SVD_LIST
    REAL(RK), DIMENSION(2*MAX(1, 3*MIN(INT(NR), INT(NC))+MAX(INT(NR), INT(NC)), 5*MIN(INT(NR), INT(NC)))) :: WORK
    INTEGER :: WORK_SIZE
    INTEGER :: INFO
    REAL(RK), DIMENSION(NR, NR) :: U_MATRIX
    REAL(RK), DIMENSION(NC, NC) :: V_MATRIX
    WORK_SIZE = SIZE(WORK)
    CALL DGESVD('N','N',NR,NC,MATRIX,NR,SVD_LIST,U_MATRIX,NR,V_MATRIX,NC,WORK,WORK_SIZE,INFO)
    WHERE (SVD_LIST < SVD_LEVEL) SVD_LIST = 0.0_RK
  END SUBROUTINE SVD_LIST_
  ! ############################################################################################################################# !
  ! TRUNCATED SVD (ARPACK)
  ! SVD_TRUNCATED_(<NR>,<NC>,<NS>,<MATRIX>(<NR>,<NC>),<LIST>(<NS>),<RVEC>(<NC>,<NS>),<LVEC>(<NR>,<NS>))
  ! <NR>                   -- (IN)     NUMBER OF ROWS (IK)
  ! <NC>                   -- (IN)     NUMBER OF COLS (IK)
  ! <MATRIX>               -- (IN)     INPUT MATRIX(<NR>, <NC>) (RK)
  ! <LIST>                 -- (OUT)    LIST OF SINGULAR VALUES (<NS>) (RK)
  ! <RVEC>                 -- (OUT)    L-SINGULAR VECTORS (<NC>, <NS>) (RK)
  ! <LVEC>                 -- (OUT)    R-SINGULAR VECTORS (<NR>, <NS>) (RK)
  MODULE SUBROUTINE SVD_TRUNCATED_(NR, NC, NS, MATRIX, LIST, RVEC, LVEC)
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
    CHARACTER(LEN=1) :: MAT
    CHARACTER(LEN=2) :: MODE
    INTEGER(IK) :: IDO, NEV, NCV, WORK, INFO, IERR, NCONV
    REAL(RK) :: TOL, SIGMA
    ! ARPACK
    NEV    = NS                    ! # OF SINGULAR VALUES TO COMPUTE, NEV < N
    NCV    = MIN(SVD_LENGTH, NC)    ! LENGTH OF ARNOLDI FACTORIZATION
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
    LIST = SQRT(ABS(S(NEV:1_IK:-1_IK,1_IK)))
    ! TRIM
    WHERE (LIST < SVD_LEVEL) LIST = 0.0_RK
    ! REVERSE AND SET RIGHT VECTORS
    RVEC(:,1_IK:NEV:1_IK) = V(1_IK:NC,NEV:1_IK:-1_IK)
    ! COMPUTE LEFT VECTOR, U = A.V.S^-1, DIRECT COMPUTATION
    DIAG = 0.0_RK
    DO I = 1_IK , NEV, 1_IK
      IF(LIST(I) > SVD_LEVEL) DIAG(I, I) = 1.0_RK/LIST(I)
    END DO
    CALL DGEMM('N','N',NC,NS,NS,1.0_RK,RVEC,NC,DIAG,NS,0.0_RK,COPY,NC)
    CALL DGEMM('N','N',NR,NS,NC,1.0_RK,MATRIX,NR,COPY,NC,0.0_RK,LVEC,NR)
  END SUBROUTINE SVD_TRUNCATED_
  ! ############################################################################################################################# !
END SUBMODULE SVD
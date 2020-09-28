
#include "signal.inc"

SUBMODULE (SIGNAL) OPTIMIZATION
  IMPLICIT NONE
  CONTAINS
  ! ############################################################################################################################# !
  ! LEAST SQUARES (SVD)
  ! (SUBROUTINE) LEAST_SQUARES_(<NR>, <NC>, <MATRIX>(<NR>, <NC>), <VECTOR>(<NR>), <SOLUTION>(<NC>))
  ! <NR>                   -- (IN)     NUMBER OF ROWS (IK)
  ! <NC>                   -- (IN)     NUMBER OF COLS (IK)
  ! <MATRIX>               -- (IN)     INPUT DATA MATRIX (<NR>, <NC>) (RK)
  ! <VECTOR>               -- (IN)     INPUT VECTOR (<NR>) (RK)
  ! <SOLUTION>             -- (OUT)    LS SOLUTION (<NC>) (RK)
  MODULE SUBROUTINE LEAST_SQUARES_(NR, NC, MATRIX, VECTOR, SOLUTION)
    INTEGER(IK), INTENT(IN) :: NR
    INTEGER(IK), INTENT(IN) :: NC
    REAL(RK), DIMENSION(NR, NC), INTENT(IN) :: MATRIX
    REAL(RK), DIMENSION(NR), INTENT(IN) :: VECTOR
    REAL(RK), DIMENSION(NC), INTENT(OUT) :: SOLUTION
    REAL(RK), DIMENSION(MIN(NR, NC)) :: SVD_LIST
    REAL(RK), DIMENSION(NR, NR) :: U_MATRIX
    REAL(RK), DIMENSION(NC, NC) :: V_MATRIX
    REAL(RK), DIMENSION(NR, NC) :: COPY
    INTEGER(IK) :: I
    REAL(RK) :: V1(NR), V2(NR), V3(NC), V4(NC)
    CALL SVD_(NR, NC, MATRIX, SVD_LIST, U_MATRIX, V_MATRIX)
    COPY = 0.0_RK
    DO I = 1_IK, INT(SIZE(SVD_LIST), IK), 1_IK
      IF (SVD_LIST(I) >= SVD_LEVEL) COPY(I, I) = 1.0_RK/SVD_LIST(I)
    END DO
    V1 = VECTOR
    CALL DGEMV('T',NR,NR,1.0_RK,U_MATRIX,NR,V1,1,0.0_RK,V2,1)
    CALL DGEMV('T',NR,NC,1.0_RK,COPY,NR,V2,1,0.0_RK,V3,1)
    CALL DGEMV('N',NC,NC,1.0_RK,V_MATRIX,NC,V3,1,0.0_RK,V4,1)
    SOLUTION = REAL(V4, RK)
  END SUBROUTINE LEAST_SQUARES_
  ! ############################################################################################################################# !
  ! FIT (HARMONIC SIGNAL)
  ! (SUBROUTINE) FIT_(<LENGTH>, <SEQUENCE>, <LOOP>, <FREQUENCY>, <MEAN>, <COS_AMP>, <SIN_AMP>, <ERROR>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK), POWER OF TWO
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  ! <LOOP>                 -- (IN)     NUMBER OF HARMONICS (IK)
  ! <FREQUENCY>            -- (IN)     FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <MEAN>                 -- (OUT)    MEAN VALUE
  ! <COS_AMP>              -- (OUT)    COS AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <SIN_AMP>              -- (OUT)    SIN AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <ERROR>                -- (OUT)    ERROR
  ! <PV>                   -- (OUT)    P-VALUE (RK)
  ! void    fit_(int*, double*, int*, double*, double*, double*, double*, double*) ;
  MODULE SUBROUTINE FIT_(LENGTH, SEQUENCE, LOOP, FREQUENCY, MEAN, COS_AMP, SIN_AMP, ERROR) &
    BIND(C, NAME = "fit_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(LENGTH), INTENT(IN) :: SEQUENCE
    INTEGER(IK), INTENT(IN) :: LOOP
    REAL(RK), DIMENSION(LOOP), INTENT(IN) :: FREQUENCY
    REAL(RK), INTENT(OUT) :: MEAN
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: COS_AMP
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: SIN_AMP
    REAL(RK), INTENT(OUT) :: ERROR
    REAL(RK), DIMENSION(LENGTH) :: RANGE
    REAL(RK), DIMENSION(LENGTH, 2_IK*LOOP+1_IK) :: MATRIX
    REAL(RK), DIMENSION(2_IK*LOOP+1_IK) :: SOLUTION
    REAL(RK), DIMENSION(LENGTH) :: FIT
    INTEGER(IK) :: I
    RANGE = TWO_PI*REAL([(I, I = 1_IK, LENGTH, 1_IK)], RK)
    MATRIX = 1.0_RK
    DO I = 1_IK, LOOP, 1_IK
      MATRIX(:,2_IK*I)      = COS(FREQUENCY(I)*RANGE)
      MATRIX(:,2_IK*I+1_IK) = SIN(FREQUENCY(I)*RANGE)
    END DO
    CALL LEAST_SQUARES_(LENGTH, 2_IK*LOOP+1_IK, MATRIX, SEQUENCE, SOLUTION)
    MEAN = SOLUTION(1_IK)
    COS_AMP = SOLUTION(2_IK:2_IK*LOOP+1_IK:2_IK)
    SIN_AMP = SOLUTION(3_IK:2_IK*LOOP+1_IK:2_IK)
    FIT = MEAN
    DO I = 1_IK, LOOP, 1_IK
      FIT = FIT + COS_AMP(I)*COS(FREQUENCY(I)*RANGE) + SIN_AMP(I)*SIN(FREQUENCY(I)*RANGE)
    END DO
    ERROR = NORM2(SEQUENCE-FIT)
  END SUBROUTINE FIT_
  ! ############################################################################################################################# !
  ! BINARY SEARCH MAXIMIZATION
  ! (FUNCTION) BINARY_(<FUN>, <GUESS>, <INTERVAL>, <LIMIT>, <TOLERANCE>)
  ! <FUN>                  -- (IN)     FUNCTION TO MAXIMIZE (RK) -> (RK)
  ! <GUESS>                -- (IN)     INITIAL GUESS VALUE (RK)
  ! <INTERVAL>             -- (IN)     SEARCH INTERVAL (RK), GUESS IS IN THE MIDLE
  ! <LIMIT>                -- (IN)     MAXIMUM NUMBER OF ITERATIONS (IK)
  ! <TOLERANCE>            -- (IN)     MAXIMUM TOLERANCE (RK)
  ! <BINARY_>              -- (OUT)    MAXIMUM POSITION
  MODULE REAL(RK) FUNCTION BINARY_(FUN, GUESS, INTERVAL, LIMIT, TOLERANCE)
    INTERFACE
      REAL(RK) FUNCTION FUN(ARG)
        IMPORT :: RK
        REAL(RK), INTENT(IN) :: ARG
      END FUNCTION FUN
    END INTERFACE
    REAL(RK), INTENT(IN) :: GUESS
    REAL(RK), INTENT(IN) :: INTERVAL
    INTEGER(IK), INTENT(IN) :: LIMIT
    REAL(RK), INTENT(IN) :: TOLERANCE
    REAL(RK) :: DELTA
    REAL(RK) :: XL, XR, XX
    REAL(RK) :: FL, FR
    INTEGER(IK) :: I
    DELTA = INTERVAL/2.0_RK
    XL = GUESS-DELTA
    XR = GUESS+DELTA
    FL = FUN(FL)
    FR = FUN(FR)
    DO I = 1_IK, LIMIT, 1_IK
      IF (FL > FR) THEN
        XX = XL
      ELSE
        XX = XR
      END IF
      XL = XX-DELTA
      XR = XX+DELTA
      DELTA = DELTA/2.0_RK
      FL = FUN(XL)
      FR = FUN(XR)
      IF (ABS(FL-FR) < TOLERANCE) EXIT
    END DO
    IF (FL > FR) THEN
      BINARY_ = XL
    ELSE
      BINARY_ = XR
    END IF
  END FUNCTION BINARY_
  ! ############################################################################################################################# !
  ! GOLDEN SEARCH MAXIMIZATION
  ! (FUNCTION) GOLDEN_(<FUN>, <GUESS>, <INTERVAL>, <LIMIT>, <TOLERANCE>)
  ! <FUN>                  -- (IN)     FUNCTION TO MAXIMIZE (RK) -> (RK)
  ! <GUESS>                -- (IN)     INITIAL GUESS VALUE (RK)
  ! <INTERVAL>             -- (IN)     SEARCH INTERVAL (RK), GUESS IS IN THE MIDLE
  ! <LIMIT>                -- (IN)     MAXIMUM NUMBER OF ITERATIONS (IK)
  ! <TOLERANCE>            -- (IN)     MAXIMUM TOLERANCE (RK)
  ! <GOLDEN_>              -- (OUT)    MAXIMUM POSITION
  MODULE REAL(RK) FUNCTION GOLDEN_(FUN, GUESS, INTERVAL, LIMIT, TOLERANCE)
    INTERFACE
      REAL(RK) FUNCTION FUN(ARG)
        IMPORT :: RK
        REAL(RK), INTENT(IN) :: ARG
      END FUNCTION FUN
    END INTERFACE
    REAL(RK), INTENT(IN) :: GUESS
    REAL(RK), INTENT(IN) :: INTERVAL
    INTEGER(IK), INTENT(IN) :: LIMIT
    REAL(RK), INTENT(IN) :: TOLERANCE
    REAL(RK), PARAMETER :: GOLDEN = (1.0_RK+SQRT(5.0_RK))/2.0_RK
    REAL(RK), PARAMETER :: PSI = 1.0_RK-1.0_RK/GOLDEN
    REAL(RK), PARAMETER :: PHI = 1.0_RK/GOLDEN
    REAL(RK) :: DELTA
    REAL(RK) :: XL, XR
    REAL(RK) :: FL, FR
    REAL(RK) :: M1, M2
    INTEGER(IK) :: I
    DELTA = INTERVAL/2.0_RK
    XL = GUESS-DELTA
    XR = GUESS+DELTA
    M1 = XL+PSI*(XR-XL)
    M2 = XL+PHI*(XR-XL)
    FL = FUN(M1)
    FR = FUN(M2)
    DO I = 1_IK, LIMIT, 1_IK
      IF (I > LIMIT) EXIT
      IF (ABS(FL-FR) < TOLERANCE) EXIT
      IF (FL > FR) THEN
        XR = M2
        M2 = M1
        M1 = PSI*XL+PHI*M2
        FR = FL
        FL = FUN(M1)
      ELSE
        XL = M1
        M1 = M2
        M2 = PHI*M1+PSI*XR
        FL = FR
        FR = FUN(M2)
      END IF
    END DO
    IF (FL > FR) THEN
      GOLDEN_ = XL
    ELSE
      GOLDEN_ = XR
    END IF
  END FUNCTION GOLDEN_
  ! ############################################################################################################################# !
  ! REFINE FREQUENCY ESTIMATION (BINARY SEARCH)
  ! (FUNCTION) BINARY_AMPLITUDE_(<FLAG>, <LENGTH>, <TOTAL>, <WINDOW>, <SEQUENCE>, <GUESS>, <INTERVAL>, <LIMIT>, <TOLERANCE>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <GUESS>                -- (IN)     INITIAL GUESS VALUE (RK)
  ! <INTERVAL>             -- (IN)     SEARCH INTERVAL (RK), GUESS IS IN THE MIDLE
  ! <LIMIT>                -- (IN)     MAXIMUM NUMBER OF ITERATIONS (IK)
  ! <TOLERANCE>            -- (IN)     MAXIMUM TOLERANCE (RK)
  ! <BINARY_AMPLITUDE_>    -- (OUT)    REFINED FREQUENCY
  ! double  binary_amplitude_(int*, int*, double*, double*, double*, double*, double*, int*, double*) ;
  MODULE REAL(RK) FUNCTION BINARY_AMPLITUDE_(FLAG, LENGTH, TOTAL, WINDOW, SEQUENCE, GUESS, INTERVAL, LIMIT, TOLERANCE) &
    BIND(C, NAME = "binary_amplitude_")
    INTEGER(IK), INTENT(IN):: FLAG
    INTEGER(IK), INTENT(IN):: LENGTH
    REAL(RK), INTENT(IN) :: TOTAL
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
    REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    REAL(RK), INTENT(IN) :: GUESS
    REAL(RK), INTENT(IN) :: INTERVAL
    INTEGER(IK), INTENT(IN) :: LIMIT
    REAL(RK), INTENT(IN) :: TOLERANCE
    BINARY_AMPLITUDE_ = BINARY_(SEARCH_, GUESS, INTERVAL, LIMIT, TOLERANCE)
  CONTAINS
    REAL(RK) FUNCTION SEARCH_(FREQUENCY)
      REAL(RK), INTENT(IN) :: FREQUENCY
      REAL(RK) :: COS_AMP
      REAL(RK) :: SIN_AMP
      REAL(RK) :: AMPLITUDE
      CALL AMPLITUDE_(FLAG, LENGTH, TOTAL, WINDOW, SEQUENCE, FREQUENCY, COS_AMP, SIN_AMP, AMPLITUDE)
      SEARCH_ = AMPLITUDE
    END FUNCTION SEARCH_
  END FUNCTION
  ! ############################################################################################################################# !
  ! REFINE FREQUENCY ESTIMATION (GOLDEN SEARCH)
  ! (FUNCTION) GOLDEN_AMPLITUDE_(<FLAG>, <LENGTH>, <TOTAL>, <WINDOW>, <SEQUENCE>, <GUESS>, <INTERVAL>, <LIMIT>, <TOLERANCE>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <GUESS>                -- (IN)     INITIAL GUESS VALUE (RK)
  ! <INTERVAL>             -- (IN)     SEARCH INTERVAL (RK), GUESS IS IN THE MIDLE
  ! <LIMIT>                -- (IN)     MAXIMUM NUMBER OF ITERATIONS (IK)
  ! <TOLERANCE>            -- (IN)     MAXIMUM TOLERANCE (RK)
  ! <GOLDEN_AMPLITUDE_>    -- (OUT)    REFINED FREQUENCY
  ! double  golden_amplitude_(int*, int*, double*, double*, double*, double*, double*, int*, double*) ;
  MODULE REAL(RK) FUNCTION GOLDEN_AMPLITUDE_(FLAG, LENGTH, TOTAL, WINDOW, SEQUENCE, GUESS, INTERVAL, LIMIT, TOLERANCE) &
    BIND(C, NAME = "golden_amplitude_")
    INTEGER(IK), INTENT(IN):: FLAG
    INTEGER(IK), INTENT(IN):: LENGTH
    REAL(RK), INTENT(IN) :: TOTAL
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
    REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    REAL(RK), INTENT(IN) :: GUESS
    REAL(RK), INTENT(IN) :: INTERVAL
    INTEGER(IK), INTENT(IN) :: LIMIT
    REAL(RK), INTENT(IN) :: TOLERANCE
    GOLDEN_AMPLITUDE_ = GOLDEN_(SEARCH_, GUESS, INTERVAL, LIMIT, TOLERANCE)
  CONTAINS
    REAL(RK) FUNCTION SEARCH_(FREQUENCY)
      REAL(RK), INTENT(IN) :: FREQUENCY
      REAL(RK) :: COS_AMP
      REAL(RK) :: SIN_AMP
      REAL(RK) :: AMPLITUDE
      CALL AMPLITUDE_(FLAG, LENGTH, TOTAL, WINDOW, SEQUENCE, FREQUENCY, COS_AMP, SIN_AMP, AMPLITUDE)
      SEARCH_ = AMPLITUDE
    END FUNCTION SEARCH_
  END FUNCTION
  ! ############################################################################################################################# !
END SUBMODULE OPTIMIZATION
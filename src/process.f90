
#include "signal.inc"

SUBMODULE (SIGNAL) PROCESS
  IMPLICIT NONE
  CONTAINS
  ! ############################################################################################################################# !
  ! CONVERT INPUT SEQUENCE (REAL)
  ! (SUBROUTINE) CONVERT_REAL_(<LENGTH>, <R_PART>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <R_PART>               -- (IN)     INPUT SEQUENCE R-PART (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (OUT)    SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...] AND SI_I=0.0_RK FOR ALL I
  ! void    convert_real_(int* length, double* r_part, double* sequence) ;
  MODULE SUBROUTINE CONVERT_REAL_(LENGTH, R_PART, SEQUENCE) &
    BIND(C, NAME = "convert_real_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: R_PART
    REAL(RK), INTENT(OUT), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    SEQUENCE = 0.0_RK
    SEQUENCE(1_IK::2_IK) = R_PART
  END SUBROUTINE CONVERT_REAL_
  ! ############################################################################################################################# !
  ! CONVERT INPUT SEQUENCE (COMPLEX)
  ! (SUBROUTINE) CONVERT_COMPLEX_(<LENGTH>, <R_PART>, <I_PART>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <R_PART>               -- (IN)     INPUT SEQUENCE R-PART (RK ARRAY OF LENGTH = <LENGTH>)
  ! <I_PART>               -- (IN)     INPUT SEQUENCE I-PART (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (OUT)    SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! void    convert_complex_(int* length, double* r_part, double* i_part, double* sequence) ;
  MODULE SUBROUTINE CONVERT_COMPLEX_(LENGTH, R_PART, I_PART, SEQUENCE) &
    BIND(C, NAME = "convert_complex_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: R_PART
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: I_PART
    REAL(RK), INTENT(OUT), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    SEQUENCE = 0.0_RK
    SEQUENCE(1_IK::2_IK) = R_PART
    SEQUENCE(2_IK::2_IK) = I_PART
  END SUBROUTINE CONVERT_COMPLEX_
  ! ############################################################################################################################# !
  ! ROUND UP (ROUND UP TO THE NEXT POWER OF TWO)
  ! (FUNCTION) ROUND_UP_(<NUMBER>)
  ! <NUMBER>               -- (IN)     NUMBER (IK)
  ! <ROUND_UP>             -- (OUT)    NEXT POWER OF TWO NUMBER (IK)
  ! int     round_up_(int* number) ;
  MODULE INTEGER(IK) FUNCTION ROUND_UP_(NUMBER) &
    BIND(C, NAME = "round_up_")
    INTEGER(IK), INTENT(IN) :: NUMBER
    ROUND_UP_ = 2_IK**CEILING(LOG(REAL(NUMBER, RK))/LOG(2.0_RK), KIND=IK)
  END FUNCTION ROUND_UP_
  ! ############################################################################################################################# !
  ! ZERO PADDING
  ! (SUBROUTINE) PAD_(<LI>, <LO>, <INPUT>, <OUTPUT>)
  ! <LI>                   -- (IN)     INPUT SEQUENCE LENGTH (IK)
  ! <LO>                   -- (IN)     OUTPUT SEQUENCE LENGTH (IK)
  ! <INPUT>                -- (IN)     INPUT SEQUENCE (RK) OF LENGTH = 2*<LI>
  ! <OUTPUT>               -- (IN)     PADDED SEQUENCE (RK) OF LENGTH = 2*<LO>
  ! void    pad_(int* linput, int* loutput, double* input, double* output) ;
  MODULE SUBROUTINE PAD_(LI, LO, INPUT, OUTPUT) &
    BIND(C, NAME = "pad_")
    INTEGER(IK), INTENT(IN) :: LI
    INTEGER(IK), INTENT(IN) :: LO
    REAL(RK), DIMENSION(2_IK*LI), INTENT(IN) :: INPUT
    REAL(RK), DIMENSION(2_IK*LO), INTENT(OUT) :: OUTPUT
    REAL(RK) :: LZERO(FLOOR(REAL(LO-LI,RK)/2_RK))
    REAL(RK) :: RZERO(CEILING(REAL(LO-LI,RK)/2_RK))
    INTEGER(IK) :: I
    LZERO = [(0.0_RK, I = 1_IK, SIZE(LZERO), 1_IK)]
    RZERO = [(0.0_RK, I = 1_IK, SIZE(RZERO), 1_IK)]
    OUTPUT(1_IK::2_IK) = [LZERO, INPUT(1_IK::2_IK), RZERO]
    OUTPUT(2_IK::2_IK) = [LZERO, INPUT(2_IK::2_IK), RZERO]
  END  SUBROUTINE PAD_
  ! ############################################################################################################################# !
  ! REMOVE MEAN
  ! (SUBROUTINE) REMOVE_MEAN_(<LENGTH>, <INPUT>, <OUTPUT> )
  ! <LENGTH>               -- (IN)     INPUT SEQUENCE LENGTH (IK)
  ! <INPUT>                -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <OUTPUT>               -- (OUT)    OUTPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! void    remove_mean_(int* length, double* input, double* output) ;
  MODULE SUBROUTINE REMOVE_MEAN_(LENGTH, INPUT, OUTPUT) &
    BIND(C, NAME = "remove_mean_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(IN) :: INPUT
    REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(OUT) :: OUTPUT
    OUTPUT(1_IK::2_IK)  = INPUT(1_IK::2_IK)-SUM(INPUT(1_IK::2_IK))/REAL(LENGTH, RK)
    OUTPUT(2_IK::2_IK)  = INPUT(2_IK::2_IK)-SUM(INPUT(2_IK::2_IK))/REAL(LENGTH, RK)
  END SUBROUTINE REMOVE_MEAN_
  ! ############################################################################################################################# !
  ! REMOVE WINDOW MEAN
  ! (SUBROUTINE) REMOVE_WINDOW_MEAN_(<LENGTH>, <TOTAL>, <WINDOW>, <INPUT>, <OUTPUT> )
  ! <LENGTH>               -- (IN)     INPUT SEQUENCE LENGTH (IK)
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <INPUT>                -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <OUTPUT>               -- (OUT)    OUTPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! void    remove_window_mean_(int* length, double* total, double* window, double* input, double* output) ;
  MODULE SUBROUTINE REMOVE_WINDOW_MEAN_(LENGTH, TOTAL, WINDOW, INPUT, OUTPUT) &
    BIND(C, NAME = "remove_window_mean_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), INTENT(IN) :: TOTAL
    REAL(RK), DIMENSION(LENGTH), INTENT(IN) :: WINDOW
    REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(IN) :: INPUT
    REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(OUT) :: OUTPUT
    OUTPUT(1_IK::2_IK) = (INPUT(1_IK::2_IK)-SUM(WINDOW*INPUT(1_IK::2_IK))/TOTAL)
    OUTPUT(2_IK::2_IK) = (INPUT(2_IK::2_IK)-SUM(WINDOW*INPUT(2_IK::2_IK))/TOTAL)
  END SUBROUTINE REMOVE_WINDOW_MEAN_
  ! ############################################################################################################################# !
  ! APPLY WINDOW
  ! (SUBROUTINE) APPLY_WINDOW_(<LENGTH>, <WINDOW>, <INPUT>, <OUTPUT> )
  ! <LENGTH>               -- (IN)     INPUT SEQUENCE LENGTH (IK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <INPUT>                -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <OUTPUT>               -- (OUT)    OUTPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! void    apply_window_(int* length, double* window, double* input, double* output) ;
  MODULE SUBROUTINE APPLY_WINDOW_(LENGTH, WINDOW, INPUT, OUTPUT) &
    BIND(C, NAME = "apply_window_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(LENGTH), INTENT(IN) :: WINDOW
    REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(IN) :: INPUT
    REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(OUT) :: OUTPUT
    OUTPUT(1_IK::2_IK) = INPUT(1_IK::2_IK)*WINDOW
    OUTPUT(2_IK::2_IK) = INPUT(2_IK::2_IK)*WINDOW
  END SUBROUTINE APPLY_WINDOW_
  ! ############################################################################################################################# !
  ! MATRIX (GENERATE MATRIX FROM SEQUENCE)
  ! (SUBROUTINE) MATRIX_(<LENGTH>, <SEQUENCE>, <MATRIX>)
  ! <LENGTH>               -- (IN)     INPUT SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK)
  ! <MATRIX>               -- (OUT)    MATRIX (<LENGTH>/2+1, <LENGTH>/2) (RK)
  MODULE SUBROUTINE MATRIX_(LENGTH, SEQUENCE, MATRIX)
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(LENGTH), INTENT(IN) :: SEQUENCE
    REAL(RK), DIMENSION(LENGTH/2_IK+1_IK, LENGTH/2_IK), INTENT(OUT) :: MATRIX
    INTEGER(IK) :: I
    DO I = 1_IK, LENGTH/2_IK+1_IK, 1_IK
      MATRIX(I, :) = SEQUENCE(I:I-1_IK+LENGTH/2_IK)
    END DO
  END SUBROUTINE MATRIX_
  ! ############################################################################################################################# !
  ! SEQUENCE (ROW) (GENERATE SEQUENCE FROM MATRIX USING 1ST AND LAST ROWS)
  ! (SUBROUTINE) SEQUENCE_ROW_(<LENGTH>, <SEQUENCE>, <MATRIX>)
  ! <LENGTH>               -- (IN)     INPUT SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (OUT)    INPUT SEQUENCE (RK)
  ! <MATRIX>               -- (IN)     MATRIX (<LENGTH>/2+1, <LENGTH>/2) (RK)
  MODULE SUBROUTINE SEQUENCE_ROW_(LENGTH, SEQUENCE, MATRIX)
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(LENGTH), INTENT(OUT) :: SEQUENCE
    REAL(RK), DIMENSION(LENGTH/2_IK+1_IK, LENGTH/2_IK), INTENT(IN) :: MATRIX
    SEQUENCE = [MATRIX(1_IK, :), MATRIX(LENGTH/2_IK+1_IK, :)]
  END SUBROUTINE SEQUENCE_ROW_
  ! ############################################################################################################################# !
  ! SEQUENCE (SUM) (GENERATE SEQUENCE FROM MATRIX USING SUMS OF SKEW DIAGONALS)
  ! (SUBROUTINE) SEQUENCE_ROW_(<LENGTH>, <SEQUENCE>, <MATRIX>)
  ! <LENGTH>               -- (IN)     INPUT SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (OUT)    INPUT SEQUENCE (RK)
  ! <MATRIX>               -- (IN)     MATRIX (<LENGTH>/2+1, <LENGTH>/2) (RK)
  MODULE SUBROUTINE SEQUENCE_SUM_(LENGTH, SEQUENCE, MATRIX)
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(LENGTH), INTENT(OUT) :: SEQUENCE
    REAL(RK), DIMENSION(LENGTH/2_IK+1_IK, LENGTH/2_IK), INTENT(IN) :: MATRIX
    INTEGER(IK) :: ROW
    INTEGER(IK) :: COL
    REAL(RK), DIMENSION((LENGTH/2_IK+1_IK)*(LENGTH/2_IK)) :: ARRAY
    INTEGER(IK) :: SHIFT
    INTEGER(IK) :: START
    INTEGER(IK) :: COUNT
    INTEGER(IK) :: Q, P
    ROW = LENGTH/2_IK+1_IK
    COL = LENGTH/2_IK
    ARRAY = RESHAPE(MATRIX, SHAPE(ARRAY))
    SHIFT = 1_IK ;
    SEQUENCE = 0.0_RK
    DO Q = 1_IK, 2_IK*COL, 1_IK
      START = Q
      IF (Q > ROW) THEN
        START = Q + COL*SHIFT
        SHIFT = SHIFT + 1_IK
      END IF
      COUNT = 0_IK
      DO P = 0_IK, MIN(Q, COL)-SHIFT, 1_IK
        SEQUENCE(Q) = SEQUENCE(Q)+ARRAY(START+P*COL)
        COUNT = COUNT + 1_IK
      END DO
      SEQUENCE(Q) = SEQUENCE(Q)/REAL(COUNT, RK) ;
    END DO
  END SUBROUTINE SEQUENCE_SUM_
  ! ############################################################################################################################# !
  ! FILTER
  ! (SUBROUTINE) FILTER(<LENGTH>, <SEQUENCE>, <LIMIT>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <SEQUENCE>             -- (INOUT)  SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  ! <LIMIT>                -- (IN)     NUMBER OF SINGULAR VALUES TO KEEP (IK)
  ! <SVD_LIST>             -- (OUT)    LIST OF SINGULAR VALUES
  ! void    filter_(int* length, double* sequence, int* limit, double* svd_list) ;
  MODULE SUBROUTINE FILTER_(LENGTH, SEQUENCE, LIMIT, SVD_LIST) &
    BIND(C, NAME = "filter_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(LENGTH), INTENT(INOUT) :: SEQUENCE
    INTEGER(IK), INTENT(IN) :: LIMIT
    REAL(RK), DIMENSION(LIMIT), INTENT(OUT) :: SVD_LIST
    REAL(RK), DIMENSION(LENGTH/2_IK+1_IK, LENGTH/2_IK) :: MATRIX
    REAL(RK), DIMENSION(LENGTH/2_IK+1_IK, LIMIT) :: U_MATRIX
    REAL(RK), DIMENSION(LENGTH/2_IK, LIMIT) :: V_MATRIX
    REAL(RK), DIMENSION(LIMIT, LIMIT) :: DIAG
    REAL(RK), DIMENSION(LIMIT, LENGTH/2_IK) :: COPY
    INTEGER(IK) :: I
    INTEGER :: NR, NC, NS
    CALL MATRIX_(LENGTH, SEQUENCE, MATRIX)
    CALL SVD_TRUNCATED_(LENGTH/2_IK+1_IK, LENGTH/2_IK, LIMIT, MATRIX, SVD_LIST, V_MATRIX, U_MATRIX)
    MATRIX = 0.0_RK
    DO I = 1_IK, LIMIT, 1_IK
      DIAG(I, I) = SVD_LIST(I)
    END DO
    NR = INT(LENGTH)/2_IK + 1_IK
    NC = INT(LENGTH)/2_IK
    NS = LIMIT
    CALL DGEMM('N','T',NS,NC,NS,1.0_RK,DIAG,SIZE(DIAG,1_IK),V_MATRIX,SIZE(V_MATRIX,1_IK),0.0_RK,COPY,SIZE(COPY,1_IK))
    CALL DGEMM('N','N',NR,NC,NS,1.0_RK,U_MATRIX,SIZE(U_MATRIX,1_IK),COPY,SIZE(COPY,1_IK),0.0_RK,MATRIX,SIZE(MATRIX,1_IK))
    CALL __SEQUENCE__(LENGTH, SEQUENCE, MATRIX)
  END SUBROUTINE FILTER_
  ! ############################################################################################################################# !
END SUBMODULE PROCESS
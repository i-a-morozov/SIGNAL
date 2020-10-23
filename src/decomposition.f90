
#include "signal.inc"

SUBMODULE (SIGNAL) DECOMPOSITION
  IMPLICIT NONE
  CONTAINS
  ! ############################################################################################################################# !
  ! ESTIMATE AMPLITUDE FOR GIVEN FREQUENCY
  ! (SUBROUTINE) AMPLITUDE_(<FLAG>, <LENGTH>, <TOTAL>, <WINDOW>, <SEQUENCE>, <FREQUENCY>, <COS_AMP>, <SIN_AMP>, <AMP>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <FREQUENCY>            -- (IN)     FREQUENCY (RK)
  ! <COS_AMP>              -- (OUT)    COS AMPLITUDE (RK)
  ! <SIN_AMP>              -- (OUT)    SIN AMPLITUDE (RK)
  ! <AMP>                  -- (OUT)    ABS AMPLITUDE (RK)
  ! void    amplitude_(int* flag, int* length, double* total, double* window, double* sequence, double* frequency, double* cos_amp, double* sin_amp, double* amp) ;
  MODULE SUBROUTINE AMPLITUDE_(FLAG, LENGTH, TOTAL, WINDOW, SEQUENCE, FREQUENCY, COS_AMP, SIN_AMP, AMP) &
    BIND(C, NAME = "amplutude_")
    INTEGER(IK), INTENT(IN):: FLAG
    INTEGER(IK), INTENT(IN):: LENGTH
    REAL(RK), INTENT(IN) :: TOTAL
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
    REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    REAL(RK), INTENT(IN) :: FREQUENCY
    REAL(RK), INTENT(OUT) :: COS_AMP
    REAL(RK), INTENT(OUT) :: SIN_AMP
    REAL(RK), INTENT(OUT) :: AMP
    REAL(RK), DIMENSION(2_IK*LENGTH) :: LOCAL
    REAL(RK), DIMENSION(LENGTH) :: LIST
    REAL(RK), DIMENSION(2_IK*LENGTH) :: DELTA
    INTEGER(IK) :: I
    CALL REMOVE_WINDOW_MEAN_(LENGTH, TOTAL, WINDOW, SEQUENCE, LOCAL)
    LIST = TWO_PI*REAL([(I, I=1_IK, LENGTH, 1_IK)], RK)
    DELTA(1_IK::2_IK) = +COS(FREQUENCY*LIST)
    DELTA(2_IK::2_IK) = -SIN(FREQUENCY*LIST)
    COS_AMP = SUM(LOCAL(1_IK::2_IK)*DELTA(1_IK::2_IK)*WINDOW)+SUM(LOCAL(2_IK::2_IK)*DELTA(2_IK::2_IK)*WINDOW)
    SIN_AMP = SUM(LOCAL(2_IK::2_IK)*DELTA(1_IK::2_IK)*WINDOW)-SUM(LOCAL(1_IK::2_IK)*DELTA(2_IK::2_IK)*WINDOW)
    COS_AMP = 1.0_RK/TOTAL*COS_AMP
    SIN_AMP = 1.0_RK/TOTAL*SIN_AMP
    IF (FLAG == 0_IK) THEN
      COS_AMP = 2.0_RK*COS_AMP
      SIN_AMP = 2.0_RK*SIN_AMP
    END IF
    AMP = SQRT(COS_AMP**2_IK+SIN_AMP**2_IK)
  END SUBROUTINE AMPLITUDE_
  ! ############################################################################################################################# !
  ! SIGNAL DECOMPOSITION
  ! (SUBROUTINE) DECOMPOSITION_(<FLAG>, <RANGE_MIN>, <RANGE_MAX>, <METHOD>, <MODE>, <LENGTH>, <PAD>, <TOTAL>, <WINDOW>, <SEQUENCE>, <LOOP>, <FREQUENCY>, <COS_AMP>, <SIN_AMP>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <RANGE_MIN>            -- (IN)     (MIN) FREQUENCY RANGE (RK)
  ! <RANGE_MAX>            -- (IN)     (MAX) FREQUENCY RANGE (RK)
  ! <METHOD>               -- (IN)     FREQUENCY APPROXIMATION METHOD (IK), FREQUENCY_FFT = 0_IK, FREQUENCY_FFRFT = 1_IK, FREQUECY_PARABOLA = 2_IK
  ! <MODE>                 -- (IN)     DECOMPOSTION MODE (IK), <MODE> = DECOMPOSITION_SUBTRACT = 0 OR <MODE> = DECOMPOSITION_PEAK = 1
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <PAD>                  -- (IN)     PADDED SEQUENCE LENGTH (IK), IF PAD > LENGTH, INPUT SEQUENCE IS PADDED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <LOOP>                 -- (IN)     NUMBER OF ITERATIONS/PEAKS (IK)
  ! <FREQUENCY>            -- (OUT)    FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <COS_AMP>              -- (OUT)    COS AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <SIN_AMP>              -- (OUT)    SIN AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! void    decomposition_(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, double* sequence, int* loop, double* frequency, double* cos_amp, double* sin_amp) ;
  MODULE SUBROUTINE DECOMPOSITION_(FLAG, RANGE_MIN, RANGE_MAX, &
    METHOD, MODE, LENGTH, PAD, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP) &
    BIND(C, NAME = "decomposition_")
    INTEGER(IK), INTENT(IN):: FLAG
    REAL(RK), INTENT(IN) :: RANGE_MIN
    REAL(RK), INTENT(IN) :: RANGE_MAX
    INTEGER(IK), INTENT(IN):: METHOD
    INTEGER(IK), INTENT(IN):: MODE
    INTEGER(IK), INTENT(IN):: LENGTH
    INTEGER(IK), INTENT(IN):: PAD
    REAL(RK), INTENT(IN) :: TOTAL
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
    REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    INTEGER(IK), INTENT(IN) :: LOOP
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: FREQUENCY
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: COS_AMP
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: SIN_AMP
    REAL(RK), DIMENSION(2_IK*LENGTH) :: LOCAL
    INTEGER(IK) :: I
    REAL(RK), DIMENSION(LENGTH) :: LIST
    REAL(RK), DIMENSION(2_IK*LENGTH) :: DELTA
    INTEGER(IK), DIMENSION(LOOP) :: ORDERING
    REAL(RK), DIMENSION(LOOP) :: AMPLITUDE
    LOCAL = SEQUENCE
    LIST = TWO_PI*REAL([(I, I=1_IK, LENGTH, 1_IK)], RK)
    DO I = 1_IK, LOOP, 1_IK
      IF (MODE == DECOMPOSITION_SUBTRACT) THEN
        FREQUENCY(I) = FREQUENCY_(FLAG, RANGE_MIN, RANGE_MAX, 0_IK, METHOD, LENGTH, PAD, TOTAL, WINDOW, LOCAL)
      ELSE
        FREQUENCY(I) = FREQUENCY_(FLAG, RANGE_MIN, RANGE_MAX, I, METHOD, LENGTH, PAD, TOTAL, WINDOW, SEQUENCE)
      END IF
      DELTA(1_IK::2_IK) = +COS(FREQUENCY(I)*LIST)
      DELTA(2_IK::2_IK) = -SIN(FREQUENCY(I)*LIST)
      COS_AMP(I) = SUM(LOCAL(1_IK::2_IK)*DELTA(1_IK::2_IK)*WINDOW)+SUM(LOCAL(2_IK::2_IK)*DELTA(2_IK::2_IK)*WINDOW)
      SIN_AMP(I) = SUM(LOCAL(2_IK::2_IK)*DELTA(1_IK::2_IK)*WINDOW)-SUM(LOCAL(1_IK::2_IK)*DELTA(2_IK::2_IK)*WINDOW)
      COS_AMP(I) = 1.0_RK/TOTAL*COS_AMP(I)
      SIN_AMP(I) = 1.0_RK/TOTAL*SIN_AMP(I)
      LOCAL(1_IK::2_IK) = LOCAL(1_IK::2_IK)-COS_AMP(I)*DELTA(1_IK::2_IK)+SIN_AMP(I)*DELTA(2_IK::2_IK)
      LOCAL(2_IK::2_IK) = LOCAL(2_IK::2_IK)-COS_AMP(I)*DELTA(2_IK::2_IK)-SIN_AMP(I)*DELTA(1_IK::2_IK)
    END DO
    IF (FLAG == 0_IK) THEN
      COS_AMP = 2.0_RK*COS_AMP
      SIN_AMP = 2.0_RK*SIN_AMP
    END IF
    AMPLITUDE = COS_AMP**2_IK+SIN_AMP**2_IK
    DO I = 1_IK, LOOP, 1_IK
      ORDERING(I) = __MAXLOC__(AMPLITUDE, 1_IK)
      AMPLITUDE(ORDERING(I)) = 0.0_RK
    END DO
    FREQUENCY = FREQUENCY([ORDERING])
    COS_AMP = COS_AMP([ORDERING])
    SIN_AMP = SIN_AMP([ORDERING])
  END SUBROUTINE DECOMPOSITION_
  ! ############################################################################################################################# !
  ! SIGNAL DECOMPOSITION (MEMORIZATION)
  ! (SUBROUTINE) DECOMPOSITION__(<FLAG>, <RANGE_MIN>, <RANGE_MAX>, <METHOD>, <MODE>, <LENGTH>, <PAD>, <TOTAL>, <WINDOW>, <SEQUENCE>, <LOOP>, <FREQUENCY>, <COS_AMP>, <SIN_AMP>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <RANGE_MIN>            -- (IN)     (MIN) FREQUENCY RANGE (RK)
  ! <RANGE_MAX>            -- (IN)     (MAX) FREQUENCY RANGE (RK)
  ! <METHOD>               -- (IN)     FREQUENCY APPROXIMATION METHOD (IK), FREQUENCY_FFT = 0_IK, FREQUENCY_FFRFT = 1_IK, FREQUECY_PARABOLA = 2_IK
  ! <MODE>                 -- (IN)     DECOMPOSTION MODE (IK), <MODE> = DECOMPOSITION_SUBTRACT = 0 OR <MODE> = DECOMPOSITION_PEAK = 1
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <PAD>                  -- (IN)     PADDED SEQUENCE LENGTH (IK), IF PAD > LENGTH, INPUT SEQUENCE IS PADDED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <LOOP>                 -- (IN)     NUMBER OF ITERATIONS/PEAKS (IK)
  ! <FREQUENCY>            -- (OUT)    FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <COS_AMP>              -- (OUT)    COS AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <SIN_AMP>              -- (OUT)    SIN AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! void    decomposition__(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, double* sequence, int* loop, double* frequency, double* cos_amp, double* sin_amp) ;
  MODULE SUBROUTINE DECOMPOSITION__(FLAG, RANGE_MIN, RANGE_MAX, &
    METHOD, MODE, LENGTH, PAD, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP) &
    BIND(C, NAME = "decomposition__")
    INTEGER(IK), INTENT(IN):: FLAG
    REAL(RK), INTENT(IN) :: RANGE_MIN
    REAL(RK), INTENT(IN) :: RANGE_MAX
    INTEGER(IK), INTENT(IN):: METHOD
    INTEGER(IK), INTENT(IN):: MODE
    INTEGER(IK), INTENT(IN):: LENGTH
    INTEGER(IK), INTENT(IN):: PAD
    REAL(RK), INTENT(IN) :: TOTAL
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
    REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    INTEGER(IK), INTENT(IN) :: LOOP
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: FREQUENCY
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: COS_AMP
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: SIN_AMP
    REAL(RK), DIMENSION(2_IK*LENGTH) :: LOCAL
    INTEGER(IK) :: I
    REAL(RK), DIMENSION(LENGTH) :: LIST
    REAL(RK), DIMENSION(2_IK*LENGTH) :: DELTA
    INTEGER(IK), DIMENSION(LOOP) :: ORDERING
    REAL(RK), DIMENSION(LOOP) :: AMPLITUDE
    LOCAL = SEQUENCE
    LIST = TWO_PI*REAL([(I, I=1_IK, LENGTH, 1_IK)], RK)
    DO I = 1_IK, LOOP, 1_IK
      IF (MODE == DECOMPOSITION_SUBTRACT) THEN
        FREQUENCY(I) = FREQUENCY__(FLAG, RANGE_MIN, RANGE_MAX, 0_IK, METHOD, LENGTH, PAD, TOTAL, WINDOW, LOCAL)
      ELSE
        FREQUENCY(I) = FREQUENCY__(FLAG, RANGE_MIN, RANGE_MAX, I, METHOD, LENGTH, PAD, TOTAL, WINDOW, SEQUENCE)
      END IF
      DELTA(1_IK::2_IK) = +COS(FREQUENCY(I)*LIST)
      DELTA(2_IK::2_IK) = -SIN(FREQUENCY(I)*LIST)
      COS_AMP(I) = SUM(LOCAL(1_IK::2_IK)*DELTA(1_IK::2_IK)*WINDOW)+SUM(LOCAL(2_IK::2_IK)*DELTA(2_IK::2_IK)*WINDOW)
      SIN_AMP(I) = SUM(LOCAL(2_IK::2_IK)*DELTA(1_IK::2_IK)*WINDOW)-SUM(LOCAL(1_IK::2_IK)*DELTA(2_IK::2_IK)*WINDOW)
      COS_AMP(I) = 1.0_RK/TOTAL*COS_AMP(I)
      SIN_AMP(I) = 1.0_RK/TOTAL*SIN_AMP(I)
      LOCAL(1_IK::2_IK) = LOCAL(1_IK::2_IK)-COS_AMP(I)*DELTA(1_IK::2_IK)+SIN_AMP(I)*DELTA(2_IK::2_IK)
      LOCAL(2_IK::2_IK) = LOCAL(2_IK::2_IK)-COS_AMP(I)*DELTA(2_IK::2_IK)-SIN_AMP(I)*DELTA(1_IK::2_IK)
    END DO
    IF (FLAG == 0_IK) THEN
      COS_AMP = 2.0_RK*COS_AMP
      SIN_AMP = 2.0_RK*SIN_AMP
    END IF
    AMPLITUDE = COS_AMP**2_IK+SIN_AMP**2_IK
    DO I = 1_IK, LOOP, 1_IK
      ORDERING(I) = __MAXLOC__(AMPLITUDE, 1_IK)
      AMPLITUDE(ORDERING(I)) = 0.0_RK
    END DO
    FREQUENCY = FREQUENCY([ORDERING])
    COS_AMP = COS_AMP([ORDERING])
    SIN_AMP = SIN_AMP([ORDERING])
  END SUBROUTINE DECOMPOSITION__
  ! ############################################################################################################################# !
  ! FREQUENCY LIST (PERFORM DECOMPOSITION AND RETURN LIST OF FREQUENCIES)
  ! (SUBROUTINE) FREQUENCY_LIST_(<FLAG>, <RANGE_MIN>, <RANGE_MAX>, <METHOD>, <MODE>, <LENGTH>, <PAD>, <TOTAL>, <WINDOW>, <SEQUENCE>, <LOOP>, <FREQUENCY>, <COS_AMP>, <SIN_AMP>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <RANGE_MIN>            -- (IN)     (MIN) FREQUENCY RANGE (RK)
  ! <RANGE_MAX>            -- (IN)     (MAX) FREQUENCY RANGE (RK)
  ! <METHOD>               -- (IN)     FREQUENCY APPROXIMATION METHOD (IK), FREQUENCY_FFT = 0_IK, FREQUENCY_FFRFT = 1_IK, FREQUECY_PARABOLA = 2_IK
  ! <MODE>                 -- (IN)     DECOMPOSTION MODE (IK), <MODE> = DECOMPOSITION_SUBTRACT = 0 OR <MODE> = DECOMPOSITION_PEAK = 1
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <PAD>                  -- (IN)     PADDED SEQUENCE LENGTH (IK), IF PAD > LENGTH, INPUT SEQUENCE IS PADDED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <LOOP>                 -- (IN)     NUMBER OF ITERATIONS/PEAKS (IK)
  ! <FREQUENCY>            -- (OUT)    FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! void    frequency_list_(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, double* sequence, int* loop, double* frequency) ;
  MODULE SUBROUTINE FREQUENCY_LIST_(FLAG, RANGE_MIN, RANGE_MAX, &
    METHOD, MODE, LENGTH, PAD, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY) &
    BIND(C, NAME = "frequency_list_")
    INTEGER(IK), INTENT(IN):: FLAG
    REAL(RK), INTENT(IN) :: RANGE_MIN
    REAL(RK), INTENT(IN) :: RANGE_MAX
    INTEGER(IK), INTENT(IN):: METHOD
    INTEGER(IK), INTENT(IN):: MODE
    INTEGER(IK), INTENT(IN):: LENGTH
    INTEGER(IK), INTENT(IN):: PAD
    REAL(RK), INTENT(IN) :: TOTAL
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
    REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    INTEGER(IK), INTENT(IN) :: LOOP
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: FREQUENCY
    REAL(RK), DIMENSION(LOOP) :: COS_AMP
    REAL(RK), DIMENSION(LOOP) :: SIN_AMP
    CALL DECOMPOSITION_(FLAG, RANGE_MIN, RANGE_MAX, &
    METHOD, MODE, LENGTH, PAD, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP)
  END SUBROUTINE FREQUENCY_LIST_
  ! ############################################################################################################################# !
  ! FREQUENCY LIST (PERFORM DECOMPOSITION AND RETURN LIST OF FREQUENCIES) (MEMORIZATION)
  ! (SUBROUTINE) FREQUENCY_LIST__(<FLAG>, <RANGE_MIN>, <RANGE_MAX>, <METHOD>, <MODE>, <LENGTH>, <PAD>, <TOTAL>, <WINDOW>, <SEQUENCE>, <LOOP>, <FREQUENCY>, <COS_AMP>, <SIN_AMP>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <RANGE_MIN>            -- (IN)     (MIN) FREQUENCY RANGE (RK)
  ! <RANGE_MAX>            -- (IN)     (MAX) FREQUENCY RANGE (RK)
  ! <METHOD>               -- (IN)     FREQUENCY APPROXIMATION METHOD (IK), FREQUENCY_FFT = 0_IK, FREQUENCY_FFRFT = 1_IK, FREQUECY_PARABOLA = 2_IK
  ! <MODE>                 -- (IN)     DECOMPOSTION MODE (IK), <MODE> = DECOMPOSITION_SUBTRACT = 0 OR <MODE> = DECOMPOSITION_PEAK = 1
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <PAD>                  -- (IN)     PADDED SEQUENCE LENGTH (IK), IF PAD > LENGTH, INPUT SEQUENCE IS PADDED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <LOOP>                 -- (IN)     NUMBER OF ITERATIONS/PEAKS (IK)
  ! <FREQUENCY>            -- (OUT)    FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! void    frequency_list__(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, double* sequence, int* loop, double* frequency) ;
  MODULE SUBROUTINE FREQUENCY_LIST__(FLAG, RANGE_MIN, RANGE_MAX, &
    METHOD, MODE, LENGTH, PAD, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY) &
    BIND(C, NAME = "frequency_list__")
    INTEGER(IK), INTENT(IN):: FLAG
    REAL(RK), INTENT(IN) :: RANGE_MIN
    REAL(RK), INTENT(IN) :: RANGE_MAX
    INTEGER(IK), INTENT(IN):: METHOD
    INTEGER(IK), INTENT(IN):: MODE
    INTEGER(IK), INTENT(IN):: LENGTH
    INTEGER(IK), INTENT(IN):: PAD
    REAL(RK), INTENT(IN) :: TOTAL
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
    REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    INTEGER(IK), INTENT(IN) :: LOOP
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: FREQUENCY
    REAL(RK), DIMENSION(LOOP) :: COS_AMP
    REAL(RK), DIMENSION(LOOP) :: SIN_AMP
    CALL DECOMPOSITION__(FLAG, RANGE_MIN, RANGE_MAX, &
    METHOD, MODE, LENGTH, PAD, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP)
  END SUBROUTINE FREQUENCY_LIST__
  ! ############################################################################################################################# !
  ! AMPLITUDE LIST (COMPUTE AMPLITUDES FOR LIST OF GIVEN FREQUENCIES)
  ! (SUBROUTINE) AMPLITUDE_LIST_(<FLAG>, <LENGTH>, <TOTAL>, <WINDOW>, <SEQUENCE>, <LOOP>, <FREQUENCY>, <COS_AMP>, <SIN_AMP>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <LOOP>                 -- (IN)     NUMBER OF ITERATIONS (IK)
  ! <FREQUENCY>            -- (IN)     FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <COS_AMP>              -- (OUT)    COS AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <SIN_AMP>              -- (OUT)    SIN AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! void    amplitude_list_(int* flag, int* length, double* total, double* window, double* sequence, int* loop, double* frequency, double* cos_amp, double* sin_amp) ;
  MODULE SUBROUTINE AMPLITUDE_LIST_(FLAG, LENGTH, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP) &
    BIND(C, NAME = "amplitude_list_")
    INTEGER(IK), INTENT(IN):: FLAG
    INTEGER(IK), INTENT(IN):: LENGTH
    REAL(RK), INTENT(IN) :: TOTAL
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
    REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    INTEGER(IK), INTENT(IN) :: LOOP
    REAL(RK), DIMENSION(LOOP), INTENT(IN) :: FREQUENCY
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: COS_AMP
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: SIN_AMP
    INTEGER(IK) :: I
    REAL(RK), DIMENSION(2_IK*LENGTH) :: LOCAL
    REAL(RK), DIMENSION(LENGTH) :: LIST
    REAL(RK), DIMENSION(2_IK*LENGTH) :: DELTA
    CALL REMOVE_WINDOW_MEAN_(LENGTH, TOTAL, WINDOW, SEQUENCE, LOCAL)
    LIST = TWO_PI*REAL([(I, I=1_IK, LENGTH, 1_IK)], RK)
    DO I = 1_IK, LOOP, 1_IK
      DELTA(1_IK::2_IK) = +COS(FREQUENCY(I)*LIST)
      DELTA(2_IK::2_IK) = -SIN(FREQUENCY(I)*LIST)
      COS_AMP(I) = SUM(LOCAL(1_IK::2_IK)*DELTA(1_IK::2_IK)*WINDOW)+SUM(LOCAL(2_IK::2_IK)*DELTA(2_IK::2_IK)*WINDOW)
      SIN_AMP(I) = SUM(LOCAL(2_IK::2_IK)*DELTA(1_IK::2_IK)*WINDOW)-SUM(LOCAL(1_IK::2_IK)*DELTA(2_IK::2_IK)*WINDOW)
      COS_AMP(I) = 1.0_RK/REAL(LENGTH, RK)*COS_AMP(I)
      SIN_AMP(I) = 1.0_RK/REAL(LENGTH, RK)*SIN_AMP(I)
      LOCAL(1_IK::2_IK) = LOCAL(1_IK::2_IK)-COS_AMP(I)*DELTA(1_IK::2_IK)+SIN_AMP(I)*DELTA(2_IK::2_IK)
      LOCAL(2_IK::2_IK) = LOCAL(2_IK::2_IK)-COS_AMP(I)*DELTA(2_IK::2_IK)-SIN_AMP(I)*DELTA(1_IK::2_IK)
    END DO
    IF (FLAG == 0_IK) THEN
      COS_AMP = 2.0_RK*COS_AMP
      SIN_AMP = 2.0_RK*SIN_AMP
    END IF
  END SUBROUTINE AMPLITUDE_LIST_
  ! ############################################################################################################################# !
  ! FREQUENCY CORRECTION
  ! (SUBROUTINE) FREQUENCY_CORRECTION_(<FLAG>, <RANGE_MIN>, <RANGE_MAX>, <METHOD>, <MODE>, <LENGTH>, <PAD>, <TOTAL>, <WINDOW>, <LOOP>, <FREQUENCY>, <COS_AMP>, <SIN_AMP>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <RANGE_MIN>            -- (IN)     (MIN) FREQUENCY RANGE (RK)
  ! <RANGE_MAX>            -- (IN)     (MAX) FREQUENCY RANGE (RK)
  ! <METHOD>               -- (IN)     FREQUENCY APPROXIMATION METHOD (IK), FREQUENCY_FFT = 0_IK, FREQUENCY_FFRFT = 1_IK, FREQUECY_PARABOLA = 2_IK
  ! <MODE>                 -- (IN)     DECOMPOSTION MODE (IK), <MODE> = DECOMPOSITION_SUBTRACT = 0 OR <MODE> = DECOMPOSITION_PEAK = 1
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <PAD>                  -- (IN)     PADDED SEQUENCE LENGTH (IK), IF PAD > LENGTH, INPUT SEQUENCE IS PADDED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <LOOP>                 -- (IN)     NUMBER OF ITERATIONS/PEAKS (IK)
  ! <FREQUENCY>            -- (INOUT)  FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <COS_AMP>              -- (INOUT)  COS AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <SIN_AMP>              -- (INOUT)  SIN AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! void    frequency_correction_(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, int* loop, double* frequency, double* cos_amp, double* sin_amp) ;
  MODULE SUBROUTINE FREQUENCY_CORRECTION_(FLAG, RANGE_MIN, RANGE_MAX, &
    METHOD, MODE, LENGTH, PAD, TOTAL, WINDOW, LOOP, FREQUENCY, COS_AMP, SIN_AMP) &
    BIND(C, NAME = "frequency_correction_")
    INTEGER(IK), INTENT(IN):: FLAG
    REAL(RK), INTENT(IN) :: RANGE_MIN
    REAL(RK), INTENT(IN) :: RANGE_MAX
    INTEGER(IK), INTENT(IN):: METHOD
    INTEGER(IK), INTENT(IN):: MODE
    INTEGER(IK), INTENT(IN):: LENGTH
    INTEGER(IK), INTENT(IN):: PAD
    REAL(RK), INTENT(IN) :: TOTAL
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
    INTEGER(IK), INTENT(IN) :: LOOP
    REAL(RK), DIMENSION(LOOP), INTENT(INOUT) :: FREQUENCY
    REAL(RK), DIMENSION(LOOP), INTENT(INOUT) :: COS_AMP
    REAL(RK), DIMENSION(LOOP), INTENT(INOUT) :: SIN_AMP
    REAL(RK), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    REAL(RK), DIMENSION(LOOP) :: F, C, S
    F = FREQUENCY
    C = COS_AMP
    S = SIN_AMP
    CALL GENERATE_SIGNAL_(FLAG, LENGTH, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP)
    CALL DECOMPOSITION_(FLAG, RANGE_MIN, RANGE_MAX, &
    METHOD, MODE, LENGTH, PAD, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP)
    FREQUENCY = 2.0_RK*F - FREQUENCY
    COS_AMP = 2.0_RK*C - COS_AMP
    SIN_AMP = 2.0_RK*S - SIN_AMP
  END SUBROUTINE FREQUENCY_CORRECTION_
  ! ############################################################################################################################# !
  ! FREQUENCY CORRECTION
  ! (SUBROUTINE) FREQUENCY_CORRECTION_(<FLAG>, <RANGE_MIN>, <RANGE_MAX>, <METHOD>, <MODE>, <LENGTH>, <PAD>, <TOTAL>, <WINDOW>, <LOOP>, <FREQUENCY>, <COS_AMP>, <SIN_AMP>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <RANGE_MIN>            -- (IN)     (MIN) FREQUENCY RANGE (RK)
  ! <RANGE_MAX>            -- (IN)     (MAX) FREQUENCY RANGE (RK)
  ! <METHOD>               -- (IN)     FREQUENCY APPROXIMATION METHOD (IK), FREQUENCY_FFT = 0_IK, FREQUENCY_FFRFT = 1_IK, FREQUECY_PARABOLA = 2_IK
  ! <MODE>                 -- (IN)     DECOMPOSTION MODE (IK), <MODE> = DECOMPOSITION_SUBTRACT = 0 OR <MODE> = DECOMPOSITION_PEAK = 1
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <PAD>                  -- (IN)     PADDED SEQUENCE LENGTH (IK), IF PAD > LENGTH, INPUT SEQUENCE IS PADDED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <LOOP>                 -- (IN)     NUMBER OF ITERATIONS/PEAKS (IK)
  ! <FREQUENCY>            -- (INOUT)  FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <COS_AMP>              -- (INOUT)  COS AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <SIN_AMP>              -- (INOUT)  SIN AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! void    frequency_correction__(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, int* loop, double* frequency, double* cos_amp, double* sin_amp) ;
  MODULE SUBROUTINE FREQUENCY_CORRECTION__(FLAG, RANGE_MIN, RANGE_MAX, &
    METHOD, MODE, LENGTH, PAD, TOTAL, WINDOW, LOOP, FREQUENCY, COS_AMP, SIN_AMP) &
    BIND(C, NAME = "frequency_correction__")
    INTEGER(IK), INTENT(IN):: FLAG
    REAL(RK), INTENT(IN) :: RANGE_MIN
    REAL(RK), INTENT(IN) :: RANGE_MAX
    INTEGER(IK), INTENT(IN):: METHOD
    INTEGER(IK), INTENT(IN):: MODE
    INTEGER(IK), INTENT(IN):: LENGTH
    INTEGER(IK), INTENT(IN):: PAD
    REAL(RK), INTENT(IN) :: TOTAL
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
    INTEGER(IK), INTENT(IN) :: LOOP
    REAL(RK), DIMENSION(LOOP), INTENT(INOUT) :: FREQUENCY
    REAL(RK), DIMENSION(LOOP), INTENT(INOUT) :: COS_AMP
    REAL(RK), DIMENSION(LOOP), INTENT(INOUT) :: SIN_AMP
    REAL(RK), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    REAL(RK), DIMENSION(LOOP) :: F, C, S
    F = FREQUENCY
    C = COS_AMP
    S = SIN_AMP
    CALL GENERATE_SIGNAL_(FLAG, LENGTH, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP)
    CALL DECOMPOSITION__(FLAG, RANGE_MIN, RANGE_MAX, &
    METHOD, MODE, LENGTH, PAD, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP)
    FREQUENCY = 2.0_RK*F - FREQUENCY
    COS_AMP = 2.0_RK*C - COS_AMP
    SIN_AMP = 2.0_RK*S - SIN_AMP
  END SUBROUTINE FREQUENCY_CORRECTION__
  ! ############################################################################################################################# !
END SUBMODULE DECOMPOSITION
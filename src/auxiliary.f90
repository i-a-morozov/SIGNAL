
#include "signal.inc"

SUBMODULE (SIGNAL) AUXILIARY
  IMPLICIT NONE
  CONTAINS
  ! ############################################################################################################################# !
  ! FACTORIAL
  ! (FUNCTION) FACTORIAL_(<NUMBER>)
  ! <NUMBER>               -- (IN)     NUMBER (IK)
  ! <FACTORIAL_>           -- (OUT)    FACTORIAL OF <N> (RK REAL)
  MODULE REAL(RK) FUNCTION FACTORIAL_(NUMBER)
    INTEGER(IK), INTENT(IN) :: NUMBER
    INTEGER(IK) :: I
    FACTORIAL_ = 1.0_RK
    DO I = 1_IK, NUMBER, 1_IK
      FACTORIAL_ = REAL(I, RK)*FACTORIAL_
    END DO
  END FUNCTION FACTORIAL_
  ! ############################################################################################################################# !
  ! GAMMA REGULARIZED
  ! (FUNCTION) GAMMA_REGULARIZED_(<A>, <X>, <Y>)
  ! <A>                    -- (IN)     A (RK)
  ! <X>                    -- (IN)     X (RK)
  ! <Y>                    -- (IN)     Y (RK)
  ! <FACTORIAL_>           -- (OUT)    GAMMA REGULARIZED
  MODULE REAL(RK) FUNCTION GAMMA_REGULARIZED_(A, X, Y)
    REAL(RK), INTENT(IN) :: A
    REAL(RK), INTENT(IN) :: X
    REAL(RK), INTENT(IN) :: Y
    GAMMA_REGULARIZED_ = (GAMMA_(A, X)-GAMMA_(A, Y))/GAMMA_(A)
  END FUNCTION GAMMA_REGULARIZED_
  ! ############################################################################################################################# !
  ! MINLOC
  ! (FUNCTION) MINLOC_(<SEQUENCE>)
  ! <SEQUENCE>             -- (IN)     SEQUENCE (RK ARRAY)
  ! <MINLOC_>              -- (OUT)    MINIMUM LOCATION (IK)
  MODULE INTEGER(IK) FUNCTION MINLOC_(SEQUENCE, EMPTY)
    REAL(RK), DIMENSION(:), CONTIGUOUS, INTENT(IN) :: SEQUENCE
    INTEGER, INTENT(IN) :: EMPTY
    REAL(RK) :: ELEMENT
    INTEGER(IK) :: I
    MINLOC_ = 0_IK
    ELEMENT = MINVAL(SEQUENCE,EMPTY)
    DO I = 1_IK, SIZE(SEQUENCE, KIND = IK), 1_IK
      IF (SEQUENCE(I) == ELEMENT) THEN
        MINLOC_ = I
        RETURN
      END IF
    END DO
  END FUNCTION MINLOC_
  ! ############################################################################################################################# !
  ! MAXLOC
  ! (FUNCTION) MAXLOC_(<SEQUENCE>)
  ! <SEQUENCE>             -- (IN)     SEQUENCE (RK ARRAY)
  ! <MAXLOC_>              -- (OUT)    MAXIMUM LOCATION (IK)
  MODULE INTEGER(IK) FUNCTION MAXLOC_(SEQUENCE, EMPTY)
    REAL(RK), DIMENSION(:), CONTIGUOUS, INTENT(IN) :: SEQUENCE
    INTEGER, INTENT(IN) :: EMPTY
    REAL(RK) :: ELEMENT
    INTEGER(IK) :: I
    MAXLOC_ = 0_IK
    ELEMENT = MAXVAL(SEQUENCE, EMPTY)
    DO I = 1_IK, SIZE(SEQUENCE, KIND = IK), 1_IK
      IF (SEQUENCE(I) == ELEMENT) THEN
        MAXLOC_ = I
        RETURN
      END IF
    END DO
  END FUNCTION MAXLOC_
  ! ############################################################################################################################# !
  ! SORT (BUBBLE, DESCENDING)
  ! (SUBROUTINE) SORT_BUBBLE_(<LENGTH>, <SEQUENCE>, <FST>, <LST>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (IN)     (UNSORTED) SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (OUT)    (SORTED, DESCENDING) SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  MODULE SUBROUTINE SORT_BUBBLE_(LENGTH, SEQUENCE, FST, LST)
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(LENGTH), INTENT(INOUT) :: SEQUENCE
    INTEGER(IK), INTENT(IN) :: FST
    INTEGER(IK), INTENT(IN) :: LST
    INTEGER(IK) :: I, J
    LOGICAL :: SWAPPED
    IF (FST > LST) RETURN
    DO I = LENGTH-1_IK, 1_IK, -1_IK
      SWAPPED = .FALSE.
      DO J = 1_IK, I, 1_IK
        IF (SEQUENCE(J) < SEQUENCE(J+1_IK)) THEN
          SEQUENCE(J:J+1_IK) = SEQUENCE([J+1_IK, J])
          SWAPPED = .TRUE.
        END IF
      END DO
      IF (.NOT. SWAPPED) EXIT
    END DO
  END SUBROUTINE SORT_BUBBLE_
  ! ############################################################################################################################# !
  ! SORT (QUICK, DESCENDING)
  ! (SUBROUTINE) SORT_QUICK_(<LENGTH>, <SEQUENCE>, <FST>, <LST>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (IN)     (UNSORTED) SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (OUT)    (SORTED, DESCENDING) SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  MODULE RECURSIVE SUBROUTINE SORT_QUICK_(LENGTH, SEQUENCE, FST, LST)
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(:), INTENT(INOUT) :: SEQUENCE
    INTEGER(IK), INTENT(IN) :: FST
    INTEGER(IK), INTENT(IN) :: LST
    REAL(RK) :: PIVOT, COPY
    INTEGER(IK) :: I, J
    PIVOT = SEQUENCE((FST+LST)/2_IK)
    I = FST
    J = LST
    DO
       DO WHILE(SEQUENCE(I) > PIVOT)
          I = I+1_IK
       END DO
       DO WHILE(PIVOT > SEQUENCE(J))
          J = J-1_IK
       END DO
       IF (I >= J) EXIT
       COPY = SEQUENCE(I)
       SEQUENCE(I) = SEQUENCE(J)
       SEQUENCE(J) = COPY
       I = I+1_IK
       J = J-1_IK
    END DO
    IF (FST < I-1_IK) CALL SORT_QUICK_(LENGTH, SEQUENCE, FST, I-1_IK)
    IF (J+1_IK < LST) CALL SORT_QUICK_(LENGTH, SEQUENCE, J+1_IK, LST)
  END SUBROUTINE SORT_QUICK_
  ! ############################################################################################################################# !
  ! GENERATE HARMONIC SIGNAL
  ! (SUBROUTINE) GENERATE_SIGNAL_(<FLAG>, <LENGTH>, <SEQUENCE>, <LOOP>, <FREQUENCY>, <MEAN>, <COS_AMP>, <SIN_AMP>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX SEQUENCE
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK), POWER OF TWO
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  ! <LOOP>                 -- (IN)     NUMBER OF HARMONICS (IK)
  ! <FREQUENCY>            -- (IN)     FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <COS_AMP>              -- (OUT)    COS AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <SIN_AMP>              -- (OUT)    SIN AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! void    generate_signal_(int*, int*, double*, int*, double*, double*, double*, double*) ;
  MODULE SUBROUTINE GENERATE_SIGNAL_(FLAG, LENGTH, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP) &
    BIND(C, NAME = "generate_signal_")
    INTEGER(IK), INTENT(IN) :: FLAG
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(OUT) :: SEQUENCE
    INTEGER(IK), INTENT(IN) :: LOOP
    REAL(RK), DIMENSION(LOOP), INTENT(IN) :: FREQUENCY
    REAL(RK), DIMENSION(LOOP), INTENT(IN) :: COS_AMP
    REAL(RK), DIMENSION(LOOP), INTENT(IN) :: SIN_AMP
    REAL(RK), DIMENSION(LENGTH) :: RANGE
    REAL(RK), DIMENSION(LENGTH) :: FC, FS
    INTEGER(IK) :: I
    RANGE = TWO_PI*REAL([(I, I = 1_IK, LENGTH, 1_IK)], RK)
    SEQUENCE = 0.0_RK
    DO I = 1_IK, LOOP, 1_IK
      FC = COS(FREQUENCY(I)*RANGE)
      FS = SIN(FREQUENCY(I)*RANGE)
      IF (FLAG == 0_IK) THEN
        SEQUENCE(1_IK::2_IK) = SEQUENCE(1_IK::2_IK)+COS_AMP(I)*FC+SIN_AMP(I)*FS
      ELSE
        SEQUENCE(1_IK::2_IK) = SEQUENCE(1_IK::2_IK)+COS_AMP(I)*FC+SIN_AMP(I)*FS
        SEQUENCE(2_IK::2_IK) = SEQUENCE(2_IK::2_IK)+SIN_AMP(I)*FC-COS_AMP(I)*FS
      END IF  
    END DO
  END SUBROUTINE GENERATE_SIGNAL_
  ! ############################################################################################################################# !
END SUBMODULE AUXILIARY
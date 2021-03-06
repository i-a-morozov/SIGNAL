! EXAMPLE-08: FREQUENCY ESTIMATION (SIGNAL WITH NOISE)
PROGRAM EXAMPLE

  USE SIGNAL

  IMPLICIT NONE

  INTEGER(IK), PARAMETER           :: LENGTH = 2_IK**10_IK  ! INPUT SIGNAL LENGTH
  INTEGER(IK)                      :: FLAG                  ! COMPLEX FLAG (0/1)
  REAL(RK)                         :: RANGE_MIN             ! (MIN) FREQUENCY RANGE
  REAL(RK)                         :: RANGE_MAX             ! (MAX) FREQUENCY RANGE
  INTEGER(IK)                      :: PEAK                  ! FOURIE SPECTRA PEAK ID
  INTEGER(IK)                      :: METHOD                ! FREQUENCY ESTIMATION METHOD
  INTEGER(IK)                      :: ORDER                 ! COSINE WINDOW ORDER
  REAL(RK), DIMENSION(LENGTH)      :: WINDOW                ! COSINE WINDOW DATA
  REAL(RK)                         :: TOTAL                 ! SUM OF WINDOW ELEMENTS
  REAL(RK)                         :: FREQUENCY             ! EXACT SIGNAL FREQUENCY
  REAL(RK), DIMENSION(LENGTH)      :: SIGNAL_R, SIGNAL_I    ! SIGNAL REAL AND COMPLEX PARTS
  REAL(RK), DIMENSION(2_IK*LENGTH) :: SIGNAL                ! INPUT SIGNAL, [..., SR_I, SI_I, ...]
  REAL(RK)                         :: RESULT                ! FREQUENCY ESTIMATION
  REAL(RK), DIMENSION(LENGTH/2_IK) :: N, M                  ! NOISE
  INTEGER(IK), PARAMETER           :: LIMIT = 4_IK          ! NUMBER OF SINGULAR VALUES TO KEEP
  REAL(RK), DIMENSION(LIMIT)       :: LIST                  ! SVD LIST
  REAL(RK)                         :: COS_AMP, SIN_AMP, AMP

  INTEGER(IK)                      :: I

  FLAG = 0_IK
  RANGE_MIN = 0.00_RK
  RANGE_MAX = 0.49_RK

  ! SET SIGNAL EXACT FREQUENCY, REAL AND IMAGINARY PARTS
  FREQUENCY = 0.123456789_RK
  DO I = 1_IK, LENGTH, 1_IK
    SIGNAL_R(I) = 1.0_RK*SIN(TWO_PI*FREQUENCY*I)+0.05*SIN(TWO_PI*2.0_RK*FREQUENCY*I)
    SIGNAL_I(I) = 0.0_RK
  END DO

  ! ADD NOISE
  CALL RANDOM_NUMBER(N)
  CALL RANDOM_NUMBER(M)

  SIGNAL_R(1_IK:LENGTH:2_IK) = SIGNAL_R(1_IK:LENGTH:2_IK) + 0.05*SQRT(-2.0_RK*LOG(N))*COS(TWO_PI*M)
  SIGNAL_R(2_IK:LENGTH:2_IK) = SIGNAL_R(2_IK:LENGTH:2_IK) + 0.05*SQRT(-2.0_RK*LOG(N))*SIN(TWO_PI*M)

  ! FORMAT TEST SIGNAL
  SIGNAL = 0.0_RK
  IF(FLAG == 0_IK) THEN
    CALL CONVERT_(LENGTH, SIGNAL_R, SIGNAL)
  ELSE
    CALL CONVERT_(LENGTH, SIGNAL_R, SIGNAL_I, SIGNAL)
  END IF

  ! SET WINDOW AND WINDOW SUM
  ORDER = 1_IK
  CALL WINDOW_(LENGTH, ORDER, WINDOW)
  TOTAL = SUM(WINDOW)

  ! FREQUENCY ESTIMATION METHOD
  METHOD = FREQUENCY_PARABOLA

  ! ESTIMATE FREQUENCY
  BLOCK
    WRITE(*, '(A)') "FILTER(-)"
    PEAK = +1_IK
    RESULT = FREQUENCY_(FLAG, RANGE_MIN, RANGE_MAX, PEAK, METHOD, LENGTH, LENGTH, TOTAL, WINDOW, SIGNAL)
    CALL AMPLITUDE_(FLAG, LENGTH, TOTAL, WINDOW, SIGNAL, RESULT, COS_AMP, SIN_AMP, AMP)
    WRITE(*, '(A)') "1ST"
    WRITE(*,'(A,E32.16,A,E32.16,A,2E32.16)')" FREQUENCY",RESULT," ERROR",ABS(RESULT-1.0_RK*FREQUENCY)," AMPLITUDE",COS_AMP,SIN_AMP
    PEAK = +2_IK
    RESULT = FREQUENCY_(FLAG, RANGE_MIN, RANGE_MAX, PEAK, METHOD, LENGTH, LENGTH, TOTAL, WINDOW, SIGNAL)
    CALL AMPLITUDE_(FLAG, LENGTH, TOTAL, WINDOW, SIGNAL, RESULT, COS_AMP, SIN_AMP, AMP)
    WRITE(*, '(A)') "2ND"
    WRITE(*,'(A,E32.16,A,E32.16,A,2E32.16)')" FREQUENCY",RESULT," ERROR",ABS(RESULT-2.0_RK*FREQUENCY)," AMPLITUDE",COS_AMP,SIN_AMP
    WRITE(*, *)
  END BLOCK

  ! FILTER AND ESTIMATE FREQUENCY
  BLOCK
    CALL FILTER_(LENGTH, SIGNAL_R, LIMIT, LIST)
    CALL CONVERT_(LENGTH, SIGNAL_R, SIGNAL)
    WRITE(*, '(A)') "FILTER(+)"
    PEAK = +1_IK
    RESULT = FREQUENCY_(FLAG, RANGE_MIN, RANGE_MAX, PEAK, METHOD, LENGTH, LENGTH, TOTAL, WINDOW, SIGNAL)
    CALL AMPLITUDE_(FLAG, LENGTH, TOTAL, WINDOW, SIGNAL, RESULT, COS_AMP, SIN_AMP, AMP)
    WRITE(*, '(A)') "1ST"
    WRITE(*,'(A,E32.16,A,E32.16,A,2E32.16)')" FREQUENCY",RESULT," ERROR",ABS(RESULT-1.0_RK*FREQUENCY)," AMPLITUDE",COS_AMP,SIN_AMP
    PEAK = +2_IK
    RESULT = FREQUENCY_(FLAG, RANGE_MIN, RANGE_MAX, PEAK, METHOD, LENGTH, LENGTH, TOTAL, WINDOW, SIGNAL)
    CALL AMPLITUDE_(FLAG, LENGTH, TOTAL, WINDOW, SIGNAL, RESULT, COS_AMP, SIN_AMP, AMP)
    WRITE(*, '(A)') "2ND"
    WRITE(*,'(A,E32.16,A,E32.16,A,2E32.16)')" FREQUENCY",RESULT," ERROR",ABS(RESULT-2.0_RK*FREQUENCY)," AMPLITUDE",COS_AMP,SIN_AMP
    WRITE(*, *)
  END BLOCK

  ! IN GENERAL APPLICATION OF FILTER MIGHT IMPROVE FREQUENCY ESTIMATION ACCURACY (STATISTICALLY)

END PROGRAM EXAMPLE

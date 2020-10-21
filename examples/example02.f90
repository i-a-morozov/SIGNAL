! ulimit -s unlimited
! EXAMPLE-02: FREQUENCY ESTIMATION (OPENMP, __FFT__ DIRECTIVE SHOULD BE FFT_RADIX_*)
PROGRAM EXAMPLE

  USE SIGNAL

  IMPLICIT NONE

  INTEGER(IK), PARAMETER                 :: LENGTH = 2_IK**10_IK  ! INPUT SIGNAL LENGTH
  INTEGER(IK)                            :: FLAG                  ! COMPLEX FLAG (0/1)
  REAL(RK)                               :: RANGE_MIN             ! (MIN) FREQUENCY RANGE
  REAL(RK)                               :: RANGE_MAX             ! (MAX) FREQUENCY RANGE
  INTEGER(IK)                            :: PEAK                  ! FOURIE SPECTRA PEAK ID
  INTEGER(IK)                            :: METHOD                ! FREQUENCY ESTIMATION METHOD
  INTEGER(IK)                            :: ORDER                 ! COSINE WINDOW ORDER
  REAL(RK), DIMENSION(LENGTH)            :: WINDOW                ! COSINE WINDOW DATA
  REAL(RK)                               :: TOTAL                 ! SUM OF WINDOW ELEMENTS
  REAL(RK)                               :: FREQUENCY             ! EXACT SIGNAL FREQUENCY
  REAL(RK), DIMENSION(LENGTH)            :: SIGNAL_R, SIGNAL_I    ! SIGNAL REAL AND COMPLEX PARTS
  REAL(RK), DIMENSION(2_IK*LENGTH)       :: SIGNAL                ! INPUT SIGNAL, [..., SR_I, SI_I, ...]
  INTEGER(IK), PARAMETER                 :: LIMIT = 128_IK        ! NUMBER OF SIGNALS
  REAL(RK), DIMENSION(LIMIT,2_IK*LENGTH) :: DATA                  ! MATRIX OF SIGNALS
  REAL(RK), DIMENSION(LIMIT)             :: LIST                  ! LIST OF EXACT FREQUENCIES
  REAL(RK), DIMENSION(LIMIT)             :: OUTPUT                ! ESTIMATED FREQUENCIES

  INTEGER(IK)                            :: I, J

  ! SET COMPLEX FLAG (0/1 FOR REAL/COMPLEX SIGNAL)
  FLAG = 1_IK
  RANGE_MIN = 0.00_RK
  RANGE_MAX = 0.99_RK

  ! SET TEST DATA
  DO I = 1_IK, LIMIT, 1_IK
    ! SET SIGNAL EXACT FREQUENCY, REAL AND IMAGINARY PARTS
    FREQUENCY = 0.623456789_RK + REAL(I, RK)/REAL(LIMIT, RK)/10.0_RK
    LIST(I) = FREQUENCY
    DO J = 1_IK, LENGTH, 1_IK
      SIGNAL_R(J) = +COS(TWO_PI*FREQUENCY*REAL(J, RK))
      SIGNAL_I(J) = -SIN(TWO_PI*FREQUENCY*REAL(J, RK))
    END DO
    ! FORMAT TEST SIGNAL
    IF(FLAG == 0_IK) THEN
      CALL CONVERT_(LENGTH, SIGNAL_R, SIGNAL)
    ELSE
      CALL CONVERT_(LENGTH, SIGNAL_R, SIGNAL_I, SIGNAL)
    END IF
    DATA(I,:) = SIGNAL
  END DO

  ! SET WINDOW AND WINDOW SUM
  ORDER = 2_IK
  CALL WINDOW_(LENGTH, ORDER, WINDOW)
  TOTAL = SUM(WINDOW)

  ! ESTIMATE FREQUENCY
  PEAK = 0_IK
  METHOD = FREQUENCY_PARABOLA
  !$OMP PARALLEL DO PRIVATE(SIGNAL)
  DO I = 1_IK, LIMIT, 1_IK
    OUTPUT(I) = FREQUENCY_(FLAG, RANGE_MIN, RANGE_MAX, PEAK, METHOD, LENGTH, LENGTH, TOTAL, WINDOW, DATA(I,:))
  END DO
  !$OMP END PARALLEL DO

  ! RESULT
  DO I = 1_IK, LIMIT, 1_IK
    WRITE(*, '(3E32.16)') LIST(I), OUTPUT(I), ABS(LIST(I)-OUTPUT(I))
  END DO

END PROGRAM EXAMPLE
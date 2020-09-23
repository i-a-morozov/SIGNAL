! gfortran -c -cpp -std=f2018 -Wall -pedantic -O3 -ffast-math -march=native -Wno-unused-function signal.f90
! ar rcs libsignal.a signal.o
! gfortran -o example01 -cpp -std=f2018 -Wall -pedantic -O3 -ffast-math -march=native -L. example01.f90 -lsignal  -llapack -lblas -lm -lfftw3 -lgsl -lgslcblas -lgfortran

! EXAMPLE-08: FREQUENCY ESTIMATION (SIGNAL WITH NOISE)
PROGRAM EXAMPLE

  USE SIGNAL

  IMPLICIT NONE

  INTEGER(IK), PARAMETER           :: LENGTH = 2_IK**10_IK  ! INPUT SIGNAL LENGTH
  INTEGER(IK)                      :: FLAG                  ! COMPLEX FLAG (0/1)
  INTEGER(IK)                      :: PEAK                  ! FOURIE SPECTRA PEAK ID
  INTEGER(IK)                      :: ORDER                 ! COSINE WINDOW ORDER
  REAL(RK), DIMENSION(LENGTH)      :: WINDOW                ! COSINE WINDOW DATA
  REAL(RK)                         :: TOTAL                 ! SUM OF WINDOW ELEMENTS
  REAL(RK)                         :: FREQUENCY             ! EXACT SIGNAL FREQUENCY
  REAL(RK), DIMENSION(LENGTH)      :: SIGNAL_R, SIGNAL_I    ! SIGNAL REAL AND COMPLEX PARTS
  REAL(RK), DIMENSION(2_IK*LENGTH) :: SIGNAL                ! INPUT SIGNAL, [..., SR_I, SI_I, ...]
  REAL(RK)                         :: RESULT                ! FREQUENCY ESTIMATION
  REAL(RK), DIMENSION(LENGTH)      :: NOISE                 ! NOISE
  INTEGER(IK), PARAMETER           :: LIMIT = 6_IK          ! NUMBER OF SINGULAR VALUES TO KEEP
  REAL(RK), DIMENSION(LIMIT)       :: LIST                  ! SVD LIST
  REAL(RK)                         :: COS_AMP, SIN_AMP, AMP

  INTEGER(IK)                      :: I

  FLAG = 0_IK

  ! SET SIGNAL EXACT FREQUENCY, REAL AND IMAGINARY PARTS
  FREQUENCY = 0.123456789_RK
  DO I = 1_IK, LENGTH, 1_IK
    SIGNAL_R(I) = 1.0_RK*SIN(TWO_PI*FREQUENCY*I)+0.01*SIN(TWO_PI*2.0_RK*FREQUENCY*I)
    SIGNAL_I(I) = 0.0_RK
  END DO

  ! ADD NOISE
  CALL RANDOM_NUMBER(NOISE)
  NOISE = 0.05_RK*2.0_RK*(0.5_RK-NOISE)
  SIGNAL_R = SIGNAL_R + NOISE

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

  ! ESTIMATE FREQUENCY
  WRITE(*, *) "FILTER(-)"
  PEAK = +1_IK
  RESULT = FREQUENCY_(FLAG, PEAK, LENGTH, TOTAL, WINDOW, SIGNAL)
  CALL AMPLITUDE_(LENGTH, TOTAL, WINDOW, SIGNAL, RESULT, COS_AMP, SIN_AMP, AMP)
  WRITE(*, *) "1ST"
  WRITE(*, *) "FREQUENCY", RESULT, "ERROR", ABS(RESULT-FREQUENCY), "AMPLITUDE", 2.0_RK*COS_AMP, 2.0_RK*SIN_AMP
  PEAK = +2_IK
  RESULT = FREQUENCY_(FLAG, PEAK, LENGTH, TOTAL, WINDOW, SIGNAL)
  CALL AMPLITUDE_(LENGTH, TOTAL, WINDOW, SIGNAL, RESULT, COS_AMP, SIN_AMP, AMP)
  WRITE(*, *) "2ND"
  WRITE(*, *) "FREQUENCY", RESULT, "ERROR", ABS(RESULT-2.0_RK*FREQUENCY), "AMPLITUDE", 2.0_RK*COS_AMP, 2.0_RK*SIN_AMP
  WRITE(*, *)


  ! FILTER AND ESTIMATE FREQUENCY
  WRITE(*, *) "FILTER(+)"
  CALL FILTER_(LENGTH, SIGNAL_R, LIMIT, LIST)
  CALL CONVERT_(LENGTH, SIGNAL_R, SIGNAL)
  PEAK = +1_IK
  RESULT = FREQUENCY_(FLAG, PEAK, LENGTH, TOTAL, WINDOW, SIGNAL)
  CALL AMPLITUDE_(LENGTH, TOTAL, WINDOW, SIGNAL, RESULT, COS_AMP, SIN_AMP, AMP)
  WRITE(*, *) "1ST"
  WRITE(*, *) "FREQUENCY", RESULT, "ERROR", ABS(RESULT-FREQUENCY), "AMPLITUDE", 2.0_RK*COS_AMP, 2.0_RK*SIN_AMP
  PEAK = +2_IK
  RESULT = FREQUENCY_(FLAG, PEAK, LENGTH, TOTAL, WINDOW, SIGNAL)
  CALL AMPLITUDE_(LENGTH, TOTAL, WINDOW, SIGNAL, RESULT, COS_AMP, SIN_AMP, AMP)
  WRITE(*, *) "2ND"
  WRITE(*, *) "FREQUENCY", RESULT, "ERROR", ABS(RESULT-2.0_RK*FREQUENCY), "AMPLITUDE", 2.0_RK*COS_AMP, 2.0_RK*SIN_AMP
  WRITE(*, *)

  WRITE(*, *) "SINGULAR VALUES"
  DO I = 1_IK, LIMIT, 1_IK
    WRITE(*, *) LIST(I)
  END DO

  ! IN GENERAL APPLICATION OF FILTER MIGHT NOT IMPROVE 1ST FREQUENCY ESTIMATION ACCURACY (BETTER AMPLITUDES)
  ! FOR HARMONICS WITH SMALL AMPLITUDES SOME IMPROVEMENT MIGHT BE OBTAINED WITH FILTERING

END PROGRAM EXAMPLE
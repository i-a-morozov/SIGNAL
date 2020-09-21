! gfortran -c -cpp -std=f2018 -Wall -pedantic -O3 -ffast-math -march=native -Wno-unused-function signal.f90
! ar rcs libsignal.a signal.o
! gfortran -o example01 -cpp -std=f2018 -Wall -pedantic -O3 -ffast-math -march=native -L. example01.f90 -lsignal  -llapack -lblas -lm -lfftw3 -lgsl -lgslcblas -lgfortran

! EXAMPLE-08: FREQUENCY ESTIMATION (SIGNAL WITH NOISE)
PROGRAM EXAMPLE

  USE SIGNAL

  IMPLICIT NONE

  INTEGER(IK), PARAMETER           :: LENGTH = 2**10_IK     ! INPUT SIGNAL LENGTH
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
  REAL(RK), DIMENSION(LENGTH/2_IK) :: LIST                  ! SVD LIST

  INTEGER(IK)                      :: I

  WRITE(*, *) IK_SIZE
  WRITE(*, *) RK_SIZE

  FLAG = 0_IK

  ! SET SIGNAL EXACT FREQUENCY, REAL AND IMAGINARY PARTS
  FREQUENCY = 0.123_RK
  DO I = 1_IK, LENGTH, 1_IK
    SIGNAL_R(I) = SIN(TWO_PI*FREQUENCY*I)
    SIGNAL_I(I) = 0.0_RK
  END DO

  ! ADD NOISE
  CALL RANDOM_NUMBER(NOISE)
  SIGNAL_R = SIGNAL_R + 0.05_RK*NOISE

  ! FORMAT TEST SIGNAL
  SIGNAL = 0.0_RK
  IF(FLAG == 0_IK) THEN
    CALL CONVERT_(LENGTH, SIGNAL_R, SIGNAL)
  ELSE
    CALL CONVERT_(LENGTH, SIGNAL_R, SIGNAL_I, SIGNAL)
  END IF

  ! SET WINDOW AND WINDOW SUM
  ORDER = 2_IK
  CALL WINDOW_(LENGTH, ORDER, WINDOW)
  TOTAL = SUM(WINDOW)

  ! ESTIMATE FREQUENCY
  PEAK = +1_IK
  RESULT = FREQUENCY_(FLAG, PEAK, LENGTH, TOTAL, WINDOW, SIGNAL)
  WRITE(*, *) RESULT, ABS(RESULT-FREQUENCY)

  ! FILTER AND ESTIMATE FREQUENCY
  CALL FILTER_(LENGTH, SIGNAL_R, 4_IK, LIST)
  CALL CONVERT_(LENGTH, SIGNAL_R, SIGNAL)
  RESULT = FREQUENCY_(FLAG, PEAK, LENGTH, TOTAL, WINDOW, SIGNAL)
  WRITE(*, *) RESULT, ABS(RESULT-FREQUENCY)

END PROGRAM EXAMPLE
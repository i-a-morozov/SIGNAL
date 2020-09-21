! SIGNAL, 2018-2020, I.A.MOROZOV@INP.NSK.SU

! STATIC LIB
! gfortran -c -cpp -fPIC -std=f2018 -Wall -pedantic -O3 -ffast-math -march=native -Wno-unused-function signal.f90
! ar rcs libsignal.a signal.o

#define __MINLOC__   MINLOC
#define __MAXLOC__   MAXLOC
#define __SORT__     SORT_QUICK_
#define __FFT__      FFT_RADIX_EIGHT_

MODULE SIGNAL

  USE, INTRINSIC :: ISO_C_BINDING,   ONLY: IK => C_INT, RK => C_DOUBLE, C_SIZEOF
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: IK
  PUBLIC :: RK
  INTEGER(IK), PUBLIC, PARAMETER :: IK_SIZE = C_SIZEOF(IK)
  INTEGER(IK), PUBLIC, PARAMETER :: RK_SIZE = C_SIZEOF(RK)

  ! GLOBAL PARAMETERS
  REAL(RK),    PUBLIC, PARAMETER :: ONE_PI                 = 2.0_RK*ACOS(0.0_RK) ! 1.0*PI
  REAL(RK),    PUBLIC, PARAMETER :: TWO_PI                 = 2.0_RK*ONE_PI       ! 2.0*PI
  INTEGER(IK), PUBLIC, PARAMETER :: FFT_FORWARD            = +1_IK               ! FORWARD FFT
  INTEGER(IK), PUBLIC, PARAMETER :: FFT_INVERSE            = -1_IK               ! INVERSE FFT
  INTEGER(IK), PUBLIC, PARAMETER :: PEAK_WIDTH             = 2_IK                ! PEAK WIDTH
  REAL(RK),    PUBLIC, PARAMETER :: PEAK_LEVEL             = -10.0_RK            ! PEAK THRESHOLD LEVEL
  INTEGER(IK), PUBLIC, PARAMETER :: DECOMPOSITION_MAX      = 0_IK                ! DECOMPOSITION WITH MAXIMUM BIN
  INTEGER(IK), PUBLIC, PARAMETER :: DECOMPOSITION_PEAK     = +1_IK               ! DECOMPOSITION WITH PEAKS
  INTEGER(IK), PUBLIC, PARAMETER :: DECOMPOSITION_PEAK_FFT = -1_IK               ! DECOMPOSITION WITH PEAKS, FFT ESTIMATION
  REAL(RK),    PUBLIC, PARAMETER :: SVD_LEVEL              = 1.0E-6_RK           ! SINGULAR VALUE THRESHOLD LEVEL

  ! GLOBAL PROCEDURES WITH C BINDING (T) AND WITHOUT (F)
  PUBLIC :: CONVERT_REAL_     ! (SUB) (T) CONVERT_REAL_(<LENGTH>,<R_PART>,<SEQUENCE>)
  PUBLIC :: CONVERT_COMPLEX_  ! (SUB) (T) CONVERT_COMPLEX_(<LENGTH>,<R_PART>,<I_PART>,<SEQUENCE>)
  PUBLIC :: CONVERT_          ! (SUB) (F) CONVERT_(<LENGTH>,<R_PART>[,<I_PART>],<SEQUENCE>)
  PUBLIC :: FFT_RADIX_TWO_    ! (SUB) (T) FFT_RADIX_TWO_(<LENGTH>,<DIRECTION>,<SEQUENCE>)
  PUBLIC :: FFT_RADIX_EIGHT_  ! (SUB) (T) FFT_RADIX_EIGHT_(<LENGTH>,<DIRECTION>,<SEQUENCE>)
  PUBLIC :: FFT_EXTERNAL_     ! (SUB) (T) FFT_EXTERNAL_(<LENGTH>,<DIRECTION>,<SEQUENCE>)
  PUBLIC :: FFRFT_            ! (SUB) (T) FFRFT_(<LENGTH>,<ARGUMENT>,<SEQUENCE>)
  PUBLIC :: WINDOW_           ! (SUB) (T) WINDOW_(<LENGTH>,<ORDER>,<SEQUENCE>)
  PUBLIC :: PEAK_             ! (FUN) (T) PEAK_(<LENGTH>,<SEQUENCE>,<PEAK_ID>)
  PUBLIC :: FREQUENCY_        ! (FUN) (T) FREQUENCY_(<FLAG>,<PEAK>,<LENGTH>,<TOTAL>,<WINDOW>,<SEQUENCE>)
  PUBLIC :: DECOMPOSITION_    ! (SUB) (T) DECOMPOSITION(<FLAG>,<METHOD>,<LENGTH>,<TOTAL>,<WINDOW>,<SEQUENCE>,<LOOP>,<FREQUENCY>,<COS_AMP>,<SIN_AMP>)
  PUBLIC :: FREQUENCY_LIST_   ! (SUB) (T) FREQUENCY_LIST_(<FLAG>,<METHOD>,<LENGTH>,<TOTAL>,<WINDOW>,<SEQUENCE>,<LOOP>,<FREQUENCY>)
  PUBLIC :: AMPLITUDE_LIST_   ! (SUB) (T) AMPLITUDE_LIST_(<FLAG>,<LENGTH>,<TOTAL>,<WINDOW>,<SEQUENCE>,<LOOP>,<FREQUENCY>,<COS_AMP>,<SIN_AMP>)
  PUBLIC :: AMPLITUDE_        ! (SUB) (T) AMPLITUDE(<LENGTH>,<TOTAL>,<WINDOW>,<SEQUENCE>,<FREQUENCY>,<COS_AMP>,<SIN_AMP>,<AMP>)
  PUBLIC :: SVD_              ! (SUB) (F) SVD_(<NR>,<NC>,<MATRIX>(<NR>,<NC>),<SVD_LIST>(MIN(<NR>,<NC>)),<U_MATRIX>(<NR>,<NR>),<V_MATRIX>(<NC>,<NC>))
  PUBLIC :: SVD_LIST_         ! (SUB) (F) SVD_LIST_(<NR>,<NC>,<MATRIX>(<NR>,<NC>),<SVD_LIST>(MIN(<NR>,<NC>)))
  PUBLIC :: LEAST_SQUARES_    ! (SUB) (F) LEAST_SQUARES_(<NR>,<NC>,<MATRIX>(<NR>,<NC>),<VECTOR>(<NR>),<SOLUTION>(<NC>))
  PUBLIC :: FIT_              ! (SUB) (T) FIT_(<LENGTH>,<SEQUENCE>,<LOOP>,<FREQUENCY>,<MEAN>,<COS_AMP>,<SIN_AMP>,<ERROR>)
  PUBLIC :: FILTER_           ! (SUB) (T) FILTER(<LENGTH>,<SEQUENCE>,<LIMIT>)

  ! EXTERNAL PROCEDURES
  EXTERNAL :: DGESVD
  EXTERNAL :: DGEMV
  EXTERNAL :: DGEMM

  ! INTERFACES
  INTERFACE CONVERT_
    MODULE PROCEDURE CONVERT_REAL_
    MODULE PROCEDURE CONVERT_COMPLEX_
  END INTERFACE

  ! MEMORIZATION
  TYPE TABLE
    INTEGER(IK), DIMENSION(:), ALLOCATABLE :: BIT_FFT
    INTEGER(IK), DIMENSION(:), ALLOCATABLE :: BIT_FFRFT
    REAL(RK), DIMENSION(:), ALLOCATABLE :: TRIG_FFT
    REAL(RK), DIMENSION(:), ALLOCATABLE :: TRIG_FFRFT
    REAL(RK), DIMENSION(:), ALLOCATABLE :: COS_FST
    REAL(RK), DIMENSION(:), ALLOCATABLE :: SIN_FST
    REAL(RK), DIMENSION(:), ALLOCATABLE :: COS_LST
    REAL(RK), DIMENSION(:), ALLOCATABLE :: SIN_LST
  END TYPE

  TYPE(TABLE), PROTECTED :: BANK

  PUBLIC :: COMPUTE_TABLE_
  PUBLIC :: DESTROY_TABLE_

  PUBLIC :: FREQUENCY__
  PUBLIC :: DECOMPOSITION__
  PUBLIC :: FREQUENCY_LIST__

  CONTAINS

  ! CONVERT INPUT SEQUENCE (REAL)
  ! (SUBROUTINE) CONVERT_REAL_(<LENGTH>, <R_PART>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <R_PART>               -- (IN)     INPUT SEQUENCE R-PART (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (OUT)    SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...] AND SI_I=0.0_RK FOR ALL I
  ! void    convert_real_(int*, double*, double*) ;
  SUBROUTINE CONVERT_REAL_(LENGTH, R_PART, SEQUENCE) &
    BIND(C, NAME = "convert_real_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: R_PART
    REAL(RK), INTENT(OUT), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    SEQUENCE = 0.0_RK
    SEQUENCE(1_IK:2_IK*LENGTH:2_IK) = R_PART
  END SUBROUTINE CONVERT_REAL_

  ! CONVERT INPUT SEQUENCE (COMPLEX)
  ! (SUBROUTINE) CONVERT_COMPLEX_(<LENGTH>, <R_PART>, <I_PART>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <R_PART>               -- (IN)     INPUT SEQUENCE R-PART (RK ARRAY OF LENGTH = <LENGTH>)
  ! <I_PART>               -- (IN)     INPUT SEQUENCE I-PART (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (OUT)    SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! void    convert_complex_(int*, double*, double*, double*) ;
  SUBROUTINE CONVERT_COMPLEX_(LENGTH, R_PART, I_PART, SEQUENCE) &
    BIND(C, NAME = "convert_complex_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: R_PART
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: I_PART
    REAL(RK), INTENT(OUT), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    SEQUENCE = 0.0_RK
    SEQUENCE(1_IK:2_IK*LENGTH:2_IK) = R_PART
    SEQUENCE(2_IK:2_IK*LENGTH:2_IK) = I_PART
  END SUBROUTINE CONVERT_COMPLEX_

  ! (NRF77) COMPLEX DISCRETE FOURIER TRANSFORM (POWER OF TWO INPUT LENGTH INPUT)
  ! (SUBROUTINE) FFT_RADIX_TWO_(<LENGTH>, <DIRECTION>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <DIRECTION>            -- (IN)     DIRECTION (IK), FFT_FORWARD = +1_IK OR FFT_INVERSE = -1_IK
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <SEQUENCE>             -- (OUT)    CDFT (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., FR_I, FI_I, ...]
  ! void    fft_radix_two_(int*, int*, double*) ;
  SUBROUTINE FFT_RADIX_TWO_(LENGTH, DIRECTION, SEQUENCE) &
    BIND(C, NAME = "fft_radix_two_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    INTEGER(IK), INTENT(IN) :: DIRECTION
    REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(INOUT) :: SEQUENCE
    INTEGER(IK) :: I, J
    INTEGER(IK) :: N, M
    INTEGER(IK) :: LIMIT, STEP
    REAL(RK) :: PIM, PRE
    REAL(RK) :: ANG
    REAL(RK) :: WR, WI, WPR, WPI, WCOPY
    REAL(RK) :: FACTOR
    N = 2_IK*LENGTH
    J = 1_IK
    FACTOR = REAL(DIRECTION, RK)*TWO_PI
    DO I = 1_IK, N, 2_IK
      IF (J > I) THEN
        PRE = SEQUENCE(J)
        PIM = SEQUENCE(J+1_IK)
        SEQUENCE(J) = SEQUENCE(I)
        SEQUENCE(J+1_IK) = SEQUENCE(I+1_IK)
        SEQUENCE(I) = PRE
        SEQUENCE(I+1_IK) = PIM
      END IF
      M = N/2_IK
      DO
        IF ((M > 2_IK) .AND. (J > M)) THEN
          J = J-M
          M = M/2_IK
        ELSE
          EXIT
        END IF
      END DO
      J = J+M
    END DO
    LIMIT = 2_IK
    DO
      IF (N > LIMIT) THEN
        STEP = 2_IK*LIMIT
        ANG = FACTOR/REAL(LIMIT, RK)
        WPI = SIN(ANG)
        ANG = SIN(0.5_RK*ANG)
        WPR = -2.0_RK*ANG*ANG
        WR = 1.0_RK
        WI = 0.0_RK
        DO M = 1_IK, LIMIT, 2_IK
          DO I = M, N, STEP
            J = I+LIMIT
            PRE = WR*SEQUENCE(J)-WI*SEQUENCE(J+1_IK)
            PIM = WR*SEQUENCE(J+1_IK)+WI*SEQUENCE(J)
            SEQUENCE(J) = SEQUENCE(I)-PRE
            SEQUENCE(J+1_IK) = SEQUENCE(I+1_IK)-PIM
            SEQUENCE(I) = SEQUENCE(I)+PRE
            SEQUENCE(I+1_IK) = SEQUENCE(I+1_IK)+PIM
          END DO
          WCOPY = WR
          WR = WR*WPR-WI*WPI+WR
          WI = WI*WPR+WCOPY*WPI+WI
        END DO
        LIMIT = STEP
      ELSE
        EXIT
      END IF
    END DO
  END SUBROUTINE FFT_RADIX_TWO_

  ! (TAKUYA OOURA) COMPLEX DISCRETE FOURIER TRANSFORM (POWER OF TWO INPUT LENGTH INPUT)
  ! (SUBROUTINE) FFT_RADIX_EIGHT_(<LENGTH>, <DIRECTION>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <DIRECTION>            -- (IN)     DIRECTION (IK), FFT_FORWARD = +1_IK OR FFT_INVERSE = -1_IK
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <SEQUENCE>             -- (OUT)    CDFT (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., FR_I, FI_I, ...]
  ! void    fft_radix_eight_(int*, int*, double*) ;
  SUBROUTINE FFT_RADIX_EIGHT_(LENGTH, DIRECTION, SEQUENCE) &
    BIND(C, NAME = "fft_radix_eight_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    INTEGER(IK), INTENT(IN) :: DIRECTION
    REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(INOUT) :: SEQUENCE
    INTEGER(IK), DIMENSION(2_IK+INT(SQRT(REAL(LENGTH/2_IK, RK)), IK)) :: IP
    REAL(RK), DIMENSION(LENGTH/2_IK) :: WORK
    IP = 0_IK
    WORK = 0.0_RK
    CALL CDFT_(2_IK*LENGTH, DIRECTION, SEQUENCE, IP, WORK)
  END SUBROUTINE FFT_RADIX_EIGHT_

  ! (TAKUYA OOURA) CDFT_
  SUBROUTINE CDFT_(LENGTH, DIRECTION, SEQUENCE, IP, WORK)
    INTEGER(IK), INTENT(IN) :: LENGTH
    INTEGER(IK), INTENT(IN) :: DIRECTION
    REAL(RK), INTENT(INOUT) :: SEQUENCE(0_IK : *)
    INTEGER(IK), INTENT(INOUT) :: IP(0_IK : *)
    REAL(RK), INTENT(INOUT) :: WORK(0_IK : *)
    CALL MAKE_FFT_TABLE_(LENGTH/4_IK, IP, WORK)
    IF (DIRECTION == FFT_FORWARD) THEN
      CALL BIT_REVERSE_(LENGTH, IP(2_IK), SEQUENCE)
      CALL CFT_FORWARD_(LENGTH, SEQUENCE, WORK)
    ELSE IF(DIRECTION == FFT_INVERSE) THEN
      CALL BIT_REVERSE_CONJUGATE_(LENGTH, IP(2_IK), SEQUENCE)
      CALL CFT_INVERSE_(LENGTH, SEQUENCE, WORK)
    END IF
  END SUBROUTINE CDFT_

  SUBROUTINE BIT_REVERSE_(N, IP, A)
    INTEGER(IK), INTENT(IN) :: N
    INTEGER(IK), INTENT(INOUT) :: IP(0_IK : *)
    REAL(RK), INTENT(INOUT) :: A(0_IK : N - 1_IK)
    INTEGER(IK) :: J, J1, K, K1, L, M, M2
    REAL(RK) :: XR, XI, YR, YI
    IP(0_IK) = 0_IK
    L = N
    M = 1_IK
    DO WHILE (8_IK * M < L)
      L = L / 2_IK
      DO J = 0_IK, M - 1_IK
        IP(M + J) = IP(J) + L
      END DO
      M = M * 2_IK
    END DO
    M2 = 2_IK * M
    IF (8_IK * M == L) THEN
      DO K = 0_IK, M - 1_IK
        DO J = 0_IK, K - 1_IK
          J1 = 2_IK * J + IP(K)
          K1 = 2_IK * K + IP(J)
          XR = A(J1)
          XI = A(J1 + 1_IK)
          YR = A(K1)
          YI = A(K1 + 1_IK)
          A(J1) = YR
          A(J1 + 1_IK) = YI
          A(K1) = XR
          A(K1 + 1_IK) = XI
          J1 = J1 + M2
          K1 = K1 + 2_IK * M2
          XR = A(J1)
          XI = A(J1 + 1_IK)
          YR = A(K1)
          YI = A(K1 + 1_IK)
          A(J1) = YR
          A(J1 + 1_IK) = YI
          A(K1) = XR
          A(K1 + 1_IK) = XI
          J1 = J1 + M2
          K1 = K1 - M2
          XR = A(J1)
          XI = A(J1 + 1_IK)
          YR = A(K1)
          YI = A(K1 + 1_IK)
          A(J1) = YR
          A(J1 + 1_IK) = YI
          A(K1) = XR
          A(K1 + 1_IK) = XI
          J1 = J1 + M2
          K1 = K1 + 2_IK * M2
          XR = A(J1)
          XI = A(J1 + 1_IK)
          YR = A(K1)
          YI = A(K1 + 1_IK)
          A(J1) = YR
          A(J1 + 1_IK) = YI
          A(K1) = XR
          A(K1 + 1_IK) = XI
        END DO
        J1 = 2_IK * K + M2 + IP(K)
        K1 = J1 + M2
        XR = A(J1)
        XI = A(J1 + 1_IK)
        YR = A(K1)
        YI = A(K1 + 1_IK)
        A(J1) = YR
        A(J1 + 1_IK) = YI
        A(K1) = XR
        A(K1 + 1_IK) = XI
      END DO
    ELSE
      DO K = 1_IK, M - 1_IK
        DO J = 0_IK, K - 1_IK
          J1 = 2_IK * J + IP(K)
          K1 = 2_IK * K + IP(J)
          XR = A(J1)
          XI = A(J1 + 1_IK)
          YR = A(K1)
          YI = A(K1 + 1_IK)
          A(J1) = YR
          A(J1 + 1_IK) = YI
          A(K1) = XR
          A(K1 + 1_IK) = XI
          J1 = J1 + M2
          K1 = K1 + M2
          XR = A(J1)
          XI = A(J1 + 1_IK)
          YR = A(K1)
          YI = A(K1 + 1_IK)
          A(J1) = YR
          A(J1 + 1_IK) = YI
          A(K1) = XR
          A(K1 + 1_IK) = XI
        END DO
      END DO
    END IF
  END SUBROUTINE BIT_REVERSE_

  SUBROUTINE BIT_REVERSE_CONJUGATE_(N, IP, A)
    INTEGER(IK), INTENT(IN) :: N
    INTEGER(IK), INTENT(INOUT) :: IP(0_IK : *)
    REAL(RK), INTENT(INOUT) :: A(0_IK : N - 1_IK)
    INTEGER(IK) :: J, J1, K, K1, L, M, M2
    REAL(RK) :: XR, XI, YR, YI
    IP(0_IK) = 0_IK
    L = N
    M = 1_IK
    DO WHILE (8_IK * M < L)
      L = L / 2_IK
      DO J = 0_IK, M - 1_IK
        IP(M + J) = IP(J) + L
      END DO
      M = M * 2_IK
    END DO
    M2 = 2_IK * M
    IF (8_IK * M == L) THEN
      DO K = 0_IK, M - 1_IK
        DO J = 0_IK, K - 1_IK
          J1 = 2_IK * J + IP(K)
          K1 = 2_IK * K + IP(J)
          XR = A(J1)
          XI = -A(J1 + 1_IK)
          YR = A(K1)
          YI = -A(K1 + 1_IK)
          A(J1) = YR
          A(J1 + 1_IK) = YI
          A(K1) = XR
          A(K1 + 1_IK) = XI
          J1 = J1 + M2
          K1 = K1 + 2_IK * M2
          XR = A(J1)
          XI = -A(J1 + 1_IK)
          YR = A(K1)
          YI = -A(K1 + 1_IK)
          A(J1) = YR
          A(J1 + 1_IK) = YI
          A(K1) = XR
          A(K1 + 1_IK) = XI
          J1 = J1 + M2
          K1 = K1 - M2
          XR = A(J1)
          XI = -A(J1 + 1_IK)
          YR = A(K1)
          YI = -A(K1 + 1_IK)
          A(J1) = YR
          A(J1 + 1_IK) = YI
          A(K1) = XR
          A(K1 + 1_IK) = XI
          J1 = J1 + M2
          K1 = K1 + 2_IK * M2
          XR = A(J1)
          XI = -A(J1 + 1_IK)
          YR = A(K1)
          YI = -A(K1 + 1_IK)
          A(J1) = YR
          A(J1 + 1_IK) = YI
          A(K1) = XR
          A(K1 + 1_IK) = XI
        END DO
        K1 = 2_IK * K + IP(K)
        A(K1 + 1_IK) = -A(K1 + 1_IK)
        J1 = K1 + M2
        K1 = J1 + M2
        XR = A(J1)
        XI = -A(J1 + 1_IK)
        YR = A(K1)
        YI = -A(K1 + 1_IK)
        A(J1) = YR
        A(J1 + 1_IK) = YI
        A(K1) = XR
        A(K1 + 1_IK) = XI
        K1 = K1 + M2
        A(K1 + 1_IK) = -A(K1 + 1_IK)
      END DO
    ELSE
      A(1_IK) = -A(1_IK)
      A(M2 + 1_IK) = -A(M2 + 1_IK)
      DO K = 1_IK, M - 1_IK
        DO J = 0_IK, K - 1_IK
          J1 = 2_IK * J + IP(K)
          K1 = 2_IK * K + IP(J)
          XR = A(J1)
          XI = -A(J1 + 1_IK)
          YR = A(K1)
          YI = -A(K1 + 1_IK)
          A(J1) = YR
          A(J1 + 1_IK) = YI
          A(K1) = XR
          A(K1 + 1_IK) = XI
          J1 = J1 + M2
          K1 = K1 + M2
          XR = A(J1)
          XI = -A(J1 + 1_IK)
          YR = A(K1)
          YI = -A(K1 + 1_IK)
          A(J1) = YR
          A(J1 + 1_IK) = YI
          A(K1) = XR
          A(K1 + 1_IK) = XI
        END DO
        K1 = 2_IK * K + IP(K)
        A(K1 + 1_IK) = -A(K1 + 1_IK)
        A(K1 + M2 + 1_IK) = -A(K1 + M2 + 1_IK)
      END DO
    END IF
  END SUBROUTINE BIT_REVERSE_CONJUGATE_

  SUBROUTINE MAKE_FFT_TABLE_(NW, IP, W)
    INTEGER(IK), INTENT(IN) :: NW
    INTEGER(IK), INTENT(INOUT) :: IP(0_IK : *)
    REAL(RK), INTENT(INOUT) :: W(0_IK : NW - 1_IK)
    INTEGER(IK) :: J, NWH
    REAL(RK) :: DELTA, X, Y
    NWH = NW / 2_IK
    DELTA = ONE_PI / (4.0_RK * REAL(NWH, RK))
    W(0_IK) = 1_IK
    W(1_IK) = 0_IK
    W(NWH) = COS(DELTA * REAL(NWH, RK))
    W(NWH + 1_IK) = W(NWH)
    IF (NWH > 2_IK) THEN
      DO J = 2_IK, NWH - 2_IK, 2_IK
        X = COS(DELTA * REAL(J,RK))
        Y = SIN(DELTA * REAL(J,RK))
        W(J) = X
        W(J + 1) = Y
        W(NW - J) = Y
        W(NW - J + 1_IK) = X
      END DO
      DO J = NWH - 2_IK, 2_IK, -2_IK
        X = W(2_IK * J)
        Y = W(2_IK * J + 1_IK)
        W(NWH + J) = X
        W(NWH + J + 1_IK) = Y
      END DO
      IP(0_IK) = NW
      IP(1_IK) = 1_IK
      CALL BIT_REVERSE_(NW, IP(2_IK), W)
    END IF
  END SUBROUTINE MAKE_FFT_TABLE_

  SUBROUTINE CFT_FORWARD_(N, A, W)
    INTEGER(IK), INTENT(IN) :: N
    REAL(RK), INTENT(INOUT) :: A(0_IK : N - 1_IK)
    REAL(RK), INTENT(IN) :: W(0_IK : *)
    INTEGER(IK) :: J, J1, J2, J3, L
    REAL(RK) :: X0R, X0I, X1R, X1I, X2R, X2I, X3R, X3I
    L = 2_IK
    IF (N >= 16_IK) THEN
      CALL CFT_1ST_(N, A, W)
      L = 16_IK
      DO WHILE (8_IK * L <= N)
        CALL CFT_MDL_(N, L, A, W)
        L = 8_IK * L
      END DO
    END IF
    IF (2_IK * L < N) THEN
      DO J = 0_IK, L - 2_IK, 2_IK
        J1 = J + L
        J2 = J1 + L
        J3 = J2 + L
        X0R = A(J) + A(J1)
        X0I = A(J + 1_IK) + A(J1 + 1_IK)
        X1R = A(J) - A(J1)
        X1I = A(J + 1_IK) - A(J1 + 1_IK)
        X2R = A(J2) + A(J3)
        X2I = A(J2 + 1_IK) + A(J3 + 1_IK)
        X3R = A(J2) - A(J3)
        X3I = A(J2 + 1_IK) - A(J3 + 1_IK)
        A(J) = X0R + X2R
        A(J + 1_IK) = X0I + X2I
        A(J2) = X0R - X2R
        A(J2 + 1_IK) = X0I - X2I
        A(J1) = X1R - X3I
        A(J1 + 1_IK) = X1I + X3R
        A(J3) = X1R + X3I
        A(J3 + 1_IK) = X1I - X3R
      END DO
    ELSE IF (2_IK * L == N) THEN
      DO J = 0_IK, L - 2_IK, 2_IK
        J1 = J + L
        X0R = A(J) - A(J1)
        X0I = A(J + 1_IK) - A(J1 + 1_IK)
        A(J) = A(J) + A(J1)
        A(J + 1_IK) = A(J + 1_IK) + A(J1 + 1_IK)
        A(J1) = X0R
        A(J1 + 1_IK) = X0I
      END DO
    END IF
  END SUBROUTINE CFT_FORWARD_

  SUBROUTINE CFT_INVERSE_(N, A, W)
    INTEGER(IK), INTENT(IN) :: N
    REAL(RK), INTENT(INOUT) :: A(0_IK : N - 1_IK)
    REAL(RK), INTENT(IN) :: W(0_IK : *)
    INTEGER(IK) :: J, J1, J2, J3, J4, J5, J6, J7, L
    REAL(RK) ::  WN4R, X0R, X0I, X1R, X1I, X2R, X2I, X3R, X3I
    REAL(RK) :: Y0R, Y0I, Y1R, Y1I, Y2R, Y2I, Y3R, Y3I
    REAL(RK) :: Y4R, Y4I, Y5R, Y5I, Y6R, Y6I, Y7R, Y7I
    L = 2_IK
    IF (N > 16_IK) THEN
      CALL CFT_1ST_(N, A, W)
      L = 16_IK
      DO WHILE (8_IK * L < N)
        CALL CFT_MDL_(N, L, A, W)
        L = 8_IK * L
      END DO
    END IF
    IF (4_IK * L < N) THEN
      WN4R = W(2_IK)
      DO J = 0_IK, L - 2_IK, 2_IK
        J1 = J + L
        J2 = J1 + L
        J3 = J2 + L
        J4 = J3 + L
        J5 = J4 + L
        J6 = J5 + L
        J7 = J6 + L
        X0R = A(J) + A(J1)
        X0I = -A(J + 1_IK) - A(J1 + 1_IK)
        X1R = A(J) - A(J1)
        X1I = -A(J + 1_IK) + A(J1 + 1_IK)
        X2R = A(J2) + A(J3)
        X2I = A(J2 + 1_IK) + A(J3 + 1_IK)
        X3R = A(J2) - A(J3)
        X3I = A(J2 + 1_IK) - A(J3 + 1_IK)
        Y0R = X0R + X2R
        Y0I = X0I - X2I
        Y2R = X0R - X2R
        Y2I = X0I + X2I
        Y1R = X1R - X3I
        Y1I = X1I - X3R
        Y3R = X1R + X3I
        Y3I = X1I + X3R
        X0R = A(J4) + A(J5)
        X0I = A(J4 + 1_IK) + A(J5 + 1_IK)
        X1R = A(J4) - A(J5)
        X1I = A(J4 + 1_IK) - A(J5 + 1_IK)
        X2R = A(J6) + A(J7)
        X2I = A(J6 + 1_IK) + A(J7 + 1_IK)
        X3R = A(J6) - A(J7)
        X3I = A(J6 + 1_IK) - A(J7 + 1_IK)
        Y4R = X0R + X2R
        Y4I = X0I + X2I
        Y6R = X0R - X2R
        Y6I = X0I - X2I
        X0R = X1R - X3I
        X0I = X1I + X3R
        X2R = X1R + X3I
        X2I = X1I - X3R
        Y5R = WN4R * (X0R - X0I)
        Y5I = WN4R * (X0R + X0I)
        Y7R = WN4R * (X2R - X2I)
        Y7I = WN4R * (X2R + X2I)
        A(J1) = Y1R + Y5R
        A(J1 + 1_IK) = Y1I - Y5I
        A(J5) = Y1R - Y5R
        A(J5 + 1_IK) = Y1I + Y5I
        A(J3) = Y3R - Y7I
        A(J3 + 1_IK) = Y3I - Y7R
        A(J7) = Y3R + Y7I
        A(J7 + 1_IK) = Y3I + Y7R
        A(J) = Y0R + Y4R
        A(J + 1_IK) = Y0I - Y4I
        A(J4) = Y0R - Y4R
        A(J4 + 1_IK) = Y0I + Y4I
        A(J2) = Y2R - Y6I
        A(J2 + 1_IK) = Y2I - Y6R
        A(J6) = Y2R + Y6I
        A(J6 + 1_IK) = Y2I + Y6R
      END DO
    ELSE IF (4_IK * L == N) THEN
      DO J = 0_IK, L - 2_IK, 2_IK
        J1 = J + L
        J2 = J1 + L
        J3 = J2 + L
        X0R = A(J) + A(J1)
        X0I = -A(J + 1_IK) - A(J1 + 1_IK)
        X1R = A(J) - A(J1)
        X1I = -A(J + 1_IK) + A(J1 + 1_IK)
        X2R = A(J2) + A(J3)
        X2I = A(J2 + 1_IK) + A(J3 + 1_IK)
        X3R = A(J2) - A(J3)
        X3I = A(J2 + 1_IK) - A(J3 + 1_IK)
        A(J) = X0R + X2R
        A(J + 1_IK) = X0I - X2I
        A(J2) = X0R - X2R
        A(J2 + 1_IK) = X0I + X2I
        A(J1) = X1R - X3I
        A(J1 + 1_IK) = X1I - X3R
        A(J3) = X1R + X3I
        A(J3 + 1_IK) = X1I + X3R
      END DO
    ELSE
      DO J = 0_IK, L - 2_IK, 2_IK
        J1 = J + L
        X0R = A(J) - A(J1)
        X0I = -A(J + 1_IK) + A(J1 + 1_IK)
        A(J) = A(J) + A(J1)
        A(J + 1_IK) = -A(J + 1_IK) - A(J1 + 1_IK)
        A(J1) = X0R
        A(J1 + 1_IK) = X0I
      END DO
    END IF
  END SUBROUTINE CFT_INVERSE_

  SUBROUTINE CFT_1ST_(N, A, W)
    INTEGER(IK), INTENT(IN) :: N
    REAL(RK), INTENT(INOUT) :: A(0_IK : N - 1_IK)
    REAL(RK), INTENT(IN) :: W(0_IK : *)
    INTEGER(IK) :: J, K1
    REAL(RK) ::  WN4R, WTMP, WK1R, WK1I, WK2R, WK2I, WK3R, WK3I
    REAL(RK) ::  WK4R, WK4I, WK5R, WK5I, WK6R, WK6I, WK7R, WK7I
    REAL(RK) ::  X0R, X0I, X1R, X1I, X2R, X2I, X3R, X3I
    REAL(RK) ::  Y0R, Y0I, Y1R, Y1I, Y2R, Y2I, Y3R, Y3I
    REAL(RK) :: Y4R, Y4I, Y5R, Y5I, Y6R, Y6I, Y7R, Y7I
    WN4R = W(2_IK)
    X0R = A(0_IK) + A(2_IK)
    X0I = A(1_IK) + A(3_IK)
    X1R = A(0_IK) - A(2_IK)
    X1I = A(1_IK) - A(3_IK)
    X2R = A(4_IK) + A(6_IK)
    X2I = A(5_IK) + A(7_IK)
    X3R = A(4_IK) - A(6_IK)
    X3I = A(5_IK) - A(7_IK)
    Y0R = X0R + X2R
    Y0I = X0I + X2I
    Y2R = X0R - X2R
    Y2I = X0I - X2I
    Y1R = X1R - X3I
    Y1I = X1I + X3R
    Y3R = X1R + X3I
    Y3I = X1I - X3R
    X0R = A(8_IK) + A(10_IK)
    X0I = A(9_IK) + A(11_IK)
    X1R = A(8_IK) - A(10_IK)
    X1I = A(9_IK) - A(11_IK)
    X2R = A(12_IK) + A(14_IK)
    X2I = A(13_IK) + A(15_IK)
    X3R = A(12_IK) - A(14_IK)
    X3I = A(13_IK) - A(15_IK)
    Y4R = X0R + X2R
    Y4I = X0I + X2I
    Y6R = X0R - X2R
    Y6I = X0I - X2I
    X0R = X1R - X3I
    X0I = X1I + X3R
    X2R = X1R + X3I
    X2I = X1I - X3R
    Y5R = WN4R * (X0R - X0I)
    Y5I = WN4R * (X0R + X0I)
    Y7R = WN4R * (X2R - X2I)
    Y7I = WN4R * (X2R + X2I)
    A(2_IK) = Y1R + Y5R
    A(3_IK) = Y1I + Y5I
    A(10_IK) = Y1R - Y5R
    A(11_IK) = Y1I - Y5I
    A(6_IK) = Y3R - Y7I
    A(7_IK) = Y3I + Y7R
    A(14_IK) = Y3R + Y7I
    A(15_IK) = Y3I - Y7R
    A(0_IK) = Y0R + Y4R
    A(1_IK) = Y0I + Y4I
    A(8_IK) = Y0R - Y4R
    A(9_IK) = Y0I - Y4I
    A(4_IK) = Y2R - Y6I
    A(5_IK) = Y2I + Y6R
    A(12_IK) = Y2R + Y6I
    A(13_IK) = Y2I - Y6R
    IF (N > 16_IK) THEN
      WK1R = W(4_IK)
      WK1I = W(5_IK)
      X0R = A(16_IK) + A(18_IK)
      X0I = A(17_IK) + A(19_IK)
      X1R = A(16_IK) - A(18_IK)
      X1I = A(17_IK) - A(19_IK)
      X2R = A(20_IK) + A(22_IK)
      X2I = A(21_IK) + A(23_IK)
      X3R = A(20_IK) - A(22_IK)
      X3I = A(21_IK) - A(23_IK)
      Y0R = X0R + X2R
      Y0I = X0I + X2I
      Y2R = X0R - X2R
      Y2I = X0I - X2I
      Y1R = X1R - X3I
      Y1I = X1I + X3R
      Y3R = X1R + X3I
      Y3I = X1I - X3R
      X0R = A(24_IK) + A(26_IK)
      X0I = A(25_IK) + A(27_IK)
      X1R = A(24_IK) - A(26_IK)
      X1I = A(25_IK) - A(27_IK)
      X2R = A(28_IK) + A(30_IK)
      X2I = A(29_IK) + A(31_IK)
      X3R = A(28_IK) - A(30_IK)
      X3I = A(29_IK) - A(31_IK)
      Y4R = X0R + X2R
      Y4I = X0I + X2I
      Y6R = X0R - X2R
      Y6I = X0I - X2I
      X0R = X1R - X3I
      X0I = X1I + X3R
      X2R = X1R + X3I
      X2I = X3R - X1I
      Y5R = WK1I * X0R - WK1R * X0I
      Y5I = WK1I * X0I + WK1R * X0R
      Y7R = WK1R * X2R + WK1I * X2I
      Y7I = WK1R * X2I - WK1I * X2R
      X0R = WK1R * Y1R - WK1I * Y1I
      X0I = WK1R * Y1I + WK1I * Y1R
      A(18_IK) = X0R + Y5R
      A(19_IK) = X0I + Y5I
      A(26_IK) = Y5I - X0I
      A(27_IK) = X0R - Y5R
      X0R = WK1I * Y3R - WK1R * Y3I
      X0I = WK1I * Y3I + WK1R * Y3R
      A(22_IK) = X0R - Y7R
      A(23_IK) = X0I + Y7I
      A(30_IK) = Y7I - X0I
      A(31_IK) = X0R + Y7R
      A(16_IK) = Y0R + Y4R
      A(17_IK) = Y0I + Y4I
      A(24_IK) = Y4I - Y0I
      A(25_IK) = Y0R - Y4R
      X0R = Y2R - Y6I
      X0I = Y2I + Y6R
      A(20_IK) = WN4R * (X0R - X0I)
      A(21_IK) = WN4R * (X0I + X0R)
      X0R = Y6R - Y2I
      X0I = Y2R + Y6I
      A(28_IK) = WN4R * (X0R - X0I)
      A(29_IK) = WN4R * (X0I + X0R)
      K1 = 4_IK
      DO J = 32_IK, N - 16_IK, 16_IK
        K1 = K1 + 4_IK
        WK1R = W(K1)
        WK1I = W(K1 + 1_IK)
        WK2R = W(K1 + 2_IK)
        WK2I = W(K1 + 3_IK)
        WTMP = 2_IK * WK2I
        WK3R = WK1R - WTMP * WK1I
        WK3I = WTMP * WK1R - WK1I
        WK4R = 1_IK - WTMP * WK2I
        WK4I = WTMP * WK2R
        WTMP = 2_IK * WK4I
        WK5R = WK3R - WTMP * WK1I
        WK5I = WTMP * WK1R - WK3I
        WK6R = WK2R - WTMP * WK2I
        WK6I = WTMP * WK2R - WK2I
        WK7R = WK1R - WTMP * WK3I
        WK7I = WTMP * WK3R - WK1I
        X0R = A(J) + A(J + 2_IK)
        X0I = A(J + 1_IK) + A(J + 3_IK)
        X1R = A(J) - A(J + 2_IK)
        X1I = A(J + 1_IK) - A(J + 3_IK)
        X2R = A(J + 4_IK) + A(J + 6_IK)
        X2I = A(J + 5_IK) + A(J + 7_IK)
        X3R = A(J + 4_IK) - A(J + 6_IK)
        X3I = A(J + 5_IK) - A(J + 7_IK)
        Y0R = X0R + X2R
        Y0I = X0I + X2I
        Y2R = X0R - X2R
        Y2I = X0I - X2I
        Y1R = X1R - X3I
        Y1I = X1I + X3R
        Y3R = X1R + X3I
        Y3I = X1I - X3R
        X0R = A(J + 8_IK) + A(J + 10_IK)
        X0I = A(J + 9_IK) + A(J + 11_IK)
        X1R = A(J + 8_IK) - A(J + 10_IK)
        X1I = A(J + 9_IK) - A(J + 11_IK)
        X2R = A(J + 12_IK) + A(J + 14_IK)
        X2I = A(J + 13_IK) + A(J + 15_IK)
        X3R = A(J + 12_IK) - A(J + 14_IK)
        X3I = A(J + 13_IK) - A(J + 15_IK)
        Y4R = X0R + X2R
        Y4I = X0I + X2I
        Y6R = X0R - X2R
        Y6I = X0I - X2I
        X0R = X1R - X3I
        X0I = X1I + X3R
        X2R = X1R + X3I
        X2I = X1I - X3R
        Y5R = WN4R * (X0R - X0I)
        Y5I = WN4R * (X0R + X0I)
        Y7R = WN4R * (X2R - X2I)
        Y7I = WN4R * (X2R + X2I)
        X0R = Y1R + Y5R
        X0I = Y1I + Y5I
        A(J + 2_IK) = WK1R * X0R - WK1I * X0I
        A(J + 3_IK) = WK1R * X0I + WK1I * X0R
        X0R = Y1R - Y5R
        X0I = Y1I - Y5I
        A(J + 10_IK) = WK5R * X0R - WK5I * X0I
        A(J + 11_IK) = WK5R * X0I + WK5I * X0R
        X0R = Y3R - Y7I
        X0I = Y3I + Y7R
        A(J + 6_IK) = WK3R * X0R - WK3I * X0I
        A(J + 7_IK) = WK3R * X0I + WK3I * X0R
        X0R = Y3R + Y7I
        X0I = Y3I - Y7R
        A(J + 14_IK) = WK7R * X0R - WK7I * X0I
        A(J + 15_IK) = WK7R * X0I + WK7I * X0R
        A(J) = Y0R + Y4R
        A(J + 1_IK) = Y0I + Y4I
        X0R = Y0R - Y4R
        X0I = Y0I - Y4I
        A(J + 8_IK) = WK4R * X0R - WK4I * X0I
        A(J + 9_IK) = WK4R * X0I + WK4I * X0R
        X0R = Y2R - Y6I
        X0I = Y2I + Y6R
        A(J + 4_IK) = WK2R * X0R - WK2I * X0I
        A(J + 5_IK) = WK2R * X0I + WK2I * X0R
        X0R = Y2R + Y6I
        X0I = Y2I - Y6R
        A(J + 12_IK) = WK6R * X0R - WK6I * X0I
        A(J + 13_IK) = WK6R * X0I + WK6I * X0R
      END DO
    END IF
  END SUBROUTINE CFT_1ST_

  SUBROUTINE CFT_MDL_(N, L, A, W)
    INTEGER(IK), INTENT(IN) :: N
    INTEGER(IK), INTENT(IN) :: L
    REAL(RK), INTENT(INOUT) :: A(0_IK : N - 1_IK)
    REAL(RK), INTENT(IN) :: W(0_IK : *)
    INTEGER(IK) :: J, J1, J2, J3, J4, J5, J6, J7, K, K1, M
    REAL(RK) :: WN4R, WTMP, WK1R, WK1I, WK2R, WK2I, WK3R, WK3I
    REAL(RK) :: WK4R, WK4I, WK5R, WK5I, WK6R, WK6I, WK7R, WK7I
    REAL(RK) :: X0R, X0I, X1R, X1I, X2R, X2I, X3R, X3I
    REAL(RK) :: Y0R, Y0I, Y1R, Y1I, Y2R, Y2I, Y3R, Y3I
    REAL(RK) :: Y4R, Y4I, Y5R, Y5I, Y6R, Y6I, Y7R, Y7I
    M = 8_IK * L
    WN4R = W(2_IK)
    DO J = 0_IK, L - 2_IK, 2_IK
      J1 = J + L
      J2 = J1 + L
      J3 = J2 + L
      J4 = J3 + L
      J5 = J4 + L
      J6 = J5 + L
      J7 = J6 + L
      X0R = A(J) + A(J1)
      X0I = A(J + 1_IK) + A(J1 + 1_IK)
      X1R = A(J) - A(J1)
      X1I = A(J + 1_IK) - A(J1 + 1_IK)
      X2R = A(J2) + A(J3)
      X2I = A(J2 + 1_IK) + A(J3 + 1_IK)
      X3R = A(J2) - A(J3)
      X3I = A(J2 + 1_IK) - A(J3 + 1_IK)
      Y0R = X0R + X2R
      Y0I = X0I + X2I
      Y2R = X0R - X2R
      Y2I = X0I - X2I
      Y1R = X1R - X3I
      Y1I = X1I + X3R
      Y3R = X1R + X3I
      Y3I = X1I - X3R
      X0R = A(J4) + A(J5)
      X0I = A(J4 + 1_IK) + A(J5 + 1_IK)
      X1R = A(J4) - A(J5)
      X1I = A(J4 + 1_IK) - A(J5 + 1_IK)
      X2R = A(J6) + A(J7)
      X2I = A(J6 + 1_IK) + A(J7 + 1_IK)
      X3R = A(J6) - A(J7)
      X3I = A(J6 + 1_IK) - A(J7 + 1_IK)
      Y4R = X0R + X2R
      Y4I = X0I + X2I
      Y6R = X0R - X2R
      Y6I = X0I - X2I
      X0R = X1R - X3I
      X0I = X1I + X3R
      X2R = X1R + X3I
      X2I = X1I - X3R
      Y5R = WN4R * (X0R - X0I)
      Y5I = WN4R * (X0R + X0I)
      Y7R = WN4R * (X2R - X2I)
      Y7I = WN4R * (X2R + X2I)
      A(J1) = Y1R + Y5R
      A(J1 + 1_IK) = Y1I + Y5I
      A(J5) = Y1R - Y5R
      A(J5 + 1_IK) = Y1I - Y5I
      A(J3) = Y3R - Y7I
      A(J3 + 1_IK) = Y3I + Y7R
      A(J7) = Y3R + Y7I
      A(J7 + 1_IK) = Y3I - Y7R
      A(J) = Y0R + Y4R
      A(J + 1_IK) = Y0I + Y4I
      A(J4) = Y0R - Y4R
      A(J4 + 1_IK) = Y0I - Y4I
      A(J2) = Y2R - Y6I
      A(J2 + 1_IK) = Y2I + Y6R
      A(J6) = Y2R + Y6I
      A(J6 + 1_IK) = Y2I - Y6R
    END DO
    IF (M < N) THEN
      WK1R = W(4_IK)
      WK1I = W(5_IK)
      DO J = M, L + M - 2_IK, 2_IK
        J1 = J + L
        J2 = J1 + L
        J3 = J2 + L
        J4 = J3 + L
        J5 = J4 + L
        J6 = J5 + L
        J7 = J6 + L
        X0R = A(J) + A(J1)
        X0I = A(J + 1_IK) + A(J1 + 1_IK)
        X1R = A(J) - A(J1)
        X1I = A(J + 1_IK) - A(J1 + 1_IK)
        X2R = A(J2) + A(J3)
        X2I = A(J2 + 1_IK) + A(J3 + 1_IK)
        X3R = A(J2) - A(J3)
        X3I = A(J2 + 1_IK) - A(J3 + 1_IK)
        Y0R = X0R + X2R
        Y0I = X0I + X2I
        Y2R = X0R - X2R
        Y2I = X0I - X2I
        Y1R = X1R - X3I
        Y1I = X1I + X3R
        Y3R = X1R + X3I
        Y3I = X1I - X3R
        X0R = A(J4) + A(J5)
        X0I = A(J4 + 1_IK) + A(J5 + 1_IK)
        X1R = A(J4) - A(J5)
        X1I = A(J4 + 1_IK) - A(J5 + 1_IK)
        X2R = A(J6) + A(J7)
        X2I = A(J6 + 1_IK) + A(J7 + 1_IK)
        X3R = A(J6) - A(J7)
        X3I = A(J6 + 1_IK) - A(J7 + 1_IK)
        Y4R = X0R + X2R
        Y4I = X0I + X2I
        Y6R = X0R - X2R
        Y6I = X0I - X2I
        X0R = X1R - X3I
        X0I = X1I + X3R
        X2R = X1R + X3I
        X2I = X3R - X1I
        Y5R = WK1I * X0R - WK1R * X0I
        Y5I = WK1I * X0I + WK1R * X0R
        Y7R = WK1R * X2R + WK1I * X2I
        Y7I = WK1R * X2I - WK1I * X2R
        X0R = WK1R * Y1R - WK1I * Y1I
        X0I = WK1R * Y1I + WK1I * Y1R
        A(J1) = X0R + Y5R
        A(J1 + 1_IK) = X0I + Y5I
        A(J5) = Y5I - X0I
        A(J5 + 1_IK) = X0R - Y5R
        X0R = WK1I * Y3R - WK1R * Y3I
        X0I = WK1I * Y3I + WK1R * Y3R
        A(J3) = X0R - Y7R
        A(J3 + 1_IK) = X0I + Y7I
        A(J7) = Y7I - X0I
        A(J7 + 1_IK) = X0R + Y7R
        A(J) = Y0R + Y4R
        A(J + 1_IK) = Y0I + Y4I
        A(J4) = Y4I - Y0I
        A(J4 + 1_IK) = Y0R - Y4R
        X0R = Y2R - Y6I
        X0I = Y2I + Y6R
        A(J2) = WN4R * (X0R - X0I)
        A(J2 + 1_IK) = WN4R * (X0I + X0R)
        X0R = Y6R - Y2I
        X0I = Y2R + Y6I
        A(J6) = WN4R * (X0R - X0I)
        A(J6 + 1_IK) = WN4R * (X0I + X0R)
      END DO
      K1 = 4_IK
      DO K = 2_IK * M, N - M, M
        K1 = K1 + 4_IK
        WK1R = W(K1)
        WK1I = W(K1 + 1_IK)
        WK2R = W(K1 + 2_IK)
        WK2I = W(K1 + 3_IK)
        WTMP = 2_IK * WK2I
        WK3R = WK1R - WTMP * WK1I
        WK3I = WTMP * WK1R - WK1I
        WK4R = 1_IK - WTMP * WK2I
        WK4I = WTMP * WK2R
        WTMP = 2_IK * WK4I
        WK5R = WK3R - WTMP * WK1I
        WK5I = WTMP * WK1R - WK3I
        WK6R = WK2R - WTMP * WK2I
        WK6I = WTMP * WK2R - WK2I
        WK7R = WK1R - WTMP * WK3I
        WK7I = WTMP * WK3R - WK1I
        DO J = K, L + K - 2_IK, 2_IK
          J1 = J + L
          J2 = J1 + L
          J3 = J2 + L
          J4 = J3 + L
          J5 = J4 + L
          J6 = J5 + L
          J7 = J6 + L
          X0R = A(J) + A(J1)
          X0I = A(J + 1_IK) + A(J1 + 1_IK)
          X1R = A(J) - A(J1)
          X1I = A(J + 1_IK) - A(J1 + 1_IK)
          X2R = A(J2) + A(J3)
          X2I = A(J2 + 1_IK) + A(J3 + 1_IK)
          X3R = A(J2) - A(J3)
          X3I = A(J2 + 1_IK) - A(J3 + 1_IK)
          Y0R = X0R + X2R
          Y0I = X0I + X2I
          Y2R = X0R - X2R
          Y2I = X0I - X2I
          Y1R = X1R - X3I
          Y1I = X1I + X3R
          Y3R = X1R + X3I
          Y3I = X1I - X3R
          X0R = A(J4) + A(J5)
          X0I = A(J4 + 1_IK) + A(J5 + 1_IK)
          X1R = A(J4) - A(J5)
          X1I = A(J4 + 1_IK) - A(J5 + 1_IK)
          X2R = A(J6) + A(J7)
          X2I = A(J6 + 1_IK) + A(J7 + 1)
          X3R = A(J6) - A(J7)
          X3I = A(J6 + 1_IK) - A(J7 + 1_IK)
          Y4R = X0R + X2R
          Y4I = X0I + X2I
          Y6R = X0R - X2R
          Y6I = X0I - X2I
          X0R = X1R - X3I
          X0I = X1I + X3R
          X2R = X1R + X3I
          X2I = X1I - X3R
          Y5R = WN4R * (X0R - X0I)
          Y5I = WN4R * (X0R + X0I)
          Y7R = WN4R * (X2R - X2I)
          Y7I = WN4R * (X2R + X2I)
          X0R = Y1R + Y5R
          X0I = Y1I + Y5I
          A(J1) = WK1R * X0R - WK1I * X0I
          A(J1 + 1_IK) = WK1R * X0I + WK1I * X0R
          X0R = Y1R - Y5R
          X0I = Y1I - Y5I
          A(J5) = WK5R * X0R - WK5I * X0I
          A(J5 + 1_IK) = WK5R * X0I + WK5I * X0R
          X0R = Y3R - Y7I
          X0I = Y3I + Y7R
          A(J3) = WK3R * X0R - WK3I * X0I
          A(J3 + 1_IK) = WK3R * X0I + WK3I * X0R
          X0R = Y3R + Y7I
          X0I = Y3I - Y7R
          A(J7) = WK7R * X0R - WK7I * X0I
          A(J7 + 1_IK) = WK7R * X0I + WK7I * X0R
          A(J) = Y0R + Y4R
          A(J + 1_IK) = Y0I + Y4I
          X0R = Y0R - Y4R
          X0I = Y0I - Y4I
          A(J4) = WK4R * X0R - WK4I * X0I
          A(J4 + 1_IK) = WK4R * X0I + WK4I * X0R
          X0R = Y2R - Y6I
          X0I = Y2I + Y6R
          A(J2) = WK2R * X0R - WK2I * X0I
          A(J2 + 1_IK) = WK2R * X0I + WK2I * X0R
          X0R = Y2R + Y6I
          X0I = Y2I - Y6R
          A(J6) = WK6R * X0R - WK6I * X0I
          A(J6 + 1_IK) = WK6R * X0I + WK6I * X0R
        END DO
      END DO
    END IF
  END SUBROUTINE CFT_MDL_

  ! (FFTW) COMPLEX DISCRETE FOURIER TRANSFORM (POWER OF TWO INPUT LENGTH INPUT)
  ! (SUBROUTINE) FFT_EXTERNAL_(<LENGTH>, <DIRECTION>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <DIRECTION>            -- (IN)     DIRECTION (IK), FFT_FORWARD = +1_IK OR FFT_INVERSE = -1_IK
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <SEQUENCE>             -- (OUT)    CDFT (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., FR_I, FI_I, ...]
  ! void    fft_external_(int*, int*, double*) ;
  SUBROUTINE FFT_EXTERNAL_(LENGTH, DIRECTION, SEQUENCE) &
    BIND(C, NAME = "fft_external_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    INTEGER(IK), INTENT(IN) :: DIRECTION
    REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(INOUT) :: SEQUENCE
    COMPLEX(RK), DIMENSION(LENGTH) :: IN, OUT
    INTEGER(IK) :: PLAN
    IN%RE = SEQUENCE(1_IK:LENGTH:2_IK)
    IN%IM = SEQUENCE(2_IK:LENGTH:2_IK)
    CALL DFFTW_PLAN_DFT_1D(PLAN, LENGTH, IN, OUT, DIRECTION, 64_IK)
    CALL DFFTW_EXECUTE_DFT(PLAN, IN, OUT)
    CALL DFFTW_DESTROY_PLAN(PLAN)
    SEQUENCE(1_IK:LENGTH:2_IK) = OUT%RE
    SEQUENCE(2_IK:LENGTH:2_IK) = OUT%IM
  END SUBROUTINE FFT_EXTERNAL_

  ! (LINEAR) FRACTIONAL COMPLEX DISCRETE FOURIER TRANSFORM (POWER OF TWO INPUT LENGTH INPUT)
  ! (SUBROUTINE) FFRFT_(<LENGTH>, <ARGUMENT>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <ARGUMENT>             -- (IN)     PARAMETER (RK)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <SEQUENCE>             -- (OUT)    CDFT (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., FR_I, FI_I, ...]
  ! void    ffrft_(int*, double*, double*) ;
  SUBROUTINE FFRFT_(LENGTH, ARGUMENT, SEQUENCE) &
    BIND(C, NAME = "ffrft_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), INTENT(IN) :: ARGUMENT
    REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(INOUT) :: SEQUENCE
    INTEGER(IK) :: I
    REAL(RK) :: FACTOR
    REAL(RK), DIMENSION(LENGTH)   :: MUL, COS_MUL, SIN_MUL
    REAL(RK), DIMENSION(4_IK*LENGTH) :: ONE, TWO, TRE
    REAL(RK), DIMENSION(2_IK*LENGTH) :: COPY
    FACTOR = ARGUMENT*ONE_PI/REAL(LENGTH, RK)
    MUL = FACTOR*REAL([(I, I = 0_IK, LENGTH-1_IK, 1_IK)], RK)**2_IK
    COS_MUL = COS(MUL)
    SIN_MUL = SIN(MUL)
    ONE = 0.0_RK
    ONE(1_IK:2_IK*LENGTH:2_IK) = SEQUENCE(1_IK:2_IK*LENGTH:2_IK)*COS_MUL-SEQUENCE(2_IK:2_IK*LENGTH:2_IK)*SIN_MUL
    ONE(2_IK:2_IK*LENGTH:2_IK) = SEQUENCE(1_IK:2_IK*LENGTH:2_IK)*SIN_MUL+SEQUENCE(2_IK:2_IK*LENGTH:2_IK)*COS_MUL
    TWO = 0.0_RK
    TWO(1_IK:2_IK*LENGTH:2_IK) = +COS_MUL
    TWO(2_IK:2_IK*LENGTH:2_IK) = -SIN_MUL
    MUL = -FACTOR*(REAL([(I, I = LENGTH+1_IK, 2_IK*LENGTH, 1_IK)], RK)-1.0_RK-2.0_RK*REAL(LENGTH, RK))**2_IK
    TWO(2_IK*LENGTH+1_IK:4_IK*LENGTH:2_IK) = COS(MUL)
    TWO(2_IK*LENGTH+2_IK:4_IK*LENGTH:2_IK) = SIN(MUL)
    CALL __FFT__(2_IK*LENGTH, FFT_FORWARD, ONE)
    CALL __FFT__(2_IK*LENGTH, FFT_FORWARD, TWO)
    TRE = ONE
    ONE(1_IK:4_IK*LENGTH:2_IK) = TRE(1_IK:4_IK*LENGTH:2_IK)*TWO(1_IK:4_IK*LENGTH:2_IK)-&
      TRE(2_IK:4_IK*LENGTH:2_IK)*TWO(2_IK:4_IK*LENGTH:2_IK)
    ONE(2_IK:4_IK*LENGTH:2_IK) = TRE(1_IK:4_IK*LENGTH:2_IK)*TWO(2_IK:4_IK*LENGTH:2_IK)+&
      TRE(2_IK:4_IK*LENGTH:2_IK)*TWO(1_IK:4_IK*LENGTH:2_IK)
    CALL __FFT__(2_IK*LENGTH, FFT_INVERSE, ONE)
    COPY = 1.0_RK/REAL(2_IK*LENGTH, RK)*ONE(1_IK:2_IK*LENGTH:1_IK)
    SEQUENCE(1_IK:2_IK*LENGTH:2_IK) = COPY(1_IK:2_IK*LENGTH:2_IK)*COS_MUL-COPY(2_IK:2_IK*LENGTH:2_IK)*SIN_MUL
    SEQUENCE(2_IK:2_IK*LENGTH:2_IK) = COPY(1_IK:2_IK*LENGTH:2_IK)*SIN_MUL+COPY(2_IK:2_IK*LENGTH:2_IK)*COS_MUL
  END SUBROUTINE FFRFT_

  ! FACTORIAL
  ! (FUNCTION) FACTORIAL_(<NUMBER>)
  ! <NUMBER>               -- (IN)     NUMBER (IK)
  ! <FACTORIAL_>           -- (OUT)    FACTORIAL OF <N> (RK REAL)
  REAL(RK) FUNCTION FACTORIAL_(NUMBER)
    INTEGER(IK), INTENT(IN) :: NUMBER
    INTEGER(IK) :: I
    FACTORIAL_ = 1.0_RK
    DO I = 1_IK, NUMBER, 1_IK
      FACTORIAL_ = REAL(I, RK)*FACTORIAL_
    END DO
  END FUNCTION FACTORIAL_

  ! WINDOW (COSINE)
  ! (SUBROUTINE) WINDOW_(<LENGTH>, <ORDER>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <ORDER>                -- (IN)     WINDOW ORDER (IK)
  ! <WINDOW>               -- (OUT)    WINDOW (RK ARRAY OF LENGTH = <LENGTH>)
  ! void    window_(int*, int*, double*) ;
  SUBROUTINE WINDOW_(LENGTH, ORDER, WINDOW) &
    BIND(C, NAME = "window_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    INTEGER(IK), INTENT(IN) :: ORDER
    REAL(RK), INTENT(OUT), DIMENSION(LENGTH) :: WINDOW
    INTEGER(IK) :: I
    REAL(RK) :: FACTOR
    FACTOR = 2.0_RK**ORDER*FACTORIAL_(ORDER)**2_IK/FACTORIAL_(2_IK*ORDER)
    WINDOW = FACTOR*(1.0_RK+COS(TWO_PI*(1.0_RK/REAL(LENGTH,RK)*REAL([(I, I = 0_IK, LENGTH-1_IK, 1_IK)], RK)-0.5_RK)))**ORDER
  END SUBROUTINE WINDOW_

  ! MINLOC
  ! (FUNCTION) MINLOC_(<SEQUENCE>)
  INTEGER(IK) FUNCTION MINLOC_(SEQUENCE, EMPTY)
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

  ! MAXLOC
  ! (FUNCTION) MAXLOC_(<SEQUENCE>)
  INTEGER(IK) FUNCTION MAXLOC_(SEQUENCE, EMPTY)
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

  ! SORT (BUBBLE, DESCENDING)
  ! (SUBROUTINE) SORT_BUBBLE_(<LENGTH>, <SEQUENCE>, <FST>, <LST>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (IN)     (UNSORTED) SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (OUT)    (SORTED, DESCENDING) SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  SUBROUTINE SORT_BUBBLE_(LENGTH, SEQUENCE, FST, LST)
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

  ! SORT (QUICK, DESCENDING)
  ! (SUBROUTINE) SORT_QUICK_(<LENGTH>, <SEQUENCE>, <FST>, <LST>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (IN)     (UNSORTED) SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (OUT)    (SORTED, DESCENDING) SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  RECURSIVE SUBROUTINE SORT_QUICK_(LENGTH, SEQUENCE, FST, LST)
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

  ! PEAK LIST
  ! (SUBROUTINE) PEAK_LIST_(<LENGTH>, <SEQUENCE>, <PEAK_LIST>)
  ! PEAK_WIDTH             -- (GLOBAL) PEAK WIDTH (IK)
  ! PEAK_LEVEL             -- (GLOBAL) PEAK LEVEL THRESHOLD (RK)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (IN)     SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  ! <PEAK_LIST>            -- (OUT)    PEAK LIST (IK ARRAY OF LENGTH = <LENGTH>), VALUES OF ONE CORRESPOND TO PEAK LOCATIONS
  SUBROUTINE PEAK_LIST_(LENGTH, SEQUENCE, PEAK_LIST)
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(LENGTH), INTENT(IN) :: SEQUENCE
    INTEGER(IK), DIMENSION(LENGTH), INTENT(OUT) :: PEAK_LIST
    INTEGER(IK) :: I
    INTEGER(IK), DIMENSION(LENGTH-(1_IK+2_IK*PEAK_WIDTH)) :: FST
    INTEGER(IK), DIMENSION(LENGTH-2_IK) :: LST
    REAL(RK) :: TOTAL
    REAL(RK) :: LAST
    REAL(RK), DIMENSION(LENGTH-2_IK*PEAK_WIDTH) :: LOCAL
    LOGICAL :: THIS, NEXT
    FST = 0_IK
    TOTAL = SUM(SEQUENCE(1_IK:2_IK*PEAK_WIDTH))
    LAST = 0.0_RK
    LOCAL = SEQUENCE(1_IK+PEAK_WIDTH:LENGTH-PEAK_WIDTH) - PEAK_LEVEL
    DO I = 0_IK, LENGTH-(1_IK+2_IK*PEAK_WIDTH)-1_IK, 1_IK
      TOTAL = TOTAL+SEQUENCE(I+1_IK+2_IK*PEAK_WIDTH)-LAST
      LAST = SEQUENCE(I+1_IK)
      LOCAL(I+1_IK) = LOCAL(I+1_IK)-TOTAL/REAL(1_IK+2_IK*PEAK_WIDTH,RK)
      IF (LOCAL(I+1_IK) >= 0.0_RK) THEN
        FST(I+1_IK) = 1_IK
      END IF
    END DO
    PEAK_LIST = [[(0_IK, I = 1_IK, PEAK_WIDTH, 1_IK)], FST, [(0_IK, I = 1_IK, PEAK_WIDTH, 1_IK)]]
    LST = 0_IK
    DO I = 1_IK, LENGTH-2_IK, 1_IK
      THIS = INT(SIGN(1.0_RK, SEQUENCE(I+1_IK)-SEQUENCE(I)), IK) == +1_IK
      NEXT = INT(SIGN(1.0_RK, SEQUENCE(I+2_IK)-SEQUENCE(I+1_IK)), IK) == -1_IK
      IF (THIS .AND. NEXT) THEN
        LST(I) = 1_IK
      END IF
    END DO
    PEAK_LIST = PEAK_LIST*[0_IK, LST, 0_IK]
  END SUBROUTINE PEAK_LIST_

  ! TOTAL NUMBER OF PEAKS
  ! (FUNCTION) PEAK_COUNT_(<LENGTH>, <PEAK_LIST>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <PEAK_LIST>            -- (IN)     PEAK LIST (IK ARRAY OF LENGTH <LENGTH>)
  ! <PEAK_COUNT_>          -- (OUT)    TOTAL NUMBER OF PEAKS (IK)
  INTEGER(IK) FUNCTION PEAK_COUNT_(LENGTH, PEAK_LIST)
    INTEGER(IK), INTENT(IN) :: LENGTH
    INTEGER(IK), DIMENSION(LENGTH), INTENT(IN) :: PEAK_LIST
    PEAK_COUNT_ = SUM(PEAK_LIST)
  END FUNCTION PEAK_COUNT_

  ! DETECT SEVERAL PEAKS (LIST OF ORDERED PEAK POSITIONS)
  ! (SUBROUTINE) PEAK_DETECT_(<LENGTH>, <SEQUENCE>, <PEAK_LENGTH>, <PEAK_ORDERED_LIST>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (IN)     SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  ! <PEAK_LENGTH>          -- (IN)     NUMBER OF PEAKS TO FIND (IK)
  ! <PEAK_ORDERED_LIST>    -- (OUT)    PEAK POSITIONS (IK ARRAY OF LENGTH = <PEAK_LENGTH>)
  SUBROUTINE PEAK_DETECT_(LENGTH, SEQUENCE, PEAK_LENGTH, PEAK_ORDERED_LIST)
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(LENGTH), INTENT(IN) :: SEQUENCE
    INTEGER(IK), INTENT(IN) :: PEAK_LENGTH
    INTEGER(IK), DIMENSION(PEAK_LENGTH), INTENT(OUT) :: PEAK_ORDERED_LIST
    INTEGER(IK), DIMENSION(LENGTH) :: PEAK_LIST
    INTEGER(IK) :: LIMIT
    INTEGER(IK) :: I, J
    INTEGER(IK) :: POSITION
    REAL(RK), DIMENSION(LENGTH) :: LOCAL
    PEAK_LIST = 0_IK
    CALL PEAK_LIST_(LENGTH, SEQUENCE, PEAK_LIST)
    PEAK_ORDERED_LIST = 1_IK
    LIMIT = PEAK_COUNT_(LENGTH, PEAK_LIST)
    POSITION = 1_IK
    DO I = 1_IK, LENGTH, 1_IK
      IF (PEAK_LIST(I) == 1_IK) THEN
        LOCAL(POSITION) = SEQUENCE(I)
        POSITION = POSITION + 1_IK
      END IF
    END DO
    CALL __SORT__(LIMIT, LOCAL, 1_IK, LIMIT)
    DO I = 1_IK, MIN(PEAK_LENGTH, LIMIT), 1_IK
      DO J = 1_IK, LENGTH, 1_IK
        IF (LOCAL(I) == SEQUENCE(J)) PEAK_ORDERED_LIST(I) = J
      END DO
    END DO
  END SUBROUTINE PEAK_DETECT_

  ! PEAK (RANKED)
  ! (FUNCTION) PEAK_(<LENGTH>, <SEQUENCE>, <PEAK_ID>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (IN)     SEQUENCE (RK ARRAY OF LENGTH <LENGTH>)
  ! <PEAK_ID>              -- (IN)     PEAK RANK (IK)
  ! <PEAK_>                -- (OUT)    PEAK POSITION (IK)
  ! int     peak_(int*, double*, int*) ;
  INTEGER(IK) FUNCTION PEAK_(LENGTH, SEQUENCE, PEAK_ID) &
    BIND(C, NAME = "peak_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(LENGTH), INTENT(IN) :: SEQUENCE
    INTEGER(IK), INTENT(IN) :: PEAK_ID
    INTEGER(IK), DIMENSION(LENGTH) :: PEAK_LIST
    CALL PEAK_DETECT_(LENGTH, SEQUENCE, LENGTH, PEAK_LIST)
    PEAK_ = PEAK_LIST(PEAK_ID)
  END FUNCTION PEAK_

  ! FREQUENCY ESTIMATION
  ! (FUNCTION) FREQUENCY_(<FLAG>, <PEAK>, <LENGTH>, <TOTAL>, <WINDOW>, <SEQUENCE>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <PEAK>                 -- (IN)     PEAK NUMBER TO USE (IK), <PEAK> = 0 USE MAXIMUM BIN, <PEAK> = +N USE N'TH PEAK, <PEAK> = -N USE N'TH PEAK AND FFT APPROXIMATION
  ! <LENGTH>               -- (IN)     INPUT SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <FREQUENCY_>           -- (OUT)    FREQUENCY ESTIMATION (RK)
  ! double  frequency_(int*, int*, int*, double*, double*, double*) ;
  REAL(RK) FUNCTION FREQUENCY_(FLAG, PEAK, LENGTH, TOTAL, WINDOW, SEQUENCE) &
    BIND(C, NAME = "frequency_")
    INTEGER(IK), INTENT(IN):: FLAG
    INTEGER(IK), INTENT(IN) :: PEAK
    INTEGER(IK), INTENT(IN):: LENGTH
    REAL(RK), INTENT(IN) :: TOTAL
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
    REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    REAL(RK), DIMENSION(2_IK*LENGTH) :: LOCAL, FOURIER
    INTEGER(IK) :: FST, CND
    REAL(RK) :: FACTOR
    REAL(RK), DIMENSION(LENGTH) :: MUL, COS_MUL, SIN_MUL
    INTEGER(IK) :: I
    ! REMOVE (WEIGHTED) MEAN AND APPLY WINDOW
    LOCAL(1_IK:2_IK*LENGTH:2_IK) = (SEQUENCE(1_IK:2_IK*LENGTH:2_IK)-SUM(WINDOW*SEQUENCE(1_IK:2_IK*LENGTH:2_IK))/TOTAL)*WINDOW
    LOCAL(2_IK:2_IK*LENGTH:2_IK) = (SEQUENCE(2_IK:2_IK*LENGTH:2_IK)-SUM(WINDOW*SEQUENCE(2_IK:2_IK*LENGTH:2_IK))/TOTAL)*WINDOW
    ! FFT APPROXIMATION
    FOURIER = LOCAL
    CALL __FFT__(LENGTH, FFT_FORWARD, FOURIER)
    MUL = LOG10(SQRT(FOURIER(1_IK:2_IK*LENGTH:2_IK)**2_IK+FOURIER(2_IK:2_IK*LENGTH:2_IK)**2_IK)+1.E-16_RK)
    IF (PEAK == 0_IK) THEN
      FST = INT(__MAXLOC__(MUL(1_IK:LENGTH/(2_IK-FLAG):1_IK), 1_IK), IK)
    ELSE
      FST = PEAK_(LENGTH/(2_IK-FLAG), MUL(1_IK:LENGTH/(2_IK-FLAG):1_IK), ABS(PEAK))
    END IF
    IF (PEAK < 0_IK) THEN
      FREQUENCY_ = REAL(FST-1_IK, RK)/REAL(LENGTH, RK)
      RETURN
    END IF
    ! MODULATE SEQUENCE
    FACTOR = TWO_PI*REAL(FST-2_IK, RK)/REAL(LENGTH, RK)
    MUL = FACTOR*REAL([(I, I = 0_IK, LENGTH-1_IK, 1_IK)], RK)
    COS_MUL = COS(MUL)
    SIN_MUL = SIN(MUL)
    FOURIER = LOCAL
    LOCAL(1_IK:2_IK*LENGTH:2_IK) = FOURIER(1_IK:2_IK*LENGTH:2_IK)*COS_MUL-FOURIER(2_IK:2_IK*LENGTH:2_IK)*SIN_MUL
    LOCAL(2_IK:2_IK*LENGTH:2_IK) = FOURIER(1_IK:2_IK*LENGTH:2_IK)*SIN_MUL+FOURIER(2_IK:2_IK*LENGTH:2_IK)*COS_MUL
    ! FFRFT APPROXIMATION
    CALL FFRFT_(LENGTH, 2.0_RK/REAL(LENGTH, RK), LOCAL)
    MUL = LOG10(SQRT(LOCAL(1_IK:2_IK*LENGTH:2_IK)**2_IK+LOCAL(2_IK:2_IK*LENGTH:2_IK)**2_IK)+1.E-16_RK)
    CND = INT(__MAXLOC__(MUL, 1), IK)
    ! PARABOLA APPROXIMATION
    FREQUENCY_ = REAL(CND, RK)-0.5_RK+(MUL(-1_IK+CND)-MUL(CND))/(MUL(-1_IK+CND)-2.0_RK*MUL(CND)+MUL(1_IK+CND))
    ! RESULT
    FREQUENCY_ = (REAL(FST, RK)-2.0_RK+2.0_RK*(FREQUENCY_-1.0_RK)/REAL(LENGTH, RK))/REAL(LENGTH, RK)
  END FUNCTION FREQUENCY_

  ! SIGNAL DECOMPOSITION (WITHOUT CONSTANT PART)
  ! (SUBROUTINE) DECOMPOSITION(<FLAG>, <METHOD>, <LENGTH>, <TOTAL>, <WINDOW>, <SEQUENCE>, <LOOP>, <FREQUENCY>, <COS_AMP>, <SIN_AMP>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <METHOD>               -- (IN)     METHOD TO USE (IK), <METHOD> = 0 USE MAXIMUM BIN, <METHOD> = +1 USE PEAKS, <METHOD> = -1 USE PEAKS AND FFT APPROXIMATION
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <LOOP>                 -- (IN)     NUMBER OF ITERATIONS (IK)
  ! <FREQUENCY>            -- (OUT)    FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <COS_AMP>              -- (OUT)    COS AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <SIN_AMP>              -- (OUT)    SIN AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! void    decomposition_(int*, int*, int*, double*, double*, double*, int*, double*, double*, double*) ;
  SUBROUTINE DECOMPOSITION_(FLAG, METHOD, LENGTH, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP) &
    BIND(C, NAME = "decomposition_")
    INTEGER(IK), INTENT(IN):: FLAG
    INTEGER(IK), INTENT(IN):: METHOD
    INTEGER(IK), INTENT(IN):: LENGTH
    REAL(RK), INTENT(IN) :: TOTAL
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
    REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    INTEGER(IK), INTENT(IN) :: LOOP
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: FREQUENCY
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: COS_AMP
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: SIN_AMP
    INTEGER(IK) :: I
    REAL(RK), DIMENSION(2_IK*LENGTH) :: LOCAL
    REAL(RK), DIMENSION(LENGTH) :: LIST
    REAL(RK), DIMENSION(2_IK*LENGTH) :: DELTA
    LOCAL(1_IK:2_IK*LENGTH:2_IK) = (SEQUENCE(1_IK:2_IK*LENGTH:2_IK)-SUM(WINDOW*SEQUENCE(1_IK:2_IK*LENGTH:2_IK))/TOTAL)
    LOCAL(2_IK:2_IK*LENGTH:2_IK) = (SEQUENCE(2_IK:2_IK*LENGTH:2_IK)-SUM(WINDOW*SEQUENCE(2_IK:2_IK*LENGTH:2_IK))/TOTAL)
    LIST = TWO_PI*REAL([(I, I=1_IK, LENGTH, 1_IK)], RK)
    DO I = 1_IK, LOOP, 1_IK
      IF (METHOD == DECOMPOSITION_MAX) THEN
        FREQUENCY(I) = FREQUENCY_(FLAG, 0_IK, LENGTH, TOTAL, WINDOW, LOCAL)
      ELSE IF(METHOD == DECOMPOSITION_PEAK) THEN
        FREQUENCY(I) = FREQUENCY_(FLAG, +I, LENGTH, TOTAL, WINDOW, SEQUENCE)
      ELSE IF(METHOD == DECOMPOSITION_PEAK_FFT) THEN
        FREQUENCY(I) = FREQUENCY_(FLAG, -I, LENGTH, TOTAL, WINDOW, SEQUENCE)
      END IF
      DELTA(1_IK:2_IK*LENGTH:2_IK) = +COS(FREQUENCY(I)*LIST)
      DELTA(2_IK:2_IK*LENGTH:2_IK) = -SIN(FREQUENCY(I)*LIST)
      COS_AMP(I) = SUM(SEQUENCE(1_IK:2_IK*LENGTH:2_IK)*DELTA(1_IK:2_IK*LENGTH:2_IK)*WINDOW)+&
        SUM(SEQUENCE(2_IK:2_IK*LENGTH:2_IK)*DELTA(2_IK:2_IK*LENGTH:2_IK)*WINDOW)
      SIN_AMP(I) = SUM(SEQUENCE(2_IK:2_IK*LENGTH:2_IK)*DELTA(1_IK:2_IK*LENGTH:2_IK)*WINDOW)-&
        SUM(SEQUENCE(1_IK:2_IK*LENGTH:2_IK)*DELTA(2_IK:2_IK*LENGTH:2_IK)*WINDOW)
      COS_AMP(I) = 1.0_RK/TOTAL*COS_AMP(I)
      SIN_AMP(I) = 1.0_RK/TOTAL*SIN_AMP(I)
      LOCAL(1_IK:2_IK*LENGTH:2_IK) = LOCAL(1_IK:2_IK*LENGTH:2_IK)-&
        COS_AMP(I)*DELTA(1_IK:2_IK*LENGTH:2_IK)+SIN_AMP(I)*DELTA(2_IK:2_IK*LENGTH:2_IK)
      LOCAL(2_IK:2_IK*LENGTH:2_IK) = LOCAL(2_IK:2_IK*LENGTH:2_IK)-&
        COS_AMP(I)*DELTA(2_IK:2_IK*LENGTH:2_IK)-SIN_AMP(I)*DELTA(1_IK:2_IK*LENGTH:2_IK)
    END DO
    IF (FLAG == 0_IK) THEN
      COS_AMP = 2.0_RK*COS_AMP
      SIN_AMP = 2.0_RK*SIN_AMP
    END IF
  END SUBROUTINE DECOMPOSITION_

  ! FREQUENCY LIST (PERFORM DECOMPOSITION AND RETURN LIST OF FREQUENCIES)
  ! (SUBROUTINE) FREQUENCY_LIST_(<FLAG>, <METHOD>, <LENGTH>, <TOTAL>, <WINDOW>, <SEQUENCE>, <LOOP>, <FREQUENCY>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <METHOD>               -- (IN)     METHOD TO USE (IK), <METHOD> = 0 USE MAXIMUM BIN, <METHOD> = +1 USE PEAKS, <METHOD> = -1 USE PEAKS AND FFT APPROXIMATION
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK), POWER OF TWO
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <LOOP>                 -- (IN)     NUMBER OF ITERATIONS (IK)
  ! <FREQUENCY>            -- (OUT)    FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! void    frequency_list_(int*, int*, int*, double*, double*, double*, int*, double*) ;
  SUBROUTINE FREQUENCY_LIST_(FLAG, METHOD, LENGTH, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY) &
    BIND(C, NAME = "frequency_list_")
    INTEGER(IK), INTENT(IN):: FLAG
    INTEGER(IK), INTENT(IN) :: METHOD
    INTEGER(IK), INTENT(IN):: LENGTH
    REAL(RK), INTENT(IN) :: TOTAL
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
    REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    INTEGER(IK), INTENT(IN) :: LOOP
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: FREQUENCY
    REAL(RK), DIMENSION(LOOP) :: COS_AMP
    REAL(RK), DIMENSION(LOOP) :: SIN_AMP
    CALL DECOMPOSITION_(FLAG, METHOD, LENGTH, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP)
  END SUBROUTINE FREQUENCY_LIST_

  ! AMPLITUDE LIST (COMPUTE AMPLITUDES FOR LIST OF GIVEN FREQUENCIES)
  ! (SUBROUTINE) AMPLITUDE_LIST_(<FLAG>, <LENGTH>, <TOTAL>, <WINDOW>, <SEQUENCE>, <LOOP>, <FREQUENCY>, <COS_AMP>, <SIN_AMP>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK), POWER OF TWO
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <LOOP>                 -- (IN)     NUMBER OF ITERATIONS (IK)
  ! <FREQUENCY>            -- (IN)     FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <COS_AMP>              -- (OUT)    COS AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <SIN_AMP>              -- (OUT)    SIN AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! void    amplitude_list_(int*, int*, double*, double*, double*, int*, double*, double*, double*) ;
  SUBROUTINE AMPLITUDE_LIST_(FLAG, LENGTH, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP) &
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
    LOCAL(1_IK:2_IK*LENGTH:2_IK) = (SEQUENCE(1_IK:2_IK*LENGTH:2_IK)-SUM(WINDOW*SEQUENCE(1_IK:2_IK*LENGTH:2_IK))/TOTAL)
    LOCAL(2_IK:2_IK*LENGTH:2_IK) = (SEQUENCE(2_IK:2_IK*LENGTH:2_IK)-SUM(WINDOW*SEQUENCE(2_IK:2_IK*LENGTH:2_IK))/TOTAL)
    LIST = TWO_PI*REAL([(I, I=1_IK, LENGTH, 1_IK)], RK)
    DO I = 1_IK, LOOP, 1_IK
      DELTA(1_IK:2_IK*LENGTH:2_IK) = +COS(FREQUENCY(I)*LIST)
      DELTA(2_IK:2_IK*LENGTH:2_IK) = -SIN(FREQUENCY(I)*LIST)
      COS_AMP(I) = SUM(SEQUENCE(1_IK:2_IK*LENGTH:2_IK)*DELTA(1_IK:2_IK*LENGTH:2_IK)*WINDOW)+&
        SUM(SEQUENCE(2_IK:2_IK*LENGTH:2_IK)*DELTA(2_IK:2_IK*LENGTH:2_IK)*WINDOW)
      SIN_AMP(I) = SUM(SEQUENCE(2_IK:2_IK*LENGTH:2_IK)*DELTA(1_IK:2_IK*LENGTH:2_IK)*WINDOW)-&
          SUM(SEQUENCE(1_IK:2_IK*LENGTH:2_IK)*DELTA(2_IK:2_IK*LENGTH:2_IK)*WINDOW)
      COS_AMP(I) = 1.0_RK/REAL(LENGTH, RK)*COS_AMP(I)
      SIN_AMP(I) = 1.0_RK/REAL(LENGTH, RK)*SIN_AMP(I)
      LOCAL(1_IK:2_IK*LENGTH:2_IK) = LOCAL(1_IK:2_IK*LENGTH:2_IK)-&
        COS_AMP(I)*DELTA(1_IK:2_IK*LENGTH:2_IK)+SIN_AMP(I)*DELTA(2_IK:2_IK*LENGTH:2_IK)
      LOCAL(2_IK:2_IK*LENGTH:2_IK) = LOCAL(2_IK:2_IK*LENGTH:2_IK)-&
        COS_AMP(I)*DELTA(2_IK:2_IK*LENGTH:2_IK)-SIN_AMP(I)*DELTA(1_IK:2_IK*LENGTH:2_IK)
    END DO
    IF (FLAG == 0_IK) THEN
      COS_AMP = 2.0_RK*COS_AMP
      SIN_AMP = 2.0_RK*SIN_AMP
    END IF
  END SUBROUTINE AMPLITUDE_LIST_

  ! ESTIMATE AMPLITUDE FOR GIVEN FREQUENCY (DTFT SPECTRA COMPUTATION)
  ! (SUBROUTINE) AMPLITUDE(<LENGTH>, <TOTAL>, <WINDOW>, <SEQUENCE>, <FREQUENCY>, <COS_AMP>, <SIN_AMP>, <AMP>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <FREQUENCY>            -- (IN)     FREQUENCY (RK)
  ! <COS_AMP>              -- (OUT)    COS AMPLITUDE (RK)
  ! <SIN_AMP>              -- (OUT)    SIN AMPLITUDE (RK)
  ! <AMP>                  -- (OUT)    ABS AMPLITUDE (RK)
  ! void    amplitude_(int*, double*, double*, double*, double*, double*, double*, double*) ;  
  SUBROUTINE AMPLITUDE_(LENGTH, TOTAL, WINDOW, SEQUENCE, FREQUENCY, COS_AMP, SIN_AMP, AMP) &
    BIND(C, NAME = "amplutude_")
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
    LOCAL(1_IK:2_IK*LENGTH:2_IK) = (SEQUENCE(1_IK:2_IK*LENGTH:2_IK)-SUM(WINDOW*SEQUENCE(1_IK:2_IK*LENGTH:2_IK))/TOTAL)
    LOCAL(2_IK:2_IK*LENGTH:2_IK) = (SEQUENCE(2_IK:2_IK*LENGTH:2_IK)-SUM(WINDOW*SEQUENCE(2_IK:2_IK*LENGTH:2_IK))/TOTAL)
    LIST = TWO_PI*REAL([(I, I=1_IK, LENGTH, 1_IK)], RK)
    DELTA(1_IK:2_IK*LENGTH:2_IK) = +COS(FREQUENCY*LIST)
    DELTA(2_IK:2_IK*LENGTH:2_IK) = -SIN(FREQUENCY*LIST)
    COS_AMP = SUM(SEQUENCE(1_IK:2_IK*LENGTH:2_IK)*DELTA(1_IK:2_IK*LENGTH:2_IK)*WINDOW)+&
      SUM(SEQUENCE(2_IK:2_IK*LENGTH:2_IK)*DELTA(2_IK:2_IK*LENGTH:2_IK)*WINDOW)
    SIN_AMP = SUM(SEQUENCE(2_IK:2_IK*LENGTH:2_IK)*DELTA(1_IK:2_IK*LENGTH:2_IK)*WINDOW)-&
      SUM(SEQUENCE(1_IK:2_IK*LENGTH:2_IK)*DELTA(2_IK:2_IK*LENGTH:2_IK)*WINDOW)
    COS_AMP = 1.0_RK/TOTAL*COS_AMP
    SIN_AMP = 1.0_RK/TOTAL*SIN_AMP
    AMP = SQRT(COS_AMP**2_IK+SIN_AMP**2_IK)
  END SUBROUTINE AMPLITUDE_

  ! SVD (DGESVD)
  ! (SUBROUTINE) SVD_(<NR>, <NC>, <MATRIX>(<NR>, <NC>), <SVD_LIST>(MIN(<NR>, <NC>)), <U_MATRIX>(<NR>, <NR>), <V_MATRIX>(<NC>, <NC>))
  SUBROUTINE SVD_(NR, NC, MATRIX, SVD_LIST, U_MATRIX, V_MATRIX)
    INTEGER(IK), INTENT(IN) :: NR
    INTEGER(IK), INTENT(IN) :: NC
    REAL(RK), DIMENSION(NR, NC), INTENT(IN) :: MATRIX
    REAL(RK), DIMENSION(MIN(NR, NC)), INTENT(OUT) :: SVD_LIST
    REAL(RK), DIMENSION(NR, NR), INTENT(OUT) :: U_MATRIX
    REAL(RK), DIMENSION(NC, NC), INTENT(OUT) :: V_MATRIX
    REAL(8), DIMENSION(2*MAX(1, 3*MIN(INT(NR), INT(NC))+MAX(INT(NR), INT(NC)), 5*MIN(INT(NR), INT(NC)))) :: WORK
    INTEGER :: WORK_SIZE
    INTEGER :: INFO
    REAL(8), DIMENSION(NR, NC) :: COPY
    REAL(8), DIMENSION(MIN(NR, NC)) :: D_SVD_LIST
    REAL(8), DIMENSION(NR, NR) :: D_U_MATRIX
    REAL(8), DIMENSION(NC, NC) :: D_V_MATRIX
    COPY = REAL(MATRIX, 8)
    WORK_SIZE = SIZE(WORK)
    CALL DGESVD('A','A',INT(NR),INT(NC),COPY,INT(NR),D_SVD_LIST,D_U_MATRIX,INT(NR),D_V_MATRIX,INT(NC),WORK,WORK_SIZE,INFO)
    SVD_LIST = REAL(D_SVD_LIST, RK)
    U_MATRIX = REAL(D_U_MATRIX, RK)
    V_MATRIX = REAL(D_V_MATRIX, RK)
    V_MATRIX = TRANSPOSE(V_MATRIX)
    WHERE (SVD_LIST < SVD_LEVEL) SVD_LIST = 0.0_RK
  END SUBROUTINE SVD_

  ! SVD LIST (DGESVD)
  ! (SUBROUTINE) SVD_LIST_(<NR>, <NC>, <MATRIX>(<NR>, <NC>), <SVD_LIST>(MIN(<NR>, <NC>)))
  SUBROUTINE SVD_LIST_(NR, NC, MATRIX, SVD_LIST)
    INTEGER(IK), INTENT(IN) :: NR
    INTEGER(IK), INTENT(IN) :: NC
    REAL(RK), DIMENSION(NR, NC), INTENT(IN) :: MATRIX
    REAL(RK), DIMENSION(MIN(NR, NC)), INTENT(OUT) :: SVD_LIST
    REAL(8), DIMENSION(2*MAX(1, 3*MIN(INT(NR), INT(NC))+MAX(INT(NR), INT(NC)), 5*MIN(INT(NR), INT(NC)))) :: WORK
    INTEGER :: WORK_SIZE
    INTEGER :: INFO
    REAL(8), DIMENSION(NR, NC) :: COPY
    REAL(8), DIMENSION(MIN(NR, NC)) :: D_SVD_LIST
    REAL(8), DIMENSION(NR, NR) :: D_U_MATRIX
    REAL(8), DIMENSION(NC, NC) :: D_V_MATRIX
    COPY = REAL(MATRIX, 8)
    WORK_SIZE = SIZE(WORK)
    CALL DGESVD('N','N',INT(NR),INT(NC),COPY,INT(NR),D_SVD_LIST,D_U_MATRIX,INT(NR),D_V_MATRIX,INT(NC),WORK,WORK_SIZE,INFO)
    SVD_LIST = REAL(D_SVD_LIST, RK)
    WHERE (SVD_LIST < SVD_LEVEL) SVD_LIST = 0.0_RK
  END SUBROUTINE SVD_LIST_

  ! LEAST SQUARES
  ! (SUBROUTINE) LEAST_SQUARES_(<NR>, <NC>, <MATRIX>(<NR>, <NC>), <VECTOR>(<NR>), <SOLUTION>(<NC>))
  SUBROUTINE LEAST_SQUARES_(NR, NC, MATRIX, VECTOR, SOLUTION)
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
    REAL(8) :: D_U_MATRIX(NR, NR), D_V_MATRIX(NC, NC), D_COPY(NR, NC)
    REAL(8) :: V1(NR), V2(NR), V3(NC), V4(NC)
    CALL SVD_(NR, NC, MATRIX, SVD_LIST, U_MATRIX, V_MATRIX)
    COPY = 0.0_RK
    DO I = 1_IK, INT(SIZE(SVD_LIST), IK), 1_IK
      IF (SVD_LIST(I) >= SVD_LEVEL) COPY(I, I) = 1.0_RK/SVD_LIST(I)
    END DO
    V1 = REAL(VECTOR, 8)
    D_U_MATRIX = REAL(U_MATRIX, 8)
    D_V_MATRIX = REAL(V_MATRIX, 8)
    D_COPY = REAL(COPY, 8)
    CALL DGEMV('T',INT(NR),INT(NR),1.0_8,D_U_MATRIX,INT(NR),V1,1,0.0_8,V2,1)
    CALL DGEMV('T',INT(NR),INT(NC),1.0_8,D_COPY,INT(NR),V2,1,0.0_8,V3,1)
    CALL DGEMV('N',INT(NC),INT(NC),1.0_8,D_V_MATRIX,INT(NC),V3,1,0.0_8,V4,1)
    SOLUTION = REAL(V4, RK)
  END SUBROUTINE LEAST_SQUARES_

  ! FIT
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
  SUBROUTINE FIT_(LENGTH, SEQUENCE, LOOP, FREQUENCY, MEAN, COS_AMP, SIN_AMP, ERROR) &
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

  ! MATRIX
  SUBROUTINE MATRIX_(LENGTH, SEQUENCE, MATRIX)
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(LENGTH), INTENT(IN) :: SEQUENCE
    REAL(RK), DIMENSION(LENGTH/2_IK+1_IK, LENGTH/2_IK), INTENT(OUT) :: MATRIX
    INTEGER(IK) :: I
    DO I = 1_IK, LENGTH/2_IK+1_IK, 1_IK
      MATRIX(I, :) = SEQUENCE(I:I-1_IK+LENGTH/2_IK)
    END DO
  END SUBROUTINE MATRIX_

  ! SEQUENCE
  SUBROUTINE SEQUENCE_(LENGTH, SEQUENCE, MATRIX)
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(LENGTH), INTENT(OUT) :: SEQUENCE
    REAL(RK), DIMENSION(LENGTH/2_IK+1_IK, LENGTH/2_IK), INTENT(IN) :: MATRIX
    SEQUENCE = [MATRIX(1_IK, :), MATRIX(LENGTH/2_IK+1_IK, :)]
  END SUBROUTINE SEQUENCE_

  ! FILTER
  ! (SUBROUTINE) FILTER(<LENGTH>, <SEQUENCE>, <LIMIT>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <SEQUENCE>             -- (INOUT)  SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  ! <LIMIT>                -- (IN)     NUMBER OF SINGULAR VALUES TO KEEP (IK)
  ! void    filter_(int*, double*, int*) ;
  SUBROUTINE FILTER_(LENGTH, SEQUENCE, LIMIT, SVD_LIST) &
    BIND(C, NAME = "filter_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(LENGTH), INTENT(INOUT) :: SEQUENCE
    REAL(RK), DIMENSION(LENGTH/2_IK), INTENT(OUT) :: SVD_LIST
    INTEGER(IK), INTENT(IN) :: LIMIT
    REAL(RK), DIMENSION(LENGTH/2_IK+1_IK, LENGTH/2_IK) :: MATRIX
    REAL(RK), DIMENSION(LENGTH/2_IK+1_IK, LENGTH/2_IK+1_IK) :: U_MATRIX
    REAL(RK), DIMENSION(LENGTH/2_IK, LENGTH/2_IK) :: V_MATRIX
    INTEGER(IK) :: I
    CALL MATRIX_(LENGTH, SEQUENCE, MATRIX)
    CALL SVD_(LENGTH/2_IK+1_IK, LENGTH/2_IK, MATRIX, SVD_LIST, U_MATRIX, V_MATRIX)
    MATRIX = 0.0_RK
    DO I = 1_IK, LIMIT, 1_IK
      MATRIX(I, I) = SVD_LIST(I)
    END DO
    MATRIX = MATMUL(U_MATRIX, MATMUL(MATRIX, TRANSPOSE(V_MATRIX)))
    CALL SEQUENCE_(LENGTH, SEQUENCE, MATRIX)
  END SUBROUTINE FILTER_

  ! COMPUTE DATA TABLE
  ! (SUBROUTINE) COMPUTE_TABLE_(<LENGTH>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! void    compute_table_(int*) ;
  SUBROUTINE COMPUTE_TABLE_(LENGTH) &
    BIND(C, NAME = "compute_table_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    ALLOCATE(BANK%BIT_FFT(2_IK+INT(SQRT(REAL(LENGTH/2_IK, RK)), IK)))
    ALLOCATE(BANK%BIT_FFRFT(2_IK+INT(SQRT(REAL(LENGTH/1_IK, RK)), IK)))
    ALLOCATE(BANK%TRIG_FFT(LENGTH/2_IK))
    ALLOCATE(BANK%TRIG_FFRFT(LENGTH/1_IK))
    ALLOCATE(BANK%COS_FST(LENGTH))
    ALLOCATE(BANK%SIN_FST(LENGTH))
    ALLOCATE(BANK%COS_LST(LENGTH))
    ALLOCATE(BANK%SIN_LST(LENGTH))
    CALL MAKE_FFT_DATA__(LENGTH, BANK%BIT_FFT, BANK%TRIG_FFT)
    CALL MAKE_FFT_DATA__(2_IK*LENGTH, BANK%BIT_FFRFT, BANK%TRIG_FFRFT)
    CALL MAKE_FFRFT_DATA__(LENGTH, BANK%COS_FST, BANK%SIN_FST, BANK%COS_LST, BANK%SIN_LST)    
  END SUBROUTINE COMPUTE_TABLE_

  ! DESTROY DATA TABLE
  ! (SUBROUTINE) DESTROY_TABLE_()
  ! void    destroy_table_(int*) ;
  SUBROUTINE DESTROY_TABLE_() &
    BIND(C, NAME = "destroy_table_")
    DEALLOCATE(BANK%BIT_FFT)
    DEALLOCATE(BANK%BIT_FFRFT)
    DEALLOCATE(BANK%TRIG_FFT)
    DEALLOCATE(BANK%TRIG_FFRFT)
    DEALLOCATE(BANK%COS_FST)
    DEALLOCATE(BANK%SIN_FST)
    DEALLOCATE(BANK%COS_LST)
    DEALLOCATE(BANK%SIN_LST)
  END SUBROUTINE DESTROY_TABLE_  

  SUBROUTINE MAKE_FFT_DATA__(LENGTH, IP, WORK)
    INTEGER(IK), INTENT(IN) :: LENGTH
    INTEGER(IK), DIMENSION(0_IK : 1_IK+INT(SQRT(REAL(LENGTH/2_IK, RK)), IK)), INTENT(OUT) :: IP
    REAL(RK), DIMENSION(0_IK : LENGTH/2_IK - 1_IK), INTENT(OUT) :: WORK
    IP = 0_IK
    WORK = 0.0_RK
    CALL MAKE_FFT_TABLE_(LENGTH/2_IK, IP, WORK)
  END SUBROUTINE MAKE_FFT_DATA__

  SUBROUTINE CDFT__(LENGTH, DIRECTION, SEQUENCE, IP, WORK)
    INTEGER(IK), INTENT(IN) :: LENGTH
    INTEGER(IK), INTENT(IN) :: DIRECTION
    REAL(RK), INTENT(INOUT) :: SEQUENCE(0_IK : *)
    INTEGER(IK), INTENT(INOUT) :: IP(0_IK : *)
    REAL(RK), INTENT(INOUT) :: WORK(0_IK : *)
    IF (DIRECTION >= 0_IK) THEN
      CALL BIT_REVERSE_(LENGTH, IP(2_IK), SEQUENCE)
      CALL CFT_FORWARD_(LENGTH, SEQUENCE, WORK)
    ELSE
      CALL BIT_REVERSE_CONJUGATE_(LENGTH, IP(2_IK), SEQUENCE)
      CALL CFT_INVERSE_(LENGTH, SEQUENCE, WORK)
    END IF
  END SUBROUTINE CDFT__

  SUBROUTINE FFT_RADIX_EIGHT__(LENGTH, DIRECTION, SEQUENCE, IP, WORK)
    INTEGER(IK), INTENT(IN) :: LENGTH
    INTEGER(IK), INTENT(IN) :: DIRECTION
    REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(INOUT) :: SEQUENCE
    INTEGER(IK), DIMENSION(0_IK : 1_IK+INT(SQRT(REAL(LENGTH/2_IK, RK)), IK)), INTENT(IN) :: IP
    REAL(RK), DIMENSION(0_IK : LENGTH/2_IK - 1_IK), INTENT(IN) :: WORK
    INTEGER(IK), DIMENSION(0_IK : 1_IK+INT(SQRT(REAL(LENGTH/2_IK, RK)), IK)) :: IP_COPY
    REAL(RK), DIMENSION(0_IK : LENGTH/2_IK - 1_IK) :: WORK_COPY
    IP_COPY = IP
    WORK_COPY = WORK
    CALL CDFT__(2_IK*LENGTH, DIRECTION, SEQUENCE, IP_COPY, WORK_COPY)
  END SUBROUTINE FFT_RADIX_EIGHT__

  SUBROUTINE MAKE_FFRFT_DATA__(LENGTH, C1, S1, C2, S2)
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(LENGTH), INTENT(OUT)   :: C1
    REAL(RK), DIMENSION(LENGTH), INTENT(OUT)   :: S1
    REAL(RK), DIMENSION(LENGTH), INTENT(OUT)   :: C2
    REAL(RK), DIMENSION(LENGTH), INTENT(OUT)   :: S2
    INTEGER(IK) :: I
    REAL(RK) :: FACTOR
    REAL(RK), DIMENSION(LENGTH) :: MUL
    FACTOR = TWO_PI/REAL(LENGTH, RK)**2_IK
    MUL = FACTOR*REAL([(I, I = 0_IK, LENGTH-1_IK, 1_IK)], RK)**2_IK
    C1 = COS(MUL)
    S1 = SIN(MUL)
    MUL = -FACTOR*(REAL([(I, I = LENGTH+1_IK, 2_IK*LENGTH, 1_IK)], RK)-1.0_RK-2.0_RK*REAL(LENGTH, RK))**2_IK
    C2 = COS(MUL)
    S2 = SIN(MUL)
  END SUBROUTINE MAKE_FFRFT_DATA__

  SUBROUTINE FFRFT__(LENGTH, SEQUENCE, IP, WORK, C1, S1, C2, S2)
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(INOUT) :: SEQUENCE
    INTEGER(IK), DIMENSION(0_IK : 1_IK+INT(SQRT(REAL(LENGTH, RK)), IK)), INTENT(IN) :: IP
    REAL(RK), DIMENSION(0_IK : LENGTH - 1_IK), INTENT(IN) :: WORK
    REAL(RK), DIMENSION(LENGTH), INTENT(IN)   :: C1
    REAL(RK), DIMENSION(LENGTH), INTENT(IN)   :: S1
    REAL(RK), DIMENSION(LENGTH), INTENT(IN)   :: C2
    REAL(RK), DIMENSION(LENGTH), INTENT(IN)   :: S2
    REAL(RK), DIMENSION(4_IK*LENGTH) :: ONE, TWO, TRE
    REAL(RK), DIMENSION(2_IK*LENGTH) :: COPY
    ONE = 0.0_RK
    ONE(1_IK:2_IK*LENGTH:2_IK) = SEQUENCE(1_IK:2_IK*LENGTH:2_IK)*C1-SEQUENCE(2_IK:2_IK*LENGTH:2_IK)*S1
    ONE(2_IK:2_IK*LENGTH:2_IK) = SEQUENCE(1_IK:2_IK*LENGTH:2_IK)*S1+SEQUENCE(2_IK:2_IK*LENGTH:2_IK)*C1
    TWO = 0.0_RK
    TWO(1_IK:2_IK*LENGTH:2_IK) = +C1
    TWO(2_IK:2_IK*LENGTH:2_IK) = -S1
    TWO(2_IK*LENGTH+1_IK:4_IK*LENGTH:2_IK) = C2
    TWO(2_IK*LENGTH+2_IK:4_IK*LENGTH:2_IK) = S2
    CALL FFT_RADIX_EIGHT__(2_IK*LENGTH, FFT_FORWARD, ONE, IP, WORK)
    CALL FFT_RADIX_EIGHT__(2_IK*LENGTH, FFT_FORWARD, TWO, IP, WORK)
    TRE = ONE
    ONE(1_IK:4_IK*LENGTH:2_IK) = TRE(1_IK:4_IK*LENGTH:2_IK)*TWO(1_IK:4_IK*LENGTH:2_IK)-&
      TRE(2_IK:4_IK*LENGTH:2_IK)*TWO(2_IK:4_IK*LENGTH:2_IK)
    ONE(2_IK:4_IK*LENGTH:2_IK) = TRE(1_IK:4_IK*LENGTH:2_IK)*TWO(2_IK:4_IK*LENGTH:2_IK)+&
      TRE(2_IK:4_IK*LENGTH:2_IK)*TWO(1_IK:4_IK*LENGTH:2_IK)
    CALL FFT_RADIX_EIGHT__(2_IK*LENGTH, FFT_INVERSE, ONE, IP, WORK)
    COPY = 1.0_RK/REAL(2_IK*LENGTH, RK)*ONE(1_IK:2_IK*LENGTH:1_IK)
    SEQUENCE(1_IK:2_IK*LENGTH:2_IK) = COPY(1_IK:2_IK*LENGTH:2_IK)*C1-COPY(2_IK:2_IK*LENGTH:2_IK)*S1
    SEQUENCE(2_IK:2_IK*LENGTH:2_IK) = COPY(1_IK:2_IK*LENGTH:2_IK)*S1+COPY(2_IK:2_IK*LENGTH:2_IK)*C1
  END SUBROUTINE FFRFT__

  REAL(RK) FUNCTION FREQUENCY__(FLAG, PEAK, LENGTH, TOTAL, WINDOW, SEQUENCE) &
    BIND(C, NAME = "frequency__")
    INTEGER(IK), INTENT(IN):: FLAG
    INTEGER(IK), INTENT(IN) :: PEAK
    INTEGER(IK), INTENT(IN):: LENGTH
    REAL(RK), INTENT(IN) :: TOTAL
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
    REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    REAL(RK), DIMENSION(2_IK*LENGTH) :: LOCAL, FOURIER
    INTEGER(IK) :: FST, CND
    REAL(RK) :: FACTOR
    REAL(RK), DIMENSION(LENGTH) :: MUL, COS_MUL, SIN_MUL
    INTEGER(IK) :: I
    LOCAL(1_IK:2_IK*LENGTH:2_IK) = (SEQUENCE(1_IK:2_IK*LENGTH:2_IK)-SUM(WINDOW*SEQUENCE(1_IK:2_IK*LENGTH:2_IK))/TOTAL)*WINDOW
    LOCAL(2_IK:2_IK*LENGTH:2_IK) = (SEQUENCE(2_IK:2_IK*LENGTH:2_IK)-SUM(WINDOW*SEQUENCE(2_IK:2_IK*LENGTH:2_IK))/TOTAL)*WINDOW
    FOURIER = LOCAL
    CALL FFT_RADIX_EIGHT__(LENGTH, FFT_FORWARD, FOURIER, BANK%BIT_FFT, BANK%TRIG_FFT)
    MUL = LOG10(SQRT(FOURIER(1_IK:2_IK*LENGTH:2_IK)**2_IK+FOURIER(2_IK:2_IK*LENGTH:2_IK)**2_IK)+1.E-16_RK)
    IF (PEAK == 0_IK) THEN
      FST = INT(__MAXLOC__(MUL(1_IK:LENGTH/(2_IK-FLAG):1_IK), 1_IK), IK)
    ELSE
      FST = PEAK_(LENGTH/(2_IK-FLAG), MUL(1_IK:LENGTH/(2_IK-FLAG):1_IK), ABS(PEAK))
    END IF
    IF (PEAK < 0_IK) THEN
      FREQUENCY__ = REAL(FST-1_IK, RK)/REAL(LENGTH, RK)
      RETURN
    END IF
    FACTOR = TWO_PI*REAL(FST-2_IK, RK)/REAL(LENGTH, RK)
    MUL = FACTOR*REAL([(I, I = 0_IK, LENGTH-1_IK, 1_IK)], RK)
    COS_MUL = COS(MUL)
    SIN_MUL = SIN(MUL)
    FOURIER = LOCAL
    LOCAL(1_IK:2_IK*LENGTH:2_IK) = FOURIER(1_IK:2_IK*LENGTH:2_IK)*COS_MUL-FOURIER(2_IK:2_IK*LENGTH:2_IK)*SIN_MUL
    LOCAL(2_IK:2_IK*LENGTH:2_IK) = FOURIER(1_IK:2_IK*LENGTH:2_IK)*SIN_MUL+FOURIER(2_IK:2_IK*LENGTH:2_IK)*COS_MUL
    CALL FFRFT__(LENGTH, LOCAL, BANK%BIT_FFRFT, BANK%TRIG_FFRFT, BANK%COS_FST, BANK%SIN_FST, BANK%COS_LST, BANK%SIN_LST)
    MUL = LOG10(SQRT(LOCAL(1_IK:2_IK*LENGTH:2_IK)**2_IK+LOCAL(2_IK:2_IK*LENGTH:2_IK)**2_IK)+1.E-16_RK)
    CND = INT(__MAXLOC__(MUL, 1), IK)
    FREQUENCY__ = REAL(CND, RK)-0.5_RK+(MUL(-1_IK+CND)-MUL(CND))/(MUL(-1_IK+CND)-2.0_RK*MUL(CND)+MUL(1_IK+CND))
    FREQUENCY__ = (REAL(FST, RK)-2.0_RK+2.0_RK*(FREQUENCY__-1.0_RK)/REAL(LENGTH, RK))/REAL(LENGTH, RK)
  END FUNCTION FREQUENCY__

  SUBROUTINE DECOMPOSITION__(FLAG, METHOD, LENGTH, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP) &
    BIND(C, NAME = "decomposition__")
    INTEGER(IK), INTENT(IN):: FLAG
    INTEGER(IK), INTENT(IN):: METHOD
    INTEGER(IK), INTENT(IN):: LENGTH
    REAL(RK), INTENT(IN) :: TOTAL
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
    REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    INTEGER(IK), INTENT(IN) :: LOOP
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: FREQUENCY
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: COS_AMP
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: SIN_AMP
    INTEGER(IK) :: I
    REAL(RK), DIMENSION(2_IK*LENGTH) :: LOCAL
    REAL(RK), DIMENSION(LENGTH) :: LIST
    REAL(RK), DIMENSION(2_IK*LENGTH) :: DELTA
    LOCAL = SEQUENCE-TOTAL
    LIST = TWO_PI*REAL([(I, I=1_IK, LENGTH, 1_IK)], RK)
    DO I = 1_IK, LOOP, 1_IK
      IF (METHOD == DECOMPOSITION_MAX) THEN
        FREQUENCY(I) = FREQUENCY__(FLAG, 0_IK, LENGTH, TOTAL, WINDOW, LOCAL)
      ELSE IF(METHOD == DECOMPOSITION_PEAK) THEN
        FREQUENCY(I) = FREQUENCY__(FLAG, +I, LENGTH, TOTAL, WINDOW, SEQUENCE)
      ELSE IF(METHOD == DECOMPOSITION_PEAK_FFT) THEN
        FREQUENCY(I) = FREQUENCY__(FLAG, -I, LENGTH, TOTAL, WINDOW, SEQUENCE)
      END IF
      DELTA(1_IK:2_IK*LENGTH:2_IK) = +COS(FREQUENCY(I)*LIST)
      DELTA(2_IK:2_IK*LENGTH:2_IK) = -SIN(FREQUENCY(I)*LIST)
      COS_AMP(I) = SUM(SEQUENCE(1_IK:2_IK*LENGTH:2_IK)*DELTA(1_IK:2_IK*LENGTH:2_IK)*WINDOW)+&
        SUM(SEQUENCE(2_IK:2_IK*LENGTH:2_IK)*DELTA(2_IK:2_IK*LENGTH:2_IK)*WINDOW)
      SIN_AMP(I) = SUM(SEQUENCE(2_IK:2_IK*LENGTH:2_IK)*DELTA(1_IK:2_IK*LENGTH:2_IK)*WINDOW)-&
        SUM(SEQUENCE(1_IK:2_IK*LENGTH:2_IK)*DELTA(2_IK:2_IK*LENGTH:2_IK)*WINDOW)
      COS_AMP(I) = 1.0_RK/REAL(LENGTH, RK)*COS_AMP(I)
      SIN_AMP(I) = 1.0_RK/REAL(LENGTH, RK)*SIN_AMP(I)
      LOCAL(1_IK:2_IK*LENGTH:2_IK) = LOCAL(1_IK:2_IK*LENGTH:2_IK)-&
        COS_AMP(I)*DELTA(1_IK:2_IK*LENGTH:2_IK)+SIN_AMP(I)*DELTA(2_IK:2_IK*LENGTH:2_IK)
      LOCAL(2_IK:2_IK*LENGTH:2_IK) = LOCAL(2_IK:2_IK*LENGTH:2_IK)-&
        COS_AMP(I)*DELTA(2_IK:2_IK*LENGTH:2_IK)-SIN_AMP(I)*DELTA(1_IK:2_IK*LENGTH:2_IK)
    END DO
    IF (FLAG == 0_IK) THEN
      COS_AMP = 2.0_RK*COS_AMP
      SIN_AMP = 2.0_RK*SIN_AMP
    END IF
  END SUBROUTINE DECOMPOSITION__

  SUBROUTINE FREQUENCY_LIST__(FLAG, METHOD, LENGTH, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY) &
    BIND(C, NAME = "frequency_list__")
    INTEGER(IK), INTENT(IN):: FLAG
    INTEGER(IK), INTENT(IN) :: METHOD
    INTEGER(IK), INTENT(IN):: LENGTH
    REAL(RK), INTENT(IN) :: TOTAL
    REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
    REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    INTEGER(IK), INTENT(IN) :: LOOP
    REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: FREQUENCY
    REAL(RK), DIMENSION(LOOP) :: COS_AMP
    REAL(RK), DIMENSION(LOOP) :: SIN_AMP
    CALL DECOMPOSITION__(FLAG, METHOD, LENGTH, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP)
  END SUBROUTINE FREQUENCY_LIST__

END MODULE SIGNAL

! ! BINARY SEARCH MAXIMIZATION
! REAL(RK) FUNCTION BINARY_(FUN, GUESS, INTERVAL, LIMIT, TOLERANCE)
! INTERFACE
!   REAL(RK) FUNCTION FUN(ARG)
!     IMPORT :: RK
!     REAL(RK), INTENT(IN) :: ARG
!   END FUNCTION FUN
! END INTERFACE
! REAL(RK), INTENT(IN) :: GUESS
! REAL(RK), INTENT(IN) :: INTERVAL
! INTEGER(IK), INTENT(IN) :: LIMIT
! REAL(RK), INTENT(IN) :: TOLERANCE
! REAL(RK) :: DELTA
! REAL(RK) :: XL, XR, XX
! REAL(RK) :: FL, FR
! INTEGER(IK) :: I
! DELTA = INTERVAL/2.0_RK
! XL = GUESS-DELTA
! XR = GUESS+DELTA
! FL = FUN(FL)
! FR = FUN(FR)
! DO I = 1_IK, LIMIT, 1_IK
!   IF (FL > FR) THEN
!     XX = XL
!   ELSE
!     XX = XR
!   END IF
!   XL = XX-DELTA
!   XR = XX+DELTA
!   DELTA = DELTA/2.0_RK
!   FL = FUN(XL)
!   FR = FUN(XR)
!   IF (ABS(FL-FR) < TOLERANCE) EXIT
! END DO
! IF (FL > FR) THEN
!   BINARY_ = XL
! ELSE
!   BINARY_ = XR
! END IF
! END FUNCTION BINARY_

! ! GOLDEN SEARCH MAXIMIZATION
! REAL(RK) FUNCTION GOLDEN_(FUN, GUESS, INTERVAL, LIMIT, TOLERANCE)
! INTERFACE
!   REAL(RK) FUNCTION FUN(ARG)
!     IMPORT :: RK
!     REAL(RK), INTENT(IN) :: ARG
!   END FUNCTION FUN
! END INTERFACE
! REAL(RK), INTENT(IN) :: GUESS
! REAL(RK), INTENT(IN) :: INTERVAL
! INTEGER(IK), INTENT(IN) :: LIMIT
! REAL(RK), INTENT(IN) :: TOLERANCE
! REAL(RK), PARAMETER :: GOLDEN = (1.0_RK+SQRT(5.0_RK))/2.0_RK
! REAL(RK), PARAMETER :: PSI = 1.0_RK-1.0_RK/GOLDEN
! REAL(RK), PARAMETER :: PHI = 1.0_RK/GOLDEN
! REAL(RK) :: DELTA
! REAL(RK) :: XL, XR
! REAL(RK) :: FL, FR
! REAL(RK) :: M1, M2
! INTEGER(IK) :: I
! DELTA = INTERVAL/2.0_RK
! XL = GUESS-DELTA
! XR = GUESS+DELTA
! M1 = XL+PSI*(XR-XL)
! M2 = XL+PHI*(XR-XL)
! FL = FUN(M1)
! FR = FUN(M2)
! DO I = 1_IK, LIMIT, 1_IK
!   IF (I > LIMIT) EXIT
!   IF (ABS(FL-FR) < TOLERANCE) EXIT
!   IF (FL > FR) THEN
!     XR = M2
!     M2 = M1
!     M1 = PSI*XL+PHI*M2
!     FR = FL
!     FL = FUN(M1)
!   ELSE
!     XL = M1
!     M1 = M2
!     M2 = PHI*M1+PSI*XR
!     FL = FR
!     FR = FUN(M2)
!   END IF
! END DO
! IF (FL > FR) THEN
!   GOLDEN_ = XL
! ELSE
!   GOLDEN_ = XR
! END IF
! END FUNCTION GOLDEN_

! INTERFACE GAMMA_
! REAL(RK) FUNCTION GAMMA_(A) &
!   BIND(C, NAME = "gsl_sf_gamma")
!   USE, INTRINSIC :: ISO_C_BINDING,   ONLY: RK => C_DOUBLE
!   REAL(RK), VALUE :: A
!   END FUNCTION GAMMA_
! END INTERFACE GAMMA_

! INTERFACE GAMMA_INCOMPLETE_
! REAL(RK) FUNCTION GAMMA_INCOMPLETE_(A, X) &
!   BIND(C, NAME = "gsl_sf_gamma_inc")
!   USE, INTRINSIC :: ISO_C_BINDING,   ONLY: RK => C_DOUBLE
!   REAL(RK), VALUE :: A
!   REAL(RK), VALUE :: X
! END FUNCTION GAMMA_INCOMPLETE_
! END INTERFACE GAMMA_INCOMPLETE_

! INTERFACE GAMMA_
! PROCEDURE GAMMA_INCOMPLETE_
! PROCEDURE GAMMA_REGULARIZED_
! END INTERFACE GAMMA_

! ! GAMMA REGULARIZED
! ! (FUNCTION) GAMMA_REGULARIZED_(<A>, <X>, <Y>)
! ! <A>                    -- (IN)     A (RK)
! ! <X>                    -- (IN)     X (RK)
! ! <Y>                    -- (IN)     Y (RK)
! ! <FACTORIAL_>           -- (OUT)    GAMMA REGULARIZED
! REAL(RK) FUNCTION GAMMA_REGULARIZED_(A, X, Y)
! REAL(RK), INTENT(IN) :: A
! REAL(RK), INTENT(IN) :: X
! REAL(RK), INTENT(IN) :: Y
! GAMMA_REGULARIZED_ = (GAMMA_(A, X)-GAMMA_(A, Y))/GAMMA_(A)
! END FUNCTION GAMMA_REGULARIZED_
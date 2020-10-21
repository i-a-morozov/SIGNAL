
#include "signal.inc"

SUBMODULE (SIGNAL) TRANSFORMATION
  IMPLICIT NONE
  CONTAINS
  ! ############################################################################################################################# !
  ! (LINEAR) FRACTIONAL COMPLEX DISCRETE FOURIER TRANSFORM (POWER OF TWO INPUT LENGTH INPUT)
  ! (SUBROUTINE) FFRFT_(<LENGTH>, <ARGUMENT>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <ARGUMENT>             -- (IN)     PARAMETER (RK)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <SEQUENCE>             -- (OUT)    CDFT (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., FR_I, FI_I, ...]
  ! void    ffrft_(int*, double*, double*) ;
  MODULE SUBROUTINE FFRFT_(LENGTH, ARGUMENT, SEQUENCE) &
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
  ! ############################################################################################################################# !
  ! (LINEAR) FRACTIONAL COMPLEX DISCRETE FOURIER TRANSFORM (POWER OF TWO INPUT LENGTH INPUT)
  ! (SUBROUTINE) FFRFT__(<LENGTH>, <ARGUMENT>, <SEQUENCE>, <IP>, <WORK>, <COS_FST>, <SIN_FST>, <COS_LST>, <SIN_LST>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <ARGUMENT>             -- (IN)     PARAMETER (RK)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <IP>                   -- (IN)     FFRFT BIT DATA
  ! <WORK>                 -- (IN)     FFRFT TRIG DATA
  ! <COS_FST>              -- (IN)     FIRST COS ARRAY
  ! <SIN_FST>              -- (IN)     FIRST SIN ARRAY
  ! <COS_LST>              -- (IN)     LAST COS ARRAY
  ! <SIN_LAT>              -- (IN)     LAST SIN ARRAY
  ! <SEQUENCE>             -- (OUT)    CDFT (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., FR_I, FI_I, ...]
  MODULE SUBROUTINE FFRFT__(LENGTH, SEQUENCE, IP, WORK, COS_FST, SIN_FST, COS_LST, SIN_LST)
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(INOUT) :: SEQUENCE
    INTEGER(IK), DIMENSION(0_IK : 1_IK+INT(SQRT(REAL(LENGTH, RK)), IK)), INTENT(IN) :: IP
    REAL(RK), DIMENSION(0_IK : LENGTH - 1_IK), INTENT(IN) :: WORK
    REAL(RK), DIMENSION(LENGTH), INTENT(IN)   :: COS_FST
    REAL(RK), DIMENSION(LENGTH), INTENT(IN)   :: SIN_FST
    REAL(RK), DIMENSION(LENGTH), INTENT(IN)   :: COS_LST
    REAL(RK), DIMENSION(LENGTH), INTENT(IN)   :: SIN_LST
    REAL(RK), DIMENSION(4_IK*LENGTH) :: ONE, TWO, TRE
    REAL(RK), DIMENSION(2_IK*LENGTH) :: COPY
    ONE = 0.0_RK
    ONE(1_IK:2_IK*LENGTH:2_IK) = SEQUENCE(1_IK:2_IK*LENGTH:2_IK)*COS_FST-SEQUENCE(2_IK:2_IK*LENGTH:2_IK)*SIN_FST
    ONE(2_IK:2_IK*LENGTH:2_IK) = SEQUENCE(1_IK:2_IK*LENGTH:2_IK)*SIN_FST+SEQUENCE(2_IK:2_IK*LENGTH:2_IK)*COS_FST
    TWO = 0.0_RK
    TWO(1_IK:2_IK*LENGTH:2_IK) = +COS_FST
    TWO(2_IK:2_IK*LENGTH:2_IK) = -SIN_FST
    TWO(2_IK*LENGTH+1_IK:4_IK*LENGTH:2_IK) = COS_LST
    TWO(2_IK*LENGTH+2_IK:4_IK*LENGTH:2_IK) = SIN_LST
    CALL FFT_RADIX_EIGHT__(2_IK*LENGTH, FFT_FORWARD, ONE, IP, WORK)
    CALL FFT_RADIX_EIGHT__(2_IK*LENGTH, FFT_FORWARD, TWO, IP, WORK)
    TRE = ONE
    ONE(1_IK:4_IK*LENGTH:2_IK) = TRE(1_IK:4_IK*LENGTH:2_IK)*TWO(1_IK:4_IK*LENGTH:2_IK)-&
      TRE(2_IK:4_IK*LENGTH:2_IK)*TWO(2_IK:4_IK*LENGTH:2_IK)
    ONE(2_IK:4_IK*LENGTH:2_IK) = TRE(1_IK:4_IK*LENGTH:2_IK)*TWO(2_IK:4_IK*LENGTH:2_IK)+&
      TRE(2_IK:4_IK*LENGTH:2_IK)*TWO(1_IK:4_IK*LENGTH:2_IK)
    CALL FFT_RADIX_EIGHT__(2_IK*LENGTH, FFT_INVERSE, ONE, IP, WORK)
    COPY = 1.0_RK/REAL(2_IK*LENGTH, RK)*ONE(1_IK:2_IK*LENGTH:1_IK)
    SEQUENCE(1_IK:2_IK*LENGTH:2_IK) = COPY(1_IK:2_IK*LENGTH:2_IK)*COS_FST-COPY(2_IK:2_IK*LENGTH:2_IK)*SIN_FST
    SEQUENCE(2_IK:2_IK*LENGTH:2_IK) = COPY(1_IK:2_IK*LENGTH:2_IK)*SIN_FST+COPY(2_IK:2_IK*LENGTH:2_IK)*COS_FST
  END SUBROUTINE FFRFT__
  ! ############################################################################################################################# !
  ! (FFTW) COMPLEX DISCRETE FOURIER TRANSFORM (POWER OF TWO INPUT LENGTH INPUT)
  ! (SUBROUTINE) FFT_EXTERNAL_(<LENGTH>, <DIRECTION>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <DIRECTION>            -- (IN)     DIRECTION (IK), FFT_FORWARD = +1_IK OR FFT_INVERSE = -1_IK
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <SEQUENCE>             -- (OUT)    CDFT (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., FR_I, FI_I, ...]
  ! void    fft_external_(int*, int*, double*) ;
  MODULE SUBROUTINE FFT_EXTERNAL_(LENGTH, DIRECTION, SEQUENCE) &
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
  ! ############################################################################################################################# !
  ! (NRF77) COMPLEX DISCRETE FOURIER TRANSFORM (POWER OF TWO INPUT LENGTH INPUT)
  ! (SUBROUTINE) FFT_RADIX_TWO_(<LENGTH>, <DIRECTION>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <DIRECTION>            -- (IN)     DIRECTION (IK), FFT_FORWARD = +1_IK OR FFT_INVERSE = -1_IK
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <SEQUENCE>             -- (OUT)    CDFT (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., FR_I, FI_I, ...]
  ! void    fft_radix_two_(int*, int*, double*) ;
  MODULE SUBROUTINE FFT_RADIX_TWO_(LENGTH, DIRECTION, SEQUENCE) &
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
  ! ############################################################################################################################# !
  ! COMPLEX DISCRETE FOURIER TRANSFORM (POWER OF TWO INPUT LENGTH INPUT)
  ! (SUBROUTINE) FFT_RADIX_EIGHT_(<LENGTH>, <DIRECTION>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <DIRECTION>            -- (IN)     DIRECTION (IK), FFT_FORWARD = +1_IK OR FFT_INVERSE = -1_IK
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <SEQUENCE>             -- (OUT)    CDFT (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., FR_I, FI_I, ...]
  ! void    fft_radix_eight_(int*, int*, double*) ;
  MODULE SUBROUTINE FFT_RADIX_EIGHT_(LENGTH, DIRECTION, SEQUENCE) &
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
  ! ############################################################################################################################# !
  ! (TAKUYA OOURA) COMPLEX DISCRETE FOURIER TRANSFORM (POWER OF TWO INPUT LENGTH INPUT)
  ! (SUBROUTINE) FFT_RADIX_EIGHT__(<LENGTH>, <DIRECTION>, <SEQUENCE>, <IP>, <WORK>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <DIRECTION>            -- (IN)     DIRECTION (IK), FFT_FORWARD = +1_IK OR FFT_INVERSE = -1_IK
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <SEQUENCE>             -- (OUT)    CDFT (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., FR_I, FI_I, ...]
  ! <IP>                   -- (IN)     FFRFT BIT DATA
  ! <WORK>                 -- (IN)     FFRFT TRIG DATA
  MODULE SUBROUTINE FFT_RADIX_EIGHT__(LENGTH, DIRECTION, SEQUENCE, IP, WORK)
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
  ! ############################################################################################################################# !
  ! COMPUTE DATA TABLE
  ! (SUBROUTINE) COMPUTE_TABLE_(<LENGTH>, <PAD>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <PAD>                  -- (IN)     PADDED LENGTH (IK)
  ! void    compute_table_(int*, int*) ;
  MODULE SUBROUTINE COMPUTE_TABLE_(LENGTH, PAD) &
    BIND(C, NAME = "compute_table_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    INTEGER(IK), INTENT(IN) :: PAD
    ALLOCATE(BANK%BIT_FFT(2_IK+INT(SQRT(REAL(PAD/2_IK, RK)), IK)))
    ALLOCATE(BANK%BIT_FFRFT(2_IK+INT(SQRT(REAL(LENGTH/1_IK, RK)), IK)))
    ALLOCATE(BANK%TRIG_FFT(PAD/2_IK))
    ALLOCATE(BANK%TRIG_FFRFT(LENGTH/1_IK))
    ALLOCATE(BANK%COS_FST(LENGTH))
    ALLOCATE(BANK%SIN_FST(LENGTH))
    ALLOCATE(BANK%COS_LST(LENGTH))
    ALLOCATE(BANK%SIN_LST(LENGTH))
    CALL MAKE_FFT_DATA__(PAD, BANK%BIT_FFT, BANK%TRIG_FFT)
    CALL MAKE_FFT_DATA__(2_IK*LENGTH, BANK%BIT_FFRFT, BANK%TRIG_FFRFT)
    CALL MAKE_FFRFT_DATA__(LENGTH, BANK%COS_FST, BANK%SIN_FST, BANK%COS_LST, BANK%SIN_LST)
  END SUBROUTINE COMPUTE_TABLE_
  ! ############################################################################################################################# !
  ! DESTROY DATA TABLE
  ! (SUBROUTINE) DESTROY_TABLE_()
  ! void    destroy_table_(int*) ;
  MODULE SUBROUTINE DESTROY_TABLE_() &
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
  ! ############################################################################################################################# !
  ! MAKE FFT DATA
  ! (SUBROUTINE) MAKE_FFT_DATA__(<LENGTH>, <IP>, <WORK>)
  MODULE SUBROUTINE MAKE_FFT_DATA__(LENGTH, IP, WORK)
    INTEGER(IK), INTENT(IN) :: LENGTH
    INTEGER(IK), DIMENSION(0_IK : 1_IK+INT(SQRT(REAL(LENGTH/2_IK, RK)), IK)), INTENT(OUT) :: IP
    REAL(RK), DIMENSION(0_IK : LENGTH/2_IK - 1_IK), INTENT(OUT) :: WORK
    IP = 0_IK
    WORK = 0.0_RK
    CALL MAKE_FFT_TABLE_(LENGTH/2_IK, IP, WORK)
  END SUBROUTINE MAKE_FFT_DATA__
  ! ############################################################################################################################# !
  ! MAKE FFRFT DATA
  ! (SUBROUTINE) MAKE_FFRFT_DATA__(LENGTH, COS_FST, SIN_FST, COS_LST, SIN_LST)
  SUBROUTINE MAKE_FFRFT_DATA__(LENGTH, COS_FST, SIN_FST, COS_LST, SIN_LST)
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(LENGTH), INTENT(OUT)   :: COS_FST
    REAL(RK), DIMENSION(LENGTH), INTENT(OUT)   :: SIN_FST
    REAL(RK), DIMENSION(LENGTH), INTENT(OUT)   :: COS_LST
    REAL(RK), DIMENSION(LENGTH), INTENT(OUT)   :: SIN_LST
    INTEGER(IK) :: I
    REAL(RK) :: FACTOR
    REAL(RK), DIMENSION(LENGTH) :: MUL
    FACTOR = TWO_PI/REAL(LENGTH, RK)**2_IK
    MUL = FACTOR*REAL([(I, I = 0_IK, LENGTH-1_IK, 1_IK)], RK)**2_IK
    COS_FST = COS(MUL)
    SIN_FST = SIN(MUL)
    MUL = -FACTOR*(REAL([(I, I = LENGTH+1_IK, 2_IK*LENGTH, 1_IK)], RK)-1.0_RK-2.0_RK*REAL(LENGTH, RK))**2_IK
    COS_LST = COS(MUL)
    SIN_LST = SIN(MUL)
  END SUBROUTINE MAKE_FFRFT_DATA__
  ! ############################################################################################################################# !
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
  ! ############################################################################################################################# !
  ! (TAKUYA OOURA) CDFT__
  MODULE SUBROUTINE CDFT__(LENGTH, DIRECTION, SEQUENCE, IP, WORK)
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
  ! ############################################################################################################################# !
  ! BIT_REVERSE_
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
  ! ############################################################################################################################# !
  ! BIT_REVERSE_CONJUGATE_
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
  ! ############################################################################################################################# !
  ! MAKE_FFT_TABLE_
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
  ! ############################################################################################################################# !
  ! CFT_FORWARD_
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
  ! ############################################################################################################################# !
  ! CFT_INVERSE_
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
  ! ############################################################################################################################# !
  ! CFT_1ST_
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
  ! ############################################################################################################################# !
  ! CFT_MDL_
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
  ! ############################################################################################################################# !
END SUBMODULE TRANSFORMATION
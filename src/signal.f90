! SIGNAL, 2018-2020, I.A.MOROZOV@INP.NSK.SU
MODULE SIGNAL
  USE, INTRINSIC :: ISO_C_BINDING,   ONLY: IK => C_INT, RK => C_DOUBLE, C_SIZEOF
  IMPLICIT NONE
  PRIVATE
  ! ############################################################################################################################# !
  ! GLOBAL
  ! ############################################################################################################################# !
  PUBLIC :: IK
  PUBLIC :: RK
  INTEGER(IK), PUBLIC, PARAMETER :: IK_SIZE                = C_SIZEOF(IK)
  INTEGER(IK), PUBLIC, PARAMETER :: RK_SIZE                = C_SIZEOF(RK)
  REAL(RK),    PUBLIC, PARAMETER :: ONE_PI                 = 2.0_RK*ACOS(0.0_RK)
  REAL(RK),    PUBLIC, PARAMETER :: TWO_PI                 = 2.0_RK*ONE_PI
  ! ############################################################################################################################# !
  ! EXTERNAL
  ! ############################################################################################################################# !
  EXTERNAL :: DGEMV               ! (BLAS)
  EXTERNAL :: DGEMM               ! (BLAS)
  EXTERNAL :: DGESVD              ! (LAPACK)
  EXTERNAL :: DSAUPD              ! (ARPACK)
  EXTERNAL :: DSEUPD              ! (ARPACK)
  EXTERNAL :: DFFTW_PLAN_DFT_1D   ! (FFTW)
  EXTERNAL :: DFFTW_EXECUTE_DFT   ! (FFTW)
  EXTERNAL :: DFFTW_DESTROY_PLAN  ! (FFTW)
  ! ############################################################################################################################# !
  ! AUXILIARY
  ! ############################################################################################################################# !
  ! FACTORIAL
  ! (FUNCTION) FACTORIAL_(<NUMBER>)
  ! <NUMBER>               -- (IN)     NUMBER (IK)
  ! <FACTORIAL_>           -- (OUT)    FACTORIAL OF <N> (RK REAL)
  INTERFACE
    MODULE REAL(RK) FUNCTION FACTORIAL_(NUMBER)
      INTEGER(IK), INTENT(IN) :: NUMBER
    END FUNCTION FACTORIAL_
  END INTERFACE
  PUBLIC :: FACTORIAL_
  ! ############################################################################################################################# !
  ! GAMMA (GSL)
  ! (FUNCTION) GAMMA_(<NUMBER>)
  ! <NUMBER>               -- (IN)     NUMBER (RK)
  ! <GAMMA_>               -- (OUT)    GAMMA OF <N> (RK REAL)
  INTERFACE GAMMA_
    REAL(RK) FUNCTION GAMMA_(NUMBER) &
      BIND(C, NAME = "gsl_sf_gamma")
      IMPORT :: RK
      REAL(RK), VALUE :: NUMBER
      END FUNCTION GAMMA_
  END INTERFACE GAMMA_
  ! ############################################################################################################################# !
  ! GAMMA INCOMPLETE (GSL)
  ! (FUNCTION) GAMMA_INCOMPLETE_(<NUMBER>)
  ! <A>                    -- (IN)     A (RK)
  ! <X>                    -- (IN)     X (RK)
  ! <GAMMA_INCOMPLETE_>    -- (OUT)    GAMMA INCOMPLETE OF <A> AND <X> (RK REAL)
  INTERFACE GAMMA_INCOMPLETE_
    REAL(RK) FUNCTION GAMMA_INCOMPLETE_(A, X) &
      BIND(C, NAME = "gsl_sf_gamma_inc")
      IMPORT :: RK
      REAL(RK), VALUE :: A
      REAL(RK), VALUE :: X
    END FUNCTION GAMMA_INCOMPLETE_
  END INTERFACE GAMMA_INCOMPLETE_
  ! ############################################################################################################################# !
  ! GAMMA REGULARIZED
  ! (FUNCTION) GAMMA_REGULARIZED_(<NUMBER>)
  ! <A>                    -- (IN)     A (RK)
  ! <X>                    -- (IN)     X (RK)
  ! <Y>                    -- (IN)     Y (RK)
  ! <GAMMA_REGULARIZED_>   -- (OUT)    GAMMA REGULARIZED OF <A>, <X> AND <Y> (RK REAL)
  INTERFACE
    MODULE REAL(RK) FUNCTION GAMMA_REGULARIZED_(A, X, Y)
      REAL(RK), INTENT(IN) :: A
      REAL(RK), INTENT(IN) :: X
      REAL(RK), INTENT(IN) :: Y
    END FUNCTION GAMMA_REGULARIZED_
  END INTERFACE
  ! ############################################################################################################################# !
  ! GENERIC GAMMA
  INTERFACE GAMMA_
    PROCEDURE GAMMA_INCOMPLETE_
    PROCEDURE GAMMA_REGULARIZED_
  END INTERFACE GAMMA_
  PUBLIC :: GAMMA_
  ! ############################################################################################################################# !
  ! MODIFIED BESSEL I_0(X) (GSL)
  ! (FUNCTION) BESSEL_(<NUMBER>)
  ! <NUMBER>               -- (IN)     NUMBER (RK)
  ! <BESSEL_>              -- (OUT)    BESSEL I_0(<X>) (RK REAL)
  INTERFACE BESSEL_
    REAL(RK) FUNCTION BESSEL_(NUMBER) &
      BIND(C, NAME = "gsl_sf_bessel_I0")
      IMPORT :: RK
      REAL(RK), VALUE :: NUMBER
      END FUNCTION BESSEL_
  END INTERFACE BESSEL_
  PUBLIC :: BESSEL_
  ! ############################################################################################################################# !
  ! MINLOC
  ! (FUNCTION) MINLOC_(<SEQUENCE>)
  ! <SEQUENCE>             -- (IN)     SEQUENCE (RK ARRAY)
  ! <MINLOC_>              -- (OUT)    MINIMUM LOCATION (IK)
  INTERFACE
    MODULE INTEGER(IK) FUNCTION MINLOC_(SEQUENCE, EMPTY)
      REAL(RK), DIMENSION(:), CONTIGUOUS, INTENT(IN) :: SEQUENCE
      INTEGER, INTENT(IN) :: EMPTY
    END FUNCTION MINLOC_
  END INTERFACE
  ! ############################################################################################################################# !
  ! MAXLOC
  ! (FUNCTION) MAXLOC_(<SEQUENCE>)
  ! <SEQUENCE>             -- (IN)     SEQUENCE (RK ARRAY)
  ! <MAXLOC_>              -- (OUT)    MAXIMUM LOCATION (IK)
  INTERFACE
    MODULE INTEGER(IK) FUNCTION MAXLOC_(SEQUENCE, EMPTY)
      REAL(RK), DIMENSION(:), CONTIGUOUS, INTENT(IN) :: SEQUENCE
      INTEGER, INTENT(IN) :: EMPTY
    END FUNCTION MAXLOC_
  END  INTERFACE
  ! ############################################################################################################################# !
  ! SORT (BUBBLE, DESCENDING)
  ! (SUBROUTINE) SORT_BUBBLE_(<LENGTH>, <SEQUENCE>, <FST>, <LST>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (IN)     (UNSORTED) SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (OUT)    (SORTED, DESCENDING) SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  INTERFACE
    MODULE SUBROUTINE SORT_BUBBLE_(LENGTH, SEQUENCE, FST, LST)
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), DIMENSION(LENGTH), INTENT(INOUT) :: SEQUENCE
      INTEGER(IK), INTENT(IN) :: FST
      INTEGER(IK), INTENT(IN) :: LST
    END SUBROUTINE SORT_BUBBLE_
  END INTERFACE
  ! ############################################################################################################################# !
  ! SORT (QUICK, DESCENDING)
  ! (SUBROUTINE) SORT_QUICK_(<LENGTH>, <SEQUENCE>, <FST>, <LST>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (IN)     (UNSORTED) SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (OUT)    (SORTED, DESCENDING) SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  INTERFACE
    MODULE RECURSIVE SUBROUTINE SORT_QUICK_(LENGTH, SEQUENCE, FST, LST)
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), DIMENSION(:), INTENT(INOUT) :: SEQUENCE
      INTEGER(IK), INTENT(IN) :: FST
      INTEGER(IK), INTENT(IN) :: LST
    END SUBROUTINE SORT_QUICK_
  END INTERFACE
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
  INTERFACE
    MODULE SUBROUTINE GENERATE_SIGNAL_(FLAG, LENGTH, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP) &
      BIND(C, NAME = "generate_signal_")
      INTEGER(IK), INTENT(IN) :: FLAG
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(OUT) :: SEQUENCE
      INTEGER(IK), INTENT(IN) :: LOOP
      REAL(RK), DIMENSION(LOOP), INTENT(IN) :: FREQUENCY
      REAL(RK), DIMENSION(LOOP), INTENT(IN) :: COS_AMP
      REAL(RK), DIMENSION(LOOP), INTENT(IN) :: SIN_AMP
    END SUBROUTINE GENERATE_SIGNAL_
  END INTERFACE
  PUBLIC :: GENERATE_SIGNAL_
  ! ############################################################################################################################# !
  ! TRANSFORMATION
  ! ############################################################################################################################# !
  INTEGER(IK), PUBLIC, PARAMETER :: FFT_FORWARD            = +1_IK               ! FORWARD FFT
  INTEGER(IK), PUBLIC, PARAMETER :: FFT_INVERSE            = -1_IK               ! INVERSE FFT
  ! ############################################################################################################################# !
  ! FFT/FFRFT DATAMEMORIZATION TYPE
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
  ! ############################################################################################################################# !
  ! FFT/FFRFT DATAMEMORIZATION CONTANER
  TYPE(TABLE), PROTECTED :: BANK
  ! ############################################################################################################################# !
  ! (LINEAR) FRACTIONAL COMPLEX DISCRETE FOURIER TRANSFORM (POWER OF TWO INPUT LENGTH INPUT)
  ! (SUBROUTINE) FFRFT_(<LENGTH>, <ARGUMENT>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <ARGUMENT>             -- (IN)     PARAMETER (RK)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <SEQUENCE>             -- (OUT)    CDFT (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., FR_I, FI_I, ...]
  ! void    ffrft_(int*, double*, double*) ;
  INTERFACE
    MODULE SUBROUTINE FFRFT_(LENGTH, ARGUMENT, SEQUENCE) &
      BIND(C, NAME = "ffrft_")
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), INTENT(IN) :: ARGUMENT
      REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(INOUT) :: SEQUENCE
    END SUBROUTINE
  END INTERFACE
  PUBLIC :: FFRFT_
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
  INTERFACE
  MODULE SUBROUTINE FFRFT__(LENGTH, SEQUENCE, IP, WORK, COS_FST, SIN_FST, COS_LST, SIN_LST)
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(INOUT) :: SEQUENCE
    INTEGER(IK), DIMENSION(0_IK : 1_IK+INT(SQRT(REAL(LENGTH, RK)), IK)), INTENT(IN) :: IP
    REAL(RK), DIMENSION(0_IK : LENGTH - 1_IK), INTENT(IN) :: WORK
    REAL(RK), DIMENSION(LENGTH), INTENT(IN)   :: COS_FST
    REAL(RK), DIMENSION(LENGTH), INTENT(IN)   :: SIN_FST
    REAL(RK), DIMENSION(LENGTH), INTENT(IN)   :: COS_LST
    REAL(RK), DIMENSION(LENGTH), INTENT(IN)   :: SIN_LST
    END SUBROUTINE FFRFT__
  END INTERFACE
  ! ############################################################################################################################# !
  ! (FFTW) COMPLEX DISCRETE FOURIER TRANSFORM (POWER OF TWO INPUT LENGTH INPUT)
  ! (SUBROUTINE) FFT_EXTERNAL_(<LENGTH>, <DIRECTION>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <DIRECTION>            -- (IN)     DIRECTION (IK), FFT_FORWARD = +1_IK OR FFT_INVERSE = -1_IK
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <SEQUENCE>             -- (OUT)    CDFT (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., FR_I, FI_I, ...]
  ! void    fft_external_(int*, int*, double*) ;
  INTERFACE
    MODULE SUBROUTINE FFT_EXTERNAL_(LENGTH, DIRECTION, SEQUENCE) &
      BIND(C, NAME = "fft_external_")
      INTEGER(IK), INTENT(IN) :: LENGTH
      INTEGER(IK), INTENT(IN) :: DIRECTION
      REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(INOUT) :: SEQUENCE
    END SUBROUTINE FFT_EXTERNAL_
  END INTERFACE
  PUBLIC :: FFT_EXTERNAL_
  ! ############################################################################################################################# !
  ! (NRF77) COMPLEX DISCRETE FOURIER TRANSFORM (POWER OF TWO INPUT LENGTH INPUT)
  ! (SUBROUTINE) FFT_RADIX_TWO_(<LENGTH>, <DIRECTION>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <DIRECTION>            -- (IN)     DIRECTION (IK), FFT_FORWARD = +1_IK OR FFT_INVERSE = -1_IK
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <SEQUENCE>             -- (OUT)    CDFT (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., FR_I, FI_I, ...]
  ! void    fft_radix_two_(int*, int*, double*) ;
  INTERFACE
    MODULE SUBROUTINE FFT_RADIX_TWO_(LENGTH, DIRECTION, SEQUENCE) &
      BIND(C, NAME = "fft_radix_two_")
      INTEGER(IK), INTENT(IN) :: LENGTH
      INTEGER(IK), INTENT(IN) :: DIRECTION
      REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(INOUT) :: SEQUENCE
    END SUBROUTINE FFT_RADIX_TWO_
  END INTERFACE
  PUBLIC :: FFT_RADIX_TWO_
  ! ############################################################################################################################# !
  ! (TAKUYA OOURA) COMPLEX DISCRETE FOURIER TRANSFORM (POWER OF TWO INPUT LENGTH INPUT)
  ! (SUBROUTINE) FFT_RADIX_EIGHT_(<LENGTH>, <DIRECTION>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <DIRECTION>            -- (IN)     DIRECTION (IK), FFT_FORWARD = +1_IK OR FFT_INVERSE = -1_IK
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <SEQUENCE>             -- (OUT)    CDFT (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., FR_I, FI_I, ...]
  ! void    fft_radix_eight_(int*, int*, double*) ;
  INTERFACE
    MODULE SUBROUTINE FFT_RADIX_EIGHT_(LENGTH, DIRECTION, SEQUENCE) &
      BIND(C, NAME = "fft_radix_eight_")
      INTEGER(IK), INTENT(IN) :: LENGTH
      INTEGER(IK), INTENT(IN) :: DIRECTION
      REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(INOUT) :: SEQUENCE
    END SUBROUTINE FFT_RADIX_EIGHT_
  END INTERFACE
  PUBLIC :: FFT_RADIX_EIGHT_
  ! ############################################################################################################################# !
  ! (TAKUYA OOURA) COMPLEX DISCRETE FOURIER TRANSFORM (POWER OF TWO INPUT LENGTH INPUT)
  ! (SUBROUTINE) FFT_RADIX_EIGHT__(<LENGTH>, <DIRECTION>, <SEQUENCE>, <IP>, <WORK>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <DIRECTION>            -- (IN)     DIRECTION (IK), FFT_FORWARD = +1_IK OR FFT_INVERSE = -1_IK
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <SEQUENCE>             -- (OUT)    CDFT (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., FR_I, FI_I, ...]
  ! <IP>                   -- (IN)     FFRFT BIT DATA
  ! <WORK>                 -- (IN)     FFRFT TRIG DATA
  INTERFACE
    MODULE SUBROUTINE FFT_RADIX_EIGHT__(LENGTH, DIRECTION, SEQUENCE, IP, WORK)
      INTEGER(IK), INTENT(IN) :: LENGTH
      INTEGER(IK), INTENT(IN) :: DIRECTION
      REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(INOUT) :: SEQUENCE
      INTEGER(IK), DIMENSION(0_IK : 1_IK+INT(SQRT(REAL(LENGTH/2_IK, RK)), IK)), INTENT(IN) :: IP
      REAL(RK), DIMENSION(0_IK : LENGTH/2_IK - 1_IK), INTENT(IN) :: WORK
    END SUBROUTINE FFT_RADIX_EIGHT__
  END INTERFACE
  ! ############################################################################################################################# !
  ! COMPUTE DATA TABLE
  ! (SUBROUTINE) COMPUTE_TABLE_(<LENGTH>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! void    compute_table_(int*) ;
  INTERFACE
    MODULE SUBROUTINE COMPUTE_TABLE_(LENGTH) &
      BIND(C, NAME = "compute_table_")
      INTEGER(IK), INTENT(IN) :: LENGTH
    END SUBROUTINE COMPUTE_TABLE_
  END INTERFACE
  PUBLIC :: COMPUTE_TABLE_
  ! ############################################################################################################################# !
  ! DESTROY DATA TABLE
  ! (SUBROUTINE) DESTROY_TABLE_()
  ! void    destroy_table_(int*) ;
  INTERFACE
    MODULE SUBROUTINE DESTROY_TABLE_() &
      BIND(C, NAME = "destroy_table_")
    END SUBROUTINE DESTROY_TABLE_
  END INTERFACE
  PUBLIC :: DESTROY_TABLE_
  ! ############################################################################################################################# !
  ! SVD
  ! ############################################################################################################################# !
  REAL(RK),    PUBLIC, PARAMETER :: SVD_LEVEL              = 1.0E-9_RK           ! SINGULAR VALUE THRESHOLD LEVEL (ABSOLUTE VALUE)
  INTEGER(IK), PUBLIC, PARAMETER :: SVD_ROW                = 2_IK**12_IK         ! MAX NUMBER OF ROWS
  INTEGER(IK), PUBLIC, PARAMETER :: SVD_COL                = 2_IK**12_IK         ! MAX NUMBER OF ROWS
  INTEGER(IK), PUBLIC, PARAMETER :: SVD_BASIS              = 128_IK              ! MAX NUMBER OF BASIS VECTORS IN THE IMPLICITLY RESTARTED ARNOLDI PROCESS, OPTIMAL VALUE (?)
  INTEGER(IK), PUBLIC, PARAMETER :: SVD_LENGTH             = 64_IK               ! LENGTH OF ARNOLDI FACTORIZATION
  INTEGER(IK), PUBLIC, PARAMETER :: SVD_LOOP               = 256_IK              ! MAX NUMBER OF MAIN LOOP ITERAIONS
  ! ############################################################################################################################# !
  ! SVD (DGESVD)
  ! (SUBROUTINE) SVD_(<NR>, <NC>, <MATRIX>(<NR>, <NC>), <SVD_LIST>(MIN(<NR>, <NC>)), <U_MATRIX>(<NR>, <NR>), <V_MATRIX>(<NC>, <NC>))
  ! <NR>                   -- (IN)     NUMBER OF ROWS (IK)
  ! <NC>                   -- (IN)     NUMBER OF COLS (IK)
  ! <MATRIX>               -- (IN)     INPUT MATRIX(<NR>, <NC>) (RK)
  ! <SVD_LIST>             -- (OUT)    LIST OF SINGULAR VALUES WITH SIZE MIN(<NR>, <NC>) (RK)
  ! <U_MATRIX>             -- (OUT)    L-SINGULAR VECTORS (<NR>, <NR>) (RK)
  ! <V_MATRIX>             -- (OUT)    R-SINGULAR VECTORS (<NC>, <NC>) (RK)
  INTERFACE
    MODULE SUBROUTINE SVD_(NR, NC, MATRIX, SVD_LIST, U_MATRIX, V_MATRIX)
      INTEGER(IK), INTENT(IN) :: NR
      INTEGER(IK), INTENT(IN) :: NC
      REAL(RK), DIMENSION(NR, NC), INTENT(IN) :: MATRIX
      REAL(RK), DIMENSION(MIN(NR, NC)), INTENT(OUT) :: SVD_LIST
      REAL(RK), DIMENSION(NR, NR), INTENT(OUT) :: U_MATRIX
      REAL(RK), DIMENSION(NC, NC), INTENT(OUT) :: V_MATRIX
    END SUBROUTINE SVD_
  END INTERFACE
  PUBLIC :: SVD_
  ! ############################################################################################################################# !
  ! SVD LIST (DGESVD)
  ! (SUBROUTINE) SVD_LIST_(<NR>, <NC>, <MATRIX>(<NR>, <NC>), <SVD_LIST>(MIN(<NR>, <NC>)))
  ! <NR>                   -- (IN)     NUMBER OF ROWS (IK)
  ! <NC>                   -- (IN)     NUMBER OF COLS (IK)
  ! <MATRIX>               -- (IN)     INPUT MATRIX(<NR>, <NC>) (RK)
  ! <SVD_LIST>             -- (OUT)    LIST OF SINGULAR VALUES WITH SIZE MIN(<NR>, <NC>) (RK)
  INTERFACE
    MODULE SUBROUTINE SVD_LIST_(NR, NC, MATRIX, SVD_LIST)
      INTEGER(IK), INTENT(IN) :: NR
      INTEGER(IK), INTENT(IN) :: NC
      REAL(RK), DIMENSION(NR, NC), INTENT(IN) :: MATRIX
      REAL(RK), DIMENSION(MIN(NR, NC)), INTENT(OUT) :: SVD_LIST
    END SUBROUTINE SVD_LIST_
  END INTERFACE
  PUBLIC :: SVD_LIST_
  ! ############################################################################################################################# !
  ! TRUNCATED SVD (ARPACK)
  ! SVD_TRUNCATED_(<NR>,<NC>,<NS>,<MATRIX>(<NR>,<NC>),<LIST>(<NS>),<RVEC>(<NC>,<NR>),<LVEC>(<NR>,<NC>))
  ! <NR>                   -- (IN)     NUMBER OF ROWS (IK)
  ! <NC>                   -- (IN)     NUMBER OF COLS (IK)
  ! <MATRIX>               -- (IN)     INPUT MATRIX(<NR>, <NC>) (RK)
  ! <LIST>                 -- (OUT)    LIST OF SINGULAR VALUES (<NS>) (RK)
  ! <RVEC>                 -- (OUT)    L-SINGULAR VECTORS (<NR>, <NS>) (RK)
  ! <LVEC>                 -- (OUT)    R-SINGULAR VECTORS (<NC>, <NS>) (RK)
  INTERFACE
    MODULE SUBROUTINE SVD_TRUNCATED_(NR, NC, NS, MATRIX, LIST, RVEC, LVEC)
      INTEGER(IK), INTENT(IN) :: NR
      INTEGER(IK), INTENT(IN) :: NC
      INTEGER(IK), INTENT(IN) :: NS
      REAL(RK), DIMENSION(NR, NC), INTENT(IN) :: MATRIX
      REAL(RK), DIMENSION(NS), INTENT(OUT) :: LIST
      REAL(RK), DIMENSION(NC, NS), INTENT(OUT) :: RVEC
      REAL(RK), DIMENSION(NR, NS), INTENT(OUT) :: LVEC
    END SUBROUTINE SVD_TRUNCATED_
  END INTERFACE
  PUBLIC :: SVD_TRUNCATED_
  ! ############################################################################################################################# !
  ! GENERIC SVD
  INTERFACE SVD_
    MODULE PROCEDURE SVD_
    MODULE PROCEDURE SVD_LIST_
    MODULE PROCEDURE SVD_TRUNCATED_
  END INTERFACE SVD_
  ! ############################################################################################################################# !
  ! PROCESS
  ! ############################################################################################################################# !
  ! CONVERT INPUT SEQUENCE (REAL)
  ! (SUBROUTINE) CONVERT_REAL_(<LENGTH>, <R_PART>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <R_PART>               -- (IN)     INPUT SEQUENCE R-PART (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (OUT)    SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...] AND SI_I=0.0_RK FOR ALL I
  ! void    convert_real_(int*, double*, double*) ;
  INTERFACE
    MODULE SUBROUTINE CONVERT_REAL_(LENGTH, R_PART, SEQUENCE) &
      BIND(C, NAME = "convert_real_")
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: R_PART
      REAL(RK), INTENT(OUT), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    END SUBROUTINE CONVERT_REAL_
  END INTERFACE
  PUBLIC :: CONVERT_REAL_
  ! ############################################################################################################################# !
  ! CONVERT INPUT SEQUENCE (COMPLEX)
  ! (SUBROUTINE) CONVERT_COMPLEX_(<LENGTH>, <R_PART>, <I_PART>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <R_PART>               -- (IN)     INPUT SEQUENCE R-PART (RK ARRAY OF LENGTH = <LENGTH>)
  ! <I_PART>               -- (IN)     INPUT SEQUENCE I-PART (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (OUT)    SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! void    convert_complex_(int*, double*, double*, double*) ;
  INTERFACE
    MODULE SUBROUTINE CONVERT_COMPLEX_(LENGTH, R_PART, I_PART, SEQUENCE) &
      BIND(C, NAME = "convert_complex_")
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: R_PART
      REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: I_PART
      REAL(RK), INTENT(OUT), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    END SUBROUTINE CONVERT_COMPLEX_
  END INTERFACE
  PUBLIC :: CONVERT_COMPLEX_
  ! ############################################################################################################################# !
  ! CONVERT INPUT SEQUENCE (REAL OR COMPLEX)
  INTERFACE CONVERT_
    MODULE PROCEDURE CONVERT_REAL_
    MODULE PROCEDURE CONVERT_COMPLEX_
  END INTERFACE
  PUBLIC :: CONVERT_
  ! ############################################################################################################################# !
  ! ROUND UP (ROUND UP TO THE NEXT POWER OF TWO)
  ! (FUNCTION) ROUND_UP_(<NUMBER>)
  ! <NUMBER>               -- (IN)     NUMBER (IK)
  ! <ROUND_UP>             -- (OUT)    NEXT POWER OF TWO NUMBER (IK)
  ! int     round_up_(int*) ;
  INTERFACE
    MODULE INTEGER(IK) FUNCTION ROUND_UP_(NUMBER) &
      BIND(C, NAME = "round_up_")
      INTEGER(IK), INTENT(IN) :: NUMBER
    END FUNCTION ROUND_UP_
  END INTERFACE
  PUBLIC :: ROUND_UP_
  ! ############################################################################################################################# !
  ! ZERO PADDING
  ! (SUBROUTINE) PAD_(<LI>, <LO>, <INPUT>, <OUTPUT>)
  ! <LI>                   -- (IN)     INPUT SEQUENCE LENGTH (IK)
  ! <LO>                   -- (IN)     OUTPUT SEQUENCE LENGTH (IK)
  ! <INPUT>                -- (IN)     INPUT SEQUENCE (RK) OF LENGTH = 2*<LI>
  ! <OUTPUT>               -- (IN)     PADDED SEQUENCE (RK) OF LENGTH = 2*<LO>
  ! void    pad_(int*, int*, double*, double*) ;
  INTERFACE
    MODULE SUBROUTINE PAD_(LI, LO, INPUT, OUTPUT) &
      BIND(C, NAME = "pad_")
      INTEGER(IK), INTENT(IN) :: LI
      INTEGER(IK), INTENT(IN) :: LO
      REAL(RK), DIMENSION(2_IK*LI), INTENT(IN) :: INPUT
      REAL(RK), DIMENSION(2_IK*LO), INTENT(OUT) :: OUTPUT
    END  SUBROUTINE PAD_
  END INTERFACE
  PUBLIC :: PAD_
  ! ############################################################################################################################# !
  ! REMOVE MEAN
  ! (SUBROUTINE) REMOVE_MEAN_(<LENGTH>, <INPUT>, <OUTPUT> )
  ! <LENGTH>               -- (IN)     INPUT SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <INPUT>                -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <OUTPUT>               -- (OUT)    OUTPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! void    remove_mean_(int*, double*, double*) ;
  INTERFACE
    MODULE SUBROUTINE REMOVE_MEAN_(LENGTH, INPUT, OUTPUT) &
      BIND(C, NAME = "remove_mean_")
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(IN) :: INPUT
      REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(OUT) :: OUTPUT
    END SUBROUTINE REMOVE_MEAN_
  END INTERFACE
  PUBLIC :: REMOVE_MEAN_
  ! ############################################################################################################################# !
  ! REMOVE WINDOW MEAN
  ! (SUBROUTINE) REMOVE_WINDOW_MEAN_(<LENGTH>, <TOTAL>, <WINDOW>, <INPUT>, <OUTPUT> )
  ! <LENGTH>               -- (IN)     INPUT SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <INPUT>                -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <OUTPUT>               -- (OUT)    OUTPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! void    remove_window_mean_(int*, double*, double*, double*, double*) ;
  INTERFACE
    MODULE SUBROUTINE REMOVE_WINDOW_MEAN_(LENGTH, TOTAL, WINDOW, INPUT, OUTPUT) &
      BIND(C, NAME = "remove_window_mean_")
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), INTENT(IN) :: TOTAL
      REAL(RK), DIMENSION(LENGTH), INTENT(IN) :: WINDOW
      REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(IN) :: INPUT
      REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(OUT) :: OUTPUT
    END SUBROUTINE REMOVE_WINDOW_MEAN_
  END INTERFACE
  PUBLIC :: REMOVE_WINDOW_MEAN_
  ! ############################################################################################################################# !
  ! APPLY WINDOW
  ! (SUBROUTINE) APPLY_WINDOW_(<LENGTH>, <WINDOW>, <INPUT>, <OUTPUT> )
  ! <LENGTH>               -- (IN)     INPUT SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <INPUT>                -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <OUTPUT>               -- (OUT)    OUTPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! void    apply_window_(int*, double*, double*, double*) ;
  INTERFACE
    MODULE SUBROUTINE APPLY_WINDOW_(LENGTH, WINDOW, INPUT, OUTPUT) &
      BIND(C, NAME = "apply_window_")
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), DIMENSION(LENGTH), INTENT(IN) :: WINDOW
      REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(IN) :: INPUT
      REAL(RK), DIMENSION(2_IK*LENGTH), INTENT(OUT) :: OUTPUT
    END SUBROUTINE APPLY_WINDOW_
  END INTERFACE
  PUBLIC :: APPLY_WINDOW_
  ! ############################################################################################################################# !
  ! MATRIX (GENERATE MATRIX FROM SEQUENCE)
  ! (SUBROUTINE) MATRIX_(<LENGTH>, <SEQUENCE>, <MATRIX>)
  ! <LENGTH>               -- (IN)     INPUT SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK)
  ! <MATRIX>               -- (OUT)    MATRIX (<LENGTH>/2+1, <LENGTH>/2) (RK)
  INTERFACE
    MODULE SUBROUTINE MATRIX_(LENGTH, SEQUENCE, MATRIX)
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), DIMENSION(LENGTH), INTENT(IN) :: SEQUENCE
      REAL(RK), DIMENSION(LENGTH/2_IK+1_IK, LENGTH/2_IK), INTENT(OUT) :: MATRIX
    END SUBROUTINE MATRIX_
  END INTERFACE
  PUBLIC :: MATRIX_
  ! ############################################################################################################################# !
  ! SEQUENCE (ROW) (GENERATE SEQUENCE FROM MATRIX USING 1ST AND LAST ROWS)
  ! (SUBROUTINE) SEQUENCE_ROW_(<LENGTH>, <SEQUENCE>, <MATRIX>)
  ! <LENGTH>               -- (IN)     INPUT SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK)
  ! <MATRIX>               -- (OUT)    MATRIX (<LENGTH>/2+1, <LENGTH>/2) (RK)
  INTERFACE
    MODULE SUBROUTINE SEQUENCE_ROW_(LENGTH, SEQUENCE, MATRIX)
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), DIMENSION(LENGTH), INTENT(OUT) :: SEQUENCE
      REAL(RK), DIMENSION(LENGTH/2_IK+1_IK, LENGTH/2_IK), INTENT(IN) :: MATRIX
    END SUBROUTINE SEQUENCE_ROW_
  END INTERFACE
  PUBLIC :: SEQUENCE_ROW_
  ! ############################################################################################################################# !
  ! SEQUENCE (SUM) (GENERATE SEQUENCE FROM MATRIX USING SUMS OF SKEW DIAGONALS)
  ! (SUBROUTINE) SEQUENCE_ROW_(<LENGTH>, <SEQUENCE>, <MATRIX>)
  ! <LENGTH>               -- (IN)     INPUT SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK)
  ! <MATRIX>               -- (OUT)    MATRIX (<LENGTH>/2+1, <LENGTH>/2) (RK)
  INTERFACE
    MODULE SUBROUTINE SEQUENCE_SUM_(LENGTH, SEQUENCE, MATRIX)
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), DIMENSION(LENGTH), INTENT(OUT) :: SEQUENCE
      REAL(RK), DIMENSION(LENGTH/2_IK+1_IK, LENGTH/2_IK), INTENT(IN) :: MATRIX
    END SUBROUTINE SEQUENCE_SUM_
  END INTERFACE
  PUBLIC :: SEQUENCE_SUM_
  ! ############################################################################################################################# !
  ! FILTER
  ! (SUBROUTINE) FILTER(<LENGTH>, <SEQUENCE>, <LIMIT>)
  ! <LENGTH>               -- (IN)     LENGTH (IK)
  ! <SEQUENCE>             -- (INOUT)  SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  ! <LIMIT>                -- (IN)     NUMBER OF SINGULAR VALUES TO KEEP (IK)
  ! void    filter_(int*, double*, int*) ;
  INTERFACE
    MODULE SUBROUTINE FILTER_(LENGTH, SEQUENCE, LIMIT, SVD_LIST) &
      BIND(C, NAME = "filter_")
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), DIMENSION(LENGTH), INTENT(INOUT) :: SEQUENCE
      INTEGER(IK), INTENT(IN) :: LIMIT
      REAL(RK), DIMENSION(LIMIT), INTENT(OUT) :: SVD_LIST
    END SUBROUTINE FILTER_
  END INTERFACE
  PUBLIC :: FILTER_
  ! ############################################################################################################################# !
  ! PEAKDETECT
  ! ############################################################################################################################# !
  INTEGER(IK), PUBLIC, PARAMETER :: PEAK_WIDTH             = 2_IK                ! PEAK WIDTH
  REAL(RK),    PUBLIC, PARAMETER :: PEAK_LEVEL             = -10.0_RK            ! PEAK THRESHOLD LEVEL
  ! ############################################################################################################################# !
  ! PEAK LIST
  ! (SUBROUTINE) PEAK_LIST_(<LENGTH>, <SEQUENCE>, <PEAK_LIST>)
  ! PEAK_WIDTH             -- (GLOBAL) PEAK WIDTH (IK)
  ! PEAK_LEVEL             -- (GLOBAL) PEAK LEVEL THRESHOLD (RK)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (IN)     SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  ! <PEAK_LIST>            -- (OUT)    PEAK LIST (IK ARRAY OF LENGTH = <LENGTH>), VALUES OF ONE CORRESPOND TO PEAK LOCATIONS
  INTERFACE
    MODULE SUBROUTINE PEAK_LIST_(LENGTH, SEQUENCE, PEAK_LIST)
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), DIMENSION(LENGTH), INTENT(IN) :: SEQUENCE
      INTEGER(IK), DIMENSION(LENGTH), INTENT(OUT) :: PEAK_LIST
    END SUBROUTINE PEAK_LIST_
  END INTERFACE
  PUBLIC :: PEAK_LIST_
  ! ############################################################################################################################# !
  ! TOTAL NUMBER OF PEAKS
  ! (FUNCTION) PEAK_COUNT_(<LENGTH>, <PEAK_LIST>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <PEAK_LIST>            -- (IN)     PEAK LIST (IK ARRAY OF LENGTH <LENGTH>)
  ! <PEAK_COUNT_>          -- (OUT)    TOTAL NUMBER OF PEAKS (IK)
  INTERFACE
    MODULE INTEGER(IK) FUNCTION PEAK_COUNT_(LENGTH, PEAK_LIST)
      INTEGER(IK), INTENT(IN) :: LENGTH
      INTEGER(IK), DIMENSION(LENGTH), INTENT(IN) :: PEAK_LIST
    END FUNCTION PEAK_COUNT_
  END INTERFACE
  PUBLIC :: PEAK_COUNT_
  ! ############################################################################################################################# !
  ! DETECT SEVERAL PEAKS (LIST OF ORDERED PEAK POSITIONS)
  ! (SUBROUTINE) PEAK_DETECT_(<LENGTH>, <SEQUENCE>, <PEAK_LENGTH>, <PEAK_ORDERED_LIST>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (IN)     SEQUENCE (RK ARRAY OF LENGTH = <LENGTH>)
  ! <PEAK_LENGTH>          -- (IN)     NUMBER OF PEAKS TO FIND (IK)
  ! <PEAK_ORDERED_LIST>    -- (OUT)    PEAK POSITIONS (IK ARRAY OF LENGTH = <PEAK_LENGTH>)
  INTERFACE
    MODULE SUBROUTINE PEAK_DETECT_(LENGTH, SEQUENCE, PEAK_LENGTH, PEAK_ORDERED_LIST)
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), DIMENSION(LENGTH), INTENT(IN) :: SEQUENCE
      INTEGER(IK), INTENT(IN) :: PEAK_LENGTH
      INTEGER(IK), DIMENSION(PEAK_LENGTH), INTENT(OUT) :: PEAK_ORDERED_LIST
    END SUBROUTINE PEAK_DETECT_
  END INTERFACE
  PUBLIC :: PEAK_DETECT_
  ! ############################################################################################################################# !
  ! PEAK (RANKED)
  ! (FUNCTION) PEAK_(<LENGTH>, <SEQUENCE>, <PEAK_ID>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <SEQUENCE>             -- (IN)     SEQUENCE (RK ARRAY OF LENGTH <LENGTH>)
  ! <PEAK_ID>              -- (IN)     PEAK RANK (IK)
  ! <PEAK_>                -- (OUT)    PEAK POSITION (IK)
  ! int     peak_(int*, double*, int*) ;
  INTERFACE
    MODULE INTEGER(IK) FUNCTION PEAK_(LENGTH, SEQUENCE, PEAK_ID) &
      BIND(C, NAME = "peak_")
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), DIMENSION(LENGTH), INTENT(IN) :: SEQUENCE
      INTEGER(IK), INTENT(IN) :: PEAK_ID
    END FUNCTION PEAK_
  END  INTERFACE
  PUBLIC :: PEAK_
  ! ############################################################################################################################# !
  ! WINDOW
  ! ############################################################################################################################# !
  ! WINDOW (COSINE)
  ! (SUBROUTINE) WINDOW_COS_(<LENGTH>, <ORDER>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <ORDER>                -- (IN)     WINDOW ORDER (IK)
  ! <WINDOW>               -- (OUT)    WINDOW (RK ARRAY OF LENGTH = <LENGTH>)
  ! void    window_cos_(int*, int*, double*) ;
  INTERFACE
    MODULE SUBROUTINE WINDOW_COS_(LENGTH, ORDER, WINDOW) &
      BIND(C, NAME = "window_cos_")
      INTEGER(IK), INTENT(IN) :: LENGTH
      INTEGER(IK), INTENT(IN) :: ORDER
      REAL(RK), INTENT(OUT), DIMENSION(LENGTH) :: WINDOW
    END SUBROUTINE WINDOW_COS_
  END INTERFACE
  PUBLIC :: WINDOW_COS_
  ! ############################################################################################################################# !
  ! WINDOW (COSINE) (GENERIC, FRACTIONAL ORDER)
  ! (SUBROUTINE) WINDOW_COS_GENERIC_(<LENGTH>, <ORDER>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <ORDER>                -- (IN)     WINDOW ORDER (RK)
  ! <WINDOW>               -- (OUT)    WINDOW (RK ARRAY OF LENGTH = <LENGTH>)
  ! void    window_cos_generic_(int*, double*, double*) ;
  INTERFACE
    MODULE SUBROUTINE WINDOW_COS_GENERIC_(LENGTH, ORDER, WINDOW)&
      BIND(C, NAME = "window_cos_generic_")
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), INTENT(IN) :: ORDER
      REAL(RK), INTENT(OUT), DIMENSION(LENGTH) :: WINDOW
    END SUBROUTINE WINDOW_COS_GENERIC_
  END INTERFACE
  PUBLIC :: WINDOW_COS_GENERIC_
  ! ############################################################################################################################# !
  ! GENERIC WINDOW
  INTERFACE WINDOW_
    MODULE PROCEDURE WINDOW_COS_
    MODULE PROCEDURE WINDOW_COS_GENERIC_
  END INTERFACE
  PUBLIC :: WINDOW_
  ! ############################################################################################################################# !
  ! WINDOW (KAISER)
  ! (SUBROUTINE) WINDOW_KAISER_(<LENGTH>, <ORDER>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <PARAMETER>            -- (IN)     WINDOW ORDER (RK)
  ! <WINDOW>               -- (OUT)    WINDOW (RK ARRAY OF LENGTH = <LENGTH>)
  ! void    window_kaiser_(int*, double*, double*) ;
  INTERFACE
    MODULE SUBROUTINE WINDOW_KAISER_(LENGTH, ORDER, WINDOW) &
      BIND(C, NAME = "window_kaiser_")
      INTEGER(IK), INTENT(IN) :: LENGTH
      REAL(RK), INTENT(IN) :: ORDER
      REAL(RK), INTENT(OUT), DIMENSION(LENGTH) :: WINDOW
    END SUBROUTINE WINDOW_KAISER_
  END INTERFACE
  PUBLIC :: WINDOW_KAISER_
  ! ############################################################################################################################# !
  ! FREQUENCY
  ! ############################################################################################################################# !
  INTEGER(IK), PUBLIC, PARAMETER :: FREQUENCY_FFT          = 0_IK                ! FFT
  INTEGER(IK), PUBLIC, PARAMETER :: FREQUENCY_FFRFT        = 1_IK                ! FFRFT
  INTEGER(IK), PUBLIC, PARAMETER :: FREQUENCY_PARABOLA     = 2_IK                ! PARABOLA
  INTEGER(IK), PUBLIC, PARAMETER :: INITIAL_FFT            = 0_IK                ! INITIAL FREQUENCY (FFT)
  INTEGER(IK), PUBLIC, PARAMETER :: REFINED_FFRFT          = 0_IK                ! REFINED FREQUENCY (FFRFT)
  INTEGER(IK), PUBLIC, PARAMETER :: REFINED_PARABOLA       = 1_IK                ! REFINED FREQUENCY (FFRFT + PARABOLA)
  INTEGER(IK), PUBLIC, PARAMETER :: REFINED_SEARCH         = 2_IK                ! REFINED FREQUENCY (MAXIMUM SEARCH)
  ! ############################################################################################################################# !
  ! FREQUENCY ESTIMATION
  ! (FUNCTION) FREQUENCY_(<FLAG>, <PEAK>, <METHOD>, <LENGTH>, <SEQUENCE>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <PEAK>                 -- (IN)     PEAK NUMBER TO USE (IK), <PEAK> = 0 USE MAXIMUM BIN, <PEAK> = +N USE N'TH PEAK
  ! <METHOD>               -- (IN)     FREQUENCY APPROXIMATION METHOD (IK), FREQUENCY_FFT = 0_IK, FREQUENCY_FFRFT = 1_IK, FREQUECY_PARABOLA = 2_IK
  ! <LENGTH>               -- (IN)     INPUT SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (PROCESSED, I.E. WITH WINDOW AND PADDED WITH ZEROS) (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <FREQUENCY_>           -- (OUT)    FREQUENCY ESTIMATION (RK)
  ! double  frequency_(int*, int*, int*, int*, double*) ;
  INTERFACE
    MODULE REAL(RK) FUNCTION FREQUENCY_(FLAG, PEAK, METHOD, LENGTH, SEQUENCE) &
      BIND(C, NAME = "frequency_")
      INTEGER(IK), INTENT(IN):: FLAG
      INTEGER(IK), INTENT(IN) :: PEAK
      INTEGER(IK), INTENT(IN) :: METHOD
      INTEGER(IK), INTENT(IN):: LENGTH
      REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    END FUNCTION FREQUENCY_
  END INTERFACE
  PUBLIC :: FREQUENCY_
  ! ############################################################################################################################# !
  ! FREQUENCY ESTIMATION (FFT AND FFRFT DATA MEMORIZATION)
  ! (FUNCTION) FREQUENCY__(<FLAG>, <PEAK>, <METHOD>, <LENGTH>, <SEQUENCE>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <PEAK>                 -- (IN)     PEAK NUMBER TO USE (IK), <PEAK> = 0 USE MAXIMUM BIN, <PEAK> = +N USE N'TH PEAK
  ! <METHOD>               -- (IN)     FREQUENCY APPROXIMATION METHOD (IK), FREQUENCY_FFT = 0_IK, FREQUENCY_FFRFT = 1_IK, FREQUECY_PARABOLA = 2_IK
  ! <LENGTH>               -- (IN)     INPUT SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (PROCESSED, I.E. WITH WINDOW AND PADDED WITH ZEROS) (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <FREQUENCY_>           -- (OUT)    FREQUENCY ESTIMATION (RK)
  ! double  frequency_(int*, int*, int*, int*, double*) ;
  INTERFACE
    MODULE REAL(RK) FUNCTION FREQUENCY__(FLAG, PEAK, METHOD, LENGTH, SEQUENCE) &
      BIND(C, NAME = "frequency__")
      INTEGER(IK), INTENT(IN):: FLAG
      INTEGER(IK), INTENT(IN) :: PEAK
      INTEGER(IK), INTENT(IN) :: METHOD
      INTEGER(IK), INTENT(IN):: LENGTH
      REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    END FUNCTION FREQUENCY__
  END INTERFACE
  PUBLIC :: FREQUENCY__
  ! ############################################################################################################################# !
  ! INITIAL FREQUENCY ESTIMATION
  ! (FUNCTION) FREQUENCY_INITIAL_(<FLAG>, <PEAK>, <METHOD>, <LENGTH>, <SEQUENCE>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <PEAK>                 -- (IN)     PEAK NUMBER TO USE (IK), <PEAK> = 0 USE MAXIMUM BIN, <PEAK> = +N USE N'TH PEAK
  ! <METHOD>               -- (IN)     INITIAL FREQUENCY APPROXIMATION METHOD (IK)
  ! <LENGTH>               -- (IN)     INPUT SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (PROCESSED, I.E. WITH WINDOW AND PADDED WITH ZEROS) (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <FREQUENCY_>           -- (OUT)    FREQUENCY ESTIMATION (RK)
  ! double  frequency_initial_(int*, int*, int*, int*, int*, double*, double*, double*) ;
  INTERFACE
    MODULE REAL(RK) FUNCTION FREQUENCY_INITIAL_(FLAG, PEAK, METHOD, LENGTH, PAD, TOTAL, WINDOW, SEQUENCE) &
      BIND(C, NAME = "frequency_initial_")
      INTEGER(IK), INTENT(IN):: FLAG
      INTEGER(IK), INTENT(IN) :: PEAK
      INTEGER(IK), INTENT(IN) :: METHOD
      INTEGER(IK), INTENT(IN):: LENGTH
      INTEGER(IK), INTENT(IN):: PAD
      REAL(RK), INTENT(IN) :: TOTAL
      REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
      REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
    END FUNCTION FREQUENCY_INITIAL_
  END INTERFACE
  PUBLIC :: FREQUENCY_INITIAL_
  ! ############################################################################################################################# !
  ! REFINE FREQUENCY ESTIMATION (COMMON INTERFACE FOR FFRFT, PARABOLA AND SEARCH)
  ! (FUNCTION) FREQUENCY_REFINED_(<FLAG>, <LENGTH>, <TOTAL>, <WINDOW>, <SEQUENCE>, <GUESS>, <INTERVAL>, <LIMIT>, <TOLERANCE>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <METHOD>               -- (IN)     METHOD (SEE GLOBAL REFINED_* CONSTANTS)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <GUESS>                -- (IN)     INITIAL GUESS VALUE (RK)
  ! <INTERVAL>             -- (IN)     SEARCH INTERVAL (RK), GUESS IS IN THE MIDLE
  ! <LIMIT>                -- (IN)     MAXIMUM NUMBER OF ITERATIONS (IK)
  ! <TOLERANCE>            -- (IN)     MAXIMUM TOLERANCE (RK)
  ! <FREQUENCY_REFINED_>   -- (OUT)    REFINED FREQUENCY
  ! double  frequenc_refined_(int*, int*, int*, double*, double*, double*, double*, double*, int*, double*) ;
  INTERFACE
    MODULE REAL(RK) FUNCTION FREQUENCY_REFINED_(FLAG, METHOD, LENGTH, TOTAL, WINDOW, SEQUENCE, GUESS, INTERVAL, LIMIT, TOLERANCE) &
      BIND(C, NAME = "frequency_refined_")
      INTEGER(IK), INTENT(IN):: FLAG
      INTEGER(IK), INTENT(IN):: METHOD
      INTEGER(IK), INTENT(IN):: LENGTH
      REAL(RK), INTENT(IN) :: TOTAL
      REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
      REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
      REAL(RK), INTENT(IN) :: GUESS
      REAL(RK), INTENT(IN) :: INTERVAL
      INTEGER(IK), INTENT(IN) :: LIMIT
      REAL(RK), INTENT(IN) :: TOLERANCE
    END FUNCTION FREQUENCY_REFINED_
  END INTERFACE
  PUBLIC :: FREQUENCY_REFINED_
  ! ############################################################################################################################# !
  ! DECOMPOSITION
  ! ############################################################################################################################# !
  INTEGER(IK), PUBLIC, PARAMETER :: DECOMPOSITION_SUBTRACT = 0_IK                ! DECOMPOSITION BY ITERATIVE SUBTRACTION
  INTEGER(IK), PUBLIC, PARAMETER :: DECOMPOSITION_PEAK     = 1_IK                ! DECOMPOSITION BY PEAKS
  ! ############################################################################################################################# !
  ! SIGNAL DECOMPOSITION
  ! (SUBROUTINE) DECOMPOSITION_(<FLAG>, <METHOD>, <MODE>, <LENGTH>, <PAD>, <TOTAL>, <WINDOW>, <SEQUENCE>, <LOOP>, <FREQUENCY>, <COS_AMP>, <SIN_AMP>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <METHOD>               -- (IN)     FREQUENCY APPROXIMATION METHOD (IK), FREQUENCY_FFT = 0_IK, FREQUENCY_FFRFT = 1_IK, FREQUECY_PARABOLA = 2_IK
  ! <MODE>                 -- (IN)     DECOMPOSTION MODE (IK), <MODE> = DECOMPOSITION_SUBTRACT = 0 OR <MODE> = DECOMPOSITION_PEAK = 1
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <PAD>                  -- (IN)     PADDED SEQUENCE LENGTH (IK), IF PAD > LENGTH, INPUT SEQUENCE IS PADDED, POWER OF TWO, NOT CHECKED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <LOOP>                 -- (IN)     NUMBER OF ITERATIONS/PEAKS (IK)
  ! <FREQUENCY>            -- (OUT)    FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <COS_AMP>              -- (OUT)    COS AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <SIN_AMP>              -- (OUT)    SIN AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! void    decomposition_(int*, int*, int*, int*, int*, double*, double*, double*, int*, double*, double*, double*) ;
  INTERFACE
    MODULE SUBROUTINE DECOMPOSITION_(FLAG, METHOD, MODE, LENGTH, PAD, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP) &
      BIND(C, NAME = "decomposition_")
      INTEGER(IK), INTENT(IN):: FLAG
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
      END SUBROUTINE DECOMPOSITION_
  END INTERFACE
  PUBLIC :: DECOMPOSITION_
  ! ############################################################################################################################# !
  ! FREQUENCY LIST (PERFORM DECOMPOSITION AND RETURN LIST OF FREQUENCIES)
  ! (SUBROUTINE) FREQUENCY_LIST_(<FLAG>, <METHOD>, <MODE>, <LENGTH>, <PAD>, <TOTAL>, <WINDOW>, <SEQUENCE>, <LOOP>, <FREQUENCY>, <COS_AMP>, <SIN_AMP>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <METHOD>               -- (IN)     FREQUENCY APPROXIMATION METHOD (IK), FREQUENCY_FFT = 0_IK, FREQUENCY_FFRFT = 1_IK, FREQUECY_PARABOLA = 2_IK
  ! <MODE>                 -- (IN)     DECOMPOSTION MODE (IK), <MODE> = DECOMPOSITION_SUBTRACT = 0 OR <MODE> = DECOMPOSITION_PEAK = 1
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <PAD>                  -- (IN)     PADDED SEQUENCE LENGTH (IK), IF PAD > LENGTH, INPUT SEQUENCE IS PADDED, POWER OF TWO, NOT CHECKED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <LOOP>                 -- (IN)     NUMBER OF ITERATIONS/PEAKS (IK)
  ! <FREQUENCY>            -- (OUT)    FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! void    frequency_list_(int*, int*, int*, int*, int*, double*, double*, double*, int*, double*) ;
  INTERFACE
    MODULE SUBROUTINE FREQUENCY_LIST_(FLAG, METHOD, MODE, LENGTH, PAD, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY) &
      BIND(C, NAME = "frequency_list_")
      INTEGER(IK), INTENT(IN):: FLAG
      INTEGER(IK), INTENT(IN):: METHOD
      INTEGER(IK), INTENT(IN):: MODE
      INTEGER(IK), INTENT(IN):: LENGTH
      INTEGER(IK), INTENT(IN):: PAD
      REAL(RK), INTENT(IN) :: TOTAL
      REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
      REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
      INTEGER(IK), INTENT(IN) :: LOOP
      REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: FREQUENCY
    END SUBROUTINE FREQUENCY_LIST_
  END INTERFACE
  PUBLIC :: FREQUENCY_LIST_
  ! ############################################################################################################################# !
  ! ESTIMATE AMPLITUDE FOR GIVEN FREQUENCY (DTFT SPECTRA COMPUTATION)
  ! (SUBROUTINE) AMPLITUDE_(<FLAG>, <LENGTH>, <TOTAL>, <WINDOW>, <SEQUENCE>, <FREQUENCY>, <COS_AMP>, <SIN_AMP>, <AMP>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <FREQUENCY>            -- (IN)     FREQUENCY (RK)
  ! <COS_AMP>              -- (OUT)    COS AMPLITUDE (RK)
  ! <SIN_AMP>              -- (OUT)    SIN AMPLITUDE (RK)
  ! <AMP>                  -- (OUT)    ABS AMPLITUDE (RK)
  ! void    amplitude_(int*, double*, double*, double*, double*, double*, double*, double*) ;
  INTERFACE
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
    END SUBROUTINE AMPLITUDE_
  END INTERFACE
  PUBLIC :: AMPLITUDE_
  ! ############################################################################################################################# !
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
  INTERFACE
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
    END SUBROUTINE AMPLITUDE_LIST_
  END INTERFACE
  PUBLIC :: AMPLITUDE_LIST_
  ! ############################################################################################################################# !
  ! SIGNAL DECOMPOSITION (FFT AND FFRFT DATA MEMOZIZATION)
  ! (SUBROUTINE) DECOMPOSITION__(<FLAG>, <METHOD>, <MODE>, <LENGTH>, <PAD>, <TOTAL>, <WINDOW>, <SEQUENCE>, <LOOP>, <FREQUENCY>, <COS_AMP>, <SIN_AMP>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <METHOD>               -- (IN)     FREQUENCY APPROXIMATION METHOD (IK), FREQUENCY_FFT = 0_IK, FREQUENCY_FFRFT = 1_IK, FREQUECY_PARABOLA = 2_IK
  ! <MODE>                 -- (IN)     DECOMPOSTION MODE (IK), <MODE> = DECOMPOSITION_SUBTRACT = 0 OR <MODE> = DECOMPOSITION_PEAK = 1
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <PAD>                  -- (IN)     PADDED SEQUENCE LENGTH (IK), IF PAD > LENGTH, INPUT SEQUENCE IS PADDED, POWER OF TWO, NOT CHECKED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <LOOP>                 -- (IN)     NUMBER OF ITERATIONS/PEAKS (IK)
  ! <FREQUENCY>            -- (OUT)    FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <COS_AMP>              -- (OUT)    COS AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <SIN_AMP>              -- (OUT)    SIN AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! void    decomposition__(int*, int*, int*, int*, int*, double*, double*, double*, int*, double*, double*, double*) ;
  INTERFACE
    MODULE SUBROUTINE DECOMPOSITION__(FLAG, METHOD, MODE, LENGTH, PAD, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY, COS_AMP, SIN_AMP) &
      BIND(C, NAME = "decomposition__")
      INTEGER(IK), INTENT(IN):: FLAG
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
    END SUBROUTINE DECOMPOSITION__
  END INTERFACE
  PUBLIC :: DECOMPOSITION__
  ! ############################################################################################################################# !
  ! FREQUENCY LIST (FFT AND FFRFT DATA MEMORIZATION) (PERFORM DECOMPOSITION AND RETURN LIST OF FREQUENCIES)
  ! (SUBROUTINE) FREQUENCY_LIST__(<FLAG>, <METHOD>, <MODE>, <LENGTH>, <PAD>, <TOTAL>, <WINDOW>, <SEQUENCE>, <LOOP>, <FREQUENCY>, <COS_AMP>, <SIN_AMP>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <METHOD>               -- (IN)     FREQUENCY APPROXIMATION METHOD (IK), FREQUENCY_FFT = 0_IK, FREQUENCY_FFRFT = 1_IK, FREQUECY_PARABOLA = 2_IK
  ! <MODE>                 -- (IN)     DECOMPOSTION MODE (IK), <MODE> = DECOMPOSITION_SUBTRACT = 0 OR <MODE> = DECOMPOSITION_PEAK = 1
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <PAD>                  -- (IN)     PADDED SEQUENCE LENGTH (IK), IF PAD > LENGTH, INPUT SEQUENCE IS PADDED, POWER OF TWO, NOT CHECKED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <SEQUENCE>             -- (IN)     INPUT SEQUENCE (RK ARRAY OF LENGTH = 2_IK*<LENGTH>), <SEQUENCE> = [..., SR_I, SI_I, ...]
  ! <LOOP>                 -- (IN)     NUMBER OF ITERATIONS/PEAKS (IK)
  ! <FREQUENCY>            -- (OUT)    FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! void    frequency_list__(int*, int*, int*, int*, int*, double*, double*, double*, int*, double*) ;
  INTERFACE
    MODULE SUBROUTINE FREQUENCY_LIST__(FLAG, METHOD, MODE, LENGTH, PAD, TOTAL, WINDOW, SEQUENCE, LOOP, FREQUENCY) &
      BIND(C, NAME = "frequency_list__")
      INTEGER(IK), INTENT(IN):: FLAG
      INTEGER(IK), INTENT(IN):: METHOD
      INTEGER(IK), INTENT(IN):: MODE
      INTEGER(IK), INTENT(IN):: LENGTH
      INTEGER(IK), INTENT(IN):: PAD
      REAL(RK), INTENT(IN) :: TOTAL
      REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
      REAL(RK), INTENT(IN), DIMENSION(2_IK*LENGTH) :: SEQUENCE
      INTEGER(IK), INTENT(IN) :: LOOP
      REAL(RK), DIMENSION(LOOP), INTENT(OUT) :: FREQUENCY
    END SUBROUTINE FREQUENCY_LIST__
  END INTERFACE
  PUBLIC :: FREQUENCY_LIST__
  ! ############################################################################################################################# !
  ! FREQUENCY CORRECTION
  ! (SUBROUTINE) FREQUENCY_CORRECTION_(<FLAG>, <METHOD>, <MODE>, <LENGTH>, <PAD>, <TOTAL>, <WINDOW>, <LOOP>, <FREQUENCY>, <COS_AMP>, <SIN_AMP>)
  ! <FLAG>                 -- (IN)     COMPLEX FLAG (IK), 0/1 FOR REAL/COMPLEX INPUT SEQUENCE
  ! <METHOD>               -- (IN)     FREQUENCY APPROXIMATION METHOD (IK), FREQUENCY_FFT = 0_IK, FREQUENCY_FFRFT = 1_IK, FREQUECY_PARABOLA = 2_IK
  ! <MODE>                 -- (IN)     DECOMPOSTION MODE (IK), <MODE> = DECOMPOSITION_SUBTRACT = 0 OR <MODE> = DECOMPOSITION_PEAK = 1
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK), POWER OF TWO, NOT CHECKED
  ! <PAD>                  -- (IN)     PADDED SEQUENCE LENGTH (IK), IF PAD > LENGTH, INPUT SEQUENCE IS PADDED, POWER OF TWO, NOT CHECKED
  ! <TOTAL>                -- (IN)     SUM(WINDOW) (RK)
  ! <WINDOW>               -- (IN)     WINDOW ARRAY (RK ARRAY OF LENGTH = <LENGTH>)
  ! <LOOP>                 -- (IN)     NUMBER OF ITERATIONS/PEAKS (IK)
  ! <FREQUENCY>            -- (IN)     FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <COS_AMP>              -- (IN)     COS AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <SIN_AMP>              -- (IN)     SIN AMPLITUDE ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! <FREQUENCY>            -- (OUT)    CORRECTED FREQUENCY ARRAY (RK ARRAY OF LENGTH = <LOOP>)
  ! void    frequency_correction_(int*, int*, int*, int*, int*, double*, double*, int*, double*, double*, double*) ;
  INTERFACE
    MODULE SUBROUTINE FREQUENCY_CORRECTION_(FLAG, METHOD, MODE, LENGTH, TOTAL, WINDOW, LOOP, FREQUENCY, COS_AMP, SIN_AMP)
      INTEGER(IK), INTENT(IN):: FLAG
      INTEGER(IK), INTENT(IN):: METHOD
      INTEGER(IK), INTENT(IN):: MODE
      INTEGER(IK), INTENT(IN):: LENGTH
      REAL(RK), INTENT(IN) :: TOTAL
      REAL(RK), INTENT(IN), DIMENSION(LENGTH) :: WINDOW
      INTEGER(IK), INTENT(IN) :: LOOP
      REAL(RK), DIMENSION(LOOP), INTENT(INOUT) :: FREQUENCY
      REAL(RK), DIMENSION(LOOP), INTENT(INOUT) :: COS_AMP
      REAL(RK), DIMENSION(LOOP), INTENT(INOUT) :: SIN_AMP
    END SUBROUTINE FREQUENCY_CORRECTION_
  END INTERFACE
  PUBLIC :: FREQUENCY_CORRECTION_
  ! ############################################################################################################################# !
  ! OPTIMIZATION
  ! ############################################################################################################################# !
  ! LEAST SQUARES (SVD)
  ! (SUBROUTINE) LEAST_SQUARES_(<NR>, <NC>, <MATRIX>(<NR>, <NC>), <VECTOR>(<NR>), <SOLUTION>(<NC>))
  ! <NR>                   -- (IN)     NUMBER OF ROWS (IK)
  ! <NC>                   -- (IN)     NUMBER OF COLS (IK)
  ! <MATRIX>               -- (IN)     INPUT DATA MATRIX (<NR>, <NC>) (RK)
  ! <VECTOR>               -- (IN)     INPUT VECTOR (<NR>) (RK)
  ! <SOLUTION>             -- (OUT)    LS SOLUTION (<NC>) (RK)
  INTERFACE
    MODULE SUBROUTINE LEAST_SQUARES_(NR, NC, MATRIX, VECTOR, SOLUTION)
      INTEGER(IK), INTENT(IN) :: NR
      INTEGER(IK), INTENT(IN) :: NC
      REAL(RK), DIMENSION(NR, NC), INTENT(IN) :: MATRIX
      REAL(RK), DIMENSION(NR), INTENT(IN) :: VECTOR
      REAL(RK), DIMENSION(NC), INTENT(OUT) :: SOLUTION
    END SUBROUTINE LEAST_SQUARES_
  END INTERFACE
  PUBLIC :: LEAST_SQUARES_
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
  INTERFACE
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
    END SUBROUTINE FIT_
  END INTERFACE
  PUBLIC :: FIT_
  ! ############################################################################################################################# !
  ! BINARY SEARCH MAXIMIZATION
  ! (FUNCTION) BINARY_(<FUN>, <GUESS>, <INTERVAL>, <LIMIT>, <TOLERANCE>)
  ! <FUN>                  -- (IN)     FUNCTION TO MAXIMIZE (RK) -> (RK)
  ! <GUESS>                -- (IN)     INITIAL GUESS VALUE (RK)
  ! <INTERVAL>             -- (IN)     SEARCH INTERVAL (RK), GUESS IS IN THE MIDLE
  ! <LIMIT>                -- (IN)     MAXIMUM NUMBER OF ITERATIONS (IK)
  ! <TOLERANCE>            -- (IN)     MAXIMUM TOLERANCE (RK)
  ! <BINARY_>              -- (OUT)    MAXIMUM POSITION
  INTERFACE
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
    END FUNCTION BINARY_
  END INTERFACE
  PUBLIC :: BINARY_
  ! ############################################################################################################################# !
  ! GOLDEN SEARCH MAXIMIZATION
  ! (FUNCTION) GOLDEN_(<FUN>, <GUESS>, <INTERVAL>, <LIMIT>, <TOLERANCE>)
  ! <FUN>                  -- (IN)     FUNCTION TO MAXIMIZE (RK) -> (RK)
  ! <GUESS>                -- (IN)     INITIAL GUESS VALUE (RK)
  ! <INTERVAL>             -- (IN)     SEARCH INTERVAL (RK), GUESS IS IN THE MIDLE
  ! <LIMIT>                -- (IN)     MAXIMUM NUMBER OF ITERATIONS (IK)
  ! <TOLERANCE>            -- (IN)     MAXIMUM TOLERANCE (RK)
  ! <GOLDEN_>              -- (OUT)    MAXIMUM POSITION
  INTERFACE
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
    END FUNCTION GOLDEN_
  END INTERFACE
  PUBLIC :: GOLDEN_
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
  INTERFACE
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
    END FUNCTION
  END INTERFACE
  PUBLIC :: BINARY_AMPLITUDE_
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
  INTERFACE
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
    END FUNCTION
  END INTERFACE
  PUBLIC :: GOLDEN_AMPLITUDE_
  ! ############################################################################################################################# !
END MODULE SIGNAL
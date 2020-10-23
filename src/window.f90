
#include "signal.inc"

SUBMODULE (SIGNAL) WINDOW
  IMPLICIT NONE
  CONTAINS
  ! ############################################################################################################################# !
  ! WINDOW (COSINE)
  ! (SUBROUTINE) WINDOW_COS_(<LENGTH>, <ORDER>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <ORDER>                -- (IN)     WINDOW ORDER (IK)
  ! <WINDOW>               -- (OUT)    WINDOW (RK ARRAY OF LENGTH = <LENGTH>)
  ! void    window_cos_(int* length, int* order, double* window) ;
  MODULE SUBROUTINE WINDOW_COS_(LENGTH, ORDER, WINDOW) &
    BIND(C, NAME = "window_cos_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    INTEGER(IK), INTENT(IN) :: ORDER
    REAL(RK), INTENT(OUT), DIMENSION(LENGTH) :: WINDOW
    INTEGER(IK) :: I
    REAL(RK) :: FACTOR
    FACTOR = 2.0_RK**ORDER*FACTORIAL_(ORDER)**2_IK/FACTORIAL_(2_IK*ORDER)
    WINDOW = FACTOR*(1.0_RK+COS(TWO_PI*(1.0_RK/REAL(LENGTH,RK)*REAL([(I, I = 0_IK, LENGTH-1_IK, 1_IK)], RK)-0.5_RK)))**ORDER
  END SUBROUTINE WINDOW_COS_
  ! ############################################################################################################################# !
  ! WINDOW (COSINE) (GENERIC, FRACTIONAL ORDER)
  ! (SUBROUTINE) WINDOW_COS_GENERIC_(<LENGTH>, <ORDER>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <ORDER>                -- (IN)     WINDOW ORDER (RK)
  ! <WINDOW>               -- (OUT)    WINDOW (RK ARRAY OF LENGTH = <LENGTH>)
  ! void    window_cos_generic_(int* length, double* order, double* window) ;
  MODULE SUBROUTINE WINDOW_COS_GENERIC_(LENGTH, ORDER, WINDOW) &
    BIND(C, NAME = "window_cos_generic_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), INTENT(IN) :: ORDER
    REAL(RK), INTENT(OUT), DIMENSION(LENGTH) :: WINDOW
    INTEGER(IK) :: I
    REAL(RK) :: FACTOR
    FACTOR = 2.0_RK**ORDER*GAMMA_(ORDER+1.0_RK)**2_IK/GAMMA_(2.0_RK*ORDER+1.0_RK)
    WINDOW = FACTOR*(1.0_RK+COS(TWO_PI*(1.0_RK/REAL(LENGTH,RK)*REAL([(I, I = 0_IK, LENGTH-1_IK, 1_IK)], RK)-0.5_RK)))**ORDER
  END SUBROUTINE WINDOW_COS_GENERIC_
  ! ############################################################################################################################# !
  ! WINDOW (KAISER)
  ! (SUBROUTINE) WINDOW_KAISER_(<LENGTH>, <ORDER>, <SEQUENCE>)
  ! <LENGTH>               -- (IN)     SEQUENCE LENGTH (IK)
  ! <PARAMETER>            -- (IN)     WINDOW ORDER (RK)
  ! <WINDOW>               -- (OUT)    WINDOW (RK ARRAY OF LENGTH = <LENGTH>)
  ! void    window_kaiser_(int* length, double* order, double* window) ;
  MODULE SUBROUTINE WINDOW_KAISER_(LENGTH, ORDER, WINDOW) &
    BIND(C, NAME = "window_kaiser_")
    INTEGER(IK), INTENT(IN) :: LENGTH
    REAL(RK), INTENT(IN) :: ORDER
    REAL(RK), INTENT(OUT), DIMENSION(LENGTH) :: WINDOW
    INTEGER(IK) :: I
    REAL(RK) :: FACTOR
    FACTOR = 1.0_RK/BESSEL_(ONE_PI*ORDER)
    WINDOW = 1.0_RK/REAL(LENGTH, RK)*REAL([(I-1_IK, I = 1_IK, LENGTH, 1_IK)], RK)-0.5_RK
    DO I = 1_IK, LENGTH, 1_IK
      WINDOW(I) = BESSEL_(ONE_PI*ORDER*SQRT(1.0_RK-2.0_RK*WINDOW(I))*SQRT(1.0_RK+2.0_RK*WINDOW(I)))
    END DO
    WINDOW = FACTOR*WINDOW
  END SUBROUTINE WINDOW_KAISER_
  ! ############################################################################################################################# !
END SUBMODULE WINDOW
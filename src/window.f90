
#include "signal.inc"

submodule (signal) window
  implicit none
  contains
  ! ############################################################################################################################# !
  ! window (cosine)
  ! (subroutine) window_cos_(<length>, <order>, <sequence>)
  ! <length>               -- (in)     sequence length (ik)
  ! <order>                -- (in)     window order (ik)
  ! <window>               -- (out)    window (rk array of length = <length>)
  ! void    window_cos_(int* length, int* order, double* window) ;
  module subroutine window_cos_(length, order, window) &
    bind(c, name = "window_cos_")
    integer(ik), intent(in) :: length
    integer(ik), intent(in) :: order
    real(rk), intent(out), dimension(length) :: window
    integer(ik) :: i
    real(rk) :: factor
    factor = 2.0_rk**order*factorial_(order)**2_ik/factorial_(2_ik*order)
    window = factor*(1.0_rk+cos(two_pi*(1.0_rk/real(length,rk)*real([(i, i = 0_ik, length-1_ik, 1_ik)], rk)-0.5_rk)))**order
  end subroutine window_cos_
  ! ############################################################################################################################# !
  ! window (cosine) (generic, fractional order)
  ! (subroutine) window_cos_generic_(<length>, <order>, <sequence>)
  ! <length>               -- (in)     sequence length (ik)
  ! <order>                -- (in)     window order (rk)
  ! <window>               -- (out)    window (rk array of length = <length>)
  ! void    window_cos_generic_(int* length, double* order, double* window) ;
  module subroutine window_cos_generic_(length, order, window) &
    bind(c, name = "window_cos_generic_")
    integer(ik), intent(in) :: length
    real(rk), intent(in) :: order
    real(rk), intent(out), dimension(length) :: window
    integer(ik) :: i
    real(rk) :: factor
    factor = 2.0_rk**order*gamma_(order+1.0_rk)**2_ik/gamma_(2.0_rk*order+1.0_rk)
    window = factor*(1.0_rk+cos(two_pi*(1.0_rk/real(length,rk)*real([(i, i = 0_ik, length-1_ik, 1_ik)], rk)-0.5_rk)))**order
  end subroutine window_cos_generic_
  ! ############################################################################################################################# !
  ! window (kaiser)
  ! (subroutine) window_kaiser_(<length>, <order>, <sequence>)
  ! <length>               -- (in)     sequence length (ik)
  ! <parameter>            -- (in)     window order (rk)
  ! <window>               -- (out)    window (rk array of length = <length>)
  ! void    window_kaiser_(int* length, double* order, double* window) ;
  module subroutine window_kaiser_(length, order, window) &
    bind(c, name = "window_kaiser_")
    integer(ik), intent(in) :: length
    real(rk), intent(in) :: order
    real(rk), intent(out), dimension(length) :: window
    integer(ik) :: i
    real(rk) :: factor
    factor = 1.0_rk/bessel_(one_pi*order)
    window = 1.0_rk/real(length, rk)*real([(i-1_ik, i = 1_ik, length, 1_ik)], rk)-0.5_rk
    do i = 1_ik, length, 1_ik
      window(i) = bessel_(one_pi*order*sqrt(1.0_rk-2.0_rk*window(i))*sqrt(1.0_rk+2.0_rk*window(i)))
    end do
    window = factor*window
  end subroutine window_kaiser_
  ! ############################################################################################################################# !
end submodule window
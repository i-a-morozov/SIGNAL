
#include "signal.inc"

submodule (signal) optimization
  implicit none
  contains
  ! ############################################################################################################################# !
  ! least squares (svd)
  ! (subroutine) least_squares_(<nr>, <nc>, <matrix>(<nr>, <nc>), <vector>(<nr>), <solution>(<nc>))
  ! <nr>                   -- (in)     number of rows (ik)
  ! <nc>                   -- (in)     number of cols (ik)
  ! <matrix>               -- (in)     input data matrix (<nr>, <nc>) (rk)
  ! <vector>               -- (in)     input vector (<nr>) (rk)
  ! <solution>             -- (out)    ls solution (<nc>) (rk)
  module subroutine least_squares_(nr, nc, matrix, vector, solution)
    integer(ik), intent(in) :: nr
    integer(ik), intent(in) :: nc
    real(rk), dimension(nr, nc), intent(in) :: matrix
    real(rk), dimension(nr), intent(in) :: vector
    real(rk), dimension(nc), intent(out) :: solution
    real(rk), dimension(min(nr, nc)) :: svd_list
    real(rk), dimension(nr, nr) :: u_matrix
    real(rk), dimension(nc, nc) :: v_matrix
    real(rk), dimension(nr, nc) :: copy
    integer(ik) :: i
    real(rk) :: v1(nr), v2(nr), v3(nc), v4(nc)
    call svd_(nr, nc, matrix, svd_list, u_matrix, v_matrix)
    copy = 0.0_rk
    do i = 1_ik, int(size(svd_list), ik), 1_ik
      if (svd_list(i) >= svd_level) copy(i, i) = 1.0_rk/svd_list(i)
    end do
    v1 = vector
    call dgemv('t',nr,nr,1.0_rk,u_matrix,nr,v1,1_ik,0.0_rk,v2,1_ik)
    call dgemv('t',nr,nc,1.0_rk,copy,nr,v2,1_ik,0.0_rk,v3,1_ik)
    call dgemv('n',nc,nc,1.0_rk,v_matrix,nc,v3,1_ik,0.0_rk,v4,1_ik)
    solution = real(v4, rk)
  end subroutine least_squares_
  ! ############################################################################################################################# !
  ! fit (harmonic signal)
  ! (subroutine) fit_(<length>, <sequence>, <loop>, <frequency>, <mean>, <cos_amp>, <sin_amp>, <error>)
  ! <length>               -- (in)     sequence length (ik), power of two
  ! <sequence>             -- (in)     input sequence (rk array of length = <length>)
  ! <loop>                 -- (in)     number of harmonics (ik)
  ! <frequency>            -- (in)     frequency array (rk array of length = <loop>)
  ! <mean>                 -- (out)    mean value
  ! <cos_amp>              -- (out)    cos amplitude array (rk array of length = <loop>)
  ! <sin_amp>              -- (out)    sin amplitude array (rk array of length = <loop>)
  ! <error>                -- (out)    error
  ! void    fit_(int* length, double* sequence, int* loop, double* frequency, double* mean, double* cos_amp, double* sin_amp, double* error) ;
  module subroutine fit_(length, sequence, loop, frequency, mean, cos_amp, sin_amp, error) &
    bind(c, name = "fit_")
    integer(ik), intent(in) :: length
    real(rk), dimension(length), intent(in) :: sequence
    integer(ik), intent(in) :: loop
    real(rk), dimension(loop), intent(in) :: frequency
    real(rk), intent(out) :: mean
    real(rk), dimension(loop), intent(out) :: cos_amp
    real(rk), dimension(loop), intent(out) :: sin_amp
    real(rk), intent(out) :: error
    real(rk), dimension(length) :: range
    real(rk), dimension(length, 2_ik*loop+1_ik) :: matrix
    real(rk), dimension(2_ik*loop+1_ik) :: solution
    real(rk), dimension(length) :: fit
    integer(ik) :: i
    range = two_pi*real([(i, i = 1_ik, length, 1_ik)], rk)
    matrix = 1.0_rk
    do i = 1_ik, loop, 1_ik
      matrix(:,2_ik*i)      = cos(frequency(i)*range)
      matrix(:,2_ik*i+1_ik) = sin(frequency(i)*range)
    end do
    call least_squares_(length, 2_ik*loop+1_ik, matrix, sequence, solution)
    mean = solution(1_ik)
    cos_amp = solution(2_ik:2_ik*loop+1_ik:2_ik)
    sin_amp = solution(3_ik:2_ik*loop+1_ik:2_ik)
    fit = mean
    do i = 1_ik, loop, 1_ik
      fit = fit + cos_amp(i)*cos(frequency(i)*range) + sin_amp(i)*sin(frequency(i)*range)
    end do
    error = norm2(sequence-fit)
  end subroutine fit_
  ! ############################################################################################################################# !
  ! fit (parabola) y = a*x**2 + b*x + c
  ! (subroutine) fit_parabola_(<length>, <x>, <y>, <a>, <b>, <c>, <maximum>)
  ! <length>               -- (in)     sequence length (ik), power of two
  ! <x>                    -- (in)     x (rk array of length = <length>)
  ! <y>                    -- (in)     y (rk array of length = <length>)
  ! <a>                    -- (out)    a (rk)
  ! <b>                    -- (out)    b (rk)
  ! <c>                    -- (out)    c (rk)
  ! <maximum>              -- (out)    maximum (minimum) position (rk)
  ! void    fit_parabola_(int* length, double* x, double* y, double* a, double* b, double* c, double* maximum) ;
  module subroutine fit_parabola_(length, x, y, a, b, c, maximum) &
    bind(c, name = "fit_parabola_")
    integer(ik), intent(in) :: length
    real(rk), dimension(length), intent(in) :: x
    real(rk), dimension(length), intent(in) :: y
    real(rk), intent(out) :: a
    real(rk), intent(out) :: b
    real(rk), intent(out) :: c
    real(rk), intent(out) :: maximum
    real(rk), dimension(length, 3_ik) :: matrix
    real(rk), dimension(3_ik) :: solution
    matrix(:,1_ik) = x**2_ik
    matrix(:,2_ik) = x
    matrix(:,3_ik) = 1.0_rk
    call least_squares_(length, 3_ik, matrix, y, solution)
    a = solution(1_ik)
    b = solution(2_ik)
    c = solution(3_ik)
    maximum = - b/(2.0_rk*a)
  end subroutine fit_parabola_
  ! ############################################################################################################################# !
  ! binary search maximization
  ! (function) binary_(<fun>, <guess>, <interval>, <limit>, <tolerance>)
  ! <fun>                  -- (in)     function to maximize (rk) -> (rk)
  ! <guess>                -- (in)     initial guess value (rk)
  ! <interval>             -- (in)     search interval (rk), guess is in the midle
  ! <limit>                -- (in)     maximum number of iterations (ik)
  ! <tolerance>            -- (in)     maximum tolerance (rk)
  ! <binary_>              -- (out)    maximum position
  module real(rk) function binary_(fun, guess, interval, limit, tolerance)
    interface
      real(rk) function search(x)
        import :: rk
        real(rk), intent(in) :: x
      end function search
    end interface
    procedure(search) :: fun
    real(rk), intent(in) :: guess
    real(rk), intent(in) :: interval
    integer(ik), intent(in) :: limit
    real(rk), intent(in) :: tolerance
    real(rk) :: delta
    real(rk) :: xl, xr, xx
    real(rk) :: fl, fr
    integer(ik) :: i
    delta = interval/2.0_rk
    xl = guess-delta
    xr = guess+delta
    fl = fun(fl)
    fr = fun(fr)
    do i = 1_ik, limit, 1_ik
      if (fl > fr) then
        xx = xl
      else
        xx = xr
      end if
      xl = xx-delta
      xr = xx+delta
      delta = delta/2.0_rk
      fl = fun(xl)
      fr = fun(xr)
      if (abs(fl-fr) < tolerance) exit
    end do
    if (fl > fr) then
      binary_ = xl
    else
      binary_ = xr
    end if
  end function binary_
  ! ############################################################################################################################# !
  ! golden search maximization
  ! (function) golden_(<fun>, <guess>, <interval>, <limit>, <tolerance>)
  ! <fun>                  -- (in)     function to maximize (rk) -> (rk)
  ! <guess>                -- (in)     initial guess value (rk)
  ! <interval>             -- (in)     search interval (rk), guess is in the midle
  ! <limit>                -- (in)     maximum number of iterations (ik)
  ! <tolerance>            -- (in)     maximum tolerance (rk)
  ! <golden_>              -- (out)    maximum position
  module real(rk) function golden_(fun, guess, interval, limit, tolerance)
    interface
      real(rk) function search(x)
        import :: rk
        real(rk), intent(in) :: x
      end function search
    end interface
    procedure(search) :: fun
    real(rk), intent(in) :: guess
    real(rk), intent(in) :: interval
    integer(ik), intent(in) :: limit
    real(rk), intent(in) :: tolerance
    real(rk), parameter :: golden = (1.0_rk+sqrt(5.0_rk))/2.0_rk
    real(rk), parameter :: psi = 1.0_rk-1.0_rk/golden
    real(rk), parameter :: phi = 1.0_rk/golden
    real(rk) :: delta
    real(rk) :: xl, xr
    real(rk) :: fl, fr
    real(rk) :: m1, m2
    integer(ik) :: i
    delta = interval/2.0_rk
    xl = guess-delta
    xr = guess+delta
    m1 = xl+psi*(xr-xl)
    m2 = xl+phi*(xr-xl)
    fl = fun(m1)
    fr = fun(m2)
    do i = 1_ik, limit, 1_ik
      if (i > limit) exit
      if (abs(fl-fr) < tolerance) exit
      if (fl > fr) then
        xr = m2
        m2 = m1
        m1 = psi*xl+phi*m2
        fr = fl
        fl = fun(m1)
      else
        xl = m1
        m1 = m2
        m2 = phi*m1+psi*xr
        fl = fr
        fr = fun(m2)
      end if
    end do
    if (fl > fr) then
      golden_ = xl
    else
      golden_ = xr
    end if
  end function golden_
  ! ############################################################################################################################# !
end submodule optimization
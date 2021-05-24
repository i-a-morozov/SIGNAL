! example-09: frequency correction
program example

  use signal

  implicit none

  integer(ik), parameter           :: length = 2_ik**9_ik   ! input signal length
  integer(ik)                      :: flag                  ! complex flag (0/1)
  real(rk)                         :: range_min             ! (min) frequency range
  real(rk)                         :: range_max             ! (max) frequency range
  integer(ik)                      :: method                ! frequency estimation method
  integer(ik)                      :: mode                  ! decomposition mode
  integer(ik)                      :: order                 ! cosine window order
  real(rk), dimension(length)      :: window                ! cosine window data
  real(rk)                         :: total                 ! sum of window elements
  real(rk), dimension(length)      :: signal_r, signal_i    ! signal real and complex parts
  real(rk), dimension(2_ik*length) :: signal                ! input signal, [..., sr_i, si_i, ...]
  integer(ik), parameter           :: loop = 4_ik           ! number of harmonics
  real(rk), dimension(loop)        :: frequency             ! frequencies
  real(rk), dimension(loop)        :: cos_amp               ! cos-amplitudes
  real(rk), dimension(loop)        :: sin_amp               ! sin-amplitudes

  real(rk), parameter              :: f(4) = [1._rk*0.123456789_rk,2._rk*0.123456789_rk,3._rk*0.123456789_rk,4._rk*0.123456789_rk]
  real(rk), parameter              :: a(4) = [1.0_rk,1.0e-1_rk,1.e-3_rk,0.0_rk]
  real(rk), parameter              :: b(4) = [1.0e-1_rk,5.0e-2_rk,0.0_rk,1.e-4_rk]

  integer(ik)                      :: i

  ! set complex flag (0/1 for real/complex signal)
  flag = 1_ik
  range_min = 0.00_rk
  range_max = 0.99_rk

  ! set signal real and imaginary parts
  signal_r = 1.0_rk
  signal_i = 0.0_rk
  do i = 1_ik, length, 1_ik
    signal_r(i) = &
      a(1)*cos(two_pi*f(1)*real(i, rk)) + b(1)*sin(two_pi*f(1)*real(i, rk)) + &
      a(2)*cos(two_pi*f(2)*real(i, rk)) + b(2)*sin(two_pi*f(2)*real(i, rk)) + &
      a(3)*cos(two_pi*f(3)*real(i, rk)) + b(3)*sin(two_pi*f(3)*real(i, rk)) + &
      a(4)*cos(two_pi*f(4)*real(i, rk)) + b(4)*sin(two_pi*f(4)*real(i, rk))
    signal_i(i) = &
      b(1)*cos(two_pi*f(1)*real(i, rk)) - a(1)*sin(two_pi*f(1)*real(i, rk)) + &
      b(2)*cos(two_pi*f(2)*real(i, rk)) - a(2)*sin(two_pi*f(2)*real(i, rk)) + &
      b(3)*cos(two_pi*f(3)*real(i, rk)) - a(3)*sin(two_pi*f(3)*real(i, rk)) + &
      b(4)*cos(two_pi*f(4)*real(i, rk)) - a(4)*sin(two_pi*f(4)*real(i, rk))
  end do

  ! format test signal
  signal = 0.0_rk
  if(flag == 0_ik) then
    call convert_(length, signal_r, signal)
  else
    call convert_(length, signal_r, signal_i, signal)
  end if

  ! set window and window sum
  order = 1_ik
  call window_(length, order, window)
  total = sum(window)

  ! print input signal parameters
  write(*, '(a)') "input"
  write(*, '(3e32.16)') f(1), a(1), b(1)
  write(*, '(3e32.16)') f(2), a(2), b(2)
  write(*, '(3e32.16)') f(3), a(3), b(3)
  write(*, '(3e32.16)') f(4), a(4), b(4)
  write(*, *)

  ! estimate frequencies and amplitudes using maximum bin and subtraction
  ! frequencies are estimated from largest bin
  write(*, '(a)') "decomposition (subtraction)"
  method = frequency_parabola
  mode = decomposition_subtract
  frequency = 0.0_rk
  cos_amp = 0.0_rk
  sin_amp = 0.0_rk
  call decomposition_(flag,range_min,range_max,method,mode,length,length,total,window,signal,loop,frequency,cos_amp,sin_amp)
  do i = 1_ik, loop, 1_ik
    write(*, '(3e32.16)') frequency(i), cos_amp(i), sin_amp(i)
  end do
  write(*, *)

  ! frequency correction
  ! this procedure can be applied to find corrections to estimated frequencies and amplitudes
  ! first, signal (complete) decompostition is performed and parameter estimations are computed
  ! next, based on estimated parameters model signal is generated and (complete) decomposition is performed
  ! and corrections are commputed
  block
    real(rk), dimension(loop) :: c_frequency
    real(rk), dimension(loop) :: c_cos_amp
    real(rk), dimension(loop) :: c_sin_amp
    c_frequency = frequency
    c_cos_amp = cos_amp
    c_sin_amp = sin_amp
    call frequency_correction_(flag,range_min,range_max,method,mode,length,length,total,window,loop,c_frequency,c_cos_amp,c_sin_amp)
    write(*, '(a)') "original errors"
    write(*, '(a,4e32.16)') "fre", abs(frequency-f)
    write(*, '(a,4e32.16)') "cos", abs(cos_amp-a)
    write(*, '(a,4e32.16)') "sin", abs(sin_amp-b)
    write(*, *)
    write(*, '(a)') "corrected errors"
    write(*, '(a,4e32.16)') "fre", abs(c_frequency-f)
    write(*, '(a,4e32.16)') "cos", abs(c_cos_amp-a)
    write(*, '(a,4e32.16)') "sin", abs(c_sin_amp-b)
  end block
  write(*, *)

end program example

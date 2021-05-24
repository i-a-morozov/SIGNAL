! example-08: frequency estimation (signal with noise)
program example

  use signal

  implicit none

  integer(ik), parameter           :: length = 2_ik**10_ik  ! input signal length
  integer(ik)                      :: flag                  ! complex flag (0/1)
  real(rk)                         :: range_min             ! (min) frequency range
  real(rk)                         :: range_max             ! (max) frequency range
  integer(ik)                      :: peak                  ! fourie spectra peak id
  integer(ik)                      :: method                ! frequency estimation method
  integer(ik)                      :: order                 ! cosine window order
  real(rk), dimension(length)      :: window                ! cosine window data
  real(rk)                         :: total                 ! sum of window elements
  real(rk)                         :: frequency             ! exact signal frequency
  real(rk), dimension(length)      :: signal_r, signal_i    ! signal real and complex parts
  real(rk), dimension(2_ik*length) :: signal                ! input signal, [..., sr_i, si_i, ...]
  real(rk)                         :: result                ! frequency estimation
  real(rk), dimension(length/2_ik) :: n, m                  ! noise
  integer(ik), parameter           :: limit = 4_ik          ! number of singular values to keep
  real(rk), dimension(limit)       :: list                  ! svd list
  real(rk)                         :: cos_amp, sin_amp, amp

  integer(ik)                      :: i

  flag = 0_ik
  range_min = 0.00_rk
  range_max = 0.49_rk

  ! set signal exact frequency, real and imaginary parts
  frequency = 0.123456789_rk
  do i = 1_ik, length, 1_ik
    signal_r(i) = 1.0_rk*sin(two_pi*frequency*i)+0.05*sin(two_pi*2.0_rk*frequency*i)
    signal_i(i) = 0.0_rk
  end do

  ! add noise
  call random_number(n)
  call random_number(m)

  signal_r(1_ik:length:2_ik) = signal_r(1_ik:length:2_ik) + 0.05*sqrt(-2.0_rk*log(n))*cos(two_pi*m)
  signal_r(2_ik:length:2_ik) = signal_r(2_ik:length:2_ik) + 0.05*sqrt(-2.0_rk*log(n))*sin(two_pi*m)

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

  ! frequency estimation method
  method = frequency_parabola

  ! estimate frequency
  block
    write(*, '(a)') "filter(-)"
    peak = +1_ik
    result = frequency_(flag, range_min, range_max, peak, method, length, length, total, window, signal)
    call amplitude_(flag, length, total, window, signal, result, cos_amp, sin_amp, amp)
    write(*, '(a)') "1st"
    write(*,'(a,e32.16,a,e32.16,a,2e32.16)')" frequency",result," error",abs(result-1.0_rk*frequency)," amplitude",cos_amp,sin_amp
    peak = +2_ik
    result = frequency_(flag, range_min, range_max, peak, method, length, length, total, window, signal)
    call amplitude_(flag, length, total, window, signal, result, cos_amp, sin_amp, amp)
    write(*, '(a)') "2nd"
    write(*,'(a,e32.16,a,e32.16,a,2e32.16)')" frequency",result," error",abs(result-2.0_rk*frequency)," amplitude",cos_amp,sin_amp
    write(*, *)
  end block

  ! filter and estimate frequency
  block
    call filter_(length, signal_r, limit, list)
    call convert_(length, signal_r, signal)
    write(*, '(a)') "filter(+)"
    peak = +1_ik
    result = frequency_(flag, range_min, range_max, peak, method, length, length, total, window, signal)
    call amplitude_(flag, length, total, window, signal, result, cos_amp, sin_amp, amp)
    write(*, '(a)') "1st"
    write(*,'(a,e32.16,a,e32.16,a,2e32.16)')" frequency",result," error",abs(result-1.0_rk*frequency)," amplitude",cos_amp,sin_amp
    peak = +2_ik
    result = frequency_(flag, range_min, range_max, peak, method, length, length, total, window, signal)
    call amplitude_(flag, length, total, window, signal, result, cos_amp, sin_amp, amp)
    write(*, '(a)') "2nd"
    write(*,'(a,e32.16,a,e32.16,a,2e32.16)')" frequency",result," error",abs(result-2.0_rk*frequency)," amplitude",cos_amp,sin_amp
    write(*, *)
  end block

  ! in general application of filter might improve frequency estimation accuracy (statistically)

end program example

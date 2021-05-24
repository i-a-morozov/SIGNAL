! example-07: frequency estimation (fft data memorization)
program example

  use signal

  implicit none

  integer(ik), parameter                 :: length = 2_ik**10_ik  ! input signal length
  integer(ik)                            :: flag                  ! complex flag (0/1)
  real(rk)                               :: range_min             ! (min) frequency range
  real(rk)                               :: range_max             ! (max) frequency range
  integer(ik)                            :: peak                  ! fourie spectra peak id
  integer(ik)                            :: method                ! frequency estimation method
  integer(ik)                            :: order                 ! cosine window order
  real(rk), dimension(length)            :: window                ! cosine window data
  real(rk)                               :: total                 ! sum of window elements
  real(rk)                               :: frequency             ! exact signal frequency
  real(rk), dimension(length)            :: signal_r, signal_i    ! signal real and complex parts
  real(rk), dimension(2_ik*length)       :: signal                ! input signal, [..., sr_i, si_i, ...]
  integer(ik), parameter                 :: limit = 128_ik        ! number of signals
  real(rk), dimension(limit,2_ik*length) :: data                  ! matrix of signals
  real(rk), dimension(limit)             :: list                  ! list of exact frequencies
  real(rk), dimension(limit)             :: output                ! estimated frequencies

  integer(ik)                            :: i, j

  ! set complex flag (0/1 for real/complex signal)
  flag = 1_ik
  range_min = 0.00_rk
  range_max = 0.99_rk

  ! set test data
  do i = 1_ik, limit, 1_ik
    ! set signal exact frequency, real and imaginary parts
    frequency = 0.623456789_rk + real(i, rk)/real(limit, rk)/10.0_rk
    list(i) = frequency
    do j = 1_ik, length, 1_ik
      signal_r(j) = +cos(two_pi*frequency*real(j, rk))
      signal_i(j) = -sin(two_pi*frequency*real(j, rk))
    end do
    ! format test signal
    if(flag == 0_ik) then
      call convert_(length, signal_r, signal)
    else
      call convert_(length, signal_r, signal_i, signal)
    end if
    data(i,:) = signal
  end do

  ! set window and window sum
  order = 2_ik
  call window_(length, order, window)
  total = sum(window)

  ! compute data tables
  call compute_table_(length, length)

  ! estimate frequency
  peak = 0_ik
  method = frequency_parabola
  !$omp parallel do private(signal)
  do i = 1_ik, limit, 1_ik
    output(i) = frequency__(flag, range_min, range_max, peak, method, length, length, total, window, data(i,:))
  end do
  !$omp end parallel do

  ! destroy data tables
  call destroy_table_()

  ! result
  do i = 1_ik, limit, 1_ik
    write(*, '(3e32.16)') list(i), output(i), abs(list(i)-output(i))
  end do

end program example
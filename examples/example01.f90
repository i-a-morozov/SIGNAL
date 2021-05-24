! example-01: frequency estimation (signal with one harmonic)
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

  integer(ik)                      :: i

  ! set complex flag (0/1 for real/complex signal)
  ! flag is used by frequency_ function
  ! for real signals (flag = 0) 1st half of fourie amplitude spetrum is used for frequency estimation
  ! frequency resolution is [0.0,0.5]
  ! for complex signals (flag = 1) full fourie amplitude spectra is used for frequency estimation
  ! frequency resolution is [0.0,1.0]
  flag = 1_ik

  ! set frequency range (0.0 <= range_min <= range_max <= 1.0)
  ! frequency range can be used to target selected range of frequencies
  ! max range should match <flag>, i.e. should be less or equal to largest frequency defined by <flag>
  range_min = 0.0_rk
  range_max = 1.0_rk/(2.0_rk-real(flag, rk))

  ! set signal exact frequency, real and imaginary parts
  frequency = 0.2143658791_rk
  do i = 1_ik, length, 1_ik
    signal_r(i) = +cos(two_pi*frequency*real(i, rk))
    signal_i(i) = -sin(two_pi*frequency*real(i, rk))
  end do

  ! add constant component
  signal_r = signal_r + 0.1_rk

  ! format test signal
  ! signal_r = [..., sr_i, ...]
  ! signal_i = [..., si_i, ...]
  ! signal = [..., sr_i, si_i, ...]
  signal = 0.0_rk
  if(flag == 0_ik) then
    call convert_(length, signal_r, signal)
  else
    call convert_(length, signal_r, signal_i, signal)
  end if

  ! set window and window sum
  order = 2_ik
  call window_(length, order, window)
  total = sum(window)

  ! if signal length is not equal to power of two
  ! signal can be padded with zeros to the nearest power of two length
  ! in this example signal has correct length

  ! round_up_(length) returns next power of two length
  write(*, '(1a,1i10)') "length         : ", length
  write(*, '(1a,1i10)') "next length    : ", round_up_(length)

  ! estimate frequency
  ! signal should have power of two length (real and imaginary part)

  ! peak = 0 uses largest fourier amplitude spectra bin
  ! peak = n uses n'th fourier amplitude spectra peak
  ! fft is defined by __fft__ directive
  ! fft_radix_two_    -- radix-two fft (nrf77), threadsafe
  ! fft_radix_eight_  -- radix-eight fft (takuya ooura), default, threadsafe
  ! fft_external_     -- fftw, not threadsafe

  write(*, '(1a,1i10)') "flag           : ", flag
  write(*, '(1a,1i10)') "order          : ", order

  ! estimate frequency (maximum bin)
  write(*, '(1a)') "(maximum bin)"
  peak = 0_ik
  ! fft max bin or selected peak
  method = frequency_fft
  result = frequency_(flag, range_min, range_max, peak, method, length, length, total, window, signal)
  write(*, '(1a,2e32.16)') "fft            : ", result, abs(result-frequency)
  ! fft max bin or selected peak and zero padding
  method = frequency_fft
  result = frequency_(flag, range_min, range_max, peak, method, length, 2_ik*length, total, window, signal)
  write(*, '(1a,2e32.16)') "fft (padded)   : ", result, abs(result-frequency)
  ! ffrft refined spectra
  method = frequency_ffrft
  result = frequency_(flag, range_min, range_max, peak, method, length, length, total, window, signal)
  write(*, '(1a,2e32.16)') "ffrft          : ", result, abs(result-frequency)
  ! parabola interpolation of ffrft refined spectra
  method = frequency_parabola
  result = frequency_(flag, range_min, range_max, peak, method, length, length, total, window, signal)
  write(*, '(1a,2e32.16)') "parabola       : ", result, abs(result-frequency)
  ! parabola fit of ffrft refined spectra
  ! number of fit points is given by <parabola_fit_length> parameter
  method = frequency_parabola_fit
  result = frequency_(flag, range_min, range_max, peak, method, length, length, total, window, signal)
  write(*, '(1a,2e32.16)') "parabola (fit) : ", result, abs(result-frequency)
  ! maximum search of dtft spectra
  ! search function is defined by __search__ directive
  ! maximum number of iterations is controlled by <search_limit> parameter
  ! search tolerance is controlled by <search_tolerance> parameter
  method = frequency_search
  result = frequency_(flag, range_min, range_max, peak, method, length, length, total, window, signal)
  write(*, '(1a,2e32.16)') "search         : ", result, abs(result-frequency)

  ! if peak \= 0, peak detection is used to locate selected peak
  ! peaks are sorted by __sort__
  ! sort_bubble_      -- bubble sort
  ! sort_quick_       -- quick sort, default
  ! estimate frequency (1st peak)
  write(*, '(1a)') "(1st peak)"
  peak = 1_ik
  method = frequency_fft
  result = frequency_(flag, range_min, range_max, peak, method, length, length, total, window, signal)
  write(*, '(1a,2e32.16)') "fft            : ", result, abs(result-frequency)
  method = frequency_fft
  result = frequency_(flag, range_min, range_max, peak, method, length, 2_ik*length, total, window, signal)
  write(*, '(1a,2e32.16)') "fft (padded)   : ", result, abs(result-frequency)
  method = frequency_ffrft
  result = frequency_(flag, range_min, range_max, peak, method, length, length, total, window, signal)
  write(*, '(1a,2e32.16)') "ffrft          : ", result, abs(result-frequency)
  method = frequency_parabola
  result = frequency_(flag, range_min, range_max, peak, method, length, length, total, window, signal)
  write(*, '(1a,2e32.16)') "parabola       : ", result, abs(result-frequency)
  method = frequency_parabola_fit
  result = frequency_(flag, range_min, range_max, peak, method, length, length, total, window, signal)
  write(*, '(1a,2e32.16)') "parabola (fit) : ", result, abs(result-frequency)
  method = frequency_search
  result = frequency_(flag, range_min, range_max, peak, method, length, length, total, window, signal)
  write(*, '(1a,2e32.16)') "search         : ", result, abs(result-frequency)

end program example
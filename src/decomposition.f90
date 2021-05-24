
#include "signal.inc"

submodule (signal) decomposition
  implicit none
  contains
  ! ############################################################################################################################# !
  ! estimate amplitude for given frequency
  ! (subroutine) amplitude_(<flag>, <length>, <total>, <window>, <sequence>, <frequency>, <cos_amp>, <sin_amp>, <amp>)
  ! <flag>                 -- (in)     complex flag (ik), 0/1 for real/complex input sequence
  ! <length>               -- (in)     sequence length (ik)
  ! <total>                -- (in)     sum(window) (rk)
  ! <window>               -- (in)     window array (rk array of length = <length>)
  ! <sequence>             -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <frequency>            -- (in)     frequency (rk)
  ! <cos_amp>              -- (out)    cos amplitude (rk)
  ! <sin_amp>              -- (out)    sin amplitude (rk)
  ! <amp>                  -- (out)    abs amplitude (rk)
  ! void    amplitude_(int* flag, int* length, double* total, double* window, double* sequence, double* frequency, double* cos_amp, double* sin_amp, double* amp) ;
  module subroutine amplitude_(flag, length, total, window, sequence, frequency, cos_amp, sin_amp, amp) &
    bind(c, name = "amplutude_")
    integer(ik), intent(in):: flag
    integer(ik), intent(in):: length
    real(rk), intent(in) :: total
    real(rk), intent(in), dimension(length) :: window
    real(rk), intent(in), dimension(2_ik*length) :: sequence
    real(rk), intent(in) :: frequency
    real(rk), intent(out) :: cos_amp
    real(rk), intent(out) :: sin_amp
    real(rk), intent(out) :: amp
    real(rk), dimension(2_ik*length) :: local
    real(rk), dimension(length) :: list
    real(rk), dimension(2_ik*length) :: delta
    integer(ik) :: i
    call remove_window_mean_(length, total, window, sequence, local)
    list = two_pi*real([(i, i=1_ik, length, 1_ik)], rk)
    delta(1_ik::2_ik) = +cos(frequency*list)
    delta(2_ik::2_ik) = -sin(frequency*list)
    cos_amp = sum(local(1_ik::2_ik)*delta(1_ik::2_ik)*window)+sum(local(2_ik::2_ik)*delta(2_ik::2_ik)*window)
    sin_amp = sum(local(2_ik::2_ik)*delta(1_ik::2_ik)*window)-sum(local(1_ik::2_ik)*delta(2_ik::2_ik)*window)
    cos_amp = 1.0_rk/total*cos_amp
    sin_amp = 1.0_rk/total*sin_amp
    if (flag == 0_ik) then
      cos_amp = 2.0_rk*cos_amp
      sin_amp = 2.0_rk*sin_amp
    end if
    amp = sqrt(cos_amp**2_ik+sin_amp**2_ik)
  end subroutine amplitude_
  ! ############################################################################################################################# !
  ! signal decomposition
  ! (subroutine) decomposition_(<flag>, <range_min>, <range_max>, <method>, <mode>, <length>, <pad>, <total>, <window>, <sequence>, <loop>, <frequency>, <cos_amp>, <sin_amp>)
  ! <flag>                 -- (in)     complex flag (ik), 0/1 for real/complex input sequence
  ! <range_min>            -- (in)     (min) frequency range (rk)
  ! <range_max>            -- (in)     (max) frequency range (rk)
  ! <method>               -- (in)     frequency approximation method (ik), frequency_fft = 0_ik, frequency_ffrft = 1_ik, frequecy_parabola = 2_ik
  ! <mode>                 -- (in)     decompostion mode (ik), <mode> = decomposition_subtract = 0 or <mode> = decomposition_peak = 1
  ! <length>               -- (in)     sequence length (ik)
  ! <pad>                  -- (in)     padded sequence length (ik), if pad > length, input sequence is padded
  ! <total>                -- (in)     sum(window) (rk)
  ! <window>               -- (in)     window array (rk array of length = <length>)
  ! <sequence>             -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <loop>                 -- (in)     number of iterations/peaks (ik)
  ! <frequency>            -- (out)    frequency array (rk array of length = <loop>)
  ! <cos_amp>              -- (out)    cos amplitude array (rk array of length = <loop>)
  ! <sin_amp>              -- (out)    sin amplitude array (rk array of length = <loop>)
  ! void    decomposition_(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, double* sequence, int* loop, double* frequency, double* cos_amp, double* sin_amp) ;
  module subroutine decomposition_(flag, range_min, range_max, &
    method, mode, length, pad, total, window, sequence, loop, frequency, cos_amp, sin_amp) &
    bind(c, name = "decomposition_")
    integer(ik), intent(in):: flag
    real(rk), intent(in) :: range_min
    real(rk), intent(in) :: range_max
    integer(ik), intent(in):: method
    integer(ik), intent(in):: mode
    integer(ik), intent(in):: length
    integer(ik), intent(in):: pad
    real(rk), intent(in) :: total
    real(rk), intent(in), dimension(length) :: window
    real(rk), intent(in), dimension(2_ik*length) :: sequence
    integer(ik), intent(in) :: loop
    real(rk), dimension(loop), intent(out) :: frequency
    real(rk), dimension(loop), intent(out) :: cos_amp
    real(rk), dimension(loop), intent(out) :: sin_amp
    real(rk), dimension(2_ik*length) :: local
    integer(ik) :: i
    real(rk), dimension(length) :: list
    real(rk), dimension(2_ik*length) :: delta
    integer(ik), dimension(loop) :: ordering
    real(rk), dimension(loop) :: amplitude
    local = sequence
    list = two_pi*real([(i, i=1_ik, length, 1_ik)], rk)
    do i = 1_ik, loop, 1_ik
      if (mode == decomposition_subtract) then
        frequency(i) = frequency_(flag, range_min, range_max, 0_ik, method, length, pad, total, window, local)
      else
        frequency(i) = frequency_(flag, range_min, range_max, i, method, length, pad, total, window, sequence)
      end if
      delta(1_ik::2_ik) = +cos(frequency(i)*list)
      delta(2_ik::2_ik) = -sin(frequency(i)*list)
      cos_amp(i) = sum(local(1_ik::2_ik)*delta(1_ik::2_ik)*window)+sum(local(2_ik::2_ik)*delta(2_ik::2_ik)*window)
      sin_amp(i) = sum(local(2_ik::2_ik)*delta(1_ik::2_ik)*window)-sum(local(1_ik::2_ik)*delta(2_ik::2_ik)*window)
      cos_amp(i) = 1.0_rk/total*cos_amp(i)
      sin_amp(i) = 1.0_rk/total*sin_amp(i)
      local(1_ik::2_ik) = local(1_ik::2_ik)-cos_amp(i)*delta(1_ik::2_ik)+sin_amp(i)*delta(2_ik::2_ik)
      local(2_ik::2_ik) = local(2_ik::2_ik)-cos_amp(i)*delta(2_ik::2_ik)-sin_amp(i)*delta(1_ik::2_ik)
    end do
    if (flag == 0_ik) then
      cos_amp = 2.0_rk*cos_amp
      sin_amp = 2.0_rk*sin_amp
    end if
    amplitude = cos_amp**2_ik+sin_amp**2_ik
    do i = 1_ik, loop, 1_ik
      ordering(i) = __maxloc__(amplitude, 1_ik)
      amplitude(ordering(i)) = 0.0_rk
    end do
    frequency = frequency([ordering])
    cos_amp = cos_amp([ordering])
    sin_amp = sin_amp([ordering])
  end subroutine decomposition_
  ! ############################################################################################################################# !
  ! signal decomposition (memorization)
  ! (subroutine) decomposition__(<flag>, <range_min>, <range_max>, <method>, <mode>, <length>, <pad>, <total>, <window>, <sequence>, <loop>, <frequency>, <cos_amp>, <sin_amp>)
  ! <flag>                 -- (in)     complex flag (ik), 0/1 for real/complex input sequence
  ! <range_min>            -- (in)     (min) frequency range (rk)
  ! <range_max>            -- (in)     (max) frequency range (rk)
  ! <method>               -- (in)     frequency approximation method (ik), frequency_fft = 0_ik, frequency_ffrft = 1_ik, frequecy_parabola = 2_ik
  ! <mode>                 -- (in)     decompostion mode (ik), <mode> = decomposition_subtract = 0 or <mode> = decomposition_peak = 1
  ! <length>               -- (in)     sequence length (ik)
  ! <pad>                  -- (in)     padded sequence length (ik), if pad > length, input sequence is padded
  ! <total>                -- (in)     sum(window) (rk)
  ! <window>               -- (in)     window array (rk array of length = <length>)
  ! <sequence>             -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <loop>                 -- (in)     number of iterations/peaks (ik)
  ! <frequency>            -- (out)    frequency array (rk array of length = <loop>)
  ! <cos_amp>              -- (out)    cos amplitude array (rk array of length = <loop>)
  ! <sin_amp>              -- (out)    sin amplitude array (rk array of length = <loop>)
  ! void    decomposition__(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, double* sequence, int* loop, double* frequency, double* cos_amp, double* sin_amp) ;
  module subroutine decomposition__(flag, range_min, range_max, &
    method, mode, length, pad, total, window, sequence, loop, frequency, cos_amp, sin_amp) &
    bind(c, name = "decomposition__")
    integer(ik), intent(in):: flag
    real(rk), intent(in) :: range_min
    real(rk), intent(in) :: range_max
    integer(ik), intent(in):: method
    integer(ik), intent(in):: mode
    integer(ik), intent(in):: length
    integer(ik), intent(in):: pad
    real(rk), intent(in) :: total
    real(rk), intent(in), dimension(length) :: window
    real(rk), intent(in), dimension(2_ik*length) :: sequence
    integer(ik), intent(in) :: loop
    real(rk), dimension(loop), intent(out) :: frequency
    real(rk), dimension(loop), intent(out) :: cos_amp
    real(rk), dimension(loop), intent(out) :: sin_amp
    real(rk), dimension(2_ik*length) :: local
    integer(ik) :: i
    real(rk), dimension(length) :: list
    real(rk), dimension(2_ik*length) :: delta
    integer(ik), dimension(loop) :: ordering
    real(rk), dimension(loop) :: amplitude
    local = sequence
    list = two_pi*real([(i, i=1_ik, length, 1_ik)], rk)
    do i = 1_ik, loop, 1_ik
      if (mode == decomposition_subtract) then
        frequency(i) = frequency__(flag, range_min, range_max, 0_ik, method, length, pad, total, window, local)
      else
        frequency(i) = frequency__(flag, range_min, range_max, i, method, length, pad, total, window, sequence)
      end if
      delta(1_ik::2_ik) = +cos(frequency(i)*list)
      delta(2_ik::2_ik) = -sin(frequency(i)*list)
      cos_amp(i) = sum(local(1_ik::2_ik)*delta(1_ik::2_ik)*window)+sum(local(2_ik::2_ik)*delta(2_ik::2_ik)*window)
      sin_amp(i) = sum(local(2_ik::2_ik)*delta(1_ik::2_ik)*window)-sum(local(1_ik::2_ik)*delta(2_ik::2_ik)*window)
      cos_amp(i) = 1.0_rk/total*cos_amp(i)
      sin_amp(i) = 1.0_rk/total*sin_amp(i)
      local(1_ik::2_ik) = local(1_ik::2_ik)-cos_amp(i)*delta(1_ik::2_ik)+sin_amp(i)*delta(2_ik::2_ik)
      local(2_ik::2_ik) = local(2_ik::2_ik)-cos_amp(i)*delta(2_ik::2_ik)-sin_amp(i)*delta(1_ik::2_ik)
    end do
    if (flag == 0_ik) then
      cos_amp = 2.0_rk*cos_amp
      sin_amp = 2.0_rk*sin_amp
    end if
    amplitude = cos_amp**2_ik+sin_amp**2_ik
    do i = 1_ik, loop, 1_ik
      ordering(i) = __maxloc__(amplitude, 1_ik)
      amplitude(ordering(i)) = 0.0_rk
    end do
    frequency = frequency([ordering])
    cos_amp = cos_amp([ordering])
    sin_amp = sin_amp([ordering])
  end subroutine decomposition__
  ! ############################################################################################################################# !
  ! frequency list (perform decomposition and return list of frequencies)
  ! (subroutine) frequency_list_(<flag>, <range_min>, <range_max>, <method>, <mode>, <length>, <pad>, <total>, <window>, <sequence>, <loop>, <frequency>, <cos_amp>, <sin_amp>)
  ! <flag>                 -- (in)     complex flag (ik), 0/1 for real/complex input sequence
  ! <range_min>            -- (in)     (min) frequency range (rk)
  ! <range_max>            -- (in)     (max) frequency range (rk)
  ! <method>               -- (in)     frequency approximation method (ik), frequency_fft = 0_ik, frequency_ffrft = 1_ik, frequecy_parabola = 2_ik
  ! <mode>                 -- (in)     decompostion mode (ik), <mode> = decomposition_subtract = 0 or <mode> = decomposition_peak = 1
  ! <length>               -- (in)     sequence length (ik)
  ! <pad>                  -- (in)     padded sequence length (ik), if pad > length, input sequence is padded
  ! <total>                -- (in)     sum(window) (rk)
  ! <window>               -- (in)     window array (rk array of length = <length>)
  ! <sequence>             -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <loop>                 -- (in)     number of iterations/peaks (ik)
  ! <frequency>            -- (out)    frequency array (rk array of length = <loop>)
  ! void    frequency_list_(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, double* sequence, int* loop, double* frequency) ;
  module subroutine frequency_list_(flag, range_min, range_max, &
    method, mode, length, pad, total, window, sequence, loop, frequency) &
    bind(c, name = "frequency_list_")
    integer(ik), intent(in):: flag
    real(rk), intent(in) :: range_min
    real(rk), intent(in) :: range_max
    integer(ik), intent(in):: method
    integer(ik), intent(in):: mode
    integer(ik), intent(in):: length
    integer(ik), intent(in):: pad
    real(rk), intent(in) :: total
    real(rk), intent(in), dimension(length) :: window
    real(rk), intent(in), dimension(2_ik*length) :: sequence
    integer(ik), intent(in) :: loop
    real(rk), dimension(loop), intent(out) :: frequency
    real(rk), dimension(loop) :: cos_amp
    real(rk), dimension(loop) :: sin_amp
    call decomposition_(flag, range_min, range_max, &
    method, mode, length, pad, total, window, sequence, loop, frequency, cos_amp, sin_amp)
  end subroutine frequency_list_
  ! ############################################################################################################################# !
  ! frequency list (perform decomposition and return list of frequencies) (memorization)
  ! (subroutine) frequency_list__(<flag>, <range_min>, <range_max>, <method>, <mode>, <length>, <pad>, <total>, <window>, <sequence>, <loop>, <frequency>, <cos_amp>, <sin_amp>)
  ! <flag>                 -- (in)     complex flag (ik), 0/1 for real/complex input sequence
  ! <range_min>            -- (in)     (min) frequency range (rk)
  ! <range_max>            -- (in)     (max) frequency range (rk)
  ! <method>               -- (in)     frequency approximation method (ik), frequency_fft = 0_ik, frequency_ffrft = 1_ik, frequecy_parabola = 2_ik
  ! <mode>                 -- (in)     decompostion mode (ik), <mode> = decomposition_subtract = 0 or <mode> = decomposition_peak = 1
  ! <length>               -- (in)     sequence length (ik)
  ! <pad>                  -- (in)     padded sequence length (ik), if pad > length, input sequence is padded
  ! <total>                -- (in)     sum(window) (rk)
  ! <window>               -- (in)     window array (rk array of length = <length>)
  ! <sequence>             -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <loop>                 -- (in)     number of iterations/peaks (ik)
  ! <frequency>            -- (out)    frequency array (rk array of length = <loop>)
  ! void    frequency_list__(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, double* sequence, int* loop, double* frequency) ;
  module subroutine frequency_list__(flag, range_min, range_max, &
    method, mode, length, pad, total, window, sequence, loop, frequency) &
    bind(c, name = "frequency_list__")
    integer(ik), intent(in):: flag
    real(rk), intent(in) :: range_min
    real(rk), intent(in) :: range_max
    integer(ik), intent(in):: method
    integer(ik), intent(in):: mode
    integer(ik), intent(in):: length
    integer(ik), intent(in):: pad
    real(rk), intent(in) :: total
    real(rk), intent(in), dimension(length) :: window
    real(rk), intent(in), dimension(2_ik*length) :: sequence
    integer(ik), intent(in) :: loop
    real(rk), dimension(loop), intent(out) :: frequency
    real(rk), dimension(loop) :: cos_amp
    real(rk), dimension(loop) :: sin_amp
    call decomposition__(flag, range_min, range_max, &
    method, mode, length, pad, total, window, sequence, loop, frequency, cos_amp, sin_amp)
  end subroutine frequency_list__
  ! ############################################################################################################################# !
  ! amplitude list (compute amplitudes for list of given frequencies)
  ! (subroutine) amplitude_list_(<flag>, <length>, <total>, <window>, <sequence>, <loop>, <frequency>, <cos_amp>, <sin_amp>)
  ! <flag>                 -- (in)     complex flag (ik), 0/1 for real/complex input sequence
  ! <length>               -- (in)     sequence length (ik)
  ! <total>                -- (in)     sum(window) (rk)
  ! <window>               -- (in)     window array (rk array of length = <length>)
  ! <sequence>             -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <loop>                 -- (in)     number of iterations (ik)
  ! <frequency>            -- (in)     frequency array (rk array of length = <loop>)
  ! <cos_amp>              -- (out)    cos amplitude array (rk array of length = <loop>)
  ! <sin_amp>              -- (out)    sin amplitude array (rk array of length = <loop>)
  ! void    amplitude_list_(int* flag, int* length, double* total, double* window, double* sequence, int* loop, double* frequency, double* cos_amp, double* sin_amp) ;
  module subroutine amplitude_list_(flag, length, total, window, sequence, loop, frequency, cos_amp, sin_amp) &
    bind(c, name = "amplitude_list_")
    integer(ik), intent(in):: flag
    integer(ik), intent(in):: length
    real(rk), intent(in) :: total
    real(rk), intent(in), dimension(length) :: window
    real(rk), intent(in), dimension(2_ik*length) :: sequence
    integer(ik), intent(in) :: loop
    real(rk), dimension(loop), intent(in) :: frequency
    real(rk), dimension(loop), intent(out) :: cos_amp
    real(rk), dimension(loop), intent(out) :: sin_amp
    integer(ik) :: i
    real(rk), dimension(2_ik*length) :: local
    real(rk), dimension(length) :: list
    real(rk), dimension(2_ik*length) :: delta
    call remove_window_mean_(length, total, window, sequence, local)
    list = two_pi*real([(i, i=1_ik, length, 1_ik)], rk)
    do i = 1_ik, loop, 1_ik
      delta(1_ik::2_ik) = +cos(frequency(i)*list)
      delta(2_ik::2_ik) = -sin(frequency(i)*list)
      cos_amp(i) = sum(local(1_ik::2_ik)*delta(1_ik::2_ik)*window)+sum(local(2_ik::2_ik)*delta(2_ik::2_ik)*window)
      sin_amp(i) = sum(local(2_ik::2_ik)*delta(1_ik::2_ik)*window)-sum(local(1_ik::2_ik)*delta(2_ik::2_ik)*window)
      cos_amp(i) = 1.0_rk/real(length, rk)*cos_amp(i)
      sin_amp(i) = 1.0_rk/real(length, rk)*sin_amp(i)
      local(1_ik::2_ik) = local(1_ik::2_ik)-cos_amp(i)*delta(1_ik::2_ik)+sin_amp(i)*delta(2_ik::2_ik)
      local(2_ik::2_ik) = local(2_ik::2_ik)-cos_amp(i)*delta(2_ik::2_ik)-sin_amp(i)*delta(1_ik::2_ik)
    end do
    if (flag == 0_ik) then
      cos_amp = 2.0_rk*cos_amp
      sin_amp = 2.0_rk*sin_amp
    end if
  end subroutine amplitude_list_
  ! ############################################################################################################################# !
  ! frequency correction
  ! (subroutine) frequency_correction_(<flag>, <range_min>, <range_max>, <method>, <mode>, <length>, <pad>, <total>, <window>, <loop>, <frequency>, <cos_amp>, <sin_amp>)
  ! <flag>                 -- (in)     complex flag (ik), 0/1 for real/complex input sequence
  ! <range_min>            -- (in)     (min) frequency range (rk)
  ! <range_max>            -- (in)     (max) frequency range (rk)
  ! <method>               -- (in)     frequency approximation method (ik), frequency_fft = 0_ik, frequency_ffrft = 1_ik, frequecy_parabola = 2_ik
  ! <mode>                 -- (in)     decompostion mode (ik), <mode> = decomposition_subtract = 0 or <mode> = decomposition_peak = 1
  ! <length>               -- (in)     sequence length (ik)
  ! <pad>                  -- (in)     padded sequence length (ik), if pad > length, input sequence is padded
  ! <total>                -- (in)     sum(window) (rk)
  ! <window>               -- (in)     window array (rk array of length = <length>)
  ! <loop>                 -- (in)     number of iterations/peaks (ik)
  ! <frequency>            -- (inout)  frequency array (rk array of length = <loop>)
  ! <cos_amp>              -- (inout)  cos amplitude array (rk array of length = <loop>)
  ! <sin_amp>              -- (inout)  sin amplitude array (rk array of length = <loop>)
  ! void    frequency_correction_(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, int* loop, double* frequency, double* cos_amp, double* sin_amp) ;
  module subroutine frequency_correction_(flag, range_min, range_max, &
    method, mode, length, pad, total, window, loop, frequency, cos_amp, sin_amp) &
    bind(c, name = "frequency_correction_")
    integer(ik), intent(in):: flag
    real(rk), intent(in) :: range_min
    real(rk), intent(in) :: range_max
    integer(ik), intent(in):: method
    integer(ik), intent(in):: mode
    integer(ik), intent(in):: length
    integer(ik), intent(in):: pad
    real(rk), intent(in) :: total
    real(rk), intent(in), dimension(length) :: window
    integer(ik), intent(in) :: loop
    real(rk), dimension(loop), intent(inout) :: frequency
    real(rk), dimension(loop), intent(inout) :: cos_amp
    real(rk), dimension(loop), intent(inout) :: sin_amp
    real(rk), dimension(2_ik*length) :: sequence
    real(rk), dimension(loop) :: f, c, s
    f = frequency
    c = cos_amp
    s = sin_amp
    call generate_signal_(flag, length, sequence, loop, frequency, cos_amp, sin_amp)
    call decomposition_(flag, range_min, range_max, &
    method, mode, length, pad, total, window, sequence, loop, frequency, cos_amp, sin_amp)
    frequency = 2.0_rk*f - frequency
    cos_amp = 2.0_rk*c - cos_amp
    sin_amp = 2.0_rk*s - sin_amp
  end subroutine frequency_correction_
  ! ############################################################################################################################# !
  ! frequency correction
  ! (subroutine) frequency_correction_(<flag>, <range_min>, <range_max>, <method>, <mode>, <length>, <pad>, <total>, <window>, <loop>, <frequency>, <cos_amp>, <sin_amp>)
  ! <flag>                 -- (in)     complex flag (ik), 0/1 for real/complex input sequence
  ! <range_min>            -- (in)     (min) frequency range (rk)
  ! <range_max>            -- (in)     (max) frequency range (rk)
  ! <method>               -- (in)     frequency approximation method (ik), frequency_fft = 0_ik, frequency_ffrft = 1_ik, frequecy_parabola = 2_ik
  ! <mode>                 -- (in)     decompostion mode (ik), <mode> = decomposition_subtract = 0 or <mode> = decomposition_peak = 1
  ! <length>               -- (in)     sequence length (ik)
  ! <pad>                  -- (in)     padded sequence length (ik), if pad > length, input sequence is padded
  ! <total>                -- (in)     sum(window) (rk)
  ! <window>               -- (in)     window array (rk array of length = <length>)
  ! <loop>                 -- (in)     number of iterations/peaks (ik)
  ! <frequency>            -- (inout)  frequency array (rk array of length = <loop>)
  ! <cos_amp>              -- (inout)  cos amplitude array (rk array of length = <loop>)
  ! <sin_amp>              -- (inout)  sin amplitude array (rk array of length = <loop>)
  ! void    frequency_correction__(int* flag, double* range_min, double* range_max, int* method, int* mode, int* length, int* pad, double* total, double* window, int* loop, double* frequency, double* cos_amp, double* sin_amp) ;
  module subroutine frequency_correction__(flag, range_min, range_max, &
    method, mode, length, pad, total, window, loop, frequency, cos_amp, sin_amp) &
    bind(c, name = "frequency_correction__")
    integer(ik), intent(in):: flag
    real(rk), intent(in) :: range_min
    real(rk), intent(in) :: range_max
    integer(ik), intent(in):: method
    integer(ik), intent(in):: mode
    integer(ik), intent(in):: length
    integer(ik), intent(in):: pad
    real(rk), intent(in) :: total
    real(rk), intent(in), dimension(length) :: window
    integer(ik), intent(in) :: loop
    real(rk), dimension(loop), intent(inout) :: frequency
    real(rk), dimension(loop), intent(inout) :: cos_amp
    real(rk), dimension(loop), intent(inout) :: sin_amp
    real(rk), dimension(2_ik*length) :: sequence
    real(rk), dimension(loop) :: f, c, s
    f = frequency
    c = cos_amp
    s = sin_amp
    call generate_signal_(flag, length, sequence, loop, frequency, cos_amp, sin_amp)
    call decomposition__(flag, range_min, range_max, &
    method, mode, length, pad, total, window, sequence, loop, frequency, cos_amp, sin_amp)
    frequency = 2.0_rk*f - frequency
    cos_amp = 2.0_rk*c - cos_amp
    sin_amp = 2.0_rk*s - sin_amp
  end subroutine frequency_correction__
  ! ############################################################################################################################# !
end submodule decomposition
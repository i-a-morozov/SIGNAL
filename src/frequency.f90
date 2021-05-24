
#include "signal.inc"

submodule (signal) frequency
  implicit none
  contains
  ! ############################################################################################################################# !
  ! initial frequency estimation
  ! (function) frequency_initial_(<range_min>, <range_max>, <peak>, <length>, <pad>, <sequence>)
  ! <range_min>            -- (in)     (min) frequency range (rk)
  ! <range_max>            -- (in)     (max) frequency range (rk)
  ! <peak>                 -- (in)     peak number to use (ik), <peak> = 0 use maximum bin, <peak> = n > 0 use n'th peak within given frequency range
  ! <length>               -- (in)     input sequence length (ik)
  ! <pad>                  -- (in)     padded sequence length (ik), if pad > length, input sequence is padded with zeros
  ! <sequence>             -- (in)     input (processed) sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <frequency_>           -- (out)    initial frequency estimation (rk)
  ! double  frequency_initial_(double* range_min, double* range_max, int* peak, int* length, int* pad, double* sequence) ;
  module real(rk) function frequency_initial_(range_min, range_max, peak, length, pad, sequence) &
    bind(c, name = "frequency_initial_")
    real(rk), intent(in) :: range_min
    real(rk), intent(in) :: range_max
    integer(ik), intent(in) :: peak
    integer(ik), intent(in):: length
    integer(ik), intent(in):: pad
    real(rk), intent(in), dimension(2_ik*length) :: sequence
    integer(ik) :: bin_min
    integer(ik) :: bin_max
    real(rk), dimension(2_ik*pad) :: fourier
    integer(ik) :: bin
    real(rk), dimension(pad) :: amplitude
    call pad_(length, pad, sequence, fourier)
    call __fft__(pad, fft_forward, fourier)
    amplitude = log10(epsilon+sqrt(fourier(1_ik::2_ik)**2_ik+fourier(2_ik::2_ik)**2_ik))
    bin_min = int(range_min*real(pad, rk), ik) + 1_ik
    bin_max = int(range_max*real(pad, rk), ik) + 0_ik
    if (peak == 0_ik) then
      bin = bin_min-1_ik+int(__maxloc__(amplitude(bin_min:bin_max:1_ik), 1_ik), ik)
    else
      bin = bin_min-1_ik+peak_(bin_max-bin_min, amplitude(bin_min:bin_max:1_ik), peak)
    end if
    frequency_initial_ = real(bin-1_ik, rk)/real(pad, rk)
  end function frequency_initial_
  ! ############################################################################################################################# !
  ! initial frequency estimation (memorization)
  ! (function) frequency_initial__(<range_min>, <range_max>, <peak>, <length>, <pad>, <sequence>)
  ! <range_min>            -- (in)     (min) frequency range (rk)
  ! <range_max>            -- (in)     (max) frequency range (rk)
  ! <peak>                 -- (in)     peak number to use (ik), <peak> = 0 use maximum bin, <peak> = n > 0 use n'th peak within given frequency range
  ! <length>               -- (in)     input sequence length (ik)
  ! <pad>                  -- (in)     padded sequence length (ik), if pad > length, input sequence is padded with zeros
  ! <sequence>             -- (in)     input (processed) sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <frequency_>           -- (out)    initial frequency estimation (rk)
  ! double  frequency_initial__(double* range_min, double* range_max, int* peak, int* length, int* pad, double* sequence) ;
  module real(rk) function frequency_initial__(range_min, range_max, peak, length, pad, sequence) &
    bind(c, name = "frequency_initial__")
    real(rk), intent(in) :: range_min
    real(rk), intent(in) :: range_max
    integer(ik), intent(in) :: peak
    integer(ik), intent(in):: length
    integer(ik), intent(in):: pad
    real(rk), intent(in), dimension(2_ik*length) :: sequence
    integer(ik) :: bin_min
    integer(ik) :: bin_max
    real(rk), dimension(2_ik*pad) :: fourier
    integer(ik) :: bin
    real(rk), dimension(pad) :: amplitude
    call pad_(length, pad, sequence, fourier)
    call fft_radix_eight__(pad, fft_forward, fourier, bank%bit_fft, bank%trig_fft)
    amplitude = log10(epsilon+sqrt(fourier(1_ik::2_ik)**2_ik+fourier(2_ik::2_ik)**2_ik))
    bin_min = int(range_min*real(pad, rk), ik) + 1_ik
    bin_max = int(range_max*real(pad, rk), ik) + 0_ik
    if (peak == 0_ik) then
      bin = bin_min-1_ik+int(__maxloc__(amplitude(bin_min:bin_max:1_ik), 1_ik), ik)
    else
      bin = bin_min-1_ik+peak_(bin_max-bin_min, amplitude(bin_min:bin_max:1_ik), peak)
    end if
    frequency_initial__ = real(bin-1_ik, rk)/real(pad, rk)
  end function frequency_initial__
  ! ############################################################################################################################# !
  ! refine frequency estimation (ffrft)
  ! (function) frequency_refine_(<method>, <length>, <sequence>, <initial>)
  ! <method>               -- (in)     method
  ! <length>               -- (in)     sequence length (ik)
  ! <sequence>             -- (in)     input (processed) sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <initial>              -- (in)     initial frequency guess (rk)
  ! <frequency_refine_>    -- (out)    refined frequency estimation (rk)
  ! double  frequency_refine_(int* method, int* length, double* sequence, double* initial) ;
  module real(rk) function frequency_refine_(method, length, sequence, initial) &
    bind(c, name = "frequency_refine_")
    integer(ik), intent(in):: method
    integer(ik), intent(in):: length
    real(rk), intent(in), dimension(2_ik*length) :: sequence
    real(rk), intent(in) :: initial
    real(rk), dimension(2_ik*length) :: fourier
    integer(ik) :: fst, cnd
    real(rk) :: factor
    real(rk), dimension(length) :: mul, cos_mul, sin_mul
    integer :: i
    frequency_refine_ = 0.0_rk
    fst = int(real(length, rk)*initial, ik)+1_ik
    factor = two_pi*real(fst-2_ik, rk)/real(length, rk)
    mul = factor*real([(i, i = 0_ik, length-1_ik, 1_ik)], rk)
    cos_mul = cos(mul)
    sin_mul = sin(mul)
    fourier(1_ik::2_ik) = sequence(1_ik::2_ik)*cos_mul-sequence(2_ik::2_ik)*sin_mul
    fourier(2_ik::2_ik) = sequence(1_ik::2_ik)*sin_mul+sequence(2_ik::2_ik)*cos_mul
    call ffrft_(length, 2.0_rk/real(length, rk), fourier)
    mul = log10(sqrt(fourier(1_ik::2_ik)**2_ik+fourier(2_ik::2_ik)**2_ik)+epsilon)
    cnd = int(__maxloc__(mul, 1_ik), ik)
    if (method == frequency_ffrft) then
      frequency_refine_ = (real(fst, rk)-2.0_rk+2.0_rk*(real(cnd, rk)-1.0_rk)/real(length, rk))/real(length, rk)
      return
    end if
    if (method == frequency_parabola) then
      frequency_refine_ = real(cnd, rk)-0.5_rk+(mul(-1_ik+cnd)-mul(cnd))/(mul(-1_ik+cnd)-2.0_rk*mul(cnd)+mul(1_ik+cnd))
      frequency_refine_ = (real(fst, rk)-2.0_rk+2.0_rk*(frequency_refine_-1.0_rk)/real(length, rk))/real(length, rk)
      return
    end if
    if (method == frequency_parabola_fit) then
      block
        real(rk), dimension(2_ik*parabola_fit_length+1_ik) :: x
        real(rk), dimension(2_ik*parabola_fit_length+1_ik) :: y
        real(rk) :: a, b, c
        integer(ik), dimension(2_ik*parabola_fit_length+1_ik) :: index
        index = cnd + [(i, i = -parabola_fit_length, +parabola_fit_length, 1_ik)]
        x = real(index, rk)
        y = mul([index])
        call fit_parabola_(2_ik*parabola_fit_length+1_ik, x, y, a, b, c, frequency_refine_)
        frequency_refine_ = (real(fst, rk)-2.0_rk+2.0_rk*(frequency_refine_-1.0_rk)/real(length, rk))/real(length, rk)
      end block
    end if
  end function frequency_refine_
  ! ############################################################################################################################# !
  ! refine frequency estimation (ffrft) (memorization)
  ! (function) frequency_refine__(<method>, <length>, <sequence>, <initial>)
  ! <method>               -- (in)     method
  ! <length>               -- (in)     sequence length (ik)
  ! <sequence>             -- (in)     input (processed) sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <initial>              -- (in)     initial frequency guess (rk)
  ! <frequency_refine_>    -- (out)    refined frequency estimation (rk)
  ! double  frequency_refine__(int* method, int* length, double* sequence, double* initial) ;
  module real(rk) function frequency_refine__(method, length, sequence, initial) &
    bind(c, name = "frequency_refine__")
    integer(ik), intent(in):: method
    integer(ik), intent(in):: length
    real(rk), intent(in), dimension(2_ik*length) :: sequence
    real(rk), intent(in) :: initial
    real(rk), dimension(2_ik*length) :: fourier
    integer(ik) :: fst, cnd
    real(rk) :: factor
    real(rk), dimension(length) :: mul, cos_mul, sin_mul
    integer :: i
    frequency_refine__ = 0.0_rk
    fst = int(real(length, rk)*initial, ik)+1_ik
    factor = two_pi*real(fst-2_ik, rk)/real(length, rk)
    mul = factor*real([(i, i = 0_ik, length-1_ik, 1_ik)], rk)
    cos_mul = cos(mul)
    sin_mul = sin(mul)
    fourier(1_ik::2_ik) = sequence(1_ik::2_ik)*cos_mul-sequence(2_ik::2_ik)*sin_mul
    fourier(2_ik::2_ik) = sequence(1_ik::2_ik)*sin_mul+sequence(2_ik::2_ik)*cos_mul
    call ffrft__(length, fourier, bank%bit_ffrft, bank%trig_ffrft, bank%cos_fst, bank%sin_fst, bank%cos_lst, bank%sin_lst)
    mul = log10(sqrt(fourier(1_ik::2_ik)**2_ik+fourier(2_ik::2_ik)**2_ik)+epsilon)
    cnd = int(__maxloc__(mul, 1_ik), ik)
    if (method == frequency_ffrft) then
      frequency_refine__ = (real(fst, rk)-2.0_rk+2.0_rk*(real(cnd, rk)-1.0_rk)/real(length, rk))/real(length, rk)
      return
    end if
    if (method == frequency_parabola) then
      frequency_refine__ = real(cnd, rk)-0.5_rk+(mul(-1_ik+cnd)-mul(cnd))/(mul(-1_ik+cnd)-2.0_rk*mul(cnd)+mul(1_ik+cnd))
      frequency_refine__ = (real(fst, rk)-2.0_rk+2.0_rk*(frequency_refine__-1.0_rk)/real(length, rk))/real(length, rk)
      return
    end if
    if (method == frequency_parabola_fit) then
      block
        real(rk), dimension(2_ik*parabola_fit_length+1_ik) :: x
        real(rk), dimension(2_ik*parabola_fit_length+1_ik) :: y
        real(rk) :: a, b, c
        integer(ik), dimension(2_ik*parabola_fit_length+1_ik) :: index
        index = cnd + [(i, i = -parabola_fit_length, +parabola_fit_length, 1_ik)]
        x = real(index, rk)
        y = mul([index])
        call fit_parabola_(2_ik*parabola_fit_length+1_ik, x, y, a, b, c, frequency_refine__)
        frequency_refine__ = (real(fst, rk)-2.0_rk+2.0_rk*(frequency_refine__-1.0_rk)/real(length, rk))/real(length, rk)
      end block
    end if
  end function frequency_refine__
  ! ############################################################################################################################# !
  ! refine frequency estimation (binary search)
  ! (function) binary_amplitude_(<flag>, <length>, <total>, <window>, <sequence>, <initial>)
  ! <flag>                 -- (in)     complex flag (ik), 0/1 for real/complex input sequence
  ! <length>               -- (in)     sequence length (ik)
  ! <total>                -- (in)     sum(window) (rk)
  ! <window>               -- (in)     window array (rk array of length = <length>)
  ! <sequence>             -- (in)     input (unprocessed) sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <initial>              -- (in)     initial frequency guess (rk)
  ! <binary_amplitude_>    -- (out)    refined frequency (rk)
  ! double  binary_amplitude_(int* flag, int* length, double* total, double* window, double* sequence, double* initial) ;
  module real(rk) function binary_amplitude_(flag, length, total, window, sequence, initial) &
    bind(c, name = "binary_amplitude_")
    integer(ik), intent(in):: flag
    integer(ik), intent(in):: length
    real(rk), intent(in) :: total
    real(rk), intent(in), dimension(length) :: window
    real(rk), intent(in), dimension(2_ik*length) :: sequence
    real(rk), intent(in) :: initial
    binary_amplitude_ = binary_(search_, initial, 1.0_rk/real(length, rk), search_limit, search_tolerance)
  contains
    real(rk) function search_(frequency)
      real(rk), intent(in) :: frequency
      real(rk) :: cos_amp
      real(rk) :: sin_amp
      real(rk) :: amplitude
      call amplitude_(flag, length, total, window, sequence, frequency, cos_amp, sin_amp, amplitude)
      search_ = amplitude
    end function search_
  end function
  ! ############################################################################################################################# !
  ! refine frequency estimation (golden search)
  ! (function) golden_amplitude_(<flag>, <length>, <total>, <window>, <sequence>, <initial>)
  ! <flag>                 -- (in)     complex flag (ik), 0/1 for real/complex input sequence
  ! <length>               -- (in)     sequence length (ik)
  ! <total>                -- (in)     sum(window) (rk)
  ! <window>               -- (in)     window array (rk array of length = <length>)
  ! <sequence>             -- (in)     input (unprocessed) sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <initial>              -- (in)     initial frequency guess (rk)
  ! <golden_amplitude_>    -- (out)    refined frequency (rk)
  ! double  golden_amplitude_(int* flag, int* length, double* total, double* window, double* sequence, double* initial) ;
  module real(rk) function golden_amplitude_(flag, length, total, window, sequence, initial) &
    bind(c, name = "golden_amplitude_")
    integer(ik), intent(in):: flag
    integer(ik), intent(in):: length
    real(rk), intent(in) :: total
    real(rk), intent(in), dimension(length) :: window
    real(rk), intent(in), dimension(2_ik*length) :: sequence
    real(rk), intent(in) :: initial
    golden_amplitude_ = golden_(search_, initial, 1.0_rk/real(length, rk), search_limit, search_tolerance)
  contains
    real(rk) function search_(frequency)
      real(rk), intent(in) :: frequency
      real(rk) :: cos_amp
      real(rk) :: sin_amp
      real(rk) :: amplitude
      call amplitude_(flag, length, total, window, sequence, frequency, cos_amp, sin_amp, amplitude)
      search_ = amplitude
    end function search_
  end function
  ! ############################################################################################################################# !
  ! frequency estimation (gereric)
  ! (function) frequency_(<flag>, <range_min>, <range_max>, <peak>, <method>, <length>, <pad>, <total>, <window>, <sequence>)
  ! <flag>                 -- (in)     complex flag (ik), 0/1 for real/complex input sequence
  ! <range_min>            -- (in)     (min) frequency range (rk)
  ! <range_max>            -- (in)     (max) frequency range (rk)
  ! <peak>                 -- (in)     peak number to use (ik), <peak> = 0 use maximum bin, <peak> = n > 0 use n'th peak within given frequency range
  ! <method>               -- (in)     frequency estimation method (ik)
  ! <length>               -- (in)     input sequence length (ik)
  ! <pad>                  -- (in)     padded sequence length (ik), if pad > length, input sequence is padded with zeros
  ! <total>                -- (in)     sum(window) (rk)
  ! <window>               -- (in)     window array (rk array of length = <length>)
  ! <sequence>             -- (in)     input (unprocessed) sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <frequency_>           -- (out)    frequency estimation (rk)
  ! double  frequency_(int* flag, double* range_min, double* range_max, int* peak, int* method, int* length, int* pad, double* total, double* window, double* sequence) ;
  module real(rk) function frequency_(flag, range_min, range_max, peak, method, length, pad, total, window, sequence) &
    bind(c, name = "frequency_")
    integer(ik), intent(in):: flag
    real(rk), intent(in) :: range_min
    real(rk), intent(in) :: range_max
    integer(ik), intent(in) :: peak
    integer(ik), intent(in) :: method
    integer(ik), intent(in):: length
    integer(ik), intent(in):: pad
    real(rk), intent(in) :: total
    real(rk), intent(in), dimension(length) :: window
    real(rk), intent(in), dimension(2_ik*length) :: sequence
    real(rk), dimension(2_ik*length) :: copy
    real(rk), dimension(2_ik*length) :: local
    call remove_window_mean_(length, total, window, sequence, copy)
    call apply_window_(length, window, copy, local)
    frequency_ = frequency_initial_(range_min, range_max, peak, length, pad, local)
    if (method == frequency_fft) return
    if (method == frequency_ffrft .or. method == frequency_parabola .or. method == frequency_parabola_fit) then
      frequency_ = frequency_refine_(method, length, local, frequency_)
      return
    end if
    if (method == frequency_search) then
      frequency_ = __search__(flag, length, total, window, sequence, frequency_)
      return
    end if
  end function frequency_
  ! ############################################################################################################################# !
  ! frequency estimation (gereric) (memorization)
  ! (function) frequency__(<flag>, <range_min>, <range_max>, <peak>, <method>, <length>, <pad>, <total>, <window>, <sequence>)
  ! <flag>                 -- (in)     complex flag (ik), 0/1 for real/complex input sequence
  ! <range_min>            -- (in)     (min) frequency range (rk)
  ! <range_max>            -- (in)     (max) frequency range (rk)
  ! <peak>                 -- (in)     peak number to use (ik), <peak> = 0 use maximum bin, <peak> = n > 0 use n'th peak within given frequency range
  ! <method>               -- (in)     frequency estimation method (ik)
  ! <length>               -- (in)     input sequence length (ik)
  ! <pad>                  -- (in)     padded sequence length (ik), if pad > length, input sequence is padded with zeros
  ! <total>                -- (in)     sum(window) (rk)
  ! <window>               -- (in)     window array (rk array of length = <length>)
  ! <sequence>             -- (in)     input (unprocessed) sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <frequency_>           -- (out)    frequency estimation (rk)
  ! double  frequency__(int* flag, double* range_min, double* range_max, int* peak, int* method, int* length, int* pad, double* total, double* window, double* sequence) ;
  module real(rk) function frequency__(flag, range_min, range_max, peak, method, length, pad, total, window, sequence) &
    bind(c, name = "frequency__")
    integer(ik), intent(in):: flag
    real(rk), intent(in) :: range_min
    real(rk), intent(in) :: range_max
    integer(ik), intent(in) :: peak
    integer(ik), intent(in) :: method
    integer(ik), intent(in):: length
    integer(ik), intent(in):: pad
    real(rk), intent(in) :: total
    real(rk), intent(in), dimension(length) :: window
    real(rk), intent(in), dimension(2_ik*length) :: sequence
    real(rk), dimension(2_ik*length) :: copy
    real(rk), dimension(2_ik*length) :: local
    call remove_window_mean_(length, total, window, sequence, copy)
    call apply_window_(length, window, copy, local)
    frequency__ = frequency_initial__(range_min, range_max, peak, length, pad, local)
    if (method == frequency_fft) return
    if (method == frequency_ffrft .or. method == frequency_parabola .or. method == frequency_parabola_fit) then
      frequency__ = frequency_refine__(method, length, local, frequency__)
      return
    end if
    if (method == frequency_search) then
      frequency__ = __search__(flag, length, total, window, sequence, frequency__)
      return
    end if
  end function frequency__
  ! ############################################################################################################################# !
end submodule frequency

#include "signal.inc"

submodule (signal) auxiliary
  implicit none
  contains
  ! ############################################################################################################################# !
  ! factorial
  ! (function) factorial_(<number>)
  ! <number>               -- (in)     number (ik)
  ! <factorial_>           -- (out)    factorial of <n> (rk)
  module real(rk) function factorial_(number)
    integer(ik), intent(in) :: number
    integer(ik) :: i
    factorial_ = 1.0_rk
    do i = 1_ik, number, 1_ik
      factorial_ = real(i, rk)*factorial_
    end do
  end function factorial_
  ! ############################################################################################################################# !
  ! gamma regularized
  ! (function) gamma_regularized_(<a>, <x>, <y>)
  ! <a>                    -- (in)     a (rk)
  ! <x>                    -- (in)     x (rk)
  ! <y>                    -- (in)     y (rk)
  ! <gamma_regularized_>   -- (out)    gamma regularized of <a>, <x> and <y> (rk)
  module real(rk) function gamma_regularized_(a, x, y)
    real(rk), intent(in) :: a
    real(rk), intent(in) :: x
    real(rk), intent(in) :: y
    gamma_regularized_ = (gamma_(a, x)-gamma_(a, y))/gamma_(a)
  end function gamma_regularized_
  ! ############################################################################################################################# !
  ! minloc
  ! (function) minloc_(<sequence>)
  ! <sequence>             -- (in)     sequence (rk array)
  ! <minloc_>              -- (out)    minimum location (ik)
  module integer(ik) function minloc_(sequence, empty)
    real(rk), dimension(:), contiguous, intent(in) :: sequence
    integer, intent(in) :: empty
    real(rk) :: element
    integer(ik) :: i
    minloc_ = 0_ik
    element = minval(sequence,empty)
    do i = 1_ik, size(sequence, kind = ik), 1_ik
      if (sequence(i) == element) then
        minloc_ = i
        return
      end if
    end do
  end function minloc_
  ! ############################################################################################################################# !
  ! maxloc
  ! (function) maxloc_(<sequence>)
  ! <sequence>             -- (in)     sequence (rk array)
  ! <maxloc_>              -- (out)    maximum location (ik)
  module integer(ik) function maxloc_(sequence, empty)
    real(rk), dimension(:), contiguous, intent(in) :: sequence
    integer, intent(in) :: empty
    real(rk) :: element
    integer(ik) :: i
    maxloc_ = 0_ik
    element = maxval(sequence, empty)
    do i = 1_ik, size(sequence, kind = ik), 1_ik
      if (sequence(i) == element) then
        maxloc_ = i
        return
      end if
    end do
  end function maxloc_
  ! ############################################################################################################################# !
  ! sort (bubble, descending)
  ! (subroutine) sort_bubble_(<length>, <sequence>, <fst>, <lst>)
  ! <length>               -- (in)     sequence length (ik)
  ! <sequence>             -- (in)     (unsorted) sequence (rk array of length = <length>)
  ! <sequence>             -- (out)    (sorted, descending) sequence (rk array of length = <length>)
  module subroutine sort_bubble_(length, sequence, fst, lst)
    integer(ik), intent(in) :: length
    real(rk), dimension(length), intent(inout) :: sequence
    integer(ik), intent(in) :: fst
    integer(ik), intent(in) :: lst
    integer(ik) :: i, j
    logical :: swapped
    if (fst > lst) return
    do i = length-1_ik, 1_ik, -1_ik
      swapped = .false.
      do j = 1_ik, i, 1_ik
        if (sequence(j) < sequence(j+1_ik)) then
          sequence(j:j+1_ik) = sequence([j+1_ik, j])
          swapped = .true.
        end if
      end do
      if (.not. swapped) exit
    end do
  end subroutine sort_bubble_
  ! ############################################################################################################################# !
  ! sort (quick, descending)
  ! (subroutine) sort_quick_(<length>, <sequence>, <fst>, <lst>)
  ! <length>               -- (in)     sequence length (ik)
  ! <sequence>             -- (in)     (unsorted) sequence (rk array of length = <length>)
  ! <sequence>             -- (out)    (sorted, descending) sequence (rk array of length = <length>)
  module recursive subroutine sort_quick_(length, sequence, fst, lst)
    integer(ik), intent(in) :: length
    real(rk), dimension(:), intent(inout) :: sequence
    integer(ik), intent(in) :: fst
    integer(ik), intent(in) :: lst
    real(rk) :: pivot, copy
    integer(ik) :: i, j
    pivot = sequence((fst+lst)/2_ik)
    i = fst
    j = lst
    do
       do while(sequence(i) > pivot)
          i = i+1_ik
       end do
       do while(pivot > sequence(j))
          j = j-1_ik
       end do
       if (i >= j) exit
       copy = sequence(i)
       sequence(i) = sequence(j)
       sequence(j) = copy
       i = i+1_ik
       j = j-1_ik
    end do
    if (fst < i-1_ik) call sort_quick_(length, sequence, fst, i-1_ik)
    if (j+1_ik < lst) call sort_quick_(length, sequence, j+1_ik, lst)
  end subroutine sort_quick_
  ! ############################################################################################################################# !
  ! generate harmonic signal
  ! (subroutine) generate_signal_(<flag>, <length>, <sequence>, <loop>, <frequency>, <cos_amp>, <sin_amp>)
  ! <flag>                 -- (in)     complex flag (ik), 0/1 for real/complex sequence
  ! <length>               -- (in)     sequence length (ik)
  ! <sequence>             -- (out)    input sequence (rk array of length = <length>)
  ! <loop>                 -- (in)     number of harmonics (ik)
  ! <frequency>            -- (in)     frequency array (rk array of length = <loop>)
  ! <cos_amp>              -- (in)     cos amplitude array (rk array of length = <loop>)
  ! <sin_amp>              -- (in)     sin amplitude array (rk array of length = <loop>)
  ! void    generate_signal_(int* flag, int* length, double* sequence, int* loop, double* frequency, double* cos_amp, double* sin_amp) ;
  module subroutine generate_signal_(flag, length, sequence, loop, frequency, cos_amp, sin_amp) &
    bind(c, name = "generate_signal_")
    integer(ik), intent(in) :: flag
    integer(ik), intent(in) :: length
    real(rk), dimension(2_ik*length), intent(out) :: sequence
    integer(ik), intent(in) :: loop
    real(rk), dimension(loop), intent(in) :: frequency
    real(rk), dimension(loop), intent(in) :: cos_amp
    real(rk), dimension(loop), intent(in) :: sin_amp
    real(rk), dimension(length) :: range
    real(rk), dimension(length) :: fc, fs
    integer(ik) :: i
    range = two_pi*real([(i, i = 1_ik, length, 1_ik)], rk)
    sequence = 0.0_rk
    do i = 1_ik, loop, 1_ik
      fc = cos(frequency(i)*range)
      fs = sin(frequency(i)*range)
      if (flag == 0_ik) then
        sequence(1_ik::2_ik) = sequence(1_ik::2_ik)+cos_amp(i)*fc+sin_amp(i)*fs
      else
        sequence(1_ik::2_ik) = sequence(1_ik::2_ik)+cos_amp(i)*fc+sin_amp(i)*fs
        sequence(2_ik::2_ik) = sequence(2_ik::2_ik)+sin_amp(i)*fc-cos_amp(i)*fs
      end if
    end do
  end subroutine generate_signal_
  ! ############################################################################################################################# !
end submodule auxiliary
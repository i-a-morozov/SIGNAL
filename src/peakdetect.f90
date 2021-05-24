#include "signal.inc"

submodule (signal) peakdetect
  implicit none
  contains
  ! ############################################################################################################################# !
  ! peak list
  ! (subroutine) peak_list_(<length>, <sequence>, <peak_list>)
  ! peak_width             -- (global) peak width (ik)
  ! peak_level             -- (global) peak level threshold (rk)
  ! <length>               -- (in)     sequence length (ik)
  ! <sequence>             -- (in)     sequence (rk array of length = <length>)
  ! <peak_list>            -- (out)    peak list (ik array of length = <length>), value of one correspond to peak location
  module subroutine peak_list_(length, sequence, peak_list)
    integer(ik), intent(in) :: length
    real(rk), dimension(length), intent(in) :: sequence
    integer(ik), dimension(length), intent(out) :: peak_list
    integer(ik) :: i
    integer(ik), dimension(length-(1_ik+2_ik*peak_width)) :: fst
    integer(ik), dimension(length-2_ik) :: lst
    real(rk) :: total
    real(rk) :: last
    real(rk), dimension(length-2_ik*peak_width) :: local
    logical :: this, next
    fst = 0_ik
    total = sum(sequence(1_ik:2_ik*peak_width))
    last = 0.0_rk
    local = sequence(1_ik+peak_width:length-peak_width) - peak_level
    do i = 0_ik, length-(1_ik+2_ik*peak_width)-1_ik, 1_ik
      total = total+sequence(i+1_ik+2_ik*peak_width)-last
      last = sequence(i+1_ik)
      local(i+1_ik) = local(i+1_ik)-total/real(1_ik+2_ik*peak_width,rk)
      if (local(i+1_ik) >= 0.0_rk) then
        fst(i+1_ik) = 1_ik
      end if
    end do
    peak_list = [[(0_ik, i = 1_ik, peak_width, 1_ik)], fst, [(0_ik, i = 1_ik, peak_width, 1_ik)]]
    lst = 0_ik
    do i = 1_ik, length-2_ik, 1_ik
      this = int(sign(1.0_rk, sequence(i+1_ik)-sequence(i)), ik) == +1_ik
      next = int(sign(1.0_rk, sequence(i+2_ik)-sequence(i+1_ik)), ik) == -1_ik
      if (this .and. next) then
        lst(i) = 1_ik
      end if
    end do
    peak_list = peak_list*[0_ik, lst, 0_ik]
  end subroutine peak_list_
  ! ############################################################################################################################# !
  ! total number of peaks
  ! (function) peak_count_(<length>, <peak_list>)
  ! <length>               -- (in)     sequence length (ik)
  ! <peak_list>            -- (in)     peak list (ik array of length <length>)
  ! <peak_count_>          -- (out)    total number of peaks (ik)
  module integer(ik) function peak_count_(length, peak_list)
    integer(ik), intent(in) :: length
    integer(ik), dimension(length), intent(in) :: peak_list
    peak_count_ = sum(peak_list)
  end function peak_count_
  ! ############################################################################################################################# !
  ! detect several peaks (list of ordered peak positions)
  ! (subroutine) peak_detect_(<length>, <sequence>, <peak_length>, <peak_ordered_list>)
  ! <length>               -- (in)     sequence length (ik)
  ! <sequence>             -- (in)     sequence (rk array of length = <length>)
  ! <peak_length>          -- (in)     number of peaks to find (ik)
  ! <peak_ordered_list>    -- (out)    peak positions (ik array of length = <peak_length>)
  module subroutine peak_detect_(length, sequence, peak_length, peak_ordered_list)
    integer(ik), intent(in) :: length
    real(rk), dimension(length), intent(in) :: sequence
    integer(ik), intent(in) :: peak_length
    integer(ik), dimension(peak_length), intent(out) :: peak_ordered_list
    integer(ik), dimension(length) :: peak_list
    integer(ik) :: limit
    integer(ik) :: i, j
    integer(ik) :: position
    real(rk), dimension(length) :: local
    peak_list = 0_ik
    call peak_list_(length, sequence, peak_list)
    peak_ordered_list = 1_ik
    limit = peak_count_(length, peak_list)
    position = 1_ik
    do i = 1_ik, length, 1_ik
      if (peak_list(i) == 1_ik) then
        local(position) = sequence(i)
        position = position + 1_ik
      end if
    end do
    call __sort__(limit, local, 1_ik, limit)
    do i = 1_ik, min(peak_length, limit), 1_ik
      do j = 1_ik, length, 1_ik
        if (local(i) == sequence(j)) peak_ordered_list(i) = j
      end do
    end do
  end subroutine peak_detect_
  ! ############################################################################################################################# !
  ! peak (ranked)
  ! (function) peak_(<length>, <sequence>, <peak_id>)
  ! <length>               -- (in)     sequence length (ik)
  ! <sequence>             -- (in)     sequence (rk array of length <length>)
  ! <peak_id>              -- (in)     peak rank (ik)
  ! <peak_>                -- (out)    peak position (ik)
  ! int     peak_(int* length, double* sequence, int* id) ;
  module integer(ik) function peak_(length, sequence, peak_id) &
    bind(c, name = "peak_")
    integer(ik), intent(in) :: length
    real(rk), dimension(length), intent(in) :: sequence
    integer(ik), intent(in) :: peak_id
    integer(ik), dimension(length) :: peak_list
    call peak_detect_(length, sequence, length, peak_list)
    peak_ = peak_list(peak_id)
  end function peak_
  ! ############################################################################################################################# !
end submodule peakdetect

#include "signal.inc"

submodule (signal) process
  implicit none
  contains
  ! ############################################################################################################################# !
  ! convert input sequence (real)
  ! (subroutine) convert_real_(<length>, <r_part>, <sequence>)
  ! <length>               -- (in)     length (ik)
  ! <r_part>               -- (in)     input sequence r-part (rk array of length = <length>)
  ! <sequence>             -- (out)    sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...] and si_i=0.0_rk for all i
  ! void    convert_real_(int* length, double* r_part, double* sequence) ;
  module subroutine convert_real_(length, r_part, sequence) &
    bind(c, name = "convert_real_")
    integer(ik), intent(in) :: length
    real(rk), intent(in), dimension(length) :: r_part
    real(rk), intent(out), dimension(2_ik*length) :: sequence
    sequence = 0.0_rk
    sequence(1_ik::2_ik) = r_part
  end subroutine convert_real_
  ! ############################################################################################################################# !
  ! convert input sequence (complex)
  ! (subroutine) convert_complex_(<length>, <r_part>, <i_part>, <sequence>)
  ! <length>               -- (in)     length (ik)
  ! <r_part>               -- (in)     input sequence r-part (rk array of length = <length>)
  ! <i_part>               -- (in)     input sequence i-part (rk array of length = <length>)
  ! <sequence>             -- (out)    sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! void    convert_complex_(int* length, double* r_part, double* i_part, double* sequence) ;
  module subroutine convert_complex_(length, r_part, i_part, sequence) &
    bind(c, name = "convert_complex_")
    integer(ik), intent(in) :: length
    real(rk), intent(in), dimension(length) :: r_part
    real(rk), intent(in), dimension(length) :: i_part
    real(rk), intent(out), dimension(2_ik*length) :: sequence
    sequence = 0.0_rk
    sequence(1_ik::2_ik) = r_part
    sequence(2_ik::2_ik) = i_part
  end subroutine convert_complex_
  ! ############################################################################################################################# !
  ! round up (round up to the next power of two)
  ! (function) round_up_(<number>)
  ! <number>               -- (in)     number (ik)
  ! <round_up>             -- (out)    next power of two number (ik)
  ! int     round_up_(int* number) ;
  module integer(ik) function round_up_(number) &
    bind(c, name = "round_up_")
    integer(ik), intent(in) :: number
    round_up_ = 2_ik**ceiling(log(real(number, rk))/log(2.0_rk), kind=ik)
  end function round_up_
  ! ############################################################################################################################# !
  ! zero padding
  ! (subroutine) pad_(<li>, <lo>, <input>, <output>)
  ! <li>                   -- (in)     input sequence length (ik)
  ! <lo>                   -- (in)     output sequence length (ik)
  ! <input>                -- (in)     input sequence (rk) of length = 2*<li>
  ! <output>               -- (in)     padded sequence (rk) of length = 2*<lo>
  ! void    pad_(int* linput, int* loutput, double* input, double* output) ;
  module subroutine pad_(li, lo, input, output) &
    bind(c, name = "pad_")
    integer(ik), intent(in) :: li
    integer(ik), intent(in) :: lo
    real(rk), dimension(2_ik*li), intent(in) :: input
    real(rk), dimension(2_ik*lo), intent(out) :: output
    real(rk) :: lzero(floor(real(lo-li,rk)/2_rk))
    real(rk) :: rzero(ceiling(real(lo-li,rk)/2_rk))
    integer(ik) :: i
    lzero = [(0.0_rk, i = 1_ik, size(lzero), 1_ik)]
    rzero = [(0.0_rk, i = 1_ik, size(rzero), 1_ik)]
    output(1_ik::2_ik) = [lzero, input(1_ik::2_ik), rzero]
    output(2_ik::2_ik) = [lzero, input(2_ik::2_ik), rzero]
  end  subroutine pad_
  ! ############################################################################################################################# !
  ! remove mean
  ! (subroutine) remove_mean_(<length>, <input>, <output> )
  ! <length>               -- (in)     input sequence length (ik)
  ! <input>                -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <output>               -- (out)    output sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! void    remove_mean_(int* length, double* input, double* output) ;
  module subroutine remove_mean_(length, input, output) &
    bind(c, name = "remove_mean_")
    integer(ik), intent(in) :: length
    real(rk), dimension(2_ik*length), intent(in) :: input
    real(rk), dimension(2_ik*length), intent(out) :: output
    output(1_ik::2_ik)  = input(1_ik::2_ik)-sum(input(1_ik::2_ik))/real(length, rk)
    output(2_ik::2_ik)  = input(2_ik::2_ik)-sum(input(2_ik::2_ik))/real(length, rk)
  end subroutine remove_mean_
  ! ############################################################################################################################# !
  ! remove window mean
  ! (subroutine) remove_window_mean_(<length>, <total>, <window>, <input>, <output> )
  ! <length>               -- (in)     input sequence length (ik)
  ! <total>                -- (in)     sum(window) (rk)
  ! <window>               -- (in)     window array (rk array of length = <length>)
  ! <input>                -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <output>               -- (out)    output sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! void    remove_window_mean_(int* length, double* total, double* window, double* input, double* output) ;
  module subroutine remove_window_mean_(length, total, window, input, output) &
    bind(c, name = "remove_window_mean_")
    integer(ik), intent(in) :: length
    real(rk), intent(in) :: total
    real(rk), dimension(length), intent(in) :: window
    real(rk), dimension(2_ik*length), intent(in) :: input
    real(rk), dimension(2_ik*length), intent(out) :: output
    output(1_ik::2_ik) = (input(1_ik::2_ik)-sum(window*input(1_ik::2_ik))/total)
    output(2_ik::2_ik) = (input(2_ik::2_ik)-sum(window*input(2_ik::2_ik))/total)
  end subroutine remove_window_mean_
  ! ############################################################################################################################# !
  ! apply window
  ! (subroutine) apply_window_(<length>, <window>, <input>, <output> )
  ! <length>               -- (in)     input sequence length (ik)
  ! <window>               -- (in)     window array (rk array of length = <length>)
  ! <input>                -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <output>               -- (out)    output sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! void    apply_window_(int* length, double* window, double* input, double* output) ;
  module subroutine apply_window_(length, window, input, output) &
    bind(c, name = "apply_window_")
    integer(ik), intent(in) :: length
    real(rk), dimension(length), intent(in) :: window
    real(rk), dimension(2_ik*length), intent(in) :: input
    real(rk), dimension(2_ik*length), intent(out) :: output
    output(1_ik::2_ik) = input(1_ik::2_ik)*window
    output(2_ik::2_ik) = input(2_ik::2_ik)*window
  end subroutine apply_window_
  ! ############################################################################################################################# !
  ! matrix (generate matrix from sequence)
  ! (subroutine) matrix_(<length>, <sequence>, <matrix>)
  ! <length>               -- (in)     input sequence length (ik)
  ! <sequence>             -- (in)     input sequence (rk)
  ! <matrix>               -- (out)    matrix (<length>/2+1, <length>/2) (rk)
  module subroutine matrix_(length, sequence, matrix)
    integer(ik), intent(in) :: length
    real(rk), dimension(length), intent(in) :: sequence
    real(rk), dimension(length/2_ik+1_ik, length/2_ik), intent(out) :: matrix
    integer(ik) :: i
    do i = 1_ik, length/2_ik+1_ik, 1_ik
      matrix(i, :) = sequence(i:i-1_ik+length/2_ik)
    end do
  end subroutine matrix_
  ! ############################################################################################################################# !
  ! sequence (row) (generate sequence from matrix using 1st and last rows)
  ! (subroutine) sequence_row_(<length>, <sequence>, <matrix>)
  ! <length>               -- (in)     input sequence length (ik)
  ! <sequence>             -- (out)    input sequence (rk)
  ! <matrix>               -- (in)     matrix (<length>/2+1, <length>/2) (rk)
  module subroutine sequence_row_(length, sequence, matrix)
    integer(ik), intent(in) :: length
    real(rk), dimension(length), intent(out) :: sequence
    real(rk), dimension(length/2_ik+1_ik, length/2_ik), intent(in) :: matrix
    sequence = [matrix(1_ik, :), matrix(length/2_ik+1_ik, :)]
  end subroutine sequence_row_
  ! ############################################################################################################################# !
  ! sequence (sum) (generate sequence from matrix using sums of skew diagonals)
  ! (subroutine) sequence_row_(<length>, <sequence>, <matrix>)
  ! <length>               -- (in)     input sequence length (ik)
  ! <sequence>             -- (out)    input sequence (rk)
  ! <matrix>               -- (in)     matrix (<length>/2+1, <length>/2) (rk)
  module subroutine sequence_sum_(length, sequence, matrix)
    integer(ik), intent(in) :: length
    real(rk), dimension(length), intent(out) :: sequence
    real(rk), dimension(length/2_ik+1_ik, length/2_ik), intent(in) :: matrix
    integer(ik) :: row
    integer(ik) :: col
    real(rk), dimension((length/2_ik+1_ik)*(length/2_ik)) :: array
    integer(ik) :: shift
    integer(ik) :: start
    integer(ik) :: count
    integer(ik) :: q, p
    row = length/2_ik+1_ik
    col = length/2_ik
    array = reshape(matrix, shape(array))
    shift = 1_ik ;
    sequence = 0.0_rk
    do q = 1_ik, 2_ik*col, 1_ik
      start = q
      if (q > row) then
        start = q + col*shift
        shift = shift + 1_ik
      end if
      count = 0_ik
      do p = 0_ik, min(q, col)-shift, 1_ik
        sequence(q) = sequence(q)+array(start+p*col)
        count = count + 1_ik
      end do
      sequence(q) = sequence(q)/real(count, rk) ;
    end do
  end subroutine sequence_sum_
  ! ############################################################################################################################# !
  ! filter
  ! (subroutine) filter(<length>, <sequence>, <limit>)
  ! <length>               -- (in)     length (ik)
  ! <sequence>             -- (inout)  sequence (rk array of length = <length>)
  ! <limit>                -- (in)     number of singular values to keep (ik)
  ! <svd_list>             -- (out)    list of singular values
  ! void    filter_(int* length, double* sequence, int* limit, double* svd_list) ;
  module subroutine filter_(length, sequence, limit, svd_list) &
    bind(c, name = "filter_")
    integer(ik), intent(in) :: length
    real(rk), dimension(length), intent(inout) :: sequence
    integer(ik), intent(in) :: limit
    real(rk), dimension(limit), intent(out) :: svd_list
    real(rk), dimension(length/2_ik+1_ik, length/2_ik) :: matrix
    real(rk), dimension(length/2_ik+1_ik, limit) :: u_matrix
    real(rk), dimension(length/2_ik, limit) :: v_matrix
    real(rk), dimension(limit, limit) :: diag
    real(rk), dimension(limit, length/2_ik) :: copy
    integer(ik) :: i
    integer :: nr, nc, ns
    call matrix_(length, sequence, matrix)
    call svd_truncated_(length/2_ik+1_ik, length/2_ik, limit, matrix, svd_list, v_matrix, u_matrix)
    matrix = 0.0_rk
    do i = 1_ik, limit, 1_ik
      diag(i, i) = svd_list(i)
    end do
    nr = int(length)/2_ik + 1_ik
    nc = int(length)/2_ik
    ns = limit
    call dgemm('n','t',ns,nc,ns,1.0_rk,diag,size(diag,1_ik),v_matrix,size(v_matrix,1_ik),0.0_rk,copy,size(copy,1_ik))
    call dgemm('n','n',nr,nc,ns,1.0_rk,u_matrix,size(u_matrix,1_ik),copy,size(copy,1_ik),0.0_rk,matrix,size(matrix,1_ik))
    call __sequence__(length, sequence, matrix)
  end subroutine filter_
  ! ############################################################################################################################# !
end submodule process
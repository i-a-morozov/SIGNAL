
#include "signal.inc"

submodule (signal) transformation
  implicit none
  contains
  ! ############################################################################################################################# !
  ! (linear) fractional complex discrete fourier transform
  ! (subroutine) ffrft_(<length>, <argument>, <sequence>)
  ! <length>               -- (in)     length (ik)
  ! <argument>             -- (in)     parameter (rk)
  ! <sequence>             -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <sequence>             -- (out)    fcdft (rk array of length = 2_ik*<length>), <sequence> = [..., fr_i, fi_i, ...]
  ! void    ffrft_(int* length, double* argument, double* sequence) ;
  module subroutine ffrft_(length, argument, sequence) &
    bind(c, name = "ffrft_")
    integer(ik), intent(in) :: length
    real(rk), intent(in) :: argument
    real(rk), dimension(2_ik*length), intent(inout) :: sequence
    integer(ik) :: i
    real(rk) :: factor
    real(rk), dimension(length)   :: mul, cos_mul, sin_mul
    real(rk), dimension(4_ik*length) :: one, two, tre
    real(rk), dimension(2_ik*length) :: copy
    factor = argument*one_pi/real(length, rk)
    mul = factor*real([(i, i = 0_ik, length-1_ik, 1_ik)], rk)**2_ik
    cos_mul = cos(mul)
    sin_mul = sin(mul)
    one = 0.0_rk
    one(1_ik:2_ik*length:2_ik) = sequence(1_ik::2_ik)*cos_mul-sequence(2_ik::2_ik)*sin_mul
    one(2_ik:2_ik*length:2_ik) = sequence(1_ik::2_ik)*sin_mul+sequence(2_ik::2_ik)*cos_mul
    two = 0.0_rk
    two(1_ik:2_ik*length:2_ik) = +cos_mul
    two(2_ik:2_ik*length:2_ik) = -sin_mul
    mul = -factor*(real([(i, i = length+1_ik, 2_ik*length, 1_ik)], rk)-1.0_rk-2.0_rk*real(length, rk))**2_ik
    two(2_ik*length+1_ik:4_ik*length:2_ik) = cos(mul)
    two(2_ik*length+2_ik:4_ik*length:2_ik) = sin(mul)
    call __fft__(2_ik*length, fft_forward, one)
    call __fft__(2_ik*length, fft_forward, two)
    tre = one
    one(1_ik::2_ik) = tre(1_ik::2_ik)*two(1_ik::2_ik)-tre(2_ik::2_ik)*two(2_ik::2_ik)
    one(2_ik::2_ik) = tre(1_ik::2_ik)*two(2_ik::2_ik)+tre(2_ik::2_ik)*two(1_ik::2_ik)
    call __fft__(2_ik*length, fft_inverse, one)
    copy = 1.0_rk/real(2_ik*length, rk)*one(1_ik:2_ik*length:1_ik)
    sequence(1_ik::2_ik) = copy(1_ik::2_ik)*cos_mul-copy(2_ik::2_ik)*sin_mul
    sequence(2_ik::2_ik) = copy(1_ik::2_ik)*sin_mul+copy(2_ik::2_ik)*cos_mul
  end subroutine ffrft_
  ! ############################################################################################################################# !
  ! (linear) fractional complex discrete fourier transform (memorization)
  ! (subroutine) ffrft__(<length>, <sequence>, <ip>, <work>, <cos_fst>, <sin_fst>, <cos_lst>, <sin_lst>)
  ! <length>               -- (in)     length (ik)
  ! <sequence>             -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <ip>                   -- (in)     ffrft bit data
  ! <work>                 -- (in)     ffrft trig data
  ! <cos_fst>              -- (in)     first cos array
  ! <sin_fst>              -- (in)     first sin array
  ! <cos_lst>              -- (in)     last cos array
  ! <sin_lat>              -- (in)     last sin array
  ! <sequence>             -- (out)    fcdft (rk array of length = 2_ik*<length>), <sequence> = [..., fr_i, fi_i, ...]
  module subroutine ffrft__(length, sequence, ip, work, cos_fst, sin_fst, cos_lst, sin_lst)
    integer(ik), intent(in) :: length
    real(rk), dimension(2_ik*length), intent(inout) :: sequence
    integer(ik), dimension(0_ik : 1_ik+int(sqrt(real(length, rk)), ik)), intent(in) :: ip
    real(rk), dimension(0_ik : length - 1_ik), intent(in) :: work
    real(rk), dimension(length), intent(in)   :: cos_fst
    real(rk), dimension(length), intent(in)   :: sin_fst
    real(rk), dimension(length), intent(in)   :: cos_lst
    real(rk), dimension(length), intent(in)   :: sin_lst
    real(rk), dimension(4_ik*length) :: one, two, tre
    real(rk), dimension(2_ik*length) :: copy
    one = 0.0_rk
    one(1_ik:2_ik*length:2_ik) = sequence(1_ik::2_ik)*cos_fst-sequence(2_ik::2_ik)*sin_fst
    one(2_ik:2_ik*length:2_ik) = sequence(1_ik::2_ik)*sin_fst+sequence(2_ik::2_ik)*cos_fst
    two = 0.0_rk
    two(1_ik:2_ik*length:2_ik) = +cos_fst
    two(2_ik:2_ik*length:2_ik) = -sin_fst
    two(2_ik*length+1_ik:4_ik*length:2_ik) = cos_lst
    two(2_ik*length+2_ik:4_ik*length:2_ik) = sin_lst
    call fft_radix_eight__(2_ik*length, fft_forward, one, ip, work)
    call fft_radix_eight__(2_ik*length, fft_forward, two, ip, work)
    tre = one
    one(1_ik::2_ik) = tre(1_ik::2_ik)*two(1_ik::2_ik)-tre(2_ik::2_ik)*two(2_ik::2_ik)
    one(2_ik::2_ik) = tre(1_ik::2_ik)*two(2_ik::2_ik)+tre(2_ik::2_ik)*two(1_ik::2_ik)
    call fft_radix_eight__(2_ik*length, fft_inverse, one, ip, work)
    copy = 1.0_rk/real(2_ik*length, rk)*one(1_ik:2_ik*length:1_ik)
    sequence(1_ik::2_ik) = copy(1_ik::2_ik)*cos_fst-copy(2_ik::2_ik)*sin_fst
    sequence(2_ik::2_ik) = copy(1_ik::2_ik)*sin_fst+copy(2_ik::2_ik)*cos_fst
  end subroutine ffrft__
  ! ############################################################################################################################# !
  ! (fftw) complex discrete fourier transform
  ! (subroutine) fft_external_(<length>, <direction>, <sequence>)
  ! <length>               -- (in)     length (ik)
  ! <direction>            -- (in)     direction (ik), fft_forward = +1_ik or fft_inverse = -1_ik
  ! <sequence>             -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <sequence>             -- (out)    cdft (rk array of length = 2_ik*<length>), <sequence> = [..., fr_i, fi_i, ...]
  ! void    fft_external_(int* length, int* direction, double* sequence) ;
  module subroutine fft_external_(length, direction, sequence) &
    bind(c, name = "fft_external_")
    integer(ik), intent(in) :: length
    integer(ik), intent(in) :: direction
    real(rk), dimension(2_ik*length), intent(inout) :: sequence
    complex(rk), dimension(length) :: in, out
    integer(ik) :: plan
    in%re = sequence(1_ik::2_ik)
    in%im = sequence(2_ik::2_ik)
    call dfftw_plan_dft_1d(plan, length, in, out, direction, 64_ik)
    call dfftw_execute_dft(plan, in, out)
    call dfftw_destroy_plan(plan)
    sequence(1_ik::2_ik) = out%re
    sequence(2_ik::2_ik) = out%im
  end subroutine fft_external_
  ! ############################################################################################################################# !
  ! (nrf77) complex discrete fourier transform
  ! (subroutine) fft_radix_two_(<length>, <direction>, <sequence>)
  ! <length>               -- (in)     length (ik)
  ! <direction>            -- (in)     direction (ik), fft_forward = +1_ik or fft_inverse = -1_ik
  ! <sequence>             -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <sequence>             -- (out)    cdft (rk array of length = 2_ik*<length>), <sequence> = [..., fr_i, fi_i, ...]
  ! void    fft_radix_two_(int* length, int* direction, double* sequence) ;
  module subroutine fft_radix_two_(length, direction, sequence) &
    bind(c, name = "fft_radix_two_")
    integer(ik), intent(in) :: length
    integer(ik), intent(in) :: direction
    real(rk), dimension(2_ik*length), intent(inout) :: sequence
    integer(ik) :: i, j
    integer(ik) :: n, m
    integer(ik) :: limit, step
    real(rk) :: pim, pre
    real(rk) :: ang
    real(rk) :: wr, wi, wpr, wpi, wcopy
    real(rk) :: factor
    n = 2_ik*length
    j = 1_ik
    factor = real(direction, rk)*two_pi
    do i = 1_ik, n, 2_ik
      if (j > i) then
        pre = sequence(j)
        pim = sequence(j+1_ik)
        sequence(j) = sequence(i)
        sequence(j+1_ik) = sequence(i+1_ik)
        sequence(i) = pre
        sequence(i+1_ik) = pim
      end if
      m = n/2_ik
      do
        if ((m > 2_ik) .and. (j > m)) then
          j = j-m
          m = m/2_ik
        else
          exit
        end if
      end do
      j = j+m
    end do
    limit = 2_ik
    do
      if (n > limit) then
        step = 2_ik*limit
        ang = factor/real(limit, rk)
        wpi = sin(ang)
        ang = sin(0.5_rk*ang)
        wpr = -2.0_rk*ang*ang
        wr = 1.0_rk
        wi = 0.0_rk
        do m = 1_ik, limit, 2_ik
          do i = m, n, step
            j = i+limit
            pre = wr*sequence(j)-wi*sequence(j+1_ik)
            pim = wr*sequence(j+1_ik)+wi*sequence(j)
            sequence(j) = sequence(i)-pre
            sequence(j+1_ik) = sequence(i+1_ik)-pim
            sequence(i) = sequence(i)+pre
            sequence(i+1_ik) = sequence(i+1_ik)+pim
          end do
          wcopy = wr
          wr = wr*wpr-wi*wpi+wr
          wi = wi*wpr+wcopy*wpi+wi
        end do
        limit = step
      else
        exit
      end if
    end do
  end subroutine fft_radix_two_
  ! ############################################################################################################################# !
  ! (takuya ooura) complex discrete fourier transform
  ! (subroutine) fft_radix_eight_(<length>, <direction>, <sequence>)
  ! <length>               -- (in)     length (ik)
  ! <direction>            -- (in)     direction (ik), fft_forward = +1_ik or fft_inverse = -1_ik
  ! <sequence>             -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <sequence>             -- (out)    cdft (rk array of length = 2_ik*<length>), <sequence> = [..., fr_i, fi_i, ...]
  ! void    fft_radix_eight_(int* length, int* direction, double* sequence) ;
  module subroutine fft_radix_eight_(length, direction, sequence) &
    bind(c, name = "fft_radix_eight_")
    integer(ik), intent(in) :: length
    integer(ik), intent(in) :: direction
    real(rk), dimension(2_ik*length), intent(inout) :: sequence
    integer(ik), dimension(2_ik+int(sqrt(real(length/2_ik, rk)), ik)) :: ip
    real(rk), dimension(length/2_ik) :: work
    ip = 0_ik
    work = 0.0_rk
    call cdft_(2_ik*length, direction, sequence, ip, work)
  end subroutine fft_radix_eight_
  ! ############################################################################################################################# !
  ! (takuya ooura) complex discrete fourier transform
  ! (subroutine) fft_radix_eight__(<length>, <direction>, <sequence>, <ip>, <work>)
  ! <length>               -- (in)     length (ik)
  ! <direction>            -- (in)     direction (ik), fft_forward = +1_ik or fft_inverse = -1_ik
  ! <sequence>             -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <sequence>             -- (out)    cdft (rk array of length = 2_ik*<length>), <sequence> = [..., fr_i, fi_i, ...]
  ! <ip>                   -- (in)     ffrft bit data
  ! <work>                 -- (in)     ffrft trig data
  module subroutine fft_radix_eight__(length, direction, sequence, ip, work)
    integer(ik), intent(in) :: length
    integer(ik), intent(in) :: direction
    real(rk), dimension(2_ik*length), intent(inout) :: sequence
    integer(ik), dimension(0_ik : 1_ik+int(sqrt(real(length/2_ik, rk)), ik)), intent(in) :: ip
    real(rk), dimension(0_ik : length/2_ik - 1_ik), intent(in) :: work
    integer(ik), dimension(0_ik : 1_ik+int(sqrt(real(length/2_ik, rk)), ik)) :: ip_copy
    real(rk), dimension(0_ik : length/2_ik - 1_ik) :: work_copy
    ip_copy = ip
    work_copy = work
    call cdft__(2_ik*length, direction, sequence, ip_copy, work_copy)
  end subroutine fft_radix_eight__
  ! ############################################################################################################################# !
  ! compute data table (memorization)
  ! (subroutine) compute_table_(<length>, <pad>)
  ! <length>               -- (in)     length (ik)
  ! <pad>                  -- (in)     padded length (ik)
  ! void    compute_table_(int* length, int* pad) ;
  module subroutine compute_table_(length, pad) &
    bind(c, name = "compute_table_")
    integer(ik), intent(in) :: length
    integer(ik), intent(in) :: pad
    allocate(bank%bit_fft(2_ik+int(sqrt(real(pad/2_ik, rk)), ik)))
    allocate(bank%bit_ffrft(2_ik+int(sqrt(real(length/1_ik, rk)), ik)))
    allocate(bank%trig_fft(pad/2_ik))
    allocate(bank%trig_ffrft(length/1_ik))
    allocate(bank%cos_fst(length))
    allocate(bank%sin_fst(length))
    allocate(bank%cos_lst(length))
    allocate(bank%sin_lst(length))
    call make_fft_data__(pad, bank%bit_fft, bank%trig_fft)
    call make_fft_data__(2_ik*length, bank%bit_ffrft, bank%trig_ffrft)
    call make_ffrft_data__(length, bank%cos_fst, bank%sin_fst, bank%cos_lst, bank%sin_lst)
  end subroutine compute_table_
  ! ############################################################################################################################# !
  ! destroy data table
  ! (subroutine) destroy_table_()
  ! void    destroy_table_() ;
  module subroutine destroy_table_() &
    bind(c, name = "destroy_table_")
    deallocate(bank%bit_fft)
    deallocate(bank%bit_ffrft)
    deallocate(bank%trig_fft)
    deallocate(bank%trig_ffrft)
    deallocate(bank%cos_fst)
    deallocate(bank%sin_fst)
    deallocate(bank%cos_lst)
    deallocate(bank%sin_lst)
  end subroutine destroy_table_
  ! ############################################################################################################################# !
  ! make fft data
  ! (subroutine) make_fft_data__(<length>, <ip>, <work>)
  module subroutine make_fft_data__(length, ip, work)
    integer(ik), intent(in) :: length
    integer(ik), dimension(0_ik : 1_ik+int(sqrt(real(length/2_ik, rk)), ik)), intent(out) :: ip
    real(rk), dimension(0_ik : length/2_ik - 1_ik), intent(out) :: work
    ip = 0_ik
    work = 0.0_rk
    call make_fft_table_(length/2_ik, ip, work)
  end subroutine make_fft_data__
  ! ############################################################################################################################# !
  ! make ffrft data
  ! (subroutine) make_ffrft_data__(length, cos_fst, sin_fst, cos_lst, sin_lst)
  module subroutine make_ffrft_data__(length, cos_fst, sin_fst, cos_lst, sin_lst)
    integer(ik), intent(in) :: length
    real(rk), dimension(length), intent(out)   :: cos_fst
    real(rk), dimension(length), intent(out)   :: sin_fst
    real(rk), dimension(length), intent(out)   :: cos_lst
    real(rk), dimension(length), intent(out)   :: sin_lst
    integer(ik) :: i
    real(rk) :: factor
    real(rk), dimension(length) :: mul
    factor = two_pi/real(length, rk)**2_ik
    mul = factor*real([(i, i = 0_ik, length-1_ik, 1_ik)], rk)**2_ik
    cos_fst = cos(mul)
    sin_fst = sin(mul)
    mul = -factor*(real([(i, i = length+1_ik, 2_ik*length, 1_ik)], rk)-1.0_rk-2.0_rk*real(length, rk))**2_ik
    cos_lst = cos(mul)
    sin_lst = sin(mul)
  end subroutine make_ffrft_data__
  ! ############################################################################################################################# !
  ! (takuya ooura) cdft_
  module subroutine cdft_(length, direction, sequence, ip, work)
    integer(ik), intent(in) :: length
    integer(ik), intent(in) :: direction
    real(rk), intent(inout) :: sequence(0_ik : *)
    integer(ik), intent(inout) :: ip(0_ik : *)
    real(rk), intent(inout) :: work(0_ik : *)
    call make_fft_table_(length/4_ik, ip, work)
    if (direction == fft_forward) then
      call bit_reverse_(length, ip(2_ik), sequence)
      call cft_forward_(length, sequence, work)
    else if(direction == fft_inverse) then
      call bit_reverse_conjugate_(length, ip(2_ik), sequence)
      call cft_inverse_(length, sequence, work)
    end if
  end subroutine cdft_
  ! ############################################################################################################################# !
  ! (takuya ooura) cdft__
  module subroutine cdft__(length, direction, sequence, ip, work)
    integer(ik), intent(in) :: length
    integer(ik), intent(in) :: direction
    real(rk), intent(inout) :: sequence(0_ik : *)
    integer(ik), intent(inout) :: ip(0_ik : *)
    real(rk), intent(inout) :: work(0_ik : *)
    if (direction == fft_forward) then
      call bit_reverse_(length, ip(2_ik), sequence)
      call cft_forward_(length, sequence, work)
    else if(direction == fft_inverse) then
      call bit_reverse_conjugate_(length, ip(2_ik), sequence)
      call cft_inverse_(length, sequence, work)
    end if
  end subroutine cdft__
  ! ############################################################################################################################# !
  ! bit_reverse_
  module subroutine bit_reverse_(n, ip, a)
    integer(ik), intent(in) :: n
    integer(ik), intent(inout) :: ip(0_ik : *)
    real(rk), intent(inout) :: a(0_ik : n - 1_ik)
    integer(ik) :: j, j1, k, k1, l, m, m2
    real(rk) :: xr, xi, yr, yi
    ip(0_ik) = 0_ik
    l = n
    m = 1_ik
    do while (8_ik * m < l)
      l = l / 2_ik
      do j = 0_ik, m - 1_ik
        ip(m + j) = ip(j) + l
      end do
      m = m * 2_ik
    end do
    m2 = 2_ik * m
    if (8_ik * m == l) then
      do k = 0_ik, m - 1_ik
        do j = 0_ik, k - 1_ik
          j1 = 2_ik * j + ip(k)
          k1 = 2_ik * k + ip(j)
          xr = a(j1)
          xi = a(j1 + 1_ik)
          yr = a(k1)
          yi = a(k1 + 1_ik)
          a(j1) = yr
          a(j1 + 1_ik) = yi
          a(k1) = xr
          a(k1 + 1_ik) = xi
          j1 = j1 + m2
          k1 = k1 + 2_ik * m2
          xr = a(j1)
          xi = a(j1 + 1_ik)
          yr = a(k1)
          yi = a(k1 + 1_ik)
          a(j1) = yr
          a(j1 + 1_ik) = yi
          a(k1) = xr
          a(k1 + 1_ik) = xi
          j1 = j1 + m2
          k1 = k1 - m2
          xr = a(j1)
          xi = a(j1 + 1_ik)
          yr = a(k1)
          yi = a(k1 + 1_ik)
          a(j1) = yr
          a(j1 + 1_ik) = yi
          a(k1) = xr
          a(k1 + 1_ik) = xi
          j1 = j1 + m2
          k1 = k1 + 2_ik * m2
          xr = a(j1)
          xi = a(j1 + 1_ik)
          yr = a(k1)
          yi = a(k1 + 1_ik)
          a(j1) = yr
          a(j1 + 1_ik) = yi
          a(k1) = xr
          a(k1 + 1_ik) = xi
        end do
        j1 = 2_ik * k + m2 + ip(k)
        k1 = j1 + m2
        xr = a(j1)
        xi = a(j1 + 1_ik)
        yr = a(k1)
        yi = a(k1 + 1_ik)
        a(j1) = yr
        a(j1 + 1_ik) = yi
        a(k1) = xr
        a(k1 + 1_ik) = xi
      end do
    else
      do k = 1_ik, m - 1_ik
        do j = 0_ik, k - 1_ik
          j1 = 2_ik * j + ip(k)
          k1 = 2_ik * k + ip(j)
          xr = a(j1)
          xi = a(j1 + 1_ik)
          yr = a(k1)
          yi = a(k1 + 1_ik)
          a(j1) = yr
          a(j1 + 1_ik) = yi
          a(k1) = xr
          a(k1 + 1_ik) = xi
          j1 = j1 + m2
          k1 = k1 + m2
          xr = a(j1)
          xi = a(j1 + 1_ik)
          yr = a(k1)
          yi = a(k1 + 1_ik)
          a(j1) = yr
          a(j1 + 1_ik) = yi
          a(k1) = xr
          a(k1 + 1_ik) = xi
        end do
      end do
    end if
  end subroutine bit_reverse_
  ! ############################################################################################################################# !
  ! bit_reverse_conjugate_
  module subroutine bit_reverse_conjugate_(n, ip, a)
    integer(ik), intent(in) :: n
    integer(ik), intent(inout) :: ip(0_ik : *)
    real(rk), intent(inout) :: a(0_ik : n - 1_ik)
    integer(ik) :: j, j1, k, k1, l, m, m2
    real(rk) :: xr, xi, yr, yi
    ip(0_ik) = 0_ik
    l = n
    m = 1_ik
    do while (8_ik * m < l)
      l = l / 2_ik
      do j = 0_ik, m - 1_ik
        ip(m + j) = ip(j) + l
      end do
      m = m * 2_ik
    end do
    m2 = 2_ik * m
    if (8_ik * m == l) then
      do k = 0_ik, m - 1_ik
        do j = 0_ik, k - 1_ik
          j1 = 2_ik * j + ip(k)
          k1 = 2_ik * k + ip(j)
          xr = a(j1)
          xi = -a(j1 + 1_ik)
          yr = a(k1)
          yi = -a(k1 + 1_ik)
          a(j1) = yr
          a(j1 + 1_ik) = yi
          a(k1) = xr
          a(k1 + 1_ik) = xi
          j1 = j1 + m2
          k1 = k1 + 2_ik * m2
          xr = a(j1)
          xi = -a(j1 + 1_ik)
          yr = a(k1)
          yi = -a(k1 + 1_ik)
          a(j1) = yr
          a(j1 + 1_ik) = yi
          a(k1) = xr
          a(k1 + 1_ik) = xi
          j1 = j1 + m2
          k1 = k1 - m2
          xr = a(j1)
          xi = -a(j1 + 1_ik)
          yr = a(k1)
          yi = -a(k1 + 1_ik)
          a(j1) = yr
          a(j1 + 1_ik) = yi
          a(k1) = xr
          a(k1 + 1_ik) = xi
          j1 = j1 + m2
          k1 = k1 + 2_ik * m2
          xr = a(j1)
          xi = -a(j1 + 1_ik)
          yr = a(k1)
          yi = -a(k1 + 1_ik)
          a(j1) = yr
          a(j1 + 1_ik) = yi
          a(k1) = xr
          a(k1 + 1_ik) = xi
        end do
        k1 = 2_ik * k + ip(k)
        a(k1 + 1_ik) = -a(k1 + 1_ik)
        j1 = k1 + m2
        k1 = j1 + m2
        xr = a(j1)
        xi = -a(j1 + 1_ik)
        yr = a(k1)
        yi = -a(k1 + 1_ik)
        a(j1) = yr
        a(j1 + 1_ik) = yi
        a(k1) = xr
        a(k1 + 1_ik) = xi
        k1 = k1 + m2
        a(k1 + 1_ik) = -a(k1 + 1_ik)
      end do
    else
      a(1_ik) = -a(1_ik)
      a(m2 + 1_ik) = -a(m2 + 1_ik)
      do k = 1_ik, m - 1_ik
        do j = 0_ik, k - 1_ik
          j1 = 2_ik * j + ip(k)
          k1 = 2_ik * k + ip(j)
          xr = a(j1)
          xi = -a(j1 + 1_ik)
          yr = a(k1)
          yi = -a(k1 + 1_ik)
          a(j1) = yr
          a(j1 + 1_ik) = yi
          a(k1) = xr
          a(k1 + 1_ik) = xi
          j1 = j1 + m2
          k1 = k1 + m2
          xr = a(j1)
          xi = -a(j1 + 1_ik)
          yr = a(k1)
          yi = -a(k1 + 1_ik)
          a(j1) = yr
          a(j1 + 1_ik) = yi
          a(k1) = xr
          a(k1 + 1_ik) = xi
        end do
        k1 = 2_ik * k + ip(k)
        a(k1 + 1_ik) = -a(k1 + 1_ik)
        a(k1 + m2 + 1_ik) = -a(k1 + m2 + 1_ik)
      end do
    end if
  end subroutine bit_reverse_conjugate_
  ! ############################################################################################################################# !
  ! make_fft_table_
  module subroutine make_fft_table_(nw, ip, w)
    integer(ik), intent(in) :: nw
    integer(ik), intent(inout) :: ip(0_ik : *)
    real(rk), intent(inout) :: w(0_ik : nw - 1_ik)
    integer(ik) :: j, nwh
    real(rk) :: delta, x, y
    nwh = nw / 2_ik
    delta = one_pi / (4.0_rk * real(nwh, rk))
    w(0_ik) = 1_ik
    w(1_ik) = 0_ik
    w(nwh) = cos(delta * real(nwh, rk))
    w(nwh + 1_ik) = w(nwh)
    if (nwh > 2_ik) then
      do j = 2_ik, nwh - 2_ik, 2_ik
        x = cos(delta * real(j,rk))
        y = sin(delta * real(j,rk))
        w(j) = x
        w(j + 1) = y
        w(nw - j) = y
        w(nw - j + 1_ik) = x
      end do
      do j = nwh - 2_ik, 2_ik, -2_ik
        x = w(2_ik * j)
        y = w(2_ik * j + 1_ik)
        w(nwh + j) = x
        w(nwh + j + 1_ik) = y
      end do
      ip(0_ik) = nw
      ip(1_ik) = 1_ik
      call bit_reverse_(nw, ip(2_ik), w)
    end if
  end subroutine make_fft_table_
  ! ############################################################################################################################# !
  ! cft_forward_
  module subroutine cft_forward_(n, a, w)
    integer(ik), intent(in) :: n
    real(rk), intent(inout) :: a(0_ik : n - 1_ik)
    real(rk), intent(in) :: w(0_ik : *)
    integer(ik) :: j, j1, j2, j3, l
    real(rk) :: x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
    l = 2_ik
    if (n >= 16_ik) then
      call cft_1st_(n, a, w)
      l = 16_ik
      do while (8_ik * l <= n)
        call cft_mdl_(n, l, a, w)
        l = 8_ik * l
      end do
    end if
    if (2_ik * l < n) then
      do j = 0_ik, l - 2_ik, 2_ik
        j1 = j + l
        j2 = j1 + l
        j3 = j2 + l
        x0r = a(j) + a(j1)
        x0i = a(j + 1_ik) + a(j1 + 1_ik)
        x1r = a(j) - a(j1)
        x1i = a(j + 1_ik) - a(j1 + 1_ik)
        x2r = a(j2) + a(j3)
        x2i = a(j2 + 1_ik) + a(j3 + 1_ik)
        x3r = a(j2) - a(j3)
        x3i = a(j2 + 1_ik) - a(j3 + 1_ik)
        a(j) = x0r + x2r
        a(j + 1_ik) = x0i + x2i
        a(j2) = x0r - x2r
        a(j2 + 1_ik) = x0i - x2i
        a(j1) = x1r - x3i
        a(j1 + 1_ik) = x1i + x3r
        a(j3) = x1r + x3i
        a(j3 + 1_ik) = x1i - x3r
      end do
    else if (2_ik * l == n) then
      do j = 0_ik, l - 2_ik, 2_ik
        j1 = j + l
        x0r = a(j) - a(j1)
        x0i = a(j + 1_ik) - a(j1 + 1_ik)
        a(j) = a(j) + a(j1)
        a(j + 1_ik) = a(j + 1_ik) + a(j1 + 1_ik)
        a(j1) = x0r
        a(j1 + 1_ik) = x0i
      end do
    end if
  end subroutine cft_forward_
  ! ############################################################################################################################# !
  ! cft_inverse_
  module subroutine cft_inverse_(n, a, w)
    integer(ik), intent(in) :: n
    real(rk), intent(inout) :: a(0_ik : n - 1_ik)
    real(rk), intent(in) :: w(0_ik : *)
    integer(ik) :: j, j1, j2, j3, j4, j5, j6, j7, l
    real(rk) ::  wn4r, x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
    real(rk) :: y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i
    real(rk) :: y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i
    l = 2_ik
    if (n > 16_ik) then
      call cft_1st_(n, a, w)
      l = 16_ik
      do while (8_ik * l < n)
        call cft_mdl_(n, l, a, w)
        l = 8_ik * l
      end do
    end if
    if (4_ik * l < n) then
      wn4r = w(2_ik)
      do j = 0_ik, l - 2_ik, 2_ik
        j1 = j + l
        j2 = j1 + l
        j3 = j2 + l
        j4 = j3 + l
        j5 = j4 + l
        j6 = j5 + l
        j7 = j6 + l
        x0r = a(j) + a(j1)
        x0i = -a(j + 1_ik) - a(j1 + 1_ik)
        x1r = a(j) - a(j1)
        x1i = -a(j + 1_ik) + a(j1 + 1_ik)
        x2r = a(j2) + a(j3)
        x2i = a(j2 + 1_ik) + a(j3 + 1_ik)
        x3r = a(j2) - a(j3)
        x3i = a(j2 + 1_ik) - a(j3 + 1_ik)
        y0r = x0r + x2r
        y0i = x0i - x2i
        y2r = x0r - x2r
        y2i = x0i + x2i
        y1r = x1r - x3i
        y1i = x1i - x3r
        y3r = x1r + x3i
        y3i = x1i + x3r
        x0r = a(j4) + a(j5)
        x0i = a(j4 + 1_ik) + a(j5 + 1_ik)
        x1r = a(j4) - a(j5)
        x1i = a(j4 + 1_ik) - a(j5 + 1_ik)
        x2r = a(j6) + a(j7)
        x2i = a(j6 + 1_ik) + a(j7 + 1_ik)
        x3r = a(j6) - a(j7)
        x3i = a(j6 + 1_ik) - a(j7 + 1_ik)
        y4r = x0r + x2r
        y4i = x0i + x2i
        y6r = x0r - x2r
        y6i = x0i - x2i
        x0r = x1r - x3i
        x0i = x1i + x3r
        x2r = x1r + x3i
        x2i = x1i - x3r
        y5r = wn4r * (x0r - x0i)
        y5i = wn4r * (x0r + x0i)
        y7r = wn4r * (x2r - x2i)
        y7i = wn4r * (x2r + x2i)
        a(j1) = y1r + y5r
        a(j1 + 1_ik) = y1i - y5i
        a(j5) = y1r - y5r
        a(j5 + 1_ik) = y1i + y5i
        a(j3) = y3r - y7i
        a(j3 + 1_ik) = y3i - y7r
        a(j7) = y3r + y7i
        a(j7 + 1_ik) = y3i + y7r
        a(j) = y0r + y4r
        a(j + 1_ik) = y0i - y4i
        a(j4) = y0r - y4r
        a(j4 + 1_ik) = y0i + y4i
        a(j2) = y2r - y6i
        a(j2 + 1_ik) = y2i - y6r
        a(j6) = y2r + y6i
        a(j6 + 1_ik) = y2i + y6r
      end do
    else if (4_ik * l == n) then
      do j = 0_ik, l - 2_ik, 2_ik
        j1 = j + l
        j2 = j1 + l
        j3 = j2 + l
        x0r = a(j) + a(j1)
        x0i = -a(j + 1_ik) - a(j1 + 1_ik)
        x1r = a(j) - a(j1)
        x1i = -a(j + 1_ik) + a(j1 + 1_ik)
        x2r = a(j2) + a(j3)
        x2i = a(j2 + 1_ik) + a(j3 + 1_ik)
        x3r = a(j2) - a(j3)
        x3i = a(j2 + 1_ik) - a(j3 + 1_ik)
        a(j) = x0r + x2r
        a(j + 1_ik) = x0i - x2i
        a(j2) = x0r - x2r
        a(j2 + 1_ik) = x0i + x2i
        a(j1) = x1r - x3i
        a(j1 + 1_ik) = x1i - x3r
        a(j3) = x1r + x3i
        a(j3 + 1_ik) = x1i + x3r
      end do
    else
      do j = 0_ik, l - 2_ik, 2_ik
        j1 = j + l
        x0r = a(j) - a(j1)
        x0i = -a(j + 1_ik) + a(j1 + 1_ik)
        a(j) = a(j) + a(j1)
        a(j + 1_ik) = -a(j + 1_ik) - a(j1 + 1_ik)
        a(j1) = x0r
        a(j1 + 1_ik) = x0i
      end do
    end if
  end subroutine cft_inverse_
  ! ############################################################################################################################# !
  ! cft_1st_
  module subroutine cft_1st_(n, a, w)
    integer(ik), intent(in) :: n
    real(rk), intent(inout) :: a(0_ik : n - 1_ik)
    real(rk), intent(in) :: w(0_ik : *)
    integer(ik) :: j, k1
    real(rk) ::  wn4r, wtmp, wk1r, wk1i, wk2r, wk2i, wk3r, wk3i
    real(rk) ::  wk4r, wk4i, wk5r, wk5i, wk6r, wk6i, wk7r, wk7i
    real(rk) ::  x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
    real(rk) ::  y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i
    real(rk) :: y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i
    wn4r = w(2_ik)
    x0r = a(0_ik) + a(2_ik)
    x0i = a(1_ik) + a(3_ik)
    x1r = a(0_ik) - a(2_ik)
    x1i = a(1_ik) - a(3_ik)
    x2r = a(4_ik) + a(6_ik)
    x2i = a(5_ik) + a(7_ik)
    x3r = a(4_ik) - a(6_ik)
    x3i = a(5_ik) - a(7_ik)
    y0r = x0r + x2r
    y0i = x0i + x2i
    y2r = x0r - x2r
    y2i = x0i - x2i
    y1r = x1r - x3i
    y1i = x1i + x3r
    y3r = x1r + x3i
    y3i = x1i - x3r
    x0r = a(8_ik) + a(10_ik)
    x0i = a(9_ik) + a(11_ik)
    x1r = a(8_ik) - a(10_ik)
    x1i = a(9_ik) - a(11_ik)
    x2r = a(12_ik) + a(14_ik)
    x2i = a(13_ik) + a(15_ik)
    x3r = a(12_ik) - a(14_ik)
    x3i = a(13_ik) - a(15_ik)
    y4r = x0r + x2r
    y4i = x0i + x2i
    y6r = x0r - x2r
    y6i = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    x2r = x1r + x3i
    x2i = x1i - x3r
    y5r = wn4r * (x0r - x0i)
    y5i = wn4r * (x0r + x0i)
    y7r = wn4r * (x2r - x2i)
    y7i = wn4r * (x2r + x2i)
    a(2_ik) = y1r + y5r
    a(3_ik) = y1i + y5i
    a(10_ik) = y1r - y5r
    a(11_ik) = y1i - y5i
    a(6_ik) = y3r - y7i
    a(7_ik) = y3i + y7r
    a(14_ik) = y3r + y7i
    a(15_ik) = y3i - y7r
    a(0_ik) = y0r + y4r
    a(1_ik) = y0i + y4i
    a(8_ik) = y0r - y4r
    a(9_ik) = y0i - y4i
    a(4_ik) = y2r - y6i
    a(5_ik) = y2i + y6r
    a(12_ik) = y2r + y6i
    a(13_ik) = y2i - y6r
    if (n > 16_ik) then
      wk1r = w(4_ik)
      wk1i = w(5_ik)
      x0r = a(16_ik) + a(18_ik)
      x0i = a(17_ik) + a(19_ik)
      x1r = a(16_ik) - a(18_ik)
      x1i = a(17_ik) - a(19_ik)
      x2r = a(20_ik) + a(22_ik)
      x2i = a(21_ik) + a(23_ik)
      x3r = a(20_ik) - a(22_ik)
      x3i = a(21_ik) - a(23_ik)
      y0r = x0r + x2r
      y0i = x0i + x2i
      y2r = x0r - x2r
      y2i = x0i - x2i
      y1r = x1r - x3i
      y1i = x1i + x3r
      y3r = x1r + x3i
      y3i = x1i - x3r
      x0r = a(24_ik) + a(26_ik)
      x0i = a(25_ik) + a(27_ik)
      x1r = a(24_ik) - a(26_ik)
      x1i = a(25_ik) - a(27_ik)
      x2r = a(28_ik) + a(30_ik)
      x2i = a(29_ik) + a(31_ik)
      x3r = a(28_ik) - a(30_ik)
      x3i = a(29_ik) - a(31_ik)
      y4r = x0r + x2r
      y4i = x0i + x2i
      y6r = x0r - x2r
      y6i = x0i - x2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      x2r = x1r + x3i
      x2i = x3r - x1i
      y5r = wk1i * x0r - wk1r * x0i
      y5i = wk1i * x0i + wk1r * x0r
      y7r = wk1r * x2r + wk1i * x2i
      y7i = wk1r * x2i - wk1i * x2r
      x0r = wk1r * y1r - wk1i * y1i
      x0i = wk1r * y1i + wk1i * y1r
      a(18_ik) = x0r + y5r
      a(19_ik) = x0i + y5i
      a(26_ik) = y5i - x0i
      a(27_ik) = x0r - y5r
      x0r = wk1i * y3r - wk1r * y3i
      x0i = wk1i * y3i + wk1r * y3r
      a(22_ik) = x0r - y7r
      a(23_ik) = x0i + y7i
      a(30_ik) = y7i - x0i
      a(31_ik) = x0r + y7r
      a(16_ik) = y0r + y4r
      a(17_ik) = y0i + y4i
      a(24_ik) = y4i - y0i
      a(25_ik) = y0r - y4r
      x0r = y2r - y6i
      x0i = y2i + y6r
      a(20_ik) = wn4r * (x0r - x0i)
      a(21_ik) = wn4r * (x0i + x0r)
      x0r = y6r - y2i
      x0i = y2r + y6i
      a(28_ik) = wn4r * (x0r - x0i)
      a(29_ik) = wn4r * (x0i + x0r)
      k1 = 4_ik
      do j = 32_ik, n - 16_ik, 16_ik
        k1 = k1 + 4_ik
        wk1r = w(k1)
        wk1i = w(k1 + 1_ik)
        wk2r = w(k1 + 2_ik)
        wk2i = w(k1 + 3_ik)
        wtmp = 2_ik * wk2i
        wk3r = wk1r - wtmp * wk1i
        wk3i = wtmp * wk1r - wk1i
        wk4r = 1_ik - wtmp * wk2i
        wk4i = wtmp * wk2r
        wtmp = 2_ik * wk4i
        wk5r = wk3r - wtmp * wk1i
        wk5i = wtmp * wk1r - wk3i
        wk6r = wk2r - wtmp * wk2i
        wk6i = wtmp * wk2r - wk2i
        wk7r = wk1r - wtmp * wk3i
        wk7i = wtmp * wk3r - wk1i
        x0r = a(j) + a(j + 2_ik)
        x0i = a(j + 1_ik) + a(j + 3_ik)
        x1r = a(j) - a(j + 2_ik)
        x1i = a(j + 1_ik) - a(j + 3_ik)
        x2r = a(j + 4_ik) + a(j + 6_ik)
        x2i = a(j + 5_ik) + a(j + 7_ik)
        x3r = a(j + 4_ik) - a(j + 6_ik)
        x3i = a(j + 5_ik) - a(j + 7_ik)
        y0r = x0r + x2r
        y0i = x0i + x2i
        y2r = x0r - x2r
        y2i = x0i - x2i
        y1r = x1r - x3i
        y1i = x1i + x3r
        y3r = x1r + x3i
        y3i = x1i - x3r
        x0r = a(j + 8_ik) + a(j + 10_ik)
        x0i = a(j + 9_ik) + a(j + 11_ik)
        x1r = a(j + 8_ik) - a(j + 10_ik)
        x1i = a(j + 9_ik) - a(j + 11_ik)
        x2r = a(j + 12_ik) + a(j + 14_ik)
        x2i = a(j + 13_ik) + a(j + 15_ik)
        x3r = a(j + 12_ik) - a(j + 14_ik)
        x3i = a(j + 13_ik) - a(j + 15_ik)
        y4r = x0r + x2r
        y4i = x0i + x2i
        y6r = x0r - x2r
        y6i = x0i - x2i
        x0r = x1r - x3i
        x0i = x1i + x3r
        x2r = x1r + x3i
        x2i = x1i - x3r
        y5r = wn4r * (x0r - x0i)
        y5i = wn4r * (x0r + x0i)
        y7r = wn4r * (x2r - x2i)
        y7i = wn4r * (x2r + x2i)
        x0r = y1r + y5r
        x0i = y1i + y5i
        a(j + 2_ik) = wk1r * x0r - wk1i * x0i
        a(j + 3_ik) = wk1r * x0i + wk1i * x0r
        x0r = y1r - y5r
        x0i = y1i - y5i
        a(j + 10_ik) = wk5r * x0r - wk5i * x0i
        a(j + 11_ik) = wk5r * x0i + wk5i * x0r
        x0r = y3r - y7i
        x0i = y3i + y7r
        a(j + 6_ik) = wk3r * x0r - wk3i * x0i
        a(j + 7_ik) = wk3r * x0i + wk3i * x0r
        x0r = y3r + y7i
        x0i = y3i - y7r
        a(j + 14_ik) = wk7r * x0r - wk7i * x0i
        a(j + 15_ik) = wk7r * x0i + wk7i * x0r
        a(j) = y0r + y4r
        a(j + 1_ik) = y0i + y4i
        x0r = y0r - y4r
        x0i = y0i - y4i
        a(j + 8_ik) = wk4r * x0r - wk4i * x0i
        a(j + 9_ik) = wk4r * x0i + wk4i * x0r
        x0r = y2r - y6i
        x0i = y2i + y6r
        a(j + 4_ik) = wk2r * x0r - wk2i * x0i
        a(j + 5_ik) = wk2r * x0i + wk2i * x0r
        x0r = y2r + y6i
        x0i = y2i - y6r
        a(j + 12_ik) = wk6r * x0r - wk6i * x0i
        a(j + 13_ik) = wk6r * x0i + wk6i * x0r
      end do
    end if
  end subroutine cft_1st_
  ! ############################################################################################################################# !
  ! cft_mdl_
  module subroutine cft_mdl_(n, l, a, w)
    integer(ik), intent(in) :: n
    integer(ik), intent(in) :: l
    real(rk), intent(inout) :: a(0_ik : n - 1_ik)
    real(rk), intent(in) :: w(0_ik : *)
    integer(ik) :: j, j1, j2, j3, j4, j5, j6, j7, k, k1, m
    real(rk) :: wn4r, wtmp, wk1r, wk1i, wk2r, wk2i, wk3r, wk3i
    real(rk) :: wk4r, wk4i, wk5r, wk5i, wk6r, wk6i, wk7r, wk7i
    real(rk) :: x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
    real(rk) :: y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i
    real(rk) :: y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i
    m = 8_ik * l
    wn4r = w(2_ik)
    do j = 0_ik, l - 2_ik, 2_ik
      j1 = j + l
      j2 = j1 + l
      j3 = j2 + l
      j4 = j3 + l
      j5 = j4 + l
      j6 = j5 + l
      j7 = j6 + l
      x0r = a(j) + a(j1)
      x0i = a(j + 1_ik) + a(j1 + 1_ik)
      x1r = a(j) - a(j1)
      x1i = a(j + 1_ik) - a(j1 + 1_ik)
      x2r = a(j2) + a(j3)
      x2i = a(j2 + 1_ik) + a(j3 + 1_ik)
      x3r = a(j2) - a(j3)
      x3i = a(j2 + 1_ik) - a(j3 + 1_ik)
      y0r = x0r + x2r
      y0i = x0i + x2i
      y2r = x0r - x2r
      y2i = x0i - x2i
      y1r = x1r - x3i
      y1i = x1i + x3r
      y3r = x1r + x3i
      y3i = x1i - x3r
      x0r = a(j4) + a(j5)
      x0i = a(j4 + 1_ik) + a(j5 + 1_ik)
      x1r = a(j4) - a(j5)
      x1i = a(j4 + 1_ik) - a(j5 + 1_ik)
      x2r = a(j6) + a(j7)
      x2i = a(j6 + 1_ik) + a(j7 + 1_ik)
      x3r = a(j6) - a(j7)
      x3i = a(j6 + 1_ik) - a(j7 + 1_ik)
      y4r = x0r + x2r
      y4i = x0i + x2i
      y6r = x0r - x2r
      y6i = x0i - x2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      x2r = x1r + x3i
      x2i = x1i - x3r
      y5r = wn4r * (x0r - x0i)
      y5i = wn4r * (x0r + x0i)
      y7r = wn4r * (x2r - x2i)
      y7i = wn4r * (x2r + x2i)
      a(j1) = y1r + y5r
      a(j1 + 1_ik) = y1i + y5i
      a(j5) = y1r - y5r
      a(j5 + 1_ik) = y1i - y5i
      a(j3) = y3r - y7i
      a(j3 + 1_ik) = y3i + y7r
      a(j7) = y3r + y7i
      a(j7 + 1_ik) = y3i - y7r
      a(j) = y0r + y4r
      a(j + 1_ik) = y0i + y4i
      a(j4) = y0r - y4r
      a(j4 + 1_ik) = y0i - y4i
      a(j2) = y2r - y6i
      a(j2 + 1_ik) = y2i + y6r
      a(j6) = y2r + y6i
      a(j6 + 1_ik) = y2i - y6r
    end do
    if (m < n) then
      wk1r = w(4_ik)
      wk1i = w(5_ik)
      do j = m, l + m - 2_ik, 2_ik
        j1 = j + l
        j2 = j1 + l
        j3 = j2 + l
        j4 = j3 + l
        j5 = j4 + l
        j6 = j5 + l
        j7 = j6 + l
        x0r = a(j) + a(j1)
        x0i = a(j + 1_ik) + a(j1 + 1_ik)
        x1r = a(j) - a(j1)
        x1i = a(j + 1_ik) - a(j1 + 1_ik)
        x2r = a(j2) + a(j3)
        x2i = a(j2 + 1_ik) + a(j3 + 1_ik)
        x3r = a(j2) - a(j3)
        x3i = a(j2 + 1_ik) - a(j3 + 1_ik)
        y0r = x0r + x2r
        y0i = x0i + x2i
        y2r = x0r - x2r
        y2i = x0i - x2i
        y1r = x1r - x3i
        y1i = x1i + x3r
        y3r = x1r + x3i
        y3i = x1i - x3r
        x0r = a(j4) + a(j5)
        x0i = a(j4 + 1_ik) + a(j5 + 1_ik)
        x1r = a(j4) - a(j5)
        x1i = a(j4 + 1_ik) - a(j5 + 1_ik)
        x2r = a(j6) + a(j7)
        x2i = a(j6 + 1_ik) + a(j7 + 1_ik)
        x3r = a(j6) - a(j7)
        x3i = a(j6 + 1_ik) - a(j7 + 1_ik)
        y4r = x0r + x2r
        y4i = x0i + x2i
        y6r = x0r - x2r
        y6i = x0i - x2i
        x0r = x1r - x3i
        x0i = x1i + x3r
        x2r = x1r + x3i
        x2i = x3r - x1i
        y5r = wk1i * x0r - wk1r * x0i
        y5i = wk1i * x0i + wk1r * x0r
        y7r = wk1r * x2r + wk1i * x2i
        y7i = wk1r * x2i - wk1i * x2r
        x0r = wk1r * y1r - wk1i * y1i
        x0i = wk1r * y1i + wk1i * y1r
        a(j1) = x0r + y5r
        a(j1 + 1_ik) = x0i + y5i
        a(j5) = y5i - x0i
        a(j5 + 1_ik) = x0r - y5r
        x0r = wk1i * y3r - wk1r * y3i
        x0i = wk1i * y3i + wk1r * y3r
        a(j3) = x0r - y7r
        a(j3 + 1_ik) = x0i + y7i
        a(j7) = y7i - x0i
        a(j7 + 1_ik) = x0r + y7r
        a(j) = y0r + y4r
        a(j + 1_ik) = y0i + y4i
        a(j4) = y4i - y0i
        a(j4 + 1_ik) = y0r - y4r
        x0r = y2r - y6i
        x0i = y2i + y6r
        a(j2) = wn4r * (x0r - x0i)
        a(j2 + 1_ik) = wn4r * (x0i + x0r)
        x0r = y6r - y2i
        x0i = y2r + y6i
        a(j6) = wn4r * (x0r - x0i)
        a(j6 + 1_ik) = wn4r * (x0i + x0r)
      end do
      k1 = 4_ik
      do k = 2_ik * m, n - m, m
        k1 = k1 + 4_ik
        wk1r = w(k1)
        wk1i = w(k1 + 1_ik)
        wk2r = w(k1 + 2_ik)
        wk2i = w(k1 + 3_ik)
        wtmp = 2_ik * wk2i
        wk3r = wk1r - wtmp * wk1i
        wk3i = wtmp * wk1r - wk1i
        wk4r = 1_ik - wtmp * wk2i
        wk4i = wtmp * wk2r
        wtmp = 2_ik * wk4i
        wk5r = wk3r - wtmp * wk1i
        wk5i = wtmp * wk1r - wk3i
        wk6r = wk2r - wtmp * wk2i
        wk6i = wtmp * wk2r - wk2i
        wk7r = wk1r - wtmp * wk3i
        wk7i = wtmp * wk3r - wk1i
        do j = k, l + k - 2_ik, 2_ik
          j1 = j + l
          j2 = j1 + l
          j3 = j2 + l
          j4 = j3 + l
          j5 = j4 + l
          j6 = j5 + l
          j7 = j6 + l
          x0r = a(j) + a(j1)
          x0i = a(j + 1_ik) + a(j1 + 1_ik)
          x1r = a(j) - a(j1)
          x1i = a(j + 1_ik) - a(j1 + 1_ik)
          x2r = a(j2) + a(j3)
          x2i = a(j2 + 1_ik) + a(j3 + 1_ik)
          x3r = a(j2) - a(j3)
          x3i = a(j2 + 1_ik) - a(j3 + 1_ik)
          y0r = x0r + x2r
          y0i = x0i + x2i
          y2r = x0r - x2r
          y2i = x0i - x2i
          y1r = x1r - x3i
          y1i = x1i + x3r
          y3r = x1r + x3i
          y3i = x1i - x3r
          x0r = a(j4) + a(j5)
          x0i = a(j4 + 1_ik) + a(j5 + 1_ik)
          x1r = a(j4) - a(j5)
          x1i = a(j4 + 1_ik) - a(j5 + 1_ik)
          x2r = a(j6) + a(j7)
          x2i = a(j6 + 1_ik) + a(j7 + 1)
          x3r = a(j6) - a(j7)
          x3i = a(j6 + 1_ik) - a(j7 + 1_ik)
          y4r = x0r + x2r
          y4i = x0i + x2i
          y6r = x0r - x2r
          y6i = x0i - x2i
          x0r = x1r - x3i
          x0i = x1i + x3r
          x2r = x1r + x3i
          x2i = x1i - x3r
          y5r = wn4r * (x0r - x0i)
          y5i = wn4r * (x0r + x0i)
          y7r = wn4r * (x2r - x2i)
          y7i = wn4r * (x2r + x2i)
          x0r = y1r + y5r
          x0i = y1i + y5i
          a(j1) = wk1r * x0r - wk1i * x0i
          a(j1 + 1_ik) = wk1r * x0i + wk1i * x0r
          x0r = y1r - y5r
          x0i = y1i - y5i
          a(j5) = wk5r * x0r - wk5i * x0i
          a(j5 + 1_ik) = wk5r * x0i + wk5i * x0r
          x0r = y3r - y7i
          x0i = y3i + y7r
          a(j3) = wk3r * x0r - wk3i * x0i
          a(j3 + 1_ik) = wk3r * x0i + wk3i * x0r
          x0r = y3r + y7i
          x0i = y3i - y7r
          a(j7) = wk7r * x0r - wk7i * x0i
          a(j7 + 1_ik) = wk7r * x0i + wk7i * x0r
          a(j) = y0r + y4r
          a(j + 1_ik) = y0i + y4i
          x0r = y0r - y4r
          x0i = y0i - y4i
          a(j4) = wk4r * x0r - wk4i * x0i
          a(j4 + 1_ik) = wk4r * x0i + wk4i * x0r
          x0r = y2r - y6i
          x0i = y2i + y6r
          a(j2) = wk2r * x0r - wk2i * x0i
          a(j2 + 1_ik) = wk2r * x0i + wk2i * x0r
          x0r = y2r + y6i
          x0i = y2i - y6r
          a(j6) = wk6r * x0r - wk6i * x0i
          a(j6 + 1_ik) = wk6r * x0i + wk6i * x0r
        end do
      end do
    end if
  end subroutine cft_mdl_
  ! ############################################################################################################################# !
end submodule transformation
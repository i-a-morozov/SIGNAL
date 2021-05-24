! signal, 2018-2020, i.a.morozov@inp.nsk.su
module signal
  use, intrinsic :: iso_c_binding,   only: ik => c_int, rk => c_double, c_sizeof
  implicit none
  private
  ! ############################################################################################################################# !
  ! global
  ! ############################################################################################################################# !
  public :: ik
  public :: rk
  integer(ik), public, parameter :: ik_size                = c_sizeof(ik)
  integer(ik), public, parameter :: rk_size                = c_sizeof(rk)
  real(rk),    public, parameter :: one_pi                 = 2.0_rk*acos(0.0_rk)
  real(rk),    public, parameter :: two_pi                 = 2.0_rk*one_pi
  real(rk),    public, parameter :: epsilon                = 1.e-16_rk
  ! ############################################################################################################################# !
  ! external
  ! ############################################################################################################################# !
  external :: dgemv               ! (blas)
  external :: dgemm               ! (blas)
  external :: dgesvd              ! (lapack)
  external :: dsaupd              ! (arpack)
  external :: dseupd              ! (arpack)
  external :: dfftw_plan_dft_1d   ! (fftw)
  external :: dfftw_execute_dft   ! (fftw)
  external :: dfftw_destroy_plan  ! (fftw)
  ! ############################################################################################################################# !
  ! auxiliary
  ! ############################################################################################################################# !
  ! factorial
  ! (function) factorial_(<number>)
  ! <number>               -- (in)     number (ik)
  ! <factorial_>           -- (out)    factorial of <n> (rk)
  interface
    module real(rk) function factorial_(number)
      integer(ik), intent(in) :: number
    end function factorial_
  end interface
  public :: factorial_
  ! ############################################################################################################################# !
  ! gamma (gsl)
  ! (function) gamma_(<number>)
  ! <number>               -- (in)     number (rk)
  ! <gamma_>               -- (out)    gamma of <n> (rk)
  interface gamma_
    real(rk) function gamma_(number) &
      bind(c, name = "gsl_sf_gamma")
      import :: rk
      real(rk), value :: number
      end function gamma_
  end interface gamma_
  ! ############################################################################################################################# !
  ! gamma incomplete (gsl)
  ! (function) gamma_incomplete_(<a>, <x>)
  ! <a>                    -- (in)     a (rk)
  ! <x>                    -- (in)     x (rk)
  ! <gamma_incomplete_>    -- (out)    gamma incomplete of <a> and <x> (rk)
  interface gamma_incomplete_
    real(rk) function gamma_incomplete_(a, x) &
      bind(c, name = "gsl_sf_gamma_inc")
      import :: rk
      real(rk), value :: a
      real(rk), value :: x
    end function gamma_incomplete_
  end interface gamma_incomplete_
  ! ############################################################################################################################# !
  ! gamma regularized
  ! (function) gamma_regularized_(<a>, <x>, <y>)
  ! <a>                    -- (in)     a (rk)
  ! <x>                    -- (in)     x (rk)
  ! <y>                    -- (in)     y (rk)
  ! <gamma_regularized_>   -- (out)    gamma regularized of <a>, <x> and <y> (rk)
  interface
    module real(rk) function gamma_regularized_(a, x, y)
      real(rk), intent(in) :: a
      real(rk), intent(in) :: x
      real(rk), intent(in) :: y
    end function gamma_regularized_
  end interface
  ! ############################################################################################################################# !
  ! generic gamma
  interface gamma_
    procedure factorial_
    procedure gamma_incomplete_
    procedure gamma_regularized_
  end interface gamma_
  public :: gamma_
  ! ############################################################################################################################# !
  ! modified bessel i_0(x) (gsl)
  ! (function) bessel_(<number>)
  ! <number>               -- (in)     number (rk)
  ! <bessel_>              -- (out)    bessel i_0(<number>) (rk)
  interface bessel_
    real(rk) function bessel_(number) &
      bind(c, name = "gsl_sf_bessel_I0")
      import :: rk
      real(rk), value :: number
      end function bessel_
  end interface bessel_
  public :: bessel_
  ! ############################################################################################################################# !
  ! minloc
  ! (function) minloc_(<sequence>)
  ! <sequence>             -- (in)     sequence (rk array)
  ! <minloc_>              -- (out)    minimum location (ik)
  interface
    module integer(ik) function minloc_(sequence, empty)
      real(rk), dimension(:), contiguous, intent(in) :: sequence
      integer, intent(in) :: empty
    end function minloc_
  end interface
  ! ############################################################################################################################# !
  ! maxloc
  ! (function) maxloc_(<sequence>)
  ! <sequence>             -- (in)     sequence (rk array)
  ! <maxloc_>              -- (out)    maximum location (ik)
  interface
    module integer(ik) function maxloc_(sequence, empty)
      real(rk), dimension(:), contiguous, intent(in) :: sequence
      integer, intent(in) :: empty
    end function maxloc_
  end  interface
  ! ############################################################################################################################# !
  ! sort (bubble, descending)
  ! (subroutine) sort_bubble_(<length>, <sequence>, <fst>, <lst>)
  ! <length>               -- (in)     sequence length (ik)
  ! <sequence>             -- (in)     (unsorted) sequence (rk array of length = <length>)
  ! <sequence>             -- (out)    (sorted, descending) sequence (rk array of length = <length>)
  interface
    module subroutine sort_bubble_(length, sequence, fst, lst)
      integer(ik), intent(in) :: length
      real(rk), dimension(length), intent(inout) :: sequence
      integer(ik), intent(in) :: fst
      integer(ik), intent(in) :: lst
    end subroutine sort_bubble_
  end interface
  ! ############################################################################################################################# !
  ! sort (quick, descending)
  ! (subroutine) sort_quick_(<length>, <sequence>, <fst>, <lst>)
  ! <length>               -- (in)     sequence length (ik)
  ! <sequence>             -- (in)     (unsorted) sequence (rk array of length = <length>)
  ! <sequence>             -- (out)    (sorted, descending) sequence (rk array of length = <length>)
  interface
    module recursive subroutine sort_quick_(length, sequence, fst, lst)
      integer(ik), intent(in) :: length
      real(rk), dimension(:), intent(inout) :: sequence
      integer(ik), intent(in) :: fst
      integer(ik), intent(in) :: lst
    end subroutine sort_quick_
  end interface
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
  interface
    module subroutine generate_signal_(flag, length, sequence, loop, frequency, cos_amp, sin_amp) &
      bind(c, name = "generate_signal_")
      integer(ik), intent(in) :: flag
      integer(ik), intent(in) :: length
      real(rk), dimension(2_ik*length), intent(out) :: sequence
      integer(ik), intent(in) :: loop
      real(rk), dimension(loop), intent(in) :: frequency
      real(rk), dimension(loop), intent(in) :: cos_amp
      real(rk), dimension(loop), intent(in) :: sin_amp
    end subroutine generate_signal_
  end interface
  public :: generate_signal_
  ! ############################################################################################################################# !
  ! transformation
  ! ############################################################################################################################# !
  integer(ik), public, parameter :: fft_forward            = +1_ik               ! forward fft
  integer(ik), public, parameter :: fft_inverse            = -1_ik               ! inverse fft
  ! ############################################################################################################################# !
  ! fft/ffrft data memorization
  type table
    integer(ik), dimension(:), allocatable :: bit_fft
    integer(ik), dimension(:), allocatable :: bit_ffrft
    real(rk), dimension(:), allocatable :: trig_fft
    real(rk), dimension(:), allocatable :: trig_ffrft
    real(rk), dimension(:), allocatable :: cos_fst
    real(rk), dimension(:), allocatable :: sin_fst
    real(rk), dimension(:), allocatable :: cos_lst
    real(rk), dimension(:), allocatable :: sin_lst
  end type
  ! ############################################################################################################################# !
  ! fft/ffrft data memorization container
  type(table), protected :: bank
  ! ############################################################################################################################# !
  ! (linear) fractional complex discrete fourier transform
  ! (subroutine) ffrft_(<length>, <argument>, <sequence>)
  ! <length>               -- (in)     length (ik)
  ! <argument>             -- (in)     parameter (rk)
  ! <sequence>             -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <sequence>             -- (out)    fcdft (rk array of length = 2_ik*<length>), <sequence> = [..., fr_i, fi_i, ...]
  ! void    ffrft_(int* length, double* argument, double* sequence) ;
  interface
    module subroutine ffrft_(length, argument, sequence) &
      bind(c, name = "ffrft_")
      integer(ik), intent(in) :: length
      real(rk), intent(in) :: argument
      real(rk), dimension(2_ik*length), intent(inout) :: sequence
    end subroutine
  end interface
  public :: ffrft_
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
  interface
  module subroutine ffrft__(length, sequence, ip, work, cos_fst, sin_fst, cos_lst, sin_lst)
    integer(ik), intent(in) :: length
    real(rk), dimension(2_ik*length), intent(inout) :: sequence
    integer(ik), dimension(0_ik : 1_ik+int(sqrt(real(length, rk)), ik)), intent(in) :: ip
    real(rk), dimension(0_ik : length-1_ik), intent(in) :: work
    real(rk), dimension(length), intent(in)   :: cos_fst
    real(rk), dimension(length), intent(in)   :: sin_fst
    real(rk), dimension(length), intent(in)   :: cos_lst
    real(rk), dimension(length), intent(in)   :: sin_lst
    end subroutine ffrft__
  end interface
  public :: ffrft__
  ! ############################################################################################################################# !
  ! (fftw) complex discrete fourier transform
  ! (subroutine) fft_external_(<length>, <direction>, <sequence>)
  ! <length>               -- (in)     length (ik)
  ! <direction>            -- (in)     direction (ik), fft_forward = +1_ik or fft_inverse = -1_ik
  ! <sequence>             -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <sequence>             -- (out)    cdft (rk array of length = 2_ik*<length>), <sequence> = [..., fr_i, fi_i, ...]
  ! void    fft_external_(int* length, int* direction, double* sequence) ;
  interface
    module subroutine fft_external_(length, direction, sequence) &
      bind(c, name = "fft_external_")
      integer(ik), intent(in) :: length
      integer(ik), intent(in) :: direction
      real(rk), dimension(2_ik*length), intent(inout) :: sequence
    end subroutine fft_external_
  end interface
  public :: fft_external_
  ! ############################################################################################################################# !
  ! (nrf77) complex discrete fourier transform
  ! (subroutine) fft_radix_two_(<length>, <direction>, <sequence>)
  ! <length>               -- (in)     length (ik)
  ! <direction>            -- (in)     direction (ik), fft_forward = +1_ik or fft_inverse = -1_ik
  ! <sequence>             -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <sequence>             -- (out)    cdft (rk array of length = 2_ik*<length>), <sequence> = [..., fr_i, fi_i, ...]
  ! void    fft_radix_two_(int* length, int* direction, double* sequence) ;
  interface
    module subroutine fft_radix_two_(length, direction, sequence) &
      bind(c, name = "fft_radix_two_")
      integer(ik), intent(in) :: length
      integer(ik), intent(in) :: direction
      real(rk), dimension(2_ik*length), intent(inout) :: sequence
    end subroutine fft_radix_two_
  end interface
  public :: fft_radix_two_
  ! ############################################################################################################################# !
  ! (takuya ooura) complex discrete fourier transform
  ! (subroutine) fft_radix_eight_(<length>, <direction>, <sequence>)
  ! <length>               -- (in)     length (ik)
  ! <direction>            -- (in)     direction (ik), fft_forward = +1_ik or fft_inverse = -1_ik
  ! <sequence>             -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <sequence>             -- (out)    cdft (rk array of length = 2_ik*<length>), <sequence> = [..., fr_i, fi_i, ...]
  ! void    fft_radix_eight_(int* length, int* direction, double* sequence) ;
  interface
    module subroutine fft_radix_eight_(length, direction, sequence) &
      bind(c, name = "fft_radix_eight_")
      integer(ik), intent(in) :: length
      integer(ik), intent(in) :: direction
      real(rk), dimension(2_ik*length), intent(inout) :: sequence
    end subroutine fft_radix_eight_
  end interface
  public :: fft_radix_eight_
  ! ############################################################################################################################# !
  ! (takuya ooura) complex discrete fourier transform
  ! (subroutine) fft_radix_eight__(<length>, <direction>, <sequence>, <ip>, <work>)
  ! <length>               -- (in)     length (ik)
  ! <direction>            -- (in)     direction (ik), fft_forward = +1_ik or fft_inverse = -1_ik
  ! <sequence>             -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <sequence>             -- (out)    cdft (rk array of length = 2_ik*<length>), <sequence> = [..., fr_i, fi_i, ...]
  ! <ip>                   -- (in)     ffrft bit data
  ! <work>                 -- (in)     ffrft trig data
  interface
    module subroutine fft_radix_eight__(length, direction, sequence, ip, work)
      integer(ik), intent(in) :: length
      integer(ik), intent(in) :: direction
      real(rk), dimension(2_ik*length), intent(inout) :: sequence
      integer(ik), dimension(0_ik : 1_ik+int(sqrt(real(length/2_ik, rk)), ik)), intent(in) :: ip
      real(rk), dimension(0_ik : length/2_ik - 1_ik), intent(in) :: work
    end subroutine fft_radix_eight__
  end interface
  public :: fft_radix_eight__
  ! ############################################################################################################################# !
  ! compute data table (memorization)
  ! (subroutine) compute_table_(<length>, <pad>)
  ! <length>               -- (in)     length (ik)
  ! <pad>                  -- (in)     padded length (ik)
  ! void    compute_table_(int* length, int* pad) ;
  interface
    module subroutine compute_table_(length, pad) &
      bind(c, name = "compute_table_")
      integer(ik), intent(in) :: length
      integer(ik), intent(in) :: pad
    end subroutine compute_table_
  end interface
  public :: compute_table_
  ! ############################################################################################################################# !
  ! destroy data table
  ! (subroutine) destroy_table_()
  ! void    destroy_table_() ;
  interface
    module subroutine destroy_table_() &
      bind(c, name = "destroy_table_")
    end subroutine destroy_table_
  end interface
  public :: destroy_table_
  ! ############################################################################################################################# !
  ! svd
  ! ############################################################################################################################# !
  real(rk),    public, parameter :: svd_level              = 1.0e-12_rk          ! singular value threshold level (absolute value)
  integer(ik), public, parameter :: svd_row                = 2_ik**12_ik         ! max number of rows
  integer(ik), public, parameter :: svd_col                = 2_ik**12_ik         ! max number of rows
  integer(ik), public, parameter :: svd_basis              = 128_ik              ! max number of basis vectors in the implicitly restarted arnoldi process, optimal value (?)
  integer(ik), public, parameter :: svd_length             = 64_ik               ! length of arnoldi factorization
  integer(ik), public, parameter :: svd_loop               = 256_ik              ! max number of main loop iteraions
  ! ############################################################################################################################# !
  ! svd (dgesvd)
  ! (subroutine) svd_(<nr>, <nc>, <matrix>(<nr>, <nc>), <svd_list>(min(<nr>, <nc>)), <u_matrix>(<nr>, <nr>), <v_matrix>(<nc>, <nc>))
  ! <nr>                   -- (in)     number of rows (ik)
  ! <nc>                   -- (in)     number of cols (ik)
  ! <matrix>               -- (in)     input matrix(<nr>, <nc>) (rk)
  ! <svd_list>             -- (out)    list of singular values with size min(<nr>, <nc>) (rk)
  ! <u_matrix>             -- (out)    l-singular vectors (<nr>, <nr>) (rk)
  ! <v_matrix>             -- (out)    r-singular vectors (<nc>, <nc>) (rk)
  interface
    module subroutine svd_(nr, nc, matrix, svd_list, u_matrix, v_matrix)
      integer(ik), intent(in) :: nr
      integer(ik), intent(in) :: nc
      real(rk), dimension(nr, nc), intent(in) :: matrix
      real(rk), dimension(min(nr, nc)), intent(out) :: svd_list
      real(rk), dimension(nr, nr), intent(out) :: u_matrix
      real(rk), dimension(nc, nc), intent(out) :: v_matrix
    end subroutine svd_
  end interface
  public :: svd_
  ! ############################################################################################################################# !
  ! svd list (dgesvd)
  ! (subroutine) svd_list_(<nr>, <nc>, <matrix>(<nr>, <nc>), <svd_list>(min(<nr>, <nc>)))
  ! <nr>                   -- (in)     number of rows (ik)
  ! <nc>                   -- (in)     number of cols (ik)
  ! <matrix>               -- (in)     input matrix(<nr>, <nc>) (rk)
  ! <svd_list>             -- (out)    list of singular values with size min(<nr>, <nc>) (rk)
  interface
    module subroutine svd_list_(nr, nc, matrix, svd_list)
      integer(ik), intent(in) :: nr
      integer(ik), intent(in) :: nc
      real(rk), dimension(nr, nc), intent(in) :: matrix
      real(rk), dimension(min(nr, nc)), intent(out) :: svd_list
    end subroutine svd_list_
  end interface
  public :: svd_list_
  ! ############################################################################################################################# !
  ! truncated svd (arpack)
  ! svd_truncated_(<nr>, <nc>, <ns>, <matrix>(<nr>, <nc>), <list>(<ns>), <rvec>(<nc>, <ns>), <lvec>(<nr>, <ns>))
  ! <nr>                   -- (in)     number of rows (ik)
  ! <nc>                   -- (in)     number of cols (ik)
  ! <ns>                   -- (in)     number of singular values to keep
  ! <matrix>               -- (in)     input matrix(<nr>, <nc>) (rk)
  ! <list>                 -- (out)    list of singular values (<ns>) (rk)
  ! <rvec>                 -- (out)    l-singular vectors (<nc>, <ns>) (rk)
  ! <lvec>                 -- (out)    r-singular vectors (<nr>, <ns>) (rk)
  interface
    module subroutine svd_truncated_(nr, nc, ns, matrix, list, rvec, lvec)
      integer(ik), intent(in) :: nr
      integer(ik), intent(in) :: nc
      integer(ik), intent(in) :: ns
      real(rk), dimension(nr, nc), intent(in) :: matrix
      real(rk), dimension(ns), intent(out) :: list
      real(rk), dimension(nc, ns), intent(out) :: rvec
      real(rk), dimension(nr, ns), intent(out) :: lvec
    end subroutine svd_truncated_
  end interface
  public :: svd_truncated_
  ! ############################################################################################################################# !
  ! generic svd
  interface svd_
    module procedure svd_
    module procedure svd_list_
    module procedure svd_truncated_
  end interface svd_
  ! ############################################################################################################################# !
  ! process
  ! ############################################################################################################################# !
  ! convert input sequence (real)
  ! (subroutine) convert_real_(<length>, <r_part>, <sequence>)
  ! <length>               -- (in)     length (ik)
  ! <r_part>               -- (in)     input sequence r-part (rk array of length = <length>)
  ! <sequence>             -- (out)    sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...] and si_i=0.0_rk for all i
  ! void    convert_real_(int* length, double* r_part, double* sequence) ;
  interface
    module subroutine convert_real_(length, r_part, sequence) &
      bind(c, name = "convert_real_")
      integer(ik), intent(in) :: length
      real(rk), intent(in), dimension(length) :: r_part
      real(rk), intent(out), dimension(2_ik*length) :: sequence
    end subroutine convert_real_
  end interface
  public :: convert_real_
  ! ############################################################################################################################# !
  ! convert input sequence (complex)
  ! (subroutine) convert_complex_(<length>, <r_part>, <i_part>, <sequence>)
  ! <length>               -- (in)     length (ik)
  ! <r_part>               -- (in)     input sequence r-part (rk array of length = <length>)
  ! <i_part>               -- (in)     input sequence i-part (rk array of length = <length>)
  ! <sequence>             -- (out)    sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! void    convert_complex_(int* length, double* r_part, double* i_part, double* sequence) ;
  interface
    module subroutine convert_complex_(length, r_part, i_part, sequence) &
      bind(c, name = "convert_complex_")
      integer(ik), intent(in) :: length
      real(rk), intent(in), dimension(length) :: r_part
      real(rk), intent(in), dimension(length) :: i_part
      real(rk), intent(out), dimension(2_ik*length) :: sequence
    end subroutine convert_complex_
  end interface
  public :: convert_complex_
  ! ############################################################################################################################# !
  ! convert input sequence (real/complex)
  interface convert_
    module procedure convert_real_
    module procedure convert_complex_
  end interface
  public :: convert_
  ! ############################################################################################################################# !
  ! round up (round up to the next power of two)
  ! (function) round_up_(<number>)
  ! <number>               -- (in)     number (ik)
  ! <round_up>             -- (out)    next power of two number (ik)
  ! int     round_up_(int* number) ;
  interface
    module integer(ik) function round_up_(number) &
      bind(c, name = "round_up_")
      integer(ik), intent(in) :: number
    end function round_up_
  end interface
  public :: round_up_
  ! ############################################################################################################################# !
  ! zero padding
  ! (subroutine) pad_(<li>, <lo>, <input>, <output>)
  ! <li>                   -- (in)     input sequence length (ik)
  ! <lo>                   -- (in)     output sequence length (ik)
  ! <input>                -- (in)     input sequence (rk) of length = 2*<li>
  ! <output>               -- (in)     padded sequence (rk) of length = 2*<lo>
  ! void    pad_(int* linput, int* loutput, double* input, double* output) ;
  interface
    module subroutine pad_(li, lo, input, output) &
      bind(c, name = "pad_")
      integer(ik), intent(in) :: li
      integer(ik), intent(in) :: lo
      real(rk), dimension(2_ik*li), intent(in) :: input
      real(rk), dimension(2_ik*lo), intent(out) :: output
    end  subroutine pad_
  end interface
  public :: pad_
  ! ############################################################################################################################# !
  ! remove mean
  ! (subroutine) remove_mean_(<length>, <input>, <output> )
  ! <length>               -- (in)     input sequence length (ik)
  ! <input>                -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <output>               -- (out)    output sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! void    remove_mean_(int* length, double* input, double* output) ;
  interface
    module subroutine remove_mean_(length, input, output) &
      bind(c, name = "remove_mean_")
      integer(ik), intent(in) :: length
      real(rk), dimension(2_ik*length), intent(in) :: input
      real(rk), dimension(2_ik*length), intent(out) :: output
    end subroutine remove_mean_
  end interface
  public :: remove_mean_
  ! ############################################################################################################################# !
  ! remove window mean
  ! (subroutine) remove_window_mean_(<length>, <total>, <window>, <input>, <output> )
  ! <length>               -- (in)     input sequence length (ik)
  ! <total>                -- (in)     sum(window) (rk)
  ! <window>               -- (in)     window array (rk array of length = <length>)
  ! <input>                -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <output>               -- (out)    output sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! void    remove_window_mean_(int* length, double* total, double* window, double* input, double* output) ;
  interface
    module subroutine remove_window_mean_(length, total, window, input, output) &
      bind(c, name = "remove_window_mean_")
      integer(ik), intent(in) :: length
      real(rk), intent(in) :: total
      real(rk), dimension(length), intent(in) :: window
      real(rk), dimension(2_ik*length), intent(in) :: input
      real(rk), dimension(2_ik*length), intent(out) :: output
    end subroutine remove_window_mean_
  end interface
  public :: remove_window_mean_
  ! ############################################################################################################################# !
  ! apply window
  ! (subroutine) apply_window_(<length>, <window>, <input>, <output> )
  ! <length>               -- (in)     input sequence length (ik)
  ! <window>               -- (in)     window array (rk array of length = <length>)
  ! <input>                -- (in)     input sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <output>               -- (out)    output sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! void    apply_window_(int* length, double* window, double* input, double* output) ;
  interface
    module subroutine apply_window_(length, window, input, output) &
      bind(c, name = "apply_window_")
      integer(ik), intent(in) :: length
      real(rk), dimension(length), intent(in) :: window
      real(rk), dimension(2_ik*length), intent(in) :: input
      real(rk), dimension(2_ik*length), intent(out) :: output
    end subroutine apply_window_
  end interface
  public :: apply_window_
  ! ############################################################################################################################# !
  ! matrix (generate matrix from sequence)
  ! (subroutine) matrix_(<length>, <sequence>, <matrix>)
  ! <length>               -- (in)     input sequence length (ik)
  ! <sequence>             -- (in)     input sequence (rk)
  ! <matrix>               -- (out)    matrix (<length>/2+1, <length>/2) (rk)
  interface
    module subroutine matrix_(length, sequence, matrix)
      integer(ik), intent(in) :: length
      real(rk), dimension(length), intent(in) :: sequence
      real(rk), dimension(length/2_ik+1_ik, length/2_ik), intent(out) :: matrix
    end subroutine matrix_
  end interface
  public :: matrix_
  ! ############################################################################################################################# !
  ! sequence (row) (generate sequence from matrix using 1st and last rows)
  ! (subroutine) sequence_row_(<length>, <sequence>, <matrix>)
  ! <length>               -- (in)     input sequence length (ik)
  ! <sequence>             -- (out)    input sequence (rk)
  ! <matrix>               -- (in)     matrix (<length>/2+1, <length>/2) (rk)
  interface
    module subroutine sequence_row_(length, sequence, matrix)
      integer(ik), intent(in) :: length
      real(rk), dimension(length), intent(out) :: sequence
      real(rk), dimension(length/2_ik+1_ik, length/2_ik), intent(in) :: matrix
    end subroutine sequence_row_
  end interface
  public :: sequence_row_
  ! ############################################################################################################################# !
  ! sequence (sum) (generate sequence from matrix using sums of skew diagonals)
  ! (subroutine) sequence_row_(<length>, <sequence>, <matrix>)
  ! <length>               -- (in)     input sequence length (ik)
  ! <sequence>             -- (out)    input sequence (rk)
  ! <matrix>               -- (in)     matrix (<length>/2+1, <length>/2) (rk)
  interface
    module subroutine sequence_sum_(length, sequence, matrix)
      integer(ik), intent(in) :: length
      real(rk), dimension(length), intent(out) :: sequence
      real(rk), dimension(length/2_ik+1_ik, length/2_ik), intent(in) :: matrix
    end subroutine sequence_sum_
  end interface
  public :: sequence_sum_
  ! ############################################################################################################################# !
  ! filter
  ! (subroutine) filter(<length>, <sequence>, <limit>)
  ! <length>               -- (in)     length (ik)
  ! <sequence>             -- (inout)  sequence (rk array of length = <length>)
  ! <limit>                -- (in)     number of singular values to keep (ik)
  ! <svd_list>             -- (out)    list of singular values
  ! void    filter_(int* length, double* sequence, int* limit, double* svd_list) ;
  interface
    module subroutine filter_(length, sequence, limit, svd_list) &
      bind(c, name = "filter_")
      integer(ik), intent(in) :: length
      real(rk), dimension(length), intent(inout) :: sequence
      integer(ik), intent(in) :: limit
      real(rk), dimension(limit), intent(out) :: svd_list
    end subroutine filter_
  end interface
  public :: filter_
  ! ############################################################################################################################# !
  ! peakdetect
  ! ############################################################################################################################# !
  integer(ik), public, parameter :: peak_width             = 2_ik                ! peak width
  real(rk),    public, parameter :: peak_level             = -10.0_rk            ! peak threshold level
  ! ############################################################################################################################# !
  ! peak list
  ! (subroutine) peak_list_(<length>, <sequence>, <peak_list>)
  ! peak_width             -- (global) peak width (ik)
  ! peak_level             -- (global) peak level threshold (rk)
  ! <length>               -- (in)     sequence length (ik)
  ! <sequence>             -- (in)     sequence (rk array of length = <length>)
  ! <peak_list>            -- (out)    peak list (ik array of length = <length>), value of one correspond to peak location
  interface
    module subroutine peak_list_(length, sequence, peak_list)
      integer(ik), intent(in) :: length
      real(rk), dimension(length), intent(in) :: sequence
      integer(ik), dimension(length), intent(out) :: peak_list
    end subroutine peak_list_
  end interface
  public :: peak_list_
  ! ############################################################################################################################# !
  ! total number of peaks
  ! (function) peak_count_(<length>, <peak_list>)
  ! <length>               -- (in)     sequence length (ik)
  ! <peak_list>            -- (in)     peak list (ik array of length <length>)
  ! <peak_count_>          -- (out)    total number of peaks (ik)
  interface
    module integer(ik) function peak_count_(length, peak_list)
      integer(ik), intent(in) :: length
      integer(ik), dimension(length), intent(in) :: peak_list
    end function peak_count_
  end interface
  public :: peak_count_
  ! ############################################################################################################################# !
  ! detect several peaks (list of ordered peak positions)
  ! (subroutine) peak_detect_(<length>, <sequence>, <peak_length>, <peak_ordered_list>)
  ! <length>               -- (in)     sequence length (ik)
  ! <sequence>             -- (in)     sequence (rk array of length = <length>)
  ! <peak_length>          -- (in)     number of peaks to find (ik)
  ! <peak_ordered_list>    -- (out)    peak positions (ik array of length = <peak_length>)
  interface
    module subroutine peak_detect_(length, sequence, peak_length, peak_ordered_list)
      integer(ik), intent(in) :: length
      real(rk), dimension(length), intent(in) :: sequence
      integer(ik), intent(in) :: peak_length
      integer(ik), dimension(peak_length), intent(out) :: peak_ordered_list
    end subroutine peak_detect_
  end interface
  public :: peak_detect_
  ! ############################################################################################################################# !
  ! peak (ranked)
  ! (function) peak_(<length>, <sequence>, <peak_id>)
  ! <length>               -- (in)     sequence length (ik)
  ! <sequence>             -- (in)     sequence (rk array of length <length>)
  ! <peak_id>              -- (in)     peak rank (ik)
  ! <peak_>                -- (out)    peak position (ik)
  ! int     peak_(int* length, double* sequence, int* id) ;
  interface
    module integer(ik) function peak_(length, sequence, peak_id) &
      bind(c, name = "peak_")
      integer(ik), intent(in) :: length
      real(rk), dimension(length), intent(in) :: sequence
      integer(ik), intent(in) :: peak_id
    end function peak_
  end  interface
  public :: peak_
  ! ############################################################################################################################# !
  ! window
  ! ############################################################################################################################# !
  ! window (cosine)
  ! (subroutine) window_cos_(<length>, <order>, <sequence>)
  ! <length>               -- (in)     sequence length (ik)
  ! <order>                -- (in)     window order (ik)
  ! <window>               -- (out)    window (rk array of length = <length>)
  ! void    window_cos_(int* length, int* order, double* window) ;
  interface
    module subroutine window_cos_(length, order, window) &
      bind(c, name = "window_cos_")
      integer(ik), intent(in) :: length
      integer(ik), intent(in) :: order
      real(rk), intent(out), dimension(length) :: window
    end subroutine window_cos_
  end interface
  public :: window_cos_
  ! ############################################################################################################################# !
  ! window (cosine) (generic, fractional order)
  ! (subroutine) window_cos_generic_(<length>, <order>, <sequence>)
  ! <length>               -- (in)     sequence length (ik)
  ! <order>                -- (in)     window order (rk)
  ! <window>               -- (out)    window (rk array of length = <length>)
  ! void    window_cos_generic_(int* length, double* order, double* window) ;
  interface
    module subroutine window_cos_generic_(length, order, window)&
      bind(c, name = "window_cos_generic_")
      integer(ik), intent(in) :: length
      real(rk), intent(in) :: order
      real(rk), intent(out), dimension(length) :: window
    end subroutine window_cos_generic_
  end interface
  public :: window_cos_generic_
  ! ############################################################################################################################# !
  ! generic window
  interface window_
    module procedure window_cos_
    module procedure window_cos_generic_
  end interface
  public :: window_
  ! ############################################################################################################################# !
  ! window (kaiser)
  ! (subroutine) window_kaiser_(<length>, <order>, <sequence>)
  ! <length>               -- (in)     sequence length (ik)
  ! <parameter>            -- (in)     window order (rk)
  ! <window>               -- (out)    window (rk array of length = <length>)
  ! void    window_kaiser_(int* length, double* order, double* window) ;
  interface
    module subroutine window_kaiser_(length, order, window) &
      bind(c, name = "window_kaiser_")
      integer(ik), intent(in) :: length
      real(rk), intent(in) :: order
      real(rk), intent(out), dimension(length) :: window
    end subroutine window_kaiser_
  end interface
  public :: window_kaiser_
  ! ############################################################################################################################# !
  ! frequency
  ! ############################################################################################################################# !
  integer(ik), public, parameter :: flag_real              = 0_ik                ! signal flag (real)
  integer(ik), public, parameter :: flag_complex           = 1_ik                ! signal flag (complex)
  integer(ik), public, parameter :: frequency_fft          = 0_ik                ! fft
  integer(ik), public, parameter :: frequency_ffrft        = 1_ik                ! ffrft
  integer(ik), public, parameter :: frequency_parabola     = 2_ik                ! parabola
  integer(ik), public, parameter :: frequency_parabola_fit = 3_ik                ! parabola fit
  integer(ik), public, parameter :: frequency_search       = 4_ik                ! maximum search
  integer(ik), public, parameter :: parabola_fit_length    = 4_ik                ! number of parabola fit points
  integer(ik), public, parameter :: search_limit           = 128_ik              ! search limit
  real(rk)   , public, parameter :: search_tolerance       = epsilon             ! search tolerance
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
  interface
    module real(rk) function frequency_initial_(range_min, range_max, peak, length, pad, sequence) &
      bind(c, name = "frequency_initial_")
      real(rk), intent(in) :: range_min
      real(rk), intent(in) :: range_max
      integer(ik), intent(in) :: peak
      integer(ik), intent(in):: length
      integer(ik), intent(in):: pad
      real(rk), intent(in), dimension(2_ik*length) :: sequence
    end function frequency_initial_
  end interface
  public :: frequency_initial_
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
  interface
    module real(rk) function frequency_initial__(range_min, range_max, peak, length, pad, sequence) &
      bind(c, name = "frequency_initial__")
      real(rk), intent(in) :: range_min
      real(rk), intent(in) :: range_max
      integer(ik), intent(in) :: peak
      integer(ik), intent(in):: length
      integer(ik), intent(in):: pad
      real(rk), intent(in), dimension(2_ik*length) :: sequence
    end function frequency_initial__
  end interface
  public :: frequency_initial__
  ! ############################################################################################################################# !
  ! refine frequency estimation (ffrft)
  ! (function) frequency_refine_(<method>, <length>, <sequence>, <initial>)
  ! <method>               -- (in)     method
  ! <length>               -- (in)     sequence length (ik)
  ! <sequence>             -- (in)     input (processed) sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <initial>              -- (in)     initial frequency guess (rk)
  ! <frequency_refine_>    -- (out)    refined frequency estimation (rk)
  ! double  frequency_refine_(int* method, int* length, double* sequence, double* initial) ;
  interface
    module real(rk) function frequency_refine_(method, length, sequence, initial) &
      bind(c, name = "frequency_refine_")
      integer(ik), intent(in):: method
      integer(ik), intent(in):: length
      real(rk), intent(in), dimension(2_ik*length) :: sequence
      real(rk), intent(in) :: initial
    end function frequency_refine_
  end interface
  public :: frequency_refine_
  ! ############################################################################################################################# !
  ! refine frequency estimation (ffrft) (memorization)
  ! (function) frequency_refine_(<method>, <length>, <sequence>, <initial>)
  ! <method>               -- (in)     method
  ! <length>               -- (in)     sequence length (ik)
  ! <sequence>             -- (in)     input (processed) sequence (rk array of length = 2_ik*<length>), <sequence> = [..., sr_i, si_i, ...]
  ! <initial>              -- (in)     initial frequency guess (rk)
  ! <frequency_refine_>    -- (out)    refined frequency estimation (rk)
  ! double  frequency_refine__(int* method, int* length, double* sequence, double* initial) ;
  interface
    module real(rk) function frequency_refine__(method, length, sequence, initial) &
      bind(c, name = "frequency_refine__")
      integer(ik), intent(in):: method
      integer(ik), intent(in):: length
      real(rk), intent(in), dimension(2_ik*length) :: sequence
      real(rk), intent(in) :: initial
    end function frequency_refine__
  end interface
  public :: frequency_refine__
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
  interface
    module real(rk) function binary_amplitude_(flag, length, total, window, sequence, initial) &
      bind(c, name = "binary_amplitude_")
      integer(ik), intent(in):: flag
      integer(ik), intent(in):: length
      real(rk), intent(in) :: total
      real(rk), intent(in), dimension(length) :: window
      real(rk), intent(in), dimension(2_ik*length) :: sequence
      real(rk), intent(in) :: initial
    end function
  end interface
  public :: binary_amplitude_
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
  interface
    module real(rk) function golden_amplitude_(flag, length, total, window, sequence, initial) &
      bind(c, name = "golden_amplitude_")
      integer(ik), intent(in):: flag
      integer(ik), intent(in):: length
      real(rk), intent(in) :: total
      real(rk), intent(in), dimension(length) :: window
      real(rk), intent(in), dimension(2_ik*length) :: sequence
      real(rk), intent(in) :: initial
    end function
  end interface
  public :: golden_amplitude_
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
  interface
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
    end function frequency_
  end interface
  public :: frequency_
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
  interface
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
    end function frequency__
  end interface
  public :: frequency__
  ! ############################################################################################################################# !
  ! decomposition
  ! ############################################################################################################################# !
  integer(ik), public, parameter :: decomposition_subtract = 0_ik                ! decomposition by iterative subtraction
  integer(ik), public, parameter :: decomposition_peak     = 1_ik                ! decomposition by peaks
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
  interface
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
    end subroutine amplitude_
  end interface
  public :: amplitude_
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
  interface
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
    end subroutine decomposition_
  end interface
  public :: decomposition_
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
  interface
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
    end subroutine decomposition__
  end interface
  public :: decomposition__
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
  interface
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
    end subroutine frequency_list_
  end interface
  public :: frequency_list_
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
  interface
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
    end subroutine frequency_list__
  end interface
  public :: frequency_list__
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
  interface
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
    end subroutine amplitude_list_
  end interface
  public :: amplitude_list_
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
  interface
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
    end subroutine frequency_correction_
  end interface
  public :: frequency_correction_
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
  interface
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
    end subroutine frequency_correction__
  end interface
  public :: frequency_correction__
  ! ############################################################################################################################# !
  ! optimization
  ! ############################################################################################################################# !
  ! least squares (svd)
  ! (subroutine) least_squares_(<nr>, <nc>, <matrix>(<nr>, <nc>), <vector>(<nr>), <solution>(<nc>))
  ! <nr>                   -- (in)     number of rows (ik)
  ! <nc>                   -- (in)     number of cols (ik)
  ! <matrix>               -- (in)     input data matrix (<nr>, <nc>) (rk)
  ! <vector>               -- (in)     input vector (<nr>) (rk)
  ! <solution>             -- (out)    ls solution (<nc>) (rk)
  interface
    module subroutine least_squares_(nr, nc, matrix, vector, solution)
      integer(ik), intent(in) :: nr
      integer(ik), intent(in) :: nc
      real(rk), dimension(nr, nc), intent(in) :: matrix
      real(rk), dimension(nr), intent(in) :: vector
      real(rk), dimension(nc), intent(out) :: solution
    end subroutine least_squares_
  end interface
  public :: least_squares_
  ! ############################################################################################################################# !
  ! fit (harmonic signal)
  ! (subroutine) fit_(<length>, <sequence>, <loop>, <frequency>, <mean>, <cos_amp>, <sin_amp>, <error>)
  ! <length>               -- (in)     sequence length (ik), power of two
  ! <sequence>             -- (in)     input sequence (rk array of length = <length>)
  ! <loop>                 -- (in)     number of harmonics (ik)
  ! <frequency>            -- (in)     frequency array (rk array of length = <loop>)
  ! <mean>                 -- (out)    mean value
  ! <cos_amp>              -- (out)    cos amplitude array (rk array of length = <loop>)
  ! <sin_amp>              -- (out)    sin amplitude array (rk array of length = <loop>)
  ! <error>                -- (out)    error
  ! void    fit_(int* length, double* sequence, int* loop, double* frequency, double* mean, double* cos_amp, double* sin_amp, double* error) ;
  interface
    module subroutine fit_(length, sequence, loop, frequency, mean, cos_amp, sin_amp, error) &
      bind(c, name = "fit_")
      integer(ik), intent(in) :: length
      real(rk), dimension(length), intent(in) :: sequence
      integer(ik), intent(in) :: loop
      real(rk), dimension(loop), intent(in) :: frequency
      real(rk), intent(out) :: mean
      real(rk), dimension(loop), intent(out) :: cos_amp
      real(rk), dimension(loop), intent(out) :: sin_amp
      real(rk), intent(out) :: error
    end subroutine fit_
  end interface
  public :: fit_
  ! ############################################################################################################################# !
  ! fit (parabola) y = a*x**2 + b*x + c
  ! (subroutine) fit_parabola_(<length>, <x>, <y>, <a>, <b>, <c>, <maximum>)
  ! <length>               -- (in)     sequence length (ik), power of two
  ! <x>                    -- (in)     x (rk array of length = <length>)
  ! <y>                    -- (in)     y (rk array of length = <length>)
  ! <a>                    -- (out)    a (rk)
  ! <b>                    -- (out)    b (rk)
  ! <c>                    -- (out)    c (rk)
  ! <maximum>              -- (out)    maximum (minimum) position (rk)
  ! void    fit_parabola_(int* length, double* x, double* y, double* a, double* b, double* c, double* maximum) ;
  interface
    module subroutine fit_parabola_(length, x, y, a, b, c, maximum) &
      bind(c, name = "fit_parabola_")
      integer(ik), intent(in) :: length
      real(rk), dimension(length), intent(in) :: x
      real(rk), dimension(length), intent(in) :: y
      real(rk), intent(out) :: a
      real(rk), intent(out) :: b
      real(rk), intent(out) :: c
      real(rk), intent(out) :: maximum
    end subroutine fit_parabola_
  end interface
  public :: fit_parabola_
  ! ############################################################################################################################# !
  ! binary search maximization
  ! (function) binary_(<fun>, <guess>, <interval>, <limit>, <tolerance>)
  ! <fun>                  -- (in)     function to maximize (rk) -> (rk)
  ! <guess>                -- (in)     initial guess value (rk)
  ! <interval>             -- (in)     search interval (rk), guess is in the midle
  ! <limit>                -- (in)     maximum number of iterations (ik)
  ! <tolerance>            -- (in)     maximum tolerance (rk)
  ! <binary_>              -- (out)    maximum position
  interface
    module real(rk) function binary_(fun, guess, interval, limit, tolerance)
      interface
        real(rk) function fun(arg)
          import :: rk
          real(rk), intent(in) :: arg
        end function fun
      end interface
      real(rk), intent(in) :: guess
      real(rk), intent(in) :: interval
      integer(ik), intent(in) :: limit
      real(rk), intent(in) :: tolerance
    end function binary_
  end interface
  public :: binary_
  ! ############################################################################################################################# !
  ! golden search maximization
  ! (function) golden_(<fun>, <guess>, <interval>, <limit>, <tolerance>)
  ! <fun>                  -- (in)     function to maximize (rk) -> (rk)
  ! <guess>                -- (in)     initial guess value (rk)
  ! <interval>             -- (in)     search interval (rk), guess is in the midle
  ! <limit>                -- (in)     maximum number of iterations (ik)
  ! <tolerance>            -- (in)     maximum tolerance (rk)
  ! <golden_>              -- (out)    maximum position
  interface
    module real(rk) function golden_(fun, guess, interval, limit, tolerance)
      interface
        real(rk) function fun(arg)
          import :: rk
          real(rk), intent(in) :: arg
        end function fun
      end interface
      real(rk), intent(in) :: guess
      real(rk), intent(in) :: interval
      integer(ik), intent(in) :: limit
      real(rk), intent(in) :: tolerance
    end function golden_
  end interface
  public :: golden_
  ! ############################################################################################################################# !
end module signal
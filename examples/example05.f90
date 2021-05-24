! example-04: harmonic signal decomposition (least squares)
program example

  use :: signal

  implicit none

  real(rk), parameter :: a0 = 0.35_rk

  real(rk), parameter :: f1 = 1.0_rk*0.123456789_rk
  real(rk), parameter :: a1 = 1.0_rk
  real(rk), parameter :: b1 = 0.5_rk

  real(rk), parameter :: f2 = 2.0_rk*0.123456789_rk
  real(rk), parameter :: a2 = 0.05_rk
  real(rk), parameter :: b2 = 0.01_rk

  real(rk), parameter :: f3 = 3.0_rk*0.123456789_rk
  real(rk), parameter :: a3 = 0.001_rk
  real(rk), parameter :: b3 = 0.005_rk

  integer(ik), parameter :: length = 256_rk
  real(rk), dimension(length) :: signal

  integer(ik) :: limit
  real(rk) :: mean
  real(rk), dimension(:), allocatable :: frequency
  real(rk), dimension(:), allocatable :: cos_amp
  real(rk), dimension(:), allocatable :: sin_amp
  real(rk) :: error

  integer :: i

  signal = a0
  do i = 1_ik, length, 1_ik
    signal(i) = &
      a1*cos(two_pi*f1*real(i, rk)) + b1*sin(two_pi*f1*real(i, rk)) + &
      a2*cos(two_pi*f2*real(i, rk)) + b2*sin(two_pi*f2*real(i, rk)) + &
      a3*cos(two_pi*f3*real(i, rk)) + b3*sin(two_pi*f3*real(i, rk))
  end do

  ! fit (one harmonic)
  limit = 1_ik
  allocate(frequency(limit))
  allocate(cos_amp(limit))
  allocate(sin_amp(limit))
  frequency = [f1,f2]
  call fit_(length, signal, limit, frequency, mean, cos_amp, sin_amp, error)
  write(*,'(1e32.16)') error
  write(*,'(1e32.16)') mean
  do i = 1_ik, limit,1_ik
    write(*,'(2e32.16)') cos_amp(i), sin_amp(i)
  end do
  write(*,*)
  deallocate(frequency)
  deallocate(cos_amp)
  deallocate(sin_amp)

  ! fit (two harmonics)
  limit = 2_ik
  allocate(frequency(limit))
  allocate(cos_amp(limit))
  allocate(sin_amp(limit))
  frequency = [f1,f2]
  call fit_(length, signal, limit, frequency, mean, cos_amp, sin_amp, error)
  write(*,'(1e32.16)') error
  write(*,'(1e32.16)') mean
  do i = 1_ik, limit,1_ik
    write(*,'(2e32.16)') cos_amp(i), sin_amp(i)
  end do
  write(*,*)
  deallocate(frequency)
  deallocate(cos_amp)
  deallocate(sin_amp)

  ! fit (three harmonics)
  limit = 3_ik
  allocate(frequency(limit))
  allocate(cos_amp(limit))
  allocate(sin_amp(limit))
  frequency = [f1,f2,f3]
  call fit_(length, signal, limit, frequency, mean, cos_amp, sin_amp, error)
  write(*,'(1e32.16)') error
  write(*,'(1e32.16)') mean
  do i = 1_ik, limit,1_ik
    write(*,'(2e32.16)') cos_amp(i), sin_amp(i)
  end do
  write(*,*)
  deallocate(frequency)
  deallocate(cos_amp)
  deallocate(sin_amp)

end program example
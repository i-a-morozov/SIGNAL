! example-04: svd and least squares
program example

  use :: signal

  implicit none

  integer(ik) :: i

  ! svd
  block
    real(rk) :: m(2_ik, 2_ik), s(2_ik), u(2_ik, 2_ik), v(2_ik, 2_ik)
    m = real(reshape([1,1,2,2], shape(m)), rk)
    write(*, *) "input matrix"
    do i = 1_ik, 2_ik, 1_ik
      write(*, *) m(i,:)
    end do
    write(*, *)
    call svd_(2_ik, 2_ik, m, s, u, v)
    write(*, *) "singular values (svd_)"
    write(*, *) s
    write(*, *)
    m = 0.0_rk
    do i = 1_ik, size(s, kind = ik), 1_ik
      m(i, i) = s(i)
    end do
    m = matmul(u, matmul(m, transpose(v)))
    write(*, *) "output matrix"
    do i = 1_ik, size(s, kind = ik), 1_ik
      write(*, *) m(i, :)
    end do
    write(*, *)
    write(*, *) "singular values (svd_list_)"
    m = real(reshape([1,1,2,2], shape(m)), rk)
    s = 0.0_rk
    call svd_list_(2_ik, 2_ik, m, s)
    write(*, *) s
    write(*, *)
  end block

  ! svd
  block
    real(rk) :: m(4_ik, 5_ik), s(4_ik), u(4_ik, 4_ik), v(5_ik, 5_ik)
    m = real(reshape([1,4,7,8,2,5,8,9,3,6,9,10,4,7,10,11,5,8,11,12], shape(m)), rk)
    write(*, *) "input matrix"
    do i = 1_ik, 4_ik, 1_ik
      write(*, *) m(i, :)
    end do
    write(*, *)
    call svd_(4_ik, 5_ik, m, s, u, v)
    write(*, *) "singular values (svd_)"
    write(*, *) s
    write(*, *)
    m = 0.0_rk
    do i = 1_ik, size(s, kind = ik), 1_ik
      m(i, i) = s(i)
    end do
    m = matmul(u, matmul(m, transpose(v)))
    write(*, *) "output matrix"
    do i = 1_ik, size(s, kind = ik), 1_ik
      write(*, *) m(i, :)
    end do
    write(*, *)
    write(*, *) "singular values (svd_list_)"
    m = real(reshape([1,4,7,8,2,5,8,9,3,6,9,10,4,7,10,11,5,8,11,12], shape(m)), rk)
    s = 0.0_rk
    call svd_list_(4_ik, 5_ik, m, s)
    write(*, *) s
    write(*, *)
  end block

  ! least squares
  block
    real(rk) :: a(3_ik, 2_ik), b(3_ik), x(2_ik)
    write(*, *) "least_squares_"
    a = real(reshape([1,1,1,1,2,3], shape(a)), rk)
    b = real([7,7,8], rk)
    call least_squares_(3_ik, 2_ik, a, b, x)
    write(*, *) x
    write(*, *)
  end block

  ! least squares
  block
    real(rk) :: a(4_ik, 3_ik), b(4_ik), x(3_ik)
    write(*, *) "least_squares_"
    a = real(reshape([1,4,7,10,2,5,8,11,3,6,9,12], shape(a)), rk)
    b = real([1,2,4,8], rk)
    call least_squares_(4_ik, 3_ik, a, b, x)
    write(*, *) x
    write(*, *)
  end block

end program example

! (* wm *)
! m1 = {1,1,2,2} ;
! m1 = n[transpose[partition[m1,2]]] ;
! diagonal[first[rest[singularvaluedecomposition[m1]]]]
! m2 = {1,4,7,8,2,5,8,9,3,6,9,10,4,7,10,11,5,8,11,12} ;
! m2 = n[transpose[partition[m2,4]]] ;
! diagonal[first[rest[singularvaluedecomposition[m2]]]]
! a1 = {1,1,1,1,2,3} ;
! a1 = n[transpose[partition[a1,3]]] ;
! b1 = {7,7,8} ;
! leastsquares[a1,b1]
! a2 = {1,4,7,10,2,5,8,11,3,6,9,12} ;
! a2 = n[transpose[partition[a2,4]]] ;
! b2 = {1,2,4,8} ;
! leastsquares[a2,b2]
! (* {3.16228, 0.} *)
! (* {34.1299, 2.26956, 0., 0.} *)
! (* {6.33333, 0.5} *)
! (* {0.872222, 0.255556, -0.361111} *)
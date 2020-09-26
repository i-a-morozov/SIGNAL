! EXAMPLE-04: SVD AND LEAST SQUARES
PROGRAM EXAMPLE

  USE :: SIGNAL

  IMPLICIT NONE

  INTEGER(IK) :: I

  ! SVD
  BLOCK
    REAL(RK) :: M(2_IK, 2_IK), S(2_IK), U(2_IK, 2_IK), V(2_IK, 2_IK)
    M = REAL(RESHAPE([1,1,2,2], SHAPE(M)), RK)
    WRITE(*, *) "INPUT MATRIX"
    DO I = 1_IK, 2_IK, 1_IK
      WRITE(*, *) M(I,:)
    END DO
    WRITE(*, *)
    CALL SVD_(2_IK, 2_IK, M, S, U, V)
    WRITE(*, *) "SINGULAR VALUES (SVD_)"
    WRITE(*, *) S
    WRITE(*, *)
    M = 0.0_RK
    DO I = 1_IK, SIZE(S, KIND = IK), 1_IK
      M(I, I) = S(I)
    END DO
    M = MATMUL(U, MATMUL(M, TRANSPOSE(V)))
    WRITE(*, *) "OUTPUT MATRIX"
    DO I = 1_IK, SIZE(S, KIND = IK), 1_IK
      WRITE(*, *) M(I, :)
    END DO
    WRITE(*, *)
    WRITE(*, *) "SINGULAR VALUES (SVD_LIST_)"
    M = REAL(RESHAPE([1,1,2,2], SHAPE(M)), RK)
    S = 0.0_RK
    CALL SVD_LIST_(2_IK, 2_IK, M, S)
    WRITE(*, *) S
    WRITE(*, *)
  END BLOCK

  ! SVD
  BLOCK
    REAL(RK) :: M(4_IK, 5_IK), S(4_IK), U(4_IK, 4_IK), V(5_IK, 5_IK)
    M = REAL(RESHAPE([1,4,7,8,2,5,8,9,3,6,9,10,4,7,10,11,5,8,11,12], SHAPE(M)), RK)
    WRITE(*, *) "INPUT MATRIX"
    DO I = 1_IK, 4_IK, 1_IK
      WRITE(*, *) M(I, :)
    END DO
    WRITE(*, *)
    CALL SVD_(4_IK, 5_IK, M, S, U, V)
    WRITE(*, *) "SINGULAR VALUES (SVD_)"
    WRITE(*, *) S
    WRITE(*, *)
    M = 0.0_RK
    DO I = 1_IK, SIZE(S, KIND = IK), 1_IK
      M(I, I) = S(I)
    END DO
    M = MATMUL(U, MATMUL(M, TRANSPOSE(V)))
    WRITE(*, *) "OUTPUT MATRIX"
    DO I = 1_IK, SIZE(S, KIND = IK), 1_IK
      WRITE(*, *) M(I, :)
    END DO
    WRITE(*, *)
    WRITE(*, *) "SINGULAR VALUES (SVD_LIST_)"
    M = REAL(RESHAPE([1,4,7,8,2,5,8,9,3,6,9,10,4,7,10,11,5,8,11,12], SHAPE(M)), RK)
    S = 0.0_RK
    CALL SVD_LIST_(4_IK, 5_IK, M, S)
    WRITE(*, *) S
    WRITE(*, *)
  END BLOCK

  ! LEAST SQUARES
  BLOCK
    REAL(RK) :: A(3_IK, 2_IK), B(3_IK), X(2_IK)
    WRITE(*, *) "LEAST_SQUARES_"
    A = REAL(RESHAPE([1,1,1,1,2,3], SHAPE(A)), RK)
    B = REAL([7,7,8], RK)
    CALL LEAST_SQUARES_(3_IK, 2_IK, A, B, X)
    WRITE(*, *) X
    WRITE(*, *)
  END BLOCK

  ! LEAST SQUARES
  BLOCK
    REAL(RK) :: A(4_IK, 3_IK), B(4_IK), X(3_IK)
    WRITE(*, *) "LEAST_SQUARES_"
    A = REAL(RESHAPE([1,4,7,10,2,5,8,11,3,6,9,12], SHAPE(A)), RK)
    B = REAL([1,2,4,8], RK)
    CALL LEAST_SQUARES_(4_IK, 3_IK, A, B, X)
    WRITE(*, *) X
    WRITE(*, *)
  END BLOCK

END PROGRAM EXAMPLE

! (* WM *)
! m1 = {1,1,2,2} ;
! m1 = N[Transpose[Partition[m1,2]]] ;
! Diagonal[First[Rest[SingularValueDecomposition[m1]]]]
! m2 = {1,4,7,8,2,5,8,9,3,6,9,10,4,7,10,11,5,8,11,12} ;
! m2 = N[Transpose[Partition[m2,4]]] ;
! Diagonal[First[Rest[SingularValueDecomposition[m2]]]]
! a1 = {1,1,1,1,2,3} ;
! a1 = N[Transpose[Partition[a1,3]]] ;
! b1 = {7,7,8} ;
! LeastSquares[a1,b1]
! a2 = {1,4,7,10,2,5,8,11,3,6,9,12} ;
! a2 = N[Transpose[Partition[a2,4]]] ;
! b2 = {1,2,4,8} ;
! LeastSquares[a2,b2]
! (* {3.16228, 0.} *)
! (* {34.1299, 2.26956, 0., 0.} *)
! (* {6.33333, 0.5} *)
! (* {0.872222, 0.255556, -0.361111} *)
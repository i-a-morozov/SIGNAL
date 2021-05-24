#include "signal.inc"

submodule (signal) svd
  implicit none
  contains
  ! ############################################################################################################################# !
  ! svd (dgesvd)
  ! (subroutine) svd_(<nr>, <nc>, <matrix>(<nr>, <nc>), <svd_list>(min(<nr>, <nc>)), <u_matrix>(<nr>, <nr>), <v_matrix>(<nc>, <nc>))
  ! <nr>                   -- (in)     number of rows (ik)
  ! <nc>                   -- (in)     number of cols (ik)
  ! <matrix>               -- (in)     input matrix(<nr>, <nc>) (rk)
  ! <svd_list>             -- (out)    list of singular values with size min(<nr>, <nc>) (rk)
  ! <u_matrix>             -- (out)    l-singular vectors (<nr>, <nr>) (rk)
  ! <v_matrix>             -- (out)    r-singular vectors (<nc>, <nc>) (rk)
  module subroutine svd_(nr, nc, matrix, svd_list, u_matrix, v_matrix)
    integer(ik), intent(in) :: nr
    integer(ik), intent(in) :: nc
    real(rk), dimension(nr, nc), intent(in) :: matrix
    real(rk), dimension(min(nr, nc)), intent(out) :: svd_list
    real(rk), dimension(nr, nr), intent(out) :: u_matrix
    real(rk), dimension(nc, nc), intent(out) :: v_matrix
    real(rk), dimension(2_ik*max(1_ik, 3_ik*min(int(nr), int(nc))+max(int(nr), int(nc)), 5_ik*min(int(nr), int(nc)))) :: work
    integer :: work_size
    integer :: info
    work_size = size(work)
    call dgesvd('a','a',nr,nc,matrix,nr,svd_list,u_matrix,nr,v_matrix,nc,work,work_size,info)
    v_matrix = transpose(v_matrix)
    where (svd_list < svd_level) svd_list = 0.0_rk
  end subroutine svd_
  ! ############################################################################################################################# !
  ! svd list (dgesvd)
  ! (subroutine) svd_list_(<nr>, <nc>, <matrix>(<nr>, <nc>), <svd_list>(min(<nr>, <nc>)))
  ! <nr>                   -- (in)     number of rows (ik)
  ! <nc>                   -- (in)     number of cols (ik)
  ! <matrix>               -- (in)     input matrix(<nr>, <nc>) (rk)
  ! <svd_list>             -- (out)    list of singular values with size min(<nr>, <nc>) (rk)
  module subroutine svd_list_(nr, nc, matrix, svd_list)
    integer(ik), intent(in) :: nr
    integer(ik), intent(in) :: nc
    real(rk), dimension(nr, nc), intent(in) :: matrix
    real(rk), dimension(min(nr, nc)), intent(out) :: svd_list
    real(rk), dimension(2*max(1, 3*min(int(nr), int(nc))+max(int(nr), int(nc)), 5*min(int(nr), int(nc)))) :: work
    integer :: work_size
    integer :: info
    real(rk), dimension(nr, nr) :: u_matrix
    real(rk), dimension(nc, nc) :: v_matrix
    work_size = size(work)
    call dgesvd('n','n',nr,nc,matrix,nr,svd_list,u_matrix,nr,v_matrix,nc,work,work_size,info)
    where (svd_list < svd_level) svd_list = 0.0_rk
  end subroutine svd_list_
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
  module subroutine svd_truncated_(nr, nc, ns, matrix, list, rvec, lvec)
    integer(ik), intent(in) :: nr
    integer(ik), intent(in) :: nc
    integer(ik), intent(in) :: ns
    real(rk), dimension(nr, nc), intent(in) :: matrix
    real(rk), dimension(ns), intent(out) :: list
    real(rk), dimension(nc, ns), intent(out) :: rvec
    real(rk), dimension(nr, ns), intent(out) :: lvec
    real(rk), dimension(ns, ns) :: diag
    real(rk), dimension(nc, ns) :: copy
    integer(ik) :: i
    real(rk), dimension(svd_col, svd_basis) :: v = 0.0_rk
    real(rk), dimension(svd_basis*(svd_basis+8_ik)) :: wl
    real(rk), dimension(3_ik*svd_col) :: wd
    real(rk), dimension(svd_basis,2_ik) :: s = 0.0_rk
    real(rk), dimension(svd_col) :: id
    real(rk), dimension(svd_row) :: ax
    integer(ik) :: par(11_ik), pnt(11_ik)
    logical :: select(svd_basis)
    character(len=1) :: mat
    character(len=2) :: mode
    integer(ik) :: ido, nev, ncv, work, info, ierr, nconv
    real(rk) :: tol, sigma
    ! arpack
    nev    = ns                    ! # of singular values to compute, nev < n
    ncv    = min(svd_length, nc)    ! length of arnoldi factorization
    mat    = 'i'                   ! op = a^t.a
    mode   = 'lm'                  ! compute largest (magnitude) singular values, also la
    work   = ncv*(ncv+8_ik)        ! work array size
    tol    = 0.0_rk                ! tolerance
    info   = 0_ik                  ! initial error code
    ido    = 0_ik                  ! reverse communication parameter
    par(1) = 1_ik                  ! shift, 0/1
    par(3) = svd_loop              ! max number of iteraions
    par(7) = 1_ik                  ! mode
    ! main loop
    main : do
      call dsaupd(ido,mat,nc,mode,nev,tol,id,ncv,v,svd_col,par,pnt,wd,wl,work,info)
      if (.not.(ido .eq. -1_ik .or. ido .eq. 1_ik)) exit
      call dgemv('n',nr,nc,1.0_rk,matrix,nr,wd(pnt(1)),1_ik,0.0_rk,ax,1_ik)
      call dgemv('t',nr,nc,1.0_rk,matrix,nr,ax,1_ik,0.0_rk,wd(pnt(2)),1_ik)
    end do main
    ! right singular verctors
    call dseupd (.true.,'all',select,s,v,svd_col,sigma,mat,nc,mode,nev,tol,id,ncv,v,svd_col,par,pnt,wd,wl,work,ierr)
    ! scale, reverse and set singular values
    nconv =  par(5)
    list = sqrt(abs(s(nev:1_ik:-1_ik,1_ik)))
    ! trim
    where (list < svd_level) list = 0.0_rk
    ! reverse and set right vectors
    rvec(:,1_ik:nev:1_ik) = v(1_ik:nc,nev:1_ik:-1_ik)
    ! compute left vector, u = a.v.s^-1, direct computation
    diag = 0.0_rk
    do i = 1_ik , nev, 1_ik
      if(list(i) > svd_level) diag(i, i) = 1.0_rk/list(i)
    end do
    call dgemm('n','n',nc,ns,ns,1.0_rk,rvec,nc,diag,ns,0.0_rk,copy,nc)
    call dgemm('n','n',nr,ns,nc,1.0_rk,matrix,nr,copy,nc,0.0_rk,lvec,nr)
  end subroutine svd_truncated_
  ! ############################################################################################################################# !
end submodule svd
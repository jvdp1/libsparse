module modtest_crs
#if (_DP==0)
 use, intrinsic :: iso_fortran_env, only: int64, wp => real32
#else
 use, intrinsic :: iso_fortran_env, only: int64, wp => real64
#endif
 use testdrive, only: new_unittest, unittest_type, error_type, check
 use modsparse, only: coosparse, crssparse, assignment(=)
 use modtest_common, only: tol_wp, verbose, ia, ja, a, aspsd, addval => addval_coo&
                       , getmat, matcheck, printmat&
                       , iaspsdf, jaspsdf, aspsdf&
                       , iaspsdf1, jaspsdf1, aspsdf1
 implicit none
 private

 public :: collect_crs

 integer, parameter :: sparse_unit = 556
 
contains

!> Collect all exported unit tests
subroutine collect_crs(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  !how to test: getordering, print, printsquare, save, sort

  testsuite = [ &
    new_unittest("crs constructor", test_constructor) &
    , new_unittest("crs add", test_add) &
    , new_unittest("crs add nel", test_add_nel) &
    , new_unittest("crs add lupper", test_add_lupper) &
    , new_unittest("crs add vect", test_add_vect) &
    , new_unittest("crs ncol add", test_ncol_add) &
    , new_unittest("crs ncol add nel", test_ncol_add_nel) &
    , new_unittest("crs ncol add lupper", test_ncol_add_lupper) &
    , new_unittest("crs test external", test_external) &
    , new_unittest("crs test harville", test_harville) &
    , new_unittest("crs cg", test_cg) &
    , new_unittest("crs chol", test_chol) &
    , new_unittest("crs diag vect", test_diag_vect) &
    , new_unittest("crs diag vect lupper", test_diag_vect_lupper) &
    , new_unittest("crs ncol diag vect", test_ncol_diag_vect) &
    , new_unittest("crs ncol diag vect lupper", test_ncol_diag_vect_lupper) &
    , new_unittest("crs get", test_get) &
    , new_unittest("crs get nel", test_get_nel) &
    , new_unittest("crs get lupper", test_get_lupper) &
    , new_unittest("crs get lupper_sym", test_get_lupper_sym) &
    , new_unittest("crs get permutation", test_getpermutation) &
    , new_unittest("crs ncol get", test_ncol_get) &
    , new_unittest("crs ncol get nel", test_ncol_get_nel) &
    , new_unittest("crs ncol get lupper", test_ncol_get_lupper) &
    , new_unittest("crs ichol", test_ichol) &
    , new_unittest("crs ichol_failed", test_ichol_failed, should_fail = .true.) &
    , new_unittest("crs isolve", test_isolve) &
    , new_unittest("crs ldlt", test_ldlt) &
    , new_unittest("crs multbyv_n_n_n", test_multbyv_n_n_n) &
    , new_unittest("crs multbyv_n_n_y", test_multbyv_n_n_y) &
    , new_unittest("crs multbyv_n_y_n", test_multbyv_n_y_n) &
    , new_unittest("crs multbyv_n_y_y", test_multbyv_n_y_y) &
    , new_unittest("crs multbyv_y_n_n", test_multbyv_y_n_n) &
    , new_unittest("crs multbyv_y_n_y", test_multbyv_y_n_y) &
    , new_unittest("crs multbyv_y_y_n", test_multbyv_y_y_n) &
    , new_unittest("crs multbyv_y_y_y", test_multbyv_y_y_y) &
    , new_unittest("crs multbyv_sym", test_multbyv_sym) &
    , new_unittest("crs multbyv_sym1", test_multbyv_sym1) &
    , new_unittest("crs multbymat_n_n_n", test_multbymat_n_n_n) &
    , new_unittest("crs multbymat_n_n_y", test_multbymat_n_n_y) &
    , new_unittest("crs multbymat_n_y_n", test_multbymat_n_y_n) &
    , new_unittest("crs multbymat_n_y_y", test_multbymat_n_y_y) &
    , new_unittest("crs multbymat_y_n_n", test_multbymat_y_n_n) &
    , new_unittest("crs multbymat_y_n_y", test_multbymat_y_n_y) &
    , new_unittest("crs multbymat_y_y_n", test_multbymat_y_y_n) &
    , new_unittest("crs multbymat_y_y_y", test_multbymat_y_y_y) &
    , new_unittest("crs multbymat_sym", test_multbymat_sym) &
    , new_unittest("crs multbymat_sym1", test_multbymat_sym1) &
    , new_unittest("crs nonzero", test_nonzero) &
    , new_unittest("crs nonzero_sym", test_nonzero_sym) &
    , new_unittest("crs scale", test_scale) &
#if (_PARDISO==1)
    , new_unittest("crs solve", test_solve_vector) &
    , new_unittest("crs solve_perm", test_solve_vector_perm) &
    , new_unittest("crs solve_arr", test_solve_array) &
    , new_unittest("crs solve_arr_perm", test_solve_array_perm) &
#else
    , new_unittest("crs solve", test_solve_vector, should_fail = .true.) &
    , new_unittest("crs solve_perm", test_solve_vector_perm, should_fail = .true.) &
    , new_unittest("crs solve_arr", test_solve_array, should_fail = .true.) &
    , new_unittest("crs solve_arr_perm", test_solve_array_perm, should_fail = .true.) &
#endif
    , new_unittest("crs solveldlt_s", test_solveldlt_s) &
    , new_unittest("crs solveldlt_vector", test_solveldlt_vector) &
    , new_unittest("crs solveldlt_array", test_solveldlt_array) &
    , new_unittest("crs spainv", test_spainv) &
    , new_unittest("crs spainv_failed", test_spainv_failed, should_fail = .true.) &
    , new_unittest("crs spainv_1", test_spainv_1) &
    , new_unittest("crs spainv_spsd", test_spainv_spsd) &
    , new_unittest("crs spainv_spsd_1", test_spainv_spsd_1) &
    , new_unittest("crs spainv_spsd_2", test_spainv_spsd_2) &
    , new_unittest("crs submatrix_off_full_full", test_submatrix_off_full_full) &
    , new_unittest("crs submatrix_off_full_upper", test_submatrix_off_full_upper) &
    , new_unittest("crs submatrix_off_upper_full", test_submatrix_off_upper_full) &
    , new_unittest("crs submatrix_full_full", test_submatrix_full_full) &
    , new_unittest("crs submatrix_full_upper", test_submatrix_full_upper) &
    , new_unittest("crs submatrix_upper_full", test_submatrix_upper_full) &
    , new_unittest("crs submatrix_dense_upper_upper", test_submatrix_dense_upper_upper) &
    , new_unittest("crs submatrix_dense_upper_upper_sym", test_submatrix_dense_upper_upper_sym) &
    , new_unittest("crs submatrix_dense_upper_full", test_submatrix_dense_upper_full) &
    ]

  !to check: diag_mat

end subroutine collect_crs

!CONSTRUCTOR + GETDIM + ISSQUARE
subroutine test_constructor(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 integer, parameter :: nel = 9
 type(crssparse) :: crs

 crs = crssparse(nrow, nel)
 call check(error, crs%getdim(1), nrow, more = 'getdim(1)')
 call check(error, crs%getdim(2), nrow, more = 'getdim(2)')
 call check(error, crs%getdim(3) < 0, more = 'getdim(3)')
 call check(error, crs%issquare(), nrow.eq.nrow, 'isquare')
 if (allocated(error)) return

 crs = crssparse(nrow, nel, lupper = .true.)
 call check(error, crs%getdim(1), nrow, more = 'getdim(1)')
 call check(error, crs%getdim(2), nrow, more = 'getdim(2)')
 call check(error, crs%issquare(), nrow.eq.nrow, 'isquare')
 if (allocated(error)) return

 crs = crssparse(nrow, nel, n = ncol)
 call check(error, crs%getdim(1), nrow, more = 'getdim(1)')
 call check(error, crs%getdim(2), ncol, more = 'getdim(2)')
 call check(error, crs%issquare(), nrow.eq.ncol, 'isquare')
 if (allocated(error)) return

 crs = crssparse(nrow, nel, n = ncol, lupper = .true.)
 call check(error, crs%getdim(1), nrow, more = 'getdim(1)')
 call check(error, crs%getdim(2), ncol, more = 'getdim(2)')
 call check(error, crs%issquare(), nrow.eq.ncol, 'isquare')
 if (allocated(error)) return

 crs = crssparse(nrow, nel, n = ncol, lupper = .true., unlog = sparse_unit)
 call check(error, crs%getdim(1), nrow, more = 'getdim(1)')
 call check(error, crs%getdim(2), ncol, more = 'getdim(2)')
 call check(error, crs%issquare(), nrow.eq.ncol, 'isquare')
 if (allocated(error)) return


end subroutine

!ADD
subroutine test_add(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, ncol)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, unlog = sparse_unit)
 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo
 call getija_crs(crs, iat, jat, at, mat)

 call check(error, size(iat), size(pack(ia, lvalid)), more = 'size(iat)')
 call check(error, size(jat), size(pack(ja, lvalid)), more = 'size(jat)')
 call check(error, size(at), size(pack(a, lvalid)), more = 'size(at)')
 if (allocated(error)) return

 call check(error, all(iat == pack(ia, lvalid)), 'ia')
 call check(error, all(jat == pack(ja, lvalid)), 'ja')
 call check(error, all(at == pack(a, lvalid)), 'a')
 if (allocated(error)) return

 if(verbose)call printmat(mat)

end subroutine

subroutine test_add_vect(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: ncrs = 3
 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol

 integer :: i
 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, ncol)
 type(coosparse) :: coo
 type(crssparse), allocatable :: crs(:)

 coo = coosparse(nrow, unlog = sparse_unit)
 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 allocate(crs(ncrs))

 do i = 1, ncrs
  crs(i) = coo
  call getija_crs(crs(i), iat, jat, at, mat)

  call check(error, size(iat), size(pack(ia, lvalid)), more = 'size(iat)')
  call check(error, size(jat), size(pack(ja, lvalid)), more = 'size(jat)')
  call check(error, size(at), size(pack(a, lvalid)), more = 'size(at)')
  if (allocated(error)) return

  call check(error, all(iat == pack(ia, lvalid)), 'ia')
  call check(error, all(jat == pack(ja, lvalid)), 'ja')
  call check(error, all(at == pack(a, lvalid)), 'a')
  if (allocated(error)) return

  if(verbose)call printmat(mat)
 end do

 deallocate(crs)

end subroutine

subroutine test_add_nel(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, ncol)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, nel = 4_int64, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo
 call getija_crs(crs, iat, jat, at, mat)

 call check(error, size(iat), size(pack(ia, lvalid)), more = 'size(iat)')
 call check(error, size(jat), size(pack(ja, lvalid)), more = 'size(jat)')
 call check(error, size(at), size(pack(a, lvalid)), more = 'size(at)')
 if (allocated(error)) return

 call check(error, all(iat == pack(ia, lvalid)), 'ia')
 call check(error, all(jat == pack(ja, lvalid)), 'ja')
 call check(error, all(at == pack(a, lvalid)), 'a')
 if (allocated(error)) return

 if(verbose)call printmat(mat)

end subroutine

subroutine test_add_lupper(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia .le. ja

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, ncol)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo
 call getija_crs(crs, iat, jat, at, mat)

 call check(error, size(iat), size(pack(ia, lvalid)), more = 'size(iat)')
 call check(error, size(jat), size(pack(ja, lvalid)), more = 'size(jat)')
 call check(error, size(at), size(pack(a, lvalid)), more = 'size(at)')
 if (allocated(error)) return

 call check(error, all(iat == pack(ia, lvalid)), 'ia')
 call check(error, all(jat == pack(ja, lvalid)), 'ja')
 call check(error, all(at == pack(a, lvalid)), 'a')
 if (allocated(error)) return

 if(verbose)call printmat(mat)

end subroutine

subroutine test_ncol_add(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, ncol)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a)

 crs = coo
 call getija_crs(crs, iat, jat, at, mat)

 call check(error, size(iat), size(pack(ia, lvalid)), more = 'size(iat)')
 call check(error, size(jat), size(pack(ja, lvalid)), more = 'size(jat)')
 call check(error, size(at), size(pack(a, lvalid)), more = 'size(at)')
 if (allocated(error)) return

 call check(error, all(iat == pack(ia, lvalid)), 'ia')
 call check(error, all(jat == pack(ja, lvalid)), 'ja')
 call check(error, all(at == pack(a, lvalid)), 'a')
 if (allocated(error)) return

 if(verbose)call printmat(mat)

end subroutine

subroutine test_ncol_add_nel(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, ncol)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, n = ncol, nel = 2_int64,  unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a)

 crs = coo
 call getija_crs(crs, iat, jat, at, mat)

 call check(error, size(iat), size(pack(ia, lvalid)), more = 'size(iat)')
 call check(error, size(jat), size(pack(ja, lvalid)), more = 'size(jat)')
 call check(error, size(at), size(pack(a, lvalid)), more = 'size(at)')
 if (allocated(error)) return

 call check(error, all(iat == pack(ia, lvalid)), 'ia')
 call check(error, all(jat == pack(ja, lvalid)), 'ja')
 call check(error, all(at == pack(a, lvalid)), 'a')
 if (allocated(error)) return

 if(verbose)call printmat(mat)

end subroutine

subroutine test_ncol_add_lupper(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 logical, parameter :: lvalid(size(ia)) =  ia.le.nrow .and. ja.le.ncol .and. ia.le.ja

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, ncol)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, n = ncol, lupper = .true.,  unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo
 call getija_crs(crs, iat, jat, at, mat)

 call check(error, size(iat), size(pack(ia, lvalid)), more = 'size(iat)')
 call check(error, size(jat), size(pack(ja, lvalid)), more = 'size(jat)')
 call check(error, size(at), size(pack(a, lvalid)), more = 'size(at)')
 if (allocated(error)) return

 call check(error, all(iat == pack(ia, lvalid)), 'ia')
 call check(error, all(jat == pack(ja, lvalid)), 'ja')
 call check(error, all(at == pack(a, lvalid)), 'a')
 if (allocated(error)) return

 if(verbose)call printmat(mat)

end subroutine

subroutine test_external(error)
  type(error_type), allocatable, intent(out) :: error
 
  integer, parameter :: nrow = 5
  integer, parameter :: ncol = 4
  logical, parameter :: lvalid(size(ia)) =  ia.le.nrow .and. ja.le.ncol .and. ia.le.ja
 
  integer, allocatable :: iat(:), jat(:), ia1(:), ja1(:), ia2(:), ja2(:)
  real(wp), allocatable :: at(:), a1(:), a2(:)
  real(wp) :: mat(nrow, ncol)
  integer :: nz
  type(coosparse) :: coo
  type(crssparse) :: crs1, crs2
 
  coo = coosparse(nrow, n = ncol, lupper = .true.,  unlog = sparse_unit)
  call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)
  crs1 = coo
  nz = crs1%nonzero()
  
  call crs1%get_rowptr(iat)
  call crs1%get_colval(jat)
  call crs1%get_nzval(at)

  crs2 = crssparse(nrow, nz, n = ncol)
  call crs2%external(iat, jat, at)

  call getija_crs(crs1, ia1, ja1, a1, mat)
  call getija_crs(crs2, ia2, ja2, a2, mat)

  call check(error, all(ia1 == ia2), 'ia')
  call check(error, all(ja1 == ja2), 'ja')
  call check(error, all(a1 == a2), 'a')
  if (allocated(error)) return
 
  if(verbose)call printmat(mat)
 
 end subroutine

!CG
subroutine test_cg(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i

 integer, parameter :: nrow = 6
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia.le.ja

 integer :: maxiter
 real(wp) :: tol
 real(wp) :: xmat(nrow, ncol)
 real(wp) :: mat_l(nrow, ncol)
 real(wp) :: y(ncol)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, aspsd)
 call coo%add(2, 2, 202._wp) !to set the matrix SPD

 crs = coo

 xmat = 0
 do i = 1, ncol
  tol = 1.e-6
  maxiter = nrow
  y = 0
  y(i) = 1
  call crs%cg(xmat(:,i), y, maxiter, tol)
 enddo

 !mat_l: expected result 
 mat_l = reshape([(merge(1._wp, 0._wp, i/nrow.eq.mod(i,nrow)), i = 0, nrow**2 - 1)], [nrow, nrow])
 where(getmat(crs) <= tol_wp) mat_l = 0._wp

 call check(error, all(abs(matmul(getmat(crs), xmat) - mat_l) < tol_wp), 'crs cg')

end subroutine

!CHOLESKY FACTOR
subroutine test_chol(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i

 integer, parameter :: nrow = 6
 integer, parameter :: ncol = nrow
 integer, parameter :: perm(nrow) = [(i, i =nrow, 1, -1)]

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, ncol)
 real(wp) :: matchol(nrow, ncol)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, n = ncol, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

! call addval(coo, coo%getdim(1), coo%getdim(2), [ia, 2], [ja, 2]&
!             , merge([a, 22._wp] + 1000, [a, 22._wp], [ia, 2] ==  [ja, 2]))
 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, aspsd)

 crs = coo

 call getija_crs(crs, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 !Complete Cholesky decomposition
 call check(error, crs%isdecomposed(), .false., 'Chol - wrong status 0')
 if(allocated(error))return

 call crs%setpermutation(perm)
 call crs%chol()

 call check(error, crs%isdecomposed(), .true., 'Chol - wrong status 1')
 if(allocated(error))return

 deallocate(iat); deallocate(jat); deallocate(at)
 call getija_crs(crs, iat, jat, at, matchol)
 
 call check(error, all(abs(mat(perm, perm) - & 
                       matmul(transpose(matchol), matchol)) < tol_wp)&
                 , 'Chol')
 if(allocated(error))return

 !Complete Cholesky decomposition
 crs = coo

 call check(error, crs%isdecomposed(), .false., 'Chol_permf - wrong status 2')
 if(allocated(error))return

 call crs%setpermutation(permf(nrow))
 call crs%chol()

 call check(error, crs%isdecomposed(), .true., 'Chol_permf - wrong status 3')
 if(allocated(error))return

 deallocate(iat); deallocate(jat); deallocate(at)
 call getija_crs(crs, iat, jat, at, matchol)
 
 call check(error, all(abs(mat(permf(nrow), permf(nrow)) - & 
                       matmul(transpose(matchol), matchol)) < tol_wp)&
                 , 'Chol_permf')

end subroutine

!DIAG VECT
subroutine test_diag_vect(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia.eq.ja
 integer, parameter :: sizediag =min(nrow, ncol)

 real(wp) :: diagcheck(sizediag)
 real(wp), allocatable :: diag(:)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, unlog = sparse_unit)
 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo

 diag = crs%diag()

 diagcheck = 0
 diagcheck(pack(ia, lvalid)) = pack(a, lvalid)

 call check(error, size(diag), sizediag, more = 'size(diag)')
 if (allocated(error)) return
 if (allocated(error)) return

 call check(error, all(diag == diagcheck), 'diag')

end subroutine

subroutine test_diag_vect_lupper(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia.eq.ja
 integer, parameter :: sizediag = min(nrow, ncol)

 real(wp) :: diagcheck(sizediag)
 real(wp), allocatable :: diag(:)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)
 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a)

 crs = coo

 diag = crs%diag()

 diagcheck = 0
 diagcheck(pack(ia, lvalid)) = pack(a, lvalid)

 call check(error, size(diag), sizediag, more = 'size(diag)')
 if (allocated(error)) return
 if (allocated(error)) return

 call check(error, all(diag == diagcheck), 'diag')

end subroutine

subroutine test_ncol_diag_vect(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia.eq.ja
 integer, parameter :: sizediag = min(nrow, ncol)

 real(wp) :: diagcheck(sizediag)
 real(wp), allocatable :: diag(:)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, n = ncol,  unlog = sparse_unit)
 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo

 diag = crs%diag()

 diagcheck = 0
 diagcheck(pack(ia, lvalid)) = pack(a, lvalid)

 call check(error, size(diag), sizediag, more = 'size(diag)')
 if (allocated(error)) return

 call check(error, all(diag == diagcheck), 'diag')

end subroutine

subroutine test_ncol_diag_vect_lupper(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia.eq.ja
 integer, parameter :: sizediag = min(nrow, ncol)

 real(wp) :: diagcheck(sizediag)
 real(wp), allocatable :: diag(:)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, n = ncol,  unlog = sparse_unit)
 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo

 diag = crs%diag()

 diagcheck = 0
 diagcheck(pack(ia, lvalid)) = pack(a, lvalid)

 call check(error, size(diag), sizediag, more = 'size(diag)')
 if (allocated(error)) return

 call check(error, all(diag == diagcheck), 'diag')

end subroutine

!GET
subroutine test_get(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol

 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, unlog = sparse_unit)
 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a)

 crs = coo

 call check(error, all(getmat(crs) == matcheck(nrow, ncol, ia, ja, a, lvalid)), 'get')

end subroutine

subroutine test_get_nel(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol

 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, nel = 4_int64, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo

 call check(error, all(getmat(crs) == matcheck(nrow, ncol, ia, ja, a, lvalid)), 'get')

end subroutine

subroutine test_get_lupper(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia .le. ja

 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo

 call check(error, all(getmat(crs) == matcheck(nrow, ncol, ia, ja, a, lvalid)), 'get')

end subroutine

subroutine test_get_lupper_sym(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia .le. ja

 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo

 call check(error, all(getmat(crs) == &
             matcheck(nrow, ncol, ia, ja, a, lvalid) &
             + transpose(matcheck(nrow, ncol, ia, ja, a, lvalid .and. ia.lt.ja)))&
             , 'get sym')

end subroutine

subroutine test_ncol_get(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol

 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo

 call check(error, all(getmat(crs) == matcheck(nrow, ncol, ia, ja, a, lvalid)), 'get')

end subroutine

subroutine test_ncol_get_nel(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol

 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, n = ncol, nel = 2_int64,  unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo

 call check(error, all(getmat(crs) == matcheck(nrow, ncol, ia, ja, a, lvalid)), 'get')

end subroutine

subroutine test_ncol_get_lupper(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 logical, parameter :: lvalid(size(ia)) =  ia.le.nrow .and. ja.le.ncol .and. ia.le.ja

 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, n = ncol, lupper = .true.,  unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo

 call check(error, all(getmat(crs) == matcheck(nrow, ncol, ia, ja, a, lvalid)), 'get')

end subroutine

!GET PERMUTATION
subroutine test_getpermutation(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i
 integer, parameter :: nrow = 6
 integer, parameter :: perm(nrow) = [(i, i = nrow, 1, -1)]

 integer, allocatable :: gperm(:)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, aspsd)

 crs = coo

 call crs%setpermutation(perm)

 call crs%getpermutation(gperm)

 call check(error, all(gperm.eq.perm), 'get permutation')
 if(allocated(error))return

end subroutine

!HARVILLE
subroutine test_harville(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 6

 real(wp), parameter :: diaginvsol(nrow) = [ real(wp):: 1.007295172270438E-002, 0.5_wp, 3.345347748825206E-003&
               , 2.546510312665334E-003, 1.987905182703515E-003, 1.664886113502529E-003]

 real(wp), allocatable :: diaginv(:)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), [ia, 2], [ja, 2], [aspsd, 2._wp])

 crs = coo

 call crs%harville(100, 10, diaginv)

 call check(error, all(abs(int(diaginv*10**7)/10._wp**7 - int(diaginvsol*10**7)/10._wp**7) < 2.5_wp*10._wp**(-5)), 'harville')
 return

end subroutine

!ICHOL
subroutine test_ichol(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i

 integer, parameter :: nrow = 6
 integer, parameter :: perm(nrow) = [(i, i = 1, nrow)]

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, nrow)
 real(wp) :: matchol(nrow, nrow)
 real(wp) :: mat_l(nrow, nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs
 type(crssparse) :: crschol
 logical :: lvalid(nrow, nrow)

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), [ia, 2], [ja, 2], [aspsd, 0._wp])

 !ICHOL
 !Cholesky
 crschol= coo

 call check(error, crschol%isdecomposed(), .false., 'ichol - status 0')
 if(allocated(error))return

 call crschol%setpermutation(perm)
 call crschol%chol()

 call check(error, crschol%isdecomposed(), .true., 'ichol - status 1')
 if(allocated(error))return
 
 call getija_crs(crschol, iat, jat, at, matchol)

 crs = coo

 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 call crs%setpermutation(perm)
 call crs%ichol()
 
 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat_l)

 lvalid = mat(perm, perm).ne.0

 call check(error, all(abs(pack(mat_l, lvalid) - pack(matchol, lvalid)) < tol_wp), 'ichol')
 if(allocated(error))return

 !ICHOL
 !Cholesky
 crschol= coo

 call check(error, crschol%isdecomposed(), .false., 'ichol_permf - status 0')
 if(allocated(error))return

 call crschol%setpermutation(permf(nrow))
 call crschol%chol()

 call check(error, crschol%isdecomposed(), .true., 'ichol_permf - status 1')
 if(allocated(error))return
 
 call getija_crs(crschol, iat, jat, at, matchol)

 crs = coo

 call crs%setpermutation(permf(nrow))
 call crs%ichol()
 
 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat_l)

 lvalid = mat(permf(nrow), permf(nrow)).ne.0

 call check(error, all(abs(pack(mat_l, lvalid) - pack(matchol, lvalid)) < tol_wp), 'ichol_permf')

end subroutine

subroutine test_ichol_failed(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i

 integer, parameter :: nrow = 6
 integer, parameter :: perm(nrow) = [(i, i = 1, nrow)]

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, nrow)
 real(wp) :: matchol(nrow, nrow)
 real(wp) :: mat_l(nrow, nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs
 type(crssparse) :: crschol

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), [ia, 2], [ja, 2], [aspsd, 0._wp])

 !Cholesky
 crschol= coo
 call crschol%setpermutation(perm)
 call crschol%chol()
 
 call getija_crs(crschol, iat, jat, at, matchol)

 !ICHOL
 crs = coo

 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 call crs%setpermutation(perm)
 call crs%ichol()
 
 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat_l)

 call check(error, all(abs(mat_l - matchol) < tol_wp), 'should_ichol_failed')

end subroutine

!ISOLVE
subroutine test_isolve(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i
 integer, parameter :: nrow = 6
 integer, parameter :: perm(nrow) = [(i, i = nrow, 1, -1)]

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, nrow)
 real(wp) :: mat_l(nrow, nrow)
 real(wp) :: mat_d(nrow, nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 !mat must be SPD
 call addval(coo, coo%getdim(1), coo%getdim(2), [ia, 2], [ja, 2]&
             , merge([a, 22._wp] + 1000, [a, 22._wp], [ia, 2] ==  [ja, 2]))

 crs = coo

 call getija_crs(crs, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 !Cholesky
 call check(error, crs%isdecomposed(), .false., 'isolve - status 0')
 if(allocated(error))return

 call crs%setpermutation(perm)
 call crs%chol()

 call check(error, crs%isdecomposed(), .true., 'isolve - status 1')
 if(allocated(error))return
 
 do i = 1, nrow
  call crs%isolve(mat_d(:,i), mat(:,i))
 enddo

 !mat_l: expected result 
 mat_l = reshape([(merge(1._wp, 0._wp, i/nrow.eq.mod(i,nrow)), i = 0, nrow**2 - 1)], [nrow, nrow])
 where(mat <= tol_wp) mat_l = 0._wp

 call check(error, all(abs(mat_d - mat_l) < tol_wp), 'isolve')
 if(allocated(error))return

 !Cholesky
 crs = coo

 call check(error, crs%isdecomposed(), .false., 'isolve_permf - status 0')
 if(allocated(error))return

 call crs%setpermutation(permf(nrow))
 call crs%chol()

 call check(error, crs%isdecomposed(), .true., 'isolve_permf - status 0')
 if(allocated(error))return
 
 do i = 1, nrow
  call crs%isolve(mat_d(:,i), mat(:,i))
 enddo

 !mat_l: expected result 
 mat_l = reshape([(merge(1._wp, 0._wp, i/nrow.eq.mod(i,nrow)), i = 0, nrow**2 - 1)], [nrow, nrow])
 where(mat <= tol_wp) mat_l = 0._wp

 call check(error, all(abs(mat_d - mat_l) < tol_wp), 'isolve_permf')

end subroutine

!LDLt
subroutine test_ldlt(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i
 integer, parameter :: nrow = 6
 integer, parameter :: ncol = nrow
 integer, parameter :: perm(nrow) = [(i, i = nrow, 1, -1)]

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, ncol)
 real(wp) :: mat_l(nrow, nrow)
 real(wp) :: mat_d(nrow, nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, n = ncol, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, aspsd)

 crs = coo

 call getija_crs(crs, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 !LDLt
 call crs%setpermutation(perm)
 call crs%getldlt()

 deallocate(iat); deallocate(jat); deallocate(at)
 call getija_crs(crs, iat, jat, at, mat_l)
 
 mat_d = 0
 do i = 1, nrow
  mat_d(i, i) = mat_l(i, i)
  mat_l(i, i) = 1
 enddo

 call check(error, all(abs(mat(perm, perm) - & 
                       matmul(transpose(mat_l), matmul(mat_d, mat_l))) < tol_wp)&
                 , 'LDLt')

 !LDLt
 crs = coo
 call crs%setpermutation(permf(nrow))
 call crs%getldlt()

 deallocate(iat); deallocate(jat); deallocate(at)
 call getija_crs(crs, iat, jat, at, mat_l)
 
 mat_d = 0
 do i = 1, nrow
  mat_d(i, i) = mat_l(i, i)
  mat_l(i, i) = 1
 enddo

 call check(error, all(abs(mat(permf(nrow), permf(nrow)) - & 
                       matmul(transpose(mat_l), matmul(mat_d, mat_l))) < tol_wp)&
                 , 'LDLt_permf')

end subroutine

!MULT BY VECT
subroutine test_multbyv_n_n_n(error)
 type(error_type), allocatable, intent(out) :: error

 call test_multbyv_gen(error, col = 'n', trans = 'n', upper = 'n')
end subroutine

subroutine test_multbyv_n_n_y(error)
 type(error_type), allocatable, intent(out) :: error

 call test_multbyv_gen(error, col = 'n', trans = 'n', upper = 'y')
end subroutine

subroutine test_multbyv_n_y_n(error)
 type(error_type), allocatable, intent(out) :: error

 call test_multbyv_gen(error, col = 'n', trans = 'y', upper = 'n')
end subroutine

subroutine test_multbyv_n_y_y(error)
 type(error_type), allocatable, intent(out) :: error

 call test_multbyv_gen(error, col = 'n', trans = 'y', upper = 'y')
end subroutine

subroutine test_multbyv_y_n_n(error)
 type(error_type), allocatable, intent(out) :: error

 call test_multbyv_gen(error, col = 'y', trans = 'n', upper = 'n')
end subroutine

subroutine test_multbyv_y_n_y(error)
 type(error_type), allocatable, intent(out) :: error

 call test_multbyv_gen(error, col = 'y', trans = 'n', upper = 'y')
end subroutine

subroutine test_multbyv_y_y_n(error)
 type(error_type), allocatable, intent(out) :: error

 call test_multbyv_gen(error, col = 'y', trans = 'y', upper = 'n')
end subroutine

subroutine test_multbyv_y_y_y(error)
 type(error_type), allocatable, intent(out) :: error

 call test_multbyv_gen(error, col = 'y', trans = 'y', upper = 'y')
end subroutine

subroutine test_multbyv_sym(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 real(wp), parameter :: alpha = 0.3_wp
 real(wp), parameter :: val = 10000._wp
 real(wp), parameter :: x(ncol) = [(real(i, wp), i = 1, ncol)]
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia.le.ja

 real(wp) :: y(nrow), ycheck(nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo

 y = 1._wp
 call crs%mult(alpha, 'n', x, val, y)

 ycheck = 1._wp
 ycheck = ycheck * val + alpha * matmul(&
   matcheck(nrow, ncol, ia, ja, a, lvalid) &
     + transpose(matcheck(nrow, ncol, ia, ja, a, lvalid .and. ia.lt.ja))&
   , x)
 
 call check(error, all(abs(y - ycheck) < tol_wp), 'multbyv_sym')

end subroutine

subroutine test_multbyv_sym1(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 real(wp), parameter :: alpha = 0.3_wp
 real(wp), parameter :: val = 10000._wp
 real(wp), parameter :: x(ncol) = [(real(i, wp), i = 1, ncol)]
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia.le.ja

 real(wp) :: y(nrow), ycheck(nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo
 call crs%setsymmetric()

 y = 1._wp
 call crs%mult(alpha, 'n', x, val, y)

 ycheck = 1._wp
 ycheck = ycheck * val + alpha * matmul(&
   matcheck(nrow, ncol, ia, ja, a, lvalid) &
     + transpose(matcheck(nrow, ncol, ia, ja, a, lvalid .and. ia.lt.ja))&
   , x)
 
 call check(error, all(abs(y - ycheck) < tol_wp), 'multbyv_sym')

end subroutine


subroutine test_multbyv_gen(error, col, trans, upper)
 character(len=*), intent(in) :: col
 character(len=*), intent(in) :: trans
 character(len=*), intent(in) :: upper
 type(error_type), allocatable, intent(out) :: error


 integer, parameter :: nrow = 5
 real(wp), parameter :: alpha = 0.3_wp
 real(wp), parameter :: val = 10000._wp

 integer :: i, ncol
 real(wp), allocatable :: x(:), y(:), ycheck(:)
 logical :: lvalid(size(ia))
 type(coosparse) :: coo
 type(crssparse) :: crs

 ncol = merge(4, nrow, col == 'y')

 if(upper == 'y')then
  lvalid = ia.le.nrow .and. ja.le.ncol .and. ia.le.ja
 else
  lvalid = ia.le.nrow .and. ja.le.ncol
 endif

 allocate(x(merge(nrow, ncol, trans == 'y')), source = &
                              [(real(i, wp), i = 1, merge(nrow, ncol, trans == 'y'))])
 allocate(y(merge(ncol, nrow, trans == 'y')), source = 1._wp)
 allocate(ycheck(merge(ncol, nrow, trans == 'y')), source = 1._wp)


 if(col == 'y')then
  if(upper == 'y')then
   !coo = coosparse(nrow, n = ncol, lupper = .true., unlog = sparse_unit)
   call coo%init(nrow, n = ncol, lupper = .true., unlog = sparse_unit)
  else
   !coo = coosparse(nrow, n = ncol, unlog = sparse_unit)
   call coo%init(nrow, n = ncol, unlog = sparse_unit)
  endif
 else
  if(upper == 'y')then
   !coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)
   call coo%init(nrow, lupper = .true., unlog = sparse_unit)
  else
   !coo = coosparse(nrow, unlog = sparse_unit)
   call coo%init(nrow, unlog = sparse_unit)
  endif
 endif

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo

 call crs%mult(alpha, merge('t', 'n', trans == 'y'), x, val, y)

 if(trans == 'y')then
  ycheck = ycheck * val + alpha * matmul(&
   transpose(matcheck(nrow, ncol, ia, ja, a, lvalid)) &
   , x)
 else
  ycheck = ycheck * val + alpha * matmul(&
   matcheck(nrow, ncol, ia, ja, a, lvalid) &
   , x)
 endif
 
 call check(error, all(abs(y - ycheck) < tol_wp)&
                 , 'multbyv_'//col//'_'//trans//'_'//upper)

end subroutine


!MULT BY MAT
subroutine test_multbymat_n_n_n(error)
 type(error_type), allocatable, intent(out) :: error

 call test_multbymat_gen(error, col = 'n', trans = 'n', upper = 'n')
end subroutine

subroutine test_multbymat_n_n_y(error)
 type(error_type), allocatable, intent(out) :: error

 call test_multbymat_gen(error, col = 'n', trans = 'n', upper = 'y')
end subroutine

subroutine test_multbymat_n_y_n(error)
 type(error_type), allocatable, intent(out) :: error

 call test_multbymat_gen(error, col = 'n', trans = 'y', upper = 'n')
end subroutine

subroutine test_multbymat_n_y_y(error)
 type(error_type), allocatable, intent(out) :: error

 call test_multbymat_gen(error, col = 'n', trans = 'y', upper = 'y')
end subroutine

subroutine test_multbymat_y_n_n(error)
 type(error_type), allocatable, intent(out) :: error

 call test_multbymat_gen(error, col = 'y', trans = 'n', upper = 'n')
end subroutine

subroutine test_multbymat_y_n_y(error)
 type(error_type), allocatable, intent(out) :: error

 call test_multbymat_gen(error, col = 'y', trans = 'n', upper = 'y')
end subroutine

subroutine test_multbymat_y_y_n(error)
 type(error_type), allocatable, intent(out) :: error

 call test_multbymat_gen(error, col = 'y', trans = 'y', upper = 'n')
end subroutine

subroutine test_multbymat_y_y_y(error)
 type(error_type), allocatable, intent(out) :: error

 call test_multbymat_gen(error, col = 'y', trans = 'y', upper = 'y')
end subroutine

subroutine test_multbymat_sym(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrhs = 3
 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 real(wp), parameter :: alpha = 0.3_wp
 real(wp), parameter :: val = 10000._wp
 real(wp), parameter :: x(ncol, nrhs) = 1._wp
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia.le.ja

 real(wp) :: y(nrow, nrhs), ycheck(nrow, nrhs)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo

 y = 1._wp
 call crs%mult(alpha, 'n', x, val, y)

 ycheck = 1._wp
 ycheck = ycheck * val + alpha * matmul(&
   matcheck(nrow, ncol, ia, ja, a, lvalid) &
   + transpose(matcheck(nrow, ncol, ia, ja, a, lvalid .and. ia.lt.ja))&
   , x)
 
 call check(error, all(abs(y - ycheck) < tol_wp), 'multbymat_sym')

end subroutine

subroutine test_multbymat_sym1(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrhs = 3
 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 real(wp), parameter :: alpha = 0.3_wp
 real(wp), parameter :: val = 10000._wp
 real(wp), parameter :: x(ncol, nrhs) = 1._wp
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia.le.ja

 real(wp) :: y(nrow, nrhs), ycheck(nrow, nrhs)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo
 call crs%setsymmetric()

 y = 1._wp
 call crs%mult(alpha, 'n', x, val, y)

 ycheck = 1._wp
 ycheck = ycheck * val + alpha * matmul(&
   matcheck(nrow, ncol, ia, ja, a, lvalid) &
   + transpose(matcheck(nrow, ncol, ia, ja, a, lvalid .and. ia.lt.ja))&
   , x)
 
 call check(error, all(abs(y - ycheck) < tol_wp), 'multbymat_sym')

end subroutine


subroutine test_multbymat_gen(error, col, trans, upper)
 character(len=*), intent(in) :: col
 character(len=*), intent(in) :: trans
 character(len=*), intent(in) :: upper
 type(error_type), allocatable, intent(out) :: error


 integer, parameter :: nrhs = 3
 integer, parameter :: nrow = 5
 real(wp), parameter :: alpha = 0.3_wp
 real(wp), parameter :: val = 10000._wp

 integer :: ncol
 real(wp), allocatable :: x(:,:), y(:,:), ycheck(:,:)
 logical :: lvalid(size(ia))
 type(coosparse) :: coo
 type(crssparse) :: crs

 ncol = merge(4, nrow, col == 'y')

 if(upper == 'y')then
  lvalid = ia.le.nrow .and. ja.le.ncol .and. ia.le.ja
 else
  lvalid = ia.le.nrow .and. ja.le.ncol
 endif

 allocate(x(merge(nrow, ncol, trans == 'y'), nrhs), source = 1._wp)
 allocate(y(merge(ncol, nrow, trans == 'y'), nrhs), source = 1._wp)
 allocate(ycheck(merge(ncol, nrow, trans == 'y'), nrhs), source = 1._wp)

 if(col == 'y')then
  if(upper == 'y')then
   !coo = coosparse(nrow, n = ncol, lupper = .true., unlog = sparse_unit)
   call coo%init(nrow, n = ncol, lupper = .true., unlog = sparse_unit)
  else
   !coo = coosparse(nrow, n = ncol, unlog = sparse_unit)
   call coo%init(nrow, n = ncol, unlog = sparse_unit)
  endif
 else
  if(upper == 'y')then
   !coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)
   call coo%init(nrow, lupper = .true., unlog = sparse_unit)
  else
   !coo = coosparse(nrow, unlog = sparse_unit)
   call coo%init(nrow, unlog = sparse_unit)
  endif
 endif

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)
 
 crs = coo

 call crs%mult(alpha, merge('t', 'n', trans == 'y'), x, val, y)

 if(trans == 'y')then
  ycheck = ycheck * val + alpha * matmul(&
   transpose(matcheck(nrow, ncol, ia, ja, a, lvalid)) &
   , x)
 else
  ycheck = ycheck * val + alpha * matmul(&
   matcheck(nrow, ncol, ia, ja, a, lvalid) &
   , x)
 endif
 
 call check(error, all(abs(y - ycheck) < tol_wp)&
                 , 'multbymat_'//col//'_'//trans//'_'//upper)

end subroutine


!NONZERO
subroutine test_nonzero(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4

 integer(int64) :: p(1)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo

 p = count(ia.le.nrow .and. ja.le.ncol .and. ja.ne.ncol .and. a.ne.0._wp) + nrow

 call check(error, crs%nonzero(), p(1), 'nonzero')

end subroutine

subroutine test_nonzero_sym(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow

 integer(int64) :: p(1)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, n = ncol, lupper = .true., unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 crs = coo

 p = count(ia.le.nrow .and. ja.le.ncol .and. ia.lt.ja .and. a.ne.0._wp) + nrow

 call check(error, crs%nonzero(), p(1), 'nonzero_sym')

end subroutine

!SCALE
subroutine test_scale(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol

 real(wp), parameter :: scalefact = 1000._wp
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a)

 crs = coo

 call crs%scale(scalefact)

 call check(error, all(getmat(crs) == scalefact * matcheck(nrow, ncol, ia, ja, a, lvalid)), 'scale')

end subroutine

!SOLVE using PARDISO VECTOR
subroutine test_solve_vector(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 6

 integer :: i
 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, nrow)
 real(wp) :: mat_l(nrow, nrow)
 real(wp) :: mat_d(nrow, nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, aspsd)

 crs = coo

 call getija_crs(crs, iat, jat, at, mat)

 !SOLVE
 do i = 1, nrow
  call crs%solve(mat_d(:,i), mat(:,i))
 enddo

 !mat_l: expected result 
 mat_l = reshape([(merge(1._wp, 0._wp, i/nrow.eq.mod(i,nrow)), i = 0, nrow**2 - 1)], [nrow, nrow])
 where(mat <= tol_wp) mat_l = 0._wp

 call check(error, all(abs(mat_d - mat_l) < tol_wp), 'Solve vector')

end subroutine

subroutine test_solve_vector_perm(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i
 integer, parameter :: nrow = 6
 integer, parameter :: perm(nrow) = [(i, i = nrow, 1, -1)]

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, nrow)
 real(wp) :: mat_l(nrow, nrow)
 real(wp) :: mat_d(nrow, nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, aspsd)

 crs = coo

 call getija_crs(crs, iat, jat, at, mat)

 !SOLVE
 call crs%setpermutation(perm)
 
 do i = 1, nrow
  call crs%solve(mat_d(:,i), mat(:,i))
 enddo

 !mat_l: expected result 
 mat_l = reshape([(merge(1._wp, 0._wp, i/nrow.eq.mod(i,nrow)), i = 0, nrow**2 - 1)], [nrow, nrow])
 where(mat <= tol_wp) mat_l = 0._wp

 call check(error, all(abs(mat_d - mat_l) < tol_wp), 'Solve vector perm')
 if(allocated(error))return

 !SOLVE
 crs = coo

 call crs%setpermutation(permf(nrow))
 
 do i = 1, nrow
  call crs%solve(mat_d(:,i), mat(:,i))
 enddo

 call check(error, all(abs(mat_d - mat_l) < tol_wp), 'Solve vector permf')


end subroutine

subroutine test_solve_array(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 6

 integer :: i
 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, nrow)
 real(wp) :: mat_l(nrow, nrow)
 real(wp) :: mat_d(nrow, nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, aspsd)

 crs = coo

 call getija_crs(crs, iat, jat, at, mat)

 !SOLVE
 call crs%solve(mat_d, mat)

 !mat_l: expected result 
 mat_l = reshape([(merge(1._wp, 0._wp, i/nrow.eq.mod(i,nrow)), i = 0, nrow**2 - 1)], [nrow, nrow])
 where(mat <= tol_wp) mat_l = 0._wp

 call check(error, all(abs(mat_d - mat_l) < tol_wp), 'Solve vector')

end subroutine

subroutine test_solve_array_perm(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i
 integer, parameter :: nrow = 6
 integer, parameter :: perm(nrow) = [(i, i = nrow, 1, -1)]

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, nrow)
 real(wp) :: mat_l(nrow, nrow)
 real(wp) :: mat_d(nrow, nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, aspsd)

 crs = coo

 call getija_crs(crs, iat, jat, at, mat)

 !SOLVE
 call crs%setpermutation(perm)
 
 call crs%solve(mat_d, mat)

 !mat_l: expected result 
 mat_l = reshape([(merge(1._wp, 0._wp, i/nrow.eq.mod(i,nrow)), i = 0, nrow**2 - 1)], [nrow, nrow])
 where(mat <= tol_wp) mat_l = 0._wp

 call check(error, all(abs(mat_d - mat_l) < tol_wp), 'Solve array perm')
 if(allocated(error))return

 !SOLVE
 crs = coo

 call crs%setpermutation(permf(nrow))
 
 call crs%solve(mat_d, mat)

 call check(error, all(abs(mat_d - mat_l) < tol_wp), 'Solve array permf')

end subroutine

!SOLVE LDLt
subroutine test_solveldlt_s(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i
 integer, parameter :: nrow = 6
 integer, parameter :: perm(nrow) = [(i, i = nrow, 1, -1)]

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, nrow)
 real(wp) :: mat_l(nrow, nrow)
 real(wp) :: mat_d(nrow, nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, aspsd)

 crs = coo

 call getija_crs(crs, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 !LDLt
 call crs%setpermutation(perm)
 call crs%getldlt()
 
 do i = 1, nrow
  call crs%solveldlt_s(mat_d(:,i), mat(:,i))
 enddo

 !mat_l: expected result 
 mat_l = reshape([(merge(1._wp, 0._wp, i/nrow.eq.mod(i,nrow)), i = 0, nrow**2 - 1)], [nrow, nrow])
 where(mat <= tol_wp) mat_l = 0._wp

 call check(error, all(abs(mat_d - mat_l) < tol_wp), 'Solve LDLt')
 if(allocated(error))return

 !LDLt
 crs = coo

 call crs%setpermutation(permf(nrow))
 call crs%getldlt()
 
 do i = 1, nrow
  call crs%solveldlt_s(mat_d(:,i), mat(:,i))
 enddo

 call check(error, all(abs(mat_d - mat_l) < tol_wp), 'Solve LDLt permf')

end subroutine

!SOLVE LDLt WITHOUT COMPUTING LDLt EXPLICITLY
subroutine test_solveldlt_vector(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i
 integer, parameter :: nrow = 6
 integer, parameter :: perm(nrow) = [(i, i = nrow, 1, -1)]

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, nrow)
 real(wp) :: mat_l(nrow, nrow)
 real(wp) :: mat_d(nrow, nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, aspsd)

 crs = coo

 call getija_crs(crs, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 !LDLt
 call crs%setpermutation(perm)

 do i = 1, nrow
  call crs%solveldlt(mat_d(:,i), mat(:,i))
 enddo

 !mat_l: expected result
 mat_l = reshape([(merge(1._wp, 0._wp, i/nrow.eq.mod(i,nrow)), i = 0, nrow**2 - 1)], [nrow, nrow])
 where(mat <= tol_wp) mat_l = 0._wp

 call check(error, all(abs(mat_d - mat_l) < tol_wp), 'Solve LDLt')
 if(allocated(error))return

 !LDLt
 crs = coo

 call crs%setpermutation(permf(nrow))

 do i = 1, nrow
  call crs%solveldlt(mat_d(:,i), mat(:,i))
 enddo

 call check(error, all(abs(mat_d - mat_l) < tol_wp), 'Solve LDLt permf')

end subroutine

subroutine test_solveldlt_array(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i
 integer, parameter :: nrow = 6
 integer, parameter :: perm(nrow) = [(i, i = nrow, 1, -1)]

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, nrow)
 real(wp) :: mat_l(nrow, nrow)
 real(wp) :: mat_d(nrow, nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, aspsd)

 crs = coo

 call getija_crs(crs, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 !LDLt
 call crs%setpermutation(perm)

 call crs%solveldlt(mat_d(:,:), mat(:,:))

 !mat_l: expected result
 mat_l = reshape([(merge(1._wp, 0._wp, i/nrow.eq.mod(i,nrow)), i = 0, nrow**2 - 1)], [nrow, nrow])
 where(mat <= tol_wp) mat_l = 0._wp

 call check(error, all(abs(mat_d - mat_l) < tol_wp), 'Solve LDLt')
 if(allocated(error))return

 !LDLt
 crs = coo

 call crs%setpermutation(permf(nrow))

 call crs%solveldlt(mat_d(:,:), mat(:,:))

 call check(error, all(abs(mat_d - mat_l) < tol_wp), 'Solve LDLt permf')

end subroutine

!GET SPARSE INVERSE
subroutine test_spainv(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i
 integer, parameter :: nrow = 6
 integer, parameter :: perm(nrow) = [(i, i = nrow, 1, -1)]

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: vect(nrow)
 real(wp) :: mat(nrow, nrow)
 real(wp) :: matinv(nrow, nrow)
 real(wp) :: mat_l(nrow, nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs
 type(crssparse) :: crschol

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), [ia, 2], [ja, 2], [aspsd, 2._wp])

 !Cholesky
 crschol= coo
 call crschol%setpermutation(perm)
 call crschol%chol()
 
 do i = 1, nrow
  vect = 0
  vect(i) = 1
  call crschol%isolve(matinv(:,i), vect) !matinv reference
 enddo

 !SPAINV
 crs = coo

 call getija_crs(crs, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 call crs%setpermutation(perm)
 call crs%spainv()
 
 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat_l)
 if(verbose)call printmat(mat_l)

 call check(error, all(abs(pack(mat_l, mat.ne.0) - pack(matinv, mat.ne.0)) < tol_wp), 'spainv')
 if(allocated(error))return

 !SPAINV
 crs = coo

 call crs%setpermutation(permf(nrow))
 call crs%spainv()
 
 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat_l)

 call check(error, all(abs(pack(mat_l, mat.ne.0) - pack(matinv, mat.ne.0)) < tol_wp), 'spainv_permf')

end subroutine

subroutine test_spainv_1(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i
 integer, parameter :: nrow = 6
 integer, parameter :: perm(nrow) = [(i, i = nrow, 1, -1)]

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: vect(nrow)
 real(wp) :: mat(nrow, nrow)
 real(wp) :: matinv(nrow, nrow)
 real(wp) :: mat_l(nrow, nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs
 type(crssparse) :: crschol

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, aspsd)

 !Cholesky
 crschol= coo

 call getija_crs(crschol, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 call crschol%setpermutation(perm)
 call crschol%chol()
 
 do i = 1, nrow
  vect = 0
  vect(i) = 1
  call crschol%isolve(matinv(:,i), vect) !matinv reference
 enddo

 !SPAINV
 crs = coo

 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 call crs%setpermutation(perm)
 call crs%spainv()
 
 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat_l)
 if(verbose)call printmat(mat_l)

 call check(error, all(abs(pack(mat_l, mat.ne.0) - pack(matinv, mat.ne.0)) < tol_wp), 'spainv_1')
 if(allocated(error))return

 !SPAINV
 crs = coo

 call crs%setpermutation(permf(nrow))
 call crs%spainv()
 
 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat_l)

 call check(error, all(abs(pack(mat_l, mat.ne.0) - pack(matinv, mat.ne.0)) < tol_wp), 'spainv_1_permf')

end subroutine

subroutine test_spainv_spsd(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i
 integer, parameter :: nrow = 6
 integer, parameter :: perm(nrow) = [(i, i = nrow, 1, -1)]

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: vect(nrow)
 real(wp) :: mat(nrow, nrow)
 real(wp) :: matinv(nrow, nrow)
 real(wp) :: mat_l(nrow, nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs
 type(crssparse) :: crschol

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), iaspsdf, jaspsdf, aspsdf)

 !Cholesky
 crschol= coo
 call crschol%setpermutation(perm)
 
 !SOLVE
 call crschol%getldlt()

 do i = 1, nrow
  vect = 0
  vect(i) = 1
  call crschol%solveldlt_s(matinv(:,i), vect)
 enddo

 if(verbose)call printmat(matinv)

 !SPAINV
 crs = coo

 call getija_crs(crs, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 call crs%setpermutation(perm)
 call crs%spainv()
 
 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat_l)
 if(verbose)call printmat(mat_l)

 call check(error, all(abs(pack(mat_l, mat.ne.0) - pack(matinv, mat.ne.0)) < tol_wp), 'spainv_spsd')
 if(allocated(error))return


 !Cholesky
 crschol= coo
 call crschol%setpermutation(permf(nrow))
 
 !SOLVE
 call crschol%getldlt()

 do i = 1, nrow
  vect = 0
  vect(i) = 1
  call crschol%solveldlt_s(matinv(:,i), vect)
 enddo

 if(verbose)call printmat(matinv)

 !SPAINV
 crs = coo

 call getija_crs(crs, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 call crs%setpermutation(permf(nrow))
 call crs%spainv()
 
 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat_l)
 if(verbose)call printmat(mat_l)

 call check(error, all(abs(pack(mat_l, mat.ne.0) - pack(matinv, mat.ne.0)) < tol_wp), 'spainv_spsd_permf')
 if(allocated(error))return

end subroutine

subroutine test_spainv_spsd_1(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i
 integer, parameter :: nrow = 6
 integer, parameter :: perm(nrow) = [(i, i = nrow, 1, -1)]

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: vect(nrow)
 real(wp) :: mat(nrow, nrow)
 real(wp) :: matinv(nrow, nrow)
 real(wp) :: mat_l(nrow, nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs
 type(crssparse) :: crschol

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), iaspsdf1, jaspsdf1, aspsdf1)

 !Cholesky
 crschol= coo
 call crschol%setpermutation(perm)
 
 !SOLVE
 call crschol%getldlt()

 do i = 1, nrow
  vect = 0
  vect(i) = 1
  call crschol%solveldlt_s(matinv(:,i), vect)
 enddo

 if(verbose)call printmat(matinv)

 !SPAINV
 crs = coo

 call getija_crs(crs, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 call crs%setpermutation(perm)
 call crs%spainv()
 
 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat_l)
 if(verbose)call printmat(mat_l)

 call check(error, all(abs(pack(mat_l, mat.ne.0) - pack(matinv, mat.ne.0)) < tol_wp), 'spainv_spsd1')
 if(allocated(error))return


 !Cholesky
 crschol= coo
 call crschol%setpermutation(permf(nrow))
 
 !SOLVE
 call crschol%getldlt()

 do i = 1, nrow
  vect = 0
  vect(i) = 1
  call crschol%solveldlt_s(matinv(:,i), vect)
 enddo

 if(verbose)call printmat(matinv)

 !SPAINV
 crs = coo

 call getija_crs(crs, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 call crs%setpermutation(permf(nrow))
 call crs%spainv()
 
 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat_l)
 if(verbose)call printmat(mat_l)

 call check(error, all(abs(pack(mat_l, mat.ne.0) - pack(matinv, mat.ne.0)) < tol_wp), 'spainv_spsd1_permf')
 if(allocated(error))return

end subroutine

subroutine test_spainv_spsd_2(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i
 integer, parameter :: nrow = 6
 integer, parameter :: perm(nrow) = [(i, i = nrow, 1, -1)]

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: vect(nrow)
 real(wp) :: ident(nrow,nrow)
 real(wp) :: mat(nrow, nrow)
 real(wp) :: matinv(nrow, nrow)
 real(wp) :: mat_l(nrow, nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs
 type(crssparse) :: crschol

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), iaspsdf, jaspsdf, aspsdf)

 !Cholesky
 crschol= coo
 call crschol%setpermutation(perm)
 
 !SOLVE
 matinv=0
 do i = 1, nrow
  vect = 0
  vect(i) = 1
  call crschol%solveldlt(matinv(:,i), vect)
 enddo

 if(verbose)call printmat(matinv)

 !SPAINV
 crs = coo

 call getija_crs(crs, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 call crs%setpermutation(perm)
 call crs%spainv()
 
 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat_l)
 if(verbose)call printmat(mat_l)

 call check(error, all(abs(pack(mat_l, mat.ne.0) - pack(matinv, mat.ne.0)) < tol_wp), 'spainv_spsd_2')
 if(allocated(error))return


 !Cholesky
 crschol= coo
 call crschol%setpermutation(permf(nrow))
 
 !SOLVE
 ident = 0
 do i = 1, nrow
  ident(i,i) = 1
 enddo
 matinv=0
 call crschol%solveldlt(matinv, ident)

 if(verbose)call printmat(matinv)

 !SPAINV
 crs = coo

 call getija_crs(crs, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 call crs%setpermutation(permf(nrow))
 call crs%spainv()
 
 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat_l)
 if(verbose)call printmat(mat_l)

 call check(error, all(abs(pack(mat_l, mat.ne.0) - pack(matinv, mat.ne.0)) < tol_wp), 'spainv_spsd_2_permf')
 if(allocated(error))return

end subroutine

subroutine test_spainv_failed(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i
 integer, parameter :: nrow = 6
 integer, parameter :: perm(nrow) = [(i, i = nrow, 1, -1)]

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: vect(nrow)
 real(wp) :: mat(nrow, nrow)
 real(wp) :: matinv(nrow, nrow)
 real(wp) :: mat_l(nrow, nrow)
 type(coosparse) :: coo
 type(crssparse) :: crs
 type(crssparse) :: crschol

 coo = coosparse(nrow, lupper = .true.,  unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), [ia, 2], [ja, 2], [aspsd, 2._wp])

 !Cholesky
 crschol= coo
 call crschol%setpermutation(perm)
 call crschol%chol()
 
 do i = 1, nrow
  vect = 0
  vect(i) = 1
  call crschol%isolve(matinv(:,i), vect) !matinv reference
 enddo

 !SPAINV
 crs = coo

 call getija_crs(crs, iat, jat, at, mat)
 if(verbose)call printmat(mat)

 call crs%setpermutation(perm)
 call crs%spainv()
 
 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat_l)

 !this test should failed because matinv contains more non-zero elements than mat_l
 call check(error, all(abs(mat_l - matinv ) < tol_wp), 'should_spainv_failed')

end subroutine

!SUBMATRIX
subroutine test_submatrix_off_full_full(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 6
 integer, parameter :: ncol = 5
 integer, parameter :: rowstart = 1, rowend = 3, colstart = 4, colend = 5
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, ncol)
 type(coosparse) :: coo
 type(crssparse) :: crs
 type(crssparse) :: crssub

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a, iat, jat, at, mat)

 crs = coo

 crssub = crs%submatrix(rowstart, rowend, colstart, colend)

 call check(error, crssub%getdim(1), rowend - rowstart + 1, 'submatrix_off_dim1')
 if(allocated(error))return

 call check(error, crssub%getdim(2), colend - colstart + 1, 'submatrix_off_dim2')
 if(allocated(error))return

 call check(error, all(getmat(crssub) == mat(rowstart:rowend, colstart:colend)) &
            , 'submatrix_off_full_full')

end subroutine

subroutine test_submatrix_off_full_upper(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i, j

 integer, parameter :: nrow = 6
 integer, parameter :: ncol = 5
 integer, parameter :: rowstart = 1, rowend = 3, colstart = 4, colend = 5
 integer, parameter :: subrow = rowend - rowstart + 1
 integer, parameter :: subcol = colend - colstart + 1
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol
 logical, parameter :: lsubval(subrow, subcol) = reshape(merge(.true., .false.&
                        , [((i, j=1, subrow), i = 1, subcol)] .ge.&
                          [((j, j=1, subrow), i = 1, subcol)])&
                        , [subrow, subcol])

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, ncol)
 type(coosparse) :: coo
 type(crssparse) :: crs
 type(crssparse) :: crssub

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a, iat, jat, at, mat)

 crs = coo

 crssub = crs%submatrix(rowstart, rowend, colstart, colend, lupper = .true.)

 call check(error, crssub%getdim(1), rowend - rowstart + 1, 'submatrix_off_dim1')
 if(allocated(error))return

 call check(error, crssub%getdim(2), colend - colstart + 1, 'submatrix_off_dim2')
 if(allocated(error))return

 call check(error, all(pack(getmat(crssub), lsubval) == pack(&
            mat(rowstart:rowend, colstart:colend), lsubval)) &
            , 'submatrix_off_full_upper')
 if(allocated(error))return

 call check(error, all(pack(getmat(crssub), .not.lsubval) ==  0._wp) &
            , 'submatrix_off_full_upper_0')

end subroutine

subroutine test_submatrix_off_upper_full(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 6
 integer, parameter :: ncol = 5
 integer, parameter :: rowstart = 1, rowend = 3, colstart = 4, colend = 5
 integer, parameter :: subrow = rowend - rowstart + 1
 integer, parameter :: subcol = colend - colstart + 1
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia.le.ja

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, ncol)
 type(coosparse) :: coo
 type(crssparse) :: crs
 type(crssparse) :: crssub

 coo = coosparse(nrow, n = ncol, lupper = .true.,  unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a, iat, jat, at, mat)

 crs = coo
 crssub = crs%submatrix(rowstart, rowend, colstart, colend, lupper = .false.)

 call check(error, crssub%getdim(1), rowend - rowstart + 1, 'submatrix_off_dim1')
 if(allocated(error))return

 call check(error, crssub%getdim(2), colend - colstart + 1, 'submatrix_off_dim2')
 if(allocated(error))return

 call check(error, all(getmat(crssub) == &
            mat(rowstart:rowend, colstart:colend)) &
            , 'submatrix_off_upper_full')
 if(allocated(error))return

end subroutine

subroutine test_submatrix_full_full(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 6
 integer, parameter :: ncol = 5
 integer, parameter :: rowstart = 1, rowend = 4, colstart = 3, colend = 5
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, ncol)
 type(coosparse) :: coo
 type(crssparse) :: crs
 type(crssparse) :: crssub

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a, iat, jat, at, mat)

 crs = coo

 crssub = crs%submatrix(rowstart, rowend, colstart, colend)

 call check(error, crssub%getdim(1), rowend - rowstart + 1, 'submatrix_dim1')
 if(allocated(error))return

 call check(error, crssub%getdim(2), colend - colstart + 1, 'submatrix_dim2')
 if(allocated(error))return

 call check(error, all(getmat(crssub) == mat(rowstart:rowend, colstart:colend)) &
            , 'submatrix_full_full')

end subroutine

subroutine test_submatrix_full_upper(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i, j

 integer, parameter :: nrow = 6
 integer, parameter :: ncol = 5
 integer, parameter :: rowstart = 1, rowend = 4, colstart = 3, colend = 5
 integer, parameter :: subrow = rowend - rowstart + 1
 integer, parameter :: subcol = colend - colstart + 1
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol
 logical, parameter :: lsubval(subrow, subcol) = reshape(merge(.true., .false.&
                        , [((i, j=1, subrow), i = 1, subcol)] .ge.&
                          [((j, j=1, subrow), i = 1, subcol)])&
                        , [subrow, subcol])

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, ncol)
 type(coosparse) :: coo
 type(crssparse) :: crs
 type(crssparse) :: crssub

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a, iat, jat, at, mat)

 crs = coo

 crssub = crs%submatrix(rowstart, rowend, colstart, colend, lupper = .true.)

 call check(error, crssub%getdim(1), rowend - rowstart + 1, 'submatrix_dim1')
 if(allocated(error))return

 call check(error, crssub%getdim(2), colend - colstart + 1, 'submatrix_dim2')
 if(allocated(error))return

 call check(error, all(pack(getmat(crssub), lsubval) == pack(&
            mat(rowstart:rowend, colstart:colend), lsubval)) &
            , 'submatrix_full_upper')
 if(allocated(error))return

 call check(error, all(pack(getmat(crssub), .not.lsubval) ==  0._wp) &
            , 'submatrix_full_upper_0')

end subroutine

subroutine test_submatrix_upper_full(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 6
 integer, parameter :: ncol = 5
 integer, parameter :: rowstart = 1, rowend = 4, colstart = 3, colend = 5
 integer, parameter :: subrow = rowend - rowstart + 1
 integer, parameter :: subcol = colend - colstart + 1
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia.le.ja

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat(nrow, ncol)
 type(coosparse) :: coo
 type(crssparse) :: crs
 type(crssparse) :: crssub

 coo = coosparse(nrow, n = ncol, lupper = .true.,  unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a, iat, jat, at, mat)

 crs = coo

 crssub = crs%submatrix(rowstart, rowend, colstart, colend, lupper = .false.)

 call check(error, crssub%getdim(1), rowend - rowstart + 1, 'submatrix_dim1')
 if(allocated(error))return

 call check(error, crssub%getdim(2), colend - colstart + 1, 'submatrix_dim2')
 if(allocated(error))return

 call check(error, all(getmat(crssub) == mat(rowstart:rowend, colstart:colend)) &
            , 'submatrix_upper_full')

end subroutine

subroutine test_submatrix_dense_upper_upper(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 6
 integer, parameter :: ncol = nrow
 integer, parameter :: vector(*) = [1, 3, 4]
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia.le.ja

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat1(nrow, ncol)
 real(wp), allocatable :: mat(:, :)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, n = ncol, lupper = .true.,  unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a, iat, jat, at, mat)

 crs = coo

 call crs%submatrix_dense(vector, mat)!, lupper = .true.)

 mat1 = getmat(crs)

 if(verbose)call printmat(mat1(vector, vector))
 if(verbose)call printmat(mat)

 call check(error, all(mat1(vector, vector) == mat) &
            , 'submatrix_dense_upper_upper')

end subroutine

subroutine test_submatrix_dense_upper_upper_sym(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 6
 integer, parameter :: ncol = nrow
 integer, parameter :: vector(*) = [1, 3, 4]
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia.le.ja

 integer :: i
 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat1(nrow, ncol)
 real(wp), allocatable :: mat(:, :)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, n = ncol, lupper = .true.,  unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a, iat, jat, at, mat)
 call coo%setsymmetric(.true.)

 crs = coo

 call crs%submatrix_dense(vector, mat, lupper = .true.)

 mat1 = getmat(crs)
 do i = 2, size(mat1, 1)
  mat1(i, 1:i-1) = 0
 enddo

 if(verbose) call printmat(mat1(vector, vector))
 if(verbose)call printmat(mat)

 call check(error, all(mat1(vector, vector) == mat) &
            , 'submatrix_dense_upper_upper_sym')

end subroutine

subroutine test_submatrix_dense_upper_full(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 6
 integer, parameter :: ncol = nrow
 integer, parameter :: vector(*) = [1, 3, 4]
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia.le.ja

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: mat1(nrow, ncol)
 real(wp), allocatable :: mat(:, :)
 type(coosparse) :: coo
 type(crssparse) :: crs

 coo = coosparse(nrow, n = ncol, lupper = .true.,  unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a, iat, jat, at, mat)
 call coo%setsymmetric(.true.)

 crs = coo

 call crs%submatrix_dense(vector, mat, lupper = .false.)

 mat1 = getmat(crs)

 if(verbose)call printmat(mat1(vector, vector))
 if(verbose)call printmat(mat)

 call check(error, all(mat1(vector, vector) == mat) &
            , 'submatrix_dense_upper_full')

end subroutine

!INTERNAL
subroutine getija_crs(crs, iat, jat, at, mat)
 type(crssparse), intent(inout) :: crs
 integer, allocatable , intent(out), optional :: iat(:), jat(:) 
 real(wp), allocatable , intent(out), optional :: at(:)
 real(wp), intent(out), optional :: mat(:,:)

 integer :: i, j
 integer :: nrow, ncol
 real(wp) :: val

 nrow = crs%getdim(1)
 ncol = crs%getdim(2)

 if(present(mat))then
  if(nrow.ne.size(mat,1) .or. ncol.ne.size(mat,2))return
 endif

 allocate(iat(0), jat(0), at(0))
 do i = 1, nrow
  do j = 1, ncol
   val = crs%get(i, j)
   if(present(mat))then
    if(i.le.nrow.and.j.le.ncol)mat(i, j) = val
   endif
   if (val .ne. 0._wp) then
    iat = [iat, i]
    jat = [jat, j]
    at = [at, val]
   endif
  enddo
 enddo

end subroutine

pure function permf(nrow) result(perm)
 integer, intent(in) :: nrow
 integer :: perm(nrow)

 integer :: i

 perm = [(i, i = nrow, nrow/2 + 1, -1), (i, i = 1, nrow/2)]

end function

end module

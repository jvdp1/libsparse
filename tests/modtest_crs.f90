module modtest_crs
#if (_DP==0)
 use, intrinsic :: iso_fortran_env, only: int64, wp => real32
#else
 use, intrinsic :: iso_fortran_env, only: int64, wp => real64
#endif
 use testdrive, only: new_unittest, unittest_type, error_type, check
 use modsparse, only: coosparse, crssparse, assignment(=)
 use modtest_common, only: tol_wp, verbose, ia, ja, a, aspsd, addval => addval_coo&
                       , getmat, matcheck, printmat
 implicit none
 private

 public :: collect_crs

 integer, parameter :: sparse_unit = 556
 
contains

!> Collect all exported unit tests
subroutine collect_crs(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  !how to test: getordering, print, printsquare, save, sort, submatrix

  testsuite = [ &
    new_unittest("crs constructor", test_constructor) &
    , new_unittest("crs add", test_add) &
    , new_unittest("crs add nel", test_add_nel) &
    , new_unittest("crs add lupper", test_add_lupper) &
    , new_unittest("crs ncol add", test_ncol_add) &
    , new_unittest("crs ncol add nel", test_ncol_add_nel) &
    , new_unittest("crs ncol add lupper", test_ncol_add_lupper) &
#if (_SPAINV==1)
    , new_unittest("crs chol", test_chol) &
#endif
    , new_unittest("crs diag vect", test_diag_vect) &
    , new_unittest("crs diag vect lupper", test_diag_vect_lupper) &
    , new_unittest("crs ncol diag vect", test_ncol_diag_vect) &
    , new_unittest("crs ncol diag vect lupper", test_ncol_diag_vect_lupper) &
    , new_unittest("crs get", test_get) &
    , new_unittest("crs get nel", test_get_nel) &
    , new_unittest("crs get lupper", test_get_lupper) &
    , new_unittest("crs get lupper_sym", test_get_lupper_sym) &
    , new_unittest("crs ncol get", test_ncol_get) &
    , new_unittest("crs ncol get nel", test_ncol_get_nel) &
    , new_unittest("crs ncol get lupper", test_ncol_get_lupper) &
#if (_SPAINV==1)
    , new_unittest("crs ichol", test_ichol) &
    , new_unittest("crs ichol_failed", test_ichol_failed, should_fail = .true.) &
    , new_unittest("crs isolve", test_isolve) &
    , new_unittest("crs ldlt", test_ldlt) &
#endif
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
    , new_unittest("crs solveldlt", test_solve_vector) &
    , new_unittest("crs solveldlt", test_solve_vector_perm) &
    , new_unittest("crs solveldlt_arr", test_solve_array) &
    , new_unittest("crs solveldlt_arr", test_solve_array_perm) &
#else
    , new_unittest("crs solveldlt", test_solve_vector, should_fail = .true.) &
    , new_unittest("crs solveldlt", test_solve_vector_perm, should_fail = .true.) &
    , new_unittest("crs solveldlt_arr", test_solve_array, should_fail = .true.) &
    , new_unittest("crs solveldlt_arr", test_solve_array_perm, should_fail = .true.) &
#endif
#if (_SPAINV==1)
    , new_unittest("crs solveldlt", test_solveldlt) &
    , new_unittest("crs spainv", test_spainv) &
    , new_unittest("crs spainv_failed", test_spainv_failed, should_fail = .true.) &
#endif
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

#if (_SPAINV==1)
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

 !Complete Cholesky decomposition
 call crs%setpermutation(perm)
 call crs%chol()

 deallocate(iat); deallocate(jat); deallocate(at)
 call getija_crs(crs, iat, jat, at, matchol)
 
 call check(error, all(abs(mat(perm, perm) - & 
                       matmul(transpose(matchol), matchol)) < tol_wp)&
                 , 'Chol')
 if(allocated(error))return

 !Complete Cholesky decomposition
 crs = coo
 call crs%setpermutation(permf(nrow))
 call crs%chol()

 deallocate(iat); deallocate(jat); deallocate(at)
 call getija_crs(crs, iat, jat, at, matchol)
 
 call check(error, all(abs(mat(permf(nrow), permf(nrow)) - & 
                       matmul(transpose(matchol), matchol)) < tol_wp)&
                 , 'Chol_permf')

end subroutine
#endif

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

#if (_SPAINV==1)
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
 call crschol%setpermutation(perm)
 call crschol%chol()
 
 call getija_crs(crschol, iat, jat, at, matchol)

 crs = coo

 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat)

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
 call crschol%setpermutation(permf(nrow))
 call crschol%chol()
 
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

 !Cholesky
 call crs%setpermutation(perm)
 call crs%chol()
 
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

 call crs%setpermutation(permf(nrow))
 call crs%chol()
 
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
#endif

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
   coo = coosparse(nrow, n = ncol, lupper = .true., unlog = sparse_unit)
  else
   coo = coosparse(nrow, n = ncol, unlog = sparse_unit)
  endif
 else
  if(upper == 'y')then
   coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)
  else
   coo = coosparse(nrow, unlog = sparse_unit)
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
   coo = coosparse(nrow, n = ncol, lupper = .true., unlog = sparse_unit)
  else
   coo = coosparse(nrow, n = ncol, unlog = sparse_unit)
  endif
 else
  if(upper == 'y')then
   coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)
  else
   coo = coosparse(nrow, unlog = sparse_unit)
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

#if (_SPAINV==1)
!SOLVE LDLt
subroutine test_solveldlt(error)
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

 !LDLt
 call crs%setpermutation(perm)
 call crs%getldlt()
 
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
 call crs%getldlt()
 
 do i = 1, nrow
  call crs%solveldlt(mat_d(:,i), mat(:,i))
 enddo

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

 call crs%setpermutation(perm)
 call crs%spainv()
 
 deallocate(iat, jat, at)
 call getija_crs(crs, iat, jat, at, mat_l)

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

subroutine test_spainv_failed(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i
 integer, parameter :: nrow = 6
 integer, parameter :: perm(nrow) = [(i, i = nrow, 1, -1)]

 integer, allocatable :: iat(:), jat(:)
 real(wp), allocatable :: at(:)
 real(wp) :: vect(nrow)
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

 call crs%setpermutation(perm)
 call crs%spainv()
 
 call getija_crs(crs, iat, jat, at, mat_l)

 !this test should failed because matinv contains more non-zero elements than mat_l
 call check(error, all(abs(mat_l - matinv ) < tol_wp), 'should_spainv_failed')

end subroutine
#endif


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

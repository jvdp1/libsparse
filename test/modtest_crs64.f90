module modtest_crs64
#if (_DP==0)
 use, intrinsic :: iso_fortran_env, only: int64, wp => real32
#else
 use, intrinsic :: iso_fortran_env, only: int64, wp => real64
#endif
 use testdrive, only: new_unittest, unittest_type, error_type, check
 use modsparse, only: coosparse, crssparse => crssparse64, assignment(=)
 use modtest_common, only: tol_wp, verbose, ia, ja, a, aspsd, addval => addval_coo&
                       , getmat, matcheck, printmat&
                       , iaspsdf, jaspsdf, aspsdf&
                       , iaspsdf1, jaspsdf1, aspsdf1
 implicit none
 private

 public :: collect_crs64

 integer, parameter :: sparse_unit = 556
 
contains

!> Collect all exported unit tests
subroutine collect_crs64(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  !how to test: getordering, print, printsquare, save, sort

  testsuite = [ &
    new_unittest("crs64 constructor", test_constructor) &
    , new_unittest("crs64 add", test_add) &
    , new_unittest("crs64 add nel", test_add_nel) &
    , new_unittest("crs64 add lupper", test_add_lupper) &
    , new_unittest("crs64 ncol add", test_ncol_add) &
    , new_unittest("crs64 ncol add nel", test_ncol_add_nel) &
    , new_unittest("crs64 ncol add lupper", test_ncol_add_lupper) &
    , new_unittest("crs64 cg", test_cg, should_fail = .true.) &
    , new_unittest("crs64 get", test_get) &
    , new_unittest("crs64 get nel", test_get_nel) &
    , new_unittest("crs64 get lupper", test_get_lupper) &
    , new_unittest("crs64 get lupper_sym", test_get_lupper_sym) &
    , new_unittest("crs64 ncol get", test_ncol_get) &
    , new_unittest("crs64 ncol get nel", test_ncol_get_nel) &
    , new_unittest("crs64 ncol get lupper", test_ncol_get_lupper) &
    , new_unittest("crs64 nonzero", test_nonzero) &
    , new_unittest("crs64 nonzero_sym", test_nonzero_sym) &
    , new_unittest("crs64 scale", test_scale) &
#if (_PARDISO==1)
    , new_unittest("crs64 solveldlt", test_solve_vector) &
    , new_unittest("crs64 solveldlt", test_solve_vector_perm) &
    , new_unittest("crs64 solveldlt_arr", test_solve_array) &
    , new_unittest("crs64 solveldlt_arr", test_solve_array_perm) &
#else
    , new_unittest("crs64 solveldlt", test_solve_vector, should_fail = .true.) &
    , new_unittest("crs64 solveldlt", test_solve_vector_perm, should_fail = .true.) &
    , new_unittest("crs64 solveldlt_arr", test_solve_array, should_fail = .true.) &
    , new_unittest("crs64 solveldlt_arr", test_solve_array_perm, should_fail = .true.) &
#endif
    ]

  !to check: diag_mat

end subroutine collect_crs64

!CONSTRUCTOR + GETDIM + ISSQUARE
subroutine test_constructor(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 integer(int64), parameter :: nel = 9
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

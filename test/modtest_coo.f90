module modtest_coo
#if (_DP==0)
 use, intrinsic :: iso_fortran_env, only: int64, output_unit, wp => real32
#else
 use, intrinsic :: iso_fortran_env, only: int64, output_unit, wp => real64
#endif
 use testdrive, only: new_unittest, unittest_type, error_type, check
  use modsparse, only: coosparse
  use modrandom, only: setseed, rand_stdnormal
  use modtest_common, only: tol_wp, verbose, ia, ja, a, aspsd&
                        , addval => addval_coo, getmat, matcheck, printmat
 implicit none
 private

 public :: collect_coo

 integer, parameter :: sparse_unit = 555
 
contains

!> Collect all exported unit tests
subroutine collect_coo(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("coo constructor", test_constructor) &
    , new_unittest("coo add", test_add) &
    , new_unittest("coo add nel", test_add_nel) &
    , new_unittest("coo add lupper", test_add_lupper) &
    , new_unittest("coo ncol add", test_ncol_add) &
    , new_unittest("coo ncol add nel", test_ncol_add_nel) &
    , new_unittest("coo ncol add lupper", test_ncol_add_lupper) &
    , new_unittest("coo cg", test_cg) &
    , new_unittest("coo diag vect", test_diag_vect) &
    , new_unittest("coo diag vect lupper", test_diag_vect_lupper) &
    , new_unittest("coo ncol diag vect", test_ncol_diag_vect) &
    , new_unittest("coo ncol diag vect lupper", test_ncol_diag_vect_lupper) &
    , new_unittest("coo get", test_get) &
    , new_unittest("coo get nel", test_get_nel) &
    , new_unittest("coo get lupper", test_get_lupper) &
    , new_unittest("coo get lupper_sym", test_get_lupper_sym) &
    , new_unittest("coo ncol get", test_ncol_get) &
    , new_unittest("coo ncol get nel", test_ncol_get_nel) &
    , new_unittest("coo ncol get lupper", test_ncol_get_lupper) &
    , new_unittest("coo multbyv_n_n_n", test_multbyv_n_n_n) &
    , new_unittest("coo multbyv_n_n_y", test_multbyv_n_n_y) &
    , new_unittest("coo multbyv_n_y_n", test_multbyv_n_y_n) &
    , new_unittest("coo multbyv_n_y_y", test_multbyv_n_y_y) &
    , new_unittest("coo multbyv_y_n_n", test_multbyv_y_n_n) &
    , new_unittest("coo multbyv_y_n_y", test_multbyv_y_n_y) &
    , new_unittest("coo multbyv_y_y_n", test_multbyv_y_y_n) &
    , new_unittest("coo multbyv_y_y_y", test_multbyv_y_y_y) &
    , new_unittest("coo multbyv_sym", test_multbyv_sym) &
    , new_unittest("coo multbymat_n_n_n", test_multbymat_n_n_n) &
    , new_unittest("coo multbymat_n_n_y", test_multbymat_n_n_y) &
    , new_unittest("coo multbymat_n_y_n", test_multbymat_n_y_n) &
    , new_unittest("coo multbymat_n_y_y", test_multbymat_n_y_y) &
    , new_unittest("coo multbymat_y_n_n", test_multbymat_y_n_n) &
    , new_unittest("coo multbymat_y_n_y", test_multbymat_y_n_y) &
    , new_unittest("coo multbymat_y_y_n", test_multbymat_y_y_n) &
    , new_unittest("coo multbymat_y_y_y", test_multbymat_y_y_y) &
    , new_unittest("coo multbymat_sym", test_multbymat_sym) &
    , new_unittest("coo xAx", test_xAx_coo) &
   , new_unittest("coo xAy", test_xAy_coo) &
    , new_unittest("coo traceproduct", test_traceproduct_coo) &
    , new_unittest("coo nonzero", test_nonzero) &
    , new_unittest("coo nonzero_sym", test_nonzero_sym) &
    , new_unittest("coo scale", test_scale) &
    , new_unittest("coo submatrix_off_full_full", test_submatrix_off_full_full) &
    , new_unittest("coo submatrix_off_full_upper", test_submatrix_off_full_upper) &
    , new_unittest("coo submatrix_off_upper_full", test_submatrix_off_upper_full) &
    , new_unittest("coo submatrix_full_full", test_submatrix_full_full) &
    , new_unittest("coo submatrix_full_upper", test_submatrix_full_upper) &
    , new_unittest("coo submatrix_upper_full", test_submatrix_upper_full) &
    ]

end subroutine collect_coo

!CONSTRUCTOR + GETDIM + ISSQUARE
subroutine test_constructor(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 type(coosparse)::coo

 coo = coosparse(nrow)
 call check(error, coo%getdim(1), nrow, more = 'getdim(1)')
 call check(error, coo%getdim(2), nrow, more = 'getdim(2)')
 call check(error, coo%getdim(3) < 0, more = 'getdim(3)')
 call check(error, coo%issquare(), nrow.eq.nrow, 'isquare')
 if (allocated(error)) return

 coo = coosparse(nrow, nel = 4_int64)
 call check(error, coo%getdim(1), nrow, more = 'getdim(1)')
 call check(error, coo%getdim(2), nrow, more = 'getdim(2)')
 call check(error, coo%issquare(), nrow.eq.nrow, 'isquare')
 if (allocated(error)) return

 coo = coosparse(nrow, lupper = .true.)
 call check(error, coo%getdim(1), nrow, more = 'getdim(1)')
 call check(error, coo%getdim(2), nrow, more = 'getdim(2)')
 call check(error, coo%issquare(), nrow.eq.nrow, 'isquare')
 if (allocated(error)) return

 coo = coosparse(nrow, n = ncol)
 call check(error, coo%getdim(1), nrow, more = 'getdim(1)')
 call check(error, coo%getdim(2), ncol, more = 'getdim(2)')
 call check(error, coo%issquare(), nrow.eq.ncol, 'isquare')
 if (allocated(error)) return

 coo = coosparse(nrow, n = ncol, nel = 4_int64)
 call check(error, coo%getdim(1), nrow, more = 'getdim(1)')
 call check(error, coo%getdim(2), ncol, more = 'getdim(2)')
 call check(error, coo%issquare(), nrow.eq.ncol, 'isquare')
 if (allocated(error)) return

 coo = coosparse(nrow, n = ncol, nel = 4_int64, lupper = .true.)
 call check(error, coo%getdim(1), nrow, more = 'getdim(1)')
 call check(error, coo%getdim(2), ncol, more = 'getdim(2)')
 call check(error, coo%issquare(), nrow.eq.ncol, 'isquare')
 if (allocated(error)) return

 coo = coosparse(nrow, n = ncol, nel = 4_int64, lupper = .true., unlog = output_unit)
 call check(error, coo%getdim(1), nrow, more = 'getdim(1)')
 call check(error, coo%getdim(2), ncol, more = 'getdim(2)')
 call check(error, coo%issquare(), nrow.eq.ncol, 'isquare')
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

 coo = coosparse(nrow, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a, iat, jat, at, mat)

 call check(error, size(iat), size(pack(ia, lvalid)), more = 'size(iat)')
 call check(error, size(jat), size(pack(ja, lvalid)), more = 'size(jat)')
 call check(error, size(at), size(pack(a, lvalid)), more = 'size(at)')
 if (allocated(error)) return

 call check(error, all(iat == pack(ia, lvalid)), 'coo_ia')
 call check(error, all(jat == pack(ja, lvalid)), 'coo_ja')
 call check(error, all(at == pack(a, lvalid)), 'coo_a')
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

 coo = coosparse(nrow, nel = 4_int64, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a, iat, jat, at, mat)

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

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a, iat, jat, at, mat)

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

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a, iat, jat, at, mat)

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

 coo = coosparse(nrow, n = ncol, nel = 2_int64,  unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a, iat, jat, at, mat)

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

 coo = coosparse(nrow, n = ncol, lupper = .true.,  unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a, iat, jat, at, mat)


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

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, aspsd)
 call coo%add(2, 2, 202._wp) !to set the matrix SPD

 xmat = 0
 do i = 1, ncol
  tol = 1.e-6
  maxiter = nrow
  y = 0
  y(i) = 1
  call coo%cg(xmat(:,i), y, maxiter, tol)
 enddo

 !mat_l: expected result 
 mat_l = reshape([(merge(1._wp, 0._wp, i/nrow.eq.mod(i,nrow)), i = 0, nrow**2 - 1)], [nrow, nrow])
 where(getmat(coo) <= tol_wp) mat_l = 0._wp

 call check(error, all(abs(matmul(getmat(coo), xmat) - mat_l) < tol_wp), 'coo cg')

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

 coo = coosparse(nrow, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a)

 diag = coo%diag()

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

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a)

 diag = coo%diag()

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

 coo = coosparse(nrow, n = ncol,  unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 diag = coo%diag()

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

 coo = coosparse(nrow, n = ncol,  unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 diag = coo%diag()

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

 coo = coosparse(nrow, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a)

 call check(error, all(getmat(coo) == matcheck(nrow, ncol, ia, ja, a, lvalid)), 'get')

end subroutine

subroutine test_get_nel(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol

 type(coosparse) :: coo

 coo = coosparse(nrow, nel = 4_int64, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 call check(error, all(getmat(coo) == matcheck(nrow, ncol, ia, ja, a, lvalid)), 'get')

end subroutine

subroutine test_get_lupper(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia .le. ja

 type(coosparse) :: coo

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 call check(error, all(getmat(coo) == matcheck(nrow, ncol, ia, ja, a, lvalid)), 'get')

end subroutine

subroutine test_get_lupper_sym(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia .le. ja

 type(coosparse) :: coo

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 call check(error, all(getmat(coo) == &
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

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 call check(error, all(getmat(coo) == matcheck(nrow, ncol, ia, ja, a, lvalid)), 'get')

end subroutine

subroutine test_ncol_get_nel(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol

 type(coosparse) :: coo

 coo = coosparse(nrow, n = ncol, nel = 2_int64,  unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 call check(error, all(getmat(coo) == matcheck(nrow, ncol, ia, ja, a, lvalid)), 'get')

end subroutine

subroutine test_ncol_get_lupper(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 logical, parameter :: lvalid(size(ia)) =  ia.le.nrow .and. ja.le.ncol .and. ia.le.ja

 type(coosparse) :: coo

 coo = coosparse(nrow, n = ncol, lupper = .true.,  unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 call check(error, all(getmat(coo) == matcheck(nrow, ncol, ia, ja, a, lvalid)), 'get')

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

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 y = 1._wp
 call coo%mult(alpha, 'n', x, val, y)

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

 call coo%mult(alpha, merge('t', 'n', trans == 'y'), x, val, y)

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

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)
 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 y = 1._wp
 call coo%mult(alpha, 'n', x, val, y)

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

 call coo%mult(alpha, merge('t', 'n', trans == 'y'), x, val, y)

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
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol

 integer(int64), parameter :: p(1) = count(lvalid .and. a.ne.0._wp)
 type(coosparse) :: coo

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 call check(error, coo%nonzero(), p(1), 'nonzero')

end subroutine

subroutine test_nonzero_sym(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia.le.ja

 integer(int64), parameter :: p(1) = count(lvalid .and. a.ne.0._wp)
 type(coosparse) :: coo

 coo = coosparse(nrow, n = ncol, lupper = .true., unlog = sparse_unit)

 call coo%setsymmetric()

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a)

 call check(error, coo%nonzero(), p(1), 'nonzero')

end subroutine

!SCALE
subroutine test_scale(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol

 real(wp), parameter :: scalefact = 1000._wp
 type(coosparse) :: coo

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a)

 call coo%scale(scalefact)

 call check(error, all(getmat(coo) == scalefact * matcheck(nrow, ncol, ia, ja, a, lvalid)), 'scale')

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
 type(coosparse) :: coosub

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a, iat, jat, at, mat)


 coosub = coo%submatrix(rowstart, rowend, colstart, colend)

 call check(error, coosub%getdim(1), rowend - rowstart + 1, 'submatrix_off_dim1')
 if(allocated(error))return

 call check(error, coosub%getdim(2), colend - colstart + 1, 'submatrix_off_dim2')
 if(allocated(error))return

 call check(error, all(getmat(coosub) == mat(rowstart:rowend, colstart:colend)) &
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
 type(coosparse) :: coosub

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a, iat, jat, at, mat)

 coosub = coo%submatrix(rowstart, rowend, colstart, colend, lupper = .true.)

 call check(error, coosub%getdim(1), rowend - rowstart + 1, 'submatrix_off_dim1')
 if(allocated(error))return

 call check(error, coosub%getdim(2), colend - colstart + 1, 'submatrix_off_dim2')
 if(allocated(error))return

 call check(error, all(pack(getmat(coosub), lsubval) == pack(&
            mat(rowstart:rowend, colstart:colend), lsubval)) &
            , 'submatrix_off_full_upper')
 if(allocated(error))return

 call check(error, all(pack(getmat(coosub), .not.lsubval) ==  0._wp) &
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
 type(coosparse) :: coosub

 coo = coosparse(nrow, n = ncol, lupper = .true.,  unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a, iat, jat, at, mat)

 coosub = coo%submatrix(rowstart, rowend, colstart, colend, lupper = .false.)

 call check(error, coosub%getdim(1), rowend - rowstart + 1, 'submatrix_off_dim1')
 if(allocated(error))return

 call check(error, coosub%getdim(2), colend - colstart + 1, 'submatrix_off_dim2')
 if(allocated(error))return

 call check(error, all(getmat(coosub) == &
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
 type(coosparse) :: coosub

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a, iat, jat, at, mat)


 coosub = coo%submatrix(rowstart, rowend, colstart, colend)

 call check(error, coosub%getdim(1), rowend - rowstart + 1, 'submatrix_dim1')
 if(allocated(error))return

 call check(error, coosub%getdim(2), colend - colstart + 1, 'submatrix_dim2')
 if(allocated(error))return

 call check(error, all(getmat(coosub) == mat(rowstart:rowend, colstart:colend)) &
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
 type(coosparse) :: coosub

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a, iat, jat, at, mat)

 coosub = coo%submatrix(rowstart, rowend, colstart, colend, lupper = .true.)

 call check(error, coosub%getdim(1), rowend - rowstart + 1, 'submatrix_dim1')
 if(allocated(error))return

 call check(error, coosub%getdim(2), colend - colstart + 1, 'submatrix_dim2')
 if(allocated(error))return

 call check(error, all(pack(getmat(coosub), lsubval) == pack(&
            mat(rowstart:rowend, colstart:colend), lsubval)) &
            , 'submatrix_full_upper')
 if(allocated(error))return

 call check(error, all(pack(getmat(coosub), .not.lsubval) ==  0._wp) &
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
 type(coosparse) :: coosub

 coo = coosparse(nrow, n = ncol, lupper = .true.,  unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, a, iat, jat, at, mat)

 coosub = coo%submatrix(rowstart, rowend, colstart, colend, lupper = .false.)

 call check(error, coosub%getdim(1), rowend - rowstart + 1, 'submatrix_dim1')
 if(allocated(error))return

 call check(error, coosub%getdim(2), colend - colstart + 1, 'submatrix_dim2')
 if(allocated(error))return

  call check(error, all(getmat(coosub) == mat(rowstart:rowend, colstart:colend)) &
            , 'submatrix_upper_full')

end subroutine

!QUADRATIC FORM
subroutine test_xAx_coo(error)
  type(error_type), allocatable, intent(out) :: error

  integer :: i
  integer, parameter :: n = 6

  type(coosparse) :: coo
  real(wp) :: mat(n,n), x(n), expected, quad
  real(wp) :: tol

  !Known SPD matrix (aspsd + a(2,2)=202 => SPD)
  coo = coosparse(n, lupper = .true., unlog = sparse_unit)
  call coo%setsymmetric()
  call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, aspsd)
  call coo%add(2, 2, 202._wp)

  mat = getmat(coo)

  !fixed x
  x = [1._wp, 2._wp, 3._wp, 4._wp, 5._wp, 6._wp]
  quad = coo%xAx(x)
  expected = dot_product(x, matmul(mat, x))
  tol = 1.e-3_wp * max(abs(expected), 1._wp)
  call check(error, abs(quad - expected) < tol, more = 'coo xAx fixed')
  if(allocated(error))return

  !random x
  call setseed(42)
  do i = 1, n
    x(i) = rand_stdnormal()
  enddo
  quad = coo%xAx(x)
  expected = dot_product(x, matmul(mat, x))
  tol = 1.e-3_wp * max(abs(expected), 1._wp)
  call check(error, abs(quad - expected) < tol, more = 'coo xAx random')
  if(allocated(error))return

  !empty matrix -> 0
  coo = coosparse(4, unlog = sparse_unit)
  call check(error, abs(coo%xAx([1._wp, 2._wp, -1._wp, 3._wp])) < tol_wp&
           , more = 'coo xAx empty')

end subroutine test_xAx_coo

!BILINEAR FORM
subroutine test_xAy_coo(error)
  type(error_type), allocatable, intent(out) :: error

  integer :: i
  integer, parameter :: n = 6

  type(coosparse) :: coo
  real(wp) :: mat(n,n), x(n), y(n), expected, bilin
  real(wp) :: tol

  !Known SPD matrix (aspsd + a(2,2)=202 => SPD)
  coo = coosparse(n, lupper = .true., unlog = sparse_unit)
  call coo%setsymmetric()
  call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, aspsd)
  call coo%add(2, 2, 202._wp)

  mat = getmat(coo)

  !fixed x, y (x != y)
  x = [1._wp, 2._wp, 3._wp, 4._wp, 5._wp, 6._wp]
  y = [2._wp, 1._wp, 4._wp, 3._wp, 6._wp, 5._wp]
  bilin = coo%xAy(x, y)
  expected = dot_product(x, matmul(mat, y))
  tol = 1.e-3_wp * max(abs(expected), 1._wp)
  call check(error, abs(bilin - expected) < tol, more = 'coo xAy fixed')
  if(allocated(error))return

  !identity: xAx(x) == xAy(x, x)
  call check(error, abs(coo%xAx(x) - coo%xAy(x, x)) < tol_wp&
           , more = 'coo xAx identity xAy(x,x)')
  if(allocated(error))return

  !random x, y
  call setseed(42)
  do i = 1, n
    x(i) = rand_stdnormal()
    y(i) = rand_stdnormal()
  enddo
  bilin = coo%xAy(x, y)
  expected = dot_product(x, matmul(mat, y))
  tol = 1.e-3_wp * max(abs(expected), 1._wp)
  call check(error, abs(bilin - expected) < tol, more = 'coo xAy random')
  if(allocated(error))return

  !empty matrix -> 0
  coo = coosparse(4, unlog = sparse_unit)
  call check(error, abs(coo%xAy([1._wp, 2._wp, -1._wp, 3._wp], [3._wp, -1._wp, 2._wp, 1._wp])) < tol_wp&
           , more = 'coo xAy empty')

end subroutine test_xAy_coo

!TRACE OF A BLOCK PRODUCT
subroutine test_traceproduct_coo(error)
  type(error_type), allocatable, intent(out) :: error

  integer :: id
  integer, parameter :: n = 6

  type(coosparse) :: coo
  real(wp) :: mat(n,n), b(n,n), expected, result
  real(wp) :: tol

  !Known SPD matrix (aspsd + a(2,2)=202)
  coo = coosparse(n, lupper = .true., unlog = sparse_unit)
  call coo%setsymmetric()
  call addval(coo, coo%getdim(1), coo%getdim(2), ia, ja, aspsd)
  call coo%add(2, 2, 202._wp)

  mat = getmat(coo)

  !Case 1: whole-matrix block, b = identity => expected = trace(mat)
  b = 0._wp
  do id = 1, n
    b(id, id) = 1._wp
  enddo
  expected = sum([(mat(id, id), id = 1, n)])
  result = coo%traceproduct(1, n, 1, n, b)
  tol = 1.e-3_wp * max(abs(expected), 1._wp)
  call check(error, abs(result - expected) < tol, more = 'coo traceproduct whole I')
  if(allocated(error))return

  !Case 2: whole-matrix block, deterministic dense b
  b = 0._wp
  do id = 0, n*n - 1
    b(mod(id, n) + 1, id/n + 1) = real(id, kind=wp)
  enddo
  expected = sum(mat * transpose(b))
  result = coo%traceproduct(1, n, 1, n, b)
  tol = 1.e-3_wp * max(abs(expected), 1._wp)
  call check(error, abs(result - expected) < tol, more = 'coo traceproduct whole dense')
  if(allocated(error))return

  !Case 3: 2x2 block at (1:2, 1:2)
  b(1, 1) = 1._wp; b(1, 2) = 2._wp; b(2, 1) = 3._wp; b(2, 2) = 4._wp
  expected = mat(1,1)*b(1,1) + mat(1,2)*b(2,1) &
           + mat(2,1)*b(1,2) + mat(2,2)*b(2,2)
  result = coo%traceproduct(1, 2, 1, 2, b(1:2, 1:2))
  tol = 1.e-3_wp * max(abs(expected), 1._wp)
  call check(error, abs(result - expected) < tol, more = 'coo traceproduct 2x2 diag')
  if(allocated(error))return

  !Case 4: 1x1 off-diagonal at (1,3): A(1,3)=13
  expected = mat(1,3) * 7._wp
  result = coo%traceproduct(1, 1, 3, 3, reshape([7._wp], [1,1]))
  tol = 1.e-3_wp * max(abs(expected), 1._wp)
  call check(error, abs(result - expected) < tol, more = 'coo traceproduct 1x1 off')
  if(allocated(error))return

  !Case 5: empty matrix -> 0
  coo = coosparse(4, unlog = sparse_unit)
  b(1:4, 1:4) = 0._wp
  do id = 1, 4
    b(id, id) = 1._wp
  enddo
  call check(error, abs(coo%traceproduct(1, 4, 1, 4, b(1:4, 1:4))) < tol_wp&
           , more = 'coo traceproduct empty')

end subroutine test_traceproduct_coo



end module modtest_coo

module modtest_coo
 use, intrinsic :: iso_fortran_env, only: int64, real64, output_unit
 use testdrive, only: new_unittest, unittest_type, error_type, check
 use modsparse, only: coosparse
 implicit none
 private

 public :: collect_coo

 integer, parameter :: sparse_unit = 555
 logical, parameter :: verbose = .false.
 
 integer, parameter :: ia(12) = [1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 5, 6]
 integer, parameter :: ja(12) = [1, 3, 6, 5, 1, 3, 4, 4, 1, 2, 5, 6]
 real(real64), parameter :: a(12) = [real(real64):: 11, 13, 22, 25, 31, 33, 34&
                                         , 44, 51, 52, 55, 66]
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
    , new_unittest("coo diag vect", test_diag_vect) &
    , new_unittest("coo diag vect lupper", test_diag_vect_lupper) &
    , new_unittest("coo ncol diag vect", test_ncol_diag_vect) &
    , new_unittest("coo ncol diag vect lupper", test_ncol_diag_vect_lupper) &
    , new_unittest("coo get", test_get) &
    , new_unittest("coo get nel", test_get_nel) &
    , new_unittest("coo get lupper", test_get_lupper) &
    , new_unittest("coo ncol get", test_ncol_get) &
    , new_unittest("coo ncol get nel", test_ncol_get_nel) &
    , new_unittest("coo ncol get lupper", test_ncol_get_lupper) &
    , new_unittest("coo multbyv", test_multbyv) &
    , new_unittest("coo multbyv lupper", test_multbyv_lupper) &
    , new_unittest("coo ncol multbyv", test_ncol_multbyv) &
    , new_unittest("coo ncol multbyv lupper", test_ncol_multbyv_lupper) &
    ]

  !to check: diag_mat

end subroutine collect_coo

!CONSTRUCTOR + GETDIM
subroutine test_constructor(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 type(coosparse)::coo

 coo = coosparse(nrow)
 call check(error, coo%getdim(1), nrow, more = 'getdim(1)')
 call check(error, coo%getdim(2), nrow, more = 'getdim(2)')
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
 real(real64), allocatable :: at(:)
 real(real64) :: mat(nrow, ncol)
 type(coosparse) :: coo

 coo = coosparse(nrow, unlog = sparse_unit)

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

subroutine test_add_nel(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol

 integer, allocatable :: iat(:), jat(:)
 real(real64), allocatable :: at(:)
 real(real64) :: mat(nrow, ncol)
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
 real(real64), allocatable :: at(:)
 real(real64) :: mat(nrow, ncol)
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
 real(real64), allocatable :: at(:)
 real(real64) :: mat(nrow, ncol)
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
 real(real64), allocatable :: at(:)
 real(real64) :: mat(nrow, ncol)
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
 real(real64), allocatable :: at(:)
 real(real64) :: mat(nrow, ncol)
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

!DIAG VECT
subroutine test_diag_vect(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia.eq.ja
 integer, parameter :: sizediag =min(nrow, ncol)

 real(real64) :: diagcheck(sizediag)
 real(real64), allocatable :: diag(:)
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

 real(real64) :: diagcheck(sizediag)
 real(real64), allocatable :: diag(:)
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

 real(real64) :: diagcheck(sizediag)
 real(real64), allocatable :: diag(:)
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

 real(real64) :: diagcheck(sizediag)
 real(real64), allocatable :: diag(:)
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
subroutine test_multbyv(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i

 character(len=1), parameter :: trans = 'n'
 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol
 real(real64), parameter :: alpha = 0.3_real64
 real(real64), parameter :: val = 10000._real64
 real(real64), parameter :: x(merge(ncol, nrow, trans == 'n')) = &
                              [(real(i, real64), i = 1, merge(ncol, nrow, trans == 'n'))]

 real(real64) :: y(merge(nrow, ncol, trans == 'n'))
 real(real64) :: ycheck(merge(nrow, ncol, trans == 'n'))
 type(coosparse) :: coo

 coo = coosparse(nrow, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a)

 y = 1
 call coo%mult(alpha, trans, x, val, y)

 ycheck = 1
 ycheck = ycheck * val + alpha * matmul(&
   merge(matcheck(nrow, ncol, ia, ja, a, lvalid), &
         transpose(matcheck(nrow, ncol, ia, ja, a, lvalid)), trans == 'n')&
   , x)
 
 call check(error, all(y == ycheck), 'multbyv')

end subroutine

subroutine test_multbyv_lupper(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i

 character(len=1), parameter :: trans = 'n'
 integer, parameter :: nrow = 5
 integer, parameter :: ncol = nrow
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia .le. ja
 real(real64), parameter :: alpha = 0.3_real64
 real(real64), parameter :: val = 10000._real64
 real(real64), parameter :: x(merge(ncol, nrow, trans == 'n')) = &
                              [(real(i, real64), i = 1, merge(ncol, nrow, trans == 'n'))]

 real(real64) :: y(merge(nrow, ncol, trans == 'n'))
 real(real64) :: ycheck(merge(nrow, ncol, trans == 'n'))
 type(coosparse) :: coo

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a)

 y = 1
 call coo%mult(alpha, trans, x, val, y)

 ycheck = 1
 ycheck = ycheck * val + alpha * matmul(&
   merge(matcheck(nrow, ncol, ia, ja, a, lvalid), &
         transpose(matcheck(nrow, ncol, ia, ja, a, lvalid)), trans == 'n')&
   , x)
 
 call check(error, all(y == ycheck), 'multbyv')

end subroutine

subroutine test_ncol_multbyv(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i

 character(len=1), parameter :: trans = 'n'
 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol
 real(real64), parameter :: alpha = 0.3_real64
 real(real64), parameter :: val = 10000._real64
 real(real64), parameter :: x(merge(ncol, nrow, trans == 'n')) = &
                              [(real(i, real64), i = 1, merge(ncol, nrow, trans == 'n'))]

 real(real64) :: y(merge(nrow, ncol, trans == 'n'))
 real(real64) :: ycheck(merge(nrow, ncol, trans == 'n'))
 type(coosparse) :: coo

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a)

 y = 1
 call coo%mult(alpha, trans, x, val, y)

 ycheck = 1
 ycheck = ycheck * val + alpha * matmul(&
   merge(matcheck(nrow, ncol, ia, ja, a, lvalid), &
         transpose(matcheck(nrow, ncol, ia, ja, a, lvalid)), trans == 'n')&
   , x)
 
 call check(error, all(y == ycheck), 'multbyv')

end subroutine

subroutine test_ncol_multbyv_lupper(error)
 type(error_type), allocatable, intent(out) :: error

 integer :: i

 character(len=1), parameter :: trans = 'n'
 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 logical, parameter :: lvalid(size(ia)) = ia.le.nrow .and. ja.le.ncol .and. ia .le. ja
 real(real64), parameter :: alpha = 0.3_real64
 real(real64), parameter :: val = 10000._real64
 real(real64), parameter :: x(merge(ncol, nrow, trans == 'n')) = &
                              [(real(i, real64), i = 1, merge(ncol, nrow, trans == 'n'))]

 real(real64) :: y(merge(nrow, ncol, trans == 'n'))
 real(real64) :: ycheck(merge(nrow, ncol, trans == 'n'))
 type(coosparse) :: coo

 coo = coosparse(nrow, n = ncol, lupper = .true., unlog = sparse_unit)

 call addval(coo, coo%getdim(1), coo%getdim(2),  ia, ja, a)

 y = 1
 call coo%mult(alpha, trans, x, val, y)

 ycheck = 1
 ycheck = ycheck * val + alpha * matmul(&
   merge(matcheck(nrow, ncol, ia, ja, a, lvalid), &
         transpose(matcheck(nrow, ncol, ia, ja, a, lvalid)), trans == 'n')&
   , x)
 
 call check(error, all(y == ycheck), 'multbyv')

end subroutine


!INTERNAL
subroutine addval(coo, nrow, ncol, ia, ja, a, iat, jat, at, mat)
 type(coosparse), intent(inout) :: coo
 integer, intent(in) :: nrow, ncol
 integer, intent(in) :: ia(:), ja(:) 
 real(real64), intent(in) :: a(:)
 integer, allocatable , intent(out), optional :: iat(:), jat(:) 
 real(real64), allocatable , intent(out), optional :: at(:)
 real(real64), intent(out), optional :: mat(:,:)

 integer :: i, j
 real(real64) :: val

 if(present(mat))then
  if(nrow.ne.size(mat,1) .or. ncol.ne.size(mat,2))return
 endif

 do i = 1, size(ia)
  call coo%add(ia(i), ja(i), a(i))
 enddo

 if(.not.present(iat).or..not.present(jat).or..not.present(at).or..not.present(mat))return

 allocate(iat(0), jat(0), at(0))
 do i = 1, maxval(ia)
  do j = 1, maxval(ja)
   val = coo%get(i, j)
   if(present(mat))then
    if(i.le.nrow.and.j.le.ncol)mat(i, j) = val
   endif
   if (val .ne. 0._real64) then
    iat = [iat, i]
    jat = [jat, j]
    at = [at, val]
   endif
  enddo
 enddo

end subroutine

subroutine printmat(mat)
 real(real64), intent(in) :: mat(:,:)

 integer :: i

 write(output_unit, '(a)')repeat('*',size(mat, 1))

 do i = 1, size(mat, 1)
  write(output_unit, '(*(f0.2,1x))')mat(i,:)
 enddo

 write(output_unit, '(a)')repeat('*',size(mat, 1))
end subroutine

function matcheck(nrow, ncol, ia, ja, a, lvalid) result(mat)
 integer, intent(in) :: nrow, ncol, ia(:), ja(:)
 real(real64), intent(in) :: a(:)
 logical, intent(in) :: lvalid(:)

 real(real64) :: mat(nrow, ncol)
 integer :: i

 mat = 0
 do i = 1, size(ia)
  if(lvalid(i))then
   mat(ia(i), ja(i)) = a(i)
  endif
 enddo

end function

function getmat(coo) result(mat)
 type(coosparse), intent(inout) :: coo

 real(real64) :: mat(coo%getdim(1), coo%getdim(2))

 integer :: i, j
 real(real64) :: val

 mat = -1

 do i = 1, coo%getdim(1)
  do j = 1, coo%getdim(2)
   mat(i, j) = coo%get(i, j)
  enddo
 enddo

end function

end module modtest_coo

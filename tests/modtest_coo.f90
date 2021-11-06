module modtest_coo
 use, intrinsic :: iso_fortran_env, only: int64, real64, output_unit
 use testdrive, only: new_unittest, unittest_type, error_type, check
 use modsparse, only: coosparse
 implicit none
 private

 public :: collect_coo

 integer, parameter :: sparse_unit = 555
 logical, parameter :: verbose = .false.
 
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
    , new_unittest("coo ncol nel", test_ncol_add_nel) &
    , new_unittest("coo ncol lupper", test_ncol_add_lupper) &
    ]

end subroutine collect_coo

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

subroutine test_add(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ia(12) = [1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 5, 6]
 integer, parameter :: ja(12) = [1, 3, 6, 5, 1, 3, 4, 4, 1, 2, 5, 6]
 real(real64), parameter :: a(12) = [real(real64):: 11, 13, 22, 25, 31, 33, 34&
                                         , 44, 51, 52, 55, 66]
 integer :: i, j
 integer, allocatable :: iat(:), jat(:)
 real(real64), allocatable :: at(:)
 real(real64) :: val
 real(real64) :: mat(nrow, nrow)
 logical, allocatable :: lvalid(:)
 type(coosparse) :: coo

 coo = coosparse(nrow, unlog = sparse_unit)

 do i = 1, size(ia)
  call coo%add(ia(i), ja(i), a(i))
 enddo

 allocate(iat(0), jat(0), at(0))
 do i = 1, maxval(ia)
  do j = 1, maxval(ja)
   val = coo%get(i, j)
   if(i.le.nrow.and.j.le.nrow)mat(i, j) = val
   if (val .ne. 0._real64) then
    iat = [iat, i]
    jat = [jat, j]
    at = [at, val]
   endif
  enddo
 enddo

 lvalid = ia.le.nrow .and. ja.le.nrow

 call check(error, size(iat), size(pack(ia, lvalid)), more = 'size(iat)')
 call check(error, size(jat), size(pack(ja, lvalid)), more = 'size(jat)')
 call check(error, size(at), size(pack(a, lvalid)), more = 'size(at)')
 if (allocated(error)) return

 call check(error, all(iat == pack(ia, lvalid)), 'ia')
 call check(error, all(jat == pack(ja, lvalid)), 'ja')
 call check(error, all(at == pack(a, lvalid)), 'a')

 if(verbose)call printmat(mat)

end subroutine

subroutine test_add_nel(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ia(12) = [1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 5, 6]
 integer, parameter :: ja(12) = [1, 3, 6, 5, 1, 3, 4, 4, 1, 2, 5, 6]
 real(real64), parameter :: a(12) = [real(real64):: 11, 13, 22, 25, 31, 33, 34&
                                         , 44, 51, 52, 55, 66]
 integer :: i, j
 integer, allocatable :: iat(:), jat(:)
 real(real64), allocatable :: at(:)
 real(real64) :: val
 real(real64) :: mat(nrow, nrow)
 logical, allocatable :: lvalid(:)
 type(coosparse) :: coo

 coo = coosparse(nrow, nel = 4_int64, unlog = sparse_unit)

 do i = 1, size(ia)
  call coo%add(ia(i), ja(i), a(i))
 enddo

 allocate(iat(0), jat(0), at(0))
 do i = 1, maxval(ia)
  do j = 1, maxval(ja)
   val = coo%get(i, j)
   if(i.le.nrow.and.j.le.nrow)mat(i, j) = val
   if (val .ne. 0._real64) then
    iat = [iat, i]
    jat = [jat, j]
    at = [at, val]
   endif
  enddo
 enddo

 lvalid = ia.le.nrow .and. ja.le.nrow

 call check(error, size(iat), size(pack(ia, lvalid)), more = 'size(iat)')
 call check(error, size(jat), size(pack(ja, lvalid)), more = 'size(jat)')
 call check(error, size(at), size(pack(a, lvalid)), more = 'size(at)')
 if (allocated(error)) return

 call check(error, all(iat == pack(ia, lvalid)), 'ia')
 call check(error, all(jat == pack(ja, lvalid)), 'ja')
 call check(error, all(at == pack(a, lvalid)), 'a')

 if(verbose)call printmat(mat)

end subroutine

subroutine test_add_lupper(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ia(12) = [1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 5, 6]
 integer, parameter :: ja(12) = [1, 3, 6, 5, 1, 3, 4, 4, 1, 2, 5, 6]
 real(real64), parameter :: a(12) = [real(real64):: 11, 13, 22, 25, 31, 33, 34&
                                         , 44, 51, 52, 55, 66]
 integer :: i, j
 integer, allocatable :: iat(:), jat(:)
 real(real64), allocatable :: at(:)
 real(real64) :: val
 real(real64) :: mat(nrow, nrow)
 logical, allocatable :: lvalid(:)
 type(coosparse) :: coo

 coo = coosparse(nrow, lupper = .true., unlog = sparse_unit)

 do i = 1, size(ia)
  call coo%add(ia(i), ja(i), a(i))
 enddo

 allocate(iat(0), jat(0), at(0))
 do i = 1, maxval(ia)
  do j = 1, maxval(ja)
   val = coo%get(i, j)
   if(i.le.nrow.and.j.le.nrow)mat(i, j) = val
   if (val .ne. 0._real64) then
    iat = [iat, i]
    jat = [jat, j]
    at = [at, val]
   endif
  enddo
 enddo

 lvalid = ia.le.nrow .and. ja.le.nrow .and. ia .le. ja

 call check(error, size(iat), size(pack(ia, lvalid)), more = 'size(iat)')
 call check(error, size(jat), size(pack(ja, lvalid)), more = 'size(jat)')
 call check(error, size(at), size(pack(a, lvalid)), more = 'size(at)')
 if (allocated(error)) return

 call check(error, all(iat == pack(ia, lvalid)), 'ia')
 call check(error, all(jat == pack(ja, lvalid)), 'ja')
 call check(error, all(at == pack(a, lvalid)), 'a')

 if(verbose)call printmat(mat)

end subroutine

subroutine test_ncol_add(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 integer, parameter :: ia(12) = [1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 5, 6]
 integer, parameter :: ja(12) = [1, 3, 6, 5, 1, 3, 4, 4, 1, 2, 5, 6]
 real(real64), parameter :: a(12) = [real(real64):: 11, 13, 22, 25, 31, 33, 34&
                                         , 44, 51, 52, 55, 66]
 integer :: i, j
 integer, allocatable :: iat(:), jat(:)
 real(real64), allocatable :: at(:)
 real(real64) :: val
 real(real64) :: mat(nrow, ncol)
 logical, allocatable :: lvalid(:)
 type(coosparse) :: coo

 coo = coosparse(nrow, n = ncol, unlog = sparse_unit)

 do i = 1, size(ia)
  call coo%add(ia(i), ja(i), a(i))
 enddo

 allocate(iat(0), jat(0), at(0))
 do i = 1, maxval(ia)
  do j = 1, maxval(ja)
   val = coo%get(i, j)
   if(i.le.nrow.and.j.le.ncol)mat(i, j) = val
   if (val .ne. 0._real64) then
    iat = [iat, i]
    jat = [jat, j]
    at = [at, val]
   endif
  enddo
 enddo

 lvalid = ia.le.nrow .and. ja.le.ncol

 call check(error, size(iat), size(pack(ia, lvalid)), more = 'size(iat)')
 call check(error, size(jat), size(pack(ja, lvalid)), more = 'size(jat)')
 call check(error, size(at), size(pack(a, lvalid)), more = 'size(at)')
 if (allocated(error)) return

 call check(error, all(iat == pack(ia, lvalid)), 'ia')
 call check(error, all(jat == pack(ja, lvalid)), 'ja')
 call check(error, all(at == pack(a, lvalid)), 'a')

 if(verbose)call printmat(mat)

end subroutine

subroutine test_ncol_add_nel(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 integer, parameter :: ia(12) = [1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 5, 6]
 integer, parameter :: ja(12) = [1, 3, 6, 5, 1, 3, 4, 4, 1, 2, 5, 6]
 real(real64), parameter :: a(12) = [real(real64):: 11, 13, 22, 25, 31, 33, 34&
                                         , 44, 51, 52, 55, 66]
 integer :: i, j
 integer, allocatable :: iat(:), jat(:)
 real(real64), allocatable :: at(:)
 real(real64) :: val
 real(real64) :: mat(nrow, ncol)
 logical, allocatable :: lvalid(:)
 type(coosparse) :: coo

 coo = coosparse(nrow, n = ncol, nel = 2_int64,  unlog = sparse_unit)

 do i = 1, size(ia)
  call coo%add(ia(i), ja(i), a(i))
 enddo

 allocate(iat(0), jat(0), at(0))
 do i = 1, maxval(ia)
  do j = 1, maxval(ja)
   val = coo%get(i, j)
   if(i.le.nrow.and.j.le.ncol)mat(i, j) = val
   if (val .ne. 0._real64) then
    iat = [iat, i]
    jat = [jat, j]
    at = [at, val]
   endif
  enddo
 enddo

 lvalid = ia.le.nrow .and. ja.le.ncol

 call check(error, size(iat), size(pack(ia, lvalid)), more = 'size(iat)')
 call check(error, size(jat), size(pack(ja, lvalid)), more = 'size(jat)')
 call check(error, size(at), size(pack(a, lvalid)), more = 'size(at)')
 if (allocated(error)) return

 call check(error, all(iat == pack(ia, lvalid)), 'ia')
 call check(error, all(jat == pack(ja, lvalid)), 'ja')
 call check(error, all(at == pack(a, lvalid)), 'a')

 if(verbose)call printmat(mat)

end subroutine

subroutine test_ncol_add_lupper(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: nrow = 5
 integer, parameter :: ncol = 4
 integer, parameter :: ia(12) = [1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 5, 6]
 integer, parameter :: ja(12) = [1, 3, 6, 5, 1, 3, 4, 4, 1, 2, 5, 6]
 real(real64), parameter :: a(12) = [real(real64):: 11, 13, 22, 25, 31, 33, 34&
                                         , 44, 51, 52, 55, 66]
 integer :: i, j
 integer, allocatable :: iat(:), jat(:)
 real(real64), allocatable :: at(:)
 real(real64) :: val
 real(real64) :: mat(nrow, ncol)
 logical, allocatable :: lvalid(:)
 type(coosparse) :: coo

 coo = coosparse(nrow, n = ncol, lupper = .true.,  unlog = sparse_unit)

 do i = 1, size(ia)
  call coo%add(ia(i), ja(i), a(i))
 enddo

 allocate(iat(0), jat(0), at(0))
 do i = 1, maxval(ia)
  do j = 1, maxval(ja)
   val = coo%get(i, j)
   if(i.le.nrow.and.j.le.ncol)mat(i, j) = val
   if (val .ne. 0._real64) then
    iat = [iat, i]
    jat = [jat, j]
    at = [at, val]
   endif
  enddo
 enddo

 lvalid = ia.le.nrow .and. ja.le.ncol .and. ia.le.ja

 call check(error, size(iat), size(pack(ia, lvalid)), more = 'size(iat)')
 call check(error, size(jat), size(pack(ja, lvalid)), more = 'size(jat)')
 call check(error, size(at), size(pack(a, lvalid)), more = 'size(at)')
 if (allocated(error)) return

 call check(error, all(iat == pack(ia, lvalid)), 'ia')
 call check(error, all(jat == pack(ja, lvalid)), 'ja')
 call check(error, all(at == pack(a, lvalid)), 'a')

 if(verbose)call printmat(mat)

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

end module modtest_coo

module modtest_random
#if (_DP==0)
 use, intrinsic :: iso_fortran_env, only: int64, wp => real32
#else
 use, intrinsic :: iso_fortran_env, only: int64, wp => real64
#endif
 use testdrive, only: new_unittest, unittest_type, error_type, check
 use modrandom, only: setseed, rand_stdnormal
 use modtest_common, only: tol_wp, verbose
 implicit none
 private

 public :: collect_random

 integer, parameter :: sparse_unit = 556
 
contains

!> Collect all exported unit tests
subroutine collect_random(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  !how to test: getordering, print, printsquare, save, sort

  testsuite = [ &
    new_unittest("stdnormal", test_stdnormal) &
    ]

end subroutine collect_random

!CONSTRUCTOR + GETDIM + ISSQUARE
subroutine test_stdnormal(error)
 type(error_type), allocatable, intent(out) :: error

 integer, parameter :: n = 100000
 integer :: i
 real(wp) :: vector(n)
 real(wp) :: mean, sd

 call setseed(1)
 
 vector = [(rand_stdnormal(), i = 1, n)]

 mean = sum(vector) / n
 sd = sqrt(sum((vector - mean)**2)/n)

 call check(error, abs(int(mean * 1000)/1000._wp - 0._wp) < tol_wp, more = 'Mean different than 0')
 if(allocated(error))return

 call check(error, abs(int(sd * 1000)/1000._wp - 1._wp) < tol_wp, more = 'SD different than 1')
 if(allocated(error))return

end subroutine


end module

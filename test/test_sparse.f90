program test_sparse
  use, intrinsic :: iso_fortran_env, only : error_unit, output_unit
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type
  use modtest_coo, only : collect_coo
  use modtest_crs, only : collect_crs
  use modtest_random, only : collect_random
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0

  testsuites = [ &
    new_testsuite("modtest_coo", collect_coo) &
    , new_testsuite("modtest_crs", collect_crs) &
    , new_testsuite("modtest_random", collect_random) &
    ]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat, parallel = .false.)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop
  else
    write(output_unit, '(a)') 'All tests passed!'
  end if

end program test_sparse

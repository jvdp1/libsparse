program cholesky
 use, intrinsic :: iso_fortran_env, only: int32, real64
 use modsparse, only: coosparse, crssparse, assignment(=)
 implicit none
 integer :: coosize
 integer :: i, un, io, row, col
 integer, allocatable :: perm(:)
 real(real64) :: val
 character(len=100) :: cdummy
 character(len=:), allocatable :: coofile
 type(coosparse) :: coo
 type(crssparse) :: crs
 logical, parameter :: ltest = .true.

 call get_command_argument(1, cdummy)
 coofile = trim(cdummy)

 call get_command_argument(2, cdummy)
 read(cdummy,*)coosize

 coo = coosparse(coosize, lupper=.true.)
 call coo%setsymmetric()

 open(newunit=un, file=coofile, status='old', action='read')
 do
  read(un, *, iostat=io)row, col, val
  if(io.ne.0)exit
  call coo%add(row, col, val)
 enddo
 close(un)

 !CRS SPARSE
 crs = coo

#if (_METIS==1)
 perm = crs%getordering()
#else
 perm = [(i, i=1, crs%getdim(1))]
#endif
 call crs%setpermutation(perm)

 call crs%printstats()

 call crs%chol()

 call crs%setsymmetric(.false.)
 call crs%printstats()

 !print Cholesky factor and permutation vector
 call crs%printtofile(coofile//'_chol')

 open(newunit=un, file = coofile//'_perm', action = 'write')
 do i = 1, size(perm)
  write(un, '(i0)')perm(i)
 enddo
 close(un)


 !test
 if(ltest)then
  if(test_cholesky(coo, perm, crs))then
   write(*,'(a)')'Cholesky solver OK'
  else
   write(*,'(a)')'Cholesky solver NOT OK'
  endif
 endif

contains

function test_cholesky(coo, perm, crschol) result(ltest)
 integer, intent(in) :: perm(:)
 type(coosparse), intent(in) :: coo
 type(crssparse), intent(in) :: crschol
 logical :: ltest

 real(real64), parameter :: tol = 1.d-10
 real(real64), allocatable :: rhs(:)
 real(real64), allocatable :: sol1(:)
 real(real64), allocatable :: sol2(:)
 type(crssparse) :: crs

 ltest = .false.

 allocate(rhs(crschol%getdim(1)))
 call random_number(rhs)

 !cholesky
 allocate(sol1(crschol%getdim(1)))
 call crschol%isolve(sol1, rhs)

 !pardiso
 allocate(sol2(coo%getdim(1)))
 crs = coo

 call crs%setpermutation(perm)
 call crs%solve(sol2, rhs)

 !test
 ltest = all(abs(sol1-sol2) < tol )

end function

end program

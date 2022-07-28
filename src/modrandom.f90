module modrandom
#if (_DP==0)
 use, intrinsic:: iso_fortran_env, only: wp=>real32
#else
 use, intrinsic:: iso_fortran_env, only: wp=>real64
#endif
 implicit none (type, external)
 private
 public :: setseed
 public :: rand_stduniform
 public :: rand_stdnormal

 real(wp), parameter :: pi = 4._wp * atan(1._wp)

contains

subroutine setseed(iseed)
 integer, intent(in) :: iseed

 integer :: n
 integer, allocatable :: seed(:)

 call random_seed(size = n)
 allocate(seed(n))
 seed = iseed
 call random_seed(put = seed)
 deallocate(seed)

end subroutine

function rand_stduniform() result(var)
 real(wp) :: var
 
 call random_number(var)
 var = 1 - var

end function

function rand_stdnormal() result(var)
 real(wp) :: var

 real(wp) :: u1, u2

 u1 = rand_stduniform()
 u2 = rand_stduniform()
 
 var = sqrt(-2*log(u1)) * cos(2*pi*u2)

end function


end module

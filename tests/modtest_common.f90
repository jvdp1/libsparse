module modtest_common
 use, intrinsic :: iso_fortran_env, only: real64, output_unit, wp => real64
 use modsparse, only: coosparse
 implicit none
 private
 public :: addval_coo, matcheck, printmat

 real(wp), parameter, public :: tol_wp = epsilon(1._wp) * 10**4

 logical, parameter, public :: verbose = .false.
 
 integer, parameter, public :: ia(12) = [1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 5, 6]
 integer, parameter, public :: ja(12) = [1, 3, 6, 5, 1, 3, 4, 4, 1, 2, 5, 6]
 real(wp), parameter, public :: a(12) = [real(wp):: 11, 13, 22, 25, 31, 33, 34&
                                         , 44, 51, 52, 55, 66]
contains
subroutine addval_coo(coo, nrow, ncol, ia, ja, a, iat, jat, at, mat)
 type(coosparse), intent(inout) :: coo
 integer, intent(in) :: nrow, ncol
 integer, intent(in) :: ia(:), ja(:) 
 real(wp), intent(in) :: a(:)
 integer, allocatable , intent(out), optional :: iat(:), jat(:) 
 real(wp), allocatable , intent(out), optional :: at(:)
 real(wp), intent(out), optional :: mat(:,:)

 integer :: i, j
 real(wp) :: val

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
   if (val .ne. 0._wp) then
    iat = [iat, i]
    jat = [jat, j]
    at = [at, val]
   endif
  enddo
 enddo

end subroutine

pure function matcheck(nrow, ncol, ia, ja, a, lvalid) result(mat)
 integer, intent(in) :: nrow, ncol, ia(:), ja(:)
 real(wp), intent(in) :: a(:)
 logical, intent(in) :: lvalid(:)

 real(wp) :: mat(nrow, ncol)
 integer :: i

 mat = 0
 do i = 1, size(ia)
  if(lvalid(i))then
   mat(ia(i), ja(i)) = a(i)
  endif
 enddo

end function

subroutine printmat(mat)
 real(wp), intent(in) :: mat(:,:)

 integer :: i

 write(output_unit, '(a)')repeat('*',size(mat, 1))

 do i = 1, size(mat, 1)
  write(output_unit, '(*(f0.2,1x))')mat(i,:)
 enddo

 write(output_unit, '(a)')repeat('*',size(mat, 1))
end subroutine

end module

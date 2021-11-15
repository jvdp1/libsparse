module modtest_common
#if (_DP==0)
 use, intrinsic :: iso_fortran_env, only: output_unit, wp => real32
#else
 use, intrinsic :: iso_fortran_env, only: output_unit, wp => real64
#endif
 use modsparse, only: gen_sparse, coosparse
 implicit none
 private
 public :: addval_coo, getmat, matcheck, printmat

 real(wp), parameter, public :: tol_wp = epsilon(1._wp) * 10**4

 logical, parameter, public :: verbose = .false.
 
 integer, parameter, public :: ia(16) = [1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6]
 integer, parameter, public :: ja(16) = [1, 3, 4, 5, 6, 6, 5, 1, 3, 4, 4, 6, 1, 2, 5, 6]
 real(wp), parameter, public :: a(16) = [real(wp):: 11, 13, 14, 15, 16, 22, 25, 31, 33, 34&
                                         , 44, 46, 51, 52, 55, 66]
 real(wp), parameter, public :: aspsd(16) = [real(wp):: 101, 13, 14, 15, 16, 0, 0, 31, 303, 34&
                                            , 404, 46, 51, 52, 505, 606]

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

function getmat(sparse) result(mat)
 class(gen_sparse), intent(inout) :: sparse

 real(wp) :: mat(sparse%getdim(1), sparse%getdim(2))

 integer :: i, j

 mat = -1

 do i = 1, sparse%getdim(1)
  do j = 1, sparse%getdim(2)
   mat(i, j) = sparse%get(i, j)
  enddo
 enddo

end function

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
  write(output_unit, '(*(g0.5,1x))')mat(i,:)
 enddo

 write(output_unit, '(a)')repeat('*',size(mat, 1))
end subroutine

end module

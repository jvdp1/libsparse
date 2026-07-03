!> Module containing various helpers
!> Most likely all of them inefficient but should do the job

module modsparse_helpers
#if (_DP==0)
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real32
#else
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real64
#endif
 !$ use omp_lib
 implicit none(type, external)
 private
 public :: csrsymv
 public :: csrmm
 public :: csrmv
 public :: csrtrsv

 contains

! Symmetric CSR matrix-vector product: y = A*x
subroutine csrsymv(uplo, n, a, ia, ja, x, y)
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 integer(kind=int32), intent(in) :: ia(:)
 integer(kind=int32), intent(in) :: ja(:)
 real(kind=wp), intent(in) :: a(:)
 real(kind=wp), intent(in) :: x(:)
 real(kind=wp), intent(out) :: y(:)

 integer :: i, j

 y = 0._wp

 select case(uplo)
 case ('U', 'u')
  do i = 1, n
   do j = ia(i), ia(i+1)-1
    y(i) = y(i) + a(j) * x(ja(j))
    if (ja(j) /= i) y(ja(j)) = y(ja(j)) + a(j) * x(i)
   enddo
  enddo
 case ('L', 'l')
  do i = 1, n
   do j = ia(i), ia(i+1)-1
    y(i) = y(i) + a(j) * x(ja(j))
    if (ja(j) /= i) y(ja(j)) = y(ja(j)) + a(j) * x(i)
   enddo
  enddo
 case default
  error stop 'csrsymv: unsupported uplo'
 end select

end subroutine csrsymv

! CSR matrix-dense matrix product: c = beta*c + alpha*op(A)*b
subroutine csrmm(transa, m, n, k, alpha, matdescra, a, ja, pntrb, pntre, b, ldb, beta, c, ldc)
 character(len=1), intent(in) :: transa
 integer, intent(in) :: m
 integer, intent(in) :: n
 integer, intent(in) :: k
 real(kind=wp), intent(in) :: alpha
 character(len=1), intent(in) :: matdescra(6)
 real(kind=wp), intent(in) :: a(:)
 integer(kind=int32), intent(in) :: ja(:)
 integer(kind=int32), intent(in) :: pntrb(:)
 integer(kind=int32), intent(in) :: pntre(:)
 real(kind=wp), intent(in) :: b(:,:)
 integer, intent(in) :: ldb
 real(kind=wp), intent(in) :: beta
 real(kind=wp), intent(inout) :: c(:,:)
 integer, intent(in) :: ldc

 integer :: i, j, p

 if (transa == 'N' .or. transa == 'n') then

  c(1:m, 1:n) = beta * c(1:m, 1:n)

  select case (matdescra(1))
  case ('G', 'g', 'T', 't')
   do j = 1, n
    do i = 1, m
     do p = pntrb(i), pntre(i)-1
      c(i,j) = c(i,j) + alpha * a(p) * b(ja(p),j)
     end do
    end do
   end do
  case ('S', 's')
   do j = 1, n
    do i = 1, m
     do p = pntrb(i), pntre(i)-1
      c(i,j) = c(i,j) + alpha * a(p) * b(ja(p),j)
      if (ja(p) /= i) c(ja(p),j) = c(ja(p),j) + alpha * a(p) * b(i,j)
     end do
    end do
   end do
  case default
   error stop 'csrmm: unsupported matdescra(1)'
  end select

 else

  select case (matdescra(1))
  case ('S', 's')
   ! Symmetric: A^T = A
   c(1:m, 1:n) = beta * c(1:m, 1:n)
   do j = 1, n
    do i = 1, m
     do p = pntrb(i), pntre(i)-1
      c(i,j) = c(i,j) + alpha * a(p) * b(ja(p),j)
      if (ja(p) /= i) c(ja(p),j) = c(ja(p),j) + alpha * a(p) * b(i,j)
     end do
    end do
   end do
  case ('G', 'g', 'T', 't')
   c(1:k, 1:n) = beta * c(1:k, 1:n)
   do j = 1, n
    do i = 1, m
     do p = pntrb(i), pntre(i)-1
      c(ja(p),j) = c(ja(p),j) + alpha * a(p) * b(i,j)
     end do
    end do
   end do
  case default
   error stop 'csrmm: unsupported matdescra(1)'
  end select

 end if

end subroutine csrmm

! General CSR matrix-vector product: y = beta*y + alpha*op(A)*x
subroutine csrmv(transa, m, k, alpha, matdescra, val, indx, pntrb, pntre, x, beta, y)
 character(len=1), intent(in) :: transa
 integer, intent(in) :: m
 integer, intent(in) :: k
 real(kind=wp), intent(in) :: alpha
 real(kind=wp), intent(in) :: beta
 character(len=1), intent(in) :: matdescra(:)
 real(kind=wp), intent(in) :: val(:)
 real(kind=wp), intent(in) :: x(:)
 integer(kind=int32), intent(in) :: indx(:)
 integer(kind=int32), intent(in) :: pntrb(:)
 integer(kind=int32), intent(in) :: pntre(:)
 real(kind=wp), intent(inout) :: y(:)

 integer :: i, j

 if (transa == 'N' .or. transa == 'n') then
  y(1:m) = beta * y(1:m)
  select case (matdescra(1))
  case ('G', 'g', 'T', 't')
   do i = 1, m
    do j = pntrb(i), pntre(i)-1
     y(i) = y(i) + alpha * val(j) * x(indx(j))
    end do
   end do
  case ('S', 's')
   do i = 1, m
    do j = pntrb(i), pntre(i)-1
     y(i) = y(i) + alpha * val(j) * x(indx(j))
     if (indx(j) /= i) y(indx(j)) = y(indx(j)) + alpha * val(j) * x(i)
    end do
   end do
  case default
   error stop 'csrmv: unsupported matdescra(1)'
  end select
 else
  select case (matdescra(1))
  case ('G', 'g', 'T', 't')
   y(1:k) = beta * y(1:k)
   do i = 1, m
    do j = pntrb(i), pntre(i)-1
     y(indx(j)) = y(indx(j)) + alpha * val(j) * x(i)
    end do
   end do
  case ('S', 's')
   y(1:m) = beta * y(1:m)
   do i = 1, m
    do j = pntrb(i), pntre(i)-1
     y(i) = y(i) + alpha * val(j) * x(indx(j))
     if (indx(j) /= i) y(indx(j)) = y(indx(j)) + alpha * val(j) * x(i)
    end do
   end do
  case default
   error stop 'csrmv: unsupported matdescra(1)'
  end select
 end if

end subroutine csrmv

! Upper triangular CSR solve: op(A)*y = x
! Diagonal is first entry in each row (sorted CSR assumed)
subroutine csrtrsv(uplo, transa, diag, m, a, ia, ja, x, y)
 character(len=1), intent(in) :: uplo
 character(len=1), intent(in) :: transa
 character(len=1), intent(in) :: diag
 integer, intent(in) :: m
 real(kind=wp), intent(in) :: a(:)
 integer(kind=int32), intent(in) :: ia(:)
 integer(kind=int32), intent(in) :: ja(:)
 real(kind=wp), intent(in) :: x(:)
 real(kind=wp), intent(out) :: y(:)

 integer :: i, k
 logical :: ldiag

 y = x

 select case(uplo)
 case ('U', 'u')
  ldiag = (diag == 'N' .or. diag == 'n')

  if (transa == 'N' .or. transa == 'n') then
   ! Solve U*y = x; backward substitution
   do i = m, 1, -1
    do k = ia(i)+1, ia(i+1)-1
     y(i) = y(i) - a(k) * y(ja(k))
    end do
    if (ldiag) y(i) = y(i) / a(ia(i))
   end do
  else
   ! Solve U^T*y = x; forward substitution
   do i = 1, m
    if (ldiag) y(i) = y(i) / a(ia(i))
    do k = ia(i)+1, ia(i+1)-1
     y(ja(k)) = y(ja(k)) - a(k) * y(i)
    end do
   end do
  end if

 case default
  error stop 'internal error csrtrsv'
 end select

end subroutine csrtrsv

end module modsparse_helpers

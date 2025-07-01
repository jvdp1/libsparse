module modsparse_inv_int
 use, intrinsic :: iso_fortran_env, only: wp => real64
 implicit none
 private
 public :: spainverse
 public :: gsfct

 real(kind=wp), parameter :: tol = 1.e-8_wp

contains

!C#######################################################################
!C#      this is subroutine "GSFCT" as given in : 'Computer Solution    #
!C#      of Sparse Linear Systems' by A. George, J. Liu and             #
!C#      E. Ng, 1994, pp. 180-183;                                      #
!C#                                                                     #
!C# Rewritten towards Mordern Fortran                                   #
!C#                                                                     #
!C# Support of sparse inverse of SPSD matrices by implementing the      #
!C# S. D. Kachman modifications                                         #
!C# (https://www.ars.usda.gov/ARSUserFiles/80420530/MTDFREML/MTDFMan.pdf #
!C# ; Chapter 6)                                                        #
!C#                                                                     #
!C#######################################################################

!c----- subroutine gsfct
!c***************************************************************
!c***************************************************************
!c******     gsfct ..... general sparse symmetric fact     ******
!c***************************************************************
!c***************************************************************
!c
!c     purpose - this subroutine performs the symmetric
!c        factorization for a general sparse system, stored in
!c        the compressed subscript data format.
!c
!c     input parameters -
!c        neqns - number of equations.
!c        xlnz - index vector for lnz. xlnz(i) points to the
!c               start of nonzeros in column i of factor l.
!c        (xnzsub, nzsub) - the compressed subscript data
!c               structure for factor l.
!c
!c     updated parameters -
!c        lnz - on input, contains nonzeros of a, and on
!c              return, the nonzeros of l.
!c        diag - the diagonal of l overwrites that of a.
!c        iflag - the error flag. it is set to 1 if a zero or
!c               negative square root occurs during the
!c               factorization.
!c        ops (removed) - a double precision common parameter that is
!c              incremented by the number of operations
!c              performed by the subroutine.
!c
!c     working parameters -
!c        link - at step j, the list in
!c                  link(j), link(link(j)), ...........
!c               consists of those columns that will modify
!c               the column l(*,j).
!c        first - temporary vector to point to the first
!c                nonzero in each column that will be used
!c                next for modification.
!c        temp - a temporary vector to accumulate modifications.
!c
!c***************************************************************
!c
 subroutine gsfct(neqns, xlnz, lnz, xnzsub, nzsub, diag, iflag)
!c
!c***************************************************************
!c
  integer, intent(in) :: neqns
  integer, intent(in) :: nzsub(:)
  integer, intent(in) :: xlnz(:), xnzsub(:)
  integer, intent(out) :: iflag
  real(kind=wp), intent(inout) :: diag(:), lnz(:)

  real(kind=wp) :: diagj, ljk
  real(kind=wp), allocatable :: temp(:)
  integer, allocatable :: link(:), first(:)
  integer :: i, ii, istop, istrt, isub, j, k, kfirst, newk
!c
!c***************************************************************
!c
!c        ------------------------------
!c        initialize working vectors ...
!c        ------------------------------
  allocate (first(neqns), source=0)
  allocate (link(neqns), source=0)
  allocate (temp(neqns), source=0.0_wp)
  iflag = 0
!c        --------------------------------------------
!c        compute column l(*,j) for j = 1,...., neqns.
!c        --------------------------------------------
  do j = 1, neqns
!c           -------------------------------------------
!c           for each column l(*,k) that affects l(*,j).
!c           -------------------------------------------
   diagj = 0.0d0
   newk = link(j)
   k = newk
   do while (k /= 0)
    newk = link(k)
!c              ---------------------------------------
!c              outer product modification of l(*,j) by
!c              l(*,k) starting at first(k) of l(*,k).
!c              ---------------------------------------
    kfirst = first(k)
    ljk = lnz(kfirst)
    diagj = diagj + ljk*ljk
    istrt = kfirst + 1
    istop = xlnz(k + 1) - 1
    if (istop >= istrt) then
!c                 ------------------------------------------
!c                 before modification, update vectors first,
!c                 and link for future modification steps.
!c                 ------------------------------------------
     first(k) = istrt
     i = xnzsub(k) + (kfirst - xlnz(k)) + 1
     isub = nzsub(i)
     link(k) = link(isub)
     link(isub) = k
!c                 ---------------------------------------
!c                 the actual mod is saved in vector temp.
!c                 ---------------------------------------
     do ii = istrt, istop
      isub = nzsub(i)
      temp(isub) = temp(isub) + lnz(ii)*ljk
      i = i + 1
     end do
    end if
    k = newk
   end do
!c           ----------------------------------------------
!c           apply the modifications accumulated in temp to
!c           column l(*,j).
!c           ----------------------------------------------
   diagj = diag(j) - diagj
   if (diag(j) < tol) then
!c              ------------------------------------------------------
!c              error - zero diagonal element
!c              ------------------------------------------------------
    diagj = 0.d0
    diag(j) = 0.d0
    iflag = iflag + 1
   else
    if (diagj < tol*diag(j)) then
!c              ------------------------------------------------------
!c              error - zero or negative square root in factorization.
!c              ------------------------------------------------------
     diagj = 0.d0
     diag(j) = 0.d0
     iflag = iflag + 1
    else
     diagj = sqrt(diagj)
     diag(j) = diagj
     diagj = 1._wp/diagj
    end if
   end if
   istrt = xlnz(j)
   istop = xlnz(j + 1) - 1
   if (istop >= istrt) then
    first(j) = istrt
    i = xnzsub(j)
    isub = nzsub(i)
    link(j) = link(isub)
    link(isub) = j
    do ii = istrt, istop
     isub = nzsub(i)
     lnz(ii) = (lnz(ii) - temp(isub))*diagj
     temp(isub) = 0.0e0
     i = i + 1
    end do
   end if
  end do

 end subroutine

 subroutine spainverse(neqns, xlnz, lnz, xnzsub, nzsub, diag)
  integer, intent(in) :: neqns
  integer, intent(in) :: xnzsub(:), nzsub(:), xlnz(:)
  real(kind=wp), intent(inout) :: diag(:), lnz(:)

  integer :: i, j, k, j1, k1, i_, j1_
  real(kind=wp), allocatable :: lnz_(:), lnzd(:)

  allocate (lnz_(neqns), source=0._wp)
  allocate (lnzd(neqns), source=0._wp)

  do i = neqns, 1, -1
   if (diag(i) < tol) cycle

   i_ = xnzsub(i) - xlnz(i)

   !Uncompression of lnz to the dense vector lnz_
   do j = xlnz(i), xlnz(i + 1) - 1
    lnz(j) = -lnz(j)/diag(i)
   end do

   do j = xlnz(i), xlnz(i + 1) - 1
    j1 = nzsub(i_ + j)
    lnz_(j1) = lnz(j)
   end do

   !Multiplication of the dense vector by the remaining sparse vector and store it in lnzd
   diag(i) = 1._wp/(diag(i)**2)

   do j = xlnz(i), xlnz(i + 1) - 1
    j1 = nzsub(i_ + j)
    lnzd(j1) = lnz(j)*diag(j1)
   end do

   do j = xlnz(i), xlnz(i + 1) - 1
    j1 = nzsub(i_ + j)
    j1_ = xnzsub(j1) - xlnz(j1)
    do k = xlnz(j1), xlnz(j1 + 1) - 1
     k1 = nzsub(j1_ + k)
     lnzd(j1) = lnzd(j1) + lnz(k)*lnz_(k1)
     lnzd(k1) = lnzd(k1) + lnz(k)*lnz_(j1)
    end do
   end do

   !Store inverse from the dense vector in diag and lnz and reset tmp dense arrays
   do j = xlnz(i), xlnz(i + 1) - 1
    j1 = nzsub(i_ + j)
    diag(i) = diag(i) + lnz(j)*lnzd(j1)
    lnz(j) = lnzd(j1)
    lnz_(j1) = 0._wp
    lnzd(j1) = 0._wp
   end do
  end do

 end subroutine

end module

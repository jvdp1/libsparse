submodule (modsparse) modsparse_gen
 !$ use omp_lib
 implicit none

 real(wp),parameter:: tolerance = 1.e-6

contains

!DESTROY
!> @brief Subroutine to reset/destroy a generic object
module subroutine destroy_gen_gen(sparse)
 class(gen_sparse),intent(inout)::sparse

 sparse%namemat='UNKNOWN'
 sparse%dim1=-1
 sparse%dim2=-1
 sparse%unlog=6
 sparse%lsymmetric=.false.
 sparse%lupperstorage=.false.
 if(allocated(sparse%perm))deallocate(sparse%perm)

end subroutine

!**CONJUGATE GRADIENT
module subroutine cg_gen(sparse,x,y,maxiter,tol)
 !sparse*x=y
 class(gen_sparse),intent(in)::sparse
 integer(kind=int32),intent(inout),optional::maxiter
 real(kind=wp),intent(inout)::x(:)
 real(kind=wp),intent(in)::y(:)
 real(kind=wp),intent(inout),optional::tol

 integer(kind=int32)::i,maxiter_
 real(kind=wp)::r(size(x))
 real(kind=wp)::p(size(x))
 real(kind=wp)::Ap(size(x))
 real(kind=wp)::rsnew,rsold,tol_,alpha
 real(kind=wp)::ynorm

 if(.not.sparse%issquare().or..not.sparse%lsymmetric&
    .or.size(x).ne.size(y)&
    .or.size(x).ne.sparse%getdim(2)&
    .or.size(y).ne.sparse%getdim(1)&
     )then
  write(sparse%unlog,'(a)')' ERROR: one of multiple arguments are not conform'
  error stop
 endif

 maxiter_ = min(1000,size(x)-1)
 if(present(maxiter)) maxiter_ = min(maxiter,size(x)-1)

 tol_ = tolerance
 if(present(tol))tol_=tol
 ynorm = norm2(y)
 tol_ = tol_ * ynorm

 Ap=0._wp
 call sparse%mult(1._wp,'n',x,0._wp,Ap)
 r = y - Ap
 p = r
 rsold = sum(r**2)
 
 do i=1, maxiter_
  call sparse%mult(1._wp,'n',p,0._wp,Ap)
  alpha = rsold / dot_product(p,Ap)
  x = x + alpha * p
  r = r - alpha * Ap
  rsnew = sum(r**2)
  if(sqrt(rsnew) < tol_)exit
  p = r +(rsnew / rsold) * p
  rsold = rsnew
 enddo
 
 if(present(maxiter)) maxiter = i
 if(present(tol)) tol = sqrt(rsnew) / ynorm

end subroutine

!**GET ELEMENTS
module function getdim_gen(sparse,dim1) result(dimget)
 class(gen_sparse),intent(in)::sparse
 integer(kind=int32),intent(in)::dim1
 integer(kind=int32)::dimget

 select case(dim1)
  case(1)
   dimget=sparse%dim1
  case(2)
   dimget=sparse%dim2
  case default
   dimget=-1
   write(sparse%unlog,'(a)')' Warning: a sparse matrix has only 2 dimensions!'
 end select

end function

!GET MEMORY
module function getmem_gen(sparse) result(getmem)
 class(gen_sparse),intent(in)::sparse
 integer(kind=int64)::getmem

 getmem=sizeof(sparse%unlog)+sizeof(sparse%dim1)+sizeof(sparse%dim2)+sizeof(sparse%namemat)&
        +sizeof(sparse%lsymmetric)+sizeof(sparse%lupperstorage)
 if(allocated(sparse%perm))getmem=getmem+sizeof(sparse%perm)

end function

!INITIATE GEN SPARSE
module subroutine init_gen(sparse,namemat,dim1,dim2)
 class(gen_sparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::dim1,dim2
 character(len=*),intent(in)::namemat

 sparse%namemat=namemat

 sparse%dim1=dim1
 sparse%dim2=dim2

 sparse%lsorted=.false.
 sparse%lsymmetric=.false.
 sparse%lupperstorage=.false.

end subroutine

!**PRINT
module subroutine print_dim_gen(sparse)
 class(gen_sparse),intent(in)::sparse

 integer(kind=int64)::nel

 write(sparse%unlog,'(/" Type of the matrix           : ",a)')trim(sparse%namemat)
 write(sparse%unlog,'( "  Output unit                 : ",i0)')sparse%unlog
 write(sparse%unlog,'( "  Dimension of the matrix     : ",i0," x ",i0)')sparse%dim1,sparse%dim2
 write(sparse%unlog,'( "  Number of non-zero elements : ",i0)')sparse%nonzero()
 write(sparse%unlog,'( "  Sorted                      : ",l1)')sparse%issorted()
 write(sparse%unlog,'( "  Symmetrix                   : ",l1)')sparse%lsymmetric
 write(sparse%unlog,'( "  Upper storage               : ",l1)')sparse%lupperstorage
 write(sparse%unlog,'( "  Permutation array provided  : ",l1)')allocated(sparse%perm)

 select type(sparse)
  type is(coosparse)
   write(sparse%unlog,'( "  Memory (B)                  : ",i0)')sparse%getmem()
   write(sparse%unlog,'( "  Size of the array           : ",i0)')sparse%nel
  type is(crssparse)
   write(sparse%unlog,'( "  Memory (B)                  : ",i0)')sparse%getmem()
#if (_PARDISO==1)
   write(sparse%unlog,'( "  PARDISO status              : ",i0)')sparse%lpardisofirst
#endif
  class default
 end select
 write(sparse%unlog,'(a)')' '

end subroutine

module subroutine printtofile_gen(sparse,namefile,lint)
 class(gen_sparse),intent(in)::sparse
 character(len=*),intent(in)::namefile
 logical,intent(in),optional::lint

 integer(kind=int32)::un
 logical::linternal

 linternal=.true.
 if(present(lint))linternal=lint

 open(newunit=un,file=namefile,status='replace',action='write')
 call sparse%print(lint=linternal,output=un)
 close(un)

end subroutine

module subroutine printsquaretofile_gen(sparse,namefile)
 class(gen_sparse),intent(inout)::sparse
 character(len=*),intent(in)::namefile

 integer(kind=int32)::un
 logical::linternal

 open(newunit=un,file=namefile,status='replace',action='write')
 call sparse%printsquare(output=un)
 close(un)

end subroutine

!**SET OUTPUT UNIT
module subroutine setoutputunit(sparse,unlog)
 class(gen_sparse),intent(inout)::sparse
 integer(kind=int32)::unlog

 sparse%unlog=unlog

end subroutine

!** SET PERMUTATION VECTOR
module subroutine setpermutation(sparse,array)
 class(gen_sparse),intent(inout)::sparse
 integer(kind=int32)::array(:)

 if(size(array).ne.sparse%getdim(1))then
  write(sparse%unlog,'(a)')' ERROR: The permutation array has a wrong size.'
  error stop
 endif

 !Probably pointer would be better???
 if(.not.allocated(sparse%perm))allocate(sparse%perm(sparse%getdim(1)))
 sparse%perm=array

end subroutine

! SET THE STATUS SORTED
module subroutine setsorted(sparse,ll)
 class(gen_sparse),intent(inout)::sparse
 logical,intent(in)::ll

 sparse%lsorted=ll

end subroutine

! SET THE STATUS SYMMETRIC
module subroutine setsymmetric(sparse,ll)
 class(gen_sparse),intent(inout)::sparse
 logical,intent(in),optional::ll

 logical::lll

 if(.not.sparse%issquare().and..not.present(ll))then
  write(sparse%unlog,'(a)')' ERROR: the sparse matrix is not square and cannot be set to symmetric!'
  error stop
 elseif(.not.sparse%issquare().and.present(ll))then
  if(ll)then
   write(sparse%unlog,'(a)')' ERROR: the sparse matrix is not square and cannot be set to symmetric as requested!'
   error stop
  endif
 endif

 lll=.true.
 if(present(ll))lll=ll

 sparse%lsymmetric=lll

end subroutine


!**OTHER
module function issorted(sparse) result(ll)
 class(gen_sparse),intent(in)::sparse
 logical::ll

 ll=sparse%lsorted

end function

module function issquare(sparse) result(ll)
 class(gen_sparse),intent(in)::sparse
 logical::ll

 ll=.true.
 if(sparse%dim1.ne.sparse%dim2)ll=.false.

end function

!CHECKS
module function validvalue_gen(sparse,row,col) result(lvalid)
 class(gen_sparse),intent(in)::sparse
 integer(kind=int32),intent(in)::row,col
 logical::lvalid

 lvalid=.true.
 if((row.lt.1.or.row.gt.sparse%dim1).or.(col.lt.1.or.col.gt.sparse%dim2))lvalid=.false.

end function

module function validnonzero_gen(sparse,val) result(lvalid)
 class(gen_sparse),intent(in)::sparse
 real(kind=wp),intent(in)::val
 logical::lvalid

 lvalid=.true.
 if((abs(val)<epsilon(val)))lvalid=.false.

end function

module function uppervalue_gen(row,col) result(lvalid)
 integer(kind=int32),intent(in)::row,col
 logical::lvalid

 lvalid=.true.
 if(row.gt.col)lvalid=.false.

end function

end submodule

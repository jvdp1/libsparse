submodule (modsparse) modsparse_coo
#if (_DP==0)
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real32
#else
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real64
#endif
 use modhash, only:hashf,roundinguppower2
 !$ use omp_lib
 implicit none

contains

!**CONSTRUCTOR
module function constructor_coo(m,n,nel,lupper,unlog) result(sparse)
 type(coosparse)::sparse
 integer(kind=int32),intent(in)::m
 integer(kind=int32),intent(in),optional::n,unlog
 integer(kind=int64),intent(in),optional::nel
 logical,intent(in),optional::lupper

 call sparse%initialize('COO',m,m)

 if(present(n))sparse%dim2=n
 if(present(lupper))sparse%lupperstorage=lupper
 if(present(unlog))sparse%unlog=unlog

 sparse%lsymmetric=.false.

 sparse%filled=0_int64

 sparse%nel=roundinguppower2(100_int64)
 if(present(nel))sparse%nel=roundinguppower2(int(nel,int64))
 allocate(sparse%ij(2,sparse%nel),sparse%a(sparse%nel))
 sparse%ij=0
 sparse%a=0._wp

end function

module subroutine constructor_sub_coo(sparse,m,n,nel,lupper,unlog)
 class(coosparse),intent(out)::sparse
 integer(kind=int32),intent(in)::m
 integer(kind=int32),intent(in),optional::n,unlog
 integer(kind=int64),intent(in),optional::nel
 logical,intent(in),optional::lupper

 call sparse%initialize('COO',m,m)

 if(present(n))sparse%dim2=n
 if(present(lupper))sparse%lupperstorage=lupper
 if(present(unlog))sparse%unlog=unlog

 sparse%lsymmetric=.false.

 sparse%filled=0_int64

 sparse%nel=roundinguppower2(100_int64)
 if(present(nel))sparse%nel=roundinguppower2(int(nel,int64))
 allocate(sparse%ij(2,sparse%nel))
 sparse%ij=0
 allocate(sparse%a(sparse%nel))
 sparse%a=0._wp

end subroutine

!**DESTROY
module subroutine destroy_coo(sparse)
 class(coosparse),intent(inout)::sparse

 call sparse%destroy_gen_gen()

 sparse%nel=-1_int64
 sparse%filled=-1_int64
 if(allocated(sparse%ij))deallocate(sparse%ij)
 if(allocated(sparse%a))deallocate(sparse%a)

end subroutine

!**DIAGONAL ELEMENTS
module function diag_vect_coo(sparse) result(array)
 class(coosparse),intent(inout)::sparse
 real(kind=wp),allocatable::array(:)

 integer(kind=int32)::ndiag,i

 ndiag=min(sparse%dim1,sparse%dim2)

 allocate(array(ndiag))
 array=0.0_wp

 do i=1,ndiag
  array(i)=sparse%get(i,i)
 enddo

end function

!**ADD ELEMENTS
module recursive subroutine add_coo(sparse,row,col,val)
 class(coosparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::row,col
 real(kind=wp),intent(in)::val

 integer(kind=int64)::hash,i8
 real(kind=real32),parameter::maxratiofilled=maxratiofilled_par
 real(kind=real32)::ratiofilled
 type(coosparse),allocatable::sptmp

 if(.not.validvalue_gen(sparse,row,col))return
 if(.not.validnonzero_gen(sparse,val))return
 if(sparse%lupperstorage.and..not.uppervalue_gen(row,col))return

 hash=hashf(row,col,sparse%ij,sparse%nel,sparse%filled,.false.)
 ratiofilled=real(sparse%filled)/real(sparse%nel)

 if(hash.eq.-1.or.ratiofilled.gt.maxratiofilled)then
  !matrix probably full, or nothing available within the n requested searches
  !1. Copy matrix
  !sptmp=coosparse(sparse%dim1,sparse%dim2,sparse%nel*2)   !to avoid the copy through a temporary array
  allocate(sptmp);call sptmp%init(sparse%dim1,sparse%dim2,sparse%nel*2)
  do i8=1_int64,sparse%nel
   call sptmp%add(sparse%ij(1,i8),sparse%ij(2,i8),sparse%a(i8))
  enddo
  !2. reallocate matrix using move_alloc
#if(_VERBOSE>0)
  write(sparse%unlog,'(2(a,i0))')'  Current | New size COO: ',sparse%nel,' | ',sptmp%nel
#endif
  sparse%nel=sptmp%nel
  sparse%filled=sptmp%filled
  if(allocated(sparse%ij))deallocate(sparse%ij)
  if(allocated(sparse%a))deallocate(sparse%a)
  call move_alloc(sptmp%ij,sparse%ij)
  call move_alloc(sptmp%a,sparse%a)
  call sptmp%destroy()
  !3. Search for a new address in the new matrix
  hash=hashf(row,col,sparse%ij,sparse%nel,sparse%filled,.false.)
  ratiofilled=real(sparse%filled)/real(sparse%nel)
 endif

 if(hash.gt.0_int64)then!.and.ratiofilled.le.maxratiofilled)then
  sparse%a(hash)=sparse%a(hash)+val
 else
  !is it possible?
  write(sparse%unlog,*)' ERROR: unexpected'!,__FILE__,__LINE__
  stop
 endif

end subroutine

!**GET ELEMENTS
module function get_coo(sparse,row,col) result(val)
 class(coosparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::row,col
 real(kind=wp)::val

 integer(kind=int32)::trow,tcol
 integer(kind=int64)::hash

 val=0.0_wp

 trow=row
 tcol=col
 if(sparse%lupperstorage.and.sparse%lsymmetric.and.row.gt.col)then
  !swap row-col
  trow=col
  tcol=row
 endif

 hash=hashf(trow,tcol,sparse%ij,sparse%nel,sparse%filled,.true.)

 if(hash.gt.0_int64)val=sparse%a(hash)

end function

!** GET MEMORY
module function getmem_coo(sparse) result(getmem)
 class(coosparse),intent(in)::sparse
 integer(kind=int64)::getmem

 getmem=sparse%getmem_gen()+sizeof(sparse%nel)+sizeof(sparse%filled)
 if(allocated(sparse%ij))getmem=getmem+sizeof(sparse%ij)
 if(allocated(sparse%a))getmem=getmem+sizeof(sparse%a)

end function

!**EXTERNAL



end submodule

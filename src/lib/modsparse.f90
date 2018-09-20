!> Module containing three types of sparse matrices:  
!> * Linked list (llsparse)  
!> * COOrdinate storage (coosparse)  
!> * Compressed Row Storage (crssparse)  
!> @todo Use of submodules (one submodule for each type of matrix)

module modsparse
 use modkind
 use modhash
 !$ use omp_lib
 implicit none
 private
 public::llsparse,coosparse,crssparse
 public::assignment(=)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
 !> @brief Generic object containing dimensions, storage format, and output unit
 type,abstract::gen_sparse
  private
  integer(kind=int4)::unlog=6
  integer(kind=int4)::dim1,dim2
  character(len=15)::namemat='UNKNOWN'
  logical::lupperstorage
  contains
  private
  procedure(destroy_gen),public,deferred::destroy
  procedure(get_gen),public,deferred::get
  procedure(nonzero_gen),public,deferred::nonzero
  procedure(print_gen),public,deferred::print
  procedure(printsquare_gen),public,deferred::printsquare

  !> @brief Returns the dimension of the matrix; e.g., ...=mat%getdim(1)
  procedure,public::getdim=>getdim_gen
  !> @brief Returns true if square matrix; else returns false
  procedure,public::lsquare
  !> @brief Returns the number of non-zero elements
  procedure,public::printtofile=>printtofile_gen
  !> @brief Prints the sparse matrix in a rectangular/square format to the output mat\%unlog
  procedure,public::printsquaretofile=>printsquaretofile_gen
  procedure,public::printstats=>print_dim_gen
  procedure::destroy_gen_gen
 end type
 
 abstract interface
  subroutine destroy_gen(sparse)
   import::gen_sparse
   class(gen_sparse),intent(inout)::sparse
  end subroutine
  function get_gen(sparse,row,col) result(val)
   import::int4,real8,gen_sparse
   class(gen_sparse),intent(inout)::sparse
   integer(kind=int4),intent(in)::row,col
   real(kind=real8)::val
  end function
  function nonzero_gen(sparse) result(nel)
   import::int8,gen_sparse
   class(gen_sparse),intent(in)::sparse
   integer(kind=int8)::nel
  end function
  subroutine print_gen(sparse,lint,output)
   import::int4,gen_sparse
   class(gen_sparse),intent(in)::sparse
   integer(kind=int4),intent(in),optional::output
   logical,intent(in),optional::lint
  end subroutine
  subroutine printsquare_gen(sparse,output)
   import::int4,gen_sparse
   class(gen_sparse),intent(inout)::sparse
   integer(kind=int4),intent(in),optional::output
   end subroutine
 end interface


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!COO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
 !> @brief Object for COOrdinate storage
 type,extends(gen_sparse)::coosparse
  private
  integer(kind=int4),allocatable::ij(:,:)
  integer(kind=int8)::nel
  integer(kind=int8)::filled
  real(kind=real8),allocatable::a(:)
  contains
  private
  !> @brief Adds the value val to mat(row,col); e.g., call mat\%add(row,col,val)
  procedure,public::add=>add_coo 
  !> @brief Deallocates the sparse matrix and sets to default values 
  procedure,public::destroy=>destroy_scal_coo
  !> @brief Returns the value of mat(row,col); e.g., ...=mat\%get(row,col)
  procedure,public::get=>get_coo
  !> @brief Returns the number of non-zero elements
  procedure,public::nonzero=>totalnumberofelements_coo
  !> @brief Prints the sparse matrix to the output mat\%unlog
  procedure,public::print=>print_coo
  !> @brief Prints the sparse matrix in a rectangular/square format to the default output
  procedure,public::printsquare=>printsquare_coo
  !> @brief Sets an entry to a certain value (even if equal to 0); e.g., call mat\%set(row,col,val)
  procedure,public::set=>set_coo
  final::deallocate_scal_coo
 end type

 !> @brief Constructor; e.g., mat=coosparse(dim1,[dim2],[#elements],[upper_storage],[output_unit])
 interface coosparse
  module procedure constructor_coo
 end interface
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CRS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
 !> @brief Object for Compressed Row Storage
 type,extends(gen_sparse)::crssparse
  private
  integer(kind=int4),allocatable::ia(:)
  integer(kind=int4),allocatable::ja(:)
  real(kind=real8),allocatable::a(:)
  contains
  private
  !> @brief Adds the value val to mat(row,col); e.g., call mat\%add(row,col,val)
  procedure,public::add=>add_crs
  !> @brief Deallocates the sparse matrix and sets to default values 
  procedure,public::destroy=>destroy_scal_crs
  !> @brief Returns the value of mat(row,col); e.g., ...=mat\%get(row,col)
  procedure,public::get=>get_crs
  !> @brief Returns the number of non-zero elements
  procedure,public::nonzero=>totalnumberofelements_crs
  !> @brief Prints the sparse matrix to the output sparse\%unlog
  procedure,public::print=>print_crs
  !> @brief Prints the sparse matrix in a rectangular/square format to the default output
  procedure,public::printsquare=>printsquare_crs
  !> @brief Sets an entry to a certain value (even if equal to 0); condition: the entry must exist; e.g., call mat\%set(row,col,val)
  procedure,public::set=>set_crs
  !> @brief Sorts the elements in a ascending order within a row
  procedure,public::sort=>sort_crs
  final::deallocate_scal_crs
 end type

 !> @brief Constructor; e.g., mat=crsparse(dim1,#elements,[dim2],[upper_storage],[output_unit])
 interface crssparse
  module procedure constructor_crs
 end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!LINKED LIST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
 type::ptrnode
  type(node),pointer::p=>null()
  contains
  procedure,public::size=>totalnumberofelements_ptrnode
 end type 

 !> @brief Object for Linked List
 type,extends(gen_sparse)::llsparse
  private
  type(ptrnode),allocatable::heads(:)
  contains
  private
  !> @brief Adds the value val to mat(row,col); e.g., call mat\%add(row,col,val)
  procedure,public::add=>addinorder_ll
  procedure,public::addfast =>addtohead_ll
  procedure,public::addtohead =>addtohead_ll
  procedure,public::addtotail =>addtotail_ll
  !> @brief Returns the value of mat(row,col); e.g., ...=mat\%get(row,col)
  procedure,public::get=>get_ll
  !> @brief Returns the number of non-zero elements
  procedure,public::nonzero=>totalnumberofelements_ll
  !> @brief Prints the sparse matrix to the output sparse\%unlog
  procedure,public::print=>print_ll
  !> @brief Prints the sparse matrix in a rectangular/square format to the default output
  procedure,public::printsquare=>printsquare_ll
  !> @brief Deallocates the sparse matrix and sets to default values 
  procedure,public::destroy=>destroy_ll
  final::deallocate_scal_ll
 end type

 type::node
  integer(kind=int4)::col
  real(kind=real8)::val
  type(ptrnode)::next
  contains
  private
  procedure::equal_node
  generic,public::assignment(=)=>equal_node
 end type

 !> @brief Constructor; e.g., mat=llsparse(dim1,[dim2],[upper_storage],[output_unit])
 interface llsparse
  module procedure constructor_ll
 end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GENERAL INTERFACES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
 !> @brief Converts sparse matrices from one format to another one; e.g., crsmat=coomat
 interface assignment(=)
  module procedure convertfromlltocoo,convertfromlltocrs,convertfromcootocrs
 end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
!DESTROY
subroutine destroy_gen_gen(sparse)
 class(gen_sparse),intent(inout)::sparse

 sparse%namemat='UNKNOWN'
 sparse%dim1=-1
 sparse%dim2=-1
 sparse%unlog=6
 sparse%lupperstorage=.false.

end subroutine

!**GET ELEMENTS
function getdim_gen(sparse,dim1) result(dimget)
 class(gen_sparse),intent(in)::sparse
 integer(kind=int4),intent(in)::dim1
 integer(kind=int4)::dimget

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

!**PRINT
subroutine print_dim_gen(sparse)
 class(gen_sparse),intent(in)::sparse

 integer(kind=int8)::nel

 write(sparse%unlog,'(/" Type of the matrix           : ",a)')trim(sparse%namemat)
 write(sparse%unlog,'( "  Dimension of the matrix     : ",i0," x ",i0)')sparse%dim1,sparse%dim2
 write(sparse%unlog,'( "  Upper storage               : ",l)')sparse%lupperstorage
 write(sparse%unlog,'( "  Number of non-zero elements : ",i0)')sparse%nonzero()
 
 select type(sparse)
  type is(coosparse)
   write(sparse%unlog,'( "  Size of the array           : ",i0)')sparse%nel
  class default
 end select

end subroutine

subroutine printtofile_gen(sparse,namefile,lint)
 class(gen_sparse),intent(in)::sparse
 character(len=*),intent(in)::namefile
 logical,intent(in),optional::lint

 integer(kind=int4)::un
 logical::linternal

 linternal=.true.
 if(present(lint))linternal=lint

 open(newunit=un,file=namefile,status='replace',action='write')
 call sparse%print(lint=linternal,output=un)
 close(un)

end subroutine

subroutine printsquaretofile_gen(sparse,namefile)
 class(gen_sparse),intent(inout)::sparse
 character(len=*),intent(in)::namefile

 integer(kind=int4)::un
 logical::linternal

 open(newunit=un,file=namefile,status='replace',action='write')
 call sparse%printsquare(output=un)
 close(un)

end subroutine

!**OTHER
function lsquare(sparse) result(ll)
 class(gen_sparse),intent(in)::sparse
 logical::ll

 ll=.true.
 if(sparse%dim1.ne.sparse%dim2)ll=.false.

end function



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!COO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
!**CONSTRUCTOR
function constructor_coo(m,n,nel,lupper,unlog) result(sparse)
 type(coosparse)::sparse
 integer(kind=int4),intent(in)::m
 integer(kind=int4),intent(in),optional::n,unlog
 integer(kind=int8),intent(in),optional::nel
 logical,intent(in),optional::lupper

 sparse%namemat='COO'
 sparse%dim1=m
 sparse%dim2=m
 if(present(n))sparse%dim2=n

 sparse%filled=0_int8

 sparse%nel=roundinguppower2(100_int8)
 if(present(nel))sparse%nel=roundinguppower2(int(nel,int8))
 allocate(sparse%ij(2,sparse%nel),sparse%a(sparse%nel))
 sparse%ij=0
 sparse%a=0._real8

 sparse%lupperstorage=.false.
 if(present(lupper))sparse%lupperstorage=lupper

 if(present(unlog))sparse%unlog=unlog

end function

!**DESTROY
subroutine destroy_scal_coo(sparse)
 class(coosparse),intent(inout)::sparse

 call sparse%destroy_gen_gen()

 sparse%nel=-1_int8
 sparse%filled=-1_int8
 if(allocated(sparse%ij))deallocate(sparse%ij)
 if(allocated(sparse%a))deallocate(sparse%a)

end subroutine

!**ADD ELEMENTS
recursive subroutine add_coo(sparse,row,col,val)
 class(coosparse),intent(inout)::sparse
 integer(kind=int4),intent(in)::row,col
 real(kind=int8),intent(in)::val

 integer(kind=int8)::hash,i8
 real(kind=real4),parameter::maxratiofilled=0.80
 real(kind=real4)::ratiofilled
 type(coosparse)::sptmp
 
 if(.not.validvalue_gen(sparse,row,col))return
 if(.not.validnonzero_gen(sparse,val))return
 if(sparse%lupperstorage.and..not.uppervalue_gen(row,col))return

 hash=hashf(row,col,sparse%ij,sparse%nel,sparse%filled,.false.)
 ratiofilled=real(sparse%filled)/real(sparse%nel)

 if(hash.eq.-1.or.ratiofilled.gt.maxratiofilled)then
  !matrix probably full, or nothing available within the n requested searches
  !1. Copy matrix
  sptmp=coosparse(sparse%dim1,sparse%dim2,sparse%nel*2)
  do i8=1_int8,sparse%nel
   call sptmp%add(sparse%ij(1,i8),sparse%ij(2,i8),sparse%a(i8))
  enddo
  !2. reallocate matrix using move_alloc
  write(sparse%unlog,'(2(a,i0))')'  Current | New size COO: ',sparse%nel,' | ',sptmp%nel
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

 if(hash.gt.0_int8)then!.and.ratiofilled.le.maxratiofilled)then
  sparse%a(hash)=sparse%a(hash)+val
 else
  !is it possible?
  write(sparse%unlog,*)' ERROR: unexpected'!,__FILE__,__LINE__
  stop
 endif
 
end subroutine

!**GET ELEMENTS
function get_coo(sparse,row,col) result(val)
 class(coosparse),intent(inout)::sparse
 integer(kind=int4),intent(in)::row,col
 real(kind=real8)::val
 
 integer(kind=int4)::trow,tcol
 integer(kind=int8)::hash

 val=0_real8
 
 trow=row
 tcol=col
 if(sparse%lupperstorage.and.row.gt.col)then
  !swap row-col
  trow=col
  tcol=row
 endif

 hash=hashf(trow,tcol,sparse%ij,sparse%nel,sparse%filled,.true.)
 
 if(hash.gt.0)val=sparse%a(hash)

end function

!**NUMBER OF ELEMENTS
function totalnumberofelements_coo(sparse) result(nel)
 class(coosparse),intent(in)::sparse
 integer(kind=int8)::nel

 nel=sparse%filled

end function

!**PRINT
subroutine print_coo(sparse,lint,output)
 class(coosparse),intent(in)::sparse
 integer(kind=int4),intent(in),optional::output
 logical,intent(in),optional::lint

 integer(kind=int4)::un,row,col
 integer(kind=int8)::i8
 real(kind=real8)::val
 logical::linternal

 linternal=.true.
 if(present(lint))linternal=lint

 un=sparse%unlog
 if(present(output))un=output

 do i8=1,sparse%nel
  row=sparse%ij(1,i8)
  col=sparse%ij(2,i8)
  if(row.eq.0.and.col.eq.0)cycle
  val=sparse%a(i8)
  !if(.not.validvalue_gen(sparse,row,col))cycle  !it should never happen
  !if(.not.validnonzero_gen(sparse,val))cycle    !to print as internal
  write(un,'(2(i0,x),g0)')row,col,val
  if(.not.linternal.and.sparse%lupperstorage.and.row.ne.col)then
   write(un,'(2(i0,x),g0)')col,row,val
  endif
 enddo

end subroutine

subroutine printsquare_coo(sparse,output)
 class(coosparse),intent(inout)::sparse
 integer(kind=int4),intent(in),optional::output

 integer(kind=int4)::i,j,un
 real(kind=real8)::val
 real(kind=real8),allocatable::tmp(:)

 un=sparse%unlog
 if(present(output))un=output

 allocate(tmp(sparse%dim2))

 do i=1,sparse%dim1
  tmp=0.d0
  do j=1,sparse%dim2
   tmp(j)=sparse%get(i,j)
  enddo
  write(un,'(10000(f9.3,x))')tmp
 enddo

 deallocate(tmp)

end subroutine

!**SET ELEMENTS
recursive subroutine set_coo(sparse,row,col,val)
 !from add_coo
 class(coosparse),intent(inout)::sparse
 integer(kind=int4),intent(in)::row,col
 real(kind=int8),intent(in)::val

 integer(kind=int8)::hash,i8
 real(kind=real4),parameter::maxratiofilled=0.80
 real(kind=real4)::ratiofilled
 type(coosparse)::sptmp
 
 if(.not.validvalue_gen(sparse,row,col))return
 !if(.not.validnonzero_gen(sparse,val))return
 if(sparse%lupperstorage.and..not.uppervalue_gen(row,col))return

 hash=hashf(row,col,sparse%ij,sparse%nel,sparse%filled,.false.)
 ratiofilled=real(sparse%filled)/real(sparse%nel)

 if(hash.eq.-1.or.ratiofilled.gt.maxratiofilled)then
  !matrix probably full, or nothing available within the n requested searches
  !1. Copy matrix
  sptmp=coosparse(sparse%dim1,sparse%dim2,sparse%nel*2)
  do i8=1_int8,sparse%nel
   call sptmp%add(sparse%ij(1,i8),sparse%ij(2,i8),sparse%a(i8))
  enddo
  !2. reallocate matrix using move_alloc
  write(sparse%unlog,'(2(a,i0))')'  Current | New size COO: ',sparse%nel,' | ',sptmp%nel
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

 if(hash.gt.0_int8)then!.and.ratiofilled.le.maxratiofilled)then
  sparse%a(hash)=val
 else
  !is it possible?
  write(sparse%unlog,*)' ERROR: unexpected'!,__FILE__,__LINE__
  stop
 endif
 
end subroutine

!FINAL
subroutine deallocate_scal_coo(sparse)
 type(coosparse),intent(inout)::sparse

 call destroy_scal_coo(sparse)

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CRS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
!**CONSTRUCTOR
function constructor_crs(m,nel,n,lupper,unlog) result(sparse)
 type(crssparse)::sparse
 integer(kind=int4),intent(in)::m
 integer(kind=int4),intent(in)::nel
 integer(kind=int4),intent(in),optional::n,unlog
 logical,intent(in),optional::lupper

 sparse%namemat='CRS'
 sparse%dim1=m
 sparse%dim2=m
 if(present(n))sparse%dim2=n
 
 allocate(sparse%ia(sparse%dim1+1),sparse%ja(nel),sparse%a(nel))
 sparse%ia=0
 sparse%ia(sparse%dim1+1)=-nel
 sparse%ja=0
 sparse%a=0._real8

 sparse%lupperstorage=.false.
 if(present(lupper))sparse%lupperstorage=lupper

 if(present(unlog))sparse%unlog=unlog

end function

!**DESTROY
subroutine destroy_scal_crs(sparse)
 class(crssparse),intent(inout)::sparse

 print*,'destroy crs'

 call sparse%destroy_gen_gen()

 if(allocated(sparse%ia))deallocate(sparse%ia)
 if(allocated(sparse%ja))deallocate(sparse%ja)
 if(allocated(sparse%a))deallocate(sparse%a)

end subroutine

!**ADD ELEMENTS
subroutine add_crs(sparse,row,col,val,error)
 !add a value only to an existing one
 class(crssparse),intent(inout)::sparse
 integer(kind=int4),intent(in)::row,col
 integer(kind=int4),intent(out),optional::error
 real(kind=int8),intent(in)::val

 integer(kind=int4)::i
 integer(kind=int4)::ierror    !added: error=0;Not existing: error=-1;matrix not inited: error=-10
 
 if(present(error))error=0

 if(.not.validvalue_gen(sparse,row,col))return
 if(.not.validnonzero_gen(sparse,val))return
 if(sparse%lupperstorage.and..not.uppervalue_gen(row,col))return

 if(sparse%ia(row).eq.0)then
  if(present(error))error=-10
  return
 endif

 do i=sparse%ia(row),sparse%ia(row+1)-1
  if(sparse%ja(i).eq.col)then
   sparse%a(i)=sparse%a(i)+val
   if(present(error))error=0
   return
  endif
 enddo
 
 if(present(error))error=-1

end subroutine

!**GET ELEMENTS
function get_crs(sparse,row,col) result(val)
 class(crssparse),intent(inout)::sparse
 integer(kind=int4),intent(in)::row,col
 real(kind=real8)::val
 
 integer(kind=int4)::i,trow,tcol

 val=0_real8
 
 trow=row
 tcol=col
 if(sparse%lupperstorage.and.row.gt.col)then
  !swap row-col
  trow=col
  tcol=row
 endif

 do i=sparse%ia(trow),sparse%ia(trow+1)-1
  if(sparse%ja(i).eq.tcol)then
   val=sparse%a(i)
   exit
  endif
 enddo

end function

!**NUMBER OF ELEMENTS
function totalnumberofelements_crs(sparse) result(nel)
 class(crssparse),intent(in)::sparse
 integer(kind=int8)::nel

 nel=int(sparse%ia(sparse%dim1+1),int8)-1_int8

end function

!**PRINT
subroutine print_crs(sparse,lint,output)
 class(crssparse),intent(in)::sparse
 integer(kind=int4),intent(in),optional::output
 logical,intent(in),optional::lint

 integer(kind=int4)::i
 integer(kind=int4)::un,j
 logical::linternal

 linternal=.true.
 if(present(lint))linternal=lint

 un=sparse%unlog
 if(present(output))un=output

 do i=1,sparse%dim1
  do j=sparse%ia(i),sparse%ia(i+1)-1
   !write(un,'(2(i0,x),g0)')i,sparse%ja(j),sparse%a(j)
   write(un,'(2i8,x,f0.4)')i,sparse%ja(j),sparse%a(j)
   if(.not.linternal.and.sparse%lupperstorage.and.i.ne.sparse%ja(j))then
    write(un,'(2(i0,x),g0)')sparse%ja(j),i,sparse%a(j)
   endif
  enddo
 enddo

end subroutine

subroutine printsquare_crs(sparse,output)
 class(crssparse),intent(inout)::sparse
 integer(kind=int4),intent(in),optional::output

 integer(kind=int4)::i,j,un
 real(kind=real8)::val
 real(kind=real8),allocatable::tmp(:)

 un=sparse%unlog
 if(present(output))un=output

 allocate(tmp(sparse%dim2))

 do i=1,sparse%dim1
  tmp=0.d0
  !could be implemented in a more efficient way
  do j=1,sparse%dim2
   tmp(j)=sparse%get(i,j)
  enddo
  write(un,'(10000(f9.3,x))')tmp
 enddo

 deallocate(tmp)

end subroutine

!**SET ELEMENTS
subroutine set_crs(sparse,row,col,val,error)
 !add a value only to an existing one
 class(crssparse),intent(inout)::sparse
 integer(kind=int4),intent(in)::row,col
 integer(kind=int4),intent(out),optional::error
 real(kind=int8),intent(in)::val

 integer(kind=int4)::i
 integer(kind=int4)::ierror    !added: error=0;Not existing: error=-1;matrix not inited: error=-10
 
 if(present(error))error=0

 if(.not.validvalue_gen(sparse,row,col))return
 !if(.not.validnonzero_gen(sparse,val))return
 if(sparse%lupperstorage.and..not.uppervalue_gen(row,col))return

 if(sparse%ia(row).eq.0)then
  if(present(error))error=-10
  return
 endif

 do i=sparse%ia(row),sparse%ia(row+1)-1
  if(sparse%ja(i).eq.col)then
   sparse%a(i)=val
   if(present(error))error=0
   return
  endif
 enddo
 
 if(present(error))error=-1

end subroutine

!**SORT ARRAY
subroutine sort_crs(sparse)
 ! sort vectors ja and a by increasing order
 class(crssparse),intent(inout)::sparse

 integer(kind=int4)::dir,endd,i,j,k,n,start,stkpnt
 integer(kind=int4)::d1,d2,d3,dmnmx,tmp
 integer(kind=int4)::stack(2,32)
 integer(kind=int4),allocatable::d(:)
 integer(kind=int4),parameter::select=20
 real(kind=real8)::umnmx,tmpu
 real(kind=real8),allocatable::u(:)

 do k=1,sparse%dim1
  n=sparse%ia(k+1)-sparse%ia(k)
  if(n.gt.1)then
   allocate(d(n),u(n))
   !copy of the vector to be sorted
   d=sparse%ja(sparse%ia(k):sparse%ia(k+1)-1)
   u=sparse%a(sparse%ia(k):sparse%ia(k+1)-1)
   !sort the vectors
   !from dlasrt.f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Quick return if possible
   stkpnt = 1
   stack( 1, 1 ) = 1
   stack( 2, 1 ) = n
   10 start = stack( 1, stkpnt )
   endd = stack( 2, stkpnt )
   stkpnt = stkpnt - 1
   IF( endd-start <= select .AND. endd-start > 0 ) THEN
   !Do Insertion sort on D( START:ENDD )
   !Sort into increasing order
     DO i = start + 1, endd
      DO j = i, start + 1, -1
       IF( d( j ) < d( j-1 ) ) THEN
         dmnmx = d( j )
         d( j ) = d( j-1 )
         d( j-1 ) = dmnmx
         umnmx = u( j )
         u( j ) = u( j-1 )
         u( j-1 ) = umnmx
       ELSE
         CYCLE
       END IF
      END DO
     END DO
   ELSE IF( endd-start > select ) THEN
   !Partition D( START:ENDD ) and stack parts, largest one first
   !Choose partition entry as median of 3
    d1 = d( start )
    d2 = d( endd )
    i = ( start + endd ) / 2
    d3 = d( i )
    IF( d1 < d2 ) THEN
      IF( d3 < d1 ) THEN
        dmnmx = d1
      ELSE IF( d3 < d2 ) THEN
        dmnmx = d3
      ELSE
        dmnmx = d2
      END IF
    ELSE
      IF( d3 < d2 ) THEN
        dmnmx = d2
      ELSE IF( d3 < d1 ) THEN
        dmnmx = d3
      ELSE
        dmnmx = d1
      END IF
    END IF
    !Sort into increasing order
     i = start - 1
     j = endd + 1
     90 j = j - 1
     IF( d( j ) > dmnmx ) GO TO 90
     110 i = i + 1
     IF( d( i ) < dmnmx ) GO TO 110
     IF( i < j ) THEN
       tmp = d( i )
       d( i ) = d( j )
       d( j ) = tmp
       tmpu = u( i )
       u( i ) = u( j )
       u( j ) = tmpu
       GO TO 90
     END IF
     IF( j-start > endd-j-1 ) THEN
       stkpnt = stkpnt + 1
       stack( 1, stkpnt ) = start
       stack( 2, stkpnt ) = j
       stkpnt = stkpnt + 1
       stack( 1, stkpnt ) = j + 1
       stack( 2, stkpnt ) = endd
     ELSE
       stkpnt = stkpnt + 1
       stack( 1, stkpnt ) = j + 1
       stack( 2, stkpnt ) = endd
       stkpnt = stkpnt + 1
       stack( 1, stkpnt ) = start
       stack( 2, stkpnt ) = j
     END IF
   END IF
   IF( stkpnt > 0 ) GO TO 10
   !end from dlasrt.f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !copy back the sorted vectors
   sparse%ja(sparse%ia(k):sparse%ia(k+1)-1)=d
   sparse%a(sparse%ia(k):sparse%ia(k+1)-1)=u
   deallocate(d,u)
  endif
 enddo
 
end subroutine

!FINAL
subroutine deallocate_scal_crs(sparse)
 type(crssparse),intent(inout)::sparse

 call destroy_scal_crs(sparse)

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!LINKED LIST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
!**CONSTRUCTOR
function constructor_ll(m,n,lupper,unlog) result(sparse)
 type(llsparse)::sparse
 integer(kind=int4),intent(in)::m
 integer(kind=int4),intent(in),optional::n,unlog
 logical,intent(in),optional::lupper

 sparse%namemat='LINKED LIST'
 sparse%dim1=m
 sparse%dim2=m
 if(present(n))sparse%dim2=n

 allocate(sparse%heads(sparse%dim1))

 sparse%lupperstorage=.false.
 if(present(lupper))sparse%lupperstorage=lupper

 if(present(unlog))sparse%unlog=unlog

end function

!**DESTROY
subroutine destroy_scal_ptrnode(pnode)
 type(ptrnode)::pnode

 type(ptrnode)::cursor
 
 do while(associated(pnode%p))
  cursor=pnode
  pnode=pnode%p%next
  deallocate(cursor%p)
  nullify(cursor%p)
 enddo

end subroutine

subroutine destroy_ll(sparse)
 class(llsparse),intent(inout)::sparse
 integer(kind=int4)::i

 print*,'destroy ll'
 call sparse%destroy_gen_gen()

 if(allocated(sparse%heads))then
  do i=1,size(sparse%heads)
   call destroy_scal_ptrnode(sparse%heads(i))
  enddo
  deallocate(sparse%heads)
 endif

end subroutine

!EQUALITIES
subroutine equal_node(nodeout,nodein)
 class(node),intent(out)::nodeout
 class(node),intent(in)::nodein
 
 nodeout%col=nodein%col
 nodeout%val=nodein%val

end subroutine

!**ADD ELEMENTS
subroutine addtohead_ptrnode(pnode,col,val)
 type(ptrnode),intent(inout),pointer::pnode
 integer(kind=int4),intent(in)::col
 real(kind=real8),intent(in)::val

 type(ptrnode)::cursor

 allocate(cursor%p)
 cursor%p%next=pnode
 cursor%p%col=col
 cursor%p%val=val
 pnode=cursor

end subroutine

subroutine addtohead_ll(sparse,row,col,val)
 class(llsparse),intent(inout),target::sparse
 integer(kind=int4),intent(in)::row,col
 real(kind=real8),intent(in)::val

 type(ptrnode)::cursor

 if(.not.validvalue_gen(sparse,row,col))return
 if(.not.validnonzero_gen(sparse,val))return
 if(sparse%lupperstorage.and..not.uppervalue_gen(row,col))return

 allocate(cursor%p)
 cursor%p%next=sparse%heads(row)
 cursor%p%col=col
 cursor%p%val=val
 sparse%heads(row)=cursor

end subroutine

subroutine addinorder_ll(sparse,row,col,val)
 class(llsparse),intent(inout),target::sparse
 integer(kind=int4),intent(in)::row,col
 real(kind=real8),intent(in)::val

 type(ptrnode),pointer::cursor

 if(.not.validvalue_gen(sparse,row,col))return
 if(.not.validnonzero_gen(sparse,val))return
 if(sparse%lupperstorage.and..not.uppervalue_gen(row,col))return

 cursor=>sparse%heads(row)
 do while(associated(cursor%p))
  if(cursor%p%col.ge.col)then
   if(.not.col.ge.cursor%p%col)then
    call addtohead_ptrnode(cursor,col,val)
   else
    cursor%p%val=cursor%p%val+val
   endif
   return
  endif 
  cursor=>cursor%p%next
 enddo
 allocate(cursor%p)
 cursor%p%col=col
 cursor%p%val=val

end subroutine

subroutine addtotail_ll(sparse,row,col,val)
 class(llsparse),intent(inout),target::sparse
 integer(kind=int4),intent(in)::row,col
 real(kind=real8),intent(in)::val

 type(ptrnode),pointer::cursor

 if(.not.validvalue_gen(sparse,row,col))return
 if(.not.validnonzero_gen(sparse,val))return
 if(sparse%lupperstorage.and..not.uppervalue_gen(row,col))return

 cursor=>sparse%heads(row)
 do while(associated(cursor%p))
  cursor=>cursor%p%next
 enddo
 allocate(cursor%p)
 cursor%p%col=col
 cursor%p%val=val

end subroutine

!**GET ELEMENTS
function get_ll(sparse,row,col) result(val)
 class(llsparse),intent(inout)::sparse
 integer(kind=int4),intent(in)::row,col
 real(kind=real8)::val

 integer(kind=int4)::trow,tcol
 integer(kind=int4)::i,un
 type(ptrnode),pointer::cursor
 type(ptrnode),target::replacecursor

 val=0_real8
 
 trow=row
 tcol=col
 if(sparse%lupperstorage.and.row.gt.col)then
  !swap row-col
  trow=col
  tcol=row
 endif

 !cursor=>sparse%heads(trow)
 replacecursor=sparse%heads(trow)
 cursor=>replacecursor
 
 do while(associated(cursor%p))
  if(cursor%p%col.eq.tcol)then
   val=cursor%p%val
   exit
  endif
  cursor=>cursor%p%next
 enddo

end function

!**NUMBER OF ELEMENTS
function totalnumberofelements_ptrnode(pnode) result(nel)
 class(ptrnode),intent(in),target::pnode
 integer(kind=int8)::nel

 integer(kind=int4)::i
 type(ptrnode),pointer::cursor

 nel=0
 cursor=>pnode
 do while(associated(cursor%p))
  nel=nel+1
  cursor=>cursor%p%next
 enddo

end function

function totalnumberofelements_ll(sparse) result(nel)
 class(llsparse),intent(in)::sparse
 integer(kind=int8)::nel

 integer(kind=int4)::i

 nel=0
 do i=1,sparse%dim1
  nel=nel+sparse%heads(i)%size()
 enddo

end function

!**PRINT
subroutine print_ll(sparse,lint,output)
 class(llsparse),intent(in)::sparse
 integer(kind=int4),intent(in),optional::output
 logical,intent(in),optional::lint

 integer(kind=int4)::i,un
 logical::linternal
 type(ptrnode),pointer::cursor
 type(ptrnode),target::replacecursor

 linternal=.true.
 if(present(lint))linternal=lint

 un=sparse%unlog
 if(present(output))un=output

 do i=1,sparse%dim1
  !cursor=>sparse%heads(i)
  replacecursor=sparse%heads(i)
  cursor=>replacecursor
  do while(associated(cursor%p))
   write(un,'(2(i0,x),g0)')i,cursor%p%col,cursor%p%val
   if(.not.linternal.and.sparse%lupperstorage.and.cursor%p%col.ne.i)then
    write(un,'(2(i0,x),g0)')cursor%p%col,i,cursor%p%val
   endif
   cursor=>cursor%p%next
  enddo
 enddo

end subroutine

subroutine printsquare_ll(sparse,output)
 class(llsparse),intent(inout)::sparse
 integer(kind=int4),intent(in),optional::output

 integer(kind=int4)::i,j,un
 real(kind=real8)::val
 real(kind=real8),allocatable::tmp(:)

 un=sparse%unlog
 if(present(output))un=output
 
 write(un,'(a)')' Warning: this subroutine is not implemented yet!'

 allocate(tmp(sparse%dim2))

end subroutine

!FINAL
subroutine deallocate_scal_ll(sparse)
 type(llsparse),intent(inout)::sparse

 call destroy_ll(sparse)

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!OTHER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
!CHECKS
function validvalue_gen(sparse,row,col) result(lvalid)
 class(gen_sparse),intent(in)::sparse
 integer(kind=int4),intent(in)::row,col
 logical::lvalid

 lvalid=.true.
 if((row.lt.1.or.row.gt.sparse%dim1).or.(col.lt.1.or.col.gt.sparse%dim2))lvalid=.false.

end function

function validnonzero_gen(sparse,val) result(lvalid)
 class(gen_sparse),intent(in)::sparse
 real(kind=real8),intent(in)::val
 logical::lvalid

 lvalid=.true.
 if((abs(val)<epsilon(val)))lvalid=.false.

end function

function uppervalue_gen(row,col) result(lvalid)
 integer(kind=int4),intent(in)::row,col
 logical::lvalid

 lvalid=.true.
 if(row.gt.col)lvalid=.false.

end function

!CONVERSIONS
subroutine convertfromlltocoo(othersparse,sparse)
 type(coosparse),intent(out)::othersparse
 type(llsparse),intent(in),target::sparse
 
 integer(kind=int4)::i
 type(ptrnode),pointer::cursor

 othersparse=coosparse(sparse%dim1,sparse%dim2,sparse%nonzero(),sparse%lupperstorage)

 do i=1,sparse%dim1
  cursor=>sparse%heads(i)
  do while(associated(cursor%p))
   call othersparse%add(i,cursor%p%col,cursor%p%val)
   cursor=>cursor%p%next
  enddo
 enddo

! call othersparse%print()

end subroutine

subroutine convertfromlltocrs(othersparse,sparse)
 type(crssparse),intent(out)::othersparse
 type(llsparse),intent(in),target::sparse
 
 integer(kind=int4)::i,ndiag,nel,col
 integer(kind=int4),allocatable::rowpos(:)
 type(ptrnode),pointer::cursor

 if(sparse%nonzero().ge.2_int8**31)then
  write(sparse%unlog,'(a)')' ERROR: impossible conversion due a too large number of non-zero elements'
  stop
 endif

 !Condition: all diagonal elements must be present

 !Number of elements=number of diagonal elements+number off-diagonal elements
 !from sparse
 ndiag=min(sparse%dim1,sparse%dim2)

 allocate(rowpos(sparse%dim1))
 rowpos=0
 do i=1,sparse%dim1
  cursor=>sparse%heads(i)
  do while(associated(cursor%p))
   if(i.ne.cursor%p%col)rowpos(i)=rowpos(i)+1
   cursor=>cursor%p%next
  enddo
 enddo

 nel=ndiag+sum(rowpos)

 othersparse=crssparse(sparse%dim1,nel,sparse%dim2,sparse%lupperstorage)

 !if(othersparse%ia(othersparse%dim1+1).ne.0)othersparse%ia=0

 !determine the number of non-zero off-diagonal elements per row
 othersparse%ia(2:othersparse%dim1+1)=rowpos

 !accumulate the number of elements and add diagonal elements
 othersparse%ia(1)=1
 do i=1,ndiag
  othersparse%ia(i+1)=othersparse%ia(i+1)+1+othersparse%ia(i)
  othersparse%ja(othersparse%ia(i))=i
 enddo

 !accumulate the number of elements (no diag elements to be added)
 do i=ndiag+1,othersparse%dim1
  othersparse%ia(i+1)=othersparse%ia(i+1)+othersparse%ia(i)
 enddo

 !add the non-zero elements to crs (othersparse)
 !allocate(rowpos(othersparse%dim1))
 rowpos=othersparse%ia(1:othersparse%dim1)
 do i=1,sparse%dim1
  cursor=>sparse%heads(i)
  do while(associated(cursor%p))
   col=cursor%p%col
   if(i.eq.col)then
    othersparse%a(othersparse%ia(i))=cursor%p%val
   else
    rowpos(i)=rowpos(i)+1
    othersparse%ja(rowpos(i))=col
    othersparse%a(rowpos(i))=cursor%p%val
   endif
   cursor=>cursor%p%next
  enddo
 enddo

 deallocate(rowpos)

! call othersparse%print()

end subroutine

subroutine convertfromcootocrs(othersparse,sparse)
 type(crssparse),intent(out)::othersparse
 type(coosparse),intent(in)::sparse
 
 integer(kind=int4)::i,ndiag,nel,row,col
 integer(kind=int4),allocatable::rowpos(:)
 integer(kind=int8)::i8

 if(sparse%nonzero().ge.2_int8**31)then
  write(sparse%unlog,'(a)')' ERROR: impossible conversion due a too large number of non-zero elements'
  stop
 endif

 !Condition: all diagonal elements must be present

 !Number of elements=number of diagonal elements+number off-diagonal elements
 !from sparse
 ndiag=min(sparse%dim1,sparse%dim2)

 allocate(rowpos(sparse%dim1))
 rowpos=0

 do i8=1_int8,sparse%nel
  row=sparse%ij(1,i8)
  if(row.ne.0.and.row.ne.sparse%ij(2,i8))then
   rowpos(row)=rowpos(row)+1
  endif
 enddo

 nel=ndiag+sum(rowpos)

 othersparse=crssparse(sparse%dim1,nel,sparse%dim2,sparse%lupperstorage)

 !if(othersparse%ia(othersparse%dim1+1).ne.0)othersparse%ia=0

 !determine the number of non-zero off-diagonal elements per row
 othersparse%ia(2:othersparse%dim1+1)=rowpos

 !accumulate the number of elements and add diagonal elements
 othersparse%ia(1)=1
 do i=1,ndiag
  othersparse%ia(i+1)=othersparse%ia(i+1)+1+othersparse%ia(i)
  othersparse%ja(othersparse%ia(i))=i   !set diagonal element for the case it would not be present
 enddo

 !accumulate the number of elements (no diag elements to be added)
 do i=ndiag+1,othersparse%dim1
  othersparse%ia(i+1)=othersparse%ia(i+1)+othersparse%ia(i)
 enddo

 !add the non-zero elements to crs (othersparse)
 !allocate(rowpos(othersparse%dim1))
 rowpos=othersparse%ia(1:othersparse%dim1)
 do i8=1_int8,sparse%nel
  row=sparse%ij(1,i8)
  if(row.gt.0)then
   col=sparse%ij(2,i8)
   if(row.eq.col)then
    othersparse%a(othersparse%ia(row))=sparse%a(i8)
   else
    rowpos(row)=rowpos(row)+1
    othersparse%ja(rowpos(row))=col
    othersparse%a(rowpos(row))=sparse%a(i8)
   endif
  endif
 enddo

 deallocate(rowpos)

! call othersparse%print()

end subroutine



end module

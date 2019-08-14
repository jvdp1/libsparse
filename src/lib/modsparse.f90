!> Module containing three types of sparse matrices:  
!> * Linked list (llsparse)  
!> * COOrdinate storage (coosparse)  
!> * Compressed Row Storage (crssparse)  

!> @todo Use of submodules (one submodule for each type of matrix)
!> @todo Implementation of parameterized derived-type declarations to allow single- and double-precision sparse matrices
!> @todo Implementation of ordering for ll and coo

#if (_PARDISO==1)
include 'mkl_pardiso.f90'
#endif

module modsparse
#if (_DP==0)
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real32
#else
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real64
#endif
 use modhash
#if (_METIS==1)
 use modmetis
 use iso_c_binding,only:c_int,c_ptr,c_null_ptr
#endif
#if (_METIS==1 .AND. _SPAINV==1)
 use modspainv
#endif
#if (_PARDISO==1)
 use mkl_pardiso
#endif
 !$ use omp_lib
 implicit none
 private
 public::llsparse,coosparse,crssparse
 public::assignment(=)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
 integer(kind=int32),parameter::typegen=1,typecoo=10,typecrs=20,typell=30

 !> @brief Generic object containing dimensions, storage format, and output unit
 type,abstract::gen_sparse
  private
  integer(kind=int32)::unlog=output_unit
  integer(kind=int32)::dim1,dim2
  integer(kind=int32),allocatable::perm(:)  !Ap(i,:)=A(perm(i),:)
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
  !> @brief Sets the output unit to value; e.g., call mat%setouputunit(unlog)
  procedure,public::setoutputunit
  !> @brief Sets the permutation vector; e.g., call mat%setpermutation(array)
  procedure,public::setpermutation
  procedure::destroy_gen_gen
 end type
 
 abstract interface
  subroutine destroy_gen(sparse)
   import::gen_sparse
   class(gen_sparse),intent(inout)::sparse
  end subroutine
  function get_gen(sparse,row,col) result(val)
   import::int32,wp,gen_sparse
   class(gen_sparse),intent(inout)::sparse
   integer(kind=int32),intent(in)::row,col
   real(kind=wp)::val
  end function
  function nonzero_gen(sparse) result(nel)
   import::int64,gen_sparse
   class(gen_sparse),intent(in)::sparse
   integer(kind=int64)::nel
  end function
  subroutine print_gen(sparse,lint,output)
   import::int32,gen_sparse
   class(gen_sparse),intent(in)::sparse
   integer(kind=int32),intent(in),optional::output
   logical,intent(in),optional::lint
  end subroutine
  subroutine printsquare_gen(sparse,output)
   import::int32,gen_sparse
   class(gen_sparse),intent(inout)::sparse
   integer(kind=int32),intent(in),optional::output
   end subroutine
 end interface


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!COO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
 real(kind=real32),parameter::maxratiofilled_par=0.80

 !> @brief Object for COOrdinate storage
 type,extends(gen_sparse)::coosparse
  private
  integer(kind=int32),allocatable::ij(:,:)
  integer(kind=int64)::nel
  integer(kind=int64)::filled
  real(kind=wp),allocatable::a(:)
  contains
  private
  !> @brief Adds the value val to mat(row,col); e.g., call mat\%add(row,col,val)
  procedure,public::add=>add_coo 
  !> @brief Deallocates the sparse matrix and sets to default values 
  procedure,public::destroy=>destroy_scal_coo
  procedure::diag_vect_coo
  procedure::diag_mat_coo
  !> @brief Gets the (upper) diagonal elements of a matrix; e.g., array=mat%diag()  OR mat=mat%diag(10) (to extract the diagonal + 10 off-diagonals)
  generic,public::diag=>diag_vect_coo,diag_mat_coo
  !> @brief Returns the value of mat(row,col); e.g., ...=mat\%get(row,col)
  procedure,public::get=>get_coo
  !> @brief Returns the number of non-zero elements
  procedure,public::nonzero=>totalnumberofelements_coo
  !> @brief Prints the sparse matrix to the output mat\%unlog
  procedure,public::print=>print_coo
  !> @brief Prints the sparse matrix in a rectangular/square format to the default output
  procedure,public::printsquare=>printsquare_coo
  !> @brief Saves the matrix (internal format) to stream file
  procedure,public::save=>save_coo
  !> @brief Sets an entry to a certain value (even if equal to 0); e.g., call mat\%set(row,col,val)
  procedure,public::set=>set_coo
  !> @brief Gets a submatrix from a sparse matrix
  procedure,public::submatrix=>submatrix_coo
  final::deallocate_scal_coo
 end type

! !> @brief Load a COO matrix from file
! interface cooload
!  module procedure load_coo
! end interface

 !> @brief Constructor; e.g., mat=coosparse(dim1,[dim2],[#elements],[upper_storage],[output_unit]) OR mat=coosparse('file',[output_unit])
 interface coosparse
  module procedure constructor_coo,load_coo
 end interface
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CRS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
 !> @brief Object for Compressed Row Storage
 type,extends(gen_sparse)::crssparse
  private
  integer(kind=int32),allocatable::ia(:)
  integer(kind=int32),allocatable::ja(:)
  real(kind=wp),allocatable::a(:)
  contains
  private
  !> @brief Adds the value val to mat(row,col); e.g., call mat\%add(row,col,val)
  procedure,public::add=>add_crs
  !> @brief Deallocates the sparse matrix and sets to default values 
  procedure,public::destroy=>destroy_scal_crs
  procedure::diag_vect_crs
  procedure::diag_mat_crs
  !> @brief Gets the (upper) diagonal elements of a matrix; e.g., array=mat%diag()  OR mat=mat%diag(10) (to extract the diagonal + 10 off-diagonals)
  generic,public::diag=>diag_vect_crs,diag_mat_crs
  !> @brief Returns the value of mat(row,col); e.g., ...=mat\%get(row,col)
  procedure,public::get=>get_crs
  !> @brief Initiate the vectors ia,ja,and a from external vectors
  procedure,public::external=>external_crs
  !> @brief Multiplication with a vector
  procedure,public::multbyv=>multgenv_csr
  !> @brief Returns the number of non-zero elements
  procedure,public::nonzero=>totalnumberofelements_crs
#if (_METIS==1)
  !> @brief Returns the ordering array obtained from METIS
  procedure,public::getordering=>getordering_crs
#endif
  !> @brief Prints the sparse matrix to the output sparse\%unlog
  procedure,public::print=>print_crs
  !> @brief Prints the sparse matrix in a rectangular/square format to the default output
  procedure,public::printsquare=>printsquare_crs
  !> @brief Saves the matrix (internal format) to stream file
  procedure,public::save=>save_crs
  !> @brief Sets an entry to a certain value (even if equal to 0); condition: the entry must exist; e.g., call mat\%set(row,col,val)
  procedure,public::set=>set_crs
  !> @brief MKL PARDISO solver
  procedure,public::solve=>solve_crs
  !> @brief Sorts the elements in a ascending order within a row
  procedure,public::sort=>sort_crs
  !> @brief Computes and replaces by the sparse inverse
#if (_METIS==1 .AND. _SPAINV==1)
  procedure,public::spainv=>getspainv_crs
#endif
  !> @brief Gets a submatrix from a sparse matrix
  procedure,public::submatrix=>submatrix_crs
  final::deallocate_scal_crs
 end type

! !> @brief Load a CRS matrix from file
! interface crsload
!  module procedure load_crs
! end interface

 !> @brief Constructor; e.g., mat=crsparse(dim1,#elements,[dim2],[upper_storage],[output_unit]) OR mat=crssparse('file',[output_unit])
 interface crssparse
  module procedure constructor_crs,load_crs
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
  integer(kind=int32)::col
  real(kind=wp)::val
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!METIS GRAPH DATA STRUCTURE!!!!!!!!!!!!!!!!!!!!!!!!!aaa
#if (_METIS==1)
 type::metisgraph
  private
  integer(kind=int32)::unlog
  integer(kind=int32)::nvertices,medges
  integer(kind=int32),allocatable::xadj(:),adjncy(:)
  type(c_ptr)::vwgt
  type(c_ptr)::adjwgt
  contains
  private
  procedure,public::destroy=>destroy_metisgraph
  final::deallocate_scall_metisgraph
 end type

 interface metisgraph
  module procedure constructor_metisgraph
 end interface
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GENERAL INTERFACES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
 !> @brief Converts sparse matrices from one format to another one; e.g., crsmat=coomat
 interface assignment(=)
  module procedure convertfromlltocoo,convertfromlltocrs&
                  ,convertfromcootocrs,convertfromcootoll&
#if (_METIS==1)
                  ,convertfromcrstometisgraph&
#endif
                  ,convertfromcrstocoo,convertfromcrstoll
 end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
!DESTROY
!> @brief Subroutine to reset/destroy a generic object
subroutine destroy_gen_gen(sparse)
 class(gen_sparse),intent(inout)::sparse

 sparse%namemat='UNKNOWN'
 sparse%dim1=-1
 sparse%dim2=-1
 sparse%unlog=6
 sparse%lupperstorage=.false.
 if(allocated(sparse%perm))deallocate(sparse%perm)

end subroutine

!**GET ELEMENTS
function getdim_gen(sparse,dim1) result(dimget)
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

!**PRINT
subroutine print_dim_gen(sparse)
 class(gen_sparse),intent(in)::sparse

 integer(kind=int64)::nel

 write(sparse%unlog,'(/" Type of the matrix           : ",a)')trim(sparse%namemat)
 write(sparse%unlog,'( "  Output unit                 : ",i0)')sparse%unlog
 write(sparse%unlog,'( "  Dimension of the matrix     : ",i0," x ",i0)')sparse%dim1,sparse%dim2
 write(sparse%unlog,'( "  Upper storage               : ",l1)')sparse%lupperstorage
 write(sparse%unlog,'( "  Number of non-zero elements : ",i0)')sparse%nonzero()
 write(sparse%unlog,'( "  Permutation array provided  : ",l1)')allocated(sparse%perm)
 
 select type(sparse)
  type is(coosparse)
   write(sparse%unlog,'( "  Size of the array           : ",i0)')sparse%nel
  class default
 end select
 write(sparse%unlog,'(a)')' '

end subroutine

subroutine printtofile_gen(sparse,namefile,lint)
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

subroutine printsquaretofile_gen(sparse,namefile)
 class(gen_sparse),intent(inout)::sparse
 character(len=*),intent(in)::namefile

 integer(kind=int32)::un
 logical::linternal

 open(newunit=un,file=namefile,status='replace',action='write')
 call sparse%printsquare(output=un)
 close(un)

end subroutine

!**SET OUTPUT UNIT
subroutine setoutputunit(sparse,unlog)
 class(gen_sparse),intent(inout)::sparse
 integer(kind=int32)::unlog

 sparse%unlog=unlog

end subroutine

!** SET PERMUTATION VECTOR
subroutine setpermutation(sparse,array)
 class(gen_sparse),intent(inout)::sparse
 integer(kind=int32)::array(:)

 if(size(array).ne.sparse%getdim(1))then
  write(sparse%unlog,'(a)')' ERROR: The permutation array has a wrong size.'
  stop
 endif

 !Probably pointer would be better???
 if(.not.allocated(sparse%perm))allocate(sparse%perm(sparse%getdim(1)))
 sparse%perm=array

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
 integer(kind=int32),intent(in)::m
 integer(kind=int32),intent(in),optional::n,unlog
 integer(kind=int64),intent(in),optional::nel
 logical,intent(in),optional::lupper

 sparse%namemat='COO'
 sparse%dim1=m
 sparse%dim2=m
 if(present(n))sparse%dim2=n

 sparse%filled=0_int64

 sparse%nel=roundinguppower2(100_int64)
 if(present(nel))sparse%nel=roundinguppower2(int(nel,int64))
 allocate(sparse%ij(2,sparse%nel),sparse%a(sparse%nel))
 sparse%ij=0
 sparse%a=0._wp

 sparse%lupperstorage=.false.
 if(present(lupper))sparse%lupperstorage=lupper

 if(present(unlog))sparse%unlog=unlog

end function

!**DESTROY
subroutine destroy_scal_coo(sparse)
 class(coosparse),intent(inout)::sparse

 call sparse%destroy_gen_gen()

 sparse%nel=-1_int64
 sparse%filled=-1_int64
 if(allocated(sparse%ij))deallocate(sparse%ij)
 if(allocated(sparse%a))deallocate(sparse%a)

end subroutine

!**DIAGONAL ELEMENTS
function diag_vect_coo(sparse) result(array)
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

function diag_mat_coo(sparse,noff) result(diagsparse)
 class(coosparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::noff
 type(coosparse)::diagsparse

 integer(kind=int32)::ndiag,i,j

 ndiag=min(sparse%dim1,sparse%dim2)

 diagsparse=coosparse(ndiag,ndiag,int(ndiag,int64),lupper=.true.)
 
 do i=1,ndiag
  call diagsparse%add(i,i,sparse%get(i,i))
  do j=i+1,i+noff
   call diagsparse%add(i,j,sparse%get(i,j))
  enddo
 enddo

end function

!**ADD ELEMENTS
recursive subroutine add_coo(sparse,row,col,val)
 class(coosparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::row,col
 real(kind=wp),intent(in)::val

 integer(kind=int64)::hash,i8
 real(kind=real32),parameter::maxratiofilled=maxratiofilled_par
 real(kind=real32)::ratiofilled
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
  do i8=1_int64,sparse%nel
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

 if(hash.gt.0_int64)then!.and.ratiofilled.le.maxratiofilled)then
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
 integer(kind=int32),intent(in)::row,col
 real(kind=wp)::val
 
 integer(kind=int32)::trow,tcol
 integer(kind=int64)::hash

 val=0.0_wp
 
 trow=row
 tcol=col
 if(sparse%lupperstorage.and.row.gt.col)then
  !swap row-col
  trow=col
  tcol=row
 endif

 hash=hashf(trow,tcol,sparse%ij,sparse%nel,sparse%filled,.true.)
 
 if(hash.gt.0_int64)val=sparse%a(hash)

end function

!**EXTERNAL

!**LOAD
function load_coo(namefile,unlog) result(sparse)
 type(coosparse)::sparse
 character(len=*),intent(in)::namefile
 integer(kind=int32),intent(in),optional::unlog

 integer(kind=int32)::un,dim1,dim2
 integer(kind=int64)::nonzero,nel
 logical::lupperstorage

 open(newunit=un,file=namefile,action='read',status='old',access='stream')!,buffered='yes')
 read(un)dim1
 if(dim1.ne.typecoo)then
  write(*,'(a)')' ERROR: the proposed file is not a COO file'
  stop
 endif
 read(un)dim1            !int32
 read(un)dim2            !int32
 read(un)nonzero         !int64
 read(un)nel             !int64
 read(un)lupperstorage   !logical

 if(present(unlog))then
  sparse=coosparse(dim1,dim2,nel,lupperstorage,unlog)
 else
  sparse=coosparse(dim1,dim2,nel,lupperstorage)
 endif

 sparse%filled=nonzero
 read(un)sparse%ij              !int32
 read(un)sparse%a               !wp
 close(un)

end function

!**MULTIPLICATIONS

!**NUMBER OF ELEMENTS
function totalnumberofelements_coo(sparse) result(nel)
 class(coosparse),intent(in)::sparse
 integer(kind=int64)::nel

 nel=sparse%filled

end function

!**PRINT
subroutine print_coo(sparse,lint,output)
 class(coosparse),intent(in)::sparse
 integer(kind=int32),intent(in),optional::output
 logical,intent(in),optional::lint

 integer(kind=int32)::un,row,col
 integer(kind=int64)::i8
 real(kind=wp)::val
 character(len=30)::frm='(2(i0,1x),g0)'
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
  write(un,frm)row,col,val
  if(.not.linternal.and.sparse%lupperstorage.and.row.ne.col)then
   write(un,frm)col,row,val
  endif
 enddo

end subroutine

subroutine printsquare_coo(sparse,output)
 class(coosparse),intent(inout)::sparse
 integer(kind=int32),intent(in),optional::output

 integer(kind=int32)::i,j,un
 real(kind=wp)::val
 real(kind=wp),allocatable::tmp(:)

 un=sparse%unlog
 if(present(output))un=output

 allocate(tmp(sparse%dim2))

 do i=1,sparse%dim1
  tmp=0._wp
  do j=1,sparse%dim2
   tmp(j)=sparse%get(i,j)
  enddo
  write(un,'(10000(f9.3,1x))')tmp
 enddo

 deallocate(tmp)

end subroutine

!**SAVE
subroutine save_coo(sparse,namefile)
 class(coosparse),intent(in)::sparse
 character(len=*),intent(in)::namefile

 integer(kind=int32)::un

 open(newunit=un,file=namefile,action='write',status='replace',access='stream')!,buffered='yes')
 write(un)typecoo                !int32
 write(un)sparse%dim1            !int32
 write(un)sparse%dim2            !int32
 write(un)sparse%nonzero()       !int64
 write(un)sparse%nel             !int64
 write(un)sparse%lupperstorage   !logical
 write(un)sparse%ij              !int32
 write(un)sparse%a               !wp
 close(un)

end subroutine

!**SET ELEMENTS
recursive subroutine set_coo(sparse,row,col,val)
 !from add_coo
 class(coosparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::row,col
 real(kind=wp),intent(in)::val

 integer(kind=int64)::hash,i8
 real(kind=real32),parameter::maxratiofilled=maxratiofilled_par
 real(kind=real32)::ratiofilled
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
  do i8=1_int64,sparse%nel
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

 if(hash.gt.0_int64)then!.and.ratiofilled.le.maxratiofilled)then
  sparse%a(hash)=val
 else
  !is it possible?
  write(sparse%unlog,*)' ERROR: unexpected'!,__FILE__,__LINE__
  stop
 endif
 
end subroutine

!**SOLVE

!**SORT ARRAY

!**SUBMATRIX
function submatrix_coo(sparse,startdim1,enddim1,startdim2,enddim2,lupper,unlog) result(subsparse)
 !Not programmed efficiently, but it should do the job
 class(coosparse),intent(in)::sparse
 type(coosparse)::subsparse
 integer(kind=int32),intent(in)::startdim1,enddim1,startdim2,enddim2
 integer(kind=int32),intent(in),optional::unlog
 logical,intent(in),optional::lupper
 
 integer(kind=int32)::i,j,k,un
 integer(kind=int64)::i8,nel=10000
 logical::lincludediag,lupperstorage

 if(.not.validvalue_gen(sparse,startdim1,startdim2))return
 if(.not.validvalue_gen(sparse,enddim1,enddim2))return
 
 !check if the submatrix include diagonal elements of sparse
 ! if yes -> lupperstorage
 ! if no  -> .not.lupperstorage
 lincludediag=.false.
 firstdim: do i=startdim1,enddim1
  do j=startdim2,enddim2
   if(i.eq.j)then
    lincludediag=.true.
    exit firstdim
   endif
  enddo
 enddo firstdim

 lupperstorage=sparse%lupperstorage
 if(present(lupper))then
  lupperstorage=lupper
 else
  if(sparse%lupperstorage)then
   lupperstorage=.false.
   if(lincludediag)lupperstorage=.true.
  endif
 endif


 if(present(unlog))then
  subsparse=coosparse(enddim1-startdim1+1,enddim2-startdim2+1,int(nel,int64),lupperstorage,unlog)
 else
  subsparse=coosparse(enddim1-startdim1+1,enddim2-startdim2+1,int(nel,int64),lupperstorage)
 endif


 if(sparse%lupperstorage.eq.lupperstorage.or.(sparse%lupperstorage.and..not.lincludediag))then
  ! upper -> upper  ||  full -> full
  do i8=1,sparse%nel
   i=sparse%ij(1,i8)
   if(i.eq.0)cycle
   j=sparse%ij(2,i8)
   if((i.ge.startdim1.and.i.le.enddim1).and.(j.ge.startdim2.and.j.le.enddim2))then
    call subsparse%add(i-startdim1+1,j-startdim2+1,sparse%a(i8))
   endif
  enddo
 elseif(sparse%lupperstorage.and..not.lupperstorage)then
  ! upper -> full
  do i8=1,sparse%nel
   i=sparse%ij(1,i8)
   if(i.eq.0)cycle
   j=sparse%ij(2,i8)
   if((i.ge.startdim1.and.i.le.enddim1).and.(j.ge.startdim2.and.j.le.enddim2))then
    call subsparse%add(i-startdim1+1,j-startdim2+1,sparse%a(i8))
   endif
   if(i.ne.k)then
    if((j.ge.startdim1.and.j.le.enddim1).and.(i.ge.startdim2.and.i.le.enddim2))then
     call subsparse%add(j-startdim1+1,i-startdim2+1,sparse%a(i8))
    endif
   endif
  enddo
 elseif(.not.sparse%lupperstorage.and.lupperstorage)then
  ! full -> upper
  do i8=1,sparse%nel
   i=sparse%ij(1,i8)
   if(i.eq.0)cycle
   j=sparse%ij(2,i8)
print*,i,j,sparse%a(i8)
   if((j-startdim2+1.ge.i-startdim1+1).and.(i.ge.startdim1.and.i.le.enddim1).and.(j.ge.startdim2.and.j.le.enddim2))then
    call subsparse%add(i-startdim1+1,j-startdim2+1,sparse%a(i8))
   endif
  enddo
 endif

end function

!FINAL
subroutine deallocate_scal_coo(sparse)
 type(coosparse),intent(inout)::sparse

 call destroy_scal_coo(sparse)

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CRS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
!**CONSTRUCTOR
function constructor_crs(m,nel,n,lupper,unlog) result(sparse)
 type(crssparse)::sparse
 integer(kind=int32),intent(in)::m
 integer(kind=int32),intent(in)::nel
 integer(kind=int32),intent(in),optional::n,unlog
 logical,intent(in),optional::lupper

 sparse%namemat='CRS'
 sparse%dim1=m
 sparse%dim2=m
 if(present(n))sparse%dim2=n
 
 allocate(sparse%ia(sparse%dim1+1),sparse%ja(nel),sparse%a(nel))
 sparse%ia=0
 sparse%ia(sparse%dim1+1)=-nel
 sparse%ja=0
 sparse%a=0._wp

 sparse%lupperstorage=.false.
 if(present(lupper))sparse%lupperstorage=lupper

 if(present(unlog))sparse%unlog=unlog

end function

!**DESTROY
subroutine destroy_scal_crs(sparse)
 class(crssparse),intent(inout)::sparse

 call sparse%destroy_gen_gen()

 if(allocated(sparse%ia))deallocate(sparse%ia)
 if(allocated(sparse%ja))deallocate(sparse%ja)
 if(allocated(sparse%a))deallocate(sparse%a)

end subroutine

!**DIAGONAL ELEMENTS
function diag_vect_crs(sparse) result(array)
 class(crssparse),intent(inout)::sparse
 real(kind=wp),allocatable::array(:)

 integer(kind=int32)::ndiag,i

 ndiag=min(sparse%dim1,sparse%dim2)

 allocate(array(ndiag))
 array=0.0_wp
 
 do i=1,ndiag
  array(i)=sparse%get(i,i)
 enddo

end function

function diag_mat_crs(sparse,noff) result(diagsparse)
 class(crssparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::noff
 type(crssparse)::diagsparse

 integer(kind=int32)::ndiag,i,j,k,startoff,endoff,nel
 integer(kind=int32),allocatable::rowpos(:)

 ndiag=min(sparse%dim1,sparse%dim2)

 allocate(rowpos(ndiag)) 
 rowpos=0

 !determine the number of elements per row
 do i=1,ndiag
  startoff=i+1
  endoff=i+noff
  do j=sparse%ia(i),sparse%ia(i+1)-1
   if(sparse%ja(j).ge.startoff.and.sparse%ja(j).le.endoff)then
    rowpos(i)=rowpos(i)+1
   endif
  enddo
 enddo

 nel=ndiag+sum(rowpos)

 diagsparse=crssparse(ndiag,nel,ndiag,lupper=.true.)

 !determine the number of non-zero off-diagonal elements per row
 diagsparse%ia(2:diagsparse%dim1+1)=rowpos

 !accumulate the number of elements and add diagonal elements
 diagsparse%ia(1)=1
 do i=1,ndiag
  diagsparse%ia(i+1)=diagsparse%ia(i+1)+1+diagsparse%ia(i)
  diagsparse%ja(diagsparse%ia(i))=i   !set diagonal element for the case it would not be present
 enddo

 !add the non-zero elements to crs (diagsparse)
 rowpos=diagsparse%ia(1:diagsparse%dim1)
 do i=1,ndiag
  startoff=i+1
  endoff=i+noff
  do j=sparse%ia(i),sparse%ia(i+1)-1
   k=sparse%ja(j)
   if(i.eq.k)then
    diagsparse%a(diagsparse%ia(i))=sparse%a(j)
   elseif(k.ge.startoff.and.k.le.endoff)then
    rowpos(i)=rowpos(i)+1
    diagsparse%ja(rowpos(i))=k
    diagsparse%a(rowpos(i))=sparse%a(j)
   endif
  enddo
 enddo

 deallocate(rowpos)

end function

!**ADD ELEMENTS
subroutine add_crs(sparse,row,col,val,error)
 !add a value only to an existing one
 class(crssparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::row,col
 integer(kind=int32),intent(out),optional::error
 real(kind=wp),intent(in)::val

 integer(kind=int32)::i
 integer(kind=int32)::ierror    !added: error=0;Not existing: error=-1;matrix not inited: error=-10
 
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
 integer(kind=int32),intent(in)::row,col
 real(kind=wp)::val
 
 integer(kind=int32)::i,trow,tcol

 val=0.0_wp
 
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

!**EXTERNAL
subroutine external_crs(sparse,ia,ja,a)
 class(crssparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::ia(:),ja(:)
 real(kind=wp),intent(in)::a(:)

 if(size(ia).ne.size(sparse%ia))then
  write(sparse%unlog,'(a)')' ERROR: The provided array ia is of a different size!'
  stop
 endif
 if(size(ja).ne.size(sparse%ja))then
  write(sparse%unlog,'(a)')' ERROR: The provided array ja is of a different size!'
  stop
 endif
 if(size(a).ne.size(sparse%a))then
  write(sparse%unlog,'(a)')' ERROR: The provided array a is of a different size!'
  stop
 endif
 
 sparse%ia=ia
 sparse%ja=ja
 sparse%a=a

end subroutine

!**LOAD
function load_crs(namefile,unlog)  result(sparse)
 type(crssparse)::sparse
 integer(kind=int32),intent(in),optional::unlog
 character(len=*),intent(in)::namefile

 integer(kind=int32)::un,dim1,dim2
 integer(kind=int64)::nonzero
 logical::lupperstorage

 open(newunit=un,file=namefile,action='read',status='old',access='stream')!,buffered='yes')
 read(un)dim1
 if(dim1.ne.typecrs)then
  write(*,'(a)')' ERROR: the proposed file is not a CRS file'
  stop
 endif
 read(un)dim1            !int32
 read(un)dim2            !int32
 read(un)nonzero         !int64
 read(un)lupperstorage   !logical

 if(present(unlog))then
  sparse=crssparse(dim1,int(nonzero,int32),dim2,lupperstorage,unlog)
 else
  sparse=crssparse(dim1,int(nonzero,int32),dim2,lupperstorage)
 endif

 read(un)sparse%ia              !int32
 read(un)sparse%ja              !int32
 read(un)sparse%a               !wp

 close(un)

end function

!**MULTIPLICATIONS
subroutine multgenv_csr(sparse,alpha,trans,x,val,y)
 !Computes y=val*y+alpha*sparse(tranposition)*x
 class(crssparse),intent(in)::sparse
 real(kind=wp),intent(in)::val,alpha
 real(kind=wp),intent(in)::x(:)
 real(kind=wp),intent(out)::y(:)
 character(len=1),intent(in)::trans

 character(len=1)::matdescra(6)


 if(.not.sparse%lupperstorage)then
  matdescra(1)='G'
 elseif(sparse%lupperstorage)then
  matdescra(1)='T'
  matdescra(2)='U'
  matdescra(3)='N'
 endif
 matdescra(4)='F'
 call mkl_dcsrmv(trans,sparse%dim1,sparse%dim2,alpha,matdescra,sparse%a,sparse%ja,sparse%ia(1:sparse%dim1),sparse%ia(2:sparse%dim1+1),x,val,y)
 
end subroutine

!**NUMBER OF ELEMENTS
function totalnumberofelements_crs(sparse) result(nel)
 class(crssparse),intent(in)::sparse
 integer(kind=int64)::nel

 nel=int(sparse%ia(sparse%dim1+1),int64)-1_int64

end function

!**GET ORDER
#if (_METIS==1)
function getordering_crs(sparse&
                          ,ctype,iptype,rtype,compress,ccorder&
                          ,pfactor,nseps,bglvl&
                          ) result(perm)
 class(crssparse),intent(in)::sparse
 integer(kind=int32),intent(in),optional::ctype,iptype,rtype,compress,ccorder,pfactor,nseps,bglvl
 integer(kind=int32),allocatable::perm(:)

 integer(kind=int32)::err
 integer(kind=int32)::pctype,piptype,prtype,pcompress,pccorder,ppfactor,pnseps,pbglvl
 integer(kind=int32),allocatable::options(:)
 integer(kind=int32),allocatable::iperm(:)
 type(metisgraph)::metis
#if (_VERBOSE>0)
 !$ real(kind=real64)::t1

 !$ t1=omp_get_wtime()
#endif
 
 metis=sparse

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,g0)')'METIS=CRS',': Elapsed time = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 pctype=METIS_CTYPE_RM
 if(present(ctype))pctype=ctype

 piptype=METIS_IPTYPE_EDGE
 if(present(iptype))piptype=iptype

 prtype=METIS_RTYPE_SEP2SIDED
 if(present(rtype))prtype=rtype

 pcompress=0
 if(present(compress))pcompress=compress

 pccorder=0
 if(present(ccorder))pccorder=ccorder

 ppfactor=0
 if(present(pfactor))ppfactor=pfactor

 if(sparse%getdim(1).lt.50000)then
  pnseps=1
 elseif(sparse%getdim(1).lt.200000)then
  pnseps=5
 else
  pnseps=10
 endif
 if(present(nseps))pnseps=nseps

#if (_VERBOSE>1)
 pbglvl=METIS_DBG_INFO
#endif
 if(present(bglvl))pbglvl=bglvl

 err=metis_setoptions(options&
                      ,ctype=pctype,iptype=piptype,rtype=prtype,compress=pcompress&
                      ,ccorder=pccorder,pfactor=ppfactor,nseps=pnseps,dbglvl=pbglvl&
                      )
 call metis_checkerror(err,sparse%unlog)

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,g0)')'METIS options setting',': Elapsed time = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif
 
 allocate(perm(metis%nvertices),iperm(metis%nvertices))

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,g0)')'METIS arrays allocation',': Elapsed time = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 err=metis_nodend(metis%nvertices,metis%xadj,metis%adjncy,metis%vwgt,options,perm,iperm)
 call metis_checkerror(err,sparse%unlog)

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,g0)')'METIS ordering',': Elapsed time = ',omp_get_wtime()-t1
#endif

end function
#endif 

!**GET SPARSE INVERSE
#if (_METIS==1 .AND. _SPAINV==1)
subroutine getspainv_crs(sparse)
 class(crssparse),intent(inout)::sparse

 type(metisgraph)::metis
#if (_VERBOSE>0)
 !$ real(kind=real64)::t1,t2

 !$ t1=omp_get_wtime()
 !$ t2=t1
#endif

 call sparse%sort()

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,g0)')'SPAINV CRS sorting',': Elapsed time = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 !Ordering
 if(.not.allocated(sparse%perm))call sparse%setpermutation(sparse%getordering())

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,g0)')'SPAINV CRS ordering',': Elapsed time = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 metis=sparse

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,g0)')'SPAINV METIS=CRS',': Elapsed time = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 call get_spainv(sparse%ia,sparse%ja,sparse%a,metis%xadj,metis%adjncy,sparse%perm,sparse%unlog)

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,g0)')'SPAINV CRS inversion',': Elapsed time = ',omp_get_wtime()-t1
 !$ write(sparse%unlog,'(x,a,t30,a,g0)')'SPAINV CRS inversion',': Total   time = ',omp_get_wtime()-t2
#endif

end subroutine
#endif
 
!**PRINT
subroutine print_crs(sparse,lint,output)
 class(crssparse),intent(in)::sparse
 integer(kind=int32),intent(in),optional::output
 logical,intent(in),optional::lint

 integer(kind=int32)::i
 integer(kind=int32)::un,j
 character(len=30)::frm='(2(i0,1x),g0)'
 logical::linternal

 linternal=.true.
 if(present(lint))linternal=lint

 un=sparse%unlog
 if(present(output))un=output

 do i=1,sparse%dim1
  do j=sparse%ia(i),sparse%ia(i+1)-1
   write(un,frm)i,sparse%ja(j),sparse%a(j)
   if(.not.linternal.and.sparse%lupperstorage.and.i.ne.sparse%ja(j))then
    write(un,frm)sparse%ja(j),i,sparse%a(j)
   endif
  enddo
 enddo

end subroutine

subroutine printsquare_crs(sparse,output)
 class(crssparse),intent(inout)::sparse
 integer(kind=int32),intent(in),optional::output

 integer(kind=int32)::i,j,un
 real(kind=wp)::val
 real(kind=wp),allocatable::tmp(:)

 un=sparse%unlog
 if(present(output))un=output

 allocate(tmp(sparse%dim2))

 do i=1,sparse%dim1
  tmp=0.0_wp
  !could be implemented in a more efficient way
  do j=1,sparse%dim2
   tmp(j)=sparse%get(i,j)
  enddo
  write(un,'(10000(f9.3,1x))')tmp
 enddo

 deallocate(tmp)

end subroutine

!**SAVE
subroutine save_crs(sparse,namefile)
 class(crssparse),intent(in)::sparse
 character(len=*),intent(in)::namefile

 integer(kind=int32)::un

 open(newunit=un,file=namefile,action='write',status='replace',access='stream')!,buffered='yes')
 write(un)typecrs                !int32
 write(un)sparse%dim1            !int32
 write(un)sparse%dim2            !int32
 write(un)sparse%nonzero()       !int64
 write(un)sparse%lupperstorage   !logical
 write(un)sparse%ia              !int32
 write(un)sparse%ja              !int32
 write(un)sparse%a               !wp
 close(un)

end subroutine

!**SET ELEMENTS
subroutine set_crs(sparse,row,col,val,error)
 !add a value only to an existing one
 class(crssparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::row,col
 integer(kind=int32),intent(out),optional::error
 real(kind=wp),intent(in)::val

 integer(kind=int32)::i
 integer(kind=int32)::ierror    !added: error=0;Not existing: error=-1;matrix not inited: error=-10
 
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

!**SOLVE
#if (_PARDISO==1)
subroutine solve_crs(sparse,x,y)
 !sparse*x=y
 class(crssparse),intent(inout)::sparse
 real(kind=wp),intent(out)::x(:)
 real(kind=wp),intent(inout)::y(:)

 !Pardiso variables
 integer(kind=int32)::mtype=-2
 !integer(kind=int32)::mtype=11
 integer(kind=int32)::solver=0,error,phase,maxfct=1,mnum=1,nrhs=1
 integer(kind=int32)::idum(1)
 integer(kind=int32),save::iparm(64),msglvl=1
 real(kind=wp)::ddum(1)
 type(MKL_PARDISO_HANDLE),allocatable,save::pt(:)
 logical,save::lpardisofirst=.true.

 integer(kind=int32)::i
 !$ real(kind=real64)::t1

 if(.not.sparse%lsquare())then
  write(sparse%unlog,'(a)')' Warning: the sparse matrix is not squared!'
  return
 endif

 if(lpardisofirst)then
  !$ t1=omp_get_wtime()
  !Preparation of Cholesky of A11 with Pardiso
  !initialize pt
  allocate(pt(64))
  do i = 1, 64
    pt(i)%DUMMY=0 
  enddo
 
  !initialize iparm
  call pardisoinit(pt,mtype,iparm)
!  do i=1,64
!   write(sparse%unlog,*)'iparm',i,iparm(i)
!  enddo
  
  !Ordering and factorization
  phase=12
  iparm(2)=3
  iparm(27)=1
#if (_DP==0)
  iparm(28)=1
#else
  iparm(28)=0
#endif
  write(sparse%unlog,'(a)')' Start ordering and factorization'
  if(allocated(sparse%perm))then
   iparm(5)=1;iparm(31)=0;iparm(36)=0
   call pardiso(pt,maxfct,mnum,mtype,phase,&
                sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,&
                sparse%perm,nrhs,iparm,msglvl,ddum,ddum,error)
  else
   call pardiso(pt,maxfct,mnum,mtype,phase,&
                sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,&
                idum,nrhs,iparm,msglvl,ddum,ddum,error)
  endif
  call checkparido(phase,error) 
 
  write(sparse%unlog,'(a,i0)')' Number of nonzeros in factors  = ',iparm(18)
  write(sparse%unlog,'(a,i0)')' Number of factorization MFLOPS = ',iparm(19)
  !$ write(sparse%unlog,'(a,g0)')' Elapsed time                   = ',omp_get_wtime()-t1
 endif 

 !Solving
 phase=33
 iparm(27)=0
 call pardiso(pt,maxfct,mnum,mtype,phase,&
              sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,&
              idum,nrhs,iparm,msglvl,y,x,error)
 call checkparido(phase,error) 

 msglvl=0
 lpardisofirst=.false.

contains

 subroutine checkparido(phase,error)
  integer(kind=int32),intent(in)::phase,error
  if(error.ne.0)then
   write(sparse%unlog,'(2(a,i0))')' The following error for phase ',phase,' was detected: ',error
   stop
  endif
 end subroutine

end subroutine
#else
subroutine solve_crs(sparse,x,y)
 !sparse*x=y
 class(crssparse),intent(in)::sparse
 real(kind=wp),intent(out)::x(:)
 real(kind=wp),intent(in)::y(:)

 write(sparse%unlog,'(a)')' Warning: Pardiso is not enabled! Array returned = rhs'
 x=y

end subroutine
#endif

!**SORT ARRAY
subroutine sort_crs(sparse)
 ! sort vectors ja and a by increasing order
 class(crssparse),intent(inout)::sparse

 integer(kind=int32)::dir,endd,i,j,k,n,start,stkpnt
 integer(kind=int32)::d1,d2,d3,dmnmx,tmp
 integer(kind=int32)::stack(2,32)
 integer(kind=int32),allocatable::d(:)
 integer(kind=int32),parameter::select=20
 real(kind=wp)::umnmx,tmpu
 real(kind=wp),allocatable::u(:)

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

!**SUBMATRIX
function submatrix_crs(sparse,startdim1,enddim1,startdim2,enddim2,lupper,unlog) result(subsparse)
 !Not programmed efficiently, but it should do the job
 class(crssparse),intent(in)::sparse
 type(crssparse)::subsparse
 integer(kind=int32),intent(in)::startdim1,enddim1,startdim2,enddim2
 integer(kind=int32),intent(in),optional::unlog
 logical,intent(in),optional::lupper
 
 integer(kind=int32)::i,j,k,nel,un
 logical::lincludediag,lupperstorage
 type(coosparse)::subcoo

 if(.not.validvalue_gen(sparse,startdim1,startdim2))return
 if(.not.validvalue_gen(sparse,enddim1,enddim2))return
 
 !check if the submatrix include diagonal elements of sparse
 ! if yes -> lupperstorage
 ! if no  -> .not.lupperstorage
 lincludediag=.false.
 firstdim: do i=startdim1,enddim1
  do j=startdim2,enddim2
   if(j.gt.i)exit
   if(i.eq.j)then
    lincludediag=.true.
    exit firstdim
   endif
  enddo
 enddo firstdim

 lupperstorage=sparse%lupperstorage
 if(present(lupper))then
  lupperstorage=lupper
 else
  if(sparse%lupperstorage)then
   lupperstorage=.false.
   if(lincludediag)lupperstorage=.true.
  endif
 endif

 !determine the number of elements
 !possibilities
 nel=0
 if(sparse%lupperstorage.eq.lupperstorage.or.(sparse%lupperstorage.and..not.lincludediag))then
  ! upper -> upper  ||  full -> full
  do i=startdim1,enddim1
   do j=sparse%ia(i),sparse%ia(i+1)-1
    k=sparse%ja(j)
    if(k.ge.startdim2.and.k.le.enddim2)then
     nel=nel+1
    endif
   enddo
  enddo
 elseif(sparse%lupperstorage.and..not.lupperstorage)then
  ! upper -> full
  do i=startdim1,enddim1
   do j=sparse%ia(i),sparse%ia(i+1)-1
    k=sparse%ja(j)
    if(k.ge.startdim2.and.k.le.enddim2)then
     nel=nel+1
    endif
   enddo
  enddo
  do i=startdim2,enddim2
   do j=sparse%ia(i),sparse%ia(i+1)-1
    k=sparse%ja(j)
    if(i.ne.k.and.k.ge.startdim1.and.k.le.enddim1)then
     nel=nel+1
    endif
   enddo
  enddo
 elseif(.not.sparse%lupperstorage.and.lupperstorage)then
  ! full -> upper
  do i=startdim1,enddim1
   do j=sparse%ia(i),sparse%ia(i+1)-1
    k=sparse%ja(j)
    if((k-startdim2+1.ge.i-startdim1+1).and.k.ge.startdim2.and.k.le.enddim2)then
     nel=nel+1
    endif
   enddo
  enddo
 endif
 
 !add the elements
 if((sparse%lupperstorage.eq.lupperstorage).or.(sparse%lupperstorage.and..not.lincludediag))then
  ! upper -> upper  ||  full -> full
  if(present(unlog))then
   subsparse=crssparse(enddim1-startdim1+1,nel,enddim2-startdim2+1,lupperstorage,unlog)
  else
   subsparse=crssparse(enddim1-startdim1+1,nel,enddim2-startdim2+1,lupperstorage)
  endif
  subsparse%ia=0
  nel=1
  un=0
  do i=startdim1,enddim1
   un=un+subsparse%ia(nel)
   nel=nel+1
   do j=sparse%ia(i),sparse%ia(i+1)-1
    k=sparse%ja(j)
    if(k.ge.startdim2.and.k.le.enddim2)then
     subsparse%ia(nel)=subsparse%ia(nel)+1
     subsparse%ja(un+subsparse%ia(nel))=k-startdim2+1
     subsparse%a(un+subsparse%ia(nel))=sparse%a(j)
    endif
   enddo
  enddo
  subsparse%ia(1)=1
  do i=2,subsparse%dim1+1
   subsparse%ia(i)=subsparse%ia(i)+subsparse%ia(i-1)
  enddo
 elseif(sparse%lupperstorage.and..not.lupperstorage)then
  ! upper -> full
  if(present(unlog))then
   subcoo=coosparse(enddim1-startdim1+1,enddim2-startdim2+1,int(nel,int64),lupperstorage,unlog)
  else
   subcoo=coosparse(enddim1-startdim1+1,enddim2-startdim2+1,int(nel,int64),lupperstorage)
  endif
  do i=1,sparse%dim1
   do j=sparse%ia(i),sparse%ia(i+1)-1
    k=sparse%ja(j)
    if(k.ge.startdim2.and.k.le.enddim2)then
     call subcoo%add(i-startdim1+1,k-startdim2+1,sparse%a(j))
    endif
    if(i.ne.k)then
     if((k.ge.startdim1.and.k.le.enddim1).and.(i.ge.startdim2.and.i.le.enddim2))then
      call subcoo%add(k-startdim1+1,i-startdim2+1,sparse%a(j))
     endif
    endif
   enddo
  enddo
  subsparse=subcoo
  call subcoo%destroy()
 elseif(.not.sparse%lupperstorage.and.lupperstorage)then
  ! full -> upper
  if(present(unlog))then
   subsparse=crssparse(enddim1-startdim1+1,nel,enddim2-startdim2+1,lupperstorage,unlog)
  else
   subsparse=crssparse(enddim1-startdim1+1,nel,enddim2-startdim2+1,lupperstorage)
  endif
  subsparse%ia=0
  nel=1
  un=0
  do i=startdim1,enddim1
   un=un+subsparse%ia(nel)
   nel=nel+1
   do j=sparse%ia(i),sparse%ia(i+1)-1
    k=sparse%ja(j)
    if((k-startdim2+1.ge.i-startdim1+1).and.k.ge.startdim2.and.k.le.enddim2)then
     subsparse%ia(nel)=subsparse%ia(nel)+1
     subsparse%ja(un+subsparse%ia(nel))=k-startdim2+1
     subsparse%a(un+subsparse%ia(nel))=sparse%a(j)
    endif
   enddo
  enddo
  subsparse%ia(1)=1
  do i=2,subsparse%dim1+1
   subsparse%ia(i)=subsparse%ia(i)+subsparse%ia(i-1)
  enddo
 endif

end function

!FINAL
subroutine deallocate_scal_crs(sparse)
 type(crssparse),intent(inout)::sparse

 call destroy_scal_crs(sparse)

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!LINKED LIST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
!**CONSTRUCTOR
function constructor_ll(m,n,lupper,unlog) result(sparse)
 type(llsparse)::sparse
 integer(kind=int32),intent(in)::m
 integer(kind=int32),intent(in),optional::n,unlog
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
 integer(kind=int32)::i

 call sparse%destroy_gen_gen()

 if(allocated(sparse%heads))then
  do i=1,size(sparse%heads)
   call destroy_scal_ptrnode(sparse%heads(i))
  enddo
  deallocate(sparse%heads)
 endif

end subroutine

!**DIAGONAL ELEMENTS

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
 integer(kind=int32),intent(in)::col
 real(kind=wp),intent(in)::val

 type(ptrnode)::cursor

 allocate(cursor%p)
 cursor%p%next=pnode
 cursor%p%col=col
 cursor%p%val=val
 pnode=cursor

end subroutine

subroutine addtohead_ll(sparse,row,col,val)
 class(llsparse),intent(inout),target::sparse
 integer(kind=int32),intent(in)::row,col
 real(kind=wp),intent(in)::val

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
 integer(kind=int32),intent(in)::row,col
 real(kind=wp),intent(in)::val

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
 integer(kind=int32),intent(in)::row,col
 real(kind=wp),intent(in)::val

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
 integer(kind=int32),intent(in)::row,col
 real(kind=wp)::val

 integer(kind=int32)::trow,tcol
 integer(kind=int32)::i,un
 type(ptrnode),pointer::cursor
 type(ptrnode),target::replacecursor

 val=0.0_wp
 
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

!**EXTERNAL

!**LOAD

!**MULTIPLICATIONS

!**NUMBER OF ELEMENTS
function totalnumberofelements_ptrnode(pnode) result(nel)
 class(ptrnode),intent(in),target::pnode
 integer(kind=int64)::nel

 integer(kind=int32)::i
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
 integer(kind=int64)::nel

 integer(kind=int32)::i

 nel=0
 do i=1,sparse%dim1
  nel=nel+sparse%heads(i)%size()
 enddo

end function

!**SAVE

!**SET ELEMENTS

!**SOLVE

!**SORT ARRAY

!**SUBMATRIX

!**PRINT
subroutine print_ll(sparse,lint,output)
 class(llsparse),intent(in)::sparse
 integer(kind=int32),intent(in),optional::output
 logical,intent(in),optional::lint

 integer(kind=int32)::i,un
 character(len=20)::frm='(2(i0,1x),g0)'
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
   write(un,frm)i,cursor%p%col,cursor%p%val
   if(.not.linternal.and.sparse%lupperstorage.and.cursor%p%col.ne.i)then
    write(un,frm)cursor%p%col,i,cursor%p%val
   endif
   cursor=>cursor%p%next
  enddo
 enddo

end subroutine

subroutine printsquare_ll(sparse,output)
 class(llsparse),intent(inout)::sparse
 integer(kind=int32),intent(in),optional::output

 integer(kind=int32)::i,j,un
 real(kind=wp)::val
 real(kind=wp),allocatable::tmp(:)

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!METIS GRAPH DATA STRUCTURE!!!!!!!!!!!!!!!!!!!!!!!!!aaa
#if (_METIS==1)
!**CONSTRUCTOR
function constructor_metisgraph(n,m,unlog) result(metis)
 type(metisgraph)::metis
 integer(kind=int32),intent(in)::n,m
 integer(kind=int32),intent(in),optional::unlog

 metis%unlog=6
 if(present(unlog))metis%unlog=unlog

 metis%nvertices=n
 metis%medges=m
 allocate(metis%xadj(metis%nvertices+1),metis%adjncy(2*metis%medges))
 metis%xadj=0
 metis%adjncy=0
 metis%vwgt=c_null_ptr
 metis%adjwgt=c_null_ptr

end function

!**DESTROY
subroutine destroy_metisgraph(metis)
 class(metisgraph),intent(inout)::metis

 metis%unlog=6
 metis%nvertices=-1
 metis%medges=-1
 metis%vwgt=c_null_ptr
 metis%adjwgt=c_null_ptr
 if(allocated(metis%xadj))deallocate(metis%xadj)
 if(allocated(metis%adjncy))deallocate(metis%adjncy)

end subroutine

!FINAL
subroutine deallocate_scall_metisgraph(metis)
 type(metisgraph),intent(inout)::metis
 
 call destroy_metisgraph(metis)

end subroutine
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!OTHER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
!CHECKS
function validvalue_gen(sparse,row,col) result(lvalid)
 class(gen_sparse),intent(in)::sparse
 integer(kind=int32),intent(in)::row,col
 logical::lvalid

 lvalid=.true.
 if((row.lt.1.or.row.gt.sparse%dim1).or.(col.lt.1.or.col.gt.sparse%dim2))lvalid=.false.

end function

function validnonzero_gen(sparse,val) result(lvalid)
 class(gen_sparse),intent(in)::sparse
 real(kind=wp),intent(in)::val
 logical::lvalid

 lvalid=.true.
 if((abs(val)<epsilon(val)))lvalid=.false.

end function

function uppervalue_gen(row,col) result(lvalid)
 integer(kind=int32),intent(in)::row,col
 logical::lvalid

 lvalid=.true.
 if(row.gt.col)lvalid=.false.

end function

!CONVERSIONS
subroutine convertfromlltocoo(othersparse,sparse)
 type(coosparse),intent(out)::othersparse
 type(llsparse),intent(in),target::sparse
 
 integer(kind=int32)::i
 type(ptrnode),pointer::cursor

 othersparse=coosparse(sparse%dim1,sparse%dim2,sparse%nonzero(),sparse%lupperstorage)

 if(allocated(sparse%perm))allocate(othersparse%perm,source=sparse%perm)

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
 
 integer(kind=int32)::i,ndiag,nel,col
 integer(kind=int32),allocatable::rowpos(:)
 type(ptrnode),pointer::cursor

 if(sparse%nonzero().ge.2_int64**31)then
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

 if(allocated(sparse%perm))allocate(othersparse%perm,source=sparse%perm)

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
 
 integer(kind=int32)::i,ndiag,nel,row,col
 integer(kind=int32),allocatable::rowpos(:)
 integer(kind=int64)::i8

 if(sparse%nonzero().ge.2_int64**31)then
  write(sparse%unlog,'(a)')' ERROR: impossible conversion due a too large number of non-zero elements'
  stop
 endif

 !Condition: all diagonal elements must be present

 !Number of elements=number of diagonal elements+number off-diagonal elements
 !from sparse
 ndiag=min(sparse%dim1,sparse%dim2)

 allocate(rowpos(sparse%dim1))
 rowpos=0

 do i8=1_int64,sparse%nel
  row=sparse%ij(1,i8)
  if(row.ne.0.and.row.ne.sparse%ij(2,i8))then
   rowpos(row)=rowpos(row)+1
  endif
 enddo

 nel=ndiag+sum(rowpos)

 othersparse=crssparse(sparse%dim1,nel,sparse%dim2,sparse%lupperstorage)

 if(allocated(sparse%perm))allocate(othersparse%perm,source=sparse%perm)

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
 do i8=1_int64,sparse%nel
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

subroutine convertfromcootoll(othersparse,sparse)
 type(llsparse),intent(out)::othersparse
 type(coosparse),intent(in)::sparse
 
 integer(kind=int32)::row
 integer(kind=int64)::i8

 othersparse=llsparse(sparse%dim1,sparse%dim2,sparse%lupperstorage)

 if(allocated(sparse%perm))allocate(othersparse%perm,source=sparse%perm)

 do i8=1_int64,sparse%nel
  row=sparse%ij(1,i8)
  if(row.ne.0)then
   call othersparse%add(row,sparse%ij(2,i8),sparse%a(i8))
  endif
 enddo

! call othersparse%print()

end subroutine

subroutine convertfromcrstocoo(othersparse,sparse)
 type(coosparse),intent(out)::othersparse
 type(crssparse),intent(in)::sparse
 
 integer(kind=int32)::i,j

 othersparse=coosparse(sparse%dim1,sparse%dim2,sparse%nonzero(),sparse%lupperstorage)

 if(allocated(sparse%perm))allocate(othersparse%perm,source=sparse%perm)

 do i=1,sparse%dim1
  do j=sparse%ia(i),sparse%ia(i+1)-1
   call othersparse%add(i,sparse%ja(j),sparse%a(j))
  enddo
 enddo

! call othersparse%print()

end subroutine

subroutine convertfromcrstoll(othersparse,sparse)
 type(llsparse),intent(out)::othersparse
 type(crssparse),intent(in)::sparse
 
 integer(kind=int32)::i,j

 othersparse=llsparse(sparse%dim1,sparse%dim2,sparse%lupperstorage)

 if(allocated(sparse%perm))allocate(othersparse%perm,source=sparse%perm)

 do i=1,sparse%dim1
  do j=sparse%ia(i),sparse%ia(i+1)-1
   call othersparse%add(i,sparse%ja(j),sparse%a(j))
  enddo
 enddo

! call othersparse%print()

end subroutine

#if (_METIS==1)
subroutine convertfromcrstometisgraph(metis,sparse)
 type(metisgraph),intent(out)::metis
 type(crssparse),intent(in)::sparse

 integer(kind=int32)::i,j,k
 integer(kind=int32)::n,nvertices,medges
 integer(kind=int32),allocatable::rowpos(:)
 
 if(.not.sparse%lupperstorage.or..not.sparse%lsquare())then
  write(sparse%unlog,'(a)')' ERROR: the CRS matrix must be square and upper triangular stored'
  stop
 endif

 n=sparse%getdim(1)
 nvertices=n

 medges=sparse%nonzero()-n  !number of off-diagonal elements
 
 metis=metisgraph(nvertices,medges,unlog=sparse%unlog)

 allocate(rowpos(n))
 rowpos=0
 do i=1,n
  do j=sparse%ia(i),sparse%ia(i+1)-1
   k=sparse%ja(j)
   if(i.ne.k)then
    rowpos(i)=rowpos(i)+1
    rowpos(k)=rowpos(k)+1
   endif
  enddo
 enddo

 metis%xadj(1)=1
 do i=1,n
  metis%xadj(i+1)=metis%xadj(i)+rowpos(i)
 enddo

 rowpos=0
 do i=1,n
  do j=sparse%ia(i),sparse%ia(i+1)-1
   k=sparse%ja(j)
   if(i.ne.k)then
     metis%adjncy(metis%xadj(i)+rowpos(i))=k
     metis%adjncy(metis%xadj(k)+rowpos(k))=i
     rowpos(i)=rowpos(i)+1
     rowpos(k)=rowpos(k)+1
   endif
  enddo
 enddo

end subroutine
#endif
end module

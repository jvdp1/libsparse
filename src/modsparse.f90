!> Module containing three types of sparse matrices:
!> * Linked list (llsparse)
!> * COOrdinate storage (coosparse)
!> * Compressed Row Storage (crssparse)

!> @todo Implementation of parameterized derived-type declarations to allow single- and double-precision sparse matrices
!> @todo Implementation of ordering for ll and coo
!> @todo Generalization of the solve_crs method (for different solvers)

module modsparse
#if (_DP==0)
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real32
#else
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real64
#endif
 use iso_c_binding,only:c_ptr,c_null_ptr
#if (_PARDISO==1)
 use modvariablepardiso, only: pardiso_variable, pardiso_variable_64
#endif
 !$ use omp_lib
 implicit none
 private
 public::gen_sparse
 public::llsparse,coosparse,crssparse,crssparse64
 public::assignment(=)
 public::deallocate_sparse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!GEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
 integer(kind=int32),parameter::typegen=1,typecoo=10,typecrs=20,typell=30
 integer(kind=int32),parameter::typecrs64=21

 real(kind=wp),parameter::tol=1.e-10_wp

 !> @brief Generic object containing dimensions, storage format, and output unit
 type,abstract::gen_sparse
  private
  integer(kind=int32)::unlog=output_unit
  integer(kind=int32)::dim1=-1,dim2=-1
  integer(kind=int32),allocatable::perm(:)  !Ap(i,:)=A(perm(i),:)
  integer(kind=int64),allocatable::perm64(:)  !Ap(i,:)=A(perm(i),:)
  character(len=15)::namemat='UNKNOWN'
  logical::lsorted=.false.
  logical::lsymmetric=.false.
  logical::lupperstorage=.false.
  contains
  private
  procedure(destroy_gen),public,deferred::destroy
  procedure(get_gen),public,deferred::get
  procedure(multbyv_gen),deferred::multbyv
  procedure(multbym_gen),deferred::multbym
  procedure(nonzero_gen),public,deferred::nonzero
  procedure(print_gen),public,deferred::print_all
  procedure(print_gen_idx),public,deferred::print_idx
  generic, public :: print => print_all, print_idx
  procedure(printsquare_gen),public,deferred::printsquare
  procedure(scale_gen),public,deferred::scale

 
  !> @brief Iterative solver with the conjugate gradient method
  procedure,public::cg=>cg_gen
  !> @brief Returns the dimension of the matrix; e.g., ...=mat%getdim(1)
  procedure,public::getdim=>getdim_gen
  !> @brief Returns the permutation vector; e.g., call mat%getpermutation(array)
  procedure :: getpermutation32
  procedure :: getpermutation64
  generic,public::getpermutation => getpermutation32, getpermutation64
  !> @brief Gets the output unit; e.g., unit = mat%getouputunit()
  procedure,public::getoutputunit
  !> @brief Initializes the values
  procedure,public::initialize=>init_gen
  !> @brief Returns true if the matrix is sorted; else returns false
  procedure,public::issorted
  !> @brief Returns true if square matrix; else returns false
  procedure,public::issquare
  !> @brief Multiplication with a vector or matrix
  generic,public::mult=>multbyv,multbym
  !> @brief Prints the sparse matrix to a file
  procedure,public::printtofile=>printtofile_gen
  !> @brief Prints the sparse matrix in a rectangular/square format to the output mat\%unlog
  procedure,public::printsquaretofile=>printsquaretofile_gen
  !> @brief Prints the current status of a sparse matrix to the output mat\%unlog ; e.g., call mat%printstats()
  procedure,public::printstats=>print_dim_gen
  !> @brief Sets the output unit to value; e.g., call mat%setouputunit(unlog)
  procedure,public::setoutputunit
  !> @brief Sets the permutation vector; e.g., call mat%setpermutation(array)
  procedure :: setpermutation32
  procedure :: setpermutation64
  generic,public::setpermutation => setpermutation32, setpermutation64
  procedure::setsorted
  !> @brief Sets the assumption that the matrix is symmetric, if the matrix is square (there is no other check)
  procedure,public::setsymmetric
  procedure::destroy_gen_gen
  procedure::getmem_gen
 end type

 abstract interface
  subroutine scale_gen(sparse, val)
   import::gen_sparse,wp
   class(gen_sparse),intent(inout)::sparse
   real(kind=wp),intent(in)::val
  end subroutine
  impure elemental subroutine destroy_gen(sparse)
   import::gen_sparse
   class(gen_sparse),intent(inout)::sparse
  end subroutine
  pure function get_gen(sparse,row,col) result(val)
   import::int32,wp,gen_sparse
   class(gen_sparse),intent(in)::sparse
   integer(kind=int32),intent(in)::row,col
   real(kind=wp)::val
  end function
  subroutine multbyv_gen(sparse,alpha,trans,x,val,y)
   import::wp,gen_sparse
   class(gen_sparse),intent(in)::sparse
   real(kind=wp),intent(in)::val,alpha
   real(kind=wp),intent(in)::x(:)
   real(kind=wp),intent(out)::y(:)
   character(len=1),intent(in)::trans
  end subroutine
  subroutine multbym_gen(sparse,alpha,trans,x,val,y)
   import::wp,gen_sparse
   class(gen_sparse),intent(in)::sparse
   real(kind=wp),intent(in)::val,alpha
   real(kind=wp),intent(in)::x(:,:)
   real(kind=wp),intent(out)::y(:,:)
   character(len=1),intent(in)::trans
  end subroutine
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
  subroutine print_gen_idx(sparse,lidx,lint,output)
   import::int32,gen_sparse
   class(gen_sparse),intent(in)::sparse
   logical,intent(in)::lidx(:) !!if .true. for both row and column the element will be printed; otherwise not
   integer(kind=int32),intent(in),optional::output
   logical,intent(in),optional::lint
  end subroutine
  subroutine printsquare_gen(sparse,output)
   import::int32,gen_sparse
   class(gen_sparse),intent(inout)::sparse
   integer(kind=int32),intent(in),optional::output
   end subroutine
 end interface

 interface
  !DESTROY
  !> @brief Subroutine to reset/destroy a generic object
  module elemental subroutine destroy_gen_gen(sparse)
   class(gen_sparse),intent(inout)::sparse
  end subroutine
  !**CONJUGATE GRADIENT
  module subroutine cg_gen(sparse,x,y,maxiter,tol)
   !sparse*x=y
   class(gen_sparse),intent(in)::sparse
   integer(kind=int32),intent(inout),optional::maxiter
   real(kind=wp),intent(inout)::x(:)
   real(kind=wp),intent(in)::y(:)
   real(kind=wp),intent(inout),optional::tol
  end subroutine
  !**GET ELEMENTS
  pure module function getdim_gen(sparse,dim1) result(dimget)
   class(gen_sparse),intent(in)::sparse
   integer(kind=int32),intent(in)::dim1
   integer(kind=int32)::dimget
  end function
  !GET MEMORY
  module function getmem_gen(sparse) result(getmem)
   class(gen_sparse),intent(in)::sparse
   integer(kind=int64)::getmem
  end function
  !** GET PERMUTATION VECTOR
  module subroutine getpermutation32(sparse,array)
   class(gen_sparse),intent(in)::sparse
   integer(kind=int32),allocatable,intent(out)::array(:)
  end subroutine
  module subroutine getpermutation64(sparse,array)
   class(gen_sparse),intent(in)::sparse
   integer(kind=int64),allocatable,intent(out)::array(:)
  end subroutine
  !**GET OUTPUT UNIT
  pure module function getoutputunit(sparse) result(val)
   class(gen_sparse),intent(in)::sparse
   integer(kind=int32)::val
  end function
  !INITIATE GEN SPARSE
  module subroutine init_gen(sparse,namemat,dim1,dim2)
   class(gen_sparse),intent(inout)::sparse
   integer(kind=int32),intent(in)::dim1,dim2
   character(len=*),intent(in)::namemat
  end subroutine
  !**PRINT
  module subroutine print_dim_gen(sparse)
   class(gen_sparse),intent(in)::sparse
  end subroutine
  module subroutine printtofile_gen(sparse,namefile,lint)
   class(gen_sparse),intent(in)::sparse
   character(len=*),intent(in)::namefile
   logical,intent(in),optional::lint
  end subroutine
  module subroutine printsquaretofile_gen(sparse,namefile)
   class(gen_sparse),intent(inout)::sparse
   character(len=*),intent(in)::namefile
  end subroutine
  !**SET OUTPUT UNIT
  pure module subroutine setoutputunit(sparse,unlog)
   class(gen_sparse),intent(inout)::sparse
   integer(kind=int32),intent(in)::unlog
  end subroutine
  !** SET PERMUTATION VECTOR
  module subroutine setpermutation32(sparse,array)
   class(gen_sparse),intent(inout)::sparse
   integer(kind=int32)::array(:)
  end subroutine
  module subroutine setpermutation64(sparse,array)
   class(gen_sparse),intent(inout)::sparse
   integer(kind=int64)::array(:)
  end subroutine
  ! SET THE STATUS SORTED
  pure module subroutine setsorted(sparse,ll)
   class(gen_sparse),intent(inout)::sparse
   logical,intent(in)::ll
  end subroutine
  ! SET THE STATUS SYMMETRIC
  module subroutine setsymmetric(sparse,ll)
   class(gen_sparse),intent(inout)::sparse
   logical,intent(in),optional::ll
  end subroutine
  !**OTHER
  pure module function issorted(sparse) result(ll)
   class(gen_sparse),intent(in)::sparse
   logical::ll
  end function
  pure module function issquare(sparse) result(ll)
   class(gen_sparse),intent(in)::sparse
   logical::ll
  end function
 end interface

 interface
  !CHECKS
  pure module function validvalue_gen(sparse,row,col) result(lvalid)
   class(gen_sparse),intent(in)::sparse
   integer(kind=int32),intent(in)::row,col
   logical::lvalid
  end function
  pure module function validnonzero_gen(sparse,val) result(lvalid)
   class(gen_sparse),intent(in)::sparse
   real(kind=wp),intent(in)::val
   logical::lvalid
  end function
  pure module function uppervalue_gen(row,col) result(lvalid)
   integer(kind=int32),intent(in)::row,col
   logical::lvalid
  end function
 end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!COO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa

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
  procedure,public::destroy=>destroy_coo
  procedure::diag_vect_coo
  procedure::diag_mat_coo
  !> @brief Gets the (upper) diagonal elements of a matrix; e.g., array=mat%diag() 
  !! OR mat=mat%diag(10) (to extract the diagonal + 10 off-diagonals)
  generic,public::diag=>diag_vect_coo,diag_mat_coo
  !> @brief Returns the value of mat(row,col); e.g., ...=mat\%get(row,col)
  procedure,public::get=>get_coo
  !> @brief Gets memory used
  procedure,public::getmem=>getmem_coo
  !> @brief Initializes coosparse
  procedure,public::init=>constructor_sub_coo
  !> @brief Multiplication with a vector (not fully implemented)
  procedure::multbyv=>multgenv_coo
  !> @brief Multiplication with a matrix (not fully implemented)
  procedure::multbym=>multgenm_coo
  !> @brief Returns the number of non-zero elements
  procedure,public::nonzero=>totalnumberofelements_coo
  !> @brief Prints the sparse matrix to the output mat\%unlog
  procedure,public::print_all=>print_coo
  !> @brief Prints the selected entries of the sparse matrix to the output mat\%unlog
  procedure,public::print_idx=>print_idx_coo
  !> @brief Prints the sparse matrix in a rectangular/square format to the default output
  procedure,public::printsquare=>printsquare_coo
  !> @brief Saves the matrix (internal format) to stream file
  procedure,public::save=>save_coo
  !> @brief Scales all entries of mat by real scalar val; e.g., call mat\%scale(val)
  procedure,public::scale=>scale_coo
  !> @brief Sets an entry to a certain value (even if equal to 0); e.g., call mat\%set(row,col,val)
  procedure,public::set=>set_coo
  !> @brief Gets a submatrix from a sparse matrix
  procedure,public::submatrix=>submatrix_coo
  !> @brief Gets a submatrix from a sparse matrix based on an index vector
  procedure,public::submatrix_index=>submatrix_index_coo
  final::deallocate_scal_coo
 end type

 interface
  !**CONSTRUCTOR
  module function constructor_coo(m,n,nel,lupper,unlog) result(sparse)
   type(coosparse)::sparse
   integer(kind=int32),intent(in)::m
   integer(kind=int32),intent(in),optional::n,unlog
   integer(kind=int64),intent(in),optional::nel
   logical,intent(in),optional::lupper
  end function
  module subroutine constructor_sub_coo(sparse,m,n,nel,lupper,unlog)
   class(coosparse),intent(out)::sparse
   integer(kind=int32),intent(in)::m
   integer(kind=int32),intent(in),optional::n,unlog
   integer(kind=int64),intent(in),optional::nel
   logical,intent(in),optional::lupper
  end subroutine
  !**DESTROY
  module elemental subroutine destroy_coo(sparse)
   class(coosparse),intent(inout)::sparse
  end subroutine
  !**DIAGONAL ELEMENTS
  module function diag_vect_coo(sparse) result(array)
   class(coosparse),intent(inout)::sparse
   real(kind=wp),allocatable::array(:)
  end function
  !**ADD ELEMENTS
  recursive module subroutine add_coo(sparse,row,col,val)
   class(coosparse),intent(inout)::sparse
   integer(kind=int32),intent(in)::row,col
   real(kind=wp),intent(in)::val
  end subroutine
  !**GET ELEMENTS
  pure module function get_coo(sparse,row,col) result(val)
   class(coosparse),intent(in)::sparse
   integer(kind=int32),intent(in)::row,col
   real(kind=wp)::val
  end function
  !** GET MEMORY
  module function getmem_coo(sparse) result(getmem)
   class(coosparse),intent(in)::sparse
   integer(kind=int64)::getmem
  end function
  !**EXTERNAL
  !**LOAD
  module function load_coo(namefile,unlog) result(sparse)
   type(coosparse)::sparse
   character(len=*),intent(in)::namefile
   integer(kind=int32),intent(in),optional::unlog
  end function
  !**MULTIPLICATIONS
  module subroutine multgenv_coo(sparse,alpha,trans,x,val,y)
   !Computes y=val*y+alpha*sparse(tranposition)*x
   class(coosparse),intent(in)::sparse
   real(kind=wp),intent(in)::val,alpha
   real(kind=wp),intent(in)::x(:)
   real(kind=wp),intent(out)::y(:)
   character(len=1),intent(in)::trans
  end subroutine
  module subroutine multgenm_coo(sparse,alpha,trans,x,val,y)
   !Computes y=val*y+alpha*sparse(tranposition)*x
   class(coosparse),intent(in)::sparse
   real(kind=wp),intent(in)::val,alpha
   real(kind=wp),intent(in)::x(:,:)
   real(kind=wp),intent(out)::y(:,:)
   character(len=1),intent(in)::trans
  end subroutine
  !**NUMBER OF ELEMENTS
  module function totalnumberofelements_coo(sparse) result(nel)
   class(coosparse),intent(in)::sparse
   integer(kind=int64)::nel
  end function
  !**PRINT
  module subroutine print_coo(sparse,lint,output)
   class(coosparse),intent(in)::sparse
   integer(kind=int32),intent(in),optional::output
   logical,intent(in),optional::lint
  end subroutine
  module subroutine print_idx_coo(sparse,lidx,lint,output)
   class(coosparse),intent(in)::sparse
   logical,intent(in)::lidx(:)
   integer(kind=int32),intent(in),optional::output
   logical,intent(in),optional::lint
  end subroutine
  module subroutine printsquare_coo(sparse,output)
   class(coosparse),intent(inout)::sparse
   integer(kind=int32),intent(in),optional::output
  end subroutine
  !**SAVE
  module subroutine save_coo(sparse,namefile)
   class(coosparse),intent(in)::sparse
   character(len=*),intent(in)::namefile
  end subroutine
  !**SCALE ALL ENTRIES
  module subroutine scale_coo(sparse,val)
   class(coosparse),intent(inout)::sparse
   real(kind=wp),intent(in)::val
  end subroutine
  !**SET ELEMENTS
  recursive module subroutine set_coo(sparse,row,col,val)
   !from add_coo
   class(coosparse),intent(inout)::sparse
   integer(kind=int32),intent(in)::row,col
   real(kind=wp),intent(in)::val
  end subroutine
  !**SOLVE
  !**SORT ARRAY
  !**SUBMATRIX
  module function submatrix_coo(sparse,startdim1,enddim1,startdim2,enddim2,lupper,unlog) result(subsparse)
   !Not programmed efficiently, but it should do the job
   class(coosparse),intent(in)::sparse
   type(coosparse)::subsparse
   integer(kind=int32),intent(in)::startdim1,enddim1,startdim2,enddim2
   integer(kind=int32),intent(in),optional::unlog
   logical,intent(in),optional::lupper
  end function
  module subroutine submatrix_index_coo(sparse,subsparse,indvector,sizeblock,unlog)
   !Not programmed efficiently, but it should do the job
   class(coosparse),intent(in)::sparse
   type(coosparse),intent(out)::subsparse
   integer(kind=int32),intent(in)::indvector(:)
   integer(kind=int32),intent(in),optional::sizeblock
   integer(kind=int32),intent(in),optional::unlog
  end subroutine
 end interface

! !> @brief Load a COO matrix from file
! interface cooload
!  module procedure load_coo
! end interface

 !> @brief Constructor; e.g., mat=coosparse(dim1,[dim2],[#elements],[upper_storage],[output_unit])
 !! OR mat=coosparse('file',[output_unit])
 interface coosparse
  module procedure constructor_coo,load_coo
 end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!CRS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
 !> @brief Object for Compressed Row Storage
 type,extends(gen_sparse)::crssparse
  private
  integer(kind=int32),allocatable::ia(:)
  integer(kind=int32),allocatable::ja(:)
  real(kind=wp),allocatable::a(:)
  logical::loriginal=.true.
#if (_PARDISO==1)
  type(pardiso_variable)::pardisovar
#endif
  contains
  private
  !> @brief Adds the value val to mat(row,col); e.g., call mat\%add(row,col,val)
  procedure,public::add=>add_crs
  !> @brief Computes and replaces the sparse matrix by the (complete) Cholesky factor
  procedure,public::chol=>getchol_crs
  !> @brief Deallocates the sparse matrix and sets to default values
  procedure,public::destroy=>destroy_crs
  procedure::diag_vect_crs
  procedure::diag_mat_crs
  !> @brief Gets the (upper) diagonal elements of a matrix; e.g., array=mat%diag()
  !! OR mat=mat%diag(10) (to extract the diagonal + 10 off-diagonals)
  generic,public::diag=>diag_vect_crs,diag_mat_crs
  !> @brief Returns the value of mat(row,col); e.g., ...=mat\%get(row,col)
  procedure,public::get=>get_crs
  !> @brief Computes and replaces the sparse matrix by the (complete) LDLt (L is stored in the upper triangle and D in the diagonal)
  procedure,public::getldlt=>getldlt_crs
  !> @brief Gets memory used
  procedure,public::getmem=>getmem_crs
  !> @brief Initializes the vectors ia,ja,and a from external vectors
  procedure,public::external=>external_crs
  !> @brief Get function for the internal vector ia of row pointers
  procedure,public::get_rowptr=>get_rowptr_crs
  !> @brief Get function for the internal vector ja of column values
  procedure,public::get_colval=>get_colval_crs
  !> @brief Get function for the internal vector a of non-zero values
  procedure,public::get_nzval=>get_nzval_crs
  !> @brief Get diagonal elements of an approximate inverse using Harville (1999)
  procedure,public::harville=>harville_crs
  !> @brief Computes and replaces the sparse matrix by an incomplete Cholesky factor
  procedure,public::ichol=>getichol_crs
  !> @brief Iniates crssparse
  procedure,public::init=>constructor_sub_crs
  !> @brief Solver with a triangular factor (e.g., a Cholesky factor or an incomplete Cholesky factor)
  procedure,public::isolve=>isolve_crs
  !> @brief Solver using pre-computed LDLt decomposition
  procedure,public::solveldlt_s=>solveldlt_s_crs
  !> @brief Solver using LDLt decomposition (after its computation if needed
  procedure,private::solveldlt_crs_vector
  procedure,private::solveldlt_crs_array
  generic,public::solveldlt=>solveldlt_crs_vector,solveldlt_crs_array
  !> @brief Multiplication with a vector
  procedure::multbyv=>multgenv_csr
  !> @brief Multiplication with a matrix
  procedure::multbym=>multgenm_csr
  !> @brief Returns the number of non-zero elements
  procedure,public::nonzero=>totalnumberofelements_crs
#if (_METIS==1)
  !> @brief Returns the ordering array obtained from METIS
  procedure,public::getordering=>getordering_crs
#endif
  !> @brief Releases Pardiso memory if possible
  procedure,public::isdecomposed=>isdecomposed_crs
  procedure,public::setdecomposed=>setdecomposed_crs
#if (_PARDISO==1)
  procedure,public::resetpardiso=>reset_pardiso_memory_crs
#endif
  !> @brief Prints the sparse matrix to the output sparse\%unlog
  procedure,public::print_all=>print_crs
  !> @brief Prints the selected entries of the sparse matrix to the output mat\%unlog
  procedure,public::print_idx=>print_idx_crs
  !> @brief Prints the sparse matrix in a rectangular/square format to the default output
  procedure,public::printsquare=>printsquare_crs
  !> @brief Saves the matrix (internal format) to stream file
  procedure,public::save=>save_crs
  !> @brief Scales all entries of mat by real scalar val; e.g., call mat\%scale(val)
  procedure,public::scale=>scale_crs
  !> @brief Sets an entry to a certain value (even if equal to 0); condition: the entry must exist; e.g., call mat\%set(row,col,val)
  procedure,public::set=>set_crs
  !> @brief MKL PARDISO solver
  procedure,private::solve_crs_vector
  procedure,private::solve_crs_array
  generic,public::solve=>solve_crs_vector,solve_crs_array
  !> @brief Sorts the elements in a ascending order within a row
  procedure,public::sort=>sort_crs
  !> @brief Computes and replaces by the sparse inverse
  procedure,public::spainv=>getspainv_crs
  !> @brief Gets a submatrix from a sparse matrix
  procedure,public::submatrix=>submatrix_crs
  !> @brief Gets a dense submatrix from a sparse matrix
  procedure,public::submatrix_dense=>submatrix_dense_crs
  final::deallocate_scal_crs
 end type

 interface
  !**CONSTRUCTOR
  module function constructor_crs(m,nel,n,lupper,unlog) result(sparse)
   type(crssparse)::sparse
   integer(kind=int32),intent(in)::m
   integer(kind=int32),intent(in)::nel
   integer(kind=int32),intent(in),optional::n,unlog
   logical,intent(in),optional::lupper
  end function
  module subroutine constructor_sub_crs(sparse,m,nel,n,lupper,unlog)
   class(crssparse),intent(out)::sparse
   integer(kind=int32),intent(in)::m
   integer(kind=int32),intent(in)::nel
   integer(kind=int32),intent(in),optional::n,unlog
   logical,intent(in),optional::lupper
  end subroutine
  !**DESTROY
  module impure elemental subroutine destroy_crs(sparse)
   class(crssparse),intent(inout)::sparse
  end subroutine
  !**DIAGONAL ELEMENTS
  module function diag_vect_crs(sparse) result(array)
   class(crssparse),intent(inout)::sparse
   real(kind=wp),allocatable::array(:)
  end function
  !**ADD ELEMENTS
  module subroutine add_crs(sparse,row,col,val,error)
   !add a value only to an existing one
   class(crssparse),intent(inout)::sparse
   integer(kind=int32),intent(in)::row,col
   integer(kind=int32),intent(out),optional::error
   real(kind=wp),intent(in)::val
  end subroutine
  !**GET ELEMENTS
  pure module function get_crs(sparse,row,col) result(val)
   class(crssparse),intent(in)::sparse
   integer(kind=int32),intent(in)::row,col
   real(kind=wp)::val
  end function
  !** GET MEMORY
  module function getmem_crs(sparse) result(getmem)
   class(crssparse),intent(in)::sparse
   integer(kind=int64)::getmem
  end function
  !**EXTERNAL
  module subroutine external_crs(sparse,ia,ja,a)
   class(crssparse),intent(inout)::sparse
   integer(kind=int32),intent(in)::ia(:),ja(:)
   real(kind=wp),intent(in)::a(:)
  end subroutine
  !**HARVILLE
  module subroutine harville_crs(sparse, ngibbs, nburn, diaginv, seed)
   class(crssparse), intent(inout)::sparse
   integer, intent(in) :: ngibbs, nburn
   real(wp), intent(out), allocatable :: diaginv(:)
   integer, intent(in), optional :: seed
  end subroutine
  !**ROWPTR
  pure module subroutine get_rowptr_crs(sparse,ia)
    class(crssparse),intent(in)::sparse
    integer(kind=int32),allocatable,intent(out)::ia(:)
  end subroutine
  !**COLVAL
  pure module subroutine get_colval_crs(sparse,ja)
    class(crssparse),intent(in)::sparse
    integer(kind=int32),allocatable,intent(out)::ja(:)
   end subroutine
  !**NZVAL
  pure module subroutine get_nzval_crs(sparse,a)
    class(crssparse),intent(in)::sparse
    real(kind=wp),allocatable,intent(out)::a(:)
  end subroutine  
  !**MULTIPLICATIONS
  module subroutine multgenv_csr(sparse,alpha,trans,x,val,y)
   !Computes y=val*y+alpha*sparse(tranposition)*x
   class(crssparse),intent(in)::sparse
   real(kind=wp),intent(in)::val,alpha
   real(kind=wp),intent(in)::x(:)
   real(kind=wp),intent(out)::y(:)
   character(len=1),intent(in)::trans
  end subroutine
  module subroutine multgenm_csr(sparse,alpha,trans,x,val,y)
   !Computes y=val*y+alpha*sparse(tranposition)*x
   class(crssparse),intent(in)::sparse
   real(kind=wp),intent(in)::val,alpha
   real(kind=wp),intent(in)::x(:,:)
   real(kind=wp),intent(out)::y(:,:)
   character(len=1),intent(in)::trans
  end subroutine
  !**NUMBER OF ELEMENTS
  module function totalnumberofelements_crs(sparse) result(nel)
   class(crssparse),intent(in)::sparse
   integer(kind=int64)::nel
  end function
  !**GET (COMPLETE) CHOLESKY FACTOR
  module subroutine getchol_crs(sparse,minsizenode)
   class(crssparse),intent(inout)::sparse
   integer(kind=int32),intent(in),optional::minsizenode
  end subroutine
  !**GET LDLt DECOMPOSITION
  module subroutine getldlt_crs(sparse,minsizenode)
   class(crssparse),intent(inout)::sparse
   integer(kind=int32),intent(in),optional::minsizenode
  end subroutine
#if (_METIS==1)
  !**GET ORDER
  module function getordering_crs(sparse&
                            ,ctype,iptype,rtype,compress,ccorder&
                            ,pfactor,nseps,bglvl&
                            ) result(perm)
   class(crssparse),intent(in)::sparse
   integer(kind=int32),intent(in),optional::ctype,iptype,rtype,compress,ccorder,pfactor,nseps,bglvl
   integer(kind=int32),allocatable::perm(:)
  end function
#endif
  !**GET INCOMPLETE CHOLESKY FACTOR
  module subroutine getichol_crs(sparse,minsizenode)
   class(crssparse),intent(inout)::sparse
   integer(kind=int32),intent(in),optional::minsizenode
  end subroutine
  !**GET SPARSE INVERSE
  module subroutine getspainv_crs(sparse,minsizenode)
   class(crssparse),intent(inout)::sparse
   integer(kind=int32),intent(in),optional::minsizenode
  end subroutine
  pure module function isdecomposed_crs(sparse) result(ll)
   class(crssparse),intent(in)::sparse
   logical::ll
  end function
  pure module subroutine setdecomposed_crs(sparse,ll)
   class(crssparse),intent(inout)::sparse
   logical,intent(in)::ll
  end subroutine
#if (_PARDISO==1)
  !**RESET PARDISO MEMORY
  module subroutine reset_pardiso_memory_crs(sparse)
   !sparse*x=y
   class(crssparse),intent(inout)::sparse
  end subroutine
#endif
  !**PRINT
  module subroutine print_crs(sparse,lint,output)
   class(crssparse),intent(in)::sparse
   integer(kind=int32),intent(in),optional::output
   logical,intent(in),optional::lint
  end subroutine
  module subroutine print_idx_crs(sparse,lidx,lint,output)
   class(crssparse),intent(in)::sparse
   logical,intent(in)::lidx(:)
   integer(kind=int32),intent(in),optional::output
   logical,intent(in),optional::lint
  end subroutine
  module subroutine printsquare_crs(sparse,output)
   class(crssparse),intent(inout)::sparse
   integer(kind=int32),intent(in),optional::output
  end subroutine
  !**SAVE
  module subroutine save_crs(sparse,namefile)
   class(crssparse),intent(in)::sparse
   character(len=*),intent(in)::namefile
  end subroutine
  !**SCALE ALL ENTRIES
  module subroutine scale_crs(sparse,val)
   class(crssparse),intent(inout)::sparse
   real(kind=wp),intent(in)::val
  end subroutine
  !**SET ELEMENTS
  module subroutine set_crs(sparse,row,col,val,error)
   !add a value only to an existing one
   class(crssparse),intent(inout)::sparse
   integer(kind=int32),intent(in)::row,col
   integer(kind=int32),intent(out),optional::error
   real(kind=wp),intent(in)::val
  end subroutine
  !**SOLVE
  module subroutine solve_crs_vector(sparse,x,y,msglvl)
   !sparse*x=y
   class(crssparse),intent(inout)::sparse
   real(kind=wp),intent(out),contiguous::x(:)
   real(kind=wp),intent(inout),contiguous::y(:)
   integer(kind=int32),intent(in),optional::msglvl
  end subroutine
  module subroutine solve_crs_array(sparse,x,y,msglvl)
   !sparse*x=y
   class(crssparse),intent(inout)::sparse
   real(kind=wp),intent(out),contiguous::x(:,:)
   real(kind=wp),intent(inout),contiguous::y(:,:)
   integer(kind=int32),intent(in),optional::msglvl
  end subroutine
  !**SOLVE WITH A TRIANGULAR FACTOR
  module subroutine isolve_crs(sparse,x,y)
   !sparse*x=y
   class(crssparse),intent(in)::sparse
   real(kind=wp),intent(out)::x(:)
   real(kind=wp),intent(in)::y(:)
  end subroutine
  !**SOLVE WITH LDLt DECOMPOSITION
  module subroutine solveldlt_s_crs(sparse,x,y)
   !sparse*x=y
   class(crssparse),intent(in)::sparse
   real(kind=wp),intent(out)::x(:)
   real(kind=wp),intent(in)::y(:)
  end subroutine
  !**SOLVE WITH LDLt DECOMPOSITION (AND COMPUTE IT IF NEEDED)
  module subroutine solveldlt_crs_vector(sparse,x,y)
   !sparse*x=y
   class(crssparse),intent(inout)::sparse
   real(kind=wp),intent(out),contiguous::x(:)
   real(kind=wp),intent(inout),contiguous::y(:)
  end subroutine
  module subroutine solveldlt_crs_array(sparse,x,y)
   !sparse*x=y
   class(crssparse),intent(inout)::sparse
   real(kind=wp),intent(out),contiguous::x(:,:)
   real(kind=wp),intent(inout),contiguous::y(:,:)
  end subroutine
  !**SORT ARRAY
  module subroutine sort_crs(sparse)
   ! sort vectors ja and a by increasing order
   class(crssparse),intent(inout)::sparse
  end subroutine
 end interface

! !> @brief Load a CRS matrix from file
! interface crsload
!  module procedure load_crs
! end interface

 !> @brief Constructor; e.g., mat=crsparse(dim1,#elements,[dim2],[upper_storage],[output_unit])
 !! OR mat=crssparse('file',[output_unit])
 interface crssparse
  module procedure constructor_crs,load_crs
 end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!CRS64!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
 !> @brief Object for Compressed Row Storage with Int64 indices
 type,extends(gen_sparse)::crssparse64
  private
  integer(kind=int64),allocatable::ia(:)
  integer(kind=int64),allocatable::ja(:)
  real(kind=wp),allocatable::a(:)
  logical::loriginal=.true.
#if (_PARDISO==1)
  type(pardiso_variable_64)::pardisovar
#endif
  contains
  private
  !> @brief Adds the value val to mat(row,col); e.g., call mat\%add(row,col,val)
  procedure,public::add=>add_crs64
!#if (_SPAINV==1)
!  !> @brief Computes and replaces the sparse matrix by the (complete) Cholesky factor
!  procedure,public::chol=>getchol_crs64
!#endif
  !> @brief Deallocates the sparse matrix and sets to default values
  procedure,public::destroy=>destroy_crs64
!  procedure::diag_vect_crs64
!  procedure::diag_mat_crs64
!  !> @brief Gets the (upper) diagonal elements of a matrix; e.g., array=mat%diag()
!  !! OR mat=mat%diag(10) (to extract the diagonal + 10 off-diagonals)
!  generic,public::diag=>diag_vect_crs64,diag_mat_crs64
  !> @brief Returns the value of mat(row,col); e.g., ...=mat\%get(row,col)
  procedure,public::get=>get_crs64
!#if (_SPAINV==1)
!  !> @brief Computes and replaces the sparse matrix by the (complete) LDLt
!  !!(L is stored in the upper triangle and D in the diagonal)
!  procedure,public::getldlt=>getldlt_crs64
!#endif
  !> @brief Gets memory used
  procedure,public::getmem=>getmem_crs64
!  !> @brief Initializes the vectors ia,ja,and a from external vectors
!  procedure,public::external=>external_crs64
  !> @brief Get function for the internal vector ia of row pointers
  procedure,public::get_rowptr=>get_rowptr_crs64
  !> @brief Get function for the internal vector ja of column values
  procedure,public::get_colval=>get_colval_crs64
  !> @brief Get function for the internal vector a of non-zero values
  procedure,public::get_nzval=>get_nzval_crs64
!  !> @brief Get diagonal elements of an approximate inverse using Harville (1999)
!  procedure,public::harville=>harville_crs64
!#if (_SPAINV==1)
!  !> @brief Computes and replaces the sparse matrix by an incomplete Cholesky factor
!  procedure,public::ichol=>getichol_crs64
!#endif
  !> @brief Iniates crssparse64
  procedure,public::init=>constructor_sub_crs64
!  !> @brief Solver with a triangular factor (e.g., a Cholesky factor or an incomplete Cholesky factor)
!  procedure,public::isolve=>isolve_crs64
!  !> @brief Solver using pre-computed LDLt decomposition
!  procedure,public::solveldlt_s=>solveldlt_s_crs64
  !> @brief Multiplication with a vector
  procedure::multbyv=>multgenv_csr64
  !> @brief Multiplication with a matrix
  procedure::multbym=>multgenm_csr64
  !> @brief Returns the number of non-zero elements
  procedure,public::nonzero=>totalnumberofelements_crs64
!#if (_METIS==1)
!  !> @brief Returns the ordering array obtained from METIS
!  procedure,public::getordering=>getordering_crs64
!#endif
!  !> @brief Releases Pardiso memory if possible
  procedure,public::isdecomposed=>isdecomposed_crs64
  procedure,public::setdecomposed=>setdecomposed_crs64
#if (_PARDISO==1)
  procedure,public::resetpardiso=>reset_pardiso_memory_crs64
#endif
  !> @brief Prints the sparse matrix to the output sparse\%unlog
  procedure,public::print_all=>print_crs64
  !> @brief Prints the selected entries of the sparse matrix to the output mat\%unlog
  procedure,public::print_idx=>print_idx_crs64
  !> @brief Prints the sparse matrix in a rectangular/square format to the default output
  procedure,public::printsquare=>printsquare_crs64
  !> @brief Saves the matrix (internal format) to stream file
  procedure,public::save=>save_crs64
  !> @brief Scales all entries of mat by real scalar val; e.g., call mat\%scale(val)
  procedure,public::scale=>scale_crs64
  !> @brief Sets an entry to a certain value (even if equal to 0); condition: the entry must exist; e.g., call mat\%set(row,col,val)
  procedure,public::set=>set_crs64
  !> @brief MKL PARDISO solver
  procedure,private::solve_crs64_vector
  procedure,private::solve_crs64_array
  generic,public::solve=>solve_crs64_vector,solve_crs64_array
  !> @brief Sorts the elements in a ascending order within a row
  procedure,public::sort=>sort_crs64
!#if (_SPAINV==1)
!  !> @brief Computes and replaces by the sparse inverse
!  procedure,public::spainv=>getspainv_crs64
!#endif
!  !> @brief Gets a submatrix from a sparse matrix
!  procedure,public::submatrix=>submatrix_crs64
  final::deallocate_scal_crs64
 end type

 interface
  !**CONSTRUCTOR
  module function constructor_crs64(m,nel,n,lupper,unlog) result(sparse)
   type(crssparse64)::sparse
   integer(kind=int32),intent(in)::m
   integer(kind=int64),intent(in)::nel
   integer(kind=int32),intent(in),optional::n,unlog
   logical,intent(in),optional::lupper
  end function
  module subroutine constructor_sub_crs64(sparse,m,nel,n,lupper,unlog)
   class(crssparse64),intent(out)::sparse
   integer(kind=int32),intent(in)::m
   integer(kind=int64),intent(in)::nel
   integer(kind=int32),intent(in),optional::n,unlog
   logical,intent(in),optional::lupper
  end subroutine
  !**DESTROY
  module impure elemental subroutine destroy_crs64(sparse)
   class(crssparse64),intent(inout)::sparse
  end subroutine
!  !**DIAGONAL ELEMENTS
!  module function diag_vect_crs64(sparse) result(array)
!   class(crssparse64),intent(inout)::sparse
!   real(kind=wp),allocatable::array(:)
!  end function
  !**ADD ELEMENTS
  module subroutine add_crs64(sparse,row,col,val,error)
   !add a value only to an existing one
   class(crssparse64),intent(inout)::sparse
   integer(kind=int32),intent(in)::row,col
   integer(kind=int32),intent(out),optional::error
   real(kind=wp),intent(in)::val
  end subroutine
  !**GET ELEMENTS
  pure module function get_crs64(sparse,row,col) result(val)
   class(crssparse64),intent(in)::sparse
   integer(kind=int32),intent(in)::row,col
   real(kind=wp)::val
  end function
  !** GET MEMORY
  module function getmem_crs64(sparse) result(getmem)
   class(crssparse64),intent(in)::sparse
   integer(kind=int64)::getmem
  end function
!  !**EXTERNAL
!  module subroutine external_crs64(sparse,ia,ja,a)
!   class(crssparse64),intent(inout)::sparse
!   integer(kind=int32),intent(in)::ia(:),ja(:)
!   real(kind=wp),intent(in)::a(:)
!  end subroutine
!  !**HARVILLE
!  module subroutine harville_crs64(sparse, ngibbs, nburn, diaginv, seed)
!   class(crssparse64), intent(inout)::sparse
!   integer, intent(in) :: ngibbs, nburn
!   real(wp), intent(out), allocatable :: diaginv(:)
!   integer, intent(in), optional :: seed
!  end subroutine
  !**ROWPTR
  module subroutine get_rowptr_crs64(sparse,ia)
    class(crssparse64),intent(in)::sparse
    integer(kind=int64),allocatable,intent(out)::ia(:)
  end subroutine
  !**COLVAL
  module subroutine get_colval_crs64(sparse,ja)
    class(crssparse64),intent(in)::sparse
    integer(kind=int64),allocatable,intent(out)::ja(:)
   end subroutine
  !**NZVAL
  module subroutine get_nzval_crs64(sparse,a)
    class(crssparse64),intent(in)::sparse
    real(kind=wp),allocatable,intent(out)::a(:)
  end subroutine  
  !**MULTIPLICATIONS
  module subroutine multgenv_csr64(sparse,alpha,trans,x,val,y)
   !Computes y=val*y+alpha*sparse(tranposition)*x
   class(crssparse64),intent(in)::sparse
   real(kind=wp),intent(in)::val,alpha
   real(kind=wp),intent(in)::x(:)
   real(kind=wp),intent(out)::y(:)
   character(len=1),intent(in)::trans
  end subroutine
  module subroutine multgenm_csr64(sparse,alpha,trans,x,val,y)
   !Computes y=val*y+alpha*sparse(tranposition)*x
   class(crssparse64),intent(in)::sparse
   real(kind=wp),intent(in)::val,alpha
   real(kind=wp),intent(in)::x(:,:)
   real(kind=wp),intent(out)::y(:,:)
   character(len=1),intent(in)::trans
  end subroutine
  !**NUMBER OF ELEMENTS
  module function totalnumberofelements_crs64(sparse) result(nel)
   class(crssparse64),intent(in)::sparse
   integer(kind=int64)::nel
  end function
!#if (_SPAINV==1)
!  !**GET (COMPLETE) CHOLESKY FACTOR
!  module subroutine getchol_crs64(sparse,minsizenode)
!   class(crssparse64),intent(inout)::sparse
!   integer(kind=int32),intent(in),optional::minsizenode
!  end subroutine
!  !**GET LDLt DECOMPOSITION
!  module subroutine getldlt_crs64(sparse,minsizenode)
!   class(crssparse64),intent(inout)::sparse
!   integer(kind=int32),intent(in),optional::minsizenode
!  end subroutine
!#endif
!#if (_METIS==1)
!  !**GET ORDER
!  module function getordering_crs64(sparse&
!                            ,ctype,iptype,rtype,compress,ccorder&
!                            ,pfactor,nseps,bglvl&
!                            ) result(perm)
!   class(crssparse64),intent(in)::sparse
!   integer(kind=int32),intent(in),optional::ctype,iptype,rtype,compress,ccorder,pfactor,nseps,bglvl
!   integer(kind=int32),allocatable::perm(:)
!  end function
!#endif
!#if (_SPAINV==1)
!  !**GET INCOMPLETE CHOLESKY FACTOR
!  module subroutine getichol_crs64(sparse,minsizenode)
!   class(crssparse64),intent(inout)::sparse
!   integer(kind=int32),intent(in),optional::minsizenode
!  end subroutine
!  !**GET SPARSE INVERSE
!  module subroutine getspainv_crs64(sparse,minsizenode)
!   class(crssparse64),intent(inout)::sparse
!   integer(kind=int32),intent(in),optional::minsizenode
!  end subroutine
!#endif
  pure module function isdecomposed_crs64(sparse) result(ll)
   class(crssparse64),intent(in)::sparse
   logical::ll
  end function
  pure module subroutine setdecomposed_crs64(sparse,ll)
   class(crssparse64),intent(inout)::sparse
   logical,intent(in)::ll
  end subroutine
#if (_PARDISO==1)
  !**RESET PARDISO MEMORY
  module subroutine reset_pardiso_memory_crs64(sparse)
   !sparse*x=y
   class(crssparse64),intent(inout)::sparse
  end subroutine
#endif
  !**PRINT
  module subroutine print_crs64(sparse,lint,output)
   class(crssparse64),intent(in)::sparse
   integer(kind=int32),intent(in),optional::output
   logical,intent(in),optional::lint
  end subroutine
  module subroutine print_idx_crs64(sparse,lidx,lint,output)
   class(crssparse64),intent(in)::sparse
   logical,intent(in)::lidx(:)
   integer(kind=int32),intent(in),optional::output
   logical,intent(in),optional::lint
  end subroutine
  module subroutine printsquare_crs64(sparse,output)
   class(crssparse64),intent(inout)::sparse
   integer(kind=int32),intent(in),optional::output
  end subroutine
  !**SAVE
  module subroutine save_crs64(sparse,namefile)
   class(crssparse64),intent(in)::sparse
   character(len=*),intent(in)::namefile
  end subroutine
  !**SCALE ALL ENTRIES
  module subroutine scale_crs64(sparse,val)
   class(crssparse64),intent(inout)::sparse
   real(kind=wp),intent(in)::val
  end subroutine
  !**SET ELEMENTS
  module subroutine set_crs64(sparse,row,col,val,error)
   !add a value only to an existing one
   class(crssparse64),intent(inout)::sparse
   integer(kind=int32),intent(in)::row,col
   integer(kind=int32),intent(out),optional::error
   real(kind=wp),intent(in)::val
  end subroutine
  !**SOLVE
  module subroutine solve_crs64_vector(sparse,x,y,msglvl)
   !sparse*x=y
   class(crssparse64),intent(inout)::sparse
   real(kind=wp),intent(out),contiguous::x(:)
   real(kind=wp),intent(inout),contiguous::y(:)
   integer(kind=int64),intent(in),optional::msglvl
  end subroutine
  module subroutine solve_crs64_array(sparse,x,y,msglvl)
   !sparse*x=y
   class(crssparse64),intent(inout)::sparse
   real(kind=wp),intent(out),contiguous::x(:,:)
   real(kind=wp),intent(inout),contiguous::y(:,:)
   integer(kind=int64),intent(in),optional::msglvl
  end subroutine
!  !**SOLVE WITH A TRIANGULAR FACTOR
!  module subroutine isolve_crs64(sparse,x,y)
!   !sparse*x=y
!   class(crssparse64),intent(in)::sparse
!   real(kind=wp),intent(out)::x(:)
!   real(kind=wp),intent(in)::y(:)
!  end subroutine
!  !**SOLVE WITH LDLt DECOMPOSITION
!  module subroutine solveldlt_s_crs64(sparse,x,y)
!   !sparse*x=y
!   class(crssparse64),intent(in)::sparse
!   real(kind=wp),intent(out)::x(:)
!   real(kind=wp),intent(in)::y(:)
!  end subroutine
  !**SORT ARRAY
  module subroutine sort_crs64(sparse)
   ! sort vectors ja and a by increasing order
   class(crssparse64),intent(inout)::sparse
  end subroutine
 end interface

 !> @brief Constructor; e.g., mat=crssparse64(dim1,#elements,[dim2],[upper_storage],[output_unit])
 !! OR mat=crssparse64('file',[output_unit])
 interface crssparse64
  module procedure constructor_crs64,load_crs64
 end interface

!!!!!!!!!!!!!!!!!!!!!!!LINKED LIST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
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
  !> @brief Initializes llsparse
  procedure,public::init=>constructor_sub_ll
  !> @brief Multiplication with a vector
  procedure::multbyv=>multgenv_ll
  !> @brief Multiplication with a matrix
  procedure::multbym=>multgenm_ll
  !> @brief Returns the number of non-zero elements
  procedure,public::nonzero=>totalnumberofelements_ll
  !> @brief Prints the sparse matrix to the output sparse\%unlog
  procedure,public::print_all=>print_ll
  !> @brief Prints the selected entries of the sparse matrix to the output mat\%unlog
  procedure,public::print_idx=>print_idx_ll
  !> @brief Prints the sparse matrix in a rectangular/square format to the default output
  procedure,public::printsquare=>printsquare_ll
  !> @brief Scales all entries of mat by real scalar val; e.g., call mat\%scale(val)
  procedure,public::scale=>scale_ll
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

 interface
  !**CONSTRUCTOR
  module subroutine constructor_sub_ll(sparse,m,n,lupper,unlog)
   class(llsparse),intent(out)::sparse
   integer(kind=int32),intent(in)::m
   integer(kind=int32),intent(in),optional::n,unlog
   logical,intent(in),optional::lupper
  end subroutine
  !**DESTROY
  module elemental subroutine destroy_scal_ptrnode(pnode)
   type(ptrnode), intent(inout)::pnode
  end subroutine
  module elemental subroutine destroy_ll(sparse)
   class(llsparse),intent(inout)::sparse
  end subroutine
  !**DIAGONAL ELEMENTS
  !EQUALITIES
  module subroutine equal_node(nodeout,nodein)
   class(node),intent(out)::nodeout
   class(node),intent(in)::nodein
  end subroutine
  !**ADD ELEMENTS
  module subroutine addtohead_ptrnode(pnode,col,val)
   type(ptrnode),intent(inout),pointer::pnode
   integer(kind=int32),intent(in)::col
   real(kind=wp),intent(in)::val
  end subroutine
  module subroutine addtohead_ll(sparse,row,col,val)
   class(llsparse),intent(inout),target::sparse
   integer(kind=int32),intent(in)::row,col
   real(kind=wp),intent(in)::val
  end subroutine
  module subroutine addinorder_ll(sparse,row,col,val)
   class(llsparse),intent(inout),target::sparse
   integer(kind=int32),intent(in)::row,col
   real(kind=wp),intent(in)::val
  end subroutine
  module subroutine addtotail_ll(sparse,row,col,val)
   class(llsparse),intent(inout),target::sparse
   integer(kind=int32),intent(in)::row,col
   real(kind=wp),intent(in)::val
  end subroutine
  !**GET ELEMENTS
  pure module function get_ll(sparse,row,col) result(val)
   class(llsparse),intent(in)::sparse
   integer(kind=int32),intent(in)::row,col
   real(kind=wp)::val
  end function
  !**EXTERNAL
  !**LOAD
  !**MULTIPLICATIONS
  module subroutine multgenv_ll(sparse,alpha,trans,x,val,y)
   !Computes y=val*y+alpha*sparse(tranposition)*x
   class(llsparse),intent(in)::sparse
   real(kind=wp),intent(in)::val,alpha
   real(kind=wp),intent(in)::x(:)
   real(kind=wp),intent(out)::y(:)
   character(len=1),intent(in)::trans
  end subroutine
  module subroutine multgenm_ll(sparse,alpha,trans,x,val,y)
   !Computes y=val*y+alpha*sparse(tranposition)*x
   class(llsparse),intent(in)::sparse
   real(kind=wp),intent(in)::val,alpha
   real(kind=wp),intent(in)::x(:,:)
   real(kind=wp),intent(out)::y(:,:)
   character(len=1),intent(in)::trans
  end subroutine
  !**NUMBER OF ELEMENTS
  module function totalnumberofelements_ptrnode(pnode) result(nel)
   class(ptrnode),intent(in),target::pnode
   integer(kind=int64)::nel
  end function
  module function totalnumberofelements_ll(sparse) result(nel)
   class(llsparse),intent(in)::sparse
   integer(kind=int64)::nel
  end function
  !**SAVE
  !**SCALE ALL ENTRIES
  module subroutine scale_ll(sparse,val)
   class(llsparse),intent(inout)::sparse
   real(kind=wp),intent(in)::val
  end subroutine
  !**SET ELEMENTS
  !**SOLVE
  !**SORT ARRAY
  !**SUBMATRIX
  !**PRINT
  module subroutine print_ll(sparse,lint,output)
   class(llsparse),intent(in)::sparse
   integer(kind=int32),intent(in),optional::output
   logical,intent(in),optional::lint
  end subroutine
  module subroutine print_idx_ll(sparse,lidx,lint,output)
   class(llsparse),intent(in)::sparse
   logical,intent(in)::lidx(:)
   integer(kind=int32),intent(in),optional::output
   logical,intent(in),optional::lint
  end subroutine
  module subroutine printsquare_ll(sparse,output)
   class(llsparse),intent(inout)::sparse
   integer(kind=int32),intent(in),optional::output
  end subroutine
 end interface

 !> @brief Constructor; e.g., mat=llsparse(dim1,[dim2],[upper_storage],[output_unit])
 interface llsparse
  module procedure constructor_ll
 end interface


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!METIS GRAPH DATA STRUCTURE!!!!!!!!!!!!!!!!!!!!!!!!!aaa
 type::metisgraph
  private
  integer(kind=int32)::unlog
  integer(kind=int32)::nvertices,medges
  integer(kind=int32),allocatable::xadj(:),adjncy(:)
  type(c_ptr)::vwgt
  type(c_ptr)::adjwgt
  contains
  private
  procedure,public::init=>constructor_sub_metisgraph
  procedure,public::destroy=>destroy_metisgraph
  procedure,public::getmem=>getmem_metisgraph
  final::deallocate_scal_metisgraph
 end type

 interface
  !**CONSTRUCTOR
  module subroutine constructor_sub_metisgraph(metis,n,m,unlog)
   class(metisgraph),intent(out)::metis
   integer(kind=int32),intent(in)::n,m
   integer(kind=int32),intent(in),optional::unlog
  end subroutine
  !** GET MEMORY
  pure module function getmem_metisgraph(metis) result(getmem)
   class(metisgraph),intent(in)::metis
   integer(kind=int64)::getmem
  end function
 end interface

 interface metisgraph
  module procedure constructor_metisgraph
 end interface

!!!!!!!!!!!!!GENERAL INTERFACES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
 !> @brief Converts sparse matrices from one format to another one; e.g., crsmat=coomat
 interface assignment(=)
  module procedure convertfromlltocoo,convertfromlltocrs&
                  ,convertfromcootocrs,convertfromcootoll&
                  ,convertfromcootocrs64&
                  ,convertfromcrstometisgraph&
                  ,convertfromcrstocoo,convertfromcrstoll&
                  ,convertfromcrs64tocoo
 end interface

 interface deallocate_sparse
  module procedure deallocate_scal_coo &
                   , deallocate_scal_crs &
                   , deallocate_scal_crs64 &
                   , deallocate_scal_ll &
                   , deallocate_scal_metisgraph
 end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!COO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
!**DIAGONAL ELEMENTS
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

!FINAL
impure elemental subroutine deallocate_scal_coo(sparse)
 type(coosparse),intent(inout)::sparse

 call sparse%destroy()

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!CRS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
!**DIAGONAL ELEMENTS
function diag_mat_crs(sparse,noff) result(diagsparse)
 class(crssparse),intent(in)::sparse
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
  error stop
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
 if(sparse%lupperstorage.eqv.lupperstorage.or.(sparse%lupperstorage.and..not.lincludediag))then
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
!  do i=startdim2,enddim2
!   do j=sparse%ia(i),sparse%ia(i+1)-1
!    k=sparse%ja(j)
!    if(i.ne.k.and.k.ge.startdim1.and.k.le.enddim1)then
!     nel=nel+1
!    endif
!   enddo
!  enddo
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
 if((sparse%lupperstorage.eqv.lupperstorage).or.(sparse%lupperstorage.and..not.lincludediag))then
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
   !subcoo=coosparse(enddim1-startdim1+1,enddim2-startdim2+1,int(nel,int64),lupperstorage,unlog)
   call subcoo%init(enddim1-startdim1+1,enddim2-startdim2+1,int(nel,int64),lupperstorage,unlog)
  else
   !subcoo=coosparse(enddim1-startdim1+1,enddim2-startdim2+1,int(nel,int64),lupperstorage)
   call subcoo%init(enddim1-startdim1+1,enddim2-startdim2+1,int(nel,int64),lupperstorage)
  endif
  do i=1,sparse%dim1
   do j=sparse%ia(i),sparse%ia(i+1)-1
    k=sparse%ja(j)
    if(k.ge.startdim2.and.k.le.enddim2)then
     call subcoo%add(i-startdim1+1,k-startdim2+1,sparse%a(j))
    endif
!    if(i.ne.k)then
!     if((k.ge.startdim1.and.k.le.enddim1).and.(i.ge.startdim2.and.i.le.enddim2))then
!      call subcoo%add(k-startdim1+1,i-startdim2+1,sparse%a(j))
!     endif
!    endif
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

pure subroutine submatrix_dense_crs(sparse,indx,dense,lupper,unlog)
 !Not programmed efficiently, but it should do the job
 class(crssparse),intent(in)::sparse
 integer(kind=int32),intent(in)::indx(:)
 real(kind=real64),intent(out),allocatable::dense(:,:)
 integer(kind=int32),intent(in),optional::unlog
 logical,intent(in),optional::lupper

 integer(kind=int32)::i,ii,j,k,nel
 integer(kind=int32) :: pos(1)
 logical::lupper_

 nel = size(indx)

 lupper_ = sparse%lupperstorage
 if(present(lupper))lupper_ = lupper


 allocate(dense(nel, nel), source = 0._real64)

 do ii = 1, size(indx)
  i = indx(ii)
  do j = sparse%ia(i), sparse%ia(i+1) - 1
   k = sparse%ja(j)
   pos = findloc(indx, value = k)
   if(pos(1).gt.0)then
    dense(ii, pos(1)) = sparse%a(j)
    if(.not.lupper_ .and. sparse%lupperstorage .and. sparse%lsymmetric)then
     dense(pos(1), ii) = sparse%a(j)
    endif
   endif
  enddo
 enddo

end subroutine

!FINAL
impure elemental subroutine deallocate_scal_crs(sparse)
 type(crssparse),intent(inout)::sparse

 call sparse%destroy()

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!CRS64!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
!**LOAD
function load_crs64(namefile,unlog)  result(sparse)
 type(crssparse64)::sparse
 integer(kind=int32),intent(in),optional::unlog
 character(len=*),intent(in)::namefile

 integer(kind=int32)::un,dim1,dim2
 integer(kind=int64)::nonzero
 logical::lupperstorage

 open(newunit=un,file=namefile,action='read',status='old',access='stream')!,buffered='yes')
 read(un)dim1
 if(dim1.ne.typecrs64)then
  write(*,'(a)')' ERROR: the proposed file is not a CRS64 file'
  error stop
 endif
 read(un)dim1            !int32
 read(un)dim2            !int32
 read(un)nonzero         !int64
 read(un)lupperstorage   !logical

 if(present(unlog))then
  sparse=crssparse64(dim1,nonzero,dim2,lupperstorage,unlog)
 else
  sparse=crssparse64(dim1,nonzero,dim2,lupperstorage)
 endif

 read(un)sparse%ia              !int64
 read(un)sparse%ja              !int64
 read(un)sparse%a               !wp

 close(un)

end function


!FINAL
impure elemental subroutine deallocate_scal_crs64(sparse)
 type(crssparse64),intent(inout)::sparse

 call sparse%destroy()

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!LINKED LIST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
!**CONSTRUCTOR
function constructor_ll(m,n,lupper,unlog) result(sparse)
 type(llsparse)::sparse
 integer(kind=int32),intent(in)::m
 integer(kind=int32),intent(in),optional::n,unlog
 logical,intent(in),optional::lupper

 call sparse%initialize('LINKED LIST',m,m)

 if(present(n))sparse%dim2=n
 if(present(lupper))sparse%lupperstorage=lupper
 if(present(unlog))sparse%unlog=unlog

 sparse%lsymmetric=.false.

 allocate(sparse%heads(sparse%dim1))

end function

!FINAL
impure elemental subroutine deallocate_scal_ll(sparse)
 type(llsparse),intent(inout)::sparse

 call sparse%destroy()

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!METIS GRAPH DATA STRUCTURE!!!!!!!!!!!!!!!!!!!!!!!!!aaa
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
elemental subroutine destroy_metisgraph(metis)
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
elemental subroutine deallocate_scal_metisgraph(metis)
 type(metisgraph),intent(inout)::metis

 call metis%destroy()

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!OTHER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaa
!CONVERSIONS
subroutine convertfromlltocoo(othersparse,sparse)
 type(coosparse),intent(out)::othersparse
 type(llsparse),intent(in),target::sparse

 integer(kind=int32)::i
 type(ptrnode),pointer::cursor

 othersparse=coosparse(sparse%dim1,sparse%dim2,sparse%nonzero(),sparse%lupperstorage)

 call othersparse%setsymmetric(sparse%lsymmetric)

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

 integer(kind=int32)::i,nel,col
 integer(kind=int32),allocatable::rowpos(:)
 type(ptrnode),pointer::cursor
 logical::lsquare

 if(sparse%nonzero().ge.2_int64**31)then
  write(sparse%unlog,'(a)')' ERROR: impossible conversion due a too large number of non-zero elements'
  error stop
 endif

 !Condition: all rows contain at least one element (diagonal element if square or one dummy entry in the last column if needed)

 !Number of elements=number of rows+number off-diagonal elements
 !from sparse
 lsquare=sparse%issquare()

 allocate(rowpos(sparse%dim1))
 rowpos=0
 if(lsquare)then
  do i=1,sparse%dim1
   cursor=>sparse%heads(i)
   do while(associated(cursor%p))
    if(i.ne.cursor%p%col)rowpos(i)=rowpos(i)+1
    cursor=>cursor%p%next
   enddo
  enddo
 else
  do i=1,sparse%dim1
   cursor=>sparse%heads(i)
   do while(associated(cursor%p))
    if(cursor%p%col.ne.sparse%dim2)rowpos(i)=rowpos(i)+1
    cursor=>cursor%p%next
   enddo
  enddo
 endif

 nel=sparse%dim1+sum(rowpos)

 othersparse=crssparse(sparse%dim1,nel,sparse%dim2,sparse%lupperstorage)

 call othersparse%setsymmetric(sparse%lsymmetric)

 if(allocated(sparse%perm))allocate(othersparse%perm,source=sparse%perm)

 !if(othersparse%ia(othersparse%dim1+1).ne.0)othersparse%ia=0

 !determine the number of non-zero off-diagonal elements per row
 othersparse%ia(2:othersparse%dim1+1)=rowpos

 !accumulate the number of elements and add diagonal elements
 othersparse%ia(1)=1
 if(lsquare)then
  do i=1,sparse%dim1
   othersparse%ia(i+1)=othersparse%ia(i+1)+1+othersparse%ia(i)
   othersparse%ja(othersparse%ia(i))=i
  enddo
 else
  do i=1,sparse%dim1
   othersparse%ia(i+1)=othersparse%ia(i+1)+1+othersparse%ia(i)
   othersparse%ja(othersparse%ia(i))=sparse%dim2
  enddo
 endif

 !add the non-zero elements to crs (othersparse)
 !allocate(rowpos(othersparse%dim1))
 rowpos=othersparse%ia(1:othersparse%dim1)
 if(lsquare)then
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
 else
  do i=1,sparse%dim1
   cursor=>sparse%heads(i)
   do while(associated(cursor%p))
    col=cursor%p%col
    if(col.eq.sparse%dim2)then
     othersparse%a(othersparse%ia(i))=cursor%p%val
    else
     rowpos(i)=rowpos(i)+1
     othersparse%ja(rowpos(i))=col
     othersparse%a(rowpos(i))=cursor%p%val
    endif
    cursor=>cursor%p%next
   enddo
  enddo
 endif

 deallocate(rowpos)

! call othersparse%print()

end subroutine

subroutine convertfromcootocrs(othersparse,sparse)
 type(crssparse),intent(out)::othersparse
 type(coosparse),intent(in)::sparse

 integer(kind=int32)::i,nel,row,col
 integer(kind=int32),allocatable::rowpos(:)
 integer(kind=int64)::i8
 logical::lsquare

 if(sparse%nonzero().ge.2_int64**31)then
  write(sparse%unlog,'(a)')' ERROR: impossible conversion due to a too large number of non-zero elements'
  error stop
 endif

 !Condition: all rows contain at least one element (diagonal element if square or one dummy entry in the last column if needed)

 !Number of elements=number of rows+number off-diagonal elements
 !from sparse
 lsquare=sparse%issquare()

 allocate(rowpos(sparse%dim1))
 rowpos=0

 if(lsquare)then
  do i8=1_int64,sparse%nel
   row=sparse%ij(1,i8)
   if(row.ne.0.and.row.ne.sparse%ij(2,i8))then
    rowpos(row)=rowpos(row)+1
   endif
  enddo
 else
  do i8=1_int64,sparse%nel
   row=sparse%ij(1,i8)
   if(row.ne.0.and.sparse%ij(2,i8).ne.sparse%dim2)then
    rowpos(row)=rowpos(row)+1
   endif
  enddo
 endif

 nel=sparse%dim1+sum(rowpos)

 !othersparse=crssparse(sparse%dim1,nel,sparse%dim2,sparse%lupperstorage)
 call othersparse%init(sparse%dim1,nel,sparse%dim2,sparse%lupperstorage)

 call othersparse%setsymmetric(sparse%lsymmetric)

 if(allocated(sparse%perm))allocate(othersparse%perm,source=sparse%perm)

 !if(othersparse%ia(othersparse%dim1+1).ne.0)othersparse%ia=0

 !determine the number of non-zero off-diagonal elements per row
 othersparse%ia(2:othersparse%dim1+1)=rowpos

 !accumulate the number of elements and add diagonal elements
 othersparse%ia(1)=1
 if(lsquare)then
  do i=1,sparse%dim1
   othersparse%ia(i+1)=othersparse%ia(i+1)+1+othersparse%ia(i)
   othersparse%ja(othersparse%ia(i))=i   !set diagonal element for the case it would not be present
  enddo
 else
  do i=1,sparse%dim1
   othersparse%ia(i+1)=othersparse%ia(i+1)+1+othersparse%ia(i)
   othersparse%ja(othersparse%ia(i))=sparse%dim2   !set last element for the case it would not be present
  enddo
 endif

 !add the non-zero elements to crs (othersparse)
 !allocate(rowpos(othersparse%dim1))
 rowpos=othersparse%ia(1:othersparse%dim1)
 if(lsquare)then
  do i8=1_int64,sparse%nel
   row=sparse%ij(1,i8)
   if(row.gt.0)then
    col=sparse%ij(2,i8)
    if(col.eq.row)then
     othersparse%a(othersparse%ia(row))=sparse%a(i8)
    else
     rowpos(row)=rowpos(row)+1
     othersparse%ja(rowpos(row))=col
     othersparse%a(rowpos(row))=sparse%a(i8)
    endif
   endif
  enddo
 else
  do i8=1_int64,sparse%nel
   row=sparse%ij(1,i8)
   if(row.gt.0)then
    col=sparse%ij(2,i8)
    if(col.eq.sparse%dim2)then
     othersparse%a(othersparse%ia(row))=sparse%a(i8)
    else
     rowpos(row)=rowpos(row)+1
     othersparse%ja(rowpos(row))=col
     othersparse%a(rowpos(row))=sparse%a(i8)
    endif
   endif
  enddo
 endif

 deallocate(rowpos)

! call othersparse%print()

end subroutine

subroutine convertfromcootocrs64(othersparse,sparse)
 type(crssparse64),intent(out)::othersparse
 type(coosparse),intent(in)::sparse

 integer(kind=int64)::i,nel,row,col
 integer(kind=int64),allocatable::rowpos(:)
 integer(kind=int64)::i8
 logical::lsquare


 !Condition: all rows contain at least one element (diagonal element if square or one dummy entry in the last column if needed)

 !Number of elements=number of rows+number off-diagonal elements
 !from sparse
 lsquare=sparse%issquare()

 allocate(rowpos(sparse%dim1))
 rowpos=0

 if(lsquare)then
  do i8=1_int64,sparse%nel
   row=sparse%ij(1,i8)
   if(row.ne.0.and.row.ne.sparse%ij(2,i8))then
    rowpos(row)=rowpos(row)+1
   endif
  enddo
 else
  do i8=1_int64,sparse%nel
   row=sparse%ij(1,i8)
   if(row.ne.0.and.sparse%ij(2,i8).ne.sparse%dim2)then
    rowpos(row)=rowpos(row)+1
   endif
  enddo
 endif

 nel=sparse%dim1+sum(rowpos)

 !othersparse=crssparse(sparse%dim1,nel,sparse%dim2,sparse%lupperstorage)
 call othersparse%init(sparse%dim1,nel,sparse%dim2,sparse%lupperstorage)

 call othersparse%setsymmetric(sparse%lsymmetric)

 if(allocated(sparse%perm))allocate(othersparse%perm,source=sparse%perm)
 if(allocated(sparse%perm64))allocate(othersparse%perm64,source=sparse%perm64)

 !if(othersparse%ia(othersparse%dim1+1).ne.0)othersparse%ia=0

 !determine the number of non-zero off-diagonal elements per row
 othersparse%ia(2:othersparse%dim1+1)=rowpos

 !accumulate the number of elements and add diagonal elements
 othersparse%ia(1)=1
 if(lsquare)then
  do i=1,sparse%dim1
   othersparse%ia(i+1)=othersparse%ia(i+1)+1+othersparse%ia(i)
   othersparse%ja(othersparse%ia(i))=i   !set diagonal element for the case it would not be present
  enddo
 else
  do i=1,sparse%dim1
   othersparse%ia(i+1)=othersparse%ia(i+1)+1+othersparse%ia(i)
   othersparse%ja(othersparse%ia(i))=sparse%dim2   !set last element for the case it would not be present
  enddo
 endif

 !add the non-zero elements to crs (othersparse)
 !allocate(rowpos(othersparse%dim1))
 rowpos=othersparse%ia(1:othersparse%dim1)
 if(lsquare)then
  do i8=1_int64,sparse%nel
   row=sparse%ij(1,i8)
   if(row.gt.0)then
    col=sparse%ij(2,i8)
    if(col.eq.row)then
     othersparse%a(othersparse%ia(row))=sparse%a(i8)
    else
     rowpos(row)=rowpos(row)+1
     othersparse%ja(rowpos(row))=col
     othersparse%a(rowpos(row))=sparse%a(i8)
    endif
   endif
  enddo
 else
  do i8=1_int64,sparse%nel
   row=sparse%ij(1,i8)
   if(row.gt.0)then
    col=sparse%ij(2,i8)
    if(col.eq.sparse%dim2)then
     othersparse%a(othersparse%ia(row))=sparse%a(i8)
    else
     rowpos(row)=rowpos(row)+1
     othersparse%ja(rowpos(row))=col
     othersparse%a(rowpos(row))=sparse%a(i8)
    endif
   endif
  enddo
 endif

 deallocate(rowpos)

! call othersparse%print()

end subroutine

subroutine convertfromcootoll(othersparse,sparse)
 type(llsparse),intent(out)::othersparse
 type(coosparse),intent(in)::sparse

 integer(kind=int32)::row
 integer(kind=int64)::i8

 othersparse=llsparse(sparse%dim1,sparse%dim2,sparse%lupperstorage)

 call othersparse%setsymmetric(sparse%lsymmetric)

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

 call othersparse%setsymmetric(sparse%lsymmetric)

 if(allocated(sparse%perm))allocate(othersparse%perm,source=sparse%perm)

 do i=1,sparse%dim1
  do j=sparse%ia(i),sparse%ia(i+1)-1
   call othersparse%add(i,sparse%ja(j),sparse%a(j))
  enddo
 enddo

! call othersparse%print()

end subroutine

subroutine convertfromcrs64tocoo(othersparse,sparse)
 type(coosparse),intent(out)::othersparse
 type(crssparse64),intent(in)::sparse

 integer(kind=int32)::i
 integer(kind=int64)::j

 othersparse=coosparse(sparse%dim1,sparse%dim2,sparse%nonzero(),sparse%lupperstorage)

 call othersparse%setsymmetric(sparse%lsymmetric)

 if(allocated(sparse%perm))allocate(othersparse%perm,source=sparse%perm)
 if(allocated(sparse%perm64))allocate(othersparse%perm64,source=sparse%perm64)

 do i=1,sparse%dim1
  do j=sparse%ia(i),sparse%ia(i+1)-1
   call othersparse%add(i,int(sparse%ja(j),int32),sparse%a(j))
  enddo
 enddo

! call othersparse%print()

end subroutine

subroutine convertfromcrstoll(othersparse,sparse)
 type(llsparse),intent(out)::othersparse
 type(crssparse),intent(in)::sparse

 integer(kind=int32)::i,j

 othersparse=llsparse(sparse%dim1,sparse%dim2,sparse%lupperstorage)

 call othersparse%setsymmetric(sparse%lsymmetric)

 if(allocated(sparse%perm))allocate(othersparse%perm,source=sparse%perm)

 do i=1,sparse%dim1
  do j=sparse%ia(i),sparse%ia(i+1)-1
   call othersparse%add(i,sparse%ja(j),sparse%a(j))
  enddo
 enddo

! call othersparse%print()

end subroutine

subroutine convertfromcrstometisgraph(metis,sparse)
 type(metisgraph),intent(out)::metis
 type(crssparse),intent(in)::sparse

 integer(kind=int32)::i,j,k
 integer(kind=int32)::n,nvertices,medges
 integer(kind=int32),allocatable::rowpos(:)

 if(.not.sparse%lupperstorage.or..not.sparse%issquare())then
  write(sparse%unlog,'(a)')' ERROR: the CRS matrix must be square and upper triangular stored'
  error stop
 endif

 n=sparse%getdim(1)
 nvertices=n

 medges=sparse%nonzero()-n  !number of off-diagonal elements

 !metis=metisgraph(nvertices,medges,unlog=sparse%unlog)
 call metis%init(nvertices,medges,unlog=sparse%unlog)

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
end module

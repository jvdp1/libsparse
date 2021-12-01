!> Module containing a derived type for MKL Pardiso variables

!> @todo Not happy with the phase variable and how it is handled in modsparse

module modvariablepardiso
 use modsparse_mkl, only: mkl_pardiso_handle
#if (_DP==0)
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real32
#else
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real64
#endif
 !$ use omp_lib
 implicit none
 private
 public::checkpardiso
 public::pardiso_variable

! type::pardiso_pt
!  integer(kind=8)::pt
! end type

 type::pardiso_variable
  !Default Pardiso variables
  integer(kind=int32)::idum(1)
  real(kind=wp)::ddum(1)
  integer(kind=int32)::maxfct
  integer(kind=int32)::mnum
  integer(kind=int32)::mtype
  integer(kind=int32)::msglvl
  integer(kind=int32)::phase
  integer(kind=int32)::solver
  integer(kind=int32),allocatable::iparm(:)
  type(mkl_pardiso_handle),allocatable::pt(:)
  contains
  final::reset_pardiso_variable
 end type

 interface pardiso_variable
  module procedure constructor_pardiso_variable
 end interface

contains

!**CONSTRUCTOR
function constructor_pardiso_variable(maxfct,mnum,mtype,solver,msglvl) result(this)
 type(pardiso_variable)::this
 integer(kind=int32),intent(in),optional::maxfct
 integer(kind=int32),intent(in),optional::mnum
 integer(kind=int32),intent(in),optional::mtype
 integer(kind=int32),intent(in),optional::solver
 integer(kind=int32),intent(in),optional::msglvl

 integer(kind=int32)::i

 this%phase=-9

 this%maxfct=1
 if(present(maxfct))this%maxfct=maxfct

 this%mnum=1
 if(present(mnum))this%mnum=mnum

 this%mtype=-2
 if(present(mtype))this%mtype=mtype

 this%solver=0
 if(present(solver))this%solver=solver

 this%msglvl=1
 if(present(msglvl))this%msglvl=msglvl

 allocate(this%iparm(64))
 this%iparm=0

 allocate(this%pt(64))
 do i=1,64
  this%pt(i)%dummy=0
 enddo
 
end function

!**OTHER
subroutine checkpardiso(phase,error,un)
 integer(kind=int32),intent(in)::phase,error
 integer(kind=int32),intent(in),optional::un

 integer(kind=int32)::unlog

 unlog=output_unit
 if(present(un))unlog=un

 if(error.ne.0)then
  write(unlog,'(2(a,i0))')' The following error for phase ',phase,' was detected: ',error
  stop
 endif

end subroutine

!FINAL
subroutine reset_pardiso_variable(this)
 type(pardiso_variable),intent(inout)::this

 this%maxfct=1
 this%mnum=1
 this%mtype=-2
 this%solver=0
 this%idum=0
 this%ddum=0._wp
 this%msglvl=1
 this%phase=-9
 if(allocated(this%iparm))deallocate(this%iparm)
 if(allocated(this%pt))deallocate(this%pt)

end subroutine

end module

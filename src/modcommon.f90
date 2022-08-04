!> Module containing subroutines and functions common to all other modules

module modcommon
#if (_DP==0)
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real32
#else
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real64
#endif
 implicit none
 private
 public::progress

contains

!PUBLIC
subroutine progress(i,n,un)
 integer(kind=int32),intent(in)::i,n
 integer(kind=int32),intent(in),optional::un
 
 integer(kind=int32)::unlog
 integer(kind=int32)::val
 integer(kind=int32)::step
 integer(kind=int32)::lastval=0

 unlog=output_unit
 if(present(un))unlog=un

 val=int(real(i)/real(n)*100.)

 step=10
 
 if(lastval+1.le.val/step)then
  write(unlog,'(2x,i0,"%")',advance='no')val
  lastval=val/step
 endif
 if(val.eq.100)then
  write(unlog,'(a/)',advance='yes')
  lastval=0
 endif

end subroutine

end module

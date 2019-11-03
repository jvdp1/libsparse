!> Module containing subroutines and functions common to all other modules

module modcommon
#if (_DP==0)
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real32
#else
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real64
#endif
 implicit none
 private
 public::alloc_err,progress

contains

!PUBLIC
subroutine alloc_err(line,namefile)
 integer(kind=int32),optional::line
 character(len=*),optional::namefile

 integer(kind=int32)::i
 character(len=30)::nfile

 i=-9
 if(present(line))i=line

 nfile='unknown'
 if(present(namefile))nfile=namefile
 
 write(*,'(a,i0,3a)')' ERROR (',i,',',namefile,'): failed allocation'
 error stop

end subroutine 

subroutine progress(i,n,un)
 integer(kind=int32),intent(in)::i,n
 integer(kind=int32),intent(in),optional::un
 
 integer(kind=int32)::unlog
 integer(kind=int32)::val
 integer(kind=int32)::lastval=0

 unlog=output_unit
 if(present(un))unlog=un

 val=int(real(i)/real(n)*100.)
 
 if(lastval+1.le.val/10)then
  write(unlog,'(2x,i0,"%")',advance='no')val
  lastval=val/10
 endif
 if(val.eq.100)then
  write(unlog,'(a/)',advance='yes')
  lastval=0
 endif

end subroutine

end module

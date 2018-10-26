program test5
#if (_DP==0)
 use iso_fortran_env,only:int32,int64,wp=>real32
#else
 use iso_fortran_env,only:int32,int64,wp=>real64
#endif
 use modsparse
 implicit none
 integer(kind=int32)::nrow
 integer(kind=int32)::row
 integer(kind=int32)::col
 integer(kind=int32)::iunit, istat
 integer(kind=int32),allocatable::iperm(:)
 real(kind=wp)::val
 logical::lup=.false.
 type(coosparse)::coo
 type(crssparse)::crs

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 !COO UPPER
 open(newunit=iunit,file='crsinput.ascii',status='old',action='read')
 read(iunit,*) nrow

 coo=coosparse(nrow,nel=4_int64,lupper=.true.)

 do
  read(iunit,*,iostat=istat) row,col,val
  if(istat.ne.0)exit
  call coo%add(row,col,val)
 end do
 close(iunit)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 !CSR UPPER
 crs=coo
 call crs%printstats()
 call crs%print(lint=.true.)

 iperm=crs%getordering()
 
 print*,'iperm',iperm
 



end program

program test1
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
 real(kind=wp)::val
 logical::lup=.false.
 type(coosparse)::coo
 type(crssparse)::crs,crs1

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 open(newunit=iunit,file='crsinput.ascii',status='old',action='read')
 read(iunit,*) nrow

 coo=coosparse(8,5)

 do
  read(iunit,*,iostat=istat) row,col,val
  if(istat.ne.0)exit
  call coo%add(row,col,val)
 end do
 close(iunit)

 call coo%printstats()
 call coo%print()
 call coo%printsquare()
 !call coo%print(lint=.false.)

 crs=coo
 call crs%printstats()
 call crs%print()
 call crs%printsquare()

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 open(newunit=iunit,file='crsinput.ascii',status='old',action='read')
 read(iunit,*) nrow

 coo=coosparse(8,5)

 do
  read(iunit,*,iostat=istat) row,col,val
  if(istat.ne.0)exit
  call coo%add(row+2,col+2,val)
 end do
 close(iunit)

 call coo%printstats()
 call coo%print()
 call coo%printsquare()
 !call coo%print(lint=.false.)

 crs=coo
 call crs%printstats()
 call crs%print()
 call crs%printsquare()


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 open(newunit=iunit,file='crsinput.ascii',status='old',action='read')
 read(iunit,*) nrow

 coo=coosparse(8)

 do
  read(iunit,*,iostat=istat) row,col,val
  if(istat.ne.0)exit
  call coo%add(row,col,val)
  call coo%add(col,row,val)
 end do
 close(iunit)

 call coo%printstats()
 call coo%print()
 call coo%printsquare()
 !call coo%print(lint=.false.)

 crs=coo
 call crs%printstats()
 call crs%print()
 call crs%printsquare()


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 open(newunit=iunit,file='crsinput.ascii',status='old',action='read')
 read(iunit,*) nrow

 coo=coosparse(8)

 do
  read(iunit,*,iostat=istat) row,col,val
  if(istat.ne.0)exit
  call coo%add(row+2,col+2,val)
  call coo%add(col+2,row+2,val)
 end do
 close(iunit)

 call coo%printstats()
 call coo%print()
 call coo%printsquare()
 !call coo%print(lint=.false.)

 crs=coo
 call crs%printstats()
 call crs%print()
 call crs%printsquare()


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 open(newunit=iunit,file='crsinput.ascii',status='old',action='read')
 read(iunit,*) nrow

 coo=coosparse(5,8)

 do
  read(iunit,*,iostat=istat) row,col,val
  if(istat.ne.0)exit
  call coo%add(row,col,val)
 end do
 close(iunit)

 call coo%printstats()
 call coo%print()
 call coo%printsquare()
 !call coo%print(lint=.false.)

 crs=coo
 call crs%printstats()
 call crs%print()
 call crs%printsquare()


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 open(newunit=iunit,file='crsinput.ascii',status='old',action='read')
 read(iunit,*) nrow

 coo=coosparse(5,8)

 do
  read(iunit,*,iostat=istat) row,col,val
  if(istat.ne.0)exit
  call coo%add(row+2,col+2,val)
  call coo%add(col+2,row+2,val)
 end do
 close(iunit)

 call coo%printstats()
 call coo%print()
 call coo%printsquare()
 !call coo%print(lint=.false.)

 crs=coo
 call crs%printstats()
 call crs%print()
 call crs%printsquare()




end program

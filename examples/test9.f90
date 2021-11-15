program test9
#if (_DP==0)
 use iso_fortran_env,only:int32,int64,wp=>real32
#else
 use iso_fortran_env,only:int32,int64,wp=>real64
#endif
 use modsparse
 implicit none
 integer(kind=int32)::i
 integer(kind=int32)::nrow
 integer(kind=int32)::row
 integer(kind=int32)::col
 integer(kind=int32)::iunit, istat
 real(kind=wp)::val
 real(kind=wp),allocatable::x(:),y(:),y1(:)
 real(kind=wp),allocatable::xa(:,:),ya(:,:),ya1(:,:)
 logical::lup=.false.
 type(coosparse)::coo
 type(crssparse)::crs

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 open(newunit=iunit,file='crsinput.ascii',status='old',action='read')
 read(iunit,*) nrow

 coo=coosparse(nrow)

 do
  read(iunit,*,iostat=istat) row,col,val
  if(istat.ne.0)exit
  call coo%add(row,col,val)
 end do
 close(iunit)


 allocate(x(coo%getdim(2)),y(coo%getdim(1)),y1(coo%getdim(1)))
 allocate(xa(coo%getdim(2),3),ya(coo%getdim(1),3),ya1(coo%getdim(1),3))

 x(:) = [(i,i=1,size(x,1))]
 xa = 3
 xa(:,1) = [(i,i=1,size(x,1))]

 call coo%printstats()
 call coo%print()
 call coo%printsquare()
 !call coo%print(lint=.false.)

 crs=coo
 call crs%printstats()
 call crs%print()
 call crs%printsquare()

 call coo%mult(1._wp,'n',x,2.5_wp,y)
 call crs%mult(1._wp,'n',x,2.5_wp,y1)
 print*,'sum ',sum(abs(y-y1))

 call coo%mult(1._wp,'n',xa,2.5_wp,ya)
 call crs%mult(1._wp,'n',xa,2.5_wp,ya1)
 print*,'sum ',sum(abs(ya-ya1))

 deallocate(x,xa,y,ya,y1,ya1)
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 open(newunit=iunit,file='crsinput.ascii',status='old',action='read')
 read(iunit,*) nrow

 coo=coosparse(nrow,lupper=.true.)

 do
  read(iunit,*,iostat=istat) row,col,val
  if(istat.ne.0)exit
  call coo%add(row,col,val)
 end do
 close(iunit)


 allocate(x(coo%getdim(2)),y(coo%getdim(1)))
 allocate(xa(coo%getdim(2),3),ya(coo%getdim(1),3))

 x(:) = [(i,i=1,size(x,1))]
 xa = 3
 xa(:,1) = [(i,i=1,size(x,1))]

 call coo%printstats()
 call coo%print()
 call coo%printsquare()
 !call coo%print(lint=.false.)

 crs=coo
 call crs%sort()
 call crs%printstats()
 call crs%print()
 call crs%printsquare()

 y=1.7
 y1=y
 call coo%mult(1._wp,'n',x,2.5_wp,y)
 call crs%mult(1._wp,'n',x,2.5_wp,y1)
 print*,'sum ',sum(abs(y-y1))

 ya=1.7
 ya1=ya
 call coo%mult(1._wp,'n',xa,2.5_wp,ya)
 call crs%mult(1._wp,'n',xa,2.5_wp,ya1)
 print*,'sum ',sum(abs(ya-ya1))


 deallocate(x,xa,y,ya,y1,ya1)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 open(newunit=iunit,file='crsinput.ascii',status='old',action='read')
 read(iunit,*) nrow

 coo=coosparse(nrow,lupper=.true.)

 call coo%setsymmetric()

 do
  read(iunit,*,iostat=istat) row,col,val
  if(istat.ne.0)exit
  call coo%add(row,col,val)
 end do
 close(iunit)


 allocate(x(coo%getdim(2)),y(coo%getdim(1)))
 allocate(xa(coo%getdim(2),3),ya(coo%getdim(1),3))

 x(:) = [(i,i=1,size(x,1))]
 xa = 3
 xa(:,1) = [(i,i=1,size(x,1))]

 call coo%printstats()
 call coo%print()
 call coo%printsquare()
 !call coo%print(lint=.false.)

 crs=coo
 call crs%sort()
 call crs%printstats()
 call crs%print()
 call crs%printsquare()

 y=1.7
 y1=y
 call coo%mult(1._wp,'n',x,2.5_wp,y)
 call crs%mult(1._wp,'n',x,2.5_wp,y1)
 print*,'sum ',sum(abs(y-y1))

 ya=1.7
 ya1=ya
 call coo%mult(1._wp,'n',xa,2.5_wp,ya)
 call crs%mult(1._wp,'n',xa,2.5_wp,ya1)
 print*,'sum ',sum(abs(ya-ya1))

end program

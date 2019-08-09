program test5
#if (_DP==0)
 use iso_fortran_env,only:int32,int64,wp=>real32
#else
 use iso_fortran_env,only:int32,int64,wp=>real64
#endif
 use modmetis,only:METIS_CTYPE_SHEM
 use modsparse
 implicit none
 integer(kind=int32)::nrow
 integer(kind=int32)::row
 integer(kind=int32)::col
 integer(kind=int32)::iunit, istat
 integer(kind=int32)::i,j
 integer(kind=int32),allocatable::iarray(:,:)
 integer(kind=int32),allocatable::perm(:)
 real(kind=wp)::val
 real(kind=wp),allocatable::x(:),y(:)
 logical::lup=.false.
 type(coosparse)::coo
 type(crssparse)::crs

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 !COO UPPER
 open(newunit=iunit,file='matrixija.ascii',status='old',action='read')
 !open(newunit=iunit,file='matkm.dat',status='old',action='read')
 read(iunit,*) nrow

 coo=coosparse(nrow,lupper=.true.)

 do
  read(iunit,*,iostat=istat) row,col,val
  if(istat.ne.0)exit
  call coo%add(row,col,val)
!  if(row.ne.col)call coo%add(col,row,val)
 end do
 close(iunit)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 !CSR UPPER
 crs=coo

 call crs%printstats()

 call crs%printsquare()

 call crs%setpermutation(crs%getordering(bglvl=0))
 call crs%printstats()

!!!!!!!!!!!!!!!!
 write(*,*)'Matrix'
 allocate(iarray(crs%getdim(1),crs%getdim(2)))
 iarray=0
 do i=1,crs%getdim(1)
  do j=1,crs%getdim(2)
   val=crs%get(i,j)
   if(val.ne.0.)iarray(i,j)=1
  enddo
 enddo

 do i=1,crs%getdim(1)
  write(*,'(10000(i2))')iarray(i,:)
 enddo

 write(*,*)'Permuted matrix'
 
 perm=crs%getordering(bglvl=0)
 iarray=0
 do i=1,crs%getdim(1)
  do j=1,crs%getdim(2)
   val=crs%get(perm(i),perm(j))
   if(val.ne.0.)iarray(i,j)=1
  enddo
 enddo

 do i=1,crs%getdim(1)
  write(*,'(10000(i2))')iarray(i,:)
 enddo

!!!!!!

 call crs%spainv()

 call crs%printsquare()


!!!!!!
 crs=coo
  
 call crs%sort()
 
 call crs%setpermutation(crs%getordering(bglvl=0))
 call crs%printstats()

 allocate(x(crs%getdim(1)),y(crs%getdim(1)))
 x=1.;y=1.

 call crs%solve(x,y)


end program

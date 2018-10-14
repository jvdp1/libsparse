program test1_ll
 use iso_fortran_env
 use modsparse
 implicit none
 integer(kind=int32)::nrow,i,j
 integer(kind=int32)::row
 integer(kind=int32)::col
 integer(kind=int32)::iunit, istat
 real(kind=real64)::val
 logical::lup=.false.
 type(llsparse)::ll
 type(coosparse)::coo
 type(crssparse)::crs

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 !LINKED LINK
 open(newunit=iunit,file='crsinput.ascii',status='old',action='read')
 read(iunit,*) nrow

 ll=llsparse(nrow,lupper=lup)

 do
  read(iunit,*,iostat=istat) row,col,val
  if(istat.ne.0)exit
  call ll%add(row,col,val)
 end do
 close(iunit)

 call ll%printstats()

 call ll%print()
 !call coo%print(lint=.false.)
 call ll%printtofile('ll.dat')

 print*,'size ll: ',ll%nonzero()

 do i=1,ll%getdim(1)
  do j=1,ll%getdim(2)
   print*,i,j,ll%get(i,j)
  enddo
 enddo

end program

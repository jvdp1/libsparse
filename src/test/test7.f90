program test7
#if (_DP==0)
 use iso_fortran_env,only:int32,int64,wp=>real32
#else
 use iso_fortran_env,only:int32,int64,wp=>real64
#endif
#if (_METIS==1)
 use modmetis,only:METIS_CTYPE_SHEM
#endif
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
 real(kind=wp),allocatable::xx(:)
 logical::lup=.false.
 type(coosparse)::coo
 type(crssparse)::crs
 type(crssparse)::crs1

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 !COO UPPER
 open(newunit=iunit,file='crsinput.ascii',status='old',action='read')
 read(iunit,*) nrow

 coo=coosparse(5,lupper=.true.)

 do
  read(iunit,*,iostat=istat) row,col,val
  if(istat.ne.0)exit
  call coo%add(row,col,val)
 end do
 close(iunit)

 call coo%add(5,5,1._wp)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 !CSR UPPER
 crs=coo

 call crs%printsquare()

 !perm=crs%getordering()
#if (_METIS==1)
 perm=crs%getordering(compress=1,ctype=METIS_CTYPE_SHEM)
#else
 perm=(/(i,i=1,crs%getdim(1))/)
#endif

 call crs%printstats()

 call crs%setpermutation(perm)
 call crs%printstats()

 allocate(x(crs%getdim(1)),y(crs%getdim(1)))
 do i=1,crs%getdim(1)
  y(i)=i
 enddo
 print*,'y=[',y,' ]'

 x=0._wp
 call crs%solve(x,y)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 !CSR UPPER
 call coo%add(5,5,10._wp)
 call coo%add(2,5,11._wp)
 call coo%add(5,2,11._wp)
 call coo%add(1,5,11._wp)
 call coo%add(5,1,11._wp)

 do i=1,coo%getdim(1)
  call coo%add(i,i,1._wp)
 enddo

 crs1=coo

 call crs1%printsquare()

#if (_METIS==1)
 perm=crs1%getordering(compress=1,ctype=METIS_CTYPE_SHEM)
#else
 perm=(/(i,i=1,crs1%getdim(1))/)
#endif

 call crs1%printstats()

 call crs1%setpermutation(perm)
 call crs1%printstats()

 allocate(xx(crs1%getdim(1)))
 do i=1,crs1%getdim(1)
  y(i)=i
 enddo

 xx=0._wp
 call crs1%solve(xx,y)

 do i=1,crs1%getdim(1)
  write(*,*)x(i),xx(i)
 enddo

end program

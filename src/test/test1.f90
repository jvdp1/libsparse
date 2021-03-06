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
 type(llsparse)::ll
 type(coosparse)::coo
 type(crssparse)::crs,crs1

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

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 !CSR FROM LL
 crs=ll
 call crs%printstats
 call crs%print()
 
 call crs%sort()
 call crs%printtofile('crsll.dat')

 call ll%destroy()

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 !COO
 open(newunit=iunit,file='crsinput.ascii',status='old',action='read')
 read(iunit,*) nrow

 coo=coosparse(nrow,nel=4_int64,lupper=lup)

 do
  read(iunit,*,iostat=istat) row,col,val
  if(istat.ne.0)exit
  call coo%add(row,col,val)
 end do
 close(iunit)

 call coo%printstats()
 !call coo%print(lint=.false.)
 call coo%print()
 call coo%printtofile('coo.dat')

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 !CSR
 crs=coo
 call crs%printstats
 call crs%print()
 
 call crs%sort()
 call crs%printtofile('crscoo.dat')


 call coo%add(1,1,150._wp)

! crs=coo
! call crs%printstats
! !call crs%print(lint=.false.)
! print*,'*****************not sorted'
! call crs%print()
! print*,'*****************sorted'
! call crs%sort()
! call crs%print()
! call crs%printtofile('crs.dat')
! call crs%printtofile('crs_int.dat',lint=.false.)
!
!
 print*,'aaaaaaaaaaaaaaCOOaaaaaaaaaaaaaaaa'
 print*,'4 3',coo%get(4,3)
 print*,'3 4',coo%get(3,4)
 print*,'3 3',coo%get(3,3)
 print*,'1 4',coo%get(1,4)
 print*,'1 1',coo%get(1,1)
 
!print*,'aaaaaaaaaaaaaaCSRaaaaaaaaaaaaaaaa'
! print*,'4 3',crs%get(4,3)
! print*,'3 4',crs%get(3,4)
! print*,'3 3',crs%get(3,3)
! print*,'1 4',crs%get(1,4)
! print*,'1 1',crs%get(1,1)
! 
 print*,'aaaaaaaaaaaaaaCOOaaaaaaaaaaaaaaaa'
 call coo%printsquare()
!print*,'aaaaaaaaaaaaaaCOOaaaaaaaaaaaaaaaa'
! call crs%printsquare()
! call crs%printsquaretofile('crsquarse.dat')


 print*,'aaaaaaaaaaaaaaCOO modifaaaaaaaaaaaaaaaa'

 call coo%set(4,4,-1._wp)
 call coo%set(3,4,-1._wp)
 call coo%set(4,3,0._wp)


 call coo%print()

 call coo%printsquare()

 print*,'aaaaaaaaaaaaaaCRS modifaaaaaaaaaaaaaaaa'
 crs=coo
 call crs%sort()
 call crs%print()


 print*,'aaaaaaaaaaaaa print COO to stream aaaaaaaaaa'
 call coo%save('coostream.dat')
 call coo%printstats()

 call coo%destroy()

 print*,'aaaaaaaaaaaaa read COO to stream aaaaaaaaaa'
 coo=coosparse('coostream.dat')!,444)
 call coo%printstats()

 call coo%destroy()

 print*,'aaaaaaaaaaaaa add to CRS aaaaaaaaaa'
 call crs%set(4,4,100._wp,row)
 print*,'error: ',row

 call crs%set(1,4,100._wp,row)
 print*,'error: ',row
 call crs%print()

 print*,'dim1',crs%getdim(1)
 print*,'dim2',crs%getdim(2)

 print*,'aaaaaaaaaaaaaaa test load/save aaaaaaaaaaaaaaa'
 call crs%save('crsstream.dat')
 
 crs1=crssparse('crsstream.dat')
 call crs1%print()
 
 print*,'aaaaaaaaaaaaaaa test copy aaaaaaaaaaaaaaa'
 call crs%destroy()

 crs=crs1

 print*,'aaaaaaaaaaaaaaaCRSaaaaaaaaaaaaaaa'
 call crs%add(1,1,130._wp)
 call crs%printstats()
 call crs%print()

 print*,'aaaaaaaaaaaaaaaCRS1aaaaaaaaaaaaaaa'
 call crs1%printstats()
 call crs1%print()


 print*,'aaaaaaaaaaaaaaa DIAGONAL COO aaaaaaaaaaaaaaaa'
 coo=crs1

 call coo%printsquare()

 print*,'aaaa',coo%diag()

 crs=coo%diag(-1)

 call crs%printsquare()

 call coo%destroy()

 print*,'aaaaaaaaaaaaaaa DIAGONAL CRS aaaaaaaaaaaaaaaa'
 coo=crs1
 crs=coo
 call crs%printsquare()

 print*,'aaaa',crs%diag()
 crs1=crs%diag(10)
 call crs1%printsquare()

end program

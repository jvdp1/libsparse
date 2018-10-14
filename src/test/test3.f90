program test3
 use iso_fortran_env
 use modsparse
 implicit none
 integer(kind=int32)::nrow
 integer(kind=int32)::row
 integer(kind=int32)::col
 integer(kind=int32)::iunit, istat
 real(kind=real64)::val
 logical::lup=.false.
 type(coosparse)::coo
 type(crssparse)::crs,subcrs

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 !COO  FULL
 open(newunit=iunit,file='crsinput.ascii',status='old',action='read')
 read(iunit,*) nrow

 coo=coosparse(nrow,nel=4_int64,lupper=.false.)

 do
  read(iunit,*,iostat=istat) row,col,val
  if(istat.ne.0)exit
  call coo%add(row,col,val)
 end do
 close(iunit)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 !CSR FULL
 crs=coo
 call crs%printstats
 call crs%print(lint=.true.)
 
 print*,'FULL => FULL xxxxxxxxxxxxxxxxxxxxxxxxxx'
 call crs%printsquare()
 print*,'SELECT UPPER BLOCK xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 1 2 3 4'
 subcrs=crs%submatrix(1,2,3,4)
 call subcrs%printstats()
 call subcrs%printsquare()

 print*,'SELECT LOWER BLOCK xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 3 4 1 2'
 subcrs=crs%submatrix(3,4,1,2)
 call subcrs%printstats()
 call subcrs%printsquare()

 print*,'SELECT BLOCK + DIAG xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 1 2 1 3'
 subcrs=crs%submatrix(1,2,1,3)
 call subcrs%printstats()
 call subcrs%printsquare()
 
 print*,''
 print*,'FULL => UPPER xxxxxxxxxxxxxxxxxxxxxxxxxx'
 call crs%printsquare()
 print*,'SELECT UPPER BLOCK xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 1 2 3 4'
 subcrs=crs%submatrix(1,2,3,4,.true.)
 call subcrs%printstats()
 call subcrs%printsquare()

 print*,'SELECT LOWER BLOCK a xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 3 4 1 2'
 subcrs=crs%submatrix(3,4,1,2,.true.)
 call subcrs%printstats()
 call subcrs%printsquare()

 print*,'SELECT BLOCK + DIAG xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 1 2 1 3'
 subcrs=crs%submatrix(1,2,1,3,.true.)
 call subcrs%printstats()
 call subcrs%printsquare()



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

 print*,'UPPER => FULL xxxxxxxxxxxxxxxxxxxxxxxxxx'
 call crs%printsquare()
 print*,'SELECT UPPER BLOCK xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 1 2 3 4'
 subcrs=crs%submatrix(1,2,3,4)
 call subcrs%printstats()
 call subcrs%printsquare()

 print*,'SELECT LOWER BLOCK xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 3 4 1 2'
 subcrs=crs%submatrix(3,4,1,2)
 call subcrs%printstats()
 call subcrs%printsquare()

 print*,'SELECT BLOCK + DIAG xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 1 2 1 3'
 subcrs=crs%submatrix(1,2,1,3)
 call subcrs%printstats()
 call subcrs%printsquare()
 
 print*,''
 print*,'UPPER => UPPER xxxxxxxxxxxxxxxxxxxxxxxxxx'
 call crs%printsquare()
 print*,'SELECT UPPER BLOCK xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 1 2 3 4'
 subcrs=crs%submatrix(1,2,3,4,.true.)
 call subcrs%printstats()
 call subcrs%printsquare()

 print*,'SELECT LOWER BLOCK xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 3 4 1 2'
 subcrs=crs%submatrix(3,4,1,2,.true.)
 call subcrs%printstats()
 call subcrs%printsquare()

 print*,'SELECT BLOCK + DIAG xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 1 2 1 3'
 subcrs=crs%submatrix(1,2,1,3,.true.)
 call subcrs%printstats()
 call subcrs%printsquare()






end program

program test4
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
 type(coosparse)::coo,subcoo

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


 print*,'FULL => FULL xxxxxxxxxxxxxxxxxxxxxxxxxx'
 call coo%printsquare()
 print*,'SELECT UPPER BLOCK xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 1 2 3 4'
 subcoo=coo%submatrix(1,2,3,4)
 call subcoo%printstats()
 call subcoo%printsquare()

 print*,'SELECT LOWER BLOCK xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 3 4 1 2'
 subcoo=coo%submatrix(3,4,1,2)
 call subcoo%printstats()
 call subcoo%printsquare()

 print*,'SELECT BLOCK + DIAG xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 1 2 1 3'
 subcoo=coo%submatrix(1,2,1,3)
 call subcoo%printstats()
 call subcoo%printsquare()

 print*,''
 print*,'FULL => UPPER xxxxxxxxxxxxxxxxxxxxxxxxxx'
 call coo%printsquare()
 print*,'SELECT UPPER BLOCK xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 1 2 3 4'
 subcoo=coo%submatrix(1,2,3,4,.true.)
 call subcoo%printstats()
 call subcoo%printsquare()

 print*,'SELECT LOWER BLOCK a xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 3 4 1 2'
 subcoo=coo%submatrix(3,4,1,2,.true.)
 call subcoo%printstats()
 call subcoo%printsquare()

 print*,'SELECT BLOCK + DIAG xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 1 2 1 3'
 subcoo=coo%submatrix(1,2,1,3,.true.)
 call subcoo%printstats()
 call subcoo%printsquare()



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
 !COO UPPER
 print*,'***COO UPPER'
 call coo%printstats()
 call coo%print(lint=.true.)

 print*,'UPPER => FULL xxxxxxxxxxxxxxxxxxxxxxxxxx'
 call coo%printsquare()
 print*,'SELECT UPPER BLOCK xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 1 2 3 4'
 subcoo=coo%submatrix(1,2,3,4)
 call subcoo%printstats()
 call subcoo%printsquare()

 print*,'SELECT LOWER BLOCK xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 3 4 1 2'
 subcoo=coo%submatrix(3,4,1,2)
 call subcoo%printstats()
 call subcoo%printsquare()

 print*,'SELECT BLOCK + DIAG xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 1 2 1 3'
 subcoo=coo%submatrix(1,2,1,3)
 call subcoo%printstats()
 call subcoo%printsquare()
 
 print*,''
 print*,'UPPER => UPPER xxxxxxxxxxxxxxxxxxxxxxxxxx'
 call coo%printsquare()
 print*,'SELECT UPPER BLOCK xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 1 2 3 4'
 subcoo=coo%submatrix(1,2,3,4,.true.)
 call subcoo%printstats()
 call subcoo%printsquare()

 print*,'SELECT LOWER BLOCK xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 3 4 1 2'
 subcoo=coo%submatrix(3,4,1,2,.true.)
 call subcoo%printstats()
 call subcoo%printsquare()

 print*,'SELECT BLOCK + DIAG xxxxxxxxxxxxxxxxxxxxxxxxxx'
 print*,'aaaa 1 2 1 3'
 subcoo=coo%submatrix(1,2,1,3,.true.)
 call subcoo%printstats()
 call subcoo%printsquare()



end program

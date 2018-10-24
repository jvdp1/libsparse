program test2
#if (_DP==0)
 use iso_fortran_env,only:int32,int64,real64,wp=>real32
#else
 use iso_fortran_env,only:int32,int64,real64,wp=>real64
#endif
 !$ use omp_lib
 use modsparse
 implicit none
 integer(kind=int32)::nrow
 integer(kind=int32)::row,col
 integer(kind=int32)::i,j
 integer(kind=int64)::nel
 real(kind=wp)::val
 !$ real(kind=real64)::t1
 !$ real(kind=real64)::t2
 type(llsparse)::sparse
 type(coosparse)::coo
 type(crssparse)::crs

 nrow=10000
 val=1._wp
 nel=1000000_int64
! nel=67108865_int64

 sparse=llsparse(nrow,lupper=.true.)
 !$ t2=0.d0
 do i=1,nrow!,2
  do j=1,nrow!,3
   val=real(i+j,wp)
   !$ t1=omp_get_wtime()
   call sparse%addtohead(i,j,val)
   !$ t2=t2+omp_get_wtime()-t1
  enddo
 enddo
 !!$ write(*,'(/a,f0.3)')' Elapsed time ll: ',omp_get_wtime()-t1
 !$ write(*,'(/a,f0.3)')' Elapsed time coo: ',t2
 call sparse%printstats()

 !$ t1=omp_get_wtime()
 crs=sparse
 !$ write(*,'(/a,f0.3)')' Elapsed time crsll: ',omp_get_wtime()-t1
 call crs%sort()
 call crs%printtofile('crsll.dat')
 call sparse%destroy()

 coo=coosparse(nrow,nel=nel,lupper=.true.)
 !$ t2=0.d0
 !$ t1=omp_get_wtime()
 do i=1,nrow!,2
  do j=1,nrow!,3
   val=real(i+j,wp)
   !$ t1=omp_get_wtime()
   call coo%add(i,j,val)
   !$ t2=t2+omp_get_wtime()-t1
  enddo
 enddo
 !!$ write(*,'(/a,f0.3)')' Elapsed time coo: ',omp_get_wtime()-t1
 !$ write(*,'(/a,f0.3)')' Elapsed time coo: ',t2
 call coo%printstats()

 !$ t1=omp_get_wtime()
 crs=coo
 !$ write(*,'(/a,f0.3)')' Elapsed time crs: ',omp_get_wtime()-t1
 call crs%printstats()

 call crs%sort()

 call crs%printtofile('crs.dat')

! call crs%print()
!
! call crs%add(1,2,5._wp,i)
! print*,'aaaaaa',i
! call crs%add(2,2,5._wp,i)
! print*,'aaaaaa',i
! call crs%print()
end program


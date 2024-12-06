submodule (modsparse) modsparse_crs64
 use modsparse_mkl, only: pardisoinit, pardiso => pardiso_64

#if (_PARDISO==1)
 use modvariablepardiso, only: checkpardiso => checkpardiso_64, pardiso_variable_64
#endif
  !$ use omp_lib
 implicit none


contains

!**CONSTRUCTOR
module function constructor_crs64(m,nel,n,lupper,unlog) result(sparse)
 type(crssparse64)::sparse
 integer(kind=int32),intent(in)::m
 integer(kind=int64),intent(in)::nel
 integer(kind=int32),intent(in),optional::n,unlog
 logical,intent(in),optional::lupper

 call sparse%initialize('CRS64',m,m)

 if(present(n))sparse%dim2=n
 if(present(lupper))sparse%lupperstorage=lupper
 if(present(unlog))sparse%unlog=unlog

 sparse%lsymmetric=.false.

 allocate(sparse%ia(sparse%dim1+1),sparse%ja(nel),sparse%a(nel))
 sparse%ia=0
 sparse%ia(sparse%dim1+1)=-nel
 sparse%ja=0
 sparse%a=0._wp

#if (_PARDISO==1)
 sparse%lpardisofirst=.true.
#endif

end function

module subroutine constructor_sub_crs64(sparse,m,nel,n,lupper,unlog)
 class(crssparse64),intent(out)::sparse
 integer(kind=int32),intent(in)::m
 integer(kind=int64),intent(in)::nel
 integer(kind=int32),intent(in),optional::n,unlog
 logical,intent(in),optional::lupper

 call sparse%initialize('CRS64',m,m)

 if(present(n))sparse%dim2=n
 if(present(lupper))sparse%lupperstorage=lupper
 if(present(unlog))sparse%unlog=unlog

 sparse%lsymmetric=.false.

 allocate(sparse%ia(sparse%dim1+1),sparse%ja(nel),sparse%a(nel))
 sparse%ia=0
 sparse%ia(sparse%dim1+1)=-nel
 sparse%ja=0
 sparse%a=0._wp

#if (_PARDISO==1)
 sparse%lpardisofirst=.true.
#endif

end subroutine

!**DESTROY
module impure elemental subroutine destroy_crs64(sparse)
 class(crssparse64),intent(inout)::sparse

#if(_PARDISO==1)
 call sparse%resetpardiso()
#endif

 call sparse%destroy_gen_gen()

 if(allocated(sparse%ia))deallocate(sparse%ia)
 if(allocated(sparse%ja))deallocate(sparse%ja)
 if(allocated(sparse%a))deallocate(sparse%a)

end subroutine

!**ADD ELEMENTS
module subroutine add_crs64(sparse,row,col,val,error)
 !add a value only to an existing one
 class(crssparse64),intent(inout)::sparse
 integer(kind=int32),intent(in)::row,col
 !added: error=0;Not existing: error=-1;matrix not initialized: error=-10
 integer(kind=int32),intent(out),optional::error
 real(kind=wp),intent(in)::val

 integer(kind=int64)::i

 if(present(error))error=0

 if(.not.validvalue_gen(sparse,row,col))return
 if(.not.validnonzero_gen(sparse,val))return
 if(sparse%lupperstorage.and..not.uppervalue_gen(row,col))return

 if(sparse%ia(row).eq.0)then
  if(present(error))error=-10
  return
 endif

 do i=sparse%ia(row),sparse%ia(row+1)-1
  if(sparse%ja(i).eq.col)then
   sparse%a(i)=sparse%a(i)+val
   if(present(error))error=0
   return
  endif
 enddo

 if(present(error))error=-1

end subroutine

!**GET ELEMENTS
pure module function get_crs64(sparse,row,col) result(val)
 class(crssparse64),intent(in)::sparse
 integer(kind=int32),intent(in)::row,col
 real(kind=wp)::val

 integer(kind=int64)::i
 integer(kind=int32)::trow,tcol

 val=0.0_wp

 trow=row
 tcol=col
 if(sparse%lupperstorage.and.sparse%lsymmetric.and.row.gt.col)then
  !swap row-col
  trow=col
  tcol=row
 endif

 do i=sparse%ia(trow),sparse%ia(trow+1)-1
  if(sparse%ja(i).eq.tcol)then
   val=sparse%a(i)
   exit
  endif
 enddo

end function

!** GET MEMORY
module function getmem_crs64(sparse) result(getmem)
 class(crssparse64),intent(in)::sparse
 integer(kind=int64)::getmem

 getmem=sparse%getmem_gen()
 if(allocated(sparse%ia))getmem=getmem+sizeof(sparse%ia)
 if(allocated(sparse%ja))getmem=getmem+sizeof(sparse%ja)
 if(allocated(sparse%a))getmem=getmem+sizeof(sparse%a)
#if (_PARDISO==1)
 getmem=getmem+sizeof(sparse%lpardisofirst)
#endif

end function

!**rowptr_crs64
module subroutine get_rowptr_crs64(sparse,ia)
  class(crssparse64),intent(in)::sparse
  integer(kind=int64),allocatable,intent(out)::ia(:)
 
  allocate(ia(size(sparse%ia, kind=int64))) 
  ia=sparse%ia

end subroutine
 
!**colval_crs64
module subroutine get_colval_crs64(sparse,ja)
  class(crssparse64),intent(in)::sparse
  integer(kind=int64),allocatable,intent(out)::ja(:)
 
  allocate(ja(size(sparse%ja, kind=int64))) 
  ja=sparse%ja

end subroutine
 
!**nzval_crs64
module subroutine get_nzval_crs64(sparse,a)
  class(crssparse64),intent(in)::sparse
  real(kind=wp),allocatable,intent(out)::a(:)
  
  allocate(a(size(sparse%a, kind=int64))) 
  a=sparse%a

end subroutine

!**MULTIPLICATIONS
module subroutine multgenv_csr64(sparse,alpha,trans,x,val,y)
 !Computes y=val*y+alpha*sparse(tranposition)*x
 class(crssparse64),intent(in)::sparse
 real(kind=wp),intent(in)::val,alpha
 real(kind=wp),intent(in)::x(:)
 real(kind=wp),intent(out)::y(:)
 character(len=1),intent(in)::trans

 write(sparse%unlog,'(a)')' Warning: mult not supported by crssparse64! Array returned = x'
 y = x

end subroutine

module subroutine multgenm_csr64(sparse,alpha,trans,x,val,y)
 !Computes y=val*y+alpha*sparse(tranposition)*x
 class(crssparse64),intent(in)::sparse
 real(kind=wp),intent(in)::val,alpha
 real(kind=wp),intent(in)::x(:,:)
 real(kind=wp),intent(out)::y(:,:)
 character(len=1),intent(in)::trans

 write(sparse%unlog,'(a)')' Warning: mult not supported by crssparse64! Array returned = x'
 y = x

end subroutine

!**NUMBER OF ELEMENTS
module function totalnumberofelements_crs64(sparse) result(nel)
 class(crssparse64),intent(in)::sparse
 integer(kind=int64)::nel

 nel=-1
 if(allocated(sparse%ia))&
 nel=int(sparse%ia(sparse%dim1+1),int64)-1_int64

end function

#if (_PARDISO==1)
!**RESET PARDISO MEMORY
module subroutine reset_pardiso_memory_crs64(sparse)
 !sparse*x=y
 class(crssparse64),intent(inout)::sparse

 !Pardiso variables
 integer(kind=int64)::error,idummy(1)
 integer(kind=int64)::nrhs

 if(.not.sparse%issquare())then
#if (_VERBOSE>0)
  write(sparse%unlog,'(a)')' Warning: the sparse matrix is not squared!'
  write(sparse%unlog,'(1x,a,1x,i0)')__FILE__,__LINE__
#endif
  return
 endif

 associate(parvar=>sparse%pardisovar)

 if(.not.allocated(parvar%pt))return

 sparse%lpardisofirst=.true.

 nrhs=1

 !Reset Pardiso memory
 parvar%phase=-1
 call pardiso(parvar%pt,parvar%maxfct,parvar%mnum,parvar%mtype,parvar%phase,&
              int(sparse%getdim(1),int64),sparse%a,sparse%ia,sparse%ja,&
              idummy,nrhs,parvar%iparm,parvar%msglvl,parvar%ddum,parvar%ddum,error)
 call checkpardiso(parvar%phase,error,sparse%unlog)

 end associate

 if(allocated(sparse%pardisovar%pt))deallocate(sparse%pardisovar%pt)

end subroutine
#endif

!**PRINT
module subroutine print_crs64(sparse,lint,output)
 class(crssparse64),intent(in)::sparse
 integer(kind=int32),intent(in),optional::output
 logical,intent(in),optional::lint

 integer(kind=int32)::i
 integer(kind=int32)::un
 integer(kind=int64)::j
 character(len=30),parameter::frm='(2(i0,1x),g0)'
 logical::linternal

 linternal=.true.
 if(present(lint))linternal=lint

 un=sparse%unlog
 if(present(output))un=output

 do i=1,sparse%dim1
  do j=sparse%ia(i),sparse%ia(i+1)-1
   write(un,frm)i,sparse%ja(j),sparse%a(j)
   if(.not.linternal.and.sparse%lupperstorage.and.sparse%lsymmetric.and.i.ne.sparse%ja(j))then
    write(un,frm)sparse%ja(j),i,sparse%a(j)
   endif
  enddo
 enddo

end subroutine

module subroutine printsquare_crs64(sparse,output)
 class(crssparse64),intent(inout)::sparse
 integer(kind=int32),intent(in),optional::output

 integer(kind=int32)::i,j,un
 real(kind=wp),allocatable::tmp(:)

 un=sparse%unlog
 if(present(output))un=output

 allocate(tmp(sparse%dim2))

 do i=1,sparse%dim1
  tmp=0.0_wp
  !could be implemented in a more efficient way
  do j=1,sparse%dim2
   tmp(j)=sparse%get(i,j)
  enddo
  write(un,'(*(g0.6,1x))')tmp
 enddo

 deallocate(tmp)

end subroutine

!**SAVE
module subroutine save_crs64(sparse,namefile)
 class(crssparse64),intent(in)::sparse
 character(len=*),intent(in)::namefile

 integer(kind=int32)::un

 open(newunit=un,file=namefile,action='write',status='replace',access='stream')!,buffered='yes')
 write(un)typecrs64              !int32
 write(un)sparse%dim1            !int32
 write(un)sparse%dim2            !int32
 write(un)sparse%nonzero()       !int64
 write(un)sparse%lupperstorage   !logical
 write(un)sparse%ia              !int64
 write(un)sparse%ja              !int64
 write(un)sparse%a               !wp
 close(un)

end subroutine

!**SCALE ALL ENTRIES
module subroutine scale_crs64(sparse,val)
 class(crssparse64),intent(inout)::sparse
 real(kind=wp),intent(in)::val
 sparse%a = sparse%a * val
end subroutine

!**SET ELEMENTS
module subroutine set_crs64(sparse,row,col,val,error)
 !add a value only to an existing one
 class(crssparse64),intent(inout)::sparse
 integer(kind=int32),intent(in)::row,col
 !added: error=0;Not existing: error=-1;matrix not initialized: error=-10
 integer(kind=int32),intent(out),optional::error
 real(kind=wp),intent(in)::val

 integer(kind=int64)::i

 if(present(error))error=0

 if(.not.validvalue_gen(sparse,row,col))return
 !if(.not.validnonzero_gen(sparse,val))return
 if(sparse%lupperstorage.and..not.uppervalue_gen(row,col))return

 if(sparse%ia(row).eq.0)then
  if(present(error))error=-10
  return
 endif

 do i=sparse%ia(row),sparse%ia(row+1)-1
  if(sparse%ja(i).eq.col)then
   sparse%a(i)=val
   if(present(error))error=0
   return
  endif
 enddo

 if(present(error))error=-1

end subroutine

!**SOLVE
#if (_PARDISO==1)
module subroutine solve_crs64_vector(sparse,x,y,msglvl)
 !sparse*x=y
 class(crssparse64),intent(inout)::sparse
 real(kind=wp),intent(out),contiguous::x(:)
 real(kind=wp),intent(inout),contiguous::y(:)
 integer(kind=int64),intent(in),optional::msglvl

 !Pardiso variables
 integer(kind=int32)::iparm(64)
 integer(kind=int64)::error
 integer(kind=int64)::nrhs
 integer(kind=int64)::msglvl_opt

 !$ real(kind=real64)::t1

 ! default value if not present
 msglvl_opt=1
 if(present(msglvl)) msglvl_opt=msglvl


 if(.not.sparse%issquare())then
  write(sparse%unlog,'(a)')' Warning: the sparse matrix is not squared!'
#if (_VERBOSE>0)
  write(sparse%unlog,'(1x,a,1x,i0)')__FILE__,__LINE__
#endif
  return
 endif

! associate(parvar=>sparse%pardisovar)

 nrhs=1 !always 1 since y is a vector

 if(sparse%lpardisofirst)then
  !$ t1=omp_get_wtime()
  !Sort the matrix
  call sparse%sort()

  !Preparation of Cholesky of A with Pardiso
  !parvar=pardiso_variable_64(maxfct=1,mnum=1,mtype=-2,solver=0,msglvl=1)   !mtype=11
  !to avoid ifort 2021.4.0 bug
  sparse%pardisovar=pardiso_variable_64(maxfct=1_int64,mnum=1_int64,mtype=-2_int64,solver=0_int64,msglvl=msglvl_opt)   !mtype=11

  !initialize iparm
  call pardisoinit(sparse%pardisovar%pt,int(sparse%pardisovar%mtype,int32),iparm)
  sparse%pardisovar%iparm=int(iparm,int64)
!  do i=1,64
!   write(sparse%unlog,*)'iparm',i,sparse%pardisovar%iparm(i)
!  enddo

  !Ordering and factorization
  sparse%pardisovar%phase=12
  sparse%pardisovar%iparm(2)=3
!  sparse%pardisovar%iparm(8)=1
#if (_VERBOSE == 1)
  sparse%pardisovar%iparm(19)=-1
#endif
  sparse%pardisovar%iparm(24)=0
  sparse%pardisovar%iparm(27)=1
#if (_DP==0)
  sparse%pardisovar%iparm(28)=1
#else
  sparse%pardisovar%iparm(28)=0
#endif
  write(sparse%unlog,'(a)')' Start ordering and factorization'
  if(allocated(sparse%perm64))then
   sparse%pardisovar%iparm(5)=1
   sparse%pardisovar%iparm(31)=0
   sparse%pardisovar%iparm(36)=0
   call pardiso(sparse%pardisovar%pt,sparse%pardisovar%maxfct,sparse%pardisovar%mnum,&
                sparse%pardisovar%mtype,sparse%pardisovar%phase,&
                int(sparse%getdim(1),int64),sparse%a,sparse%ia,sparse%ja,&
                sparse%perm64,nrhs,sparse%pardisovar%iparm,sparse%pardisovar%msglvl,&
                sparse%pardisovar%ddum,sparse%pardisovar%ddum,error)
  else
   call pardiso(sparse%pardisovar%pt,sparse%pardisovar%maxfct,sparse%pardisovar%mnum,&
                sparse%pardisovar%mtype,sparse%pardisovar%phase,&
                int(sparse%getdim(1),int64),sparse%a,sparse%ia,sparse%ja,&
                sparse%pardisovar%idum,nrhs,sparse%pardisovar%iparm,&
                sparse%pardisovar%msglvl,sparse%pardisovar%ddum,sparse%pardisovar%ddum,&
                error)
  endif
  call checkpardiso(sparse%pardisovar%phase,error,sparse%unlog)

  write(sparse%unlog,'(a,i0)')' Number of nonzeros in factors  = ',&
                              sparse%pardisovar%iparm(18)
#if (_VERBOSE == 1)
  write(sparse%unlog,'(a,i0)')' Number of factorization MFLOPS = ',&
                              sparse%pardisovar%iparm(19)
#endif
  !$ write(sparse%unlog,'(a,f0.5)')' Elapsed time (s)               = ',omp_get_wtime()-t1

  sparse%pardisovar%iparm(27)=0 !disable Pardiso checker
  sparse%lpardisofirst=.false.

 endif

 !Solving
 sparse%pardisovar%phase=33
 call pardiso(sparse%pardisovar%pt,sparse%pardisovar%maxfct,sparse%pardisovar%mnum,&
              sparse%pardisovar%mtype,sparse%pardisovar%phase,&
              int(sparse%getdim(1),int64),sparse%a,sparse%ia,sparse%ja,&
              sparse%pardisovar%idum,nrhs,sparse%pardisovar%iparm,&
              sparse%pardisovar%msglvl,y,x,error)
 call checkpardiso(sparse%pardisovar%phase,error,sparse%unlog)

#if (_VERBOSE>0)
 sparse%pardisovar%msglvl=1
#else
 sparse%pardisovar%msglvl=0
#endif

! end associate

end subroutine

module subroutine solve_crs64_array(sparse,x,y,msglvl)
 !sparse*x=y
 class(crssparse64),intent(inout)::sparse
 real(kind=wp),intent(out),contiguous::x(:,:)
 real(kind=wp),intent(inout),contiguous::y(:,:)
 integer(kind=int64),intent(in),optional::msglvl

 !Pardiso variables
 integer(kind=int32)::iparm(64)
 integer(kind=int64)::error
 integer(kind=int64)::nrhs
 integer(kind=int64)::msglvl_opt

 !$ real(kind=real64)::t1

 ! default value if not present
 msglvl_opt=1
 if(present(msglvl)) msglvl_opt=msglvl


 if(.not.sparse%issquare())then
  write(sparse%unlog,'(a)')' Warning: the sparse matrix is not squared!'
#if (_VERBOSE>0)
  write(sparse%unlog,'(1x,a,1x,i0)')__FILE__,__LINE__
#endif
  return
 endif

! associate(parvar=>sparse%pardisovar)

 nrhs=size(x,2)
 if(nrhs.ne.size(y,2))then
  write(sparse%unlog,'(a)')' ERROR: the number of colums of x and y provided to solve are different!'
  error stop
 endif

 if(sparse%lpardisofirst)then
  !$ t1=omp_get_wtime()
  !Sort the matrix
  call sparse%sort()

  !Preparation of Cholesky of A with Pardiso
  !parvar=pardiso_variable_64(maxfct=1,mnum=1,mtype=-2,solver=0,msglvl=1)   !mtype=11
  !to avoid ifort 2021.4.0 bug
  sparse%pardisovar=pardiso_variable_64(maxfct=1_int64,mnum=1_int64,mtype=-2_int64,solver=0_int64,msglvl=msglvl_opt)   !mtype=11

  !initialize iparm
  call pardisoinit(sparse%pardisovar%pt,int(sparse%pardisovar%mtype,int32),iparm)
  sparse%pardisovar%iparm=int(iparm,int64)
!  do i=1,64
!   write(sparse%unlog,*)'iparm',i,sparse%pardisovar%iparm(i)
!  enddo

  !Ordering and factorization
  sparse%pardisovar%phase=12
  sparse%pardisovar%iparm(2)=3
!  sparse%pardisovar%iparm(8)=1
#if (_VERBOSE == 1)
  sparse%pardisovar%iparm(19)=-1
#endif
  sparse%pardisovar%iparm(24)=0
  sparse%pardisovar%iparm(27)=1
#if (_DP==0)
  sparse%pardisovar%iparm(28)=1
#else
  sparse%pardisovar%iparm(28)=0
#endif
  write(sparse%unlog,'(a)')' Start ordering and factorization'
  if(allocated(sparse%perm64))then
   sparse%pardisovar%iparm(5)=1
   sparse%pardisovar%iparm(31)=0
   sparse%pardisovar%iparm(36)=0
   call pardiso(sparse%pardisovar%pt,sparse%pardisovar%maxfct,sparse%pardisovar%mnum,&
                sparse%pardisovar%mtype,sparse%pardisovar%phase,&
                int(sparse%getdim(1),int64),sparse%a,sparse%ia,sparse%ja,&
                sparse%perm64,nrhs,sparse%pardisovar%iparm,sparse%pardisovar%msglvl,&
                sparse%pardisovar%ddum,sparse%pardisovar%ddum,error)
  else
   call pardiso(sparse%pardisovar%pt,sparse%pardisovar%maxfct,sparse%pardisovar%mnum,&
                sparse%pardisovar%mtype,sparse%pardisovar%phase,&
                int(sparse%getdim(1),int64),sparse%a,sparse%ia,sparse%ja,&
                sparse%pardisovar%idum,nrhs,sparse%pardisovar%iparm,&
                sparse%pardisovar%msglvl,sparse%pardisovar%ddum,sparse%pardisovar%ddum,&
                error)
  endif
  call checkpardiso(sparse%pardisovar%phase,error,sparse%unlog)

  write(sparse%unlog,'(a,i0)')' Number of nonzeros in factors  = ',&
                               sparse%pardisovar%iparm(18)
#if (_VERBOSE == 1)
  write(sparse%unlog,'(a,i0)')' Number of factorization MFLOPS = ',&
                               sparse%pardisovar%iparm(19)
#endif
  !$ write(sparse%unlog,'(a,f0.5)')' Elapsed time (s)               = ',omp_get_wtime()-t1

  sparse%pardisovar%iparm(27)=0 !disable Pardiso checker
  sparse%lpardisofirst=.false.

 endif

 !Solving
 sparse%pardisovar%phase=33
 call pardiso(sparse%pardisovar%pt,sparse%pardisovar%maxfct,sparse%pardisovar%mnum,&
              sparse%pardisovar%mtype,sparse%pardisovar%phase,&
              int(sparse%getdim(1),int64),sparse%a,sparse%ia,sparse%ja,&
              sparse%pardisovar%idum,nrhs,sparse%pardisovar%iparm,&
              sparse%pardisovar%msglvl,y,x,error)
 call checkpardiso(sparse%pardisovar%phase,error,sparse%unlog)

#if (_VERBOSE>0)
 sparse%pardisovar%msglvl=1
#else
 sparse%pardisovar%msglvl=0
#endif

! end associate

end subroutine
#else
module subroutine solve_crs64_vector(sparse,x,y,msglvl)
 !sparse*x=y
 class(crssparse64),intent(inout)::sparse
 real(kind=wp),intent(out),contiguous::x(:)
 real(kind=wp),intent(inout),contiguous::y(:)
 integer(kind=int64),intent(in),optional::msglvl

 write(sparse%unlog,'(a)')' Warning: Pardiso is not enabled! Array returned = rhs'
 x=y

end subroutine

module subroutine solve_crs64_array(sparse,x,y,msglvl)
 !sparse*x=y
 class(crssparse64),intent(inout)::sparse
 real(kind=wp),intent(out),contiguous::x(:,:)
 real(kind=wp),intent(inout),contiguous::y(:,:)
 integer(kind=int64),intent(in),optional::msglvl

 write(sparse%unlog,'(a)')' Warning: Pardiso is not enabled! Array returned = rhs'
 x=y

end subroutine
#endif

!**SORT ARRAY
module subroutine sort_crs64(sparse)
 ! sort vectors ja and a by increasing order
 class(crssparse64),intent(inout)::sparse

 integer(kind=int32)::k
 integer(kind=int64)::endd,i,j,n,start,stkpnt
 integer(kind=int64)::d1,d2,d3,dmnmx,tmp
 integer(kind=int64)::stack(2,32)
 integer(kind=int64),allocatable::d(:)
 integer(kind=int64),parameter::select=20
 real(kind=wp)::umnmx,tmpu
 real(kind=wp),allocatable::u(:)

 if(sparse%issorted())then
  return
 endif

 do k=1,sparse%dim1
  n=sparse%ia(k+1)-sparse%ia(k)
  if(n.gt.1)then
   allocate(d(n),u(n))
   !copy of the vector to be sorted
   d=sparse%ja(sparse%ia(k):sparse%ia(k+1)-1)
   u=sparse%a(sparse%ia(k):sparse%ia(k+1)-1)
   !sort the vectors
   !from dlasrt.f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Quick return if possible
   stkpnt = 1
   stack( 1, 1 ) = 1
   stack( 2, 1 ) = n
   10 start = stack( 1, stkpnt )
   endd = stack( 2, stkpnt )
   stkpnt = stkpnt - 1
   IF( endd-start <= select .AND. endd-start > 0 ) THEN
   !Do Insertion sort on D( START:ENDD )
   !Sort into increasing order
     DO i = start + 1, endd
      DO j = i, start + 1, -1
       IF( d( j ) < d( j-1 ) ) THEN
         dmnmx = d( j )
         d( j ) = d( j-1 )
         d( j-1 ) = dmnmx
         umnmx = u( j )
         u( j ) = u( j-1 )
         u( j-1 ) = umnmx
       ELSE
         CYCLE
       END IF
      END DO
     END DO
   ELSE IF( endd-start > select ) THEN
   !Partition D( START:ENDD ) and stack parts, largest one first
   !Choose partition entry as median of 3
    d1 = d( start )
    d2 = d( endd )
    i = ( start + endd ) / 2
    d3 = d( i )
    IF( d1 < d2 ) THEN
      IF( d3 < d1 ) THEN
        dmnmx = d1
      ELSE IF( d3 < d2 ) THEN
        dmnmx = d3
      ELSE
        dmnmx = d2
      END IF
    ELSE
      IF( d3 < d2 ) THEN
        dmnmx = d2
      ELSE IF( d3 < d1 ) THEN
        dmnmx = d3
      ELSE
        dmnmx = d1
      END IF
    END IF
    !Sort into increasing order
     i = start - 1
     j = endd + 1
     90 j = j - 1
     IF( d( j ) > dmnmx ) GO TO 90
     110 i = i + 1
     IF( d( i ) < dmnmx ) GO TO 110
     IF( i < j ) THEN
       tmp = d( i )
       d( i ) = d( j )
       d( j ) = tmp
       tmpu = u( i )
       u( i ) = u( j )
       u( j ) = tmpu
       GO TO 90
     END IF
     IF( j-start > endd-j-1 ) THEN
       stkpnt = stkpnt + 1
       stack( 1, stkpnt ) = start
       stack( 2, stkpnt ) = j
       stkpnt = stkpnt + 1
       stack( 1, stkpnt ) = j + 1
       stack( 2, stkpnt ) = endd
     ELSE
       stkpnt = stkpnt + 1
       stack( 1, stkpnt ) = j + 1
       stack( 2, stkpnt ) = endd
       stkpnt = stkpnt + 1
       stack( 1, stkpnt ) = start
       stack( 2, stkpnt ) = j
     END IF
   END IF
   IF( stkpnt > 0 ) GO TO 10
   !end from dlasrt.f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !copy back the sorted vectors
   sparse%ja(sparse%ia(k):sparse%ia(k+1)-1)=d
   sparse%a(sparse%ia(k):sparse%ia(k+1)-1)=u
   deallocate(d,u)
  endif
 enddo

 call sparse%setsorted(.true.)

end subroutine

!**STATUS DECOMPOSITION
pure module function isdecomposed_crs(sparse) result(ll)
 class(crssparse),intent(in)::sparse
 logical::ll

 ll=(.not.sparse%lpardisofirst)

end function

module subroutine setdecomposed_crs(sparse,ll)
 class(crssparse),intent(inout)::sparse
 logical,intent(in)::ll

 sparse%lpardisofirst=(.not.ll)

end subroutine

end submodule

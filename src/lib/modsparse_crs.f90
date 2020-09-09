submodule (modsparse) modsparse_crs
#if (_DP==0)
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real32
#else
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real64
#endif
#if (_SPAINV==1)
 use modspainv
#endif
#if (_PARDISO==1)
 use mkl_pardiso
 use modvariablepardiso, only: checkpardiso, pardiso_variable
#endif
#if (_METIS==1)
 use modmetis, only: metis_nodend,metis_setoptions,metis_checkerror&
    ,METIS_CTYPE_RM,  METIS_IPTYPE_EDGE,  METIS_RTYPE_SEP2SIDED,  METIS_DBG_INFO
#endif
  !$ use omp_lib
 implicit none


contains

!**CONSTRUCTOR
module function constructor_crs(m,nel,n,lupper,unlog) result(sparse)
 type(crssparse)::sparse
 integer(kind=int32),intent(in)::m
 integer(kind=int32),intent(in)::nel
 integer(kind=int32),intent(in),optional::n,unlog
 logical,intent(in),optional::lupper

 call sparse%initialize('CRS',m,m)

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

!**DIAGONAL ELEMENTS
module function diag_vect_crs(sparse) result(array)
 class(crssparse),intent(inout)::sparse
 real(kind=wp),allocatable::array(:)

 integer(kind=int32)::ndiag,i

 ndiag=min(sparse%dim1,sparse%dim2)

 allocate(array(ndiag))
 array=0.0_wp

 do i=1,ndiag
  array(i)=sparse%get(i,i)
 enddo

end function

#if (_SPAINV==1)
!**GET (COMPLETE) CHOLESKY FACTOR
module subroutine getchol_crs(sparse,minsizenode)
 class(crssparse),intent(inout)::sparse
 integer(kind=int32),intent(in),optional::minsizenode

 type(metisgraph)::metis
#if (_VERBOSE>0)
 !$ real(kind=real64)::t1,t2

 !$ t1=omp_get_wtime()
 !$ t2=t1
#endif

 call sparse%sort()

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'CHOL CRS sorting',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 !Ordering
 if(.not.allocated(sparse%perm))then
#if (_METIS==1)
  call sparse%setpermutation(sparse%getordering())
#else
  write(sparse%unlog,'(a)')' ERROR: A permutation vector must be set before calling chol'
  error stop
#endif
 endif

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'CHOL CRS ordering',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 metis=sparse

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'CHOL METIS=CRS',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 if(present(minsizenode))then
  call get_chol(sparse%ia,sparse%ja,sparse%a,metis%xadj,metis%adjncy,sparse%perm,minsizenode=minsizenode,un=sparse%unlog)
 else
  call get_chol(sparse%ia,sparse%ja,sparse%a,metis%xadj,metis%adjncy,sparse%perm,un=sparse%unlog)
 endif

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'CHOL CRS Chol. fact.',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 call sparse%setsorted(.false.)
 call sparse%sort()

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'CHOL CRS sorting fact.',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'CHOL CRS',': Total   time (s) = ',omp_get_wtime()-t2
#endif

end subroutine

!**GET LDLt DECOMPOSITION
module subroutine getldlt_crs(sparse,minsizenode)
 class(crssparse),intent(inout)::sparse
 integer(kind=int32),intent(in),optional::minsizenode

 integer(kind=int32)::i,j
 real(kind=wp)::s
! real(kind=wp),allocatable::s(:)
#if (_VERBOSE>0)
 !$ real(kind=real64)::t1,t2

 !$ t1=omp_get_wtime()
 !$ t2=t1
#endif

 if(present(minsizenode))then
  call sparse%chol(minsizenode)
 else
  call sparse%chol()
 endif

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'LDLt Cholesky fact.',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 do i=1,sparse%getdim(1)
  s=0._wp
  if(sparse%a(sparse%ia(i)).gt.0._wp)s=1._wp/sparse%a(sparse%ia(i))   !aaaa use a tolerance factor
  sparse%a(sparse%ia(i))=sparse%a(sparse%ia(i))**2
  do j=sparse%ia(i)+1,sparse%ia(i+1)-1
   sparse%a(j)=s*sparse%a(j)
  enddo
 enddo

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'LDLt CRS Decomposition',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'LDLt CRS',': Total   time (s) = ',omp_get_wtime()-t2
#endif

end subroutine
#endif

#if (_METIS==1)
!**GET ORDER
module function getordering_crs(sparse&
                          ,ctype,iptype,rtype,compress,ccorder&
                          ,pfactor,nseps,bglvl&
                          ) result(perm)
 class(crssparse),intent(in)::sparse
 integer(kind=int32),intent(in),optional::ctype,iptype,rtype,compress,ccorder,pfactor,nseps,bglvl
 integer(kind=int32),allocatable::perm(:)

 integer(kind=int32)::err
 integer(kind=int32)::pctype,piptype,prtype,pcompress,pccorder,ppfactor,pnseps,pbglvl
 integer(kind=int32),allocatable::options(:)
 integer(kind=int32),allocatable::iperm(:)
 type(metisgraph)::metis
#if (_VERBOSE>0)
 !$ real(kind=real64)::t1

 !$ t1=omp_get_wtime()
#endif

 metis=sparse

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'METIS=CRS',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 pctype=METIS_CTYPE_RM
 if(present(ctype))pctype=ctype

 piptype=METIS_IPTYPE_EDGE
 if(present(iptype))piptype=iptype

 prtype=METIS_RTYPE_SEP2SIDED
 if(present(rtype))prtype=rtype

 pcompress=0
 if(present(compress))pcompress=compress

 pccorder=0
 if(present(ccorder))pccorder=ccorder

 ppfactor=0
 if(present(pfactor))ppfactor=pfactor

 if(sparse%getdim(1).lt.50000)then
  pnseps=1
 elseif(sparse%getdim(1).lt.200000)then
  pnseps=5
 else
  pnseps=10
 endif
 if(present(nseps))pnseps=nseps

#if (_VERBOSE>1)
 pbglvl=METIS_DBG_INFO
#endif
 if(present(bglvl))pbglvl=bglvl

 err=metis_setoptions(options&
                      ,ctype=pctype,iptype=piptype,rtype=prtype,compress=pcompress&
                      ,ccorder=pccorder,pfactor=ppfactor,nseps=pnseps,dbglvl=pbglvl&
                      )
 call metis_checkerror(err,sparse%unlog)

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'METIS options setting',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 allocate(perm(metis%nvertices),iperm(metis%nvertices))

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'METIS arrays allocation',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 err=metis_nodend(metis%nvertices,metis%xadj,metis%adjncy,metis%vwgt,options,perm,iperm)
 call metis_checkerror(err,sparse%unlog)

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'METIS ordering',': Elapsed time (s) = ',omp_get_wtime()-t1
#endif

end function
#endif

#if (_SPAINV==1)
!**GET INCOMPLETE CHOLESKY FACTOR
module subroutine getichol_crs(sparse,minsizenode)
 class(crssparse),intent(inout)::sparse
 integer(kind=int32),intent(in),optional::minsizenode

 type(metisgraph)::metis
#if (_VERBOSE>0)
 !$ real(kind=real64)::t1,t2

 !$ t1=omp_get_wtime()
 !$ t2=t1
#endif

 call sparse%sort()

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'ICHOL CRS sorting',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 !Ordering
 if(.not.allocated(sparse%perm))then
#if (_METIS==1)
  call sparse%setpermutation(sparse%getordering())
#else
  write(sparse%unlog,'(a)')' ERROR: A permutation vector must be set before calling ichol'
  error stop
#endif
 endif

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'ICHOL CRS ordering',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 metis=sparse

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'ICHOL METIS=CRS',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 if(present(minsizenode))then
  call get_ichol(sparse%ia,sparse%ja,sparse%a,metis%xadj,metis%adjncy,sparse%perm,minsizenode=minsizenode,un=sparse%unlog)
 else
  call get_ichol(sparse%ia,sparse%ja,sparse%a,metis%xadj,metis%adjncy,sparse%perm,un=sparse%unlog)
 endif

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'ICHOL CRS Chol. fact.',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'ICHOL CRS',': Total   time (s) = ',omp_get_wtime()-t2
#endif

end subroutine

!**GET SPARSE INVERSE
module subroutine getspainv_crs(sparse,minsizenode)
 class(crssparse),intent(inout)::sparse
 integer(kind=int32),intent(in),optional::minsizenode

 type(metisgraph)::metis
#if (_VERBOSE>0)
 !$ real(kind=real64)::t1,t2

 !$ t1=omp_get_wtime()
 !$ t2=t1
#endif

 call sparse%sort()

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'SPAINV CRS sorting',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 !Ordering
 if(.not.allocated(sparse%perm))then
#if (_METIS==1)
  call sparse%setpermutation(sparse%getordering())
#else
  write(sparse%unlog,'(a)')' ERROR: A permutation vector must be set before calling spainv'
  error stop
#endif
 endif

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'SPAINV CRS ordering',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 metis=sparse

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'SPAINV METIS=CRS',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 if(present(minsizenode))then
  call get_spainv(sparse%ia,sparse%ja,sparse%a,metis%xadj,metis%adjncy,sparse%perm,minsizenode=minsizenode,un=sparse%unlog)
 else
  call get_spainv(sparse%ia,sparse%ja,sparse%a,metis%xadj,metis%adjncy,sparse%perm,un=sparse%unlog)
 endif

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'SPAINV CRS inversion',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'SPAINV CRS',': Total   time (s) = ',omp_get_wtime()-t2
#endif

end subroutine
#endif

#if (_PARDISO==1)
!**RESET PARDISO MEMORY
module subroutine reset_pardiso_memory_crs(sparse)
 !sparse*x=y
 class(crssparse),intent(inout)::sparse

 !Pardiso variables
 integer(kind=int32)::error,idummy(1)
 integer(kind=int32)::nrhs

 if(.not.sparse%issquare())then
#if (_VERBOSE>0)
  write(sparse%unlog,'(a)')' Warning: the sparse matrix is not squared!'
  write(sparse%unlog,'(x,a,x,i0)')__FILE__,__LINE__
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
              sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,&
              idummy,nrhs,parvar%iparm,parvar%msglvl,parvar%ddum,parvar%ddum,error)
 call checkpardiso(parvar%phase,error,sparse%unlog)

 end associate

end subroutine
#endif

!**SOLVE
#if (_PARDISO==1)
module subroutine solve_crs_vector(sparse,x,y)
 !sparse*x=y
 class(crssparse),intent(inout)::sparse
 real(kind=wp),intent(out)::x(:)
 real(kind=wp),intent(inout)::y(:)

 !Pardiso variables
 integer(kind=int32)::error,phase
 integer(kind=int32)::nrhs

 integer(kind=int32)::i
 !$ real(kind=real64)::t1

 if(.not.sparse%issquare())then
#if (_VERBOSE>0)
  write(sparse%unlog,'(a)')' Warning: the sparse matrix is not squared!'
  write(sparse%unlog,'(x,a,x,i0)')__FILE__,__LINE__
#endif
  return
 endif

 associate(parvar=>sparse%pardisovar)

 nrhs=1 !always 1 since y is a vector

 if(sparse%lpardisofirst)then
  !$ t1=omp_get_wtime()
  !Sort the matrix
  call sparse%sort()

  !Preparation of Cholesky of A with Pardiso
  parvar=pardiso_variable(maxfct=1,mnum=1,mtype=-2,solver=0,msglvl=1)   !mtype=11

  !initialize iparm
  call pardisoinit(parvar%pt,parvar%mtype,parvar%iparm)
!  do i=1,64
!   write(sparse%unlog,*)'iparm',i,parvar%iparm(i)
!  enddo

  !Ordering and factorization
  parvar%phase=12
  parvar%iparm(2)=3
!  parvar%iparm(8)=1
  parvar%iparm(27)=1
#if (_DP==0)
  parvar%iparm(28)=1
#else
  parvar%iparm(28)=0
#endif
  write(sparse%unlog,'(a)')' Start ordering and factorization'
  if(allocated(sparse%perm))then
   parvar%iparm(5)=1;parvar%iparm(31)=0;parvar%iparm(36)=0
   call pardiso(parvar%pt,parvar%maxfct,parvar%mnum,parvar%mtype,parvar%phase,&
                sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,&
                sparse%perm,nrhs,parvar%iparm,parvar%msglvl,parvar%ddum,parvar%ddum,error)
  else
   call pardiso(parvar%pt,parvar%maxfct,parvar%mnum,parvar%mtype,parvar%phase,&
                sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,&
                parvar%idum,nrhs,parvar%iparm,parvar%msglvl,parvar%ddum,parvar%ddum,error)
  endif
  call checkpardiso(parvar%phase,error,sparse%unlog)

  write(sparse%unlog,'(a,i0)')' Number of nonzeros in factors  = ',parvar%iparm(18)
  write(sparse%unlog,'(a,i0)')' Number of factorization MFLOPS = ',parvar%iparm(19)
  !$ write(sparse%unlog,'(a,f0.5)')' Elapsed time (s)               = ',omp_get_wtime()-t1
 endif

 !Solving
 parvar%phase=33
 parvar%iparm(27)=0
 call pardiso(parvar%pt,parvar%maxfct,parvar%mnum,parvar%mtype,parvar%phase,&
              sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,&
              parvar%idum,nrhs,parvar%iparm,parvar%msglvl,y,x,error)
 call checkpardiso(parvar%phase,error,sparse%unlog)

#if (_VERBOSE>0)
 parvar%msglvl=1
#else
 parvar%msglvl=0
#endif

 sparse%lpardisofirst=.false.

 end associate

end subroutine

module subroutine solve_crs_array(sparse,x,y)
 !sparse*x=y
 class(crssparse),intent(inout)::sparse
 real(kind=wp),intent(out)::x(:,:)
 real(kind=wp),intent(inout)::y(:,:)

 !Pardiso variables
 integer(kind=int32)::error,phase
 integer(kind=int32)::nrhs

 integer(kind=int32)::i
 !$ real(kind=real64)::t1

 if(.not.sparse%issquare())then
  write(sparse%unlog,'(a)')' Warning: the sparse matrix is not squared!'
#if (_VERBOSE>0)
  write(sparse%unlog,'(x,a,x,i0)')__FILE__,__LINE__
#endif
  return
 endif

 associate(parvar=>sparse%pardisovar)

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
  parvar=pardiso_variable(maxfct=1,mnum=1,mtype=-2,solver=0,msglvl=1)   !mtype=11

  !initialize iparm
  call pardisoinit(parvar%pt,parvar%mtype,parvar%iparm)
!  do i=1,64
!   write(sparse%unlog,*)'iparm',i,parvar%iparm(i)
!  enddo

  !Ordering and factorization
  parvar%phase=12
  parvar%iparm(2)=3
  parvar%iparm(27)=1
#if (_DP==0)
  parvar%iparm(28)=1
#else
  parvar%iparm(28)=0
#endif
  write(sparse%unlog,'(a)')' Start ordering and factorization'
  if(allocated(sparse%perm))then
   parvar%iparm(5)=1;parvar%iparm(31)=0;parvar%iparm(36)=0
   call pardiso(parvar%pt,parvar%maxfct,parvar%mnum,parvar%mtype,parvar%phase,&
                sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,&
                sparse%perm,nrhs,parvar%iparm,parvar%msglvl,parvar%ddum,parvar%ddum,error)
  else
   call pardiso(parvar%pt,parvar%maxfct,parvar%mnum,parvar%mtype,parvar%phase,&
                sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,&
                parvar%idum,nrhs,parvar%iparm,parvar%msglvl,parvar%ddum,parvar%ddum,error)
  endif
  call checkpardiso(parvar%phase,error,sparse%unlog)

  write(sparse%unlog,'(a,i0)')' Number of nonzeros in factors  = ',parvar%iparm(18)
  write(sparse%unlog,'(a,i0)')' Number of factorization MFLOPS = ',parvar%iparm(19)
  !$ write(sparse%unlog,'(a,f0.5)')' Elapsed time (s)               = ',omp_get_wtime()-t1
 endif

 !Solving
 parvar%phase=33
 parvar%iparm(27)=0
 call pardiso(parvar%pt,parvar%maxfct,parvar%mnum,parvar%mtype,parvar%phase,&
              sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,&
              parvar%idum,nrhs,parvar%iparm,parvar%msglvl,y,x,error)
 call checkpardiso(parvar%phase,error,sparse%unlog)

 parvar%msglvl=0
 sparse%lpardisofirst=.false.

 end associate

end subroutine
#else
module subroutine solve_crs_vector(sparse,x,y)
 !sparse*x=y
 class(crssparse),intent(in)::sparse
 real(kind=wp),intent(out)::x(:)
 real(kind=wp),intent(inout)::y(:)

 write(sparse%unlog,'(a)')' Warning: Pardiso is not enabled! Array returned = rhs'
 x=y

end subroutine

module subroutine solve_crs_array(sparse,x,y)
 !sparse*x=y
 class(crssparse),intent(in)::sparse
 real(kind=wp),intent(out)::x(:,:)
 real(kind=wp),intent(inout)::y(:,:)

 write(sparse%unlog,'(a)')' Warning: Pardiso is not enabled! Array returned = rhs'
 x=y

end subroutine
#endif

!**SOLVE WITH A TRIANGULAR FACTOR
module subroutine isolve_crs(sparse,x,y)
 !sparse*x=y
 class(crssparse),intent(in)::sparse
 real(kind=wp),intent(out)::x(:)
 real(kind=wp),intent(in)::y(:)

 integer(kind=int32)::i
 real(kind=wp),allocatable::x_(:)
#if (_VERBOSE>0)
 !$ real(kind=real64)::t1,t2

 !$ t1=omp_get_wtime()
 !$ t2=omp_get_wtime()
#endif

 allocate(x_(1:size(x)))

 do i=1,size(x)
  x_(i)=y(sparse%perm(i))
 enddo

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'ISOLVE CRS y permutation',': Elapsed time (s) = ',omp_get_wtime()-t2
 !$ t2=omp_get_wtime()
#endif

#if(_DP==0)
 call mkl_scsrtrsv('U','T','N',sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,x_,x)
#else
 call mkl_dcsrtrsv('U','T','N',sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,x_,x)
#endif

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'ISOLVE CRS 1st triangular solve',': Elapsed time (s) = ',omp_get_wtime()-t2
 !$ t2=omp_get_wtime()
#endif

#if(_DP==0)
 call mkl_scsrtrsv('U','N','N',sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,x,x_)
#else
 call mkl_dcsrtrsv('U','N','N',sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,x,x_)
#endif

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'ISOLVE CRS 2nd triangular solve',': Elapsed time (s) = ',omp_get_wtime()-t2
 !$ t2=omp_get_wtime()
#endif

 do i=1,size(x)
  x(sparse%perm(i))=x_(i)
 enddo

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'ISOLVE CRS x permutation',': Elapsed time (s) = ',omp_get_wtime()-t2
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'ISOLVE CRS',': Total   time (s) = ',omp_get_wtime()-t1
#endif

end subroutine

!**SOLVE WITH LDLt DECOMPOSITION
module subroutine solveldlt_crs(sparse,x,y)
 !sparse*x=y
 class(crssparse),intent(in)::sparse
 real(kind=wp),intent(out)::x(:)
 real(kind=wp),intent(in)::y(:)

 integer(kind=int32)::i
 real(kind=wp),allocatable::x_(:)
#if (_VERBOSE>0)
 !$ real(kind=real64)::t1,t2

 !$ t1=omp_get_wtime()
 !$ t2=omp_get_wtime()
#endif

 allocate(x_(1:size(x)))

 do i=1,size(x)
  x_(i)=y(sparse%perm(i))
 enddo

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'SOLVE LDLt CRS y permutation',': Elapsed time (s) = ',omp_get_wtime()-t2
 !$ t2=omp_get_wtime()
#endif

#if(_DP==0)
 call mkl_scsrtrsv('U','T','U',sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,x_,x)
#else
 call mkl_dcsrtrsv('U','T','U',sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,x_,x)
#endif
#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'SOLVE LDLt CRS 1st triangular solve',': Elapsed time (s) = ',omp_get_wtime()-t2
 !$ t2=omp_get_wtime()
#endif
 do i=1,sparse%getdim(1)
  if(sparse%a(sparse%ia(i)).gt.tol)then  !aaa must use a tol parameter
   x(i)=x(i)/sparse%a(sparse%ia(i))
  else
   x(i)=0._wp
  endif
 enddo

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'SOLVE LDLt CRS Diagonal solve',': Elapsed time (s) = ',omp_get_wtime()-t2
 !$ t2=omp_get_wtime()
#endif

#if(_DP==0)
 call mkl_scsrtrsv('U','N','U',sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,x,x_)
#else
 call mkl_dcsrtrsv('U','N','U',sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,x,x_)
#endif

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'SOLVE LDLt CRS 2nd triangular solve',': Elapsed time (s) = ',omp_get_wtime()-t2
 !$ t2=omp_get_wtime()
#endif

 do i=1,size(x)
  x(sparse%perm(i))=x_(i)
 enddo

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'SOLVE LDLt CRS x permutation',': Elapsed time (s) = ',omp_get_wtime()-t2
 !$ write(sparse%unlog,'(x,a,t30,a,f0.5)')'SOLVE LDLt CRS',': Total   time (s) = ',omp_get_wtime()-t1
#endif

end subroutine



end submodule

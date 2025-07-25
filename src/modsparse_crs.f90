submodule (modsparse) modsparse_crs
 use modrandom, only: setseed, snorm=>rand_stdnormal
 use modsparse_mkl, only: pardisoinit, pardiso &
                          , mkl_scsrmv, mkl_dcsrmv &
                          , mkl_scsrmm, mkl_dcsrmm &
                          , mkl_scsrtrsv, mkl_dcsrtrsv &
                          , mkl_scsrsymv, mkl_dcsrsymv
 use modsparse_inv, only: get_chol, get_ichol, get_spainv
#if (_PARDISO==1)
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

 sparse%loriginal=.true.

end function

module subroutine constructor_sub_crs(sparse,m,nel,n,lupper,unlog)
 class(crssparse),intent(out)::sparse
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

 sparse%loriginal=.true.

end subroutine

!**DESTROY
module impure elemental subroutine destroy_crs(sparse)
 class(crssparse),intent(inout)::sparse

#if(_PARDISO==1)
 call sparse%resetpardiso()
#endif

 call sparse%destroy_gen_gen()

 if(allocated(sparse%ia))deallocate(sparse%ia)
 if(allocated(sparse%ja))deallocate(sparse%ja)
 if(allocated(sparse%a))deallocate(sparse%a)

end subroutine

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

!**ADD ELEMENTS
module subroutine add_crs(sparse,row,col,val,error)
 !add a value only to an existing one
 class(crssparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::row,col
 !added: error=0;Not existing: error=-1;matrix not initialized: error=-10
 integer(kind=int32),intent(out),optional::error
 real(kind=wp),intent(in)::val

 integer(kind=int32)::i

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
pure module function get_crs(sparse,row,col) result(val)
 class(crssparse),intent(in)::sparse
 integer(kind=int32),intent(in)::row,col
 real(kind=wp)::val

 integer(kind=int32)::i
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
module function getmem_crs(sparse) result(getmem)
 class(crssparse),intent(in)::sparse
 integer(kind=int64)::getmem

 getmem=sparse%getmem_gen()
 if(allocated(sparse%ia))getmem=getmem+sizeof(sparse%ia)
 if(allocated(sparse%ja))getmem=getmem+sizeof(sparse%ja)
 if(allocated(sparse%a))getmem=getmem+sizeof(sparse%a)
 getmem=getmem+sizeof(sparse%loriginal)

end function

!**EXTERNAL
module subroutine external_crs(sparse,ia,ja,a)
 class(crssparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::ia(:),ja(:)
 real(kind=wp),intent(in)::a(:)

 if(size(ia).ne.size(sparse%ia))then
  write(sparse%unlog,'(a)')' ERROR: The provided array ia is of a different size!'
  error stop
 endif
 if(size(ja).ne.size(sparse%ja))then
  write(sparse%unlog,'(a)')' ERROR: The provided array ja is of a different size!'
  error stop
 endif
 if(size(a).ne.size(sparse%a))then
  write(sparse%unlog,'(a)')' ERROR: The provided array a is of a different size!'
  error stop
 endif

 sparse%ia=ia
 sparse%ja=ja
 sparse%a=a

end subroutine

!**HARVILLE
module subroutine harville_crs(sparse, ngibbs, nburn, diaginv, seed)
 ! Approximate the diagonal elements of the inverse of the sparse matrix
 ! following Harville (1999)
 !
 !Possibility to implement multiple chains by extenting x to a matrix x(n,m)
 !where m is the number of independent chains.
 class(crssparse), intent(inout)::sparse
 integer, intent(in) :: ngibbs, nburn
 real(wp), intent(out), allocatable :: diaginv(:)
 integer, intent(in), optional :: seed

 integer :: n, i, j, k
 real(wp) :: d
 real(wp), allocatable :: a(:)
 real(wp), allocatable :: w(:), xt(:), dinv(:) ,x(:)

 if(.not.sparse%lsymmetric.or..not.sparse%lupperstorage) return

 call sparse%sort()

 n = size(sparse%ia) - 1

 allocate(a, source = sparse%a)
 allocate(w(n), source = 0._wp)

 do i = 1, n
  if (a(sparse%ia(i)).ne.0._wp) then
   w(i) = 1 / a(sparse%ia(i))
   a(sparse%ia(i)) = 0
  else
   w(i) = 0
  endif
 enddo

 allocate(xt(n), source = 0._wp)
 allocate(dinv(n), source = 0._wp)
 allocate(x(n), source = 0._wp)

 if(present(seed))then
  call setseed(seed)
 else
  call setseed(1)
 endif

 do i = 1, ngibbs
#if(_DP==0)
  call mkl_scsrsymv('U', n, a, sparse%ia, sparse%ja, xt, x)
#else
  call mkl_dcsrsymv('U', n, a, sparse%ia, sparse%ja, xt, x)
#endif
  x = -1 * x
  do j=1,n
   x(j) = w(j) * x(j)
   if (i.gt.nburn) dinv(j) = dinv(j) + x(j)**2
   x(j) = x(j) + sqrt(w(j)) * snorm()
   d = x(j) - xt(j)
   do k = sparse%ia(j)+1, sparse%ia(j+1)-1
    x(sparse%ja(k)) = x(sparse%ja(k)) - a(k) * d
   enddo
  enddo
  xt = x
  x  = 0
 enddo
 deallocate(a)
 deallocate(x)
 deallocate(xt)

 allocate(diaginv(n), source = 0._wp)
 do j = 1, n
  diaginv(j) = w(j) + dinv(j)/(ngibbs-nburn)
 enddo

end subroutine

!**rowptr_crs
pure module subroutine get_rowptr_crs(sparse,ia)
  class(crssparse),intent(in)::sparse
  integer(kind=int32),allocatable,intent(out)::ia(:)
 
  allocate(ia, source = sparse%ia)

end subroutine
 
!**colval_crs
pure module subroutine get_colval_crs(sparse,ja)
  class(crssparse),intent(in)::sparse
  integer(kind=int32),allocatable,intent(out)::ja(:)
 
  allocate(ja, source = sparse%ja)

end subroutine
 
!**nzval_crs
pure module subroutine get_nzval_crs(sparse,a)
  class(crssparse),intent(in)::sparse
  real(kind=wp),allocatable,intent(out)::a(:)
  
  allocate(a, source = sparse%a)

end subroutine

!**MULTIPLICATIONS
module subroutine multgenv_csr(sparse,alpha,trans,x,val,y)
 !Computes y=val*y+alpha*sparse(tranposition)*x
 class(crssparse),intent(in)::sparse
 real(kind=wp),intent(in)::val,alpha
 real(kind=wp),intent(in)::x(:)
 real(kind=wp),intent(out)::y(:)
 character(len=1),intent(in)::trans

 character(len=1)::matdescra(6)

 matdescra=''

 if(sparse%lsymmetric.and.sparse%lupperstorage)then
  matdescra(1)='S'
 elseif(.not.sparse%lsymmetric.and.sparse%lupperstorage)then
  matdescra(1)='T'
 elseif(.not.sparse%lsymmetric.and..not.sparse%lupperstorage)then
  matdescra(1)='G'
 else
  write(sparse%unlog,'(a)')'  ERROR (multbyv): unsupported format'
  call sparse%printstats
  error stop
 endif

 if(sparse%lupperstorage)then
  matdescra(2)='U'
  matdescra(3)='N'
 endif

 matdescra(4)='F'

#if(_DP==0)
 call mkl_scsrmv(trans,sparse%dim1,sparse%dim2,alpha,matdescra&
                  ,sparse%a,sparse%ja,sparse%ia(1:sparse%dim1),sparse%ia(2:sparse%dim1+1)&
                  ,x,val,y)
#else
 call mkl_dcsrmv(trans,sparse%dim1,sparse%dim2,alpha,matdescra&
                  ,sparse%a,sparse%ja,sparse%ia(1:sparse%dim1),sparse%ia(2:sparse%dim1+1)&
                  ,x,val,y)
#endif

end subroutine

module subroutine multgenm_csr(sparse,alpha,trans,x,val,y)
 !Computes y=val*y+alpha*sparse(tranposition)*x
 class(crssparse),intent(in)::sparse
 real(kind=wp),intent(in)::val,alpha
 real(kind=wp),intent(in)::x(:,:)
 real(kind=wp),intent(out)::y(:,:)
 character(len=1),intent(in)::trans

 character(len=1)::matdescra(6)

 matdescra=''

 if(sparse%lsymmetric.and.sparse%lupperstorage)then
  matdescra(1)='S'
 elseif(.not.sparse%lsymmetric.and.sparse%lupperstorage)then
  matdescra(1)='T'
 elseif(.not.sparse%lsymmetric.and..not.sparse%lupperstorage)then
  matdescra(1)='G'
 else
  write(sparse%unlog,'(a)')'  ERROR (multbyv): unsupported format'
  call sparse%printstats
  error stop
 endif

 if(sparse%lupperstorage)then
  matdescra(2)='U'
  matdescra(3)='N'
 endif

 matdescra(4)='F'

#if(_DP==0)
 call mkl_scsrmm(trans,sparse%dim1,size(y,2),sparse%dim2,&
                 alpha,matdescra,sparse%a,sparse%ja,sparse%ia(1:sparse%dim1),sparse%ia(2:sparse%dim1+1),&
                 x,size(x,1),&
                 val,y,size(y,1))
#else
 call mkl_dcsrmm(trans,sparse%dim1,size(y,2),sparse%dim2,&
                 alpha,matdescra,sparse%a,sparse%ja,sparse%ia(1:sparse%dim1),sparse%ia(2:sparse%dim1+1),&
                 x,size(x,1),&
                 val,y,size(y,1))
#endif

end subroutine

!**NUMBER OF ELEMENTS
module function totalnumberofelements_crs(sparse) result(nel)
 class(crssparse),intent(in)::sparse
 integer(kind=int64)::nel

 nel=-1
 if(allocated(sparse%ia))nel=sparse%ia(sparse%dim1+1)-1_int64

end function

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

 if (sparse%isdecomposed())return

 call sparse%sort()

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'CHOL CRS sorting',': Elapsed time (s) = ',omp_get_wtime()-t1
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
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'CHOL CRS ordering',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 metis=sparse

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'CHOL METIS=CRS',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 if(present(minsizenode))then
  call get_chol(sparse%ia,sparse%ja,sparse%a,metis%xadj,metis%adjncy,sparse%perm,minsizenode=minsizenode,un=sparse%unlog)
 else
  call get_chol(sparse%ia,sparse%ja,sparse%a,metis%xadj,metis%adjncy,sparse%perm,un=sparse%unlog)
 endif

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'CHOL CRS Chol. fact.',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 call sparse%setsymmetric(.false.)
 call sparse%setdecomposed(.true.)

 call sparse%setsorted(.false.)
 call sparse%sort()

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'CHOL CRS sorting fact.',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'CHOL CRS',': Total   time (s) = ',omp_get_wtime()-t2
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
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'LDLt Cholesky fact.',': Elapsed time (s) = ',omp_get_wtime()-t1
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
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'LDLt CRS Decomposition',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'LDLt CRS',': Total   time (s) = ',omp_get_wtime()-t2
#endif

end subroutine

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
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'METIS=CRS',': Elapsed time (s) = ',omp_get_wtime()-t1
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
#else
 pbglvl=0
#endif
 if(present(bglvl))pbglvl=bglvl

 err=metis_setoptions(options&
                      ,ctype=pctype,iptype=piptype,rtype=prtype,compress=pcompress&
                      ,ccorder=pccorder,pfactor=ppfactor,nseps=pnseps,dbglvl=pbglvl&
                      )
 call metis_checkerror(err,sparse%unlog)

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'METIS options setting',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 allocate(perm(metis%nvertices),iperm(metis%nvertices))

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'METIS arrays allocation',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 err=metis_nodend(metis%nvertices,metis%xadj,metis%adjncy,metis%vwgt,options,perm,iperm)
 call metis_checkerror(err,sparse%unlog)

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'METIS ordering',': Elapsed time (s) = ',omp_get_wtime()-t1
#endif

end function
#endif

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
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'ICHOL CRS sorting',': Elapsed time (s) = ',omp_get_wtime()-t1
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
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'ICHOL CRS ordering',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 metis=sparse

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'ICHOL METIS=CRS',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 if(present(minsizenode))then
  call get_ichol(sparse%ia,sparse%ja,sparse%a,metis%xadj,metis%adjncy,sparse%perm,minsizenode=minsizenode,un=sparse%unlog)
 else
  call get_ichol(sparse%ia,sparse%ja,sparse%a,metis%xadj,metis%adjncy,sparse%perm,un=sparse%unlog)
 endif

 call sparse%setsymmetric(.false.)
 call sparse%setdecomposed(.true.)

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'ICHOL CRS Chol. fact.',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'ICHOL CRS',': Total   time (s) = ',omp_get_wtime()-t2
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
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'SPAINV CRS sorting',': Elapsed time (s) = ',omp_get_wtime()-t1
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
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'SPAINV CRS ordering',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 metis=sparse

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'SPAINV METIS=CRS',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#endif

 if(present(minsizenode))then
  call get_spainv(sparse%ia,sparse%ja,sparse%a,metis%xadj,metis%adjncy,sparse%perm,minsizenode=minsizenode,un=sparse%unlog)
 else
  call get_spainv(sparse%ia,sparse%ja,sparse%a,metis%xadj,metis%adjncy,sparse%perm,un=sparse%unlog)
 endif

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'SPAINV CRS inversion',': Elapsed time (s) = ',omp_get_wtime()-t1
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'SPAINV CRS',': Total   time (s) = ',omp_get_wtime()-t2
#endif

end subroutine

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
  write(sparse%unlog,'(1x,a,1x,i0)')__FILE__,__LINE__
#endif
  return
 endif

 associate(parvar=>sparse%pardisovar)

 if(.not.allocated(parvar%pt))return

 sparse%loriginal=.true.

 nrhs=1

 !Reset Pardiso memory
 parvar%phase=-1
 call pardiso(parvar%pt,parvar%maxfct,parvar%mnum,parvar%mtype,parvar%phase,&
              sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,&
              idummy,nrhs,parvar%iparm,parvar%msglvl,parvar%ddum,parvar%ddum,error)
 call checkpardiso(parvar%phase,error,sparse%unlog)

 end associate

 if(allocated(sparse%pardisovar%pt))deallocate(sparse%pardisovar%pt)

end subroutine
#endif

!**PRINT
module subroutine print_crs(sparse,lint,output)
 class(crssparse),intent(in)::sparse
 integer(kind=int32),intent(in),optional::output
 logical,intent(in),optional::lint

 integer(kind=int32)::i
 integer(kind=int32)::un
 integer(kind=int32)::j
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

module subroutine print_idx_crs(sparse,lidx,lint,output)
 class(crssparse),intent(in)::sparse
 logical,intent(in)::lidx(:)
 integer(kind=int32),intent(in),optional::output
 logical,intent(in),optional::lint

 integer(kind=int32)::i
 integer(kind=int32)::un
 integer(kind=int32)::j
 character(len=30),parameter::frm='(2(i0,1x),g0)'
 logical::linternal

 if(size(lidx).lt.sparse%getdim(1).or.size(lidx).lt.sparse%getdim(2))then
  error stop ' ERROR: the vector lidx is smaller than the smallest dimension of CRS'
 endif

 linternal=.true.
 if(present(lint))linternal=lint

 un=sparse%unlog
 if(present(output))un=output

 do i=1,sparse%dim1
  if(.not.lidx(i))cycle
  do j=sparse%ia(i),sparse%ia(i+1)-1
   if(.not.lidx(sparse%ja(j)))cycle
   write(un,frm)i,sparse%ja(j),sparse%a(j)
   if(.not.linternal.and.sparse%lupperstorage.and.sparse%lsymmetric.and.i.ne.sparse%ja(j))then
    write(un,frm)sparse%ja(j),i,sparse%a(j)
   endif
  enddo
 enddo

end subroutine

module subroutine printsquare_crs(sparse,output)
 class(crssparse),intent(inout)::sparse
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
module subroutine save_crs(sparse,namefile)
 class(crssparse),intent(in)::sparse
 character(len=*),intent(in)::namefile

 integer(kind=int32)::un

 open(newunit=un,file=namefile,action='write',status='replace',access='stream')!,buffered='yes')
 write(un)typecrs                !int32
 write(un)sparse%dim1            !int32
 write(un)sparse%dim2            !int32
 write(un)sparse%nonzero()       !int64
 write(un)sparse%lupperstorage   !logical
 write(un)sparse%ia              !int32
 write(un)sparse%ja              !int32
 write(un)sparse%a               !wp
 close(un)

end subroutine

!**SCALE ALL ENTRIES
module subroutine scale_crs(sparse,val)
 class(crssparse),intent(inout)::sparse
 real(kind=wp),intent(in)::val
 sparse%a = sparse%a * val
end subroutine

!**SET ELEMENTS
module subroutine set_crs(sparse,row,col,val,error)
 !add a value only to an existing one
 class(crssparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::row,col
 !added: error=0;Not existing: error=-1;matrix not initialized: error=-10
 integer(kind=int32),intent(out),optional::error
 real(kind=wp),intent(in)::val

 integer(kind=int32)::i

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
module subroutine solve_crs_vector(sparse,x,y,msglvl)
 !sparse*x=y
 class(crssparse),intent(inout)::sparse
 real(kind=wp),intent(out),contiguous::x(:)
 real(kind=wp),intent(inout),contiguous::y(:)
 integer(kind=int32),intent(in),optional::msglvl

 !Pardiso variables
 integer(kind=int32)::error
 integer(kind=int32)::nrhs
 integer(kind=int32)::msglvl_opt

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

 if(sparse%loriginal)then
  !$ t1=omp_get_wtime()
  !Sort the matrix
  call sparse%sort()

  !Preparation of Cholesky of A with Pardiso
  !parvar=pardiso_variable(maxfct=1,mnum=1,mtype=-2,solver=0,msglvl=1)   !mtype=11
  !to avoid ifort 2021.4.0 bug
  sparse%pardisovar=pardiso_variable(maxfct=1,mnum=1,mtype=-2,solver=0,msglvl=msglvl_opt)   !mtype=11

  !initialize iparm
  call pardisoinit(sparse%pardisovar%pt,sparse%pardisovar%mtype,sparse%pardisovar%iparm)
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
  sparse%pardisovar%iparm(24)=1
  sparse%pardisovar%iparm(27)=1
#if (_DP==0)
  sparse%pardisovar%iparm(28)=1
#else
  sparse%pardisovar%iparm(28)=0
#endif
  write(sparse%unlog,'(a)')' Start ordering and factorization'
  if(allocated(sparse%perm))then
   sparse%pardisovar%iparm(5)=1
   sparse%pardisovar%iparm(31)=0
   sparse%pardisovar%iparm(36)=0
   call pardiso(sparse%pardisovar%pt,sparse%pardisovar%maxfct,sparse%pardisovar%mnum,&
                sparse%pardisovar%mtype,sparse%pardisovar%phase,&
                sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,&
                sparse%perm,nrhs,sparse%pardisovar%iparm,sparse%pardisovar%msglvl,&
                sparse%pardisovar%ddum,sparse%pardisovar%ddum,error)
  else
   call pardiso(sparse%pardisovar%pt,sparse%pardisovar%maxfct,sparse%pardisovar%mnum,&
                sparse%pardisovar%mtype,sparse%pardisovar%phase,&
                sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,&
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
  sparse%loriginal=.false.

 endif

 !Solving
 sparse%pardisovar%phase=33
 call pardiso(sparse%pardisovar%pt,sparse%pardisovar%maxfct,sparse%pardisovar%mnum,&
              sparse%pardisovar%mtype,sparse%pardisovar%phase,&
              sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,&
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

module subroutine solve_crs_array(sparse,x,y,msglvl)
 !sparse*x=y
 class(crssparse),intent(inout)::sparse
 real(kind=wp),intent(out),contiguous::x(:,:)
 real(kind=wp),intent(inout),contiguous::y(:,:)
 integer(kind=int32),intent(in),optional::msglvl

 !Pardiso variables
 integer(kind=int32)::error
 integer(kind=int32)::nrhs
 integer(kind=int32)::msglvl_opt

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

 if(sparse%loriginal)then
  !$ t1=omp_get_wtime()
  !Sort the matrix
  call sparse%sort()

  !Preparation of Cholesky of A with Pardiso
  !parvar=pardiso_variable(maxfct=1,mnum=1,mtype=-2,solver=0,msglvl=1)   !mtype=11
  !to avoid ifort 2021.4.0 bug
  sparse%pardisovar=pardiso_variable(maxfct=1,mnum=1,mtype=-2,solver=0,msglvl=msglvl_opt)   !mtype=11

  !initialize iparm
  call pardisoinit(sparse%pardisovar%pt,sparse%pardisovar%mtype,sparse%pardisovar%iparm)
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
  sparse%pardisovar%iparm(24)=1
  sparse%pardisovar%iparm(27)=1
#if (_DP==0)
  sparse%pardisovar%iparm(28)=1
#else
  sparse%pardisovar%iparm(28)=0
#endif
  write(sparse%unlog,'(a)')' Start ordering and factorization'
  if(allocated(sparse%perm))then
   sparse%pardisovar%iparm(5)=1
   sparse%pardisovar%iparm(31)=0
   sparse%pardisovar%iparm(36)=0
   call pardiso(sparse%pardisovar%pt,sparse%pardisovar%maxfct,sparse%pardisovar%mnum,&
                sparse%pardisovar%mtype,sparse%pardisovar%phase,&
                sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,&
                sparse%perm,nrhs,sparse%pardisovar%iparm,sparse%pardisovar%msglvl,&
                sparse%pardisovar%ddum,sparse%pardisovar%ddum,error)
  else
   call pardiso(sparse%pardisovar%pt,sparse%pardisovar%maxfct,sparse%pardisovar%mnum,&
                sparse%pardisovar%mtype,sparse%pardisovar%phase,&
                sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,&
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
  sparse%loriginal=.false.

 endif

 !Solving
 sparse%pardisovar%phase=33
 call pardiso(sparse%pardisovar%pt,sparse%pardisovar%maxfct,sparse%pardisovar%mnum,&
              sparse%pardisovar%mtype,sparse%pardisovar%phase,&
              sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,&
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
module subroutine solve_crs_vector(sparse,x,y,msglvl)
 !sparse*x=y
 class(crssparse),intent(inout)::sparse
 real(kind=wp),intent(out),contiguous::x(:)
 real(kind=wp),intent(inout),contiguous::y(:)
 integer(kind=int32),intent(in),optional::msglvl

 write(sparse%unlog,'(a)')' Warning: Pardiso is not enabled! Array returned = rhs'
 x=y

end subroutine

module subroutine solve_crs_array(sparse,x,y,msglvl)
 !sparse*x=y
 class(crssparse),intent(inout)::sparse
 real(kind=wp),intent(out),contiguous::x(:,:)
 real(kind=wp),intent(inout),contiguous::y(:,:)
 integer(kind=int32),intent(in),optional::msglvl

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
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'ISOLVE CRS y permutation',': Elapsed time (s) = ',omp_get_wtime()-t2
 !$ t2=omp_get_wtime()
#endif

#if(_DP==0)
 call mkl_scsrtrsv('U','T','N',sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,x_,x)
#else
 call mkl_dcsrtrsv('U','T','N',sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,x_,x)
#endif

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'ISOLVE CRS 1st triangular solve',': Elapsed time (s) = ',omp_get_wtime()-t2
 !$ t2=omp_get_wtime()
#endif

#if(_DP==0)
 call mkl_scsrtrsv('U','N','N',sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,x,x_)
#else
 call mkl_dcsrtrsv('U','N','N',sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,x,x_)
#endif

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'ISOLVE CRS 2nd triangular solve',': Elapsed time (s) = ',omp_get_wtime()-t2
 !$ t2=omp_get_wtime()
#endif

 do i=1,size(x)
  x(sparse%perm(i))=x_(i)
 enddo

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'ISOLVE CRS x permutation',': Elapsed time (s) = ',omp_get_wtime()-t2
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'ISOLVE CRS',': Total   time (s) = ',omp_get_wtime()-t1
#endif

end subroutine

!**SOLVE WITH PRE-COMPUTED LDLt DECOMPOSITION
module subroutine solveldlt_s_crs(sparse,x,y)
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
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'SOLVE LDLt CRS y permutation',': Elapsed time (s) = ',omp_get_wtime()-t2
 !$ t2=omp_get_wtime()
#endif

#if(_DP==0)
 call mkl_scsrtrsv('U','T','U',sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,x_,x)
#else
 call mkl_dcsrtrsv('U','T','U',sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,x_,x)
#endif
#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'SOLVE LDLt CRS 1st triangular solve',': Elapsed time (s) = ',omp_get_wtime()-t2
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
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'SOLVE LDLt CRS Diagonal solve',': Elapsed time (s) = ',omp_get_wtime()-t2
 !$ t2=omp_get_wtime()
#endif

#if(_DP==0)
 call mkl_scsrtrsv('U','N','U',sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,x,x_)
#else
 call mkl_dcsrtrsv('U','N','U',sparse%getdim(1),sparse%a,sparse%ia,sparse%ja,x,x_)
#endif

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'SOLVE LDLt CRS 2nd triangular solve',': Elapsed time (s) = ',omp_get_wtime()-t2
 !$ t2=omp_get_wtime()
#endif

 do i=1,size(x)
  x(sparse%perm(i))=x_(i)
 enddo

#if (_VERBOSE>0)
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'SOLVE LDLt CRS x permutation',': Elapsed time (s) = ',omp_get_wtime()-t2
 !$ write(sparse%unlog,'(1x,a,t30,a,f0.5)')'SOLVE LDLt CRS',': Total   time (s) = ',omp_get_wtime()-t1
#endif

end subroutine

!**SOLVE WITH LDLt DECOMPOSITION (AND COMPUTE IT IF NEEDED)
module subroutine solveldlt_crs_vector(sparse,x,y)
 !sparse*x=y
 class(crssparse),intent(inout)::sparse
 real(kind=wp),intent(out),contiguous::x(:)
 real(kind=wp),intent(inout),contiguous::y(:)

 !$ real(kind=real64)::t1

 if(.not.sparse%issquare())then
  write(sparse%unlog,'(a)')' Warning: the sparse matrix is not squared!'
#if (_VERBOSE>0)
  write(sparse%unlog,'(1x,a,1x,i0)')__FILE__,__LINE__
#endif
  return
 endif

 if(sparse%loriginal)then
  !$ t1=omp_get_wtime()
  !Sort the matrix
  call sparse%sort()

  !Ordering and factorization
  write(sparse%unlog,'(a)')' Start ordering and factorization'

 !Preparation of Cholesky of A
!   call sparse%chol()
 call sparse%getldlt()
 call sparse%sort()

  !$ write(sparse%unlog,'(a,f0.5)')' Elapsed time (s)               = ',omp_get_wtime()-t1

  sparse%loriginal=.false.

 endif

 !Solving
!  call sparse%isolve(x(:),y(:))
 call sparse%solveldlt_s(x(:),y(:))

end subroutine

module subroutine solveldlt_crs_array(sparse,x,y)
 !sparse*x=y
 class(crssparse),intent(inout)::sparse
 real(kind=wp),intent(out),contiguous::x(:,:)
 real(kind=wp),intent(inout),contiguous::y(:,:)

 integer(kind=int32)::i
 !$ real(kind=real64)::t1

 if(.not.sparse%issquare())then
  write(sparse%unlog,'(a)')' Warning: the sparse matrix is not squared!'
#if (_VERBOSE>0)
  write(sparse%unlog,'(1x,a,1x,i0)')__FILE__,__LINE__
#endif
  return
 endif

 if(size(x,2).ne.size(y,2))then
  write(sparse%unlog,'(a)')' ERROR: the number of colums of x and y provided to solve are different!'
  error stop
 endif

 if(sparse%loriginal)then
  !$ t1=omp_get_wtime()
  !Sort the matrix
  call sparse%sort()

  !Ordering and factorization
  write(sparse%unlog,'(a)')' Start ordering and factorization'

 !Preparation of Cholesky of A
!   call sparse%chol()
 call sparse%getldlt()
 call sparse%sort()

  !$ write(sparse%unlog,'(a,f0.5)')' Elapsed time (s)               = ',omp_get_wtime()-t1

  sparse%loriginal=.false.

 endif

 !Solving
 do i=1,size(x,2)
!   call sparse%isolve(x(:,i),y(:,i))
  call sparse%solveldlt_s(x(:,i),y(:,i))
 end do

end subroutine

!**SORT ARRAY
module subroutine sort_crs(sparse)
 ! sort vectors ja and a by increasing order
 class(crssparse),intent(inout)::sparse

 integer(kind=int32)::k
 integer(kind=int32)::endd,i,j,n,start,stkpnt
 integer(kind=int32)::d1,d2,d3,dmnmx,tmp
 integer(kind=int32)::stack(2,32)
 integer(kind=int32),allocatable::d(:)
 integer(kind=int32),parameter::select=20
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

 ll=(.not.sparse%loriginal)

end function

pure module subroutine setdecomposed_crs(sparse,ll)
 class(crssparse),intent(inout)::sparse
 logical,intent(in)::ll

 sparse%loriginal=(.not.ll)

end subroutine

end submodule

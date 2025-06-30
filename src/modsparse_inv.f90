!> Module containing subroutines for sparse inverses of symmetric positive definite matrices

!> @todo Still raw and not very efficient

module modsparse_inv
#if (_DP==0)
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real32
#else
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real64
#endif
 use modspainv, only: super_nodes, super_gsfct, super_sparsinv
 !$ use omp_lib
 implicit none
 private
 public::get_chol,get_ichol,get_spainv

 integer(kind=int32),parameter::minsizesupernode=0  !values other than 0 (e.g., 256) may give troubles

 interface get_chol
  module procedure get_chol_crs
 end interface

 interface get_ichol
  module procedure get_ichol_crs
 end interface

 interface get_spainv
  module procedure get_spainv_crs
 end interface

 interface
  subroutine smbfct(neqns,xadj,adjncy,perm,invp,&
              xlnz,maxlnz,xnzsub,nzsub,maxsub,&
              rchlnk,mrglnk,marker,flag&
              )
   integer::neqns
   integer::maxsub,flag,maxlnz
   integer::xadj(*),adjncy(*)
   integer::perm(*),invp(*)
   integer::xlnz(*),xnzsub(*),nzsub(*),rchlnk(*),mrglnk(*),marker(*)
  end subroutine
 end interface

contains

!PUBLIC
!Should return the Cholesky factor in the permuted order
subroutine get_ichol_crs(ia,ja,a,xadj,adjncy,perm,minsizenode,un)
 integer(kind=int32),intent(inout)::ia(:)
 integer(kind=int32),intent(inout)::ja(:)
 integer(kind=int32),intent(in)::xadj(:),adjncy(:)
 integer(kind=int32),intent(in)::perm(:)  !Ap(i,:)=A(perm(i),:)
 integer(kind=int32),intent(in),optional::minsizenode
 integer(kind=int32),intent(in),optional::un
 real(kind=wp),intent(inout)::a(:)
 
 integer(kind=int32)::unlog,neqns
 integer(kind=int32)::mssn
 integer(kind=int32),allocatable::xlnz(:),xnzsub(:),nzsub(:)
 real(kind=wp),allocatable::xspars(:),diag(:)
 !$ real(kind=real64)::t1
 real(kind=real64)::time(6)

 mssn=minsizesupernode
 if(present(minsizenode))mssn=minsizenode

 unlog=output_unit
 if(present(un))unlog=un

 time=0._real64

 neqns=size(ia)-1

 call get_ichol_spainv_crs(neqns,ia,ja,a,xadj,adjncy,perm,.false.,xlnz,xspars,xnzsub,nzsub,diag,mssn,time)

 ! Cholesky factorization
 !Convert to ija
 !$ t1=omp_get_wtime()
 call converttoija_noperm(neqns,xlnz,xspars,xnzsub,nzsub,diag,ia,ja,a,perm)
 !$ time(6)=omp_get_wtime()-t1

 call writetime(unlog,time,'SPARSE CHOL FACT.')

end subroutine

subroutine get_chol_crs(ia,ja,a,xadj,adjncy,perm,minsizenode,un)
 integer(kind=int32),intent(inout),allocatable::ia(:)
 integer(kind=int32),intent(inout),allocatable::ja(:)
 integer(kind=int32),intent(in)::xadj(:),adjncy(:)
 integer(kind=int32),intent(in)::perm(:)  !Ap(i,:)=A(perm(i),:)
 integer(kind=int32),intent(in),optional::minsizenode
 integer(kind=int32),intent(in),optional::un
 real(kind=wp),intent(inout),allocatable::a(:)
 
 integer(kind=int32)::unlog,neqns
 integer(kind=int32)::mssn
 integer(kind=int32),allocatable::xlnz(:),xnzsub(:),nzsub(:)
 real(kind=wp),allocatable::xspars(:),diag(:)
 !$ real(kind=real64)::t1
 real(kind=real64)::time(6)

 mssn=minsizesupernode
 if(present(minsizenode))mssn=minsizenode

 unlog=output_unit
 if(present(un))unlog=un

 time=0._real64

 neqns=size(ia)-1

 call get_ichol_spainv_crs(neqns,ia,ja,a,xadj,adjncy,perm,.false.,xlnz,xspars,xnzsub,nzsub,diag,mssn,time)

 !Convert to ija
 !$ t1=omp_get_wtime()
 call convertfactortoija(neqns,xlnz,xspars,xnzsub,nzsub,diag,ia,ja,a)
 !$ time(6)=omp_get_wtime()-t1

 call writetime(unlog,time,'CHOL FACT.')

end subroutine

subroutine get_spainv_crs(ia,ja,a,xadj,adjncy,perm,minsizenode,un)
 integer(kind=int32),intent(in)::ia(:)
 integer(kind=int32),intent(in)::ja(:)
 integer(kind=int32),intent(in)::xadj(:),adjncy(:)
 integer(kind=int32),intent(inout)::perm(:)  !Ap(i,:)=A(perm(i),:)
 integer(kind=int32),intent(in),optional::minsizenode
 integer(kind=int32),intent(in),optional::un
 real(kind=wp),intent(inout)::a(:)
 
 integer(kind=int32)::unlog,neqns
 integer(kind=int32)::mssn
 integer(kind=int32),allocatable::xlnz(:),xnzsub(:),nzsub(:)
 real(kind=wp),allocatable::xspars(:),diag(:)
 !$ real(kind=real64)::t1
 real(kind=real64)::time(6)

 mssn=minsizesupernode
 if(present(minsizenode))mssn=minsizenode

 unlog=output_unit
 if(present(un))unlog=un

 time=0._real64

 neqns=size(ia)-1

 call get_ichol_spainv_crs(neqns,ia,ja,a,xadj,adjncy,perm,.true.,xlnz,xspars,xnzsub,nzsub,diag,mssn,time)

 !Convert to ija
 !$ t1=omp_get_wtime()
 call converttoija(neqns,xlnz,xspars,xnzsub,nzsub,diag,ia,ja,a,perm)
 !$ time(6)=omp_get_wtime()-t1

 call writetime(unlog,time,'INVERSION')

end subroutine

subroutine get_ichol_spainv_crs(neqns,ia,ja,a,xadj,adjncy,perm,lspainv,xlnz,xspars,xnzsub,nzsub,diag,mssn,time)
 integer(kind=int32),intent(in)::neqns
 integer(kind=int32),intent(inout)::mssn
 integer(kind=int32),intent(in)::ia(:)
 integer(kind=int32),intent(in)::ja(:)
 integer(kind=int32),intent(in)::xadj(:),adjncy(:)
 integer(kind=int32),intent(in)::perm(:)  !Ap(i,:)=A(perm(i),:)
 real(kind=wp),intent(inout)::a(:)
 real(kind=real64),intent(inout)::time(:)
 logical,intent(in)::lspainv

 integer(kind=int32),allocatable,intent(out)::xlnz(:),xnzsub(:),nzsub(:)
 real(kind=wp),allocatable,intent(out)::xspars(:),diag(:)
  
#if (_VERBOSE >2)
 integer(kind=int32)::i
#endif
 integer(kind=int32)::nnode
 integer(kind=int32)::maxnode
 integer(kind=int32)::maxsub,flag,maxlnz
 integer(kind=int32)::rank
 integer(kind=int32),allocatable::inode(:)
!$ real(kind=real64)::t1

 !symbolic factorization
 !$ t1=omp_get_wtime()
 call symbolicfact(neqns,ia(neqns+1)-1,xadj,adjncy,perm,xlnz,maxlnz,xnzsub,nzsub,maxsub,flag)
 !$ time(1)=omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#if (_VERBOSE >0)
 !$ write(*,'(1x,a,t30,a,g0)')'CRS Symb. fact.',': Elapsed time = ',time(1)
#endif

 call computexsparsdiag(neqns,ia,ja,a,xlnz,nzsub,xnzsub,maxlnz,xspars,diag,perm)
 !$ time(2)=omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#if (_VERBOSE >0)
 !$ write(*,'(1x,a,t30,a,g0)')'CRS Compute xspars',': Elapsed time = ',time(2)
#endif

 allocate(inode(neqns+1))  !to allow diagonal matrices (for which the number of super-nodes is equal to the number of equations
 call super_nodes(mssn,neqns,xlnz,xnzsub,nzsub,nnode,inode,maxnode)
 !$ time(3)=omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#if (_VERBOSE >0)
 !$ write(*,'(1x,a,t30,a,g0)')'CRS Super nodes',': Elapsed time = ',time(3)
#endif

 ! Cholesky factorization
 call super_gsfct(neqns,xlnz,xspars,xnzsub,nzsub,diag,nnode,inode,rank)
 !$ time(4)=omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#if (_VERBOSE >0)
 !$ write(*,'(1x,a,t30,a,g0)')'CRS Chol. fact.',': Elapsed time = ',time(4)
#endif

 if(lspainv)then
  ! Matrix inverse
  call super_sparsinv(neqns,xlnz,xspars,xnzsub,nzsub,diag,nnode,inode)
  !$ time(5)=omp_get_wtime()-t1
  !$ t1=omp_get_wtime()
#if (_VERBOSE >0)
  !$ write(*,'(1x,a,t30,a,g0)')'CRS Inversion',': Elapsed time = ',time(5)
#endif
 endif
 
#if (_VERBOSE >0)
 write(*,'(/2x,a,i0)')'Flag symbolic factorization    : ',flag
 !write(*,'(2x,a,i0)')'Number of non-zero in the facor: ',maxlnz
 write(*,'(2x,a,i0)')'Number of super-nodes          : ',nnode
 write(*,'(2x,a,i0)')'Min size of super-nodes        : ',mssn
 write(*,'(2x,a,i0)')'Max size of super-nodes        : ',maxnode
 write(*,'(2x,a,i0/)')'Rank of the matrix             : ',rank
#endif
#if (_VERBOSE >2)
 do i=1,nnode
  write(*,'(1x,4(a,i8))') 'node',i,' from row', inode(i+1)+1,'  to row', inode(i),' size ',inode(i)-inode(i+1)
 end do
#endif

end subroutine

!PRIVATE
subroutine symbolicfact(neqns,nnzeros,xadj,adjncy,perm,xlnz,maxlnz,xnzsub,nzsub,maxsub,flag)
 integer(kind=int32),intent(in)::neqns,nnzeros
 integer(kind=int32),intent(in)::xadj(:),adjncy(:),perm(:)
 integer(kind=int32),intent(out)::maxlnz,maxsub,flag
 integer(kind=int32),allocatable,intent(out)::xlnz(:),xnzsub(:),nzsub(:)

 integer(kind=int32)::i,maxsubinit
 integer(kind=int32),allocatable::invp(:),rchlnk(:),mrglnk(:),marker(:)

 allocate(invp(size(perm)))
 do i=1,size(perm)
  invp(perm(i))=i
 enddo

 allocate(xlnz(neqns+1),xnzsub(neqns+1))
 allocate(rchlnk(neqns),mrglnk(neqns),marker(neqns))
 xlnz=0
 xnzsub=0
 rchlnk=0
 mrglnk=0
 marker=0
 !maxsub=ia(neqns+1)-1
 maxsub=nnzeros

 do
  flag=0
  maxsubinit=maxsub
  if(allocated(nzsub))deallocate(nzsub)
  allocate(nzsub(maxsub))
  nzsub=0
  call smbfct(neqns,xadj,adjncy,perm,invp,xlnz,maxlnz,xnzsub,nzsub,maxsub,&
              rchlnk,mrglnk,marker,flag&
              )
  if(maxsub.ne.maxsubinit)cycle
  if(flag.eq.0)exit
  maxsub=maxsub*2
 enddo
 deallocate(rchlnk,mrglnk,marker)

end subroutine

pure subroutine computexsparsdiag(neqns,ia,ja,a,xlnz,nzsub,xnzsub,maxlnz,xspars,diag,perm)
 integer(kind=int32),intent(in)::neqns,maxlnz
 integer(kind=int32),intent(in)::ia(:),ja(:)
 integer(kind=int32),intent(in)::perm(:)
 integer(kind=int32),intent(in)::xlnz(:),nzsub(:),xnzsub(:)
 real(kind=wp),intent(in)::a(:)
 real(kind=wp),intent(out),allocatable::xspars(:),diag(:)

 integer(kind=int32)::irow,iirow,icol,i,j,k

 allocate(xspars(maxlnz),diag(neqns))
 xspars=0._wp
 diag=0._wp

 do i=1,neqns
  irow=perm(i)
  diag(i)=a(ia(irow))
  do k=xlnz(i),xlnz(i+1)-1
   iirow=irow
   icol=perm(nzsub(xnzsub(i)+k-xlnz(i)))
   if(iirow.gt.icol)then
    iirow=icol
    icol=irow
   endif
   do j=ia(iirow)+1,ia(iirow+1)-1
    if(ja(j).eq.icol)then
     xspars(k)=a(j)
     exit
    endif
   enddo
  enddo
 enddo

end subroutine

pure subroutine converttoija(neqns,xlnz,xspars,xnzsub,ixsub,diag,ia,ja,a,perm)
 integer(kind=int32),intent(in)::neqns
 integer(kind=int32),intent(in)::ixsub(:),xlnz(:),xnzsub(:)
 integer(kind=int32),intent(in)::ia(:),ja(:),perm(:)
 real(kind=wp),intent(in)::xspars(:),diag(:)
 real(kind=wp),intent(inout):: a(:)

 integer(kind=int32)::irow,ksub,i,icol
 integer(kind=int32)::pirow,ppirow,picol,ip

 do irow = 1, neqns
  pirow=perm(irow)
  a(ia(pirow))=diag(irow)
  ksub = xnzsub(irow)
  do i = xlnz(irow), xlnz(irow+1)-1
   ppirow=pirow
   icol = ixsub(ksub)
   picol=perm(icol)
   ksub = ksub + 1
   if(ppirow.gt.picol)then
    ppirow=picol
    picol=pirow
   endif
   intloop: do ip=ia(ppirow)+1,ia(ppirow+1)-1
    if(ja(ip).eq.picol)then
     a(ip)=xspars(i)
     exit intloop
    endif
   enddo intloop
  end do
 end do

end subroutine 

pure subroutine converttoija_noperm(neqns,xlnz,xspars,xnzsub,ixsub,diag,ia,ja,a,perm)
 integer(kind=int32),intent(in)::neqns
 integer(kind=int32),intent(in)::ixsub(:),xlnz(:),xnzsub(:)
 integer(kind=int32),intent(in)::perm(:)
 integer(kind=int32),intent(inout)::ia(:),ja(:)
 real(kind=wp),intent(in)::xspars(:),diag(:)
 real(kind=wp),intent(out):: a(:)

 integer(kind=int32)::irow,ksub,i,icol
 integer(kind=int32)::pirow,ppirow,picol,ip
 integer(kind=int32)::nel
 integer(kind=int32),allocatable::iperm(:)
 integer(kind=int32),allocatable::iia(:),jja(:)

 
 !inefficient but it works
 
 !1. Permute ia and ja
 allocate(iperm,source=perm)
 do i=1,neqns
  iperm(perm(i))=i
 enddo
 nel=ia(neqns+1)-1
 allocate(iia(nel),jja(nel))
 
 do irow=1,neqns
  pirow=iperm(irow)
  do i=ia(irow),ia(irow+1)-1
   if(pirow.lt.iperm(ja(i)))then
    iia(i)=pirow
    jja(i)=iperm(ja(i))
   else
    iia(i)=iperm(ja(i))
    jja(i)=pirow
   endif
  enddo
 enddo

 ia=0
 ia(1)=1
 do i=1,nel
  ia(iia(i)+1)=ia(iia(i)+1)+1
 enddo

 ja=0
 do i=1,neqns
  ia(i+1)=ia(i+1)+ia(i)
  ja(ia(i))=i
 enddo

 iperm=1
 do i=1,nel
  if(iia(i).ne.jja(i))then
   ja(ia(iia(i))+iperm(iia(i)))=jja(i)
   iperm(iia(i))=iperm(iia(i))+1
  endif
 enddo

 !2. Replace a
 a=0._wp
 do irow = 1, neqns
  pirow=irow
  a(ia(pirow))=diag(irow)
  ksub = xnzsub(irow)
  do i = xlnz(irow), xlnz(irow+1)-1
   ppirow=pirow
   icol = ixsub(ksub)
   picol=icol
   ksub = ksub + 1
   if(ppirow.gt.picol)then
    ppirow=picol
    picol=pirow
   endif
   intloop: do ip=ia(ppirow)+1,ia(ppirow+1)-1
    if(ja(ip).eq.picol)then
     a(ip)=xspars(i)
     exit intloop
    endif
   enddo intloop
  end do
 end do

end subroutine 

pure subroutine convertfactortoija(neqns,xlnz,xspars,xnzsub,ixsub,diag,ia,ja,a)
 integer(kind=int32),intent(in)::neqns
 integer(kind=int32),intent(in)::ixsub(:),xlnz(:),xnzsub(:)
 integer(kind=int32),intent(inout)::ia(:)
 integer(kind=int32),intent(out),allocatable::ja(:)
 real(kind=wp),intent(in)::xspars(:),diag(:)
 real(kind=wp),intent(out),allocatable::a(:)

 integer(kind=int32)::j,irow,ksub,icol

 allocate(ja(size(xspars)+neqns), source=0)
 allocate(a(size(xspars)+neqns), source=0._wp)

 do irow=1,neqns
  ia(irow)=xlnz(irow)+irow-1
  ja(ia(irow))=irow
  a(ia(irow))=diag(irow)
  ksub=xnzsub(irow)
  do j=xlnz(irow),xlnz(irow+1)-1
   icol=ixsub(ksub)
   ksub=ksub+1
   ja(j+irow)=icol
   a(j+irow)=xspars(j)
  enddo
 enddo
 ia(neqns+1)=xlnz(neqns+1)+neqns

end subroutine 

subroutine writetime(unlog,time,a)
 integer(kind=int32)::unlog
 real(kind=real64),intent(in)::time(:)
 character(len=*),intent(in)::a

 write(unlog,'(/a)')' CRS MATRIX '//trim(a)
 !$ write(unlog,'(2x,a,t31,a,t33,f0.5,a)')'Symbolic factorization',':',time(1),' s'
 !$ write(unlog,'(2x,a,t31,a,t33,f0.5,a)')'Setup of tmp arrays',':',time(2),' s'
 !$ write(unlog,'(2x,a,t31,a,t33,f0.5,a)')'Node determination',':',time(3),' s'
 !$ write(unlog,'(2x,a,t31,a,t33,f0.5,a)')'Cholesky factorization',':',time(4),' s'
 !$ write(unlog,'(2x,a,t31,a,t33,f0.5,a)')'Matrix inversion',':',time(5),' s'
 !$ write(unlog,'(2x,a,t31,a,t33,f0.5,a)')'Conversion to CRS',':',time(6),' s'
 !$ write(unlog,'(2x,a,t31,a,t33,f0.5,a)')'Total time',':',sum(time),' s'

end subroutine

end module

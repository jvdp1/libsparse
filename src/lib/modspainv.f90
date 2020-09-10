!> Module containing subroutines for sparse inverses of symmetric positive definite matrices

!> @todo Still raw and not very efficient

!Based on Karin Meyer 's code !(didgeridoo.une.edu.au/womwiki/doku.php?id=fortran:fortran)
!Rewritten for my purposes

module modspainv
#if (_DP==0)
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real32
#else
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real64
#endif
 use modcommon
 !$ use omp_lib
 implicit none
 private
 public::get_chol,get_ichol,get_spainv

 integer(kind=int32),parameter::minsizesupernode=0  !values other than 0 (e.g., 256) may give troubles
 real(kind=wp),parameter::tol=1.e-10_wp

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
   integer,intent(in)::neqns
   integer::maxsub,flag,maxlnz
   integer,intent(in)::xadj(*),adjncy(*)
   integer,intent(in)::perm(*),invp(*)
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

subroutine get_chol_crs(ia,ja,a,xadj,adjncy,perm,minsizenode,lzerodiag,un)
 integer(kind=int32),intent(inout),allocatable::ia(:)
 integer(kind=int32),intent(inout),allocatable::ja(:)
 integer(kind=int32),intent(in)::xadj(:),adjncy(:)
 integer(kind=int32),intent(in)::perm(:)  !Ap(i,:)=A(perm(i),:)
 integer(kind=int32),intent(in),optional::minsizenode
 integer(kind=int32),intent(in),optional::un
 real(kind=wp),intent(inout),allocatable::a(:)
 logical,intent(in),optional::lzerodiag
 
 integer(kind=int32)::unlog,neqns
 integer(kind=int32)::mssn
 integer(kind=int32),allocatable::xlnz(:),xnzsub(:),nzsub(:)
 real(kind=wp),allocatable::xspars(:),diag(:)
 !$ real(kind=real64)::t1
 real(kind=real64)::time(6)
 logical::lzero

 mssn=minsizesupernode
 if(present(minsizenode))mssn=minsizenode

 unlog=output_unit
 if(present(un))unlog=un

 time=0._real64

 neqns=size(ia)-1

 call get_ichol_spainv_crs(neqns,ia,ja,a,xadj,adjncy,perm,.false.,xlnz,xspars,xnzsub,nzsub,diag,mssn,time)

 !Convert to ija
 !$ t1=omp_get_wtime()
 lzero=.true.
 if(present(lzerodiag))lzero=lzerodiag
print*,'wwww ',present(lzerodiag),lzero
 call convertfactortoija(neqns,xlnz,xspars,xnzsub,nzsub,diag,ia,ja,a,perm,lzero)
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
  
 integer(kind=int32)::i
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
 !$ write(*,'(x,a,t30,a,g0)')'CRS Symb. fact.',': Elapsed time = ',time(1)
#endif

 call computexsparsdiag(neqns,ia,ja,a,xlnz,nzsub,xnzsub,maxlnz,xspars,diag,perm)
 !$ time(2)=omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#if (_VERBOSE >0)
 !$ write(*,'(x,a,t30,a,g0)')'CRS Compute xspars',': Elapsed time = ',time(2)
#endif

 !super following Karin Meyer
 allocate(inode(neqns+1))  !to allow diagonal matrices (for which the number of super-nodes is equal to the number of equations
 call super_nodes(mssn,neqns,xlnz,xnzsub,nzsub,nnode,inode,maxnode)
 !$ time(3)=omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#if (_VERBOSE >0)
 !$ write(*,'(x,a,t30,a,g0)')'CRS Super nodes',': Elapsed time = ',time(3)
#endif

 ! Cholesky factorization
 call super_gsfct(neqns,xlnz,xspars,xnzsub,nzsub,diag,nnode,inode,rank)
 !$ time(4)=omp_get_wtime()-t1
 !$ t1=omp_get_wtime()
#if (_VERBOSE >0)
 !$ write(*,'(x,a,t30,a,g0)')'CRS Chol. fact.',': Elapsed time = ',time(4)
#endif

 if(lspainv)then
  ! Matrix inverse
  call super_sparsinv(neqns,xlnz,xspars,xnzsub,nzsub,diag,nnode,inode)
  !$ time(5)=omp_get_wtime()-t1
  !$ t1=omp_get_wtime()
#if (_VERBOSE >0)
  !$ write(*,'(x,a,t30,a,g0)')'CRS Inversion',': Elapsed time = ',time(5)
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
#if (_VERBOSE >4)
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

subroutine super_nodes(mssn, neqns, xlnz, xnzsub, ixsub, nnode, inode,maxnode)
 integer(kind=int32),intent(inout)::mssn
 integer(kind=int32),intent(in)::neqns
 integer(kind=int32),intent(in)::ixsub(:),xlnz(:),xnzsub(:)
 integer(kind=int32),intent(out)::nnode
 integer(kind=int32),intent(out)::maxnode
 integer(kind=int32),intent(out)::inode(:) !size=neqns+1

 integer(kind=int32)::i,ii,j,n,ilast,kk
 integer(kind=int32)::minnode
 real(kind=wp)::xx

 ! establish boundaries between diaggonal blocks 100% full
 ilast = neqns
 nnode = 0
 minnode=2**30
 maxnode=0
 do i = neqns-1, 1, -1
    ii = xnzsub(i)
    n = 0
    do j = xlnz(i), xlnz(i+1)-1
       kk = ixsub(ii)
       ii = ii + 1
       if(  kk > ilast ) exit
       n = n + 1
    end do
    xx = dble( n ) / dble( ilast - i )
    if( xx < 1._wp .and. (ilast-i).ge.mssn ) then
     nnode = nnode +1
     minnode=min(minnode,ilast-i)
     maxnode=max(maxnode,ilast-i)
     inode(nnode) = ilast
     ilast = i
    end if
 end do
 nnode = nnode +1
 inode(nnode) = ilast
 inode(nnode+1) = 0

 mssn=minnode

end subroutine 

subroutine super_gsfct(neqns,xlnz,xspars,xnzsub,ixsub,diag,nnode,inode,rank)
 integer(kind=int32),intent(in)::neqns,nnode
 integer(kind=int32),intent(in)::ixsub(:),xlnz(:),xnzsub(:),inode(:)
 integer(kind=int32),intent(out),optional::rank
 real(kind=wp),intent(inout)::xspars(:),diag(:)

 integer(kind=int32)::i,j,k,jrow,n,ksub,irow,jnode,icol1,icol2,jcol,ii,jj,mm,kk
#if (_VERBOSE>0)
 integer(kind=int32)::ninit
#endif
 integer(kind=int32)::orank
 integer(kind=int32),allocatable::jvec(:),kvec(:)
 real(kind=wp),allocatable::ttt(:,:),s21(:,:),s22(:,:)
! real(kind=wp),allocatable::ttt1(:,:)   !aaaa
 logical::lpos

 write(output_unit,'(/" Cholesky factorization...")')

 allocate(jvec(neqns),kvec(neqns),stat=ii)
 if(ii.ne.0)call alloc_err(__LINE__,__FILE__)

 jvec=0 !it must be re-initialized to 0 only if it is modified, i.e, when n.gt.0
#if (_VERBOSE>0)
 ninit=0
#endif
 orank=0

 do jnode = nnode, 1, -1
  icol1 = inode(jnode+1) + 1
  icol2 = inode(jnode)
  mm = icol2 - icol1 +1

  call progress(icol2,neqns)

  !pick out diagonal block
  allocate( ttt(icol1:icol2,icol1:icol2), stat = ii )
  if(ii.ne.0)call alloc_err(__LINE__,__FILE__)
  ttt=0._wp
  !jvec=0
  n=0
  do irow = icol1, icol2
   ttt(irow,irow) = diag(irow)
   ksub = xnzsub(irow)
   do i = xlnz(irow), xlnz(irow+1)-1
    jcol = ixsub(ksub)
    ksub = ksub + 1
    if( jcol <= icol2 ) then
     ttt(jcol,irow) = xspars(i) 
    else
     if(jvec(jcol).eq.0)then
      n=n+1
      jvec(jcol)=n
      kvec(n)=jcol
     endif
    end if
   end do
  end do

  !factorise

  lpos=.true.
  irow=mm
  if(irow.eq.1)then
   if(ttt(icol1,icol1).gt.0._wp)then
    ttt=sqrt(ttt)
   else
    lpos=.false.
    ttt=0._wp
    irow=0
   endif
  else
#if(_DP==0)
   call spotrf('L',mm,ttt,mm,ii)
#else
   call dpotrf('L',mm,ttt,mm,ii)
#endif
   if(ii.ne.0)then
    write(*,'(a,i0,a)')'Routine DPOTRF returned error code: ',ii,' (matrix is not positive definite)'
    error stop
   endif
  endif
  orank=orank+irow

  !adjust block below diagonal
  if(n.gt.0)then
         !... pick out rows
         allocate( s21(n, icol1:icol2), stat = ii )
         if(ii.ne.0)call alloc_err(__LINE__,__FILE__)
         s21 = 0._wp
         do irow = icol1, icol2
            ksub = xnzsub(irow)
            do i = xlnz(irow), xlnz(irow+1)-1
               jcol = ixsub(ksub)
               ksub = ksub + 1
               if( jcol <= icol2 ) cycle
               jj = jvec(jcol)
               s21(jj,irow) = xspars(i)
            end do
         end do

         !calculate L21
   if(lpos)then
#if(_DP==0)
    call strsm( 'R', 'L', 'T', 'N', n, mm, 1._wp, ttt, mm, s21, n )
   else
    call sgtrsm( 'R', 'L', 'T', 'N', n, mm, 1._wp, ttt, mm, s21, n )
#else
    call dtrsm( 'R', 'L', 'T', 'N', n, mm, 1._wp, ttt, mm, s21, n )
   else
    call dgtrsm( 'R', 'L', 'T', 'N', n, mm, 1._wp, ttt, mm, s21, n )
#endif
   endif
         !adjust remaining triangle to right: A22 := A22 - L21 L21'
         allocate( s22(n,n), stat = ii )
         if(ii.ne.0)call alloc_err(__LINE__,__FILE__)
#if(_DP==0)
         call ssyrk( 'L', 'N', n, mm, 1._wp, s21, n, 0._wp, s22, n )
#else
!         call dsyrk( 'L', 'N', n, mm, 1._wp, s21, n, 0._wp, s22, n )
block      !aaaaaaaaaaaaaa !to be checked
real(kind=wp)::tmp
         do j=1,n
          s22(:,j)=0._wp
          do ii=icol1,icol2
           if(s21(j,ii).ne.0._wp)then
            tmp=s21(j,ii)
            do i=j,n
             s22(i,j)=s22(i,j)+tmp*s21(i,ii)
            enddo
           endif
          enddo 
         enddo 
end block
#endif
         kk = maxval(kvec(1:n))
         do i = 1, n
             jrow = kvec(i)
             diag(jrow) = diag(jrow) - s22(i,i)
             ksub = xnzsub(jrow)
             do j = xlnz(jrow), xlnz(jrow+1)-1
                jcol = ixsub(ksub)
                if( jcol > kk ) exit
                ksub = ksub + 1
                ii = jvec(jcol)
                if( ii > 0 ) then
                 if(ii>i)then
                  xspars(j) = xspars(j) - s22(ii,i)
                 else
                  xspars(j) = xspars(j) - s22(i,ii)
                 endif
                endif
             end do
         end do
         deallocate( s22 )
     end if

     !transfer block back to sparse storage
     do irow = icol1, icol2 
       diag(irow) = ttt(irow,irow)
        ksub = xnzsub(irow)
        do i = xlnz(irow), xlnz(irow+1)-1
           jcol = ixsub(ksub)
           ksub = ksub + 1
           if( jcol <= icol2 ) then
               xspars(i) = ttt(jcol,irow)
           else
               xspars(i) = s21( jvec(jcol), irow)
           end if
        end do
     end do

  deallocate( ttt )
  if(n.gt.0)then
   deallocate(s21)
   !jvec=0
   !this should be faster than jvec=0
   do irow = icol1, icol2
    ksub = xnzsub(irow)
    do i = xlnz(irow), xlnz(irow+1)-1
     jcol = ixsub(ksub)
     ksub = ksub + 1
     if( jcol <= icol2 ) then
     else
       jvec(jcol)=0
     end if
    end do
   end do
#if (_VERBOSE>0)
   ninit=ninit+1
#endif
  endif

 enddo ! jnode

 !zeroing columns for positive semi-definite matrix
 do irow=neqns,1,-1
  if(diag(irow).lt.tol)then
   jvec(irow)=1
  endif
  ksub = xnzsub(irow)
  do i = xlnz(irow), xlnz(irow+1)-1
   jcol = ixsub(ksub)
   ksub = ksub + 1
   if(jvec(jcol).ne.0)xspars(i)=0._wp
  end do
 end do

 deallocate(jvec,kvec)
 
 if(present(rank))rank=orank

#if (_VERBOSE>0)
 write(*,'(a,i0)')' Number of zeroing a large vector (jvec): ',ninit
#endif
 
end subroutine

subroutine super_sparsinv(neqns,xlnz,xspars,xnzsub,ixsub,diag,nnode,inode)
  integer(kind=int32),intent(in)::neqns,nnode
  integer(kind=int32),intent(in)::ixsub(:),xlnz(:),xnzsub(:),inode(:)
  real(kind=wp),intent(inout) ::xspars(:),diag(:)

  integer(kind=int32)::irow,ksub,i,j,k,m,jcol,jrow, jnode, icol2, icol1, ii,jj, mm, n21, iopt
  integer(kind=int32),allocatable::kvec(:),jvec(:)
  real(kind=wp)::tt,xx
  real(kind=wp),dimension(:,:),allocatable:: ttt, s21, s22, f21
  real(kind=wp),dimension(:),allocatable:: rr, qx

  write(output_unit,'(/" Sparse inversion...")')

  allocate( jvec(neqns), kvec(neqns),stat = ii )
  if(ii.ne.0)call alloc_err(__LINE__,__FILE__)
 
  jvec=0  !placed here and at the end of the loop to avoid to re-initialize it at each iteration and when n21=0

! backwards flops: determine inverse using supernodal blocks
  do jnode = 1, nnode
     icol1 = inode(jnode+1) + 1
     icol2 = inode(jnode)
     mm = icol2 - icol1 +1

   call progress(icol2,neqns)
     
     !pick out diagonal block
     allocate( ttt(icol1:icol2,icol1:icol2), stat = ii )
     if(ii.ne.0)call alloc_err(__LINE__,__FILE__)
     ttt=0._wp
     !jvec=0
     n21=0
     do irow=icol1,icol2
      ttt(irow,irow)=diag(irow)
      ksub=xnzsub(irow)
      do i=xlnz(irow),xlnz(irow+1)-1
       jcol=ixsub(ksub)
       ksub=ksub+1
       if(jcol<=icol2)then
        ttt(jcol,irow)=xspars(i) 
       else 
        if(jvec(jcol).eq.0)then
         n21=n21+1
         jvec(jcol)=n21
         kvec(n21)=jcol
        endif
       end if
      end do
     end do

     !pick out lead columns (condensed)
     if( n21 > 0 ) then
         allocate( s21(n21, icol1:icol2), stat = ii )
         if(ii.ne.0)call alloc_err(__LINE__,__FILE__)
         allocate( f21(n21, icol1:icol2), stat = ii )
         if(ii.ne.0)call alloc_err(__LINE__,__FILE__)
         s21 = 0._wp
         f21 = 0._wp
         do irow = icol1, icol2
            ksub = xnzsub(irow)
            do i = xlnz(irow), xlnz(irow+1)-1
               jcol = ixsub(ksub)
               ksub = ksub + 1
               if( jcol <= icol2 ) cycle
               jj = jvec(jcol)
               s21(jj,irow) = xspars(i)
            end do
         end do
!        ... post-multiply with inverse Chol factor -> solve
#if(_DP==0)
         call strsm( 'R', 'L', 'N', 'N', n21, mm, 1._wp, ttt, mm, s21, n21 )
#else
         call dtrsm( 'R', 'L', 'N', 'N', n21, mm, 1._wp, ttt, mm, s21, n21 )
#endif
!        ... invert Cholesky factor
#if(_DP==0)
         call spotri( 'L',  mm, ttt, mm, ii )
#else
         call dpotri( 'L',  mm, ttt, mm, ii )
#endif
         if( ii /= 0 ) then
             write(*,*) 'Routine DPOTRI returned error code', ii
             stop
         end if
!        ... pre-multiply by already inverted submatrix
         iopt = 2
44       if( iopt == 1 ) then
            allocate( rr(icol1:icol2), qx(icol1:icol2), stat = ii )
            if(ii.ne.0)call alloc_err(__LINE__,__FILE__)
            f21 = 0._wp
            do k = 1, n21
               jrow = kvec(k)
               rr = s21(k,:)
               qx =  diag(jrow) * rr 
               ksub = xnzsub(jrow)
               do i = xlnz(jrow), xlnz(jrow+1)-1
                  jcol = ixsub(ksub)
                  ksub = ksub + 1
                  m = jvec(jcol)
                  if( m < 1 ) cycle
                  xx = xspars(i)
                  f21(m,:) = f21(m,:) + xx * rr
                  qx = qx + xx * s21(m,:) 
               end do
               f21(k,:) = f21(k,:) + qx 
            end do
            deallocate( rr, qx )
         else
            allocate( s22(n21,n21), stat = ii )
            if( ii /= 0 ) then
                iopt = 1
                go to 44
            end if
            s22 = 0._wp
            do k = 1, n21
               jrow = kvec(k)
               s22(k,k) = diag(jrow)
               ksub = xnzsub(jrow)
               do i = xlnz(jrow), xlnz(jrow+1)-1
                  jcol = ixsub(ksub)
                  ksub = ksub + 1
                  m = jvec(jcol)
                  if( m > 0 )then
                   if(m>k)then
                    s22(m,k) = xspars(i)
                   else
                    s22(k,m) = xspars(i)
                   endif
                  endif
               end do
            end do
#if(_DP==0)
            call ssymm( 'L', 'L', n21, mm, 1._wp, s22, n21, s21, n21, 0._wp,f21, n21  )
#else
            call dsymm( 'L', 'L', n21, mm, 1._wp, s22, n21, s21, n21, 0._wp,f21, n21  )
#endif
            deallocate( s22 )
         end if
!        ... adjustments to current block
#if(_DP==0)
         call sgemm( 'T', 'N', mm, mm, n21, 1._wp, f21, n21, s21, n21, 1._wp,ttt, mm )
#else
         call dgemm( 'T', 'N', mm, mm, n21, 1._wp, f21, n21, s21, n21, 1._wp,ttt, mm )
#endif
     else
#if(_DP==0)
         call spotri( 'L',  mm, ttt, mm, ii )
#else
         call dpotri( 'L',  mm, ttt, mm, ii )
#endif
         if( ii /= 0 ) then
             write(*,*) 'Routine DPOTRI returned error code', ii
             stop
         end if
     end if

!    save current block  
     do irow = icol1, icol2
        diag(irow) = ttt(irow,irow)
        ksub = xnzsub(irow)
        do i = xlnz(irow), xlnz(irow+1)-1
           jcol = ixsub(ksub)
           ksub = ksub + 1
           if( jcol <= icol2 ) then
               xspars(i) = ttt(jcol,irow) 
           else 
               xspars(i) = - f21(jvec(jcol),irow)
           end if
        end do
     end do

   deallocate( ttt )
   if(n21.gt.0)then
    deallocate(s21,f21)
    jvec=0
   endif

  enddo ! jnode

  deallocate(jvec,kvec)

end subroutine 

subroutine computexsparsdiag(neqns,ia,ja,a,xlnz,nzsub,xnzsub,maxlnz,xspars,diag,perm)
 integer(kind=int32),intent(in)::neqns,maxlnz
 integer(kind=int32),intent(in)::ia(:),ja(:)
 integer(kind=int32),intent(in)::perm(:)
 integer(kind=int32),intent(in)::xlnz(:),nzsub(:),xnzsub(:)
 real(kind=wp),intent(in)::a(:)
 real(kind=wp),intent(out),allocatable::xspars(:),diag(:)

 integer(kind=int32)::irow,iirow,icol,i,j,k,kk
 integer(kind=int32),allocatable::tmp(:)
 real(kind=wp),allocatable::rtmp(:)

 allocate(xspars(maxlnz),diag(neqns))
 xspars=0._wp
 diag=0._wp

 do i=1,neqns
  irow=perm(i)
  do j=ia(irow),ia(irow+1)-1
   if(ja(j).eq.irow)then
    diag(i)=a(j)
    exit
   endif
  enddo
  !diag(i)=a(ia(irow))
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

subroutine converttoija(neqns,xlnz,xspars,xnzsub,ixsub,diag,ia,ja,a,perm)
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

subroutine converttoija_noperm(neqns,xlnz,xspars,xnzsub,ixsub,diag,ia,ja,a,perm)
 integer(kind=int32),intent(in)::neqns
 integer(kind=int32),intent(in)::ixsub(:),xlnz(:),xnzsub(:)
 integer(kind=int32),intent(in)::perm(:)
 integer(kind=int32),intent(inout)::ia(:),ja(:)
 real(kind=wp),intent(in)::xspars(:),diag(:)
 real(kind=wp),intent(inout):: a(:)

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

subroutine convertfactortoija(neqns,xlnz,xspars,xnzsub,ixsub,diag,ia,ja,a,perm,lzerodiag)
 integer(kind=int32),intent(in)::neqns
 integer(kind=int32),intent(in)::ixsub(:),xlnz(:),xnzsub(:)
 integer(kind=int32),intent(inout)::ia(:)
 integer(kind=int32),intent(inout),allocatable::ja(:)
 integer(kind=int32),intent(in)::perm(:)
 real(kind=wp),intent(in)::xspars(:),diag(:)
 real(kind=wp),intent(inout),allocatable::a(:)
 logical,optional::lzerodiag

 integer(kind=int32)::ndiag
 integer(kind=int32)::j,irow,ksub,icol
 logical::lzero

 deallocate(ja)
 deallocate(a)

 lzero=.true.
 if(present(lzerodiag))lzero=lzerodiag

 ndiag=neqns
 if(.not.lzero)then
  lzero=.false.
  ndiag=count(diag.ne.0._wp)
 endif

 allocate(ja(size(xspars)+ndiag))
 allocate(a(size(xspars)+ndiag))
 ja=0
 a=0._wp

 ia(1:neqns+1)=xlnz(1:neqns+1)

 ndiag=0
 if(lzero)then
  do irow=1,neqns
   ia(irow)=ia(irow)+ndiag
   ndiag=ndiag+1
   ja(ia(irow))=irow
   a(ia(irow))=diag(irow)
  enddo
 else
  do irow=1,neqns
   ia(irow)=ia(irow)+ndiag
   if(diag(irow).ne.0._wp)then
    ndiag=ndiag+1
    ja(ia(irow))=irow
    a(ia(irow))=diag(irow)
   endif
  enddo
 endif

 ia(neqns+1)=xlnz(neqns+1)+ndiag

 ndiag=0
 do irow=1,neqns
  ksub=xnzsub(irow)
  if(lzero)then
   ndiag=ndiag+1
  else
   if(diag(irow).ne.0._wp)ndiag=ndiag+1
  endif
  do j=xlnz(irow),xlnz(irow+1)-1
   icol=ixsub(ksub)
   ksub=ksub+1
   ja(j+ndiag)=icol
   a(j+ndiag)=xspars(j)
  enddo
 enddo

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

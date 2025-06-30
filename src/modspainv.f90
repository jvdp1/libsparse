!> Module containing subroutines for sparse inverses of symmetric positive definite matrices

!Based on Karin Meyer 's code (didgeridoo.une.edu.au/womwiki/doku.php?id=fortran:fortran).
!The content of Karin Meyer 's wiki is licensed under CC Attribution-Share Alike 4.0 International.
!
!Support of sparse inverse of SPSD matrices by implementing the S. D. Kachman modifications
!(https://www.ars.usda.gov/ARSUserFiles/80420530/MTDFREML/MTDFMan.pdf ; Chapter 6)

!Rewritten for my purposes

module modspainv
#if (_DP==0)
 use iso_fortran_env,only:output_unit,int32,real32,real64,wp=>real32
 use modsparse_mkl, only: spotrf, spotri, strsm, ssymm, sgemm
#else
 use iso_fortran_env,only:output_unit,int32,real32,real64,wp=>real64
 use modsparse_mkl, only: dpotrf, dpotri, dtrsm, dsymm, dgemm
#endif
 use modcommon, only: progress
 implicit none(type, external)
 private
 public :: super_nodes, super_gsfct, super_sparsinv

 real(kind=wp),parameter::tol=1.e-8_wp

 interface
  SUBROUTINE dgtrsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
   !.. Scalar Arguments ..
   DOUBLE PRECISION ALPHA
   INTEGER LDA,LDB,M,N
   CHARACTER DIAG,SIDE,TRANSA,UPLO
   !.. Array Arguments ..
   DOUBLE PRECISION A(lda,*),B(ldb,*)
  end subroutine
 end interface

contains

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

 ! find boundaries between diagonal blocks 100% full
 ilast = neqns
 nnode = 0
 inode = 0
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
    xx = real( n, kind=wp ) / real( ilast - i, kind=wp )
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
 integer(kind=int32),intent(in)::neqns
 integer(kind=int32),intent(inout)::nnode
 integer(kind=int32),intent(in)::ixsub(:),xlnz(:),xnzsub(:)
 integer(kind=int32),intent(inout)::inode(:)
 integer(kind=int32),intent(out),optional::rank
 real(kind=wp),intent(inout)::xspars(:),diag(:)

 integer(kind=int32)::i,j,jrow,n,ksub,irow,jnode,icol1,icol2,jcol,ii,jj,mm,kk
#if (_VERBOSE>0)
 integer(kind=int32)::ninit
#endif
 integer(kind=int32)::orank
 integer(kind=int32) :: nnode_
 integer(kind=int32), allocatable :: inode_(:)
 integer(kind=int32),allocatable::jvec(:),kvec(:)
 real(kind=wp),allocatable::ttt(:,:),s21(:,:),s22(:,:)
! real(kind=wp),allocatable::ttt1(:,:)   !aaaa
 logical::lpos
 logical :: ldpotrf

 write(output_unit,'(/" Cholesky factorization...")')

 allocate(jvec(neqns),kvec(neqns))

 nnode_ = nnode
 allocate(inode_, source = inode)

 jvec=0 !it must be re-initialized to 0 only if it is modified, i.e, when n.gt.0
#if (_VERBOSE>0)
 ninit=0
#endif
 orank=0

 jnode = nnode_ + 1
 ldpotrf = .true.

 loopjnode: do !jnode = nnode, 1, -1
  if(.not.ldpotrf)then
   call split_inode(inode_, nnode_, jnode)
   ldpotrf = .true.
   jnode = jnode + 1
  endif

  jnode = jnode - 1
  if(jnode < 1)exit
  icol1 = inode_(jnode+1) + 1
  icol2 = inode_(jnode)
  mm = icol2 - icol1 +1

  !check all diagonal elements are > tol
  if(mm.gt.1)then
   do irow = icol1, icol2
    if(diag(irow).le.tol)then
     ldpotrf = .false.
     cycle loopjnode
    endif
   enddo
  endif

  call progress(icol2,neqns)

  !pick out diagonal block
  allocate( ttt(icol1:icol2,icol1:icol2))
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
   if(ttt(icol1,icol1).gt.tol)then
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
    ldpotrf = .false.
    deallocate(ttt)
    jvec = 0
    cycle loopjnode
!    error stop
   endif
  endif
  orank=orank+irow

  !adjust block below diagonal
  if(n.gt.0)then
         !... pick out rows
         allocate( s21(n, icol1:icol2))
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
         allocate( s22(n,n))
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

 enddo loopjnode ! jnode

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
 
 inode = inode_
 nnode = nnode_

 if(present(rank))rank=orank

#if (_VERBOSE>0)
 write(*,'(a,i0)')' Number of zeroing a large vector (jvec): ',ninit
#endif
 
end subroutine

subroutine super_sparsinv(neqns,xlnz,xspars,xnzsub,ixsub,diag,nnode,inode)
  integer(kind=int32),intent(in)::neqns,nnode
  integer(kind=int32),intent(in)::ixsub(:),xlnz(:),xnzsub(:),inode(:)
  real(kind=wp),intent(inout) ::xspars(:),diag(:)

  integer(kind=int32)::irow,ksub,i,k,m,jcol,jrow, jnode, icol2, icol1, ii,jj, mm, n21, iopt
  integer(kind=int32),allocatable::kvec(:),jvec(:)
  real(kind=wp)::xx
  real(kind=wp),dimension(:,:),allocatable:: ttt, s21, s22, f21
  real(kind=wp),dimension(:),allocatable:: rr, qx
  logical :: lpos

  write(output_unit,'(/" Sparse inversion...")')

  allocate( jvec(neqns), kvec(neqns))
 
  jvec=0  !placed here and at the end of the loop to avoid to re-initialize it at each iteration and when n21=0

! backwards flops: determine inverse using supernodal blocks
  do jnode = 1, nnode
     icol1 = inode(jnode+1) + 1
     icol2 = inode(jnode)
     mm = icol2 - icol1 +1

   call progress(icol2,neqns)
     
     !pick out diagonal block
     allocate( ttt(icol1:icol2,icol1:icol2))
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
         allocate( s21(n21, icol1:icol2))
         allocate( f21(n21, icol1:icol2))
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
         lpos = .not.(mm.eq.1.and.ttt(icol1,icol1).le.tol)

         if(lpos)then
#if(_DP==0)
         call strsm( 'R', 'L', 'N', 'N', n21, mm, 1._wp, ttt, mm, s21, n21 )
         else
         call sgrsm( 'R', 'L', 'N', 'N', n21, mm, 1._wp, ttt, mm, s21, n21 )
#else
         call dtrsm( 'R', 'L', 'N', 'N', n21, mm, 1._wp, ttt, mm, s21, n21 )
         else
         call dgtrsm( 'R', 'L', 'N', 'N', n21, mm, 1._wp, ttt, mm, s21, n21 )
#endif
         endif

!        ... invert Cholesky factor
         if(.not.lpos)then
             ttt=0._wp
         else
#if(_DP==0)
         call spotri( 'L',  mm, ttt, mm, ii )
#else
         call dpotri( 'L',  mm, ttt, mm, ii )
#endif
         if( ii /= 0 ) then
             write(*,'(a,i0)') 'Routine DPOTRI returned error code ', ii
             error stop
         end if
         endif
!        ... pre-multiply by already inverted submatrix
         iopt = 2
44       if( iopt == 1 ) then
            allocate( rr(icol1:icol2), qx(icol1:icol2))
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
         if(mm.eq.1)then
           if(ttt(icol1,icol1).le.tol)then
            ttt=0._wp
           else
            call dpotri( 'L',  mm, ttt, mm, ii )
           endif
         else
#if(_DP==0)
         call spotri( 'L',  mm, ttt, mm, ii )
#else
         call dpotri( 'L',  mm, ttt, mm, ii )
#endif
         if( ii /= 0 ) then
             write(*,'(a,i0)') 'Routine DPOTRI returned error code ', ii
             error stop
         end if
         endif
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

pure subroutine split_inode(inode, nnode, jnode)
 integer(kind=int32), intent(inout) :: inode(:)
 integer(kind=int32), intent(inout) :: nnode
 integer(kind=int32), intent(inout) :: jnode

 integer(kind=int32) :: i, n

 n = inode(jnode) - inode(jnode+1) - 1

 inode(n+jnode+1:n+nnode+1) = inode(jnode+1:nnode+1)

 inode(jnode+1:n+jnode) = [ (i, i = inode(jnode)-1, inode(jnode)-n, -1)]

 nnode = nnode + n

 jnode = jnode + n

end subroutine

end module

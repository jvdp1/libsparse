!> Module containing subroutines for sparse inverses of symmetric positive definite matrices

!> @todo Still raw and not very efficient

!Based on Karin Meyer 's code !(didgeridoo.une.edu.au/womwiki/doku.php?id=fortran:fortran)
!Slightly rewritten for my purposes

module modspainv
#if (_DP==0)
 use iso_fortran_env,only:int32,int64,real32,real64,wp=>real32
#else
 use iso_fortran_env,only:int32,int64,real32,real64,wp=>real64
#endif
 !$ use omp_lib
 implicit none
 private
 public::get_spainv

 interface get_spainv
  module procedure get_spainv_crs
 end interface

 interface
  subroutine smbfct(neqns,xadj,adjncy,perm,invp,&
              xlnz,maxlnz,xnzsub,nzsub,maxsub,&
              rchlnk,mrglnk,marker,flag&
              )
   import::int32
   integer(kind=int32),intent(in)::neqns
   integer(kind=int32)::maxsub,maxsubinit,flag,maxlnz
   integer(kind=int32),intent(in)::xadj(*),adjncy(*)
   integer(kind=int32),intent(in)::perm(*),invp(*)
   integer(kind=int32)::xlnz(*),xnzsub(*),nzsub(*),rchlnk(*),mrglnk(*),marker(*)
  end subroutine
 end interface

contains

!PUBLIC
subroutine get_spainv_crs(ia,ja,a,xadj,adjncy,perm)
 integer(kind=int32),intent(in)::ia(:)
 integer(kind=int32),intent(in)::ja(:)
 integer(kind=int32),intent(inout)::xadj(:),adjncy(:)
 integer(kind=int32),intent(inout)::perm(:)  !Ap(i,:)=A(perm(i),:)
 real(kind=wp),intent(inout)::a(:)
 
 integer(kind=int32)::i,neqns
 integer(kind=int32)::nnode
 integer(kind=int32)::maxsub,flag,maxlnz!,maxsubinit
 !integer(kind=int32),allocatable::invp(:)
 integer(kind=int32),allocatable::xlnz(:),xnzsub(:),nzsub(:)!,rchlnk(:),mrglnk(:),marker(:)
 integer(kind=int32),allocatable::inode(:)
 real(kind=wp),allocatable::xspars(:),diag(:)

 real(wp),dimension(size(ia)-1,size(ia)-1) :: hh, hh1, hh2, hh3

 print*,' start inverse aaaaaaaaaaaa'

 neqns=size(ia)-1

 !symbolic factorization

! allocate(invp,source=perm)
! 
! do i=1,size(perm)
!  invp(perm(i))=i
! enddo
!
! allocate(xlnz(neqns+1),xnzsub(neqns+1),rchlnk(neqns),mrglnk(neqns),marker(neqns))
! maxsub=ia(neqns+1)-1
!
! do
!  maxsubinit=maxsub
!  if(allocated(nzsub))deallocate(nzsub)
!  allocate(nzsub(maxsub))
!  call smbfct(neqns,xadj,adjncy,perm,invp,xlnz,maxlnz,xnzsub,nzsub,maxsub,&
!              rchlnk,mrglnk,marker,flag&
!              )
!  if(maxsub.eq.maxsubinit)exit
!  if(flag.eq.0)then
!   exit
!  else
!   write(*,*)' ERROR: smbfct ', flag
!   error stop
!  endif
! enddo
! deallocate(rchlnk,mrglnk,marker)
 
 call symbolicfact(neqns,ia(neqns+1)-1,xadj,adjncy,perm,xlnz,maxlnz,xnzsub,nzsub,maxsub,flag)
 write(*,*)' Symbolic factorization: ',flag,maxsub

 !super following Karin Meyer
 allocate(xspars(ia(neqns+1)),diag(neqns))

 call computexsparsdiag(neqns,ia,ja,a,xlnz,nzsub,xnzsub,xspars,diag,perm)

call expand( neqns, xlnz, xspars, xnzsub, nzsub, diag,hh,.true.)
write(*,*)'matrix init';call print_matrix(neqns,hh)

print*,'xlnz',xlnz
print*,'xnzsub',xnzsub
print*,'ixsub',nzsub
print*,'xspars xspars',xspars
print*,'diag',diag

call expand( neqns, xlnz, xspars, xnzsub, nzsub, diag,hh,.true.)
write(*,*)'matrix init';call print_matrix(neqns,hh)
write(*,*)'matrix init2';call print_matrix(neqns,hh-hh1)

 allocate(inode(neqns))
 call super_nodes(neqns,xlnz,xnzsub,nzsub,nnode,inode)
 write(*,*) 'no. of rows =',neqns,'  no. of super-nodes =',nnode

 do i = 1, nnode
  write(*,'(1x,3(a,i3))') 'node',i,' from row', inode(i+1)+1,'  to row', inode(i)
 end do

 ! Cholesky factorisation
print*,'ixlnz',xlnz
print*,'ixnzsub',xnzsub
print*,'iixsub',nzsub
print*,'idiag',diag
print*,'innode',nnode
print*,'iinode',inode

 call super_gsfct( neqns, xlnz, xspars, xnzsub, nzsub, diag, nnode, inode )

call expand( neqns, xlnz, xspars, xnzsub, nzsub, diag,hh1,.false.)
hh2 = matmul( hh1, transpose(hh1) )
write(*,*)'Cholesky factor'; call print_matrix(neqns, hh1 )
write(*,*)'Matrix';call print_matrix(neqns,hh2)

print*,'diff   ',sum(abs(hh2-hh))
print*,'diff   ',hh2-hh
do i=1,neqns
 write(*,'(1000(f10.4,x))')hh2(i,:)-hh(i,:)
enddo

 ! Matrix inverse
 call super_sparsinv( neqns, xlnz, xspars, xnzsub, nzsub, diag, nnode, inode )
 
call expand( neqns, xlnz, xspars, xnzsub, nzsub, diag,hh3,.true.)
write(*,*)'"Inverse" '; call print_matrix(neqns, hh3 )
 
print*,'xlnz',xlnz
print*,'xspars',xspars
print*,'xnzsub',xnzsub
print*,'nzsub',nzsub
 
 !Convert to ija
 call converttoija( neqns, xlnz, xspars, xnzsub, nzsub, diag,ia,ja,a,perm)

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
 !maxsub=ia(neqns+1)-1
 maxsub=nnzeros

 do
  maxsubinit=maxsub
  if(allocated(nzsub))deallocate(nzsub)
  allocate(nzsub(maxsub))
  call smbfct(neqns,xadj,adjncy,perm,invp,xlnz,maxlnz,xnzsub,nzsub,maxsub,&
              rchlnk,mrglnk,marker,flag&
              )
  if(maxsub.eq.maxsubinit)exit
  if(flag.eq.0)then
   exit
  else
   write(*,*)' ERROR: smbfct ', flag
   error stop
  endif
 enddo
 deallocate(rchlnk,mrglnk,marker)

end subroutine


subroutine super_nodes( neqns, xlnz, xnzsub, ixsub, nnode, inode )
 integer(kind=int32),intent(in)::neqns
 integer(kind=int32),intent(in)::ixsub(:),xlnz(:),xnzsub(:)
 integer(kind=int32),intent(out)::nnode
 integer(kind=int32),intent(out)::inode(:) !size=neqns

 integer(kind=int32)::i,ii,j,n,ilast,kk
 real(kind=wp)::xx

 ! establish boundaries between diaggonal blocks 100% full
 ilast = neqns
 nnode = 0
 do i = neqns-1, 1, -1
    ii = xnzsub(i)
    n = 0
    do j = xlnz(i), xlnz(i+1)-1
       kk = ixsub(ii); ii = ii + 1
       if(  kk > ilast ) exit
       n = n + 1
    end do
    xx = dble( n ) / dble( ilast - i )
    if( xx < 1._wp ) then
     nnode = nnode +1
     inode(nnode) = ilast
     ilast = i
    end if
 end do
 nnode = nnode +1
 inode(nnode) = ilast
 inode(nnode+1) = 0

end subroutine 

subroutine super_gsfct(neqns,xlnz,xspars,xnzsub,ixsub,diag,nnode,inode)
 integer(kind=int32),intent(in)::neqns,nnode
 integer(kind=int32),intent(in)::ixsub(:),xlnz(:),xnzsub(:),inode(:)
 real(kind=wp),intent(inout)::xspars(:),diag(:)

 integer(kind=int32)::i,j,k,jrow,n,ksub,irow,jnode,icol1,icol2,jcol,ii,jj,mm,kk
 integer(kind=int32),allocatable::jvec(:),kvec(:)
 real(kind=wp),allocatable::ttt(:,:),s21(:,:),s22(:,:)

 allocate( jvec(neqns), kvec(neqns),stat = ii )
 if( ii /= 0 ) call alloc_err

 do jnode = nnode, 1, -1
    icol1 = inode(jnode+1) + 1
    icol2 = inode(jnode)
    mm = icol2 - icol1 +1

    !pick out diaggonal block
    allocate( ttt(icol1:icol2,icol1:icol2), stat = ii )
    if( ii /= 0 ) call alloc_err
    ttt = 0.d0; jvec = 0
    do irow = icol1, icol2
       ttt(irow,irow) = diag(irow)
       ksub = xnzsub(irow)
       do i = xlnz(irow), xlnz(irow+1)-1
          jcol = ixsub(ksub); ksub = ksub + 1
          if( jcol <= icol2 ) then
              ttt(jcol,irow) = xspars(i) 
          else
              jvec(jcol) = 1
          end if
       end do
     end do

     !factorise
     call dpotrf( 'L', mm, ttt, mm, ii )
     if( ii /= 0 ) then
         write(*,*)'Dense fact: Routine DPOTRF returned error code', ii
         write(*,*)'... coefficient matrix must be positive definite'; stop
     end if

     !adjust block below diaggonal
     !.. count no. of rows
     n = 0
     do i = icol2+1, neqns
        if( jvec(i) == 0 ) cycle
        n = n + 1
        jvec(i) = n
        kvec(n) = i
     end do
     if( n > 0 ) then
         !... pick out rows
         allocate( s21(n, icol1:icol2), stat = ii )
         if( ii /= 0 ) call alloc_err
         s21 = 0.d0
         do irow = icol1, icol2
            ksub = xnzsub(irow)
            do i = xlnz(irow), xlnz(irow+1)-1
               jcol = ixsub(ksub); ksub = ksub + 1
               if( jcol <= icol2 ) cycle
               jj = jvec(jcol); s21(jj,irow) = xspars(i)
            end do
         end do

         !calculate L21
         call dtrsm( 'R', 'L', 'T', 'N', n, mm, 1.d0, ttt, mm, s21, n )

         !adjust remaining triangle to right: A22 := A22 - L21 L21'
         allocate( s22(n,n), stat = ii ); if( ii /= 0 ) call alloc_err
         call dsyrk( 'L', 'N', n, mm, 1.d0, s21, n, 0.d0, s22, n )
         kk = kvec(n)
         do i = 1, n
             jrow = kvec(i); diag(jrow) = diag(jrow) - s22(i,i)
             ksub = xnzsub(jrow)
             do j = xlnz(jrow), xlnz(jrow+1)-1
                jcol = ixsub(ksub); if( jcol > kk ) exit
                ksub = ksub + 1; ii = jvec(jcol)
                if( ii > 0 ) xspars(j) = xspars(j) - s22(ii,i)
             end do
         end do
         deallocate( s22 )
     end if

     !transfer block back to sparse storage
     do irow = icol1, icol2 
       diag(irow) = ttt(irow,irow)
        ksub = xnzsub(irow)
        do i = xlnz(irow), xlnz(irow+1)-1
           jcol = ixsub(ksub); ksub = ksub + 1
           if( jcol <= icol2 ) then
               xspars(i) = ttt(jcol,irow)
           else
               xspars(i) = s21( jvec(jcol), irow)
           end if
        end do
     end do
     deallocate( ttt )
     if( n > 0 ) deallocate( s21)

 end do ! jnode

 deallocate( jvec, kvec )

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

  allocate( jvec(neqns), kvec(neqns),stat = ii )
  if( ii /= 0 ) call alloc_err

! backwards flops: determine inverse using supernodal blocks
  do jnode = 1, nnode
     icol1 = inode(jnode+1) + 1; icol2 = inode(jnode);  mm = icol2 - icol1 +1
      
!    pick out diagonal block
     allocate( ttt(icol1:icol2,icol1:icol2), stat = ii )
     if( ii /= 0 ) call alloc_err
     ttt = 0.d0; jvec = 0
     do irow = icol1, icol2
        ttt(irow,irow) = diag(irow);  ksub = xnzsub(irow)
        do i = xlnz(irow), xlnz(irow+1)-1
           jcol = ixsub(ksub); ksub = ksub + 1
           if( jcol <= icol2 ) then
               ttt(jcol,irow) = xspars(i) 
          else 
               jvec(jcol) = 1
           end if
        end do
     end do

!    pick out lead columns (condensed)
     n21 = 0
     do i = icol2+1, neqns
        if( jvec(i) == 0 ) cycle
        n21 = n21 + 1; jvec(i) = n21; kvec(n21) = i
     end do
     if( n21 > 0 ) then
         allocate( s21(n21, icol1:icol2), stat = ii )
         if( ii /= 0 ) call alloc_err
         allocate( f21(n21, icol1:icol2), stat = ii )
         if( ii /= 0 ) call alloc_err
         s21 = 0.d0; f21 = 0.d0
         do irow = icol1, icol2
            ksub = xnzsub(irow)
            do i = xlnz(irow), xlnz(irow+1)-1
               jcol = ixsub(ksub); ksub = ksub + 1
               if( jcol <= icol2 ) cycle
               jj = jvec(jcol); s21(jj,irow) = xspars(i)
            end do
         end do
!        ... post-multiply with inverse Chol factor -> solve
         call dtrsm( 'R', 'L', 'N', 'N', n21, mm, 1.d0, ttt, mm, s21, n21 )
!        ... invert Cholesky factor
         call dpotri( 'L',  mm, ttt, mm, ii )
         if( ii /= 0 ) then
             write(*,*) 'Routine DPOTRI returned error code', ii; stop
         end if
!        ... pre-multiply by already inverted submatrix
         iopt = 2
44       if( iopt == 1 ) then
            allocate( rr(icol1:icol2), qx(icol1:icol2), stat = ii )
            if( ii /= 0 ) call alloc_err
            f21 = 0.d0
            do k = 1, n21
               jrow = kvec(k); rr = s21(k,:); qx =  diag(jrow) * rr 
               ksub = xnzsub(jrow)
               do i = xlnz(jrow), xlnz(jrow+1)-1
                  jcol = ixsub(ksub); ksub = ksub + 1
                  m = jvec(jcol); if( m < 1 ) cycle
                  xx = xspars(i); f21(m,:) = f21(m,:) + xx * rr
                  qx = qx + xx * s21(m,:) 
               end do
               f21(k,:) = f21(k,:) + qx 
            end do
            deallocate( rr, qx )
         else
            allocate( s22(n21,n21), stat = ii )
            if( ii /= 0 ) then
                iopt = 1; go to 44
            end if
            s22 = 0.d0
            do k = 1, n21
               jrow = kvec(k); s22(k,k) = diag(jrow)
               ksub = xnzsub(jrow)
               do i = xlnz(jrow), xlnz(jrow+1)-1
                  jcol = ixsub(ksub); ksub = ksub + 1
                  m = jvec(jcol); if( m > 0 ) s22(m,k) = xspars(i)
               end do
            end do
            call dsymm( 'L', 'L', n21, mm, 1.d0, s22, n21, s21, n21, 0.d0,     &
&                                                                f21, n21  )
            deallocate( s22 )
         end if
!        ... adjustments to current block
         call dgemm( 'T', 'N', mm, mm, n21, 1.d0, f21, n21, s21, n21, 1.d0,    &
&                                                                 ttt, mm )
     else
         call dpotri( 'L',  mm, ttt, mm, ii )
         if( ii /= 0 ) then
             write(*,*) 'Routine DPOTRI returned error code', ii; stop
         end if
     end if

!    save current block  
     do irow = icol1, icol2
        diag(irow) = ttt(irow,irow); ksub = xnzsub(irow)
        do i = xlnz(irow), xlnz(irow+1)-1
           jcol = ixsub(ksub); ksub = ksub + 1
           if( jcol <= icol2 ) then
               xspars(i) = ttt(jcol,irow) 
           else 
               xspars(i) = - f21(jvec(jcol),irow)
           end if
        end do
     end do
     deallocate( ttt ); if( n21 > 0 ) deallocate( s21, f21 )
  end do ! jnode

  deallocate( jvec, kvec )

end subroutine 

subroutine computexsparsdiag(n,ia,ja,a,iu,ju,iju,xspars,diag,ip)
 integer:: n,ia(:),ja(:),ip(:),i,j,k,kk
 integer::iu(:),ju(:),iju(:)
 real(kind=wp)::xspars(:)
 real(kind=wp)::a(:),diag(:)

 integer(kind=int32)::irow,iirow,icol
 integer(kind=int32),allocatable::tmp(:)
 real(kind=wp),allocatable::rtmp(:)

 do i=1,n
  irow=ip(i)
  diag(i)=a(ia(irow))
  do k=iu(i),iu(i+1)-1
   iirow=irow
   icol=ip(ju(iju(i)+k-iu(i)))
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

subroutine alloc_err()
 write(*,'(a,i0,3a)')' ERROR (',__LINE__,',',__FILE__,'): failed allocation'
 error stop
end subroutine 

!TO BE DELETED
subroutine print_matrix(neqns, hh )
  real(8), dimension(neqns,neqns), intent(in) :: hh
  integer                                     :: i,neqns

  do i = 1, neqns
     write(*,'(1000(f10.3,x))') hh(i,:)
  end do

end subroutine print_matrix

subroutine expand(neqns,xlnz,xspars,xnzsub,ixsub,diag, hhh, opt )
  real(8), intent(inout), dimension(neqns) :: diag
  real(8),  intent(inout)     :: xspars(:)
  integer,  intent(in)        :: ixsub(:) ,xlnz(:), xnzsub(:)
  integer,  intent(in)        :: neqns
 
  real(8), dimension(neqns,neqns), intent(out) :: hhh
  logical, intent(in)                          :: opt
  integer                                      :: irow, ksub, i, icol

! expand sparse-stored matrix into two-dimensional array
  hhh = 0.d0
  do irow = 1, neqns
     hhh(irow,irow) = diag(irow)
     ksub = xnzsub(irow)
     do i = xlnz(irow), xlnz(irow+1)-1
        icol = ixsub(ksub)
        ksub = ksub + 1
        hhh(icol,irow) = xspars(i); if( opt ) hhh(irow,icol) = xspars(i)
     end do
  end do

end subroutine expand


end module

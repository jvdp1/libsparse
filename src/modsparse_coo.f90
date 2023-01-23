submodule (modsparse) modsparse_coo
 use modhash, only:hashf,roundinguppower2
 !$ use omp_lib
 implicit none

 real(kind=real32),parameter::maxratiofilled_par=0.80

contains

!**CONSTRUCTOR
module function constructor_coo(m,n,nel,lupper,unlog) result(sparse)
 type(coosparse)::sparse
 integer(kind=int32),intent(in)::m
 integer(kind=int32),intent(in),optional::n,unlog
 integer(kind=int64),intent(in),optional::nel
 logical,intent(in),optional::lupper

 call sparse%initialize('COO',m,m)

 if(present(n))sparse%dim2=n
 if(present(lupper))sparse%lupperstorage=lupper
 if(present(unlog))sparse%unlog=unlog

 sparse%lsymmetric=.false.

 sparse%filled=0_int64

 sparse%nel=roundinguppower2(100_int64)
 if(present(nel))sparse%nel=roundinguppower2(int(nel,int64))
 allocate(sparse%ij(2,sparse%nel),sparse%a(sparse%nel))
 sparse%ij=0
 sparse%a=0._wp

end function

module subroutine constructor_sub_coo(sparse,m,n,nel,lupper,unlog)
 class(coosparse),intent(out)::sparse
 integer(kind=int32),intent(in)::m
 integer(kind=int32),intent(in),optional::n,unlog
 integer(kind=int64),intent(in),optional::nel
 logical,intent(in),optional::lupper

 call sparse%initialize('COO',m,m)

 if(present(n))sparse%dim2=n
 if(present(lupper))sparse%lupperstorage=lupper
 if(present(unlog))sparse%unlog=unlog

 sparse%lsymmetric=.false.

 sparse%filled=0_int64

 sparse%nel=roundinguppower2(100_int64)
 if(present(nel))sparse%nel=roundinguppower2(int(nel,int64))
 allocate(sparse%ij(2,sparse%nel))
 sparse%ij=0
 allocate(sparse%a(sparse%nel))
 sparse%a=0._wp

end subroutine

!**DESTROY
module subroutine destroy_coo(sparse)
 class(coosparse),intent(inout)::sparse

 call sparse%destroy_gen_gen()

 sparse%nel=-1_int64
 sparse%filled=-1_int64
 if(allocated(sparse%ij))deallocate(sparse%ij)
 if(allocated(sparse%a))deallocate(sparse%a)

end subroutine

!**DIAGONAL ELEMENTS
module function diag_vect_coo(sparse) result(array)
 class(coosparse),intent(inout)::sparse
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
module recursive subroutine add_coo(sparse,row,col,val)
 class(coosparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::row,col
 real(kind=wp),intent(in)::val

 integer(kind=int64)::hash,i8
 real(kind=real32),parameter::maxratiofilled=maxratiofilled_par
 real(kind=real32)::ratiofilled
 type(coosparse),allocatable::sptmp

 if(.not.validvalue_gen(sparse,row,col))return
 if(.not.validnonzero_gen(sparse,val))return
 if(sparse%lupperstorage.and..not.uppervalue_gen(row,col))return

 hash=hashf(row,col,sparse%ij,sparse%nel,sparse%filled,.false.)
 ratiofilled=real(sparse%filled)/real(sparse%nel)

 if(hash.eq.-1.or.ratiofilled.gt.maxratiofilled)then
  !matrix probably full, or nothing available within the n requested searches
  !1. Copy matrix
  !sptmp=coosparse(sparse%dim1,sparse%dim2,sparse%nel*2)   !to avoid the copy through a temporary array
  allocate(sptmp);call sptmp%init(sparse%dim1,sparse%dim2,int(sparse%nel*1.5_real64, int64))
  do i8=1_int64,sparse%nel
   call sptmp%add(sparse%ij(1,i8),sparse%ij(2,i8),sparse%a(i8))
  enddo
  !2. reallocate matrix using move_alloc
#if(_VERBOSE>0)
  write(sparse%unlog,'(2(a,i0))')'  Current | New size COO: ',sparse%nel,' | ',sptmp%nel
#endif
  sparse%nel=sptmp%nel
  sparse%filled=sptmp%filled
  if(allocated(sparse%ij))deallocate(sparse%ij)
  if(allocated(sparse%a))deallocate(sparse%a)
  call move_alloc(sptmp%ij,sparse%ij)
  call move_alloc(sptmp%a,sparse%a)
  call sptmp%destroy()
  !3. Search for a new address in the new matrix
  hash=hashf(row,col,sparse%ij,sparse%nel,sparse%filled,.false.)
  ratiofilled=real(sparse%filled)/real(sparse%nel)
 endif

 if(hash.gt.0_int64)then!.and.ratiofilled.le.maxratiofilled)then
  sparse%a(hash)=sparse%a(hash)+val
 else
  !is it possible?
  write(sparse%unlog,*)' ERROR: unexpected'!,__FILE__,__LINE__
  stop
 endif

end subroutine

!**GET ELEMENTS
module function get_coo(sparse,row,col) result(val)
 class(coosparse),intent(in)::sparse
 integer(kind=int32),intent(in)::row,col
 real(kind=wp)::val

 integer(kind=int32)::trow,tcol
 integer(kind=int64)::hash

 val=0.0_wp

 trow=row
 tcol=col
 if(sparse%lupperstorage.and.sparse%lsymmetric.and.row.gt.col)then
  !swap row-col
  trow=col
  tcol=row
 endif

 hash=hashf(trow,tcol,sparse%ij,sparse%nel)

 if(hash.gt.0_int64)val=sparse%a(hash)

end function

!** GET MEMORY
module function getmem_coo(sparse) result(getmem)
 class(coosparse),intent(in)::sparse
 integer(kind=int64)::getmem

 getmem=sparse%getmem_gen()+sizeof(sparse%nel)+sizeof(sparse%filled)
 if(allocated(sparse%ij))getmem=getmem+sizeof(sparse%ij)
 if(allocated(sparse%a))getmem=getmem+sizeof(sparse%a)

end function

!**EXTERNAL

!**LOAD
module function load_coo(namefile,unlog) result(sparse)
 type(coosparse)::sparse
 character(len=*),intent(in)::namefile
 integer(kind=int32),intent(in),optional::unlog

 integer(kind=int32)::un,dim1,dim2
 integer(kind=int64)::nonzero,nel
 logical::lupperstorage

 open(newunit=un,file=namefile,action='read',status='old',access='stream')!,buffered='yes')
 read(un)dim1
 if(dim1.ne.typecoo)then
  write(*,'(a)')' ERROR: the proposed file is not a COO file'
  stop
 endif
 read(un)dim1            !int32
 read(un)dim2            !int32
 read(un)nonzero         !int64
 read(un)nel             !int64
 read(un)lupperstorage   !logical

 if(present(unlog))then
  sparse=coosparse(dim1,dim2,nel,lupperstorage,unlog)
 else
  sparse=coosparse(dim1,dim2,nel,lupperstorage)
 endif

 sparse%filled=nonzero
 read(un)sparse%ij              !int32
 read(un)sparse%a               !wp
 close(un)

end function

!**MULTIPLICATIONS
module subroutine multgenv_coo(sparse,alpha,trans,x,val,y)
 !Computes y=val*y+alpha*sparse(tranposition)*x
 class(coosparse),intent(in)::sparse
 real(kind=wp),intent(in)::val,alpha
 real(kind=wp),intent(in)::x(:)
 real(kind=wp),intent(out)::y(:)
 character(len=1),intent(in)::trans

 integer(kind=int64)::i
 integer(kind=int32)::j,k
 character(len=1)::matdescra(6)

 if(trans.eq.'N'.or.trans.eq.'n')then
  if(size(y).ne.sparse%getdim(1).or.size(x).ne.sparse%getdim(2))then
   write(sparse%unlog,'(a)')'  ERROR (mult): wrong dimensions'
   error stop
  endif
 elseif(trans.eq.'T'.or.trans.eq.'t')then
  if(size(y).ne.sparse%getdim(2).or.size(x).ne.sparse%getdim(1))then
   write(sparse%unlog,'(a)')'  ERROR (mult): wrong dimensions'
   error stop
  endif
 else
  write(sparse%unlog,'(a)')'  ERROR (mult): wrong transposition'
  error stop
 endif

 matdescra=''

 if(sparse%lsymmetric.and.sparse%lupperstorage)then
  matdescra(1)='S'
 elseif(.not.sparse%lsymmetric.and.sparse%lupperstorage)then
  matdescra(1)='T'
 elseif(.not.sparse%lsymmetric.and..not.sparse%lupperstorage)then
  matdescra(1)='G'
 else
  write(sparse%unlog,'(a)')'  ERROR (mult): unsupported format'
  call sparse%printstats
  error stop
 endif
!
! if(sparse%lupperstorage)then
!  matdescra(2)='U'
!  matdescra(3)='N'
! endif
!
! matdescra(4)='F'

 !don't forget transposition

 y = val * y
 select case(matdescra(1))
  case('S')
   do i = 1, sparse%nel
    if(sparse%ij(1,i).eq.0)cycle
    j=sparse%ij(1,i)
    k=sparse%ij(2,i)
    y(j) = y(j) + alpha * sparse%a(i) * x(k)
    if(j.ne.k) y(k) = y(k) + alpha * sparse%a(i) * x(j)
   enddo
  case('T','G')
   if(trans.eq.'N'.or.trans.eq.'n')then
    do i = 1, sparse%nel
     if(sparse%ij(1,i).eq.0)cycle
     j=sparse%ij(1,i)
     k=sparse%ij(2,i)
     y(j) = y(j) + alpha * sparse%a(i) * x(k)
    enddo
   elseif(trans.eq.'T'.or.trans.eq.'t')then
    do i = 1, sparse%nel
     if(sparse%ij(1,i).eq.0)cycle
     j=sparse%ij(2,i)
     k=sparse%ij(1,i)
     y(j) = y(j) + alpha * sparse%a(i) * x(k)
    enddo
   endif
  case default
   write(sparse%unlog,'(a)')'  ERROR (multbyv): unsupported format'
   error stop
 end select

end subroutine

module subroutine multgenm_coo(sparse,alpha,trans,x,val,y)
 !Computes y=val*y+alpha*sparse(tranposition)*x
 class(coosparse),intent(in)::sparse
 real(kind=wp),intent(in)::val,alpha
 real(kind=wp),intent(in)::x(:,:)
 real(kind=wp),intent(out)::y(:,:)
 character(len=1),intent(in)::trans

 integer(kind=int64)::i
 integer(kind=int32)::j,k
 character(len=1)::matdescra(6)

 if(trans.eq.'N'.or.trans.eq.'n')then
  if(size(y,1).ne.sparse%getdim(1).or.size(x,1).ne.sparse%getdim(2))then
   write(sparse%unlog,'(a)')'  ERROR (mult): wrong dimensions'
   error stop
  endif
 elseif(trans.eq.'T'.or.trans.eq.'t')then
  if(size(y,1).ne.sparse%getdim(2).or.size(x,1).ne.sparse%getdim(1))then
   write(sparse%unlog,'(a)')'  ERROR (mult): wrong dimensions'
   error stop
  endif
 else
  write(sparse%unlog,'(a)')'  ERROR (mult): wrong transposition'
  error stop
 endif

 matdescra=''

 if(sparse%lsymmetric.and.sparse%lupperstorage)then
  matdescra(1)='S'
 elseif(.not.sparse%lsymmetric.and.sparse%lupperstorage)then
  matdescra(1)='T'
 elseif(.not.sparse%lsymmetric.and..not.sparse%lupperstorage)then
  matdescra(1)='G'
 else
  write(sparse%unlog,'(a)')'  ERROR (mult): unsupported format'
  call sparse%printstats
  error stop
 endif
!
! if(sparse%lupperstorage)then
!  matdescra(2)='U'
!  matdescra(3)='N'
! endif
!
! matdescra(4)='F'

 y = val * y
 select case(matdescra(1))
  case('S')
   do i = 1, sparse%nel
    if(sparse%ij(1,i).eq.0)cycle
    j=sparse%ij(1,i)
    k=sparse%ij(2,i)
    y(j,:) = y(j,:) + alpha * sparse%a(i) * x(k,:)
    if(j.ne.k) y(k,:) = y(k,:) + alpha * sparse%a(i) * x(j,:)
   enddo
  case('T','G')
   if(trans.eq.'N'.or.trans.eq.'n')then
    do i = 1, sparse%nel
     if(sparse%ij(1,i).eq.0)cycle
     j=sparse%ij(1,i)
     k=sparse%ij(2,i)
     y(j,:) = y(j,:) + alpha * sparse%a(i) * x(k,:)
    enddo
   elseif(trans.eq.'T'.or.trans.eq.'t')then
    do i = 1, sparse%nel
     if(sparse%ij(1,i).eq.0)cycle
     j=sparse%ij(2,i)
     k=sparse%ij(1,i)
     y(j,:) = y(j,:) + alpha * sparse%a(i) * x(k,:)
    enddo
   endif
  case default
   write(sparse%unlog,'(a)')'  ERROR (multbyv): unsupported format'
   error stop
 end select

end subroutine

!**NUMBER OF ELEMENTS
module function totalnumberofelements_coo(sparse) result(nel)
 class(coosparse),intent(in)::sparse
 integer(kind=int64)::nel

 nel=sparse%filled

end function

!**PRINT
module subroutine print_coo(sparse,lint,output)
 class(coosparse),intent(in)::sparse
 integer(kind=int32),intent(in),optional::output
 logical,intent(in),optional::lint

 integer(kind=int32)::un,row,col
 integer(kind=int64)::i8
 real(kind=wp)::val
 character(len=30)::frm='(2(i0,1x),g0)'
 logical::linternal

 linternal=.true.
 if(present(lint))linternal=lint

 un=sparse%unlog
 if(present(output))un=output

 do i8=1,sparse%nel
  row=sparse%ij(1,i8)
  col=sparse%ij(2,i8)
  if(row.eq.0.and.col.eq.0)cycle
  val=sparse%a(i8)
  !if(.not.validvalue_gen(sparse,row,col))cycle  !it should never happen
  !if(.not.validnonzero_gen(sparse,val))cycle    !to print as internal
  write(un,frm)row,col,val
  if(.not.linternal.and.sparse%lupperstorage.and.sparse%lsymmetric.and.row.ne.col)then
   write(un,frm)col,row,val
  endif
 enddo

end subroutine

module subroutine printsquare_coo(sparse,output)
 class(coosparse),intent(inout)::sparse
 integer(kind=int32),intent(in),optional::output

 integer(kind=int32)::i,j,un
 real(kind=wp),allocatable::tmp(:)

 un=sparse%unlog
 if(present(output))un=output

 allocate(tmp(sparse%dim2))

 do i=1,sparse%dim1
  tmp=0._wp
  do j=1,sparse%dim2
   tmp(j)=sparse%get(i,j)
  enddo
  write(un,'(*(f9.3,1x))')tmp
 enddo

 deallocate(tmp)

end subroutine

!**SAVE
module subroutine save_coo(sparse,namefile)
 class(coosparse),intent(in)::sparse
 character(len=*),intent(in)::namefile

 integer(kind=int32)::un

 open(newunit=un,file=namefile,action='write',status='replace',access='stream')!,buffered='yes')
 write(un)typecoo                !int32
 write(un)sparse%dim1            !int32
 write(un)sparse%dim2            !int32
 write(un)sparse%nonzero()       !int64
 write(un)sparse%nel             !int64
 write(un)sparse%lupperstorage   !logical
 write(un)sparse%ij              !int32
 write(un)sparse%a               !wp
 close(un)

end subroutine

!**SCALE ALL ENTRIES
module subroutine scale_coo(sparse,val)
 class(coosparse),intent(inout)::sparse
 real(kind=wp),intent(in)::val
 sparse%a = sparse%a * val
end subroutine

!**SET ELEMENTS
module recursive subroutine set_coo(sparse,row,col,val)
 !from add_coo
 class(coosparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::row,col
 real(kind=wp),intent(in)::val

 integer(kind=int64)::hash,i8
 real(kind=real32),parameter::maxratiofilled=maxratiofilled_par
 real(kind=real32)::ratiofilled
 type(coosparse)::sptmp

 if(.not.validvalue_gen(sparse,row,col))return
 !if(.not.validnonzero_gen(sparse,val))return
 if(sparse%lupperstorage.and..not.uppervalue_gen(row,col))return

 hash=hashf(row,col,sparse%ij,sparse%nel,sparse%filled,.false.)
 ratiofilled=real(sparse%filled)/real(sparse%nel)

 if(hash.eq.-1.or.ratiofilled.gt.maxratiofilled)then
  !matrix probably full, or nothing available within the n requested searches
  !1. Copy matrix
  sptmp=coosparse(sparse%dim1,sparse%dim2,int(sparse%nel*1.5_real64, int64))
  do i8=1_int64,sparse%nel
   call sptmp%add(sparse%ij(1,i8),sparse%ij(2,i8),sparse%a(i8))
  enddo
  !2. reallocate matrix using move_alloc
  write(sparse%unlog,'(2(a,i0))')'  Current | New size COO: ',sparse%nel,' | ',sptmp%nel
  sparse%nel=sptmp%nel
  sparse%filled=sptmp%filled
  if(allocated(sparse%ij))deallocate(sparse%ij)
  if(allocated(sparse%a))deallocate(sparse%a)
  call move_alloc(sptmp%ij,sparse%ij)
  call move_alloc(sptmp%a,sparse%a)
  call sptmp%destroy()
  !3. Search for a new address in the new matrix
  hash=hashf(row,col,sparse%ij,sparse%nel,sparse%filled,.false.)
  ratiofilled=real(sparse%filled)/real(sparse%nel)
 endif

 if(hash.gt.0_int64)then!.and.ratiofilled.le.maxratiofilled)then
  sparse%a(hash)=val
 else
  !is it possible?
  write(sparse%unlog,*)' ERROR: unexpected'!,__FILE__,__LINE__
  stop
 endif

end subroutine

!**SOLVE

!**SORT ARRAY

!**SUBMATRIX
module function submatrix_coo(sparse,startdim1,enddim1,startdim2,enddim2,lupper,unlog) result(subsparse)
 !Not programmed efficiently, but it should do the job
 class(coosparse),intent(in)::sparse
 type(coosparse)::subsparse
 integer(kind=int32),intent(in)::startdim1,enddim1,startdim2,enddim2
 integer(kind=int32),intent(in),optional::unlog
 logical,intent(in),optional::lupper

 integer(kind=int32)::i,j
 integer(kind=int64)::i8,nel
 logical::lincludediag,lupperstorage


 if(.not.validvalue_gen(sparse,startdim1,startdim2))return
 if(.not.validvalue_gen(sparse,enddim1,enddim2))return

 nel=10000

 !check if the submatrix include diagonal elements of sparse
 ! if yes -> lupperstorage
 ! if no  -> .not.lupperstorage
 lincludediag=.false.
 firstdim: do i=startdim1,enddim1
  do j=startdim2,enddim2
   if(i.eq.j)then
    lincludediag=.true.
    exit firstdim
   endif
  enddo
 enddo firstdim

 lupperstorage=sparse%lupperstorage
 if(present(lupper))then
  lupperstorage=lupper
 else
  if(sparse%lupperstorage)then
   lupperstorage=.false.
   if(lincludediag)lupperstorage=.true.
  endif
 endif


 if(present(unlog))then
  subsparse=coosparse(enddim1-startdim1+1,enddim2-startdim2+1,int(nel,int64),lupperstorage,unlog)
 else
  subsparse=coosparse(enddim1-startdim1+1,enddim2-startdim2+1,int(nel,int64),lupperstorage)
 endif

 if(sparse%lupperstorage.eqv.lupperstorage.or.(sparse%lupperstorage.and..not.lincludediag))then
  ! upper -> upper  ||  full -> full
  do i8=1,sparse%nel
   i=sparse%ij(1,i8)
   if(i.eq.0)cycle
   j=sparse%ij(2,i8)
   if((i.ge.startdim1.and.i.le.enddim1).and.(j.ge.startdim2.and.j.le.enddim2))then
    call subsparse%add(i-startdim1+1,j-startdim2+1,sparse%a(i8))
   endif
  enddo
 elseif(sparse%lupperstorage.and..not.lupperstorage)then
  ! upper -> full
  do i8=1,sparse%nel
   i=sparse%ij(1,i8)
   if(i.eq.0)cycle
   j=sparse%ij(2,i8)
   if((i.ge.startdim1.and.i.le.enddim1).and.(j.ge.startdim2.and.j.le.enddim2))then
    call subsparse%add(i-startdim1+1,j-startdim2+1,sparse%a(i8))
   endif
!   if(i.ne.j)then
!    if((j.ge.startdim1.and.j.le.enddim1).and.(i.ge.startdim2.and.i.le.enddim2))then
!     call subsparse%add(j-startdim1+1,i-startdim2+1,sparse%a(i8))
!    endif
!   endif
  enddo
 elseif(.not.sparse%lupperstorage.and.lupperstorage)then
  ! full -> upper
  do i8=1,sparse%nel
   i=sparse%ij(1,i8)
   if(i.eq.0)cycle
   j=sparse%ij(2,i8)
   if((j-startdim2+1.ge.i-startdim1+1).and.(i.ge.startdim1.and.i.le.enddim1).and.(j.ge.startdim2.and.j.le.enddim2))then
    call subsparse%add(i-startdim1+1,j-startdim2+1,sparse%a(i8))
   endif
  enddo
 endif

end function

module function submatrix_index_coo(sparse,indvector,sizeblock,unlog) result(subsparse)
 !Not programmed efficiently, but it should do the job
 class(coosparse),intent(in)::sparse
 type(coosparse)::subsparse
 integer(kind=int32),intent(in)::indvector(:)
 integer(kind=int32),intent(in),optional::sizeblock
 integer(kind=int32),intent(in),optional::unlog

 integer(kind=int32)::i,j
 integer(kind=int32)::i_
 integer(kind=int32)::ii,jj
 integer(kind=int64)::sizeblock_
 integer(kind=int32)::startdim1,enddim1,startdim2,enddim2
 integer(kind=int64)::nel


 if(.not.validvalue_gen(sparse,minval(indvector),minval(indvector)))return
 if(.not.validvalue_gen(sparse,maxval(indvector),maxval(indvector)))return
 if(.not.sparse%lupperstorage)return

 nel = size(indvector)

 sizeblock_ = nel
 if(present(sizeblock))sizeblock_=sizeblock


 startdim1 = minval(indvector)
 enddim1 = maxval(indvector)
 startdim2 = minval(indvector)
 enddim2 = maxval(indvector)

 if(present(unlog))then
  call subsparse%init(size(indvector), size(indvector), nel, sparse%lupperstorage, unlog)
 else
  call subsparse%init(size(indvector), size(indvector), nel, sparse%lupperstorage)
 endif

 ! upper -> upper
 do i_ = 1, nel, sizeblock_
  do i = i_, i_ + sizeblock_ - 1
   do j = i, i_ + sizeblock_ - 1
    ii = indvector(i)
    jj = indvector(j)
    if(ii.le.jj)then
     call subsparse%add(i, j, sparse%get(ii,jj))
    else
     call subsparse%add(i, j, sparse%get(jj,ii))
    endif
   enddo
  enddo
 enddo

end function

end submodule

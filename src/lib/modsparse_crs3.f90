submodule (modsparse) modsparse_crs3
#if (_DP==0)
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real32
#else
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real64
#endif
  !$ use omp_lib
 implicit none


contains

!**CONSTRUCTOR
module function constructor_crs3(m,nel,n,lupper,unlog) result(sparse)
 type(crs3sparse)::sparse
 integer(kind=int32),intent(in)::m
 integer(kind=int32),intent(in)::nel
 integer(kind=int32),intent(in),optional::n,unlog
 logical,intent(in),optional::lupper

 call sparse%initialize('CRS3',m,m)

 if(present(n))sparse%dim2=n
 if(present(lupper))sparse%lupperstorage=lupper
 if(present(unlog))sparse%unlog=unlog

 sparse%lsymmetric=.false.

 allocate(sparse%ia(sparse%dim1+1),sparse%ja(nel),sparse%a(nel))
 sparse%ia=0
 sparse%ia(sparse%dim1+1)=-nel
 sparse%ja=0
 sparse%a=0._wp

end function

module subroutine constructor_sub_crs3(sparse,m,nel,n,lupper,unlog)
 class(crs3sparse),intent(out)::sparse
 integer(kind=int32),intent(in)::m
 integer(kind=int32),intent(in)::nel
 integer(kind=int32),intent(in),optional::n,unlog
 logical,intent(in),optional::lupper

 call sparse%initialize('CRS3',m,m)

 if(present(n))sparse%dim2=n
 if(present(lupper))sparse%lupperstorage=lupper
 if(present(unlog))sparse%unlog=unlog

 sparse%lsymmetric=.false.

 allocate(sparse%ia(sparse%dim1+1),sparse%ja(nel),sparse%a(nel))
 sparse%ia=0
 sparse%ia(sparse%dim1+1)=-nel
 sparse%ja=0
 sparse%a=0._wp

end subroutine

!**ADD ELEMENTS
module subroutine add_crs3(sparse,row,col,val,error)
 !add a value only to an existing one
 class(crs3sparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::row,col
 integer(kind=int32),intent(out),optional::error
 real(kind=wp),intent(in)::val

 integer(kind=int32)::i
 integer(kind=int32)::ierror    !added: error=0;Not existing: error=-1;matrix not initialized: error=-10

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
module function get_crs3(sparse,row,col) result(val)
 class(crs3sparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::row,col
 real(kind=wp)::val

 integer(kind=int32)::i,trow,tcol

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
module function getmem_crs3(sparse) result(getmem)
 class(crs3sparse),intent(in)::sparse
 integer(kind=int64)::getmem

 getmem=sparse%getmem_gen()
 if(allocated(sparse%ia))getmem=getmem+sizeof(sparse%ia)
 if(allocated(sparse%ja))getmem=getmem+sizeof(sparse%ja)
 if(allocated(sparse%a))getmem=getmem+sizeof(sparse%a)
#if (_PARDISO==1)
 select type(sparse)
  type is(crssparse)
   getmem=getmem+sizeof(sparse%lpardisofirst)
 end select
#endif

end function


!**EXTERNAL
module subroutine external_crs3(sparse,ia,ja,a)
 class(crs3sparse),intent(inout)::sparse
 integer(kind=int32),intent(in)::ia(:),ja(:)
 real(kind=wp),intent(in)::a(:)

 if(size(ia).ne.size(sparse%ia))then
  write(sparse%unlog,'(a)')' ERROR: The provided array ia is of a different size!'
  stop
 endif
 if(size(ja).ne.size(sparse%ja))then
  write(sparse%unlog,'(a)')' ERROR: The provided array ja is of a different size!'
  stop
 endif
 if(size(a).ne.size(sparse%a))then
  write(sparse%unlog,'(a)')' ERROR: The provided array a is of a different size!'
  stop
 endif

 sparse%ia=ia
 sparse%ja=ja
 sparse%a=a

end subroutine


!**MULTIPLICATIONS
module subroutine multgenv_csr3(sparse,alpha,trans,x,val,y)
 !Computes y=val*y+alpha*sparse(tranposition)*x
 class(crs3sparse),intent(in)::sparse
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
 call mkl_scsrmv(trans,sparse%dim1,sparse%dim2,alpha,matdescra,sparse%a,sparse%ja,sparse%ia(1:sparse%dim1),sparse%ia(2:sparse%dim1+1),x,val,y)
#else
 call mkl_dcsrmv(trans,sparse%dim1,sparse%dim2,alpha,matdescra,sparse%a,sparse%ja,sparse%ia(1:sparse%dim1),sparse%ia(2:sparse%dim1+1),x,val,y)
#endif

end subroutine

module subroutine multgenm_csr3(sparse,alpha,trans,x,val,y)
 !Computes y=val*y+alpha*sparse(tranposition)*x
 class(crs3sparse),intent(in)::sparse
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
module function totalnumberofelements_crs3(sparse) result(nel)
 class(crs3sparse),intent(in)::sparse
 integer(kind=int64)::nel

 nel=int(sparse%ia(sparse%dim1+1),int64)-1_int64

end function

!**PRINT
module subroutine print_crs3(sparse,lint,output)
 class(crs3sparse),intent(in)::sparse
 integer(kind=int32),intent(in),optional::output
 logical,intent(in),optional::lint

 integer(kind=int32)::i
 integer(kind=int32)::un,j
 character(len=30)::frm='(2(i0,1x),g0)'
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

module subroutine printsquare_crs3(sparse,output)
 class(crs3sparse),intent(inout)::sparse
 integer(kind=int32),intent(in),optional::output

 integer(kind=int32)::i,j,un
 real(kind=wp)::val
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
module subroutine save_crs3(sparse,namefile)
 class(crs3sparse),intent(in)::sparse
 character(len=*),intent(in)::namefile

 integer(kind=int32)::un

 open(newunit=un,file=namefile,action='write',status='replace',access='stream')!,buffered='yes')

 select type(sparse)
  type is(crs3sparse)
   write(un)typecrs3             !int32
  type is(crssparse)
   write(un)typecrs              !int32
 end select

 write(un)sparse%dim1            !int32
 write(un)sparse%dim2            !int32
 write(un)sparse%nonzero()       !int64
 write(un)sparse%lupperstorage   !logical
 write(un)sparse%ia              !int32
 write(un)sparse%ja              !int32
 write(un)sparse%a               !wp
 close(un)

end subroutine


end submodule

submodule (modsparse) modsparse_crs3
#if (_DP==0)
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real32
#else
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real64
#endif
  !$ use omp_lib
 implicit none


contains

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

end submodule

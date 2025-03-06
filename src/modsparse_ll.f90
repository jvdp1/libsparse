submodule (modsparse) modsparse_ll
 !$ use omp_lib
 implicit none

contains

!**CONSTRUCTOR
module subroutine constructor_sub_ll(sparse,m,n,lupper,unlog)
 class(llsparse),intent(out)::sparse
 integer(kind=int32),intent(in)::m
 integer(kind=int32),intent(in),optional::n,unlog
 logical,intent(in),optional::lupper

 call sparse%initialize('LINKED LIST',m,m)

 if(present(n))sparse%dim2=n
 if(present(lupper))sparse%lupperstorage=lupper
 if(present(unlog))sparse%unlog=unlog

 sparse%lsymmetric=.false.

 allocate(sparse%heads(sparse%dim1))

end subroutine

!**DESTROY
module elemental subroutine destroy_scal_ptrnode(pnode)
 type(ptrnode), intent(inout)::pnode

 type(ptrnode)::cursor

 do while(associated(pnode%p))
  cursor=pnode
  pnode=pnode%p%next
  deallocate(cursor%p)
  nullify(cursor%p)
 enddo

end subroutine

module elemental subroutine destroy_ll(sparse)
 class(llsparse),intent(inout)::sparse
 integer(kind=int32)::i

 call sparse%destroy_gen_gen()

 if(allocated(sparse%heads))then
  do i=1,size(sparse%heads)
   call destroy_scal_ptrnode(sparse%heads(i))
  enddo
  deallocate(sparse%heads)
 endif

end subroutine

!**DIAGONAL ELEMENTS

!EQUALITIES
module subroutine equal_node(nodeout,nodein)
 class(node),intent(out)::nodeout
 class(node),intent(in)::nodein

 nodeout%col=nodein%col
 nodeout%val=nodein%val

end subroutine

!**ADD ELEMENTS
module subroutine addtohead_ptrnode(pnode,col,val)
 type(ptrnode),intent(inout),pointer::pnode
 integer(kind=int32),intent(in)::col
 real(kind=wp),intent(in)::val

 type(ptrnode)::cursor

 allocate(cursor%p)
 cursor%p%next=pnode
 cursor%p%col=col
 cursor%p%val=val
 pnode=cursor

end subroutine

module subroutine addtohead_ll(sparse,row,col,val)
 class(llsparse),intent(inout),target::sparse
 integer(kind=int32),intent(in)::row,col
 real(kind=wp),intent(in)::val

 type(ptrnode)::cursor

 if(.not.validvalue_gen(sparse,row,col))return
 if(.not.validnonzero_gen(sparse,val))return
 if(sparse%lupperstorage.and..not.uppervalue_gen(row,col))return

 allocate(cursor%p)
 cursor%p%next=sparse%heads(row)
 cursor%p%col=col
 cursor%p%val=val
 sparse%heads(row)=cursor

end subroutine

module subroutine addinorder_ll(sparse,row,col,val)
 class(llsparse),intent(inout),target::sparse
 integer(kind=int32),intent(in)::row,col
 real(kind=wp),intent(in)::val

 type(ptrnode),pointer::cursor

 if(.not.validvalue_gen(sparse,row,col))return
 if(.not.validnonzero_gen(sparse,val))return
 if(sparse%lupperstorage.and..not.uppervalue_gen(row,col))return

 cursor=>sparse%heads(row)
 do while(associated(cursor%p))
  if(cursor%p%col.ge.col)then
   if(.not.col.ge.cursor%p%col)then
    call addtohead_ptrnode(cursor,col,val)
   else
    cursor%p%val=cursor%p%val+val
   endif
   return
  endif
  cursor=>cursor%p%next
 enddo
 allocate(cursor%p)
 cursor%p%col=col
 cursor%p%val=val

end subroutine

module subroutine addtotail_ll(sparse,row,col,val)
 class(llsparse),intent(inout),target::sparse
 integer(kind=int32),intent(in)::row,col
 real(kind=wp),intent(in)::val

 type(ptrnode),pointer::cursor

 if(.not.validvalue_gen(sparse,row,col))return
 if(.not.validnonzero_gen(sparse,val))return
 if(sparse%lupperstorage.and..not.uppervalue_gen(row,col))return

 cursor=>sparse%heads(row)
 do while(associated(cursor%p))
  cursor=>cursor%p%next
 enddo
 allocate(cursor%p)
 cursor%p%col=col
 cursor%p%val=val

end subroutine

!**GET ELEMENTS
pure module function get_ll(sparse,row,col) result(val)
 class(llsparse),intent(in)::sparse
 integer(kind=int32),intent(in)::row,col
 real(kind=wp)::val

 integer(kind=int32)::trow,tcol
 type(ptrnode),pointer::cursor
 type(ptrnode),target, allocatable::replacecursor

 val=0.0_wp

 trow=row
 tcol=col
 if(sparse%lupperstorage.and.row.gt.col)then
  !swap row-col
  trow=col
  tcol=row
 endif

 !cursor=>sparse%heads(trow)
 allocate(replacecursor, source=sparse%heads(trow))
 cursor=>replacecursor

 do while(associated(cursor%p))
  if(cursor%p%col.eq.tcol)then
   val=cursor%p%val
   exit
  endif
  cursor=>cursor%p%next
 enddo

end function

!**EXTERNAL

!**LOAD

!**MULTIPLICATIONS
module subroutine multgenv_ll(sparse,alpha,trans,x,val,y)
 !Computes y=val*y+alpha*sparse(tranposition)*x
 class(llsparse),intent(in)::sparse
 real(kind=wp),intent(in)::val,alpha
 real(kind=wp),intent(in)::x(:)
 real(kind=wp),intent(out)::y(:)
 character(len=1),intent(in)::trans

 y=0._wp
 write(sparse%unlog,'(a)')' ERROR: Multiplication mult not implemented for llsparse'
 error stop

end subroutine

module subroutine multgenm_ll(sparse,alpha,trans,x,val,y)
 !Computes y=val*y+alpha*sparse(tranposition)*x
 class(llsparse),intent(in)::sparse
 real(kind=wp),intent(in)::val,alpha
 real(kind=wp),intent(in)::x(:,:)
 real(kind=wp),intent(out)::y(:,:)
 character(len=1),intent(in)::trans

 y=0._wp
 write(sparse%unlog,'(a)')' ERROR: Multiplication mult not implemented for llsparse'
 error stop

end subroutine

!**NUMBER OF ELEMENTS
module function totalnumberofelements_ptrnode(pnode) result(nel)
 class(ptrnode),intent(in),target::pnode
 integer(kind=int64)::nel

 type(ptrnode),pointer::cursor

 nel=0
 cursor=>pnode
 do while(associated(cursor%p))
  nel=nel+1
  cursor=>cursor%p%next
 enddo

end function

module function totalnumberofelements_ll(sparse) result(nel)
 class(llsparse),intent(in)::sparse
 integer(kind=int64)::nel

 integer(kind=int32)::i

 nel=0
 do i=1,sparse%dim1
  nel=nel+sparse%heads(i)%size()
 enddo

end function

!**SAVE

!**SCALE ALL ENTRIES
module subroutine scale_ll(sparse,val)
 class(llsparse),intent(inout)::sparse
 real(kind=wp),intent(in)::val

 write(sparse%unlog,'(a)')' ERROR: scale not implemented for llsparse'
 error stop

end subroutine

!**SET ELEMENTS

!**SOLVE

!**SORT ARRAY

!**SUBMATRIX

!**PRINT
module subroutine print_ll(sparse,lint,output)
 class(llsparse),intent(in)::sparse
 integer(kind=int32),intent(in),optional::output
 logical,intent(in),optional::lint

 integer(kind=int32)::i,un
 character(len=20)::frm='(2(i0,1x),g0)'
 logical::linternal
 type(ptrnode),pointer::cursor
 type(ptrnode),target::replacecursor

 linternal=.true.
 if(present(lint))linternal=lint

 un=sparse%unlog
 if(present(output))un=output

 do i=1,sparse%dim1
  !cursor=>sparse%heads(i)
  replacecursor=sparse%heads(i)
  cursor=>replacecursor
  do while(associated(cursor%p))
   write(un,frm)i,cursor%p%col,cursor%p%val
   if(.not.linternal.and.sparse%lupperstorage.and.sparse%lsymmetric.and.cursor%p%col.ne.i)then
    write(un,frm)cursor%p%col,i,cursor%p%val
   endif
   cursor=>cursor%p%next
  enddo
 enddo

end subroutine

module subroutine print_idx_ll(sparse,lidx,lint,output)
 class(llsparse),intent(in)::sparse
 logical,intent(in)::lidx(:)
 integer(kind=int32),intent(in),optional::output
 logical,intent(in),optional::lint

 integer(kind=int32)::i,un
 character(len=20)::frm='(2(i0,1x),g0)'
 logical::linternal
 type(ptrnode),pointer::cursor
 type(ptrnode),target::replacecursor

 linternal=.true.
 if(present(lint))linternal=lint

 un=sparse%unlog
 if(present(output))un=output

 do i=1,sparse%dim1
  if(.not.lidx(i))cycle
  !cursor=>sparse%heads(i)
  replacecursor=sparse%heads(i)
  cursor=>replacecursor
  do while(associated(cursor%p))
   if(.not.lidx(cursor%p%col))then
    write(un,frm)i,cursor%p%col,cursor%p%val
    if(.not.linternal.and.sparse%lupperstorage.and.sparse%lsymmetric.and.cursor%p%col.ne.i)then
     write(un,frm)cursor%p%col,i,cursor%p%val
    endif
   endif
   cursor=>cursor%p%next
  enddo
 enddo

end subroutine

module subroutine printsquare_ll(sparse,output)
 class(llsparse),intent(inout)::sparse
 integer(kind=int32),intent(in),optional::output

 integer(kind=int32)::un
 real(kind=wp),allocatable::tmp(:)

 un=sparse%unlog
 if(present(output))un=output

 write(un,'(a)')' Warning: this subroutine is not implemented yet!'

 allocate(tmp(sparse%dim2))

end subroutine


end submodule

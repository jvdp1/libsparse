module modhash
 use modkind
 implicit none
 public::hashf

contains
!Inspired by lookup3.c from Bob Jenkins (http://burtleburtle.net/bob/hash/index.html#lookup)
!Converted in Fortran by Francois Guillaume - 2011
!Simplifications made by Jeremie Vandenplas - 2018

!PUBLIC
function hashf(row,col,mat,dim2,filled,getval) result(hashvr)
 !hashvr: address (column) of mat
 !mat of size dim1 (=2) x dim2
 !filled: number of elements
 !getval .eq. .true. : search for a value and returns 0 if absent
 !getval .eq. .false.: add a value if row,col was not present before
 integer(kind=int4)::hashvr
 integer(kind=int4),intent(in)::row,col!,mode
 integer(kind=int4),intent(inout)::mat(:,:)
 integer(kind=int8),intent(in)::dim2
 integer(kind=int8),intent(inout)::filled
 logical,intent(in)::getval

 integer(kind=int4),parameter::dim1=2
 integer(kind=int4),parameter::maxiter=5000

 integer(kind=int8)::a,b,c,iaddress,plage
 integer(kind=int4)::i,j,k 
 logical::indnotzero,indnotequal
 !logical::indzero,indequal

 plage=dim2                !size of the array
 a=int(row,kind(a))        !conversion of 1st coordinate
 b=int(col,kind(b))        !conversion of 2nd coordinate
 c=305419896_int8          !default value for 3rd coordinate  

 !Hashing
 call mix(a,b,c)
  
 !Computation of the address
 iaddress=iand(c,plage-1)+1
  
 !Cycle until a free entry is found
 do k=1,maxiter
!  indnotzero=.false.;indnotequal=.false.
!  do i=1,dim1
!   j=mat(i,iaddress)
!   if(j.ne.rowcol(i))indnotequal=.true.
!   if(j.ne.0)indnotzero=.true.
!  enddo     
  indnotzero=.false.;indnotequal=.false.
  if(mat(1,iaddress).ne.row.or.mat(2,iaddress).ne.col)indnotequal=.true.
  if(mat(1,iaddress).ne.0.or.mat(2,iaddress).ne.0)indnotzero=.true.
  if(.not.indnotzero.or..not.indnotequal)then
   if(.not.getval.and..not.indnotzero)then
    mat(1,iaddress)=row
    mat(2,iaddress)=col
    filled=filled+1
   endif
   if(getval.and..not.indnotzero)then
    hashvr=0
   else
    hashvr=iaddress
   endif
   return
  endif

!  indzero=.false.;indequal=.false.
!  if(mat(1,iaddress).eq.row.and.mat(2,iaddress).eq.col)indequal=.true.
!  if(mat(1,iaddress).eq.0.or.mat(2,iaddress).eq.0)indzero=.true.
!  if(indzero.or.indequal)then
!   if(.not.getval.and.indzero)then
!    mat(1,iaddress)=row
!    mat(2,iaddress)=col
!    filled=filled+1
!   endif
!   if(getval.and.indzero)then
!    hashvr=0
!   else
!    hasvr=iaddress
!   endif
!   return
!  endif

  !Hash again
  call mix(a,b,c)
  iaddress=iand(c,plage-1)+1
 enddo

 hashvr=-1
 write(*,'(a)')' Warning: the maximum number of searches was reached!'

end function
  
!PRIVATE
function rot(i,j) result(rota)
 integer(kind=int8) :: i,j
 integer(kind=int8)::rota

 rota=ior(ishft(i,j),ishft(i,-(32-j)))

end function
    
subroutine mix(a,b,c)
 integer(kind=int8),intent(inout)::a,b,c

 a=a-c;a=ieor(a,rot(c,4));c=c+b 
 b=b-a;b=ieor(b,rot(a,6));a=a+c 
 c=c-b;c=ieor(c,rot(b,8));b=b+a 
 a=a-c;a=ieor(a,rot(c,16));c=c+b 
 b=b-a;b=ieor(b,rot(a,19));a=a+c 
 c=c-b;c=ieor(c,rot(b,4));b=b+a 

end subroutine 

end module

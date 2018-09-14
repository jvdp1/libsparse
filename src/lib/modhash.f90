module modhash
 use modkind
 implicit none
 public::hashf

contains
!Inspired by lookup3.c from Bob Jenkins (http://burtleburtle.net/bob/hash/index.html#lookup)
!Converted in Fortran by Francois Guillaume - 2011
!Simplifications made by Jeremie Vandenplas - 2018

!PUBLIC
function hashf(row,col,mat,dim2,filled,mode) result(hashvr)
 !hashvr: address (column) of mat
 !mat of size dim1 (=2) x dim2
 !filled: number of elements
 !mode = 0: search for a value and returns 0 if absent
 !mode = 1: add a value if row,col was not present before
 integer(kind=int4)::hashvr
 integer(kind=int4),intent(in)::mode,row,col
 integer(kind=int4),intent(inout)::mat(:,:)
 integer(kind=int8),intent(in)::dim2
 integer(kind=int8),intent(inout)::filled

 integer(kind=int4),parameter::dim1=2
 integer(kind=int4),parameter::maxiter=5000

 integer(kind=int8)::a,b,c,iaddress,plage
 integer(kind=int4)::rowcol(dim1) 
 integer(kind=int4)::i,j,k 
 logical::indnotzero,indnotequal

 rowcol(1)=row
 rowcol(2)=col
  
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
  indnotzero=.false.;indnotequal=.false.
  do i=1,dim1
   j=mat(i,iaddress)
   if(j.ne.rowcol(i))indnotequal=.true.
   if(j.ne.0)indnotzero=.true.
  enddo     
  if(.not.indnotzero.or..not.indnotequal)then
   if(mode.eq.1.and..not.indnotzero)then
    do i=1,dim1
     mat(i,iaddress)=rowcol(i)
    enddo
    filled=filled+1
   endif
   if(mode.eq.0.and..not.indnotzero)then
    hashvr=0
   else
    hashvr=iaddress
   endif
   return
  endif
  !Hashing again
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

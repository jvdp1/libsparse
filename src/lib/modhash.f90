module modhash
 use modkind
 implicit none
 public::hashf,roundinguppower2

contains
!Inspired by lookup3.c from Bob Jenkins (http://burtleburtle.net/bob/hash/index.html#lookup)
!Converted in Fortran by Francois Guillaume - 2011
!Simplifications made by Jeremie Vandenplas - 2018

!PUBLIC
function hashf(row,col,mat,dim2,filled,getval) result(address)
 !address: address (column) of mat
 !mat of size dim1 (=2) x dim2
 !filled: number of elements
 !getval .eq. .true. : search for a value and returns 0 if absent
 !getval .eq. .false.: add a value if row,col was not present before
 integer(kind=int8)::address
 integer(kind=int4),intent(in)::row,col
 integer(kind=int4),intent(inout)::mat(:,:)
 integer(kind=int8),intent(in)::dim2
 integer(kind=int8),intent(inout)::filled
 logical,intent(in)::getval

 integer(kind=int4),parameter::maxiter=5000

 integer(kind=int8)::a,b,c
 integer(kind=int4)::i 
 logical::indzero,indequal

 a=int(row,kind(a))        !conversion of 1st coordinate
 b=int(col,kind(b))        !conversion of 2nd coordinate
 c=305419896_int8          !default value for 3rd coordinate  

 !Cycle until a free entry is found
 do i=1,maxiter
  !Hashing
  call mix(a,b,c)
  !Computation of the address
  address=iand(c,dim2-1)+1
  !Check if the address is correct
  indzero=.false.;indequal=.false.
  if(mat(1,address).eq.row.and.mat(2,address).eq.col)indequal=.true.
  if(mat(1,address).eq.0.or.mat(2,address).eq.0)indzero=.true.
  if(indzero.or.indequal)then
   if(.not.getval.and.indzero)then
    mat(1,address)=row
    mat(2,address)=col
    filled=filled+1
    return
   endif
   if(getval.and.indzero)then
    address=0
   endif
   return
  endif
 enddo

 address=-1
 write(*,'(a)')' Warning: the maximum number of searches was reached!'

end function
  
function roundinguppower2(x) result(next)
 integer(kind=int8),intent(in)::x
 integer(kind=int8)::next

 !https://stackoverflow.com/questions/466204/rounding-up-to-next-power-of-2
 
 next=2_int8**int((ceiling(log(real(x,real8))/log(real(2,real8)))),int8)

end function

!PRIVATE
function rot(i,j) result(rota)
 integer(kind=int8) :: i,j
 integer(kind=int8)::rota

 rota=ior(ishft(i,j),ishft(i,-(32_int8-j)))

end function
    
subroutine mix(a,b,c)
 integer(kind=int8),intent(inout)::a,b,c

 a=a-c;a=ieor(a,rot(c,4_int8));c=c+b 
 b=b-a;b=ieor(b,rot(a,6_int8));a=a+c 
 c=c-b;c=ieor(c,rot(b,8_int8));b=b+a 
 a=a-c;a=ieor(a,rot(c,16_int8));c=c+b 
 b=b-a;b=ieor(b,rot(a,19_int8));a=a+c 
 c=c-b;c=ieor(c,rot(b,4_int8));b=b+a 

end subroutine 

end module

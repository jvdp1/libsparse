module modhash
 use iso_fortran_env
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
 integer(kind=int64)::address
 integer(kind=int32),intent(in)::row,col
 integer(kind=int32),intent(inout)::mat(:,:)
 integer(kind=int64),intent(in)::dim2
 integer(kind=int64),intent(inout)::filled
 logical,intent(in)::getval

 integer(kind=int32),parameter::maxiter=5000

 integer(kind=int64)::a,b,c
 integer(kind=int32)::i 
 logical::indzero,indequal

 a=int(row,kind(a))        !conversion of 1st coordinate
 b=int(col,kind(b))        !conversion of 2nd coordinate
 c=305419896_int64          !default value for 3rd coordinate  

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
 integer(kind=int64),intent(in)::x
 integer(kind=int64)::next

 !https://stackoverflow.com/questions/466204/rounding-up-to-next-power-of-2
 
 next=2_int64**int((ceiling(log(real(x,real64))/log(real(2,real64)))),int64)

end function

!PRIVATE
function rot(i,j) result(rota)
 integer(kind=int64) :: i,j
 integer(kind=int64)::rota

 rota=ior(ishft(i,j),ishft(i,-(32_int64-j)))

end function
    
subroutine mix(a,b,c)
 integer(kind=int64),intent(inout)::a,b,c

 a=a-c;a=ieor(a,rot(c,4_int64));c=c+b 
 b=b-a;b=ieor(b,rot(a,6_int64));a=a+c 
 c=c-b;c=ieor(c,rot(b,8_int64));b=b+a 
 a=a-c;a=ieor(a,rot(c,16_int64));c=c+b 
 b=b-a;b=ieor(b,rot(a,19_int64));a=a+c 
 c=c-b;c=ieor(c,rot(b,4_int64));b=b+a 

end subroutine 

end module

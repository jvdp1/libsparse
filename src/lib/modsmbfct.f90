module modsmbfct
#if (_dp==0)
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real32
#else
 use iso_fortran_env,only:output_unit,int32,int64,real32,real64,wp=>real64
#endif
 implicit none
 private
 public::smbfctf90

contains

subroutine smbfctf90(neqns, xadj, adjncy, perm, invp,&
                      xlnz, maxlnz, xnzsub, nzsub, maxsub,&
                      rchlnk, mrglnk, marker, flag)
 integer(kind=int32),intent(in)::neqns,xadj(:),adjncy(:),perm(:),invp(:)
 integer(kind=int32),intent(out)::xlnz(:),xnzsub(:),nzsub(:)
 integer(kind=int32),intent(out)::maxlnz,flag
 integer(kind=int32),intent(out)::rchlnk(neqns),mrglnk(neqns),marker(neqns)
 integer(kind=int32),intent(inout)::maxsub

 integer(kind=int32)::i,j,k,m
 integer(kind=int32)::nzbeg,nzbend,np1,knz,mrgk,mrkflg,nzend,node,jstrt,jstop,nabor,rchm,lmax,inz,kxsub
 logical::l_350,l_1200

 !initialization
 nzbeg = 1
 nzbend = 0
 xlnz(1) = 1
 mrglnk = 0
 marker = 0

 !for each column:  knz counts the number of nonzeros in column k accumulated in rchlnk.
 np1 = neqns + 1

 do k = 1, neqns
  knz = 0
  mrgk = mrglnk(k)
  mrkflg = 0
  marker(k) = k
  if (mrgk .ne. 0) marker(k) = marker(mrgk)
  xnzsub(k) = nzend
  node = perm(k)
  jstrt = xadj(node)
  jstop = xadj(node+1) - 1
  if(jstrt.le.jstop)then
   !use rchlnk to link through the structure of a(*,k) below diagonal
   rchlnk(k) = np1
   do_300: do j = jstrt, jstop
    nabor = adjncy(j)
    nabor = invp(nabor)
    if ( nabor .le. k ) exit do_300
    rchm = k
    do_200: do
     m = rchm
     rchm = rchlnk(m)
     if ( rchm .gt. nabor ) exit do_200
    enddo do_200
    knz = knz+1
    rchlnk(m) = nabor
    rchlnk(nabor) = rchm
    if ( marker(nabor) .ne. marker(k) ) mrkflg = 1
   enddo do_300
   !test for mass symbolic elimination ...
   lmax = 0
   l_350=.false.
   if(mrkflg.ne.0 .or. mrgk.eq.0)l_350=.true.
   if(.not.l_350 .and. mrglnk(mrgk).ne.0)l_350=.true.
   !link through each column i that affects l(*,k).
   if(l_350)then
    i = k
    do_400: do
     i = mrglnk(i)
     if (i.eq.0) exit do_400
     inz = xlnz(i+1) - (xlnz(i)+1)
     jstrt = xnzsub(i) +  1
     jstop = xnzsub(i) + inz
     if (inz.gt.lmax)then
      lmax = inz
      xnzsub(k) = jstrt
     endif
     !merge structure of l(*,i) in nzsub into rchlnk.
     rchm = k
     do_700: do j = jstrt, jstop
      nabor = nzsub(j)
      do_600: do
       m = rchm
       rchm = rchlnk(m)
       if(rchm.ge.nabor) exit do_600
      enddo do_600
      if (rchm.eq.nabor) exit do_700
      knz = knz+1
      rchlnk(m) = nabor
      rchlnk(nabor) = rchm
      rchm = nabor
     enddo do_700
    enddo do_400
    !check if subscripts duplicate those of another column.
    if (knz.ne.lmax) then
     !or if tail of k-1st column matches head of kth.
     if (nzbeg.le.nzend)then
      i = rchlnk(k)
      l_1200 = .false.
      do_900: do jstrt=nzbeg,nzend
       if ((nzsub(jstrt)-i).lt.0) then
        cycle do_900
       elseif ((nzsub(jstrt)-i).eq.0)then
        xnzsub(k) = jstrt
        do j=jstrt,nzend
         if (nzsub(j).ne.i)then
          l_1200=.true.
          exit do_900
         endif
         i = rchlnk(i)
         if (i.gt.neqns) exit do_900
        enddo
        nzend = jstrt - 1
       elseif ((nzsub(jstrt)-i) .gt.0)then
        l_1200=.true.
        exit do_900
       endif
      enddo do_900
     endif
     if(l_1200)then
      !copy the structure of l(*,k) from rchlnk to the data structure (xnzsub, nzsub).
      nzbeg = nzend +  1
      nzend = nzend + knz
      if (nzend.gt.maxsub)then
       !error - insufficient storage for nonzero subscripts.
       flag = 1
       return
      endif
      i = k
      do j=nzbeg,nzend
       i = rchlnk(i)
       nzsub(j) = i
       marker(i) = k
      enddo
      xnzsub(k) = nzbeg
      marker(k) = k
     endif
    endif
   else
    xnzsub(k) = xnzsub(mrgk) + 1
    knz = xlnz(mrgk+1) - (xlnz(mrgk) + 1)
   endif
   !update the vector mrglnk.  note column l(*,k) just found is required to determine column l(*,j), where l(j,k) is the first nonzero in l(*,k) below diagonal.
   if(knz.gt.1)then
    kxsub = xnzsub(k)
    i = nzsub(kxsub)
    mrglnk(k) = mrglnk(i)
    mrglnk(i) = k
   endif
  endif
  xlnz(k+1) = xlnz(k) + knz
 enddo

 maxlnz = xlnz(neqns) - 1
 maxsub = xnzsub(neqns)
 xnzsub(neqns+1) = xnzsub(neqns)
 flag = 0

end subroutine

end module

submodule (modsparse) modsparse_metisgraph
 !$ use omp_lib
 implicit none

contains
!**CONSTRUCTOR
module subroutine constructor_sub_metisgraph(metis,n,m,unlog)
 class(metisgraph),intent(out)::metis
 integer(kind=int32),intent(in)::n,m
 integer(kind=int32),intent(in),optional::unlog

 metis%unlog=6
 if(present(unlog))metis%unlog=unlog

 metis%nvertices=n
 metis%medges=m
 allocate(metis%xadj(metis%nvertices+1),metis%adjncy(2*metis%medges))
 metis%xadj=0
 metis%adjncy=0
 metis%vwgt=c_null_ptr
 metis%adjwgt=c_null_ptr

end subroutine


!** GET MEMORY
module function getmem_metisgraph(metis) result(getmem)
 class(metisgraph),intent(in)::metis
 integer(kind=int64)::getmem

 getmem=sizeof(metis%unlog)+sizeof(metis%nvertices)+sizeof(metis%medges)
 if(allocated(metis%xadj))getmem=getmem+sizeof(metis%xadj)
 if(allocated(metis%adjncy))getmem=getmem+sizeof(metis%adjncy)

 !what to do with the pointers???

end function


end submodule

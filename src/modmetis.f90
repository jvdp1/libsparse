!> Module containing interfaces for METIS 5

!> @todo Currently supports only ordering. Should be extended to the other options.

module modmetis
 !based on https://glaros.dtc.umn.edu/gkhome/node/877
 use iso_fortran_env,only:int32
 use iso_c_binding,only:c_int,c_ptr
 implicit none
 private
 public::metis_nodend,metis_setoptions,metis_checkerror

 integer(kind=c_int),parameter::METIS_NOPTIONS=40

 !options, possible values
 integer(kind=c_int),parameter,public::METIS_OPTION_PTYPE=0      ,METIS_PTYPE_RB=0,METIS_PTYPE_KWAY=1
 integer(kind=c_int),parameter,public::METIS_OPTION_OBJTYPE=1    ,METIS_OBJTYPE_CUT=1,METIS_OBJTYPE_VOL=1 
 integer(kind=c_int),parameter,public::METIS_OPTION_CTYPE=2      ,METIS_CTYPE_RM=0,METIS_CTYPE_SHEM=1
 integer(kind=c_int),parameter,public::METIS_OPTION_IPTYPE=3     ,METIS_IPTYPE_GROW=0,METIS_IPTYPE_RANDOM=1&
                                                                 ,METIS_IPTYPE_EDGE=2,METIS_IPTYPE_NODE=3
 integer(kind=c_int),parameter,public::METIS_OPTION_RTYPE=4      ,METIS_RTYPE_FM=0,METIS_RTYPE_GREEDY=1&
                                                                 ,METIS_RTYPE_SEP2SIDED=2,METIS_RTYPE_SEP1SIDED=3
 integer(kind=c_int),parameter,public::METIS_OPTION_DBGLVL=5     ,METIS_DBG_INFO=1,METIS_DBG_TIME=2&
                                                                 ,METIS_DBG_COARSEN=4,METIS_DBG_REFINE=8&
                                                                 ,METIS_DBG_IPART=16,METIS_DBG_MOVEINFO=32,METIS_DBG_SEPINFO=64&
                                                                 ,METIS_DBG_CONNINFO=128,METIS_DBG_CONTIGINFO=256
 integer(kind=c_int),parameter,public::METIS_OPTION_NITER=6      !default: 10
 integer(kind=c_int),parameter,public::METIS_OPTION_NCUTS=7      !default: 1
 integer(kind=c_int),parameter,public::METIS_OPTION_SEED=8      
 integer(kind=c_int),parameter,public::METIS_OPTION_NO2HOP=9     !0 or 1
 integer(kind=c_int),parameter,public::METIS_OPTION_MINCONN=10   !0 or 1
 integer(kind=c_int),parameter,public::METIS_OPTION_CONTIG=11    !0 or 1
 integer(kind=c_int),parameter,public::METIS_OPTION_COMPRESS=12  !0 or 1
 integer(kind=c_int),parameter,public::METIS_OPTION_CCORDER=13   !0 or 1
 integer(kind=c_int),parameter,public::METIS_OPTION_PFACTOR=14   
 integer(kind=c_int),parameter,public::METIS_OPTION_NSEPS=15     !default: 1
 integer(kind=c_int),parameter,public::METIS_OPTION_UFACTOR=16   !default: 1 or 30
 integer(kind=c_int),parameter,public::METIS_OPTION_NUMBERING=17 !0 or 1

 !error values
 integer(kind=c_int),parameter::METIS_OK=1,METIS_ERROR_INPUT=-2,METIS_ERROR_MEMORY=-3,METIS_ERROR=-4

 interface
  function metis_setdefaultoptions(options) result(err) bind(C,name='METIS_SetDefaultOptions')
   import c_int,METIS_NOPTIONS
   integer(kind=c_int),intent(inout)::options(0:METIS_NOPTIONS)
   integer(kind=c_int)::err
  end function
  !METIS_NodND
  !OTPIONS Version 4                   -> Version 5
  ![1] 1(RM) 2(HEM) 3(SHEM)            -> CTYPE
  ![2] 1(edge-based) 2(node-based)     -> ?
  ![3] 1(2-sided node) 2(1-sided node) -> RTYPE
  ![4] 0
  ![5] 0(no compress+no order) 1(compress) 2(order) 3(compress + order) -> COMPRESS +  CCORDER
  ![6] remove vertices                 -> PFACTOR
  ![7] separators                      -> NSEPS
  function metis_nodend(nvtxs,xadj,adjncy,vwgt,options,perm,iperm) result(err) bind(C,name='METIS_NodeND')
   import c_int,c_ptr,METIS_NOPTIONS
   integer(kind=c_int),intent(in)::nvtxs
   integer(kind=c_int),intent(in)::xadj(*),adjncy(*)
   type(c_ptr),intent(in),value::vwgt
   integer(kind=c_int),intent(in)::options(0:METIS_NOPTIONS)   !options is mandatory for Fortran when array start at pos 1
   integer(kind=c_int),intent(out)::perm(*),iperm(*)
   integer(kind=c_int)::err
  end function
 end interface

contains

function metis_setoptions(options&
                           ,ptype,objtype,ctype,iptype,rtype,ncuts&
                           ,nseps,niter,seed,minconn,no2hop,contig&
                           ,compress,ccorder,pfactor,ufactor,dbglvl&
                           ) result(err)
 integer(kind=c_int),allocatable,intent(out)::options(:)
 integer(kind=c_int),intent(in),optional::ptype,objtype,ctype,iptype,rtype,ncuts,&
                                           nseps,niter,seed,minconn,no2hop,&
                                           contig,compress,ccorder,pfactor,ufactor,dbglvl
 integer(kind=c_int)::err

 if(allocated(options))deallocate(options)
 allocate(options(0:METIS_NOPTIONS))

 err=metis_setdefaultoptions(options)
 !DEFAULT
 options(METIS_OPTION_NUMBERING)=1

 !OPTIONAL
 if(present(ptype))options(METIS_OPTION_PTYPE)=ptype
 if(present(objtype))options(METIS_OPTION_OBJTYPE)=objtype
 if(present(ctype))options(METIS_OPTION_CTYPE)=ctype
 if(present(iptype))options(METIS_OPTION_IPTYPE)=iptype
 if(present(rtype))options(METIS_OPTION_RTYPE)=rtype
 if(present(ncuts))options(METIS_OPTION_NCUTS)=ncuts
 if(present(nseps))options(METIS_OPTION_NSEPS)=nseps
 if(present(niter))options(METIS_OPTION_NITER)=niter
 if(present(seed))options(METIS_OPTION_SEED)=seed
 if(present(minconn))options(METIS_OPTION_MINCONN)=minconn
 if(present(no2hop))options(METIS_OPTION_NO2HOP)=no2hop
 if(present(contig))options(METIS_OPTION_CONTIG)=contig
 if(present(compress))options(METIS_OPTION_COMPRESS)=compress
 if(present(ccorder))options(METIS_OPTION_CCORDER)=ccorder
 if(present(pfactor))options(METIS_OPTION_PFACTOR)=pfactor
 if(present(ufactor))options(METIS_OPTION_UFACTOR)=ufactor
 if(present(dbglvl))options(METIS_OPTION_DBGLVL)=dbglvl

end function

subroutine metis_checkerror(err,unlog)
 integer(kind=c_int),intent(in)::err
 integer(kind=int32),intent(in),optional::unlog

 integer(kind=int32)::un

 un=6
 if(present(unlog))un=unlog
 
 select case(err)
  case(METIS_OK)
   !write(un,'(/a/)')' METIS_OK'
  case(METIS_ERROR_INPUT)
   write(un,'(/a/)')' METIS_ERROR_INPUT'
   error stop
  case(METIS_ERROR_MEMORY)
   write(un,'(/a/)')' METIS_ERROR_MEMORY'
   error stop
  case(METIS_ERROR)
   write(un,'(/a/)')' METIS_ERROR'
   error stop
  case default
   write(un,'(/a/)')' UNKNOWN METIS_ERROR'
   error stop
 end select
 
end subroutine

end module

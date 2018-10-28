module modmetis
 !based on https://glaros.dtc.umn.edu/gkhome/node/877
 use iso_fortran_env,only:int32
 use iso_c_binding,only:c_int,c_ptr
 implicit none
 private
 public::METIS_NodeND,METIS_SetOptions,METIS_CheckError

 integer(kind=c_int),parameter::METIS_NOPTIONS=40

 !options, possible values
 integer(kind=c_int),parameter,public::METIS_OPTION_CTYPE=2      ,METIS_CTYPE_RM=0,METIS_CTYPE_SHEM=1
 integer(kind=c_int),parameter,public::METIS_OPTION_RTYPE=4      ,METIS_RTYPE_FM=0,METIS_RTYPE_GREEDY=1&
                                                                 ,METIS_RTYPE_SEP2SIDED=2,METIS_RTYPE_SEP1SIDED=3
 integer(kind=c_int),parameter,public::METIS_OPTION_NO2HOP=9     !0 or 1
 integer(kind=c_int),parameter,public::METIS_OPTION_NSEPS=15     !default: 1
 integer(kind=c_int),parameter,public::METIS_OPTION_NITER=6      !default: 10
 integer(kind=c_int),parameter,public::METIS_OPTION_UFACTOR=16   !default: 1 or 30
 integer(kind=c_int),parameter,public::METIS_OPTION_COMPRESS=12  !0 or 1
 integer(kind=c_int),parameter,public::METIS_OPTION_CCORDER=13   !0 or 1
 integer(kind=c_int),parameter,public::METIS_OPTION_SEED=8      
 integer(kind=c_int),parameter,public::METIS_OPTION_PFACTOR=14   
 integer(kind=c_int),parameter,public::METIS_OPTION_NUMBERING=17 !0 or 1
 integer(kind=c_int),parameter,public::METIS_OPTION_DBGLVL=5     ,METIS_DBG_INFO=1,METIS_DBG_TIME=2&
                                                                 ,METIS_DBG_COARSEN=4,METIS_DBG_REFINE=8&
                                                                 ,METIS_DBG_IPART=16,METIS_DBG_MOVEINFO=32,METIS_DBG_SEPINFO=64&
                                                                 ,METIS_DBG_CONNINFO=128,METIS_DBG_CONTIGINFO=256

 !error values
 integer(kind=c_int),parameter::METIS_OK=1,METIS_ERROR_INPUT=-2,METIS_ERROR_MEMORY=-3,METIS_ERROR=-4

 interface
  function METIS_SetDefaultOptions(options) result(err) bind(C,name='METIS_SetDefaultOptions')
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
  function METIS_NodeND(nvtxs,xadj,adjncy,vwgt,options,perm,iperm) result(err) bind(C,name='METIS_NodeND')
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

function METIS_SetOptions(options,dbglvl) result(err)
 integer(kind=c_int),allocatable,intent(out)::options(:)
 integer(kind=c_int),intent(in),optional::dbglvl
 integer(kind=c_int)::err

 if(allocated(options))deallocate(options)
 allocate(options(0:METIS_NOPTIONS))

 err=METIS_SetDefaultOptions(options)
 !DEFAULT
 options(METIS_OPTION_NUMBERING)=1

 !OPTIONAL
 if(present(dbglvl))options(METIS_OPTION_DBGLVL)=dbglvl

end function

subroutine METIS_CheckError(err,unlog)
 integer(kind=c_int),intent(in)::err
 integer(kind=int32),intent(in),optional::unlog

 integer(kind=int32)::un

 un=6
 if(present(unlog))un=unlog
 
 select case(err)
  case(METIS_OK)
   write(un,'(/a/)')' METIS_OK'
  case(METIS_ERROR_INPUT)
   write(un,'(/a/)')' METIS_ERROR_INPUT'
   stop
  case(METIS_ERROR_MEMORY)
   write(un,'(/a/)')' METIS_ERROR_MEMORY'
   stop
  case(METIS_ERROR)
   write(un,'(/a/)')' METIS_ERROR'
   stop
  case default
   write(un,'(/a/)')' UNKNOWN METIS_ERROR'
   stop
 end select
 
end subroutine

end module

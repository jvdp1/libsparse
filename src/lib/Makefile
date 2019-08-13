# Makes libraries

include ../Makeinit

OBJ=

ifeq ($(DEBUGENABLE),1)
 DEBUG=-check all -traceback -debug extended
else
 DEBUG=
endif

ifeq ($(DPENABLE),0)
 DP=0
else
 DP=1
endif

ifeq ($(PARDISOENABLE),0)
 PARDISO=0
 F90FLAGS+=
else
 PARDISO=1
 F90FLAGS+=
endif

ifeq ($(METISENABLE),0)
 METIS=0
 LIBMETISROOT=
 LIBMETIS=
 LMETIS=
else
 METIS=1
 F90FLAGS+=-I$(LIBMETISROOT)$(LIBMETIS)
 OBJ+=modmetis.o
endif

ifeq ($(SPAINVENABLE),0)
 SPAINV=0
else
 SPAINV=1
 OBJ+=smbfct.o modspainv.o
endif


all: libsparse.a

OBJ+=modhash.o modsparse.o

libsparse.a: $(OBJ)
	ar cr  libsparse.a $(OBJ)

%.o:	%.f90
	$(f90) -c  $(F90FLAGS) $(DEBUG) -fpp -D_PARDISO=$(PARDISO) -D_DP=$(DP) -D_METIS=$(METIS) -D_SPAINV=$(SPAINV) -D_VERBOSE=$(VERBOSE) $<

%.o:	%.f
	$(f90) -c  $(F90FLAGS) $(DEBUG) -fpp -D_PARDISO=$(PARDISO) -D_DP=$(DP) -D_METIS=$(METIS) -D_SPAINV=$(SPAINV) -D_VERBOSE=$(VERBOSE) $<

clean:
	rm -fv *.o *.optrpt *.a *.mod
	
cleanlib:
	rm -fv *.a
	

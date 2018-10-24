# Makes libraries

include ../Makeinit

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


all: libsparse.a

OBJ=modhash.o modsparse.o

libsparse.a: $(OBJ)
	ar cr  libsparse.a $(OBJ)

%.o:	%.f90 
	$(f90) -c  $(F90FLAGS) $(DEBUG) -fpp -D_PARDISO=$(PARDISO) -D_DP=$(DP) $<

clean:
	rm -fv *.o *.optrpt *.a *.mod
	
cleanlib:
	rm -fv *.a
	

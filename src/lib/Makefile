# Makes libraries

include ../Makeinit

ifeq ($(DEBUGENABLE),1)
 DEBUG=-check all -traceback -debug extended
else
 DEBUG=
endif

PARDISOENABLE=1
ifeq ($(PARDISOENABLE),1)
 PARDISO=1
 F90FLAGS+=
endif

all: libsparse.a

OBJ=modkind.o modhash.o modsparse.o

libsparse.a: $(OBJ)
	ar cr  libsparse.a $(OBJ)

%.o:	%.f90 
	$(f90) -c  $(F90FLAGS) $(DEBUG) -fpp -D_PARDISO=$(PARDISO) $<

clean:
	rm -fv *.o *.optrpt *.a *.mod
	
cleanlib:
	rm -fv *.a
	

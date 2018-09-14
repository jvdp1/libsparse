# Makes libraries

include ../Makeinit

ifeq ($(DEBUGENABLE),1)
 DEBUG=-check all -traceback -debug extended
else
 DEBUG=
endif


all: libsparse.a

OBJ=modkind.o modhash.o modsparse.o

libsparse.a: $(OBJ)
	ar cr  libsparse.a $(OBJ)

%.o:	%.f90 
	$(f90) -c  $(optf90) $(DEBUG) $<

clean:
	rm -fv *.o *.optrpt *.a *.mod
	
cleanlib:
	rm -fv *.a
	

# Makes libraries

OBJ=modcommon.o

ifeq ($(PARDISOENABLE),0)
 PARDISO=0
 FFLAGS+=
else
 PARDISO=1
 FFLAGS+=
 OBJ+=modvariablepardiso.o
endif

ifeq ($(METISENABLE),0)
else
 OBJ+=modmetis.o
endif

ifeq ($(SPAINVENABLE),0)
else
 OBJ+=sgtrsm.o dgtrsm.o smbfct.o modspainv.o
endif


all: libsparse.a

OBJ+=modhash.o modsparse.o modsparse_gen.o modsparse_coo.o modsparse_crs.o modsparse_ll.o modsparse_metisgraph.o

libsparse.a: $(OBJ)
	ar cr  libsparse.a $(OBJ)

%.o:	%.f90
	$(FC) -c  $(FFLAGS) $(DEBUG) -fpp -D_PARDISO=$(PARDISO) $<

%.o:	%.f
	$(FC) -c  $(FFLAGS) $(DEBUG) -fpp -D_PARDISO=$(PARDISO) $<

clean:
	rm -fv *.o *.optrpt *.a *.mod *.smod
	
cleanlib:
	rm -fv *.a
	

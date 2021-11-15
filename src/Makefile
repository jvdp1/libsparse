# Makes libraries

ifeq ($(PARDISOENABLE),1)
 OBJPARDISO = modvariablepardiso.o
endif

ifeq ($(METISENABLE),1)
 OBJMETIS = modmetis.o
endif

ifeq ($(SPAINVENABLE),1)
 OBJSPAINV = sgtrsm.o dgtrsm.o smbfct.o modspainv.o
endif


all: libsparse.a

OBJ = modcommon.o \
      $(OBJPARDISO) $(OBJMETIS) $(OBJSPAINV) \
      modhash.o modsparse.o modsparse_gen.o modsparse_coo.o \
      modsparse_crs.o modsparse_ll.o modsparse_metisgraph.o

libsparse.a: $(OBJ)
	ar cr  libsparse.a $(OBJ)

%.o:	%.f90
	$(FC) -c  $(FFLAGS) $(DEBUG) $<

%.o:	%.f
	$(FC) -c  $(FFLAGS) $(DEBUG) $<

clean:
	rm -fv *.o *.optrpt *.a *.mod *.smod
	
cleanlib:
	rm -fv *.a


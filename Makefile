#Makefile

DEBUGENABLE=0
DPENABLE=1
METISENABLE=0
PARDISOENABLE=1
SPAINVENABLE=1
VERBOSE=0

FC = gfortran
#FC = ifort

#ifort
ifeq ($(FC), ifort)
FFLAGS=-O3 -fpp -heap-arrays -qopenmp -parallel -qopt-matmul -qopt-report=5
FFLAGS+=-I${MKLROOT}/include -I${MKLROOT}/include/intel64/lp64 -L${MKLROOT}/lib/intel64

FLIBS += -liomp5 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core 
FLIBS += -lpthread -lm -ldl

ifeq ($(DEBUGENABLE), 1)
FFLAGS += -g -check all -traceback -warn all -stand f18
endif

endif


#gfortan
ifeq ($(FC), gfortran)
FFLAGS=-O3 -cpp -fopenmp -fall-intrinsics
FFLAGS+=-I${MKLROOT}/include -I${MKLROOT}/include/intel64/lp64 -L${MKLROOT}/lib/intel64

FLIBS += -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp
FLIBS += -lpthread -lm -ldl

ifeq ($(DEBUGENABLE), 1)
FFLAGS += -g -Wall -fcheck=all -fbacktrace -std=f2008
endif

endif


ifeq ($(METISENABLE), 1)
 LIBMETISROOT=~/metis-5.1.0
 LIBMETIS=$(LIBMETISROOT)/build/Linux-x86_64/libmetis
 FLIBS += $(LIBMETIS)/libmetis.a
 METIS = 1
else
 METIS = 0
endif


ifeq ($(DPENABLE),0)
 DP=0
else
 DP=1
endif

ifeq ($(PARDISOENABLE),1)
 PARDISO=1
else
 PARDISO=0
endif


ifeq ($(SPAINVENABLE),1)
 SPAINV=1
else
 SPAINV=0
endif

FFLAGS += -D_DP=$(DP) -D_METIS=$(METIS) -D_PARDISO=$(PARDISO) -D_SPAINV=$(SPAINV) -D_VERBOSE=$(VERBOSE)

FYPPFLAGS =

export FC
export FFLAGS
export FLIBS
export FYPPFLAGS

export METISENABLE
export PARDISOENABLE
export SPAINVENABLE

.PHONY: all clean examples lib test

all: lib
	$(MAKE) --directory=test

examples:
	$(MAKE) --directory=examples

lib:
	$(MAKE) --directory=src -j

test:
	$(MAKE) --directory=test test
	@echo
	@echo "All tests passed."

clean:
	$(MAKE) clean --directory=src
	$(MAKE) clean --directory=examples
	$(MAKE) clean --directory=test

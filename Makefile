#Makefile

DEBUGENABLE=1
DPENABLE=1
METISENABLE=1
PARDISOENABLE=1
SPAINVENABLE=1
VERBOSE=0

FC = ifort
#FC = gfortran

#FFLAGS = -g -Wall -std=f2008
FFLAGS=-O3 -heap-arrays -qopenmp -liomp5 -parallel -qopt-matmul -qopt-report=5
FFLAGS+=-I${MKLROOT}/include -I${MKLROOT}/include/intel64/lp64 -L${MKLROOT}/lib/intel64

FLIBS = -lpthread -ldl
FLIBS += -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core 

ifeq ($(METISENABLE), 1)
 LIBMETISROOT=~/metis-5.1.0
 LIBMETIS=$(LIBMETISROOT)/build/Linux-x86_64/libmetis
 FFLAGS += -D_METIS=$(METISENABLE)
 FLIBS += $(LIBMETIS)/libmetis.a
endif

ifeq ($(DEBUGENABLE), 1)
FFLAGS += -g -check all -traceback -stand f08
endif

ifeq ($(DPENABLE),0)
 DP=0
else
 DP=1
endif

ifeq ($(PARDISOENABLE),0)
 PARDISO=0
else
 PARDISO=1
endif


ifeq ($(SPAINVENABLE),0)
 SPAINV=0
else
 SPAINV=1
endif

FFLAGS += -fpp -D_DP=$(DP)  -D_PARDISO=$(PARDISO) -D_SPAINV=$(SPAINV) -D_VERBOSE=$(VERBOSE)

FYPPFLAGS=

export FC
export FFLAGS
export FLIBS
export FYPPFLAGS

.PHONY: all clean test

all:
	$(MAKE) --directory=src/lib
	$(MAKE) --directory=tests

example:
	$(MAKE) --directory=src/test

test:
	$(MAKE) --directory=tests test
	@echo
	@echo "All tests passed."

clean:
	$(MAKE) clean --directory=src/lib
	$(MAKE) clean --directory=src/test
	$(MAKE) clean --directory=tests

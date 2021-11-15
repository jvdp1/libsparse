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
 FLIBS += $(LIBMETIS)/libmetis.a
 METIS = 1
else
 METIS = 0
endif

ifeq ($(DEBUGENABLE), 1)
FFLAGS += -g -check all -traceback -stand f08
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

FFLAGS += -fpp -D_DP=$(DP) -D_METIS=$(METIS) -D_PARDISO=$(PARDISO) -D_SPAINV=$(SPAINV) -D_VERBOSE=$(VERBOSE)

FYPPFLAGS =

MAKEFLAGS = METISENABLE=$(METISENABLE) PARDISOENABLE=$(PARDISOENABLE) SPAINVENABLE=$(SPAINVENABLE)

export FC
export FFLAGS
export FLIBS
export FYPPFLAGS

.PHONY: all clean test

all:
	$(MAKE) --directory=src $(MAKEFLAGS)
	$(MAKE) --directory=tests

example:
	$(MAKE) --directory=src/test $(MAKEFLAGS)

test:
	$(MAKE) --directory=tests test
	@echo
	@echo "All tests passed."

clean:
	$(MAKE) clean --directory=src $(MAKEFLAGS)
	$(MAKE) clean --directory=src/test $(MAKEFLAGS)
	$(MAKE) clean --directory=tests

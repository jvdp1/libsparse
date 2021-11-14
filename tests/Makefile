include ../Makeinit

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

ifeq ($(METISENABLE),0)
 METIS=0
 LIBMETISROOT=
 LIBMETIS=
 LMETIS=
else
 METIS=1
 F90FLAGS+=-I$(LIBMETISROOT)$(LIBMETIS)
endif

ifeq ($(SPAINVENABLE),0)
 SPAINV=0
else
 SPAINV=1
endif

PROGS_SRC = test_sparse.f90

LIBTEST = libtestrive.a
OBJTEST = testdrive.o

NAME = sparse

CPPFLAGS += -I../src/lib/ -I.

CPPFLAGS += -fpp -D_PARDISO=$(PARDISO) -D_DP=$(DP) -D_METIS=$(METIS) -D_SPAINV=$(SPAINV) -D_VERBOSE=$(VERBOSE) 

OBJS = modtest_common.o modtest_coo.o modtest_crs.o

LIB = $(patsubst %, ../src/lib/lib%.a, $(NAME))

PROGS = $(PROGS_SRC:.f90=)
TESTPROGS = $(PROGS:=TEST)

F90FLAGS += -check all -traceback -debug extended -stand f18

.PHONY: all clean test $(TESTPROGS)

all: $(PROGS)

test: $(TESTPROGS)

$(TESTPROGS):
	./$(@:TEST=)

clean:
	$(RM) $(PROGS) $(OBJS) $(CLEAN_FILES) *.o* *.mod *.a testdrive.f90

$(LIBTEST): $(OBJTEST)
	ar rcs $@ $^

%.o: %.f90 $(LIBTEST)
	$(f90) $(F90FLAGS) $(CPPFLAGS) -c $<

$(PROGS): %: $(OBJS) $(LIBTEST) $(LIB)  %.o
	$(f90) $(F90FLAGS) $(CPPFLAGS) -o $@ $^ $(LIBMETISROOT)$(LIBMETIS)$(LMETIS)

testdrive.f90: ../test-drive/src/testdrive.f90
	cp ../test-drive/src/testdrive.f90 .

test_coo.o: modtest_coo.o
modtest_coo.o: modtest_common.o
modtest_crs.o: modtest_common.o
testdrive.o: testdrive.f90

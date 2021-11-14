ifeq ($(PARDISOENABLE),0)
 PARDISO=0
else
 PARDISO=1
endif

PROGS_SRC = test_sparse.f90

LIBTEST = libtestrive.a
OBJTEST = testdrive.o

NAME = sparse

CPPFLAGS += -I../src/lib/ -I.

FFLAGS += -fpp -D_PARDISO=$(PARDISO)

OBJS = modtest_common.o modtest_coo.o modtest_crs.o

LIB = $(patsubst %, ../src/lib/lib%.a, $(NAME))

PROGS = $(PROGS_SRC:.f90=)
TESTPROGS = $(PROGS:=TEST)

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
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $<

$(PROGS): %: $(OBJS) $(LIBTEST) $(LIB)  %.o
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(FLIBS)

testdrive.f90: ../test-drive/src/testdrive.f90
	cp ../test-drive/src/testdrive.f90 .

test_coo.o: modtest_coo.o
modtest_coo.o: modtest_common.o
modtest_crs.o: modtest_common.o
testdrive.o: testdrive.f90

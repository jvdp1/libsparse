PROGS_SRC = test_sparse.f90

PATHTESST = ../test-drive/src
LIBTEST = libtestrive.a
OBJTEST = testdrive.o

NAME = sparse

CPPFLAGS += -I../src/ -I.

OBJS = modtest_common.o modtest_coo.o modtest_crs.o modtest_random.o

LIB = $(patsubst %, ../src/lib%.a, $(NAME))

PROGS = $(PROGS_SRC:.f90=)
TESTPROGS = $(PROGS:=TEST)

.PHONY: all clean test $(TESTPROGS)

all: $(PROGS)

test: $(TESTPROGS)

$(TESTPROGS):
	./$(@:TEST=)

clean:
	$(RM) $(PROGS) $(OBJS) $(CLEAN_FILES) *.o* *.*mod *.a

$(LIBTEST): $(OBJTEST)
	ar rcs $@ $^

testdrive.o: $(PATHTESST)/testdrive.F90
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $<

%.o: %.f90 $(LIBTEST)
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $<

$(PROGS): %: $(OBJS) $(LIBTEST) $(LIB)  %.o
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(FLIBS)

modtest_common.o: $(LIB)
modtest_coo.o: modtest_common.o $(LIB) testdrive.o
modtest_crs.o: modtest_common.o $(LIB) testdrive.o
modtest_random.o: $(LIB) testdrive.o
test_coo.o: modtest_coo.o modtest_crs.o modtest_random.o testdrive.o

PROGS_SRC = test_sparse.f90

PATHTESST = ../test-drive/src
LIBTEST = libtestrive.a
OBJTEST = testdrive.o

NAME = sparse

CPPFLAGS += -I../src/ -I.

OBJS = modtest_common.o modtest_coo.o modtest_crs.o

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

testdrive.o: $(PATHTESST)/testdrive.f90
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $<

%.o: %.f90 $(LIBTEST)
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $<

$(PROGS): %: $(OBJS) $(LIBTEST) $(LIB)  %.o
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(FLIBS)

test_coo.o: modtest_coo.o
modtest_coo.o: modtest_common.o
modtest_crs.o: modtest_common.o

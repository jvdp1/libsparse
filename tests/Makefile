include ../src/Makeinit
PROGS_SRC = test_sparse.f90

LIBTEST = libtestrive.a
OBJTEST = testdrive.o

NAME = sparse

CPPFLAGS += -I../src/lib/ -I.

OBJS = modtest_common.o modtest_coo.o

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
	$(RM) $(PROGS) $(OBJS) $(CLEAN_FILES) *.o* *.mod *.a

$(LIBTEST): $(OBJTEST)
	ar rcs $@ $^

%.o: %.f90 $(LIBTEST)
	$(f90) $(F90FLAGS) $(CPPFLAGS) -c $<

$(PROGS): %: $(OBJS) $(LIBTEST) $(LIB)  %.o
	$(f90) $(F90FLAGS) $(CPPFLAGS) -o $@ $^ $(LIBMETISROOT)$(LIBMETIS)$(LMETIS)

test_coo.o: modtest_coo.o
modtest_coo.o: modtest_common.o

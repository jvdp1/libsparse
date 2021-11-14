include ../../Makeinit

F90FLAGS+=-I$(LIBSPARSEROOT)$(LIBSPARSE)

EXEC=

ifeq ($(DEBUGENABLE),1)
 DEBUG=-check all -traceback -debug extended
else
 DEBUG=
endif

ifeq ($(DPENABLE),0)
 DP=0
else
 DP=1
endif

ifeq ($(METISENABLE),0)
 METIS=0
 LIBMETISROOT=
 LIBMETIS=
 LMETIS=
else
 METIS=1
 F90FLAGS+=-I$(LIBMETISROOT)$(LIBMETIS)
 EXEC+=test5
endif

ifeq ($(SPAINVENABLE),0)
 SPAINV=0
else
 SPAINV=1
 EXEC+=test6
endif

EXEC+=test1 test2 test3 test4 test1_ll test7 test8 test9 test10

all:	$(EXEC)

%: %.f90
	$(f90)  $(F90FLAGS) $(DEBUG) -fpp -D_DP=$(DP) -D_METIS=$(METIS) -D_SPAINV=$(SPAINV)\
	 $^ $(LIBSPARSEROOT)$(LIBSPARSE)$(LSPARSE) $(LIBMETISROOT)$(LIBMETIS)$(LMETIS) -o $@

clean:
	rm -fv *.dat *.o *.optrpt $(EXEC)

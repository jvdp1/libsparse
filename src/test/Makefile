include ../Makeinit

LIBSPARSEROOT=../../../
LIBSPARSE=libsparse/src/lib/
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

LIBMETISROOT=
LIBMETIS=
LMETIS=
ifeq ($(METISENABLE),0)
 METIS=0
else
 METIS=1
 LIBMETISROOT=../../../
 LIBMETIS=metis-5.1.0/build/Linux-x86_64/libmetis/
 LMETIS=libmetis.a
 F90FLAGS+=-I$(LIBMETISROOT)$(LIBMETIS)
 EXEC+=test5
endif

EXEC+=test1 test2 test3 test4 test1_ll

all:	$(EXEC)

test1:	test1.f90 
	$(f90)  $(F90FLAGS) $(DEBUG) -fpp -D_DP=$(DP) $^ $(LIBSPARSEROOT)$(LIBSPARSE)/libsparse.a $(LIBMETISROOT)$(LIBMETIS)$(LMETIS) -o $@

test1_ll:	test1_ll.f90 
	$(f90)  $(F90FLAGS) $(DEBUG) -fpp -D_DP=$(DP) $^ $(LIBSPARSEROOT)$(LIBSPARSE)/libsparse.a $(LIBMETISROOT)$(LIBMETIS)$(LMETIS) -o $@

test2:	test2.f90
	$(f90)  $(F90FLAGS) $(DEBUG) -fpp -D_DP=$(DP) $^ $(LIBSPARSEROOT)$(LIBSPARSE)/libsparse.a $(LIBMETISROOT)$(LIBMETIS)$(LMETIS) -o $@

test3:	test3.f90
	$(f90)  $(F90FLAGS) $(DEBUG) -fpp -D_DP=$(DP) $^ $(LIBSPARSEROOT)$(LIBSPARSE)/libsparse.a $(LIBMETISROOT)$(LIBMETIS)$(LMETIS) -o $@

test4:	test4.f90
	$(f90)  $(F90FLAGS) $(DEBUG) -fpp -D_DP=$(DP) $^ $(LIBSPARSEROOT)$(LIBSPARSE)/libsparse.a $(LIBMETISROOT)$(LIBMETIS)$(LMETIS) -o $@

test5:	test5.f90
	$(f90)  $(F90FLAGS) $(DEBUG) -fpp -D_DP=$(DP) $^ $(LIBSPARSEROOT)$(LIBSPARSE)/libsparse.a $(LIBMETISROOT)$(LIBMETIS)$(LMETIS) -o $@


clean:
	rm -fv *.dat *.o *.optrpt $(EXEC)

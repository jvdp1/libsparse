include ../Makeinit

ifeq ($(DEBUGENABLE),1)
 DEBUG=-check all -traceback -debug extended
else
 DEBUG=
endif

LIBSPARSEROOT=../../../
LIBSPARSE=libsparse/src/lib/

F90FLAGS+=-I$(LIBSPARSEROOT)$(LIBSPARSE)

EXEC=test1 test2 test3 test4 test1_ll

all:	$(EXEC)

test1:	test1.f90 
	$(f90)  $(F90FLAGS) $(DEBUG) $^ $(LIBSPARSEROOT)$(LIBSPARSE)/libsparse.a -o $@

test1_ll:	test1_ll.f90 
	$(f90)  $(F90FLAGS) $(DEBUG) $^ $(LIBSPARSEROOT)$(LIBSPARSE)/libsparse.a -o $@

test2:	test2.f90
	$(f90)  $(F90FLAGS) $(DEBUG) $^ $(LIBSPARSEROOT)$(LIBSPARSE)/libsparse.a -o $@

test3:	test3.f90
	$(f90)  $(F90FLAGS) $(DEBUG) $^ $(LIBSPARSEROOT)$(LIBSPARSE)/libsparse.a -o $@

test4:	test4.f90
	$(f90)  $(F90FLAGS) $(DEBUG) $^ $(LIBSPARSEROOT)$(LIBSPARSE)/libsparse.a -o $@

clean:
	rm -fv *.dat *.o *.optrpt $(EXEC)

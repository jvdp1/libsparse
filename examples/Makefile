FFLAGS+=-I../src

LIB = ../src/libsparse.a

EXEC=

ifeq ($(METISENABLE),0)
else
 EXEC+=test5
endif

ifeq ($(SPAINVENABLE),0)
else
 EXEC+=test6
endif

EXEC+=test1 test2 test3 test4 test1_ll test7 test8 test9 test10

all:	$(EXEC)

%: %.f90
	$(FC)  $(FFLAGS) $(DEBUG) $^ $(LIB) $(FLIBS) -o $@

clean:
	rm -fv *.dat *.o *.optrpt $(EXEC)

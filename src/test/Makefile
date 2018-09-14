include ../Makeinit

ifeq ($(DEBUGENABLE),1)
 DEBUG=-check all -traceback -debug extended
else
 DEBUG=
endif

EXEC=test1 test2

all:	$(EXEC)

test1:	test1.f90
	$(f90)  -I../lib $(optf90) $(DEBUG) $^ ../lib/libsparse.a -o $@

test2:	test2.f90
	$(f90)  -I../lib $(optf90) $(DEBUG) $^ ../lib/libsparse.a -o $@

clean:
	rm -fv *.dat *.o *.optrpt $(EXEC)

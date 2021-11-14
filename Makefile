# Fortran stdlib Makefile

FC = ifort
FFLAGS = -g -check all -traceback -stand f08
#FC = gfortran
#FFLAGS = -g -Wall -std=f2008
FYPPFLAGS=

export FC
export FFLAGS
export FYPPFLAGS

.PHONY: all clean test

all:
	$(MAKE) --directory=src/lib
	$(MAKE) --directory=tests

test:
	$(MAKE) --directory=tests test
	@echo
	@echo "All tests passed."

clean:
	$(MAKE) clean --directory=src/lib
	$(MAKE) clean --directory=tests

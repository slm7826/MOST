# EXECUTABLES = MOST.x phi.x SIM.x
EXECUTABLES = MOST.x PROFILE.x evaluate_RSL_I.x evaluate_RSL_I2.x evaluate_RSL_F.x  lookup2nc.x test_lookup_2D.x
TESTS       = test_integrate.x test_lookup.x lookup2nc.x
GARBAGE=$(EXECUTABLES) $(TESTS)

INCLUDE = -I. -I$(shell nf-config --includedir)
LIBS    = $(shell nf-config --flibs)

# use gfortran by default
ifndef FC
    FC := gfortran
endif

ifeq ($(FC),gfortran)
    FFLAGS = -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -Waliasing $(INCLUDE)
    # DFLAGS = -Wall -Wno-unused-function -fcheck=all -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file -pedantic  -ffpe-trap=invalid,overflow,zero -fbacktrace -finit-real=snan
    DFLAGS = -Wall -Wno-unused-function -fcheck=all -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file -pedantic  -fbacktrace -finit-real=snan
else ifeq ($(FC),ifort)
    FFLAGS = -O2 -fp-model fast -fpp -Wp,-w -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -i4 -r8 -heap-arrays -g -traceback $(INCLUDE)
    DFLAGS =
else
    $(error Not sure what flafs to use for FORTRAN compiler "${FC}". Check your environment variable FC.)
endif

all: $(EXECUTABLES) $(TESTS)

.PHONY: clean realclean
clean:
	rm -fr $(GARBAGE) *.dSYM *.dwarf *.o *.mod *.i90 a.out
realclean: clean
	rm -f input.nml *.out.txt

# dependencies
$(EXECUTABLES) test_lookup.x : integrate.o rsl_functions.o monin_obukhov_kernel.o monin_obukhov_functions.o monin_obukhov_mod.o
monin_obukhov_functions.o : rsl_functions.o integrate.o
monin_obukhov_kernel.o    : monin_obukhov_functions.o
test_integrate.x          : integrate.o
# rules
%.x : %.o
	$(FC) $(FFLAGS) $(DFLAGS) $(LIBS) -o $@ $^

%.o : %.F90
	$(FC) $(FFLAGS) $(DFLAGS) -c -o $@ $<

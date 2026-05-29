EXECUTABLES = MOST.x SIM.x SIM2.x SIMA.x

GARBAGE=$(EXECUTABLES)

FC     = gfortran

ifeq ($(FC),gfortran)
    FFLAGS = -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -Waliasing -I.
    # DFLAGS = -Wall -Wno-unused-function -fcheck=all -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file -pedantic  -ffpe-trap=invalid,overflow,zero -fbacktrace -finit-real=snan
    DFLAGS = -Wall -Wno-unused-function -fcheck=all -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file -pedantic  -fbacktrace -finit-real=snan
else ifeq ($(FC),ifort)
    FFLAGS = -O2 -fp-model fast -fpp -Wp,-w -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -i4 -r8 -heap-arrays -g -traceback -I.
    DFLAGS =
endif


all: $(EXECUTABLES) $(TESTS)

.PHONY: clean realclean
clean:
	rm -fr $(GARBAGE) *.dSYM *.dwarf *.o *.mod *.i90 a.out
realclean: clean
	rm -f input.nml *.out.txt

# dependencies
$(EXECUTABLES)        : monin_obukhov_inter.o

# rules
%.x : %.o
	$(FC) $(FFLAGS) $(DFLAGS) -o $@ $^

%.o : %.F90
	$(FC) $(FFLAGS) $(DFLAGS) -c -o $@ $<

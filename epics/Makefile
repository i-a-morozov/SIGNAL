BINDIR   = .

CC = gcc-10
FC = gfortran-10

PROGRAM  = decomposition

INCDIR += $(EPICS_BASE)/include/os/Linux
INCDIR += $(EPICS_BASE)/include
INCDIR += $(EPICS_BASE)/include/compiler/gcc
INCDIR += .

LIBDIR += $(EPICS_BASE)/lib/$(EPICS_HOST_ARCH)
LIBDIR += .

LIB    += ca Com 
LIB    += m
LIB    += signal lapack blas fftw3 gsl gslcblas gfortran

CFLAG  += -O3
CFLAG  += -std=c99 -Wall -pedantic
CFLAG  += -MMD -MP

CPPFLAG += $(INCDIR:%=-I %)

LDFLAG += $(LIBDIR:%=-L %)
LDFLAG += $(LIBDIR:%=-Wl,-rpath,%)
LDFLAG += $(LIB:%=-l %)

SRC += $(wildcard *.c)

OBJ = $(addprefix build_$(EPICS_HOST_ARCH)/,$(addsuffix .o,$(basename $(SRC))))

all: $(PROGRAM)

.PHONY: clean realclean
clean:
	rm -rf build_*
	rm libsignal.a

realclean: clean
	rm -f $(foreach file,$(PROGRAM),$(BINDIR)/$(file))
	rm libsignal.a

build_%:
	mkdir $@

$(BINDIR):
	mkdir -p $@

build: $(BINDIR)/$(PROGRAM)

$(BINDIR)/$(PROGRAM): $(BINDIR) build_$(EPICS_HOST_ARCH) build_$(EPICS_HOST_ARCH)/$(PROGRAM)
	rm -f $@
	cp build_$(EPICS_HOST_ARCH)/$(@F) $@

build_$(EPICS_HOST_ARCH)/$(PROGRAM): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAG)
	# $(CXX) -o $@ $^ $(LDFLAG)

build_$(EPICS_HOST_ARCH)/%.o: %.c libsignal.a signal.h
	$(CC) -c $(CPPFLAG) $(CFLAG) $< -o $@
	# $(CC) -c $(CPPFLAG) $(CFLAG) $< -o $@

libsignal.a: signal.f90
	$(FC) -c -cpp -Wno-unused-function -std=f2018 -Wall -pedantic -O3 -ffast-math -march=native signal.f90
	ar rcs libsignal.a signal.o
	rm -f signal.o signal.mod

-include build_$(EPICS_HOST_ARCH)/*.d

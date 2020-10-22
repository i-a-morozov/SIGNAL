BINDIR   = .

CC = gcc-10
FC = gfortran-10

PROGRAM  = decomposition

INCDIR += $(EPICS_BASE)/include/os/Linux
INCDIR += $(EPICS_BASE)/include
INCDIR += $(EPICS_BASE)/include/compiler/gcc
INCDIR += ../src

LIBDIR += $(EPICS_BASE)/lib/$(EPICS_HOST_ARCH)
LIBDIR += ../src

LIB    += ca Com
LIB    += signal arpack lapack blas fftw3 gsl gslcblas gfortran
LIB    += m

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
	rm -f decomposition
	rm -f *.mod *.smod

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

build_$(EPICS_HOST_ARCH)/%.o: %.c ../src/libsignal.a ../src/signal.h
	$(CC) -c $(CPPFLAG) $(CFLAG) $< -o $@
	# $(CC) -c $(CPPFLAG) $(CFLAG) $< -o $@

../src/libsignal.a: ../src/*.f90
	make -C ../src

-include build_$(EPICS_HOST_ARCH)/*.d

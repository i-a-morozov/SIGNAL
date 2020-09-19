lib += -lsignal
lib += -llapack -lblas
lib += -lm
lib += -lfftw3
lib += -lgsl -lgslcblas
lib += -lgfortran

fc = gfortran-10
cc = gcc-10

fcf = -cpp -std=f2018 -Wall -pedantic -O3 -ffast-math -march=native
ccf = -std=c99 -Wall -pedantic -O3 -ffast-math -march=native

all: example01
all: example02
all: example03
all: example04
all: example05
all: example06
all: example07
all: example08

example01: libsignal.a signal.mod example01.f90
	$(fc) -o example01 $(fcf) -L. -I. example01.f90 $(lib)

example02: libsignal.a signal.mod example02.f90
	$(fc) -o example02 -fopenmp $(fcf) -L. -I. example02.f90 $(lib)

example03: libsignal.a signal.mod example03.f90
	$(fc) -o example03 $(fcf) -L. -I. example03.f90 $(lib)

example04: libsignal.a signal.mod example04.f90
	$(fc) -o example04 $(fcf) -L. -I. example04.f90 $(lib)

example05: libsignal.a signal.mod example05.f90
	$(fc) -o example05 $(fcf) -L. -I. example05.f90 $(lib)

example06: libsignal.a signal.mod signal.h example06.c
	$(cc) -o example06 $(ccf) -L. -I. example06.c $(lib)

example07: libsignal.a signal.mod example07.f90
	$(fc) -o example07 -fopenmp $(fcf) -L. -I. example07.f90 $(lib)

example08: libsignal.a signal.mod example08.f90
	$(fc) -o example08 $(fcf) -L. -I. example08.f90 $(lib)

libsignal.a: signal.o
	ar rcs libsignal.a signal.o

signal.o: signal.f90
	$(fc) -c $(fcf) -Wno-unused-function signal.f90

clean:
	rm libsignal.a
	rm signal.o
	rm signal.mod
	rm example01
	rm example02
	rm example03
	rm example04
	rm example05
	rm example06
	rm example07
	rm example08

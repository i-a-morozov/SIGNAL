lib += -lsignal
lib += -llapack -lblas -larpack
lib += -lm
lib += -lfftw3
lib += -lgsl -lgslcblas
lib += -lgfortran

fc = gfortran-10
cc = gcc-10

fcf = -cpp -std=f2018 -fPIC -frecursive -Wall -pedantic -O3 -ffast-math -march=native
ccf = -std=c99 -Wall -pedantic -O3 -ffast-math -march=native

src = ../src
mod = ../src/mod

all: example01
all: example02
all: example03
all: example04
all: example05
all: example06
all: example07
all: example08

example01: example01.f90 ../src/libsignal.a
	$(fc) -o example01 $(fcf) -L$(src) -I$(src) -I$(mod) example01.f90 $(lib)

example02: example02.f90 ../src/libsignal.a
	$(fc) -o example02 -fopenmp $(fcf) -L$(src) -I$(src) -I$(mod) example02.f90 $(lib)

example03: example03.f90 ../src/libsignal.a
	$(fc) -o example03 $(fcf) -L$(src) -I$(src) -I$(mod) example03.f90 $(lib)

example04: example04.f90 ../src/libsignal.a
	$(fc) -o example04 $(fcf) -L$(src) -I$(src) -I$(mod) example04.f90 $(lib)

example05: example05.f90 ../src/libsignal.a
	$(fc) -o example05 $(fcf) -L$(src) -I$(src) -I$(mod) example05.f90 $(lib)

example06: example06.c ../src/libsignal.a ../src/signal.h
	$(cc) -o example06 $(ccf) -L$(src) -I$(src) -I$(mod) example06.c $(lib)

example07: example07.f90 ../src/libsignal.a
	$(fc) -o example07 -fopenmp $(fcf) -L$(src) -I$(src) -I$(mod) example07.f90 $(lib)

example08: example08.f90 ../src/libsignal.a
	$(fc) -o example08 $(fcf) -L$(src) -I$(src) -I$(mod) example08.f90 $(lib)

../src/libsignal.a: ../src/*.f90
	make -C ../src

clean:
	rm -f  example01
	rm -f  example02
	rm -f  example03
	rm -f  example04
	rm -f  example05
	rm -f  example06
	rm -f  example07
	rm -f  example08
	rm -f *.mod *.smod

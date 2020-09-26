MOD  := mod
OBJ  := obj

$(shell mkdir -p $(MOD) $(OBJ))

DEP  := \
	$(OBJ)/signal.o         \
	$(OBJ)/auxiliary.o      \
	$(OBJ)/transformation.o \
	$(OBJ)/svd.o            \
	$(OBJ)/process.o        \
	$(OBJ)/peakdetect.o     \
	$(OBJ)/window.o         \
	$(OBJ)/frequency.o      \
	$(OBJ)/decomposition.o  \
	$(OBJ)/optimization.o


FC   :=  gfortran-10

FCF  := -cpp
FCF  += -std=f2018 -Wall -pedantic
FCF  += -fPIC
FCF  += -frecursive
FCF  += -O3 -ffast-math -march=native
FCF  += -Wno-unused-function -Wno-unused-value
FCF  +=  -J $(MOD)



all: libsignal.a

libsignal.a: $(DEP)
	ar rcs libsignal.a $(DEP)

$(OBJ)/signal.o: signal.f90
	$(FC) $(FCF) -c signal.f90 -o $(OBJ)/signal.o

$(OBJ)/auxiliary.o: auxiliary.f90
	$(FC) $(FCF) -c auxiliary.f90 -o $(OBJ)/auxiliary.o

$(OBJ)/transformation.o: transformation.f90
	$(FC) $(FCF) -c transformation.f90 -o $(OBJ)/transformation.o

$(OBJ)/svd.o: svd.f90
	$(FC) $(FCF) -c svd.f90 -o $(OBJ)/svd.o

$(OBJ)/process.o: process.f90
	$(FC) $(FCF) -c process.f90 -o $(OBJ)/process.o

$(OBJ)/peakdetect.o: peakdetect.f90
	$(FC) $(FCF) -c peakdetect.f90 -o $(OBJ)/peakdetect.o

$(OBJ)/window.o: window.f90
	$(FC) $(FCF) -c window.f90 -o $(OBJ)/window.o

$(OBJ)/frequency.o: frequency.f90
	$(FC) $(FCF) -c frequency.f90 -o $(OBJ)/frequency.o

$(OBJ)/decomposition.o: decomposition.f90
	$(FC) $(FCF) -c decomposition.f90 -o $(OBJ)/decomposition.o	

$(OBJ)/optimization.o: optimization.f90
	$(FC) $(FCF) -c optimization.f90 -o $(OBJ)/optimization.o

clean:
	rm -fr $(MOD)
	rm -fr $(OBJ)
	rm -f libsignal.a
	rm -f *.mod *.smod

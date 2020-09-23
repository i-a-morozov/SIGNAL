# SIGNAL, 2018-2020, I.A.MOROZOV@INP.NSK.SU
# FREQUENCY IDENTIFICATION, QUASIPERIODIC DECOMPOSITION AND CHAOS INDICATORS

*Files*:
1) 'SRC_SIGNAL.WL'       WM code (use Get[] function to load or see WM examples)
2) 'signal_examples.nb'  WM examples
3) 'data'                directory with WM examples' data
4) 'signal.f90'          FORTRAN code
5) 'signal.h'            C INTERFACE
6) 'examples'            FORTRAN and C examples
7) 'epics'               EPICS

*Dependencies*:
Basic FORTRAN functions can be used without external libraries.
FFT_EXTERNAL_ uses FFTW3.
SVD related procedures use BLAS, LAPACK and ARPACK.
GSL library will be used for special functions.

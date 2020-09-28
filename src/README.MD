# SIGNAL, 2018-2020, I.A.MOROZOV@INP.NSK.SU
# FREQUENCY IDENTIFICATION, QUASIPERIODIC DECOMPOSITION AND CHAOS INDICATORS


*Files*:
01) signal.h             -- C header
02) signal.inc           -- C style directives
03) signal.f90           -- main module (interfaces and parameters)
04) auxiliary.f90        -- auxiliary procedures
05) transformation.f90   -- transformations procedures (fft, ffrft)
06) svd.f90              -- svd procedures (full, truncated)
07) process.f90          -- data processing procedures (remove mean, filtering)
08) peakdetect.f90       -- peak detection procedures
09) window.f90           -- window generation procedures
10) frequency.f90        -- frequency identification procedures
11) decomposition.f90    -- signal decomposition procedures
12) optimization.f90     -- optimization procedures (LS, maximum search)
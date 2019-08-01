# SIGNAL, 2019, I.A.MOROZOV@INP.NSK.SU
# QUASIPERIODIC DECOMPOSITION AND CHAOS INDICATORS IN WOLFRAM MATHEMATICA

1) 'SRC_SIGNAL'       main code (use Get[] function to load)
2) 'frequency.f90'    FORTAN module for frequency identification (see EXAMPLE-08 in 'examples.nb')
3) 'examples.nb'      SIGNAL examples
4) 'data'             directory with examples' data

Chaos indicators:
1) Frequency Map Analysis
2) Fractal Dimension Estimation
3) Fast Lyapunov Indicator
4) The Smaller and the Generalized Alignment Indices

List of functions:
S$SAMPLING$RATE -- signal sampling rate (integer) 
S$COMPLEX$FLAG -- complex flag (true or false) 
S$RANGE -- frequency range 
S$WINDOW$ORDER -- window order (integer or real) 
S$WINDOW$LENGTH -- window length (integer) 
S$INTERPOLATION$ORDER -- spectra peak interpolation order (integer) 
S$INTERPOLATION$POINTS -- number of points to use for spectra interpolation (integer) 
S$FIT$LENGTH -- fit sample length (integer) with randomly sampled data points 
S$FIT$METHOD -- fit method to be used with LinearModelFit, GeneralizedLinearModelFit or NonlinearModelFit 
S$COSINE$WINDOW[ORD,LEN] -- generate cosine window of order <ORD> (positive integer or real) and window length <LEN> (integer) 
S$WINDOW$DATA -- cosine window generated for given order and number of points (list)
S$SPECTRA[SIG,OPT] -- compute Fourier amplitude specta for given signal <SIG> (list) and options <OPT> 
S$FIND$PEAKS[NUM,DAT] -- compute up to <NUM> (integer) peak positions sorted by values at peaks for given data <DAT> (list) 
S$PAD$ZEROS[LEF,RIG,SIG] -- pad <LEF> (integer) zeros to the left and <RIG> (integer) zeros to the right for given signal <SIG> (list)
S$PROCESS -- signal process function 
S$PRINT$VARIABLES[] -- print global variables 
S$SET$VARIABLES[SAM,FLA,ORD,WIN,INT,POI,FIT,MET] -- initialize SIGNAL with <SAM> (integer) -- sampling rate, <FLA> (true or false) -- complex signal flag, <ORD> (integer or real) -- window order, <WIN> (integer) -- window length, <INT> (integer) -- spectra interpolation order, <POI> (integer) -- number of points to use for spectra interpolation (on each side), <FIT> (integer) -- fitting (random) sample length, <MET> -- fit method 
S$FREQUENCY[SIG] -- estimate main frequency (i.e. with largest Fourier amplitude) for input signal <SIG> (list)
S$FREQUENCY[NUM,SIG] -- estimate given peak frequency for input signal <SIG> (list) and peak number <NUM> (integer) 
S$HARMONICS[ORD,FUN] -- generate list of harmonics up to given order <ORD> (integer) and list of given fundamental frequencies <FUN> (list), returns assosiation {N_1_1,N_1_2} -> H_1
S$IDENTIFY[ORD,FUN,HAR] -- identify list of harmonics <HAR> (list) up to maximum order <ORD> (integer) given list of fundamental frequencies <FUN> (list) and options <OPT>
S$GENERATE$MODEL[FRE,OPT] -- generate signal model to be fitted for given list of known frequencies <FRE> (list) and options <OPT> 
S$FIT$MODEL[FRE,SIG,OPT] -- fit model for given list of frequencies <FRE> (list), input data <SIG> (list) and options <OPT> (only real part of <SIG> is fitted) 
S$QUASIPERIODIC$SUBTRACT[ITR,TYP,SIG] -- compute quasiperiodic decomposition based on subtraction for <ITR> (integer) number of iterations, fitting function <TYP> (function) and signal <SIG> (list of reals) 
S$QUASIPERIODIC$PEAKS[ITR,TYP,SIG] -- compute quasiperiodic decomposition based on peak detection in frequency space for <ITR> (integer) number of peaks, fitting function <TYP> (function) and signal <SIG> (list of reals) 
S$QUASIPERIODIC$HARMONICS$L2[ORD,FUN,TYP,SIG] -- compute quasiperiodic decomposition based on harmonics up to order <ORD> (integer) of fundamental frequencies <FUN> (List), fitting function <TYP> (function) and signal <SIG> (list) 
S$ORTHOGONAL$MATCHING$PURSUIT[MAT,VEC,ERR] -- orthogonal matching pursuit 
S$QUASIPERIODIC$HARMONICS$L1[ORD,BAS,FAC,ERR,MET,SIG] -- compute quasiperiodic decomposition based on harmonics using compressed sensing (FAC ~ 0.8)
S$NAFF[NUM,SIG,LEV] -- implementation of Laskar's NAFF (without orthogonalization), compute <NUM> (integer) terms of quasiperiodic decomposition for given signal <SIG> (real or complex list) with chop level <LEV> (real), S$FREQUENCY[1,<...>] is used for frequency computation 
S$PRONY[SIG,NUM,TOL] -- apply Prony's method for signal reconstruction <SIG> (list) based on <NUM> (integer) singular values and exponent tolerance <TOL> (real)
S$FREQUENCY$LIST[NUM,SIG] -- estimate <NUM> (integer) frequencies for given signal <SIG> (list)
S$DEVIATION[SIG] -- estimate main frequency deviation for given signal <SIG> (list)
S$DEVIATION[NUM,SIG] -- estimate deviation of selected frequency <NUM> (integer) for given signal <SIG> (list) 
S$SHIFT[SHI,SIG] -- main frequency as a function of window shift <SHI> (integer)
S$FMA[SIG] -- perform FMA based on the main frequency for list of signals <SIG> = {<SIG1>,<SIG2>,...}
S$FMA[NUM,SIG] -- perform FMA based on <NUM> (integer) peak for list of signals <SIG> = {<SIG1>,<SIG2>,...} 
S$FMA$SHIFT[SHI,SIGS] -- perform FMA based on the main frequency with window shift <SHI> (integer) 
S$BINARIZE$REGION[BINX,BINY,INTX,INTY,DATA] -- binarization of 3D data <DATA> = {{x1,y1,z1},...} (list) with <BINX> (integer) and <BINY> (integer) number of bins within region defined by <INTX> (list) and <INTY> (list) based on mean value (default), other fuction can be given in options <OPT> 
S$COLOR$PLOT$LEVEL -- color plot base level color 
S$COLOR$PLOT$COLOR -- color plot color scheme 
S$COLOR$PLOT$LIMIT -- color plot max resolution 
S$COLOR$PLOT$PADDING -- color plot image padding 
S$COLOR$PLOT[LEV,INT,HOR,VER,DAT,OPT] -- generate color plot
S$COLOR$LEGEND[LEV,INT,PLO,RAT] -- generate color legend 
S$RESONANCE$LINES[FRE,ORD,HNUM,{HMIN,HMAX},VNUM,{VMIN,VMAX}] -- compute resonance lines for given resonace order <ORD>, integer tune parts <HNUM> and <VNUM> within tune range {{<HMIN>,<HMAX>},{<VMIN>,<VMAX>}} 
S$RESONANCE$PLOT[FRE,ORD,HNUM,{HMIN,HMAX},VNUM,{VMIN,VMAX}] -- plot resonance lines for given resonace order <ORD>, integer tune parts <HNUM> and <VNUM> within tune range {{<HMIN>,<HMAX>},{<VMIN>,<VMAX>}} 
S$SELECT[RES,DEL,DAT,IND] -- select point close (frequency distance <DEL> (real)) to a given resonance <RES> (list of integers) for given data <DAT> with frequencies at positions <IND> (list of integers) 
S$REFINE$INITIALS[NUM,AMP,DAT] -- refine (perturb) initial conditions <DAT> with amplitude <AMP> (real) using duplication number <NUM> (integer) 
S$FILTER$DATA[{I,J},{B_I,{MIN_I,MAX_I}},{B_J,{MIN_J,MAX_J}},DAT] -- filter data <DAT> = {DAT_1,DAT_2,...} = {{X_1,X_2,...,X_I,...,X_J}_1,...} based on {I,J} entries within region {{MIN_I,MAX_I},{MIN_J,MAX_J}} using number of bins {B_I,B_J} 
S$SDDSNAFF[NUM,SIG] -- use sddsnaff to estimate <NUM> frequencies and corresponding parameters for given signal <SIG> (needs sdds tools to be installed and linux environment) 
S$DHT[SIG] -- return discrete Hilbert transform S = U + I H[V] for given signal U = <SIG> (instantaneous frequency W(N) = arg(S(N) conj(S(N-1)))) 
S$PSL[SIG] -- compute (normalized with respect to maximum value) power spectra length 
S$FDE[RAN,SIG,OPT] -- fractal dimension estimation for given list of sample length <RAN> (list of integers), signal <SIG> and options <OPT> 
S$INDEX$LEVEL -- index level for SALI/GALI indicators
S$SALI[{DV1,DV2}] -- compute SALI indicator for given normalized deviation vectors <DV1> and <DV2> 
S$FLI$FUN -- FLI folding function 
S$FLI$ONE[VEC] -- compute FLI indicator using single deviation vector <VEC> given for all iterations, <VEC> = {<VEC>_1, <VEC>_2, ...} 
S$FLI$ALL[VEC] -- compute FLI indicator using several (orthonormal) deviation vectors <VEC> given for all iterations, <VEC> = {<VEC>_1, <VEC>_2, ...}, <VEC>_I = {<VEC>_I_1, <VEC>_I_2, ...} 
S$FLI$LAS[VEC] -- compute FLI indicator using several (orthonormal) deviation vectors <VEC> given for all iterations, <VEC> = {<VEC>_1, <VEC>_2, ...}, <VEC>_I = {<VEC>_I_1, <VEC>_I_2, ...} (returns only last value) 
S$FLI[VEC] -- compute FLI indicator using several (orthonormal) deviation vectors <VEC> given for all iterations, <VEC> = {<VEC>_1, <VEC>_2, ...}, <VEC>_I = {<VEC>_I_1, <VEC>_I_2, ...} (returns only last value) 
S$FILTER$SVD[NUM,SIG] -- apply svd filter using <NUM> (integer) singular value pairs to given signal <SIG> (list) 
S$INFORMATION[] -- list of (inactive) functions and global variables

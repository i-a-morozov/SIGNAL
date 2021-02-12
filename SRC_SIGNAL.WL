(* ################################################################################################################################################################ *)
(* SIGNAL (FREQUENCY IDENTIFICATION, QUASIPERIODIC SIGNAL DECOMPOSITION AND CHAOS INDICATORS), 2018-2020, I.A.MOROZOV@INP.NSK.SU *)
(* ################################################################################################################################################################ *)
ClearSystemCache[] ;                (* -- CLEAR SYSTEM CACHE *)
Set[$HistoryLength,0] ;             (* -- NO HISTORY *)
Get["CompiledFunctionTools`"] ;     (* -- LOAD COMPILED FUNCTION PACKAGE *)
Get["CCompilerDriver`"] ;           (* -- LOAD COMPILER DRIVER PACKAGE *)
Get["MaTeX`"] ;                     (* -- LOAD MATEX PACKAGE *)
(* ################################################################################################################################################################ *)
(* SYSTEM COMPILE OPTIONS  *)
(* ################################################################################################################################################################ *)
Compiler`$CCompilerOptions = List[Rule["SystemCompileOptions","-O3  -ffast-math -march=native"]] ;
(* ################################################################################################################################################################ *)
(* COMPILE FLAG  *)
(* ################################################################################################################################################################ *)
ClearAll[S$COMPILE$FLAG] ;
S$COMPILE$FLAG::usage = "
S$COMPILE$FLAG -- flag to use compiled version of S$FREQUENCY[] function (logical), frequency estimation accuracy corresponds to \"PARABOLA\" method, S$FREQUENCY[] options will be ignored, decomposition functions that pass non-default options to S$FREQUENCY[] might produce wrong results
" ;
S$COMPILE$FLAG = False ;
(* ################################################################################################################################################################ *)
(* SIGNAL SAMPLING RATE *)
(* ################################################################################################################################################################ *)
ClearAll[S$SAMPLING$RATE] ;
S$SAMPLING$RATE::usage = "
S$SAMPLING$RATE -- signal sampling rate (integer)
" ;
S$SAMPLING$RATE = 1 ;
(* ################################################################################################################################################################ *)
(* COMPLEX SIGNAL FLAG *)
(* ################################################################################################################################################################ *)
ClearAll[S$COMPLEX$FLAG] ;
S$COMPLEX$FLAG::usage = "
S$COMPLEX$FLAG -- complex signal flag (logical)
" ;
S$COMPLEX$FLAG = False ;
(* ################################################################################################################################################################ *)
(* FREQUENCY RANGE *)
(* ################################################################################################################################################################ *)
ClearAll[S$FREQUENCY$RANGE] ;
S$FREQUENCY$RANGE::usage = "
S$FREQUENCY$RANGE -- frequency range (real), recomputed by S$SET$VARIABLES[] function
" ;
S$FREQUENCY$RANGE = N[If[S$COMPLEX$FLAG,S$SAMPLING$RATE,Divide[S$SAMPLING$RATE,Plus[1,1]]]] ;
(* ################################################################################################################################################################ *)
(* WINDOW ORDER *)
(* ################################################################################################################################################################ *)
ClearAll[S$WINDOW$ORDER] ;
S$WINDOW$ORDER::usage = "
S$WINDOW$ORDER -- cosine window order (integer)
" ;
S$WINDOW$ORDER = Plus[1,1] ;
(* ################################################################################################################################################################ *)
(* WINDOW LENGTH *)
(* ################################################################################################################################################################ *)
ClearAll[S$WINDOW$LENGTH] ;
S$WINDOW$LENGTH::usage = "
S$WINDOW$LENGTH -- window length (integer)
" ;
S$WINDOW$LENGTH = 2^10 ;
(* ################################################################################################################################################################ *)
(* PEAK INTERPOLATION/FITTING ORDER *)
(* ################################################################################################################################################################ *)
ClearAll[S$PEAK$ORDER] ;
S$PEAK$ORDER::usage = "
S$PEAK$ORDER -- fourier amplitude spectra peak interpolation or fitting order (integer)
" ;
S$PEAK$ORDER = 10 ;
(* ################################################################################################################################################################ *)
(* NUMBER OF INTERPOLATION/FITTING POINTS *)
(* ################################################################################################################################################################ *)
ClearAll[S$PEAK$POINTS] ;
S$PEAK$POINTS::usage = "
S$PEAK$POINTS -- number of points to use for fourier amplitude spectra peak interpolation or fitting (integer)
" ;
S$PEAK$POINTS = 20 ;
(* ################################################################################################################################################################ *)
(* SIGNAL FITTING METHOD *)
(* ################################################################################################################################################################ *)
ClearAll[S$FIT$METHOD] ;
S$FIT$METHOD::usage = "
S$FIT$METHOD -- signal fitting method used with NonlinearModelFit[]
" ;
S$FIT$METHOD = Automatic ;
(* ################################################################################################################################################################ *)
(* COSINE WINDOW FUNCTION *)
(* ################################################################################################################################################################ *)
ClearAll[S$COSINE$WINDOW] ;
S$COSINE$WINDOW::usage = "
S$COSINE$WINDOW[ORDER,LENGTH] -- generate cosine window data of order <ORDER> (integer) and window length <LENGTH> (integer)
" ;
S$COSINE$WINDOW[                    (* -- GENERATE COSINE WINDOW DATA (LIST) *)
  ORDER_?NumericQ,                  (* -- WINDOW ORDER (INTEGER/REAL) *)
  LENGTH_Integer                    (* -- WINDOW LENGTH (INTEGER) *)
] := If[
  SameQ[N[ORDER],N[0]],
  N[ConstantArray[1,LENGTH]],
  Times[
    Power[Plus[1,1],ORDER],
    Power[Factorial[ORDER],Plus[1,1]],
    Power[Factorial[Times[Plus[1,1],ORDER]],Subtract[0,1]],
    Power[Plus[1,Cos[Times[Plus[1,1],Pi,Subtract[Divide[Subtract[N[Range[LENGTH]],1],LENGTH],Divide[1,Plus[1,1]]]]]],ORDER]
  ]
] /; And[
  GreaterEqual[ORDER,0],
  GreaterEqual[LENGTH,1]
] ;
(* ################################################################################################################################################################ *)
(* KAISER WINDOW FUNCTION *)
(* ################################################################################################################################################################ *)
ClearAll[S$KAISER$WINDOW] ;
S$KAISER$WINDOW::usage = "S$KAISER$WINDOW[ORDER,LENGTH] -- generate kaiser window data for order parameter <ORDER> (real) and window length <LENGTH> (integer) " ;
S$KAISER$WINDOW[ORDER_?NumericQ,LENGTH_Integer] := Block[
  {DATA},
  DATA = Subtract[Divide[Subtract[N[Range[LENGTH]],1],LENGTH],Divide[1,Plus[1,1]]] ;
  DATA = BesselI[0,Pi*ORDER*Sqrt[1-2*DATA]*Sqrt[1+2*DATA]]/BesselI[0,Pi*ORDER]
] /; And[
  GreaterEqual[ORDER,0],
  GreaterEqual[LENGTH,1]
] ;
(* ################################################################################################################################################################ *)
(* WINDOW DATA *)
(* ################################################################################################################################################################ *)
ClearAll[S$WINDOW$DATA] ;
S$WINDOW$DATA::usage = "
S$WINDOW$DATA -- cosine window data generated for given order S$WINDOW$ORDER and window length S$WINDOW$LENGTH, recomputed by S$SET$VARIABLES[] function
" ;
S$WINDOW$DATA = S$COSINE$WINDOW[S$WINDOW$ORDER,S$WINDOW$LENGTH] ;
(* ################################################################################################################################################################ *)
(* PRINT GLOBAL VARIABLES *)
(* ################################################################################################################################################################ *)
ClearAll[S$PRINT$GLOBAL] ;
S$PRINT$GLOBAL::usage = "
S$PRINT$GLOBAL[] -- print global variables
" ;
S$PRINT$GLOBAL[] := Column[
  List[
    StringTemplate["S$SAMPLING$RATE         :    `1`"][S$SAMPLING$RATE],
    StringTemplate["S$COMPLEX$FLAG          :    `1`"][S$COMPLEX$FLAG],
    StringTemplate["S$FREQUENCY$RANGE       :    `1`"][S$FREQUENCY$RANGE],
    StringTemplate["S$WINDOW$ORDER          :    `1`"][S$WINDOW$ORDER],
    StringTemplate["S$WINDOW$LENGTH         :    `1`"][S$WINDOW$LENGTH],
    StringTemplate["S$PEAK$ORDER            :    `1`"][S$PEAK$ORDER],
    StringTemplate["S$PEAK$POINTS           :    `1`"][S$PEAK$POINTS],
    StringTemplate["S$FIT$METHOD            :    `1`"][S$FIT$METHOD]
  ],
  Rule[Alignment,Left]
] ;
(* ################################################################################################################################################################ *)
(* SET GLOBAL VARIABLES *)
(* ################################################################################################################################################################ *)
ClearAll[S$SET$VARIABLES] ;
S$SET$VARIABLES::usage = "
S$SET$VARIABLES[RATE,FLAG,ORDER,LENGTH,DEGREE,POINT,METHOD] -- initialize global variables with <RATE> (integer) -- sampling rate, <FLAG> (logical) -- complex signal flag, <ORDER> (integer) -- cosine window order, <LENGTH> (integer) -- window length, <DEGREE> (integer) -- amplitude spectra interpolation/fitting order, <POINT> (integer) -- number of points to use in amplitude spectra interpolation/fitting and <METHOD> -- signal fitting method
" ;
Options[S$SET$VARIABLES] = List[
  Rule["VERBOSE",False]             (* -- VERBOSE FLAG (LOGICAL) *)
] ;
S$SET$VARIABLES[                    (* -- SET GLOBAL VARIABLES *)
  RATE_Integer,                     (* -- SAMPLING RATE (INTEGER) *)
  FLAG_,                            (* -- COMPLEX FLAG (LOGICAL) *)
  ORDER_?NumericQ,                  (* -- WINDOW ORDER (INTEGER) *)
  LENGTH_Integer,                   (* -- WINDOW LENGTH (INTEGER) *)
  DEGREE_Integer,                   (* -- PEAK ORDER (INTEGER) *)
  POINT_Integer,                    (* -- PEAK POINTS (INTEGER) *)
  METHOD_:Automatic,                (* -- SIGNAL FIT METHOD *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Block[
  List[],
  S$SAMPLING$RATE = RATE ;
  S$COMPLEX$FLAG = FLAG ;
  S$FREQUENCY$RANGE = N[If[S$COMPLEX$FLAG,S$SAMPLING$RATE,Divide[S$SAMPLING$RATE,Plus[1,1]]]] ;
  S$WINDOW$ORDER = ORDER ;
  S$WINDOW$LENGTH = LENGTH ;
  S$WINDOW$DATA = S$COSINE$WINDOW[S$WINDOW$ORDER,S$WINDOW$LENGTH] ;
  S$PEAK$ORDER = DEGREE ;
  S$PEAK$POINTS = POINT ;
  S$FIT$METHOD = METHOD ;
  If[OptionValue["VERBOSE"],Print[S$PRINT$GLOBAL[]]] ;
] ;
(* ################################################################################################################################################################ *)
(* FOURIER AMPLITUDE SPECTRA *)
(* ################################################################################################################################################################ *)
ClearAll[S$SPECTRA] ;
S$SPECTRA::usage = "
S$SPECTRA[SIGNAL] -- compute fourier amplitude specta for given signal <SIGNAL> (list)
" ;
Options[S$SPECTRA] = List[
  Rule["FOURIER",List[-1,1]],       (* -- FOURIER PARAMETERS *)
  Rule["NORMALIZED",False],         (* -- FLAG TO NORMALIZE OUTPUT WITH RESPECT TO IT'S MAXIMUM VALUE, ELSE LOG10 OF CORRESPONDING AMPLITUDE SPECTRA IS RETURNED (LOGICAL) *)
  Rule["COMPLEX",False],            (* -- COMPLEX FLAG (LOGICAL) *)
  Rule["2D",False],                 (* -- RETURN 2D DATA (LOGICAL) *)
  Rule["LEVEL",-15.0]               (* -- BASE LEVEL (REAL) *)
] ;
S$SPECTRA[                          (* -- COMPUTE FOURIER AMPLITUDE SPECTA *)
  SIGNAL_List,                      (* -- SIGNAL (LIST) *)
  OPTIONS:OptionsPattern[]          (* -- OPTIONS(S) *)
] := Block[
  {LEVEL,FOURIER},
  LEVEL = 10^OptionValue["LEVEL"] ;
  FOURIER = Plus[Abs[Fourier[SIGNAL,Rule[FourierParameters,OptionValue["FOURIER"]]]],LEVEL] ;
  FOURIER = If[OptionValue["COMPLEX"],FOURIER,Take[FOURIER,List[1,Floor[Divide[Length[FOURIER],Plus[1,1]]]]]] ;
  FOURIER = If[OptionValue["NORMALIZED"],Divide[FOURIER,Max[FOURIER]],Log10[FOURIER]] ;
  If[
    OptionValue["2D"],
    Transpose[List[Times[Plus[1,1],S$FREQUENCY$RANGE,Divide[N[Take[Range[0,Subtract[Length[SIGNAL],1]],Length[FOURIER]]],Length[SIGNAL]]],FOURIER]],
    FOURIER
  ]
] ;
(* ################################################################################################################################################################ *)
(* FIND PEAKS *)
(* ################################################################################################################################################################ *)
ClearAll[S$FIND$PEAKS] ;
S$FIND$PEAKS::usage = "
S$FIND$PEAKS[NUMBER,DATA] -- compute up to <NUMBER> (integer) peak positions sorted by data values at peaks for given data <DATA> (list)
" ;
S$FIND$PEAKS[                       (* -- FIND PEAK POSITIONS (LIST) *)
  NUMBER_Integer,                   (* -- NUMBER OF PEAKS TO FIND (INTEGER) *)
  DATA_List,                        (* -- DATA (LIST) *)
  OPTIONS:OptionsPattern[FindPeaks] (* -- OPTIONS(S) *)
] := Block[
  {POSITION},
  POSITION = Reverse[SortBy[FindPeaks[DATA,0,OPTIONS,Rule[Padding,1]],List[Last]]] ;
  POSITION = Map[First,POSITION] ;
  If[Greater[Length[POSITION],NUMBER],Take[POSITION,List[1,NUMBER]],POSITION]
] /; GreaterEqual[NUMBER,1] ;
(* ################################################################################################################################################################ *)
(* PAD ZEROS *)
(* ################################################################################################################################################################ *)
ClearAll[S$PAD$ZEROS] ;
S$PAD$ZEROS::usage = "
S$PAD$ZEROS[LPAD,RPAD,SIGNAL] -- pad <LPAD> (integer) zeros to the left and <RPAD> (integer) zeros to the right for given signal <SIGNAL> (list)
" ;
S$PAD$ZEROS[                        (* -- PAD ZEROS (LIST) *)
  LPAD_Integer,                     (* -- NUMBER OF LEFT ZEROS (INTEGER) *)
  RPAD_Integer,                     (* -- NUMBER OF RIGHT ZEROS (INTEGER) *)
  SIGNAL_List                       (* -- SIGNAL (LIST) *)
] := Developer`ToPackedArray[Flatten[List[ConstantArray[N[0],LPAD],SIGNAL,ConstantArray[N[0],RPAD]]]] ;
ClearAll[S$ROUND$UP] ;
S$ROUND$UP::usage = "
S$ROUND$UP[NUMBER] -- return next power of two number greater then <NUMBER> (integer)
" ;
S$ROUND$UP[NUMBER_Integer] := 2^Ceiling[Log[N[NUMBER]]/Log[N[2.0]]] ;
ClearAll[S$PAD$LENGTH] ;
S$PAD$LENGTH::usage = "
S$PAD$LENGTH[NUMBER] -- return numbers to pad (sum) to the next power of two
" ;
S$PAD$LENGTH[NUMBER_Integer] := Block[
	{TOTAL},
	TOTAL = Subtract[S$ROUND$UP[NUMBER],NUMBER] ;
	List[Floor[Divide[TOTAL,Plus[1,1]]],Ceiling[Divide[TOTAL,Plus[1,1]]]]
] ;
S$PAD$ZEROS[SIGNAL_List] := Block[
  {LPAD, RPAD},
  {LPAD, RPAD} = S$PAD$LENGTH[Length[SIGNAL]] ;
  S$PAD$ZEROS[LPAD,RPAD,SIGNAL]
] ;
(* ################################################################################################################################################################ *)
(* WINDOW WEIGHTED MEAN *)
(* ################################################################################################################################################################ *)
ClearAll[S$WEIGHTED$MEAN] ;
S$WEIGHTED$MEAN::usage = "
S$WEIGHTED$MEAN[WINDOW,SIGNAL] -- compute window-weighted mean for given window data <WINDOW> (list) and signal <SIGNAL> (list), should have identical length
" ;
S$WEIGHTED$MEAN[                    (* -- WINDOW-WEIGHTED MEAN *)
  WINDOW_List,                      (* -- WINDOW DATA (LIST) *)
  SIGNAL_List                       (* -- INPUT SIGNAL (LIST) *)
] := Divide[
  Total[Times[SIGNAL,WINDOW],Rule[Method,"CompensatedSummation"]],
  Total[WINDOW,Rule[Method,"CompensatedSummation"]]
] /; SameQ[Length[WINDOW],Length[SIGNAL]] ;
(* ################################################################################################################################################################ *)
(* SIGNAL PROCESS FUNCTION (APPLIED BEFORE WINDOW AND/OR ZERO PADDING) *)
(* ################################################################################################################################################################ *)
ClearAll[S$PROCESS] ;
S$PROCESS::usage = "
S$PROCESS -- signal pre-processing (function to apply before window and zero padding) function used in frequency identification functions (by default window-weighted mean is removed)
" ;
S$PROCESS[                          (* -- SIGNAL PRE-PROCESSING FUNCTION *)
  SIGNAL_List                       (* -- INPUT SIGNAL (LIST) *)
] := Subtract[SIGNAL,S$WEIGHTED$MEAN[S$WINDOW$DATA,SIGNAL]] /; SameQ[Length[S$WINDOW$DATA],Length[SIGNAL]] ;
(* ################################################################################################################################################################ *)
(* MAKE MATRIX *)
(* ################################################################################################################################################################ *)
ClearAll[S$MAKE$MATRIX] ;
S$MAKE$MATRIX::usage = "
S$MAKE$MATRIX[SIGNAL] -- generate history matrix for given signal <SIGNAL> (list)
" ;
S$MAKE$MATRIX[                      (* -- MAKE MATRIX *)
  SIGNAL_List                       (* -- INPUT SIGNAL (LIST) *)
] := Block[
  {DATA,LENGTH,MATRIX,MTU,MTS,MTV,SVL},
  DATA = SIGNAL ;
  LENGTH = Length[DATA] ;
  MATRIX = Partition[DATA,Floor[Divide[LENGTH,Plus[1,1]]],1] ;
  MATRIX = Developer`ToPackedArray[MATRIX] ;
  MATRIX
] ;
(* ################################################################################################################################################################ *)
(* MAKE SIGNAL *)
(* ################################################################################################################################################################ *)
ClearAll[S$MAKE$SIGNAL$ROW] ;
S$MAKE$SIGNAL$ROW[MATRIX_] := Flatten[List[First[MATRIX],Last[MATRIX]]] ;
ClearAll[S$MAKE$SIGNAL$DIAGONAL] ;
S$MAKE$SIGNAL$$DIAGONAL[MATRIX_] := Block[
  {NR, NC, ROTATE, DIAGONAL},
  {NR, NC} = Dimensions[MATRIX] ;
  ROTATE = Reverse[Transpose[MATRIX]] ;
  DIAGONAL = Table[Diagonal[ROTATE,COUNT],{COUNT,1-NC,NR-1}] ;
  Map[Mean,DIAGONAL]
] ;
ClearAll[S$MAKE$SIGNAL$SUM] ;
S$MAKE$SIGNAL$SUM[MATRIX_] := Block[
  {ROW, COL, ARRAY, SUM, SHIFT, START, COUNT, Q, P},
  {ROW, COL} = Dimensions[MATRIX] ;
  ARRAY = Flatten[Transpose[MATRIX]]    ;
  SUM = Table[0,2*COL]  ;
  SHIFT = 1 ;
  Do[
  {
    START = Q ;
    If[
      Q > ROW,
      START = Q + COL*SHIFT ;
      SHIFT++ ;
    ] ;
    COUNT = 0 ;
    Do[
    {
      SUM[[Q]] += ARRAY[[START+P*COL]] ;
      COUNT++ ;
    },
    {P,0,Min[Q,COL]-SHIFT,1}
    ] ;
    SUM[[Q]]/=COUNT ;
  },
  {Q,1,2*COL,1}
  ] ;
  SUM
] ;
ClearAll[S$MAKE$SIGNAL] ;
S$MAKE$SIGNAL::usage = "
S$MAKE$SIGNAL[MATRIX] -- generate signal for given history matrix <MATRIX> (matrix)
" ;
S$MAKE$SIGNAL[                      (* -- MAKE SIGNAL *)
  List[MATRIX__List],               (* -- INPUT HISTORY MATRIX (MATRIX) *)
  BUILDER_:S$MAKE$SIGNAL$DIAGONAL   (* -- BUILDER FUNCTION (S$MAKE$SIGNAL$ROW, S$MAKE$SIGNAL$SUM OR S$MAKE$SIGNAL$DIAGONAL (DEFAULT)) *)
] := Block[
  {SIGNAL},
  SIGNAL = BUILDER[List[MATRIX]] ;
  SIGNAL = Developer`ToPackedArray[SIGNAL] ;
  SIGNAL
] ;
(* ################################################################################################################################################################ *)
(* SVD LIST *)
(* ################################################################################################################################################################ *)
ClearAll[S$SVD$LIST] ;
S$SVD$LIST::usage = "
S$SVD$LIST[NUMBER,SIGNAL] -- compute <NUMBER> (integer) singular values for given signal <SIGNAL> (list)
" ;
S$SVD$LIST[                         (* -- SVD LIST *)
  NUMBER_Integer,                   (* -- NUMBER OF SINGULAR VALUES TO COMPUTE (INTEGER) *)
  SIGNAL_List                       (* -- INPUT SIGNAL (LIST) *)
] := Block[
  {DATA,LENGTH,MATRIX,MTU,MTS,MTV,SVL},
  DATA = SIGNAL ;
  LENGTH = Length[DATA] ;
  MATRIX = S$MAKE$MATRIX[SIGNAL] ;
  List[MTU,MTS,MTV] = SingularValueDecomposition[MATRIX,NUMBER] ;
  Diagonal[MTS]
] ;
(* ################################################################################################################################################################ *)
(* SVD FILTER *)
(* ################################################################################################################################################################ *)
ClearAll[S$SVD$FILTER] ;
S$SVD$FILTER::usage = "
S$SVD$FILTER[NUMBER,SIGNAL] -- apply svd filter using <NUMBER> (integer) singular values to given signal <SIGNAL> (list)
" ;
S$SVD$FILTER[                       (* -- SVD FILTER *)
  NUMBER_Integer,                   (* -- NUMBER OF SINGULAR VALUE TO KEEP (INTEGER) *)
  SIGNAL_List,                      (* -- INPUT SIGNAL (LIST) *)
  BUILDER_:S$MAKE$SIGNAL$DIAGONAL   (* -- BUILDER FUNCTION (S$MAKE$SIGNAL$ROW, S$MAKE$SIGNAL$SUM OR S$MAKE$SIGNAL$DIAGONAL (DEFAULT)) *)
] := Block[
  {DATA,LENGTH,MATRIX,MTU,MTS,MTV,SVL},
  DATA = SIGNAL ;
  LENGTH = Length[DATA] ;
  MATRIX = S$MAKE$MATRIX[SIGNAL] ;
  List[MTU,MTS,MTV] = SingularValueDecomposition[MATRIX,NUMBER] ;
  MATRIX = Dot[MTU,MTS,Transpose[MTV]] ;
  S$MAKE$SIGNAL[MATRIX,BUILDER]
] ;
(* ################################################################################################################################################################ *)
(* GENERATE HARMONICS *)
(* ################################################################################################################################################################ *)
ClearAll[S$HARMONICS] ;
S$HARMONICS::usage = "
S$HARMONICS[ORDER,BASIS] -- generate list of harmonics up to given order <ORDER> (integer) for list of given basis frequencies <BASIS> (list)
" ;
S$HARMONICS[                        (* -- GENERATE HARMONICS FOR GIVEN BASIS (ASSOSIATION) *)
  ORDER_Integer,                    (* -- MAXIMUM ORDER (INTEGER) *)
  BASIS_List                        (* -- LIST OF FUNDAMENTAL FREQUENCIES (BASIS) *)
] := Block[
  {HARMONIC,RULE,LIST,TABLE},
  RULE = Rule[Array[HARMONIC,Length[BASIS]],Abs[Mod[Dot[BASIS,Array[HARMONIC,Length[BASIS]]],Times[Plus[1,1],S$FREQUENCY$RANGE],Subtract[0,S$FREQUENCY$RANGE]]]] ;
  LIST = ConstantArray[Range[Subtract[0,ORDER],ORDER],Length[BASIS]] ;
  LIST = Transpose[List[Array[HARMONIC,Length[BASIS]],LIST]] ;
  LIST = MapAt[Select[NonNegative],LIST,List[1,Plus[1,1]]] ;
  TABLE = Flatten[Apply[Table,Join[List[RULE],LIST]]] ;
  TABLE = Select[TABLE,Composition[Curry[LessEqual][ORDER],Total,Abs,First]] ;
  TABLE = Rest[DeleteCases[TABLE,Rule[List[Repeated[0],_?Negative,___],_]]] ;
  KeySortBy[TABLE,Composition[Total,Abs]]
] ;
(* ################################################################################################################################################################ *)
(* IDENTIFY HARMONICS *)
(* ################################################################################################################################################################ *)
ClearAll[S$IDENTIFY] ;
S$IDENTIFY::usage = "
S$IDENTIFY[ORDER,BASIS,LIST] -- identify list of harmonics <LIST> (list) up to maximum order <ORDER> (integer) for given frequency basis <BASIS> (list)
" ;
Options[S$IDENTIFY] = List[
  Rule["TABLE",False]               (* -- FLAG TO RETURN TABLE OF HARMONICS *)
] ;
S$IDENTIFY[                         (* -- IDENTIFY HARMONICS *)
  ORDER_Integer,                    (* -- ORDER LIMIT (INTEGER) *)
  BASIS_List,                       (* -- FUNDAMENTAL FREQUENCIES (LIST) *)
  LIST_List,                        (* -- LIST OF FREQUENCIES TO IDENTIFY (LIST)*)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Block[
  {TABLE,SELECT},
  TABLE = S$HARMONICS[ORDER,BASIS] ;
  SELECT = Map[
    Function[
      Block[
        {HARMONIC,INDEX,ERROR},
        HARMONIC = Slot[1] ;
        List[INDEX,ERROR] = Apply[List,First[Normal[Take[Sort[Abs[Subtract[TABLE,Slot[1]]]],List[1,1]]]]] ;
        List[HARMONIC,Abs[Mod[Dot[INDEX,BASIS],Times[Plus[1,1],S$FREQUENCY$RANGE],Subtract[0,S$FREQUENCY$RANGE]]],ERROR,INDEX]
      ]
    ] ,
    LIST
  ] ;
  If[OptionValue["TABLE"],List[SELECT,TABLE],SELECT]
] ;
(* ################################################################################################################################################################ *)
(* IDENTIFY FUNDAMENTAL BASIS *)
(* ################################################################################################################################################################ *)
ClearAll[S$FUNDAMENTAL] ;
S$FUNDAMENTAL::usage = "
S$FUNDAMENTAL[DIMENSION,ORDER,LIST] -- identify fundamental basis for given dimension <DIMENSION> (integer), maximum order <ORDER> (integer) and list of frequencies <LIST> (list), it is assumed that <DIMENSION> basis frequencies are present in a <LIST>, this function might give wrong result
" ;
S$FUNDAMENTAL[                      (* -- IDENTIFY FUNDAMENTAL BASIS *)
  DIMENSION_Integer,                (* -- BASIS DIMENSION (INTEGER) *)
  ORDER_Integer,                    (* -- MAXIMUM ORDER (INTEGER) *)
  LIST_List                         (* -- LIST OF FREQUENCIES (LIST) *)
] := Block[
  {BASIS,DATA},
  BASIS = Tuples[LIST,DIMENSION] ;
  BASIS = Select[BASIS,Apply[UnsameQ]] ;
  BASIS = DeleteDuplicates[Map[Sort,BASIS]] ;
  DATA = Map[Function[S$IDENTIFY[ORDER,Slot[1],LIST]],BASIS] ;
  DATA = DeleteCases[DATA,List[__,List[RepeatedNull[0],1,RepeatedNull[0]]],Infinity] ;
  DATA = Map[Composition[Total,Map[Composition[Last,Most]]],DATA] ;
  DATA = SortBy[Transpose[List[BASIS,DATA]],Last] ;
  DATA = First[DATA] ;
  First[DATA]
] ;
(* ################################################################################################################################################################ *)
(* AMPLITUDE *)
(* ################################################################################################################################################################ *)
ClearAll[S$AMPLITUDE] ;
S$AMPLITUDE::usage = "
S$AMPLITUDE[SIGNAL][FREQUENCY] -- estimate absolute amplitude for given frequency <FREQUENCY> (real) and input signal <SIGNAL> (list)
" ;
Options[S$AMPLITUDE] = List[
  Rule["ABS",True]                  (* -- RETURN ABSOLUTE VALUE *)
] ;
S$AMPLITUDE[                        (* -- ESTIMATE ABSOLUTE AMPLITYDE (REAL) *)
  SIGNAL_List,                      (* -- INPUT SIGNAL (LIST) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
][
  FREQUENCY_?NumberQ                (* -- FREQUENCY (REAL) *)
] := Block[
  {DELTA, AMPLITUDE},
  DELTA = Exp[Times[Subtract[0,Plus[1,1]],Pi,FREQUENCY,I,Range[S$WINDOW$LENGTH]]] ;
  AMPLITUDE = Divide[
    Total[Times[S$WINDOW$DATA,SIGNAL,Conjugate[DELTA]],Rule[Method,"CompensatedSummation"]],
    Total[S$WINDOW$DATA,Rule[Method,"CompensatedSummation"]]
  ] ;
  If[
    OptionValue["ABS"],
    Log10[Abs[AMPLITUDE]],
    AMPLITUDE
  ]
] ;
(* ################################################################################################################################################################ *)
(* BINARY SEARCH (MAXIMIZATION) *)
(* ################################################################################################################################################################ *)
ClearAll[S$SEARCH$BINARY] ;
S$SEARCH$BINARY::usage = "
S$SEARCH$BINARY[FUNCTION,INITIAL,STEP,LIMIT,TOLERANCE] -- perform binary search maximization
" ;
S$SEARCH$BINARY[                    (* -- BINARY SEARCH MAXIMIZATION *)
  FUNCTION_,                        (* -- FUNCTION TO MAXIMAZE (REAL) -> (REAL) *)
  INITIAL_,                         (* -- INITIAL GUESS (REAL) *)
  STEP_,                            (* -- HALF INTERVAL SIZE (REAL) *)
  LIMIT_,                           (* -- MAXIMUM NUMBER OF ITERATIONS (INTEGER) *)
  TOLERANCE_                        (* -- FREQUENCY STEP TOLERANCE (REAL) *)
] := Block[
  {DELTA, XL, XR, FL, FR, DIRECTION},
  DELTA = STEP ;
  XL = INITIAL - DELTA ;
  XR = INITIAL + DELTA ;
  FL = FUNCTION[XL] ;
  FR = FUNCTION[XR] ;
  Do[
    {
      DIRECTION = Sign[Subtract[FL,FR]] ;
      {XL,XR} = If[
        DIRECTION > 0,
        XL + {-DELTA, +DELTA} ,
        XR + {-DELTA, +DELTA}
      ] ;
      DELTA = DELTA/2;
      FL = FUNCTION[XL] ;
      FR = FUNCTION[XR] ;
      If[
        LessEqual[Abs[Subtract[FL,FR]],TOLERANCE],
        Break[]
      ] ;
    },
    {LIMIT}
  ] ;
  DIRECTION = Sign[Subtract[FL,FR]] ;
  If[DIRECTION > 0, XL, XR]
] ;
(* ################################################################################################################################################################ *)
(* GOLDEN SEARCH (MAXIMIZATION) *)
(* ################################################################################################################################################################ *)
ClearAll[S$SEARCH$GOLDEN] ;
S$SEARCH$GOLDEN::usage = "
S$SEARCH$GOLDEN[FUNCTION,INITIAL,STEP,LIMIT,TOLERANCE] -- perform golden search maximization
" ;
S$SEARCH$GOLDEN[                    (* -- GOLDEN SEARCH MAXIMIZATION *)
  FUNCTION_,                        (* -- FUNCTION TO MAXIMAZE (REAL) -> (REAL) *)
  INITIAL_,                         (* -- INITIAL GUESS (REAL) *)
  STEP_,                            (* -- HALF INTERVAL SIZE (REAL) *)
  LIMIT_,                           (* -- MAXIMUM NUMBER OF ITERATIONS (INTEGER) *)
  TOLERANCE_                        (* -- FREQUENCY STEP TOLERANCE (REAL) *)
] := Block[
  {DELTA, XL, XR, PSI, PHI, M1, M2, FL, FR, LOOP},
  DELTA = STEP ;
  XL = INITIAL - DELTA ;
  XR = INITIAL + DELTA ;
  {PSI, PHI} = N[{1 - 1/GoldenRatio, 1/GoldenRatio}] ;
  M1 = Rescale[PSI,{0,1},{XL,XR}] ;
  M2 = Rescale[PHI,{0,1},{XL,XR}] ;
  FL = FUNCTION[M1] ;
  FR = FUNCTION[M2] ;
  LOOP = 0 ;
  While[
    Abs[Subtract[FL,FR]] > TOLERANCE && LOOP < LIMIT,
    LOOP = LOOP + 1;
    If[
      FL > FR,
      XR = M2 ;
      M2 = M1 ;
      M1 = PSI*XL + PHI*M2;
      FR = FL ;
      FL = FUNCTION[M1],
      XL = M1 ;
      M1 = M2 ;
      M2 = PHI*M1 + PSI*XR;
      FL = FR ;
      FR = FUNCTION[M2] ;
    ] ;
  ] ;
  If[FL > FR, XL, XR]
] ;
(* ################################################################################################################################################################ *)
(* FREQUENCY *)
(* ################################################################################################################################################################ *)
ClearAll[S$FREQUENCY] ;
S$FREQUENCY::usage = "
S$FREQUENCY[SIGNAL] -- estimate frequency for given input signal <SIGNAL> (list)
S$FREQUENCY[SIGNAL,{MIN,MAX}] -- estimate frequency for given input signal <SIGNAL> (list) and frequency range {<MIN>,<MAX>}
S$FREQUENCY[SIGNAL,SEARCH,LIMIT,TOLERANCE] -- estimate frequency for given input signal <SIGNAL> (list), search function <SEARCH>, maximum number of iterations <LIMIT> (integer) and tolerance <TOLERANCE> (real)
" ;
Options[S$FREQUENCY] = List[
  Rule["PROCESS",True],             (* -- SIGNAL PROCESSING FLAG (LOGICAL) *)
  Rule["PADDING",0],                (* -- ZERO PADDING (INTEGER) *)
  Rule["PEAK",0],                   (* -- PEAK NUMBER (INTEGER) *)
  Rule["METHOD","PARABOLA"]         (* -- FREQUENCY ESTIMATION METHOD (STRING) ("FFT", "REFINE", "PARABOLA", "INTERPOLATION" OR "FIT") *)
] ;
S$FREQUENCY[                        (* -- ESTIMATE FREQUENCY (LIST) *)
  SIGNAL_List,                      (* -- INPUT SIGNAL (LIST) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Catch[
  Block[
    {LENGTH,SAMPLE,PROCESS,PADDING,FOURIER,BIN,SHIFT,MAX,LIST,DATA,ROOT,MODEL,COEFFICIENT,VARIABLE,FREQUENCY},
    (* SET INITIAL SAMPLE LENGTH *)
    LENGTH = S$WINDOW$LENGTH ;
    (* SET SAMPLE *)
    SAMPLE = Take[SIGNAL,LENGTH] ;
    (* APPLY PROCESS FUNCTION *)
    PROCESS = OptionValue["PROCESS"] ;
    If[
      PROCESS,
      SAMPLE = S$PROCESS[SAMPLE] ;
    ] ;
    (* COMPILED VERSION FLAG *)
    If[
      S$COMPILE$FLAG,
      SAMPLE = Developer`ToPackedArray[N[Flatten[Transpose[List[Re[SAMPLE],Im[SAMPLE]]]]]] ;
      Throw[S$COMPILED$FREQUENCY[Boole[S$COMPLEX$FLAG],S$WINDOW$LENGTH,Total[S$WINDOW$DATA],S$WINDOW$DATA,SAMPLE]] ;
    ] ;
    (* APPLY WINDOW *)
    SAMPLE = Times[S$WINDOW$DATA,SAMPLE] ;
    (* ZERO PADDING *)
    PADDING = OptionValue["PADDING"] ;
    If[
      UnsameQ[PADDING,0],
      SAMPLE = S$PAD$ZEROS[PADDING,PADDING,SAMPLE] ;
      LENGTH = Length[SAMPLE] ;
    ] ;
    (* COMPUTE FOURIER AMPLITUDE SPECTRA *)
    FOURIER = S$SPECTRA[SAMPLE,Rule["COMPLEX",S$COMPLEX$FLAG],Rule["FOURIER",List[-1,1]]] ;
    (* MAXIMUM OR SELECTED PEAK POSITION *)
    BIN = If[
      SameQ[OptionValue["PEAK"],0],
      Extract[Position[FOURIER,Max[FOURIER]],List[1,1]],
      Last[S$FIND$PEAKS[OptionValue["PEAK"],FOURIER]]
    ] ;
    (* METHOD (FFT) *)
    If[
      SameQ[OptionValue["METHOD"],"FFT"],
      FREQUENCY = N[Times[S$SAMPLING$RATE,Divide[Subtract[BIN,1],LENGTH]]] ;
      Throw[FREQUENCY]
    ] ;
    (* SHIFTED FREQUENCY *)
    SHIFT = N[Times[S$SAMPLING$RATE,Divide[Subtract[BIN,Plus[1,1]],LENGTH]]] ;
    (* MODULATE SIGNAL *)
    SAMPLE = Times[SAMPLE,Exp[Times[Plus[1,1],Pi,I,SHIFT,N[Divide[Range[0,Subtract[LENGTH,1]],S$SAMPLING$RATE]]]]] ;
    (* REFINED SPECTRA *)
    FOURIER = S$SPECTRA[SAMPLE,Rule["COMPLEX",True],Rule["FOURIER",List[0,Divide[Plus[1,1],LENGTH]]]] ;
    (* REFINED SPECTRA MAXIMUM *)
    MAX = Extract[Position[FOURIER,Max[FOURIER]],List[1,1]] ;
    (* METHOD (REFINE) *)
    If[
      SameQ[OptionValue["METHOD"],"REFINE"],
      FREQUENCY = N[Times[S$SAMPLING$RATE,Power[LENGTH,Subtract[0,1]],Plus[Times[Plus[1,1],Subtract[0,1]],BIN,Times[Plus[1,1],Power[LENGTH,Subtract[0,1]],Plus[Subtract[0,1],MAX]]]]] ;
      Throw[FREQUENCY]
    ] ;
    (* METHOD (PARABOLA) *)
    If[
      SameQ[OptionValue["METHOD"],"PARABOLA"],
      If[
        Or[SameQ[MAX,1],SameQ[MAX,Length[FOURIER]]],
        FREQUENCY = N[Times[S$SAMPLING$RATE,Power[LENGTH,Subtract[0,1]],Plus[Times[Plus[1,1],Subtract[0,1]],BIN,Times[Plus[1,1],Power[LENGTH,Subtract[0,1]],Plus[Subtract[0,1],MAX]]]]] ;
        Throw[FREQUENCY]
      ] ;
      MAX = MAX-1/2+(FOURIER[[MAX-1]]-FOURIER[[MAX]])/(FOURIER[[MAX-1]]-2*FOURIER[[MAX]]+FOURIER[[MAX+1]]) ;
      FREQUENCY = N[Times[S$SAMPLING$RATE,Power[LENGTH,Subtract[0,1]],Plus[Times[Plus[1,1],Subtract[0,1]],BIN,Times[Plus[1,1],Power[LENGTH,Subtract[0,1]],Plus[Subtract[0,1],MAX]]]]] ;
      Throw[FREQUENCY]
    ] ;
    (* INTERPOLATION/FIT POINTS *)
    LIST = S$PEAK$POINTS ;
    LIST = Flatten[List[Reverse[Subtract[MAX,Range[LIST]]],MAX,Plus[MAX,Range[LIST]]]] ;
    LIST = DeleteCases[LIST,PatternTest[Blank[],NonPositive]] ;
    LIST = TakeWhile[LIST,Function[Less[Slot[1],LENGTH]]] ;
    (* DATA TO USE IN INTERPOLATION/FITTING *)
    DATA = Transpose[List[N[LIST],Part[FOURIER,LIST]]] ;
    (* METHOD (INTERPOLATION) *)
    If[
      SameQ[OptionValue["METHOD"],"INTERPOLATION"],
      DATA = Interpolation[DATA,Rule[InterpolationOrder,S$PEAK$ORDER],Rule[Method,"Spline"]] ;
      MAX = Extract[Quiet[FindRoot[D[DATA[ROOT],ROOT],List[ROOT,MAX],Rule[Method,"Newton"],Rule[MaxIterations,Infinity],Rule[AccuracyGoal,$MachinePrecision],Rule[PrecisionGoal,$MachinePrecision]]],List[1,Plus[1,1]]] ;
      FREQUENCY = N[Times[S$SAMPLING$RATE,Power[LENGTH,Subtract[0,1]],Plus[Times[Plus[1,1],Subtract[0,1]],BIN,Times[Plus[1,1],Power[LENGTH,Subtract[0,1]],Plus[Subtract[0,1],MAX]]]]] ;
      Throw[FREQUENCY]
    ] ;
    (* METHOD (FIT) *)
    If[
      SameQ[OptionValue["METHOD"],"FIT"],
      MODEL = Dot[
        Array[COEFFICIENT,S$PEAK$ORDER],
        Power[VARIABLE,Range[S$PEAK$ORDER]]
      ] ;
      MODEL = Quiet[NonlinearModelFit[DATA,MODEL,Array[COEFFICIENT,S$PEAK$ORDER],VARIABLE,Rule[MaxIterations,Infinity],Rule[AccuracyGoal,$MachinePrecision],Rule[PrecisionGoal,$MachinePrecision]]] ;
      MODEL = MODEL["Function"] ;
      MAX = Extract[Quiet[FindRoot[D[MODEL[ROOT],ROOT],List[ROOT,MAX],Rule[Method,"Newton"],Rule[MaxIterations,Infinity],Rule[AccuracyGoal,$MachinePrecision],Rule[PrecisionGoal,$MachinePrecision]]],List[1,Plus[1,1]]] ;
      FREQUENCY = N[Times[S$SAMPLING$RATE,Power[LENGTH,Subtract[0,1]],Plus[Times[Plus[1,1],Subtract[0,1]],BIN,Times[Plus[1,1],Power[LENGTH,Subtract[0,1]],Plus[Subtract[0,1],MAX]]]]] ;
      Throw[FREQUENCY]
    ] ;
  ]
] ;
S$FREQUENCY[                        (* -- ESTIMATE FREQUENCY (LIST) *)
  SIGNAL_List,                      (* -- INPUT SIGNAL (LIST) *)
  List[FMIN_,FMAX_],                (* -- FREQUENCY RANGE *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Catch[
  Block[
    {LENGTH,SAMPLE,PROCESS,PADDING,FOURIER,BMIN,BMAX,BIN,SHIFT,MAX,LIST,DATA,ROOT,MODEL,COEFFICIENT,VARIABLE,FREQUENCY},
    (* SET INITIAL SAMPLE LENGTH *)
    LENGTH = S$WINDOW$LENGTH ;
    (* SET SAMPLE *)
    SAMPLE = Take[SIGNAL,LENGTH] ;
    (* APPLY PROCESS FUNCTION *)
    PROCESS = OptionValue["PROCESS"] ;
    If[
      PROCESS,
      SAMPLE = S$PROCESS[SAMPLE] ;
    ] ;
    (* COMPILED VERSION FLAG *)
    If[
      S$COMPILE$FLAG,
      SAMPLE = Developer`ToPackedArray[N[Flatten[Transpose[List[Re[SAMPLE],Im[SAMPLE]]]]]] ;
      Throw[S$COMPILED$FREQUENCY[Boole[S$COMPLEX$FLAG],S$WINDOW$LENGTH,Total[S$WINDOW$DATA],S$WINDOW$DATA,SAMPLE]] ;
    ] ;
    (* APPLY WINDOW *)
    SAMPLE = Times[S$WINDOW$DATA,SAMPLE] ;
    (* ZERO PADDING *)
    PADDING = OptionValue["PADDING"] ;
    If[
      UnsameQ[PADDING,0],
      SAMPLE = S$PAD$ZEROS[PADDING,PADDING,SAMPLE] ;
      LENGTH = Length[SAMPLE] ;
    ] ;
    (* COMPUTE FOURIER AMPLITUDE SPECTRA *)
    FOURIER = S$SPECTRA[SAMPLE,Rule["COMPLEX",S$COMPLEX$FLAG],Rule["FOURIER",List[-1,1]]] ;
    BMIN = Floor[FMIN*LENGTH] + 1 ;
    BMAX = Floor[FMAX*LENGTH] + 0 ;
    FOURIER = Take[FOURIER,List[BMIN,BMAX]] ;
    (* MAXIMUM OR SELECTED PEAK POSITION *)
    BIN = If[
      SameQ[OptionValue["PEAK"],0],
      Extract[Position[FOURIER,Max[FOURIER]],List[1,1]],
      Last[S$FIND$PEAKS[OptionValue["PEAK"],FOURIER]]
    ] ;
    BIN = BMIN-1+BIN ;
    (* METHOD (FFT) *)
    If[
      SameQ[OptionValue["METHOD"],"FFT"],
      FREQUENCY = N[Times[S$SAMPLING$RATE,Divide[Subtract[BIN,1],LENGTH]]] ;
      Throw[FREQUENCY]
    ] ;
    (* SHIFTED FREQUENCY *)
    SHIFT = N[Times[S$SAMPLING$RATE,Divide[Subtract[BIN,Plus[1,1]],LENGTH]]] ;
    (* MODULATE SIGNAL *)
    SAMPLE = Times[SAMPLE,Exp[Times[Plus[1,1],Pi,I,SHIFT,N[Divide[Range[0,Subtract[LENGTH,1]],S$SAMPLING$RATE]]]]] ;
    (* REFINED SPECTRA *)
    FOURIER = S$SPECTRA[SAMPLE,Rule["COMPLEX",True],Rule["FOURIER",List[0,Divide[Plus[1,1],LENGTH]]]] ;
    (* REFINED SPECTRA MAXIMUM *)
    MAX = Extract[Position[FOURIER,Max[FOURIER]],List[1,1]] ;
    (* METHOD (REFINE) *)
    If[
      SameQ[OptionValue["METHOD"],"REFINE"],
      FREQUENCY = N[Times[S$SAMPLING$RATE,Power[LENGTH,Subtract[0,1]],Plus[Times[Plus[1,1],Subtract[0,1]],BIN,Times[Plus[1,1],Power[LENGTH,Subtract[0,1]],Plus[Subtract[0,1],MAX]]]]] ;
      Throw[FREQUENCY]
    ] ;
    (* METHOD (PARABOLA) *)
    If[
      SameQ[OptionValue["METHOD"],"PARABOLA"],
      If[
        Or[SameQ[MAX,1],SameQ[MAX,Length[FOURIER]]],
        FREQUENCY = N[Times[S$SAMPLING$RATE,Power[LENGTH,Subtract[0,1]],Plus[Times[Plus[1,1],Subtract[0,1]],BIN,Times[Plus[1,1],Power[LENGTH,Subtract[0,1]],Plus[Subtract[0,1],MAX]]]]] ;
        Throw[FREQUENCY]
      ] ;
      MAX = MAX-1/2+(FOURIER[[MAX-1]]-FOURIER[[MAX]])/(FOURIER[[MAX-1]]-2*FOURIER[[MAX]]+FOURIER[[MAX+1]]) ;
      FREQUENCY = N[Times[S$SAMPLING$RATE,Power[LENGTH,Subtract[0,1]],Plus[Times[Plus[1,1],Subtract[0,1]],BIN,Times[Plus[1,1],Power[LENGTH,Subtract[0,1]],Plus[Subtract[0,1],MAX]]]]] ;
      Throw[FREQUENCY]
    ] ;
    (* INTERPOLATION/FIT POINTS *)
    LIST = S$PEAK$POINTS ;
    LIST = Flatten[List[Reverse[Subtract[MAX,Range[LIST]]],MAX,Plus[MAX,Range[LIST]]]] ;
    LIST = DeleteCases[LIST,PatternTest[Blank[],NonPositive]] ;
    LIST = TakeWhile[LIST,Function[Less[Slot[1],LENGTH]]] ;
    (* DATA TO USE IN INTERPOLATION/FITTING *)
    DATA = Transpose[List[N[LIST],Part[FOURIER,LIST]]] ;
    (* METHOD (INTERPOLATION) *)
    If[
      SameQ[OptionValue["METHOD"],"INTERPOLATION"],
      DATA = Interpolation[DATA,Rule[InterpolationOrder,S$PEAK$ORDER],Rule[Method,"Spline"]] ;
      MAX = Extract[Quiet[FindRoot[D[DATA[ROOT],ROOT],List[ROOT,MAX],Rule[Method,"Newton"],Rule[MaxIterations,Infinity],Rule[AccuracyGoal,$MachinePrecision],Rule[PrecisionGoal,$MachinePrecision]]],List[1,Plus[1,1]]] ;
      FREQUENCY = N[Times[S$SAMPLING$RATE,Power[LENGTH,Subtract[0,1]],Plus[Times[Plus[1,1],Subtract[0,1]],BIN,Times[Plus[1,1],Power[LENGTH,Subtract[0,1]],Plus[Subtract[0,1],MAX]]]]] ;
      Throw[FREQUENCY]
    ] ;
    (* METHOD (FIT) *)
    If[
      SameQ[OptionValue["METHOD"],"FIT"],
      MODEL = Dot[
        Array[COEFFICIENT,S$PEAK$ORDER],
        Power[VARIABLE,Range[S$PEAK$ORDER]]
      ] ;
      MODEL = Quiet[NonlinearModelFit[DATA,MODEL,Array[COEFFICIENT,S$PEAK$ORDER],VARIABLE,Rule[MaxIterations,Infinity],Rule[AccuracyGoal,$MachinePrecision],Rule[PrecisionGoal,$MachinePrecision]]] ;
      MODEL = MODEL["Function"] ;
      MAX = Extract[Quiet[FindRoot[D[MODEL[ROOT],ROOT],List[ROOT,MAX],Rule[Method,"Newton"],Rule[MaxIterations,Infinity],Rule[AccuracyGoal,$MachinePrecision],Rule[PrecisionGoal,$MachinePrecision]]],List[1,Plus[1,1]]] ;
      FREQUENCY = N[Times[S$SAMPLING$RATE,Power[LENGTH,Subtract[0,1]],Plus[Times[Plus[1,1],Subtract[0,1]],BIN,Times[Plus[1,1],Power[LENGTH,Subtract[0,1]],Plus[Subtract[0,1],MAX]]]]] ;
      Throw[FREQUENCY]
    ] ;
  ]
] ;
S$FREQUENCY[                        (* -- ESTIMATE FREQUENCY (DTFT MAXIMUM SEARCH) (LIST) *)
  SIGNAL_List,                      (* -- INPUT SIGNAL (LIST) *)
  SEARCH_,                          (* -- SEARCH FUNCTION (S$SEARCH$BINARY[] OR S$SEARCH$GOLDEN[]) *)
  LIMIT_Integer,                    (* -- MAXIMUM NUMBER OF ITERATIONS (INTEGER) *)
  TOLERANCE_,                       (* -- TOLERANCE (REAL) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Catch[
  Block[
    {LENGTH,SAMPLE,PROCESS,PADDING,FOURIER,BIN},
    LENGTH = S$WINDOW$LENGTH ;
    SAMPLE = Take[SIGNAL,LENGTH] ;
    PROCESS = OptionValue["PROCESS"] ;
    If[
      PROCESS,
      SAMPLE = S$PROCESS[SAMPLE] ;
    ] ;
    SAMPLE = Times[S$WINDOW$DATA,SAMPLE] ;
    PADDING = OptionValue["PADDING"] ;
    If[
      UnsameQ[PADDING,0],
      SAMPLE = S$PAD$ZEROS[PADDING,PADDING,SAMPLE] ;
      LENGTH = Length[SAMPLE] ;
    ] ;
    FOURIER = S$SPECTRA[SAMPLE,Rule["COMPLEX",S$COMPLEX$FLAG],Rule["FOURIER",List[-1,1]]] ;
    BIN = If[
      SameQ[OptionValue["PEAK"],0],
      Extract[Position[FOURIER,Max[FOURIER]],List[1,1]],
      Last[S$FIND$PEAKS[OptionValue["PEAK"],FOURIER]]
    ] ;
    SEARCH[
      S$AMPLITUDE[SAMPLE],
      N[Times[S$SAMPLING$RATE,Divide[Subtract[BIN,1],LENGTH]]],
      N[1/S$WINDOW$LENGTH],
      LIMIT,
      TOLERANCE
    ]
  ]
] ;
(* ################################################################################################################################################################ *)
(* DECOMPOSITION (FOURIER) *)
(* ################################################################################################################################################################ *)
Clear[S$DECOMPOSITION] ;
S$DECOMPOSITION::usage = "
S$DECOMPOSITION[NUMBER,SIGNAL] -- decompose signal <SIGNAL> (list) into <NUMBER> (integer) number of components sorted by amplitude
" ;
Options[S$DECOMPOSITION] = List[
  Rule["TYPE","PEAK"],              (* -- FREQUENCY ESTIMATION TYPE (STRING) ("MAXIMUM" OR "PEAK") *)
  Rule["CHOP",10.^-10],             (* -- AMPLITUDE CHOP LEVEL (REAL) *)
  Rule["PADDING",0],                (* -- ZERO PADDING (INTEGER) *)
  Rule["METHOD","PARABOLA"]         (* -- FREQUENCY ESTIMATION METHOD (STRING) ("FFT", "REFINE", "PARABOLA", "INTERPOLATION" OR "FIT") *)
] ;
S$DECOMPOSITION[                    (* -- DECOMPOSITION (FOURIER) *)
  NUMBER_Integer,                   (* -- NUMBER OF COMPONENTS (INTEGER) *)
  SIGNAL_List,                      (* -- SIGNAL (LIST) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Block[
  {MEAN,SAMPLE,DATA,FREQUENCY,DELTA,AMPLITUDE,LOOP,TIME,SIN,COS,FUNCTION},
  SAMPLE = Take[SIGNAL,S$WINDOW$LENGTH] ;
  MEAN = S$WEIGHTED$MEAN[S$WINDOW$DATA,SAMPLE] ;
  MEAN = Chop[MEAN,OptionValue["CHOP"]] ;
  SAMPLE = Subtract[SAMPLE,MEAN] ;
  DATA = Which[
    SameQ[OptionValue["TYPE"],"MAXIMUM"],
    Reap[
      Do[
        List[
          FREQUENCY = S$FREQUENCY[
            SAMPLE,
            Rule["PROCESS",False],
            Rule["PADDING",OptionValue["PADDING"]],
            Rule["METHOD",OptionValue["METHOD"]],
            Rule["PEAK",0]
          ] ;
          DELTA = Exp[Times[Subtract[0,Plus[1,1]],Pi,FREQUENCY,I,Range[S$WINDOW$LENGTH]]] ;
          AMPLITUDE = Divide[
            Total[Times[S$WINDOW$DATA,SAMPLE,Conjugate[DELTA]],Rule[Method,"CompensatedSummation"]],
            Total[S$WINDOW$DATA,Rule[Method,"CompensatedSummation"]]
          ] ;
          SAMPLE = Subtract[SAMPLE,Times[AMPLITUDE,DELTA]] ;
          Sow[List[FREQUENCY,AMPLITUDE]]
        ],
        List[LOOP,1,NUMBER]
      ] ;
    ],
    SameQ[OptionValue["TYPE"],"PEAK"],
    Reap[
      Do[
        List[
          FREQUENCY = S$FREQUENCY[
            Take[SIGNAL,S$WINDOW$LENGTH],
            Rule["PROCESS",True],
            Rule["PADDING",OptionValue["PADDING"]],
            Rule["METHOD",OptionValue["METHOD"]],
            Rule["PEAK",LOOP]
          ] ;
          DELTA = Exp[Times[Subtract[0,Plus[1,1]],Pi,FREQUENCY,I,Range[S$WINDOW$LENGTH]]] ;
          AMPLITUDE = Divide[
            Total[Times[S$WINDOW$DATA,SAMPLE,Conjugate[DELTA]],Rule[Method,"CompensatedSummation"]],
            Total[S$WINDOW$DATA,Rule[Method,"CompensatedSummation"]]
          ] ;
          SAMPLE = Subtract[SAMPLE,Times[AMPLITUDE,DELTA]] ;
          Sow[List[FREQUENCY,AMPLITUDE]]
        ],
        List[LOOP,1,NUMBER]
      ] ;
    ]
  ] ;
  DATA = First[Last[DATA]] ;
  FREQUENCY = Map[First,DATA] ;
  AMPLITUDE = Map[Last,DATA] ;
  AMPLITUDE = Map[Curry[Times][N[If[S$COMPLEX$FLAG,1,Plus[1,1]]]],AMPLITUDE] ;
  AMPLITUDE = Chop[AMPLITUDE,OptionValue["CHOP"]] ;
  COS = N[Re[AMPLITUDE]] ;
  SIN = N[Im[AMPLITUDE]] ;
  FUNCTION = If[
    S$COMPLEX$FLAG,
    (COS+I*SIN)*Cos[2*Pi*FREQUENCY*TIME]+(SIN-I*COS)*Sin[2*Pi*FREQUENCY*TIME],
    COS*Cos[2*Pi*FREQUENCY*TIME]+SIN*Sin[2*Pi*FREQUENCY*TIME]
  ] ;
  FUNCTION = Apply[Function,List[TIME,Plus[MEAN,Total[FUNCTION]]]] ;
  Association[
    Rule["MEAN",MEAN],
    Rule["DATA",Reverse[SortBy[Transpose[List[FREQUENCY,COS,SIN]],Composition[Norm,Rest]]]],
    Rule["FUNCTION",FUNCTION]
  ]
] /; GreaterEqual[NUMBER,1] ;
(* ################################################################################################################################################################ *)
(* MODEL (GENERATE) *)
(* ################################################################################################################################################################ *)
ClearAll[S$GENERATE$MODEL] ;
S$GENERATE$MODEL::usage = "
S$GENERATE$MODEL[LIST] -- generate signal model to be fitted for given list of known frequencies <LIST> (list)
" ;
Options[S$GENERATE$MODEL] = List[
  Rule["DAMPING",False]           (* -- DAMPING EXPONENTS FLAG (LOGICAL) *)
] ;
S$GENERATE$MODEL[                 (* -- GENERATE SIGNAL MODEL *)
  LIST_List,                      (* -- LIST OF FREQUENCIES (LIST)*)
  OPTIONS:OptionsPattern[]        (* -- OPTION(S) *)
] := Block[
  {LENGTH,COS,SIN,EXP,MEAN,TIME,MODEL,VARIABLE,DATA},
  LENGTH = Length[LIST] ;
  (* SET MODEL PARAMETERS *)
  COS = ToExpression[Array[StringTemplate["COS$`1`"],LENGTH]] ;
  SIN = ToExpression[Array[StringTemplate["SIN$`1`"],LENGTH]] ;
  EXP = ToExpression[Array[StringTemplate["EXP$`1`"],LENGTH]] ;
  MEAN = ToExpression["MEAN$"] ;
  MODEL = Map[Total,Transpose[List[Times[COS,Cos[Times[Plus[1,1],Pi,TIME,LIST]]],Times[SIN,Sin[Times[Plus[1,1],Pi,TIME,LIST]]]]]] ;
  MODEL = Flatten[List[MEAN,MODEL]] ;
  VARIABLE = List[MEAN,COS,SIN] ;
  DATA = List[MODEL,Flatten[VARIABLE],TIME] ;
  DATA = If[
    OptionValue["DAMPING"],
    List[First[MapAt[Composition[Total,Function[Times[Slot[1],Flatten[List[1,Exp[Times[EXP,TIME,Subtract[0,1]]]]]]]],DATA,1]],Flatten[List[VARIABLE,EXP]],TIME],
    MapAt[Total,DATA,1]
  ] ;
  VARIABLE = If[OptionValue["DAMPING"],Flatten[List[VARIABLE,EXP,TIME]],Flatten[List[VARIABLE,TIME]]] ;
  Apply[Function,List[VARIABLE,DATA]]
] ;
(* ################################################################################################################################################################ *)
(* MODEL (FIT) *)
(* ################################################################################################################################################################ *)
ClearAll[S$FIT$MODEL] ;
S$FIT$MODEL::usage = "
S$FIT$MODEL[LIST,SIGNAL] -- fit signal for given list of frequencies <LIST> (list) and input signal <SIGNAL> (list)
" ;
Options[S$FIT$MODEL] = List[
  Rule["DAMPING",False],            (* -- DAMPING EXPONENTS FLAG (LOGICAL) *)
  Rule["NONLINEAR",False],          (* -- INCLUDE FREQUENCIES INTO FITTING (LOGICAL) *)
  Rule["DECOMPOSITION",False],      (* -- USE DECOMPOSITION (LOGICAL) *)
  Rule["TYPE","PEAK"],              (* -- FREQUENCY ESTIMATION TYPE (STRING) ("MAXIMUM" OR "PEAK") *)
  Rule["CHOP",10.^-10],             (* -- AMPLITUDE CHOP LEVEL (REAL) *)
  Rule["PADDING",0],                (* -- ZERO PADDING (INTEGER) *)
  Rule["METHOD","PARABOLA"]         (* -- FREQUENCY ESTIMATION METHOD (STRING) ("FFT", "REFINE", "PARABOLA", "INTERPOLATION" OR "FIT") *)
] ;
S$FIT$MODEL[                        (* -- FIT MODEL FOR GIVEN FREQUENCY LIST *)
  LIST_List,                        (* -- LIST OF FREQUENCIES (LIST) *)
  SIGNAL_List,                      (* -- SIGNAL (LIST)*)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Block[
  {RULE,ERROR,FREQUENCY,MODEL,KNOB,INITIAL,DATA,MEAN,SIN,COS,RESULT},
  RULE = List[Rule[MaxIterations,Infinity],Rule[AccuracyGoal,Infinity],Rule[Method,S$FIT$METHOD]] ;
  ERROR = List["RSquared","AIC","BIC"] ;
  FREQUENCY = Abs[Mod[LIST,Times[Plus[1,1],S$FREQUENCY$RANGE],Subtract[0,S$FREQUENCY$RANGE]]] ;
  MODEL = S$GENERATE$MODEL[FREQUENCY,Rule["DAMPING",OptionValue["DAMPING"]]] ;
  KNOB = First[MODEL] ;
  INITIAL = ConstantArray[0,Length[KNOB]] ;
  If[
    OptionValue["DECOMPOSITION"],
    DATA = S$DECOMPOSITION[
      Length[FREQUENCY],
      SIGNAL,
      Rule["TYPE",OptionValue["TYPE"]],
      Rule["CHOP",OptionValue["CHOP"]],
      Rule["PADDING",OptionValue["PADDING"]],
      Rule["METHOD",OptionValue["METHOD"]]
    ] ;
    MEAN = DATA["MEAN"] ;
    List[COS,SIN] = Transpose[Map[Rest,DATA["DATA"]]] ;
    INITIAL = PadRight[Flatten[List[MEAN,COS,SIN]],Length[INITIAL]] ;
  ] ;
  MODEL = List[First[Apply[MODEL,Plus[KNOB,INITIAL]]],Most[KNOB],Last[KNOB]] ;
  MODEL = Apply[Sequence,MODEL] ;
  RESULT = Quiet[Apply[NonlinearModelFit,List[SIGNAL,MODEL,Apply[Sequence,RULE]]]] ;
  DATA = ReplaceAll[Most[Plus[KNOB,INITIAL]],RESULT["BestFitParameters"]] ;
  MEAN = First[DATA] ;
  DATA = Rest[DATA] ;
  DATA = Take[DATA,Times[Plus[1,1],Length[FREQUENCY]]] ;
  DATA = Partition[DATA,Length[FREQUENCY]] ;
  List[COS,SIN] = DATA ;
  COS = Re[COS] ;
  SIN = Re[SIN] ;
  If[
    OptionValue["NONLINEAR"],
    INITIAL = PadRight[Most[ReplaceAll[KNOB,RESULT["BestFitParameters"]]],Length[KNOB]] ;
    FREQUENCY = ToExpression[Array[StringTemplate["F$`1`"],Length[LIST]]] ;
    MODEL = S$GENERATE$MODEL[Plus[LIST,FREQUENCY],Rule["DAMPING",OptionValue["DAMPING"]]] ;
    MODEL = List[First[Apply[MODEL,Plus[KNOB,INITIAL]]],Flatten[List[Most[KNOB],FREQUENCY]],Last[KNOB]] ;
    MODEL = Apply[Sequence,MODEL] ;
    RESULT = Quiet[Apply[NonlinearModelFit,List[SIGNAL,MODEL,Apply[Sequence,RULE]]]] ;
    DATA = ReplaceAll[Most[Plus[KNOB,INITIAL]],RESULT["BestFitParameters"]] ;
    MEAN = First[DATA] ;
    DATA = Rest[DATA] ;
    DATA = Take[DATA,Times[Plus[1,1],Length[FREQUENCY]]] ;
    DATA = Partition[DATA,Length[FREQUENCY]] ;
    List[COS,SIN] = DATA ;
    FREQUENCY = Plus[LIST,ReplaceAll[FREQUENCY,RESULT["BestFitParameters"]]] ;
    FREQUENCY = Abs[Mod[FREQUENCY,Times[Plus[1,1],S$FREQUENCY$RANGE],Subtract[0,S$FREQUENCY$RANGE]]] ;
  ] ;
  Association[
    Rule["ERROR",Transpose[List[ERROR,Quiet[RESULT[ERROR]]]]],
    Rule["MEAN",MEAN],
    Rule["DATA",Reverse[SortBy[Transpose[List[FREQUENCY,COS,SIN]],Composition[Norm,Rest]]]],
    Rule["MODEL",RESULT],
    Rule["FUNCTION",RESULT["Function"]]
  ]
] ;
(* ################################################################################################################################################################ *)
(*  QUASIPERIODIC DECOMPOSITION (SUBTRACT) *)
(* ################################################################################################################################################################ *)
ClearAll[S$QUASIPERIODIC$SUBTRACT] ;
S$QUASIPERIODIC$SUBTRACT::usage = "
S$QUASIPERIODIC$SUBTRACT[NUMBER,SIGNAL] -- decomposition based on subtraction for <NUMBER> (integer) number of iterations and signal <SIGNAL> (list)
" ;
Options[S$QUASIPERIODIC$SUBTRACT] = List[
  Rule["DAMPING",False],            (* -- DAMPING EXPONENTS FLAG (LOGICAL) *)
  Rule["NONLINEAR",False],          (* -- INCLUDE FREQUENCIES INTO FITTING (LOGICAL), CAN IMPROVE FREQUENCY ESTIMATION IN SOME CASES, COMPARE STATISTICAL INDICATORS WITH AND WITHOUT FLAG *)
  Rule["DECOMPOSITION",False],      (* -- USE DECOMPOSITION (LOGICAL) *)
  Rule["TYPE","PEAK"],              (* -- FREQUENCY ESTIMATION TYPE (STRING) ("MAXIMUM" OR "PEAK") *)
  Rule["CHOP",10.^-10],             (* -- AMPLITUDE CHOP LEVEL (REAL) *)
  Rule["PADDING",0],                (* -- ZERO PADDING (INTEGER) *)
  Rule["METHOD","PARABOLA"]         (* -- FREQUENCY ESTIMATION METHOD (STRING) ("FFT", "REFINE", "PARABOLA", "INTERPOLATION" OR "FIT") *)
] ;
S$QUASIPERIODIC$SUBTRACT[           (* -- QUASIPERIODIC DECOMPOSITION BASED ON ITERATIVE SUBTRACTION *)
  NUMBER_Integer,                   (* -- # OF ITERATIONS (INTEGER) *)
  SIGNAL_List,                      (* -- SIGNAL (LIST) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Block[
  {SAMPLE,FIT,LIST,FREQUENCY,MODEL},
  SAMPLE = Take[SIGNAL,List[1,S$WINDOW$LENGTH]] ;
  FIT = SAMPLE ;
  LIST = List[] ;
  Do[
    List[
      FREQUENCY = S$FREQUENCY[
        SAMPLE,
        Rule["PROCESS",True],
        Rule["PADDING",OptionValue["PADDING"]],
        Rule["PEAK",0],
        Rule["METHOD",OptionValue["METHOD"]]
      ] ;
      LIST = Join[LIST,List[FREQUENCY]] ;
      MODEL = S$FIT$MODEL[
        LIST,
        FIT,
        Rule["DAMPING",OptionValue["DAMPING"]],
        Rule["NONLINEAR",OptionValue["NONLINEAR"]],
        Rule["DECOMPOSITION",OptionValue["DECOMPOSITION"]],
        Rule["TYPE",OptionValue["TYPE"]],
        Rule["CHOP",OptionValue["CHOP"]],
        Rule["PADDING",OptionValue["PADDING"]],
        Rule["METHOD",OptionValue["METHOD"]]
      ] ;
      SAMPLE = Subtract[Take[SIGNAL,List[1,S$WINDOW$LENGTH]],MODEL["MODEL"]["Function"][Divide[Range[S$WINDOW$LENGTH],S$SAMPLING$RATE]]] ;
    ],
    NUMBER
  ] ;
  MODEL
] ;
(* ################################################################################################################################################################ *)
(*  QUASIPERIODIC DECOMPOSITION (PEAK) *)
(* ################################################################################################################################################################ *)
ClearAll[S$QUASIPERIODIC$PEAK] ;
S$QUASIPERIODIC$PEAK::usage = "
S$QUASIPERIODIC$PEAK[NUMBER,SIGNAL] -- decomposition based on <NUMBER> (integer) number of peaks and signal <SIGNAL> (list)
" ;
Options[S$QUASIPERIODIC$PEAK] = List[
  Rule["DAMPING",False],            (* -- DAMPING EXPONENTS FLAG (LOGICAL) *)
  Rule["NONLINEAR",False],          (* -- INCLUDE FREQUENCIES INTO FITTING (LOGICAL), CAN IMPROVE FREQUENCY ESTIMATION IN SOME CASES, COMPARE STATISTICAL INDICATORS WITH AND WITHOUT FLAG *)
  Rule["DECOMPOSITION",False],      (* -- USE DECOMPOSITION (LOGICAL) *)
  Rule["TYPE","PEAK"],              (* -- FREQUENCY ESTIMATION TYPE (STRING) ("MAXIMUM" OR "PEAK") *)
  Rule["CHOP",10.^-10],             (* -- AMPLITUDE CHOP LEVEL (REAL) *)
  Rule["PADDING",0],                (* -- ZERO PADDING (INTEGER) *)
  Rule["METHOD","PARABOLA"]         (* -- FREQUENCY ESTIMATION METHOD (STRING) ("FFT", "REFINE", "PARABOLA", "INTERPOLATION" OR "FIT") *)
] ;
S$QUASIPERIODIC$PEAK[               (* -- QUASIPERIODIC DECOMPOSITION BASED ON ITERATIVE SUBTRACTION *)
  NUMBER_Integer,                   (* -- # OF ITERATIONS (INTEGER) *)
  SIGNAL_List,                      (* -- SIGNAL (LIST) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Block[
  {SAMPLE,FIT,LIST,LOOP,FREQUENCY,MODEL},
  SAMPLE = Take[SIGNAL,List[1,S$WINDOW$LENGTH]] ;
  FIT = SAMPLE ;
  LIST = List[] ;
  Do[
    List[
      FREQUENCY = S$FREQUENCY[
        SAMPLE,
        Rule["PROCESS",True],
        Rule["PADDING",OptionValue["PADDING"]],
        Rule["PEAK",LOOP],
        Rule["METHOD",OptionValue["METHOD"]]
      ] ;
      LIST = Join[LIST,List[FREQUENCY]] ;
    ],
    List[LOOP,1,NUMBER]
  ] ;
  MODEL = S$FIT$MODEL[
    LIST,
    FIT,
    Rule["DAMPING",OptionValue["DAMPING"]],
    Rule["NONLINEAR",OptionValue["NONLINEAR"]],
    Rule["DECOMPOSITION",OptionValue["DECOMPOSITION"]],
    Rule["TYPE",OptionValue["TYPE"]],
    Rule["CHOP",OptionValue["CHOP"]],
    Rule["PADDING",OptionValue["PADDING"]],
    Rule["METHOD",OptionValue["METHOD"]]
  ] ;
  MODEL
] ;
(* ################################################################################################################################################################ *)
(*  QUASIPERIODIC DECOMPOSITION (HARMONICS) *)
(* ################################################################################################################################################################ *)
ClearAll[S$QUASIPERIODIC$HARMONICS] ;
S$QUASIPERIODIC$HARMONICS::usage = "
S$QUASIPERIODIC$HARMONICS[ORDER,BASE,SIGNAL] -- decomposition based on harmonics up to order <ORDER> (integer) of fundamental frequencies <BASIS> (List) and signal <SIGNAL> (list)
" ;
Options[S$QUASIPERIODIC$HARMONICS] = List[
  Rule["DAMPING",False],            (* -- DAMPING EXPONENTS FLAG (LOGICAL) *)
  Rule["NONLINEAR",False],          (* -- INCLUDE FREQUENCIES INTO FITTING (LOGICAL), CAN IMPROVE FREQUENCY ESTIMATION IN SOME CASES, COMPARE STATISTICAL INDICATORS WITH AND WITHOUT FLAG *)
  Rule["DECOMPOSITION",False],      (* -- USE DECOMPOSITION (LOGICAL) *)
  Rule["TYPE","PEAK"],              (* -- FREQUENCY ESTIMATION TYPE (STRING) ("MAXIMUM" OR "PEAK") *)
  Rule["CHOP",10.^-10],             (* -- AMPLITUDE CHOP LEVEL (REAL) *)
  Rule["PADDING",0],                (* -- ZERO PADDING (INTEGER) *)
  Rule["METHOD","PARABOLA"]         (* -- FREQUENCY ESTIMATION METHOD (STRING) ("FFT", "REFINE", "PARABOLA", "INTERPOLATION" OR "FIT") *)
] ;
S$QUASIPERIODIC$HARMONICS[          (* -- QUASIPERIODIC DECOMPOSITION BASED ON HARMONICS OF FUNDAMENTAL FREQUENCIES *)
  ORDER_Integer,                    (* -- MAXIMUM ORDER OF HARMONICS (INTEGER) *)
  BASIS_List,                       (* -- LIST OF FUNDAMENTAL FREQUENCIES (LIST) *)
  SIGNAL_List,                      (* -- SIGNAL (LIST) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Block[
  {LIST,MODEL},
  LIST = Values[S$HARMONICS[ORDER,BASIS]] ;
  MODEL = S$FIT$MODEL[
    LIST,
    SIGNAL,
    Rule["DAMPING",OptionValue["DAMPING"]],
    Rule["NONLINEAR",OptionValue["NONLINEAR"]],
    Rule["DECOMPOSITION",OptionValue["DECOMPOSITION"]],
    Rule["TYPE",OptionValue["TYPE"]],
    Rule["CHOP",OptionValue["CHOP"]],
    Rule["PADDING",OptionValue["PADDING"]],
    Rule["METHOD",OptionValue["METHOD"]]
  ] ;
  MODEL
] ;
(* ################################################################################################################################################################ *)
(* ORTHOGONAL MATCHING PURSUIT *)
(* ################################################################################################################################################################ *)
ClearAll[S$ORTHOGONAL$MATCHING$PURSUIT] ;
S$ORTHOGONAL$MATCHING$PURSUIT::usage = "S$ORTHOGONAL$MATCHING$PURSUIT[MAT,VEC,ERR] -- orthogonal matching pursuit " ;
S$ORTHOGONAL$MATCHING$PURSUIT[      (* -- ORTHOGONAL MATCHING PURSUIT *)
  MATRIX_,                          (* -- MATRIX *)
  VECTOR_,                          (* -- VECTOR *)
  ERROR_,                           (* -- ERROR *)
  MAX_:1024                         (* -- MAXIMUM NUMBER OF ITERATIONS *)
] := Block[
  {LIST,RESULT},
  LIST[] = List[] ;
  LIST[ARGUMENT_] := LIST[] = Flatten[Sort[Join[LIST[],First[Position[ARGUMENT,Max[ARGUMENT]]]]]] ;
  RESULT := Function[
    Block[
      {LOCAL},
      LOCAL = Part[MATRIX,Span[1,All],Slot[1]] ;
      Subtract[VECTOR,Dot[Dot[LOCAL,PseudoInverse[LOCAL]],VECTOR]]
    ]
  ] ;
  NestWhile[Composition[RESULT,LIST,Function[Abs[Dot[Transpose[MATRIX],Slot[1]]]]],VECTOR,Function[Greater[Norm[Flatten[Slot[1]]],ERROR]],1,MAX] ;
  N[Chop[Normal[SparseArray[Thread[Function[Rule[LIST[],Dot[PseudoInverse[Part[MATRIX,Span[1,All],LIST[]]],VECTOR]]][LIST[]]],Last[Dimensions[MATRIX]]]],ERROR]]
] ;
(* ################################################################################################################################################################ *)
(* LINEAR PROGRAMMING (L1-MINIMIZATION) *)
(* ################################################################################################################################################################ *)
ClearAll[S$LINEAR$PROGRAMMING] ;
S$LINEAR$PROGRAMMING::usage = "S$LINEAR$PROGRAMMING[MAT,VEC] -- solve l1 minimization problem for given real matrix <MAT> and real vector <VEC> " ;
S$LINEAR$PROGRAMMING[               (* -- L1 MINIMIZATION WITH LINEAR PROGRAMMING *)
  MATRIX_,                          (* -- MATRIX (REAL) *)
  VECTOR_                           (* -- VECTOR (REAL) *)
] := Block[
  {M, N, TRANSPOSE, SPARSE, A, B, C, BOUND, SOLUTION},
  {M, N} = Dimensions[MATRIX] ;
  TRANSPOSE = Transpose[MATRIX] ;
  SPARSE = SparseArray[{{X_, X_} -> 1}, {M, M}] ;
  A = Join[Transpose[Join[SPARSE, -TRANSPOSE]], Transpose[Join[SPARSE, TRANSPOSE]]] ;
  B = Join[-VECTOR, VECTOR] ;
  C = Join[Table[1, {M}], Table[0, {N}]] ;
  BOUND = Table[-Infinity, {M + N}] ;
  SOLUTION = LinearProgramming[C, A, B, BOUND] ;
  Drop[SOLUTION, M]
] ;
(* ################################################################################################################################################################ *)
(*  QUASIPERIODIC DECOMPOSITION (HARMONICS/COMPRESSED) *)
(* ################################################################################################################################################################ *)
ClearAll[S$QUASIPERIODIC$COMPRESSED]
S$QUASIPERIODIC$COMPRESSED::usage = "
S$QUASIPERIODIC$COMPRESSED[ORDER,BASIS,FACTOR,ERROR,SIGNAL] -- compute quasiperiodic decomposition based on harmonics using compressed sensing
" ;
S$QUASIPERIODIC$COMPRESSED[         (* -- FIND COEFFICIENTS WITH COMPRESSED SENSING *)
  ORDER_Integer,                    (* -- MAXIMUM ORDER (INTEGER) *)
  BASIS_List,                       (* -- BASIS FREQUENCIES (LIST) *)
  FACTOR_Real,                      (* -- RATIO FACTOR (REAL IN [0,1], CHOOSE VALUE CLOSE TO 1 FOR BETTER PERFORMANCE, E.G. 0.8) *)
  ERROR_Real,                       (* -- ERROR (REAL) *)
  SIGNAL_List                       (* -- SIGNAL (LIST) *)
] := Block[
  {LIST,MODEL,TIME,COEFFICIENT,ROW,LENGTH,POSITION,MATRIX,VECTOR,SOLUTION,MEAN,FREQUENCY,COS,SIN},
  LIST = Values[S$HARMONICS[ORDER,BASIS]] ;
  MODEL = S$GENERATE$MODEL[LIST] ;
  COEFFICIENT = Most[First[MODEL]] ;
  MODEL = Last[MODEL] ;
  ROW[TIME_] := Evaluate[N[D[First[MODEL],List[COEFFICIENT]]]] ;
  LENGTH = Floor[Times[Length[COEFFICIENT],FACTOR]] ;
  POSITION = Sort[RandomSample[Range[1,Length[SIGNAL]],LENGTH]] ;
  MATRIX = Map[ROW,Divide[POSITION,S$SAMPLING$RATE]] ;
  VECTOR = Part[SIGNAL,POSITION] ;
  SOLUTION = S$ORTHOGONAL$MATCHING$PURSUIT[MATRIX,VECTOR,ERROR] ;
  MEAN = First[SOLUTION] ;
  MODEL = ReplaceAll[First[MODEL],Thread[Rule[COEFFICIENT,SOLUTION]]] ;
  MODEL = Apply[Function,List[TIME,MODEL]] ;
  SOLUTION = Thread[Rule[LIST,Transpose[Partition[Rest[SOLUTION],Divide[Length[Rest[SOLUTION]],Plus[1,1]]]]]] ;
  SOLUTION = DeleteCases[SOLUTION,Rule[_,List[0.,0.]]] ;
  SOLUTION = Reverse[SortBy[SOLUTION,Function[Norm[Last[Slot[1]]]]]] ;
  FREQUENCY = Map[First,SOLUTION] ;
  List[COS,SIN] = Transpose[Map[Last,SOLUTION]] ;
  COS = Re[COS] ;
  SIN = Re[SIN] ;
  Association[
    Rule["ERROR",Undefined],
    Rule["MEAN",MEAN],
    Rule["DATA",Reverse[SortBy[Transpose[List[FREQUENCY,COS,SIN]],Composition[Norm,Rest]]]],
    Rule["MODEL",MODEL]
  ]
] ;
(* ################################################################################################################################################################ *)
(* DECOMPOSITION (PRONY) *)
(* ################################################################################################################################################################ *)
ClearAll[S$PRONY] ;
S$PRONY::usage = "
S$PRONY[NUMBER,SIGNAL] -- decompose signal <SIGNAL> (list) using <NUMBER> (integer) number of singular values
" ;
Options[S$PRONY] = List[
  Rule["CHOP",10.0^-10]             (* -- AMPLITUDE CHOP LEVEL (REAL) *)
] ;
S$PRONY[                            (* -- DECOMPOSITION (PRONY) *)
  NUMBER_Integer,                   (* -- NUMBER OF SINGULAR VALUES (INTEGER) *)
  SIGNAL_List,                      (* -- INPUT SIGNAL (LIST) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Quiet[
  Block[
    {DATA,SIZE,MATRIX,LIST,RANGE,VECTOR,SOLUTION,ROOT,VARIABLE,EXPONENT,AMPLITUDE,FUNCTION,TIME,MEAN,FREQUENCY,COS,SIN},
    DATA = SIGNAL ;
    SIZE = Length[DATA] ;
    MATRIX = Most[Partition[DATA,Floor[Divide[SIZE,Plus[1,1]]],1]] ;
    MATRIX = Developer`ToPackedArray[MATRIX] ;
    LIST = SingularValueList[MATRIX,NUMBER,Rule[Tolerance,0]] ;
    RANGE = Length[LIST] ;
    MATRIX = Most[Partition[DATA,RANGE,1]] ;
    VECTOR = Drop[DATA,RANGE] ;
    SOLUTION = Dot[PseudoInverse[MATRIX],VECTOR] ;
    ROOT = ReplaceAll[VARIABLE,NSolve[Equal[Subtract[Power[VARIABLE,RANGE],Dot[SOLUTION,Power[VARIABLE,Range[0,Subtract[RANGE,1]]]]],0],VARIABLE]] ;
    EXPONENT = Chop[Divide[Log[ROOT],Times[Plus[1,1],Pi]],OptionValue["CHOP"]] ;
    EXPONENT = Times[S$SAMPLING$RATE,EXPONENT] ;
    MATRIX = Map[Function[Power[ROOT,Slot[1]]],Range[1,SIZE]] ;
    AMPLITUDE = Dot[PseudoInverse[MATRIX],DATA] ;
    FUNCTION = TrigExpand[ExpToTrig[Dot[AMPLITUDE,Exp[Times[Plus[1,1],Pi,EXPONENT,TIME]]]]] ;
    FUNCTION = ReplaceAll[FUNCTION,List[RuleDelayed[Cosh[TIME_],TrigToExp[Cosh[TIME]]],RuleDelayed[Sinh[TIME_],TrigToExp[Sinh[TIME]]]]] ;
    FUNCTION = Collect[FUNCTION,Exp[__]] ;
    FUNCTION = Apply[Function,List[TIME,FUNCTION]] ;
    MEAN = ReplaceAll[FUNCTION[TIME],List[Rule[Cos[__],0],Rule[Sin[__],0]]]  /. TIME -> 0 ;
    FREQUENCY = Cases[Expand[FUNCTION[TIME]],RuleDelayed[Times[AMPLITUDE_,Alternatives[Cos,Sin][FREQUENCY_*TIME]],FREQUENCY]] ;
    FREQUENCY = DeleteDuplicates[Divide[FREQUENCY,Times[Plus[1,1],Pi]]] ;
    COS = N[Re[Coefficient[Expand[FUNCTION[TIME]],Cos[2*Pi*FREQUENCY*TIME]]]] /. TIME -> 0 ;
    SIN = N[Re[Coefficient[Expand[FUNCTION[TIME]],Sin[2*Pi*FREQUENCY*TIME]]]] /. TIME -> 0 ;
    FREQUENCY = Re[FREQUENCY] ;
    Association[
      Rule["ERROR",Undefined],
      Rule["MEAN",MEAN],
      Rule["DATA",Reverse[SortBy[Transpose[List[FREQUENCY,COS,SIN]],Composition[Norm,Rest]]]],
      Rule["MODEL",FUNCTION],
      Rule["FUNCTION",FUNCTION]
    ]
  ]
] ;
(* ################################################################################################################################################################ *)
(* SEARCH FREQUNCIES (SDDSNAFF) *)
(* ################################################################################################################################################################ *)
ClearAll[S$SDDSNAFF] ;
S$SDDSNAFF::usage = "
S$SDDSNAFF[NUMBER,SIGNAL] -- use sddsnaff to estimate <NUMBER> frequencies and corresponding parameters for given signal <SIGNAL>
" ;
S$SDDSNAFF[                         (* -- SEARCH FREQUENCIES WITH SDDSNAFF *)
  NUMBER_Integer,                   (* -- NUMBER OF FREQUENCIES TO SEARCH *)
  SIGNAL_List                       (* -- INPUT SIGNAL *)
] := Block[
  {RUN,DATA},
  RUN = StringJoin[
    "csv2sdds \\\n",
    "  signal.csv -pipe=out \\\n",
    "  -columnData=name=s,type=long \\\n",
    "  -columnData=name=x,type=double |\\\n",
    "sddsnaff \\\n",
    "  -pipe \\\n",
    "  -columns=s,x \\\n",
    StringTemplate["  -terminateSearch=maxFrequencies=`1` \\\n"][NUMBER],
    "  -iterateFrequency=cycleLimit=1000,accuracyLimit=0.000000000000001 |\\\n",
    "sddsconvert \\\n",
    "  -pipe=in signal.sdds \\\n",
    "  -ascii"
  ] ;
  Export["signal.csv",Transpose[List[Range[Length[SIGNAL]],SIGNAL]]] ;
  Run[RUN] ;
  DATA = Select[
    Take[Import["signal.sdds","Table"],List[10,-1]],
    Function[UnsameQ[Slot[1],List[-1.,-1.,-1.,-1.]]]
  ] ;
  Run["rm signal.csv"] ;
  Run["rm signal.sdds"] ;
  DATA
] ;
(* ################################################################################################################################################################ *)
(* FMA *)
(* ################################################################################################################################################################ *)
ClearAll[S$FMA] ;
S$FMA::usage = "
S$FMA[MATRIX] -- estimate mean frequency deviations for given matrix of signals <MATRIX> = {...,<SIGNAL_I>,...} (matrix)
S$FMA[SIGNAL] -- estimate mean frequency deviation for given signal <SIGNAL> (list)
" ;
Options[S$FMA] = List[
  Rule["TYPE","PEAK"],              (* -- FREQUENCY ESTIMATION TYPE (STRING) ("MAXIMUM" OR "PEAK") *)
  Rule["CHOP",10.^-10],             (* -- AMPLITUDE CHOP LEVEL (REAL) *)
  Rule["PADDING",0],                (* -- ZERO PADDING (INTEGER) *)
  Rule["METHOD","PARABOLA"],        (* -- FREQUENCY ESTIMATION METHOD (STRING) ("FFT", "REFINE", "PARABOLA", "INTERPOLATION" OR "FIT") *)
  Rule["SCAN",0],                   (* -- WINDOW SCAN STEP (INTEGER) *)
  Rule["MODE","SINGLE"],            (* -- MODE (STRING) ("SINGLE" OR "MULTIPLE") *)
  Rule["PEAK",Plus[1,1]]            (* -- NUMBER OF PEAKS USED IN "MULTIPLE" MODE *)
] ;
S$FMA[                              (* -- FMA (LIST) *)
  SIGNAL_List,                      (* -- SIGNAL (LIST) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Catch[
  Block[
    {LIST,DATA,FREQUENCY,DEVIATION},
    LIST = List[] ;
    DATA = If[
      SameQ[OptionValue["SCAN"],0],
      TakeList[SIGNAL,Times[S$WINDOW$LENGTH,List[1,Subtract[0,1]]]],
      Partition[SIGNAL,S$WINDOW$LENGTH,OptionValue["SCAN"]]
    ] ;
    DATA = Which[
      SameQ[OptionValue["MODE"],"SINGLE"],
      Map[
        Composition[
          List,
          Function[
            S$FREQUENCY[
              Slot[1],
              Rule["PROCESS",True],
              Rule["PADDING",OptionValue["PADDING"]],
              Rule["PEAK",0],
              Rule["METHOD",OptionValue["METHOD"]]
            ]
          ]
        ],
        DATA
      ],
      SameQ[OptionValue["MODE"],"MULTIPLE"],
      Map[
        Composition[
          Map[First],
          Curry[Apply][List["DATA"]],
          Function[
            S$DECOMPOSITION[
              OptionValue["PEAK"],
              Slot[1],
              Rule["TYPE",OptionValue["TYPE"]],
              Rule["CHOP",OptionValue["CHOP"]],
              Rule["PADDING",OptionValue["PADDING"]],
              Rule["METHOD",OptionValue["METHOD"]]
            ]
          ]
        ],
        DATA
      ]
    ] ;
    Abs[Mod[DATA,Times[Plus[1,1],S$FREQUENCY$RANGE],Subtract[0,S$FREQUENCY$RANGE]]]
  ]
] ;
S$FMA[                              (* -- FMA (MATRIX) *)
  List[MATRIX__List],               (* -- SIGNAL (LIST) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Catch[
  Block[
     {DATA},
     DATA = Map[Function[S$FMA[Slot[1],OPTIONS]],List[MATRIX]] ;
     DATA
   ]
] ;
(* ################################################################################################################################################################ *)
(* 3D DATA BINARIZATION *)
(* ################################################################################################################################################################ *)
(* RELATION BETWEEN SPARCE ARRAY INDEXES {I_S,J_S} AND LAYER INDEXES {I_L,J_L} FOR ARRAY WITH DIMENSIONS (2*N+1)x(2*N+1) *)
(* I_S - I_L = N + 1, I_S = J_M *)
(* J_S + J_L = N + 1, J_S = I_M *)
(* {I_M,J_M} DEFINE WM INDEXED ARRAY {{X_1_1,X_1_2,...},{X_2_1,X_2_2,...}} *)
(* ################################################################################################################################################################ *)
ClearAll[S$BINARIZE$REGION] ;
S$BINARIZE$REGION::usage = "
S$BINARIZE$REGION[{HBIN,HMIN,HMAX},{VBIN,VMIN,VMAX},DATA] -- binarize 3d data <DATA> = {...,{X_I,Y_I,Z_I},...} (list) based on 3rd dimension using <HBIN> (integer) and <VBIN> (integer) number of bins within rectangular region defined by {HMIN (real),HMAX (real)} and {VMIN (real),VMAX (real)}, returns sparse array with WM indexing
S$BINARIZE$REGION[{BIN,HDEL,VDEL},DATA] -- binarize 3d data <DATA> = {...,{X_I,Y_I,Z_I},...} (list) based on 3rd dimension using <BIN> (integer) number of bins within square region defined by horizontal bin size <HDEL> (real) and verticle step size <VDEL> (real), returns sparse array with WM indexing
" ;
Options[S$BINARIZE$REGION] = List[
  Rule["ERROR",10.^-10],            (* -- BINARIZATION CELL ERROR (REAL) *)
  Rule["FUNCTION",Max],             (* -- REPEATED FUNCTION HANDLER, FUNCTION ACTING ON LIST (FUNCTION), DEFAULT VALUE IS "Max" *)
  Rule["DEFAULT",0.0],              (* -- DEFAULT VALUE (REAL) *)
  Rule["ORDERING",False]            (* -- FLAG TO RETURN INDEX TABLE (LOGICAL) *)
] ;
S$BINARIZE$REGION[                  (* -- BINARIZE 3D DATA BASED ON 3RD DIMENSION (RECTANGULAR REGION) *)
  List[HBIN_Integer,HMIN_,HMAX_],   (* -- NUMBER OF HORIZONTAL BINS <HBIN> (INTEGER) AND INTERVAL <HMIN> (REAL) AND <HMAX> (REAL) *)
  List[VBIN_Integer,VMIN_,VMAX_],   (* -- NUMBER OF VERTICLE BINS <VBIN> (INTEGER) AND INTERVAL <VMIN> (REAL) AND <VMAX> (REAL) *)
  DATA_List,                        (* -- INPUT DATA = {...,{X_I,Y_I,Z_I},...} (LIST) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Catch[
  Block[
    {FUNCTION,ERROR,DEFAULT,LOCAL,VALUE,HIND,VIND,SPARSE},
    (* SET REPEATED ENTRIES FUNCTION HANDLER *)
    FUNCTION = OptionValue["FUNCTION"] ;
    (* SET BINARIZATION CELL ERROR *)
    ERROR = N[Subtract[1,OptionValue["ERROR"]]] ;
    (* SET DEFAULT VALUE *)
    DEFAULT = OptionValue["DEFAULT"] ;
    (* SELECTA DATA POINTS INSIDE GIVEN REGION *)
    LOCAL = Select[DATA,Function[Apply[And,Thread[LessEqual[List[HMIN,VMIN],Most[Slot[1]],List[HMAX,VMAX]]]]]] ;
    (* ADD CORNER MARKERS WITH DEFAULT VALUES *)
    LOCAL = Developer`ToPackedArray[Join[LOCAL,List[List[HMIN,VMIN,DEFAULT],List[HMIN,VMAX,DEFAULT],List[HMAX,VMIN,DEFAULT],List[HMAX,VMAX,DEFAULT]]]] ;
    (* SET LIST OF VALUES *)
    VALUE = Part[LOCAL,All,Subtract[0,1]] ;
    (* GENERATE INDEXES *)
    HIND = Plus[1,Floor[ERROR*HBIN*Rescale[Part[LOCAL,All,Plus[0,1]]]]] ;
    VIND = Plus[1,Floor[ERROR*VBIN*Rescale[Part[LOCAL,All,Plus[1,1]]]]] ;
    (* CHECK ORDERING FLAG *)
    If[OptionValue["ORDERING"],Throw[Drop[Transpose[List[HIND,VIND]],Subtract[0,4]]] ; ] ;
    (* CHANGE REPEARTED ENTRIES HANDLER TO CUSTOM FUNCTION *)
    System`SetSystemOptions["SparseArrayOptions" -> {"TreatRepeatedEntries" -> Composition[FUNCTION,List]}] ;
    (* GENERATE SPARSE ARRAY *)
    SPARSE = SparseArray[Rule[Transpose[List[HIND,VIND]],VALUE],List[HBIN,VBIN],DEFAULT] ;
    (* CHANGE REPEARTED ENTRIES HANDLER TO DEFAULT FUNCTION *)
    System`SetSystemOptions["SparseArrayOptions" -> {"TreatRepeatedEntries" -> First}] ;
    (* RETURN, WM INDEXED ARRAY, (1,1) ELEMENT IS AT THE TOP-LEFT CORNER *)
    Reverse[Transpose[SPARSE]]
  ]
] ;
S$BINARIZE$REGION[                  (* -- BINARIZE 3D DATA BASED ON 3RD DIMENSION (SQUARE REGION) *)
  List[BIN_Integer,HDEL_,VDEL_],    (* -- NUMBER OF BINS <BIN> (INTEGER), HORIZONTAL BIN SIZE <HDEL> (REAL) AND VERTICLE BIN SIZE <VDEL> (REAL) *)
  DATA_List,                        (* -- INPUT DATA = {...,{X_I,Y_I,Z_I},...} (LIST) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := S$BINARIZE$REGION[BIN*List[1,-0.5*HDEL,0.5*HDEL],BIN*List[1,-0.5*VDEL,0.5*VDEL],DATA,OPTIONS] ;
(* ################################################################################################################################################################ *)
(* FIRST BOUNDARY LAYER (FMA/AUX) *)
(* ################################################################################################################################################################ *)
Clear[S$FIRST$BOUNDARY$LAYER] ;
S$FIRST$BOUNDARY$LAYER::usage = "S$FIRST$BOUNDARY$LAYER -- 1st boundary layer " ;
S$FIRST$BOUNDARY$LAYER = {{{-1,1},{0,1},{1,1}},{{1,1},{1,0},{1,-1}},{{1,-1},{0,-1},{-1,-1}},{{-1,-1},{-1,0},{-1,1}}} ;
(* ################################################################################################################################################################ *)
(* NEXT BOUNDARY LAYER (FMA/AUX) *)
(* ################################################################################################################################################################ *)
ClearAll[S$NEXT$BOUNDARY$LAYER] ;
S$NEXT$BOUNDARY$LAYER::usage = "S$NEXT$BOUNDARY$LAYER[LAYER] -- generate next boundary layer from given current layer <LAYER> (list) " ;
S$NEXT$BOUNDARY$LAYER = Compile[
  {
    {LAYER,_Integer,3}
  },
  Block[
    {UPLEFT,UPRIGHT,DOWNRIGHT,DOWNLEFT,UP,RIGHT,DOWN,LEFT},
    UPLEFT = Compile`GetElement[LAYER,1,1]+{-1,1} ;
    UP = Transpose[Transpose[Compile`GetElement[LAYER,1]]+{0,1}] ;
    UPRIGHT = Compile`GetElement[LAYER,2,1]+{1,1} ;
    RIGHT = Transpose[Transpose[Compile`GetElement[LAYER,2]]+{1,0}] ;
    DOWNRIGHT = Compile`GetElement[LAYER,3,1]+{1,-1} ;
    DOWN = Transpose[Transpose[Compile`GetElement[LAYER,3]]+{0,-1}] ;
    DOWNLEFT = Compile`GetElement[LAYER,4,1]+{-1,-1} ;
    LEFT = Transpose[Transpose[Compile`GetElement[LAYER,4]]+{-1,0}] ;
    {Join[{UPLEFT},UP,{UPRIGHT}],Join[{UPRIGHT},RIGHT,{DOWNRIGHT}],Join[{DOWNRIGHT},DOWN,{DOWNLEFT}],Join[{DOWNLEFT},LEFT,{UPLEFT}]}
  ],
  CompilationTarget -> "C",
  RuntimeOptions -> "Speed"
] ;
(* ################################################################################################################################################################ *)
(* REORDER ARRAY BASED ON REFFERENCE ARRAY (FMA/AUX) *)
(* ################################################################################################################################################################ *)
ClearAll[S$REORDER$ARRAY] ;
S$REORDER$ARRAY::usage = "S$REORDER$ARRAY[ARRAY,REFERENCE] -- reorder array <ARRAY> (list) based on refference array <REFERENCE> (list) " ;
S$REORDER$ARRAY = Compile[
  {
    {ARRAY,_Real,1},
    {REFERENCE,_Real,1}
  },
  Block[
    {LENGTH,COPY,RESULT,LOOP},
    LENGTH = Length[ARRAY] ;
    COPY = ARRAY ;
    RESULT = COPY ;
    Do[
      Block[
        {POSITION},
        POSITION = First[Ordering[Abs[COPY-Compile`GetElement[REFERENCE,LOOP]],1]] ;
        RESULT[[LOOP]] = Compile`GetElement[COPY,POSITION] ;
        COPY = Delete[COPY,POSITION] ;
      ],
      {LOOP,1,LENGTH}
    ] ;
    RESULT
  ],
  CompilationTarget -> "C",
  RuntimeOptions -> "Speed",
  RuntimeAttributes -> {Listable},
  Parallelization -> True
] ;
(* ################################################################################################################################################################ *)
(* SHIFT LAYER (FMA/AUX) *)
(* ################################################################################################################################################################ *)
ClearAll[S$SHIFT$LAYER] ;
S$SHIFT$LAYER::usage = "S$SHIFT$LAYER[CENTER,LAYER] -- shift layer <LAYER> (list) to center <CENTER> (list) " ;
S$SHIFT$LAYER = Compile[
  {
    {CENTER,_Integer,1},
    {LAYER,_Integer,2}
  },
  Map[Function[CENTER + Slot[1]],LAYER],
  CompilationTarget -> "C",
  RuntimeOptions -> "Speed",
  RuntimeAttributes -> {Listable},
  Parallelization -> True
] ;
(* ################################################################################################################################################################ *)
(* ORDER RANGE (FMA/AUX) *)
(* ################################################################################################################################################################ *)
ClearAll[S$ORDER$RANGE] ;
S$ORDER$RANGE::usage = "S$ORDER$RANGE[ORDER] -- generate order range for given order <ORDER> (integer) " ;
S$ORDER$RANGE = Compile[
  {
    {ORDER,_Integer,0}
  },
  Block[
    {START,END},
    START = Total[8*Range[0,ORDER-1]]+2 ;
    END = START+8*ORDER-1 ;
    {START,END}
  ],
  CompilationTarget -> "C",
  RuntimeOptions -> "Speed",
  RuntimeAttributes -> {Listable},
  Parallelization -> True
] ;
(* ################################################################################################################################################################ *)
(* CHECK INDEX (FMA/AUX) *)
(* ################################################################################################################################################################ *)
ClearAll[S$CHECK$INDEX] ;
S$CHECK$INDEX::usage = "S$CHECK$INDEX[MAXIMUM,CURRENT,INDEX,TABLE,DATA,DEFAULT,LIMIT] -- check data ordering assosiated with index <INDEX> (list) for given maximum order <MAXIMUM> (integer), index order <CURRENT> (integer), table of layers <TABLE> (list), data <DATA> (list), default data value <DEFAULT> (real) and maximum number of checked close indexes to use <LIMIT> (integer) " ;
S$CHECK$INDEX = Compile[
  {
    {MAXIMUM,_Integer,0},
    {CURRENT,_Integer,0},
    {INDEX,_Integer,1},
    {TABLE,_Integer,2},
    {DATA,_Real,3},
    {DEFAULT,_Real,0},
    {LIMIT,_Integer,0}
  },
  Catch[
    Block[
      {LOCAL,HIND,VIND,VALUE,ZERO,LIST,CHECK,LOOP},
      (* SET LOCAL INDEX VALUE *)
      LOCAL = INDEX ;
      (* SET HORIZONTAL AND VERTICLE INDEX VALUES *)
      {HIND,VIND} = LOCAL ;
      (* SET DATA VALUE *)
      VALUE = Compile`GetElement[DATA,MAXIMUM+1-VIND,MAXIMUM+1+HIND] ;
      (* SET DEFAULT DATA VALUE *)
      ZERO = 0.0*VALUE+DEFAULT ;
      (* CHECK 1ST ORDER INDEXES, COMPARE WITH (0,0) *)
      If[
        CURRENT == 1,
        VALUE = VALUE[[Flatten[Map[Ordering[Abs[(VALUE-#)],1]&,Compile`GetElement[DATA,MAXIMUM+1,MAXIMUM+1]]]]] ;
        Throw[VALUE] ;
      ] ;
      (* INITIALIZE LIST OF CLOSE CHECKED INDEXES, INCLUDE (0,0) INDEX *)
      CHECK = {{0,0}} ;
      (* INITIALIZE LIST OF DATA VALUES FOR VALID CLOSE INDEXES *)
      LIST = {ZERO} ;
      (* MAIN LOOP, SELECT CLOSE VALID POINTS BY INCREASING LAYER ORDER AROUND GIVEN INPUT INDEX *)
      Do[
        Block[
          {START,END,SHIFT,CLOSE,SELECT,MEAN},
          (* SET LAYER RANGE *)
          START = Total[8*Range[0,CURRENT-LOOP]]+2 ;
          END = START+8*(CURRENT-LOOP+1)-1 ;
          (* SHIFT LAYER *)
          SHIFT = Map[#+LOCAL &,Take[TABLE,{START,END}]] ;
          (* SET CHECKED TABLE RANGE *)
          If[
            LOOP == 1,
            START = 1 ;
            END = 1 ;
          ] ;
          If[
            LOOP > 1,
            START = Total[8*Range[0,LOOP-2]]+2 ;
            END = START+8*(LOOP-1)-1 ;
          ] ;
          (* SELECT AND ACCUMULATE LIST OF CHECKED INDEXES *)
          CHECK = Join[Take[TABLE,{START,END}],CHECK] ;
          (* FILTER CLOSE INDEXES (INTERSECTION OF CHECKED WITH CURRENT LAYER) *)
          CLOSE = Intersection[CHECK,SHIFT] ;
          (* CONVERT TO DATA FRAME *)
          CLOSE = Transpose[Reverse[(MAXIMUM+1)+{1,-1}*Transpose[CLOSE]]] ;
          (* EXTRACT DATA *)
          SELECT = Join[{1.0+ZERO},Map[Function[(DATA[[Slot[1][[1]],Slot[1][[2]]]])],CLOSE]] ;
          (* DELETE DEFAULT *)
          SELECT = DeleteCases[SELECT,ZERO] ;
          If[Length[SELECT]==1,Continue[]] ;
          SELECT = Rest[SELECT] ;
          (* ACCUMULATE VALID DATA POINTS *)
          LIST = Join[SELECT,LIST] ;
          (* CHECK EXIT CONDITION, NUMBER OF VALID SELECTED IS EQUAL OR LARGER THEN INPUT *)
          If[
            (Length[LIST]-1 >= LIMIT) || (LOOP == 1),
            (* COMPUTE MEAN VALUE OF CLOSE CHECKED POINTS *)
            LIST = Reverse[Most[LIST]] ;
            MEAN = Map[Mean,Transpose[Take[LIST,Min[LIMIT,Length[LIST]]]]] ;
            (* RE-ORDER AND RETURN *)
            Block[
              {LENGTH,COPY,RESULT,MOVE},
              LENGTH = Length[VALUE] ;
              COPY = VALUE ;
              RESULT = COPY ;
              Do[
                Block[
                  {POSITION},
                  POSITION = First[Ordering[Abs[COPY-Compile`GetElement[MEAN,MOVE]],1]] ;
                  RESULT[[MOVE]] = Compile`GetElement[COPY,POSITION] ;
                  COPY = Delete[COPY,POSITION] ;
                ],
                {MOVE,1,LENGTH}
              ] ;
              Throw[RESULT]
            ] ;
          ] ;
        ],
        {LOOP,CURRENT,1,-1}
      ] ;
    ]
  ],
  CompilationTarget -> "C",
  RuntimeOptions -> "Speed",
  RuntimeAttributes -> {Listable},
  Parallelization -> True
] ;
(* ################################################################################################################################################################ *)
(* CHECK LAYER (FMA/AUX) *)
(* ################################################################################################################################################################ *)
ClearAll[S$CHECK$LAYER] ;
S$CHECK$LAYER::usage = "S$CHECK$LAYER[MAXIMUM,CURRENT,TABLE,DATA,DEFAULT,LIMIT] -- check data ordering assosiated with given layer order <CURRENT> (integer) for given maximum order <MAXIMUM> (integer), table of layers <TABLE> (list), data <DATA> (list), default data value <DEFAULT> (real) and maximum number of checked close indexes to use <LIMIT> (integer) " ;
Options[S$CHECK$LAYER] = List[
  Rule["CASES",Null]                (* -- CASES TO DELETE *)
] ;
S$CHECK$LAYER[                      (* -- CHECK LAYER, RETURNS MODIFIED DATA WITH CHECKED LAYER OF SPECIFIED ORDER *)
  MAXIMUM_Integer,                  (* -- MAXIMUM ORDER (INTEGER) *)
  CURRENT_Integer,                  (* -- CURRENT ORDER (INTEGER) *)
  TABLE_List,                       (* -- LAYER TABLE (LIST) *)
  DATA_,                            (* -- DATA (MATRIX) *)
  DEFAULT_?NumericQ,                (* -- DEFAULT DATA ELEMENT VALUE (REAL) *)
  LIMIT_Integer,                    (* -- MAXIMUM NUMBER OF CHECKED CLOSE INDEXES TO USE (INTEGER) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Catch[
  Block[
    {COPY,DIMENSION,DELETE,LAYER,CHECKED},
    (* LOCAL COPY *)
    COPY = DATA ;
    (* DATA DIMENSION *)
    DIMENSION = Last[Dimensions[COPY]] ;
    (* DEFAULT DATA VALUE *)
    DELETE = ConstantArray[DEFAULT,DIMENSION] ;
    (* CURRENT LAYER INDEXES *)
    LAYER = Take[TABLE,S$ORDER$RANGE[CURRENT]] ;
    (* REMOVE DEFAULT *)
    LAYER = Pick[LAYER,Map[Composition[Count[DELETE],List],Extract[DATA,Transpose[Reverse[(MAXIMUM+1)+{1,-1}*Transpose[LAYER]]]]],0] ;
    (* EXIT IF EMPTY *)
    If[SameQ[LAYER,List[]],Throw[DATA]] ;
    (* CHECK OPTIONS *)
    If[UnsameQ[OptionValue["CASES"],Null],LAYER = DeleteCases[LAYER,OptionValue["CASES"]] ; ] ;
    (* EXIT IF EMPTY *)
    If[SameQ[LAYER,List[]],Throw[DATA]] ;
    (* CHECK SELECTED *)
    CHECKED = S$CHECK$INDEX[MAXIMUM,CURRENT,LAYER,TABLE,DATA,DEFAULT,LIMIT] ;
    (* SHIFT INDEXES *)
    LAYER = Transpose[Reverse[(MAXIMUM+1)+{1,-1}*Transpose[LAYER]]] ;
    (* CHANGE DATA *)
    MapThread[Function[(Part[COPY,#1[[1]],#1[[2]]] = #2)],{LAYER,CHECKED}] ;
    (* RETURN *)
    COPY
  ]
] ;
(* ################################################################################################################################################################ *)
(* CHECK DATA (FMA/AUX) *)
(* ################################################################################################################################################################ *)
ClearAll[S$CHECK$DATA] ;
S$CHECK$DATA::usage = "S$CHECK$DATA[MAXIMUM,ORDER,TABLE,DATA,DEFAULT,LIMIT]-- check data ordering assosiated for given maximum order <MAXIMUM> (integer), number of orders to check <ORDER> (integer), table of layers <TABLE> (list), data <DATA> (list), default data value <DEFAULT> (real) and maximum number of checked close indexes to use <LIMIT> (integer) " ;
Options[S$CHECK$DATA] = List[
  Rule["CASES",Null]                (* -- CASES TO DELETE *)
] ;
S$CHECK$DATA[                       (* -- check data,returns modified data *)
  MAXIMUM_Integer,                  (* -- MAXIMUM ORDER (INTEGER) *)
  ORDER_Integer,                    (* -- NUMBER OF ORDERS TO CHECK (INTEGER) *)
  TABLE_List,                       (* -- LAYER TABLE (LIST) *)
  DATA_,                            (* -- DATA (MATRIX) *)
  DEFAULT_?NumericQ,                (* -- DEFAULT DATA ELEMENT VALUE (REAL) *)
  LIMIT_Integer,                    (* -- MAXIMUM NUMBER OF CHECKED CLOSE INDEXES TO USE (INTEGER) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Block[
  {COPY,CURRENT},
  (* LOCAL COPY *)
  COPY = DATA ;
  (* CHECK LAYER-BY-LAYER *)
  Do[COPY = S$CHECK$LAYER[MAXIMUM,CURRENT,TABLE,COPY,DEFAULT,LIMIT,OPTIONS],{CURRENT,1,ORDER}] ;
  (* RETURN *)
  COPY
] /; ORDER <= MAXIMUM ;
(* ################################################################################################################################################################ *)
(* COLOR PLOT *)
(* ################################################################################################################################################################ *)
ClearAll[S$COLOR$PLOT] ;
S$COLOR$PLOT::usage = "
S$COLOR$PLOT[LEVEL,INTERVAL,HRANGE,VRANGE,DATA] -- generate color plot with default data level <LEVEL> (real), color interval <INTERVAL> (list), horizontal data range <HRANGE> (list), vertivle data <VRANGE> (list)  and data <DATA> (matrix)
" ;
Options[S$COLOR$PLOT] = List[
  Rule["LEVEL",LightGray],          (* -- DEFAULT DATA COLOR *)
  Rule["SCHEME","SunsetColors"],    (* -- COLOR SCHEME *)
  Rule["SIZE",400],                 (* -- IMAGE SIZE *)
  Rule["PADDING",40]                (* -- IMAGE PADDING *)
] ;
Options[S$COLOR$PLOT] = Join[
  Options[S$COLOR$PLOT],
  Options[Graphics]
] ;
S$COLOR$PLOT[                       (* -- COLOR PLOT *)
  LEVEL_,                           (* -- LEVEL VALUE *)
  INTERVAL_List,                    (* -- DATA COLOR INTERVAL *)
  HRANGE_List,                      (* -- HORIZONTAL PLOT RANGE *)
  VRANGE_List,                      (* -- VERTICLE PLOT RANGE *)
  DATA_,                            (* -- DATA *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Block[
  {MIN,MAX,COLOR,SELECT,LOCAL,PLOT},
  (* SET COLOR INTERVAL *)
  List[MIN,MAX] = N[INTERVAL] ;
  (* SET COLOR FUNCTION *)
  COLOR = With[
    {
      BACKGROUND = OptionValue["LEVEL"],
      SCHEME = ColorData[OptionValue["SCHEME"]]
    },
    Function[If[SameQ[Slot[1],N[LEVEL]],BACKGROUND,SCHEME[Rescale[Slot[1],INTERVAL]]]]
  ] ;
  (* SET SELECTION FUNCTION *)
  SELECT = Function[Which[SameQ[Slot[1],N[LEVEL]],N[LEVEL],Less[MIN,Slot[1],MAX],Slot[1],GreaterEqual[Slot[1],MAX],MAX,LessEqual[Slot[1],MIN],MIN]] ;  
  (* ADJUST DATA (DATA VALUES OUTSIDE GIVEN INTERVAL ARE SHIFTED TO END POINTS) *)
  LOCAL = Map[SELECT,DATA,List[Plus[1,1]]] ;
  (* GENERATE PLOT *)
  PLOT = ArrayPlot[
    LOCAL,
    Rule[DataRange,List[HRANGE,VRANGE]],
    FilterRules[List[OPTIONS],Options[ArrayPlot]],
    Rule[PlotLegends,False],
    Rule[ColorFunction,COLOR],
    Rule[ColorFunctionScaling,False]
  ] ;
  (* RESULT *)
  Show[
    Graphics[
      List[White,Opacity[0],Apply[Rectangle,Transpose[N[List[HRANGE,VRANGE]]]]],
      FilterRules[List[OPTIONS],Options[Graphics]],
      Rule[AspectRatio,1],
      Rule[Frame,True],
      Rule[FrameStyle,Directive[List[Black]]],
      Rule[FrameTicksStyle,Directive[List[Black,12]]],
      Rule[PlotRangePadding,None],
      Rule[ImagePadding,OptionValue["PADDING"]],
      Rule[ImageSize,OptionValue["SIZE"]]
    ],
    PLOT
  ]
] ;
(* ################################################################################################################################################################ *)
(* COLOR LEGEND *)
(* ################################################################################################################################################################ *)
ClearAll[S$COLOR$LEGEND] ;
S$COLOR$LEGEND::usage = "
S$COLOR$LEGEND[LEVEL,INTERVAL,PLOT,COLOR,SCHEME,PADDING,RATIO] -- generate color plot legend with default data level <LEVEL> (real), color interval <INTERVAL> (list), color plot <PLOT>, background color <COLOR>, color scheme <SCHEME>, plot padding <PADDING> and aspect ratio <RATIO>
" ;
S$COLOR$LEGEND[                     (* -- COLOR LEGEND *)
  LEVEL_,                           (* -- LEVEL VALUE *)
  INTERVAL_,                        (* -- DATA INTERVAL *)
  PLOT_,                            (* -- COLOR PLOT *)
  COLOR_,                           (* -- BACKGROUND COLOR *)
  SCHEME_,                          (* -- SCHEME *)
  PADDIND_,                         (* -- PADDING *)
  RATIO_,                           (* -- SIZE RATIO *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Block[
  {FUNCTION,LEGEND,SIZE},
  FUNCTION = Function[If[SameQ[Slot[1],N[LEVEL]],COLOR,ColorData[SCHEME][Rescale[Slot[1],INTERVAL]]]] ;
  LEGEND = List[Range[First[INTERVAL],Last[INTERVAL],Divide[Total[Abs[INTERVAL]],1000]]] ;
  LEGEND = Reverse[Transpose[LEGEND]] ;
  SIZE = Times[ReplaceAll[ImageSize,Options[PLOT]],ReplaceAll[AspectRatio,Options[PLOT]]] ;
  LEGEND = ArrayPlot[
    LEGEND,
    Rule[PlotLegends,False],
    Rule[Frame,False],
    Rule[ImagePadding,None],
    Rule[PlotRangePadding,None],
    Rule[ColorFunction,FUNCTION],
    Rule[ColorFunctionScaling,False],
    Rule[DataRange,List[List[0,1],INTERVAL]]
  ] ;
  (* RESULT *)
  Show[
    Graphics[
      List[White,Opacity[0],Apply[Rectangle,Transpose[N[List[List[0.0,1.0],INTERVAL]]]]],
      FilterRules[List[OPTIONS],Options[Graphics]],
      Rule[AspectRatio,RATIO],
      Rule[FrameTicks,List[List[All,All],List[None,None]]],
      Rule[Frame,True],
      Rule[FrameStyle,Directive[List[Black]]],
      Rule[FrameTicksStyle,Directive[List[Black,12]]],
      Rule[PlotRangePadding,None],
      Rule[ImagePadding,List[List[0,PADDIND],List[PADDIND,PADDIND]]],
      Rule[ImageSize,List[Automatic,SIZE]]
    ],
    LEGEND
  ]
] ;
(* ################################################################################################################################################################ *)
(* GENERATE RESONANCE DIAGRAM *)
(* ################################################################################################################################################################ *)
ClearAll[S$RESONANCE$DIAGRAM] ;
S$RESONANCE$DIAGRAM::usage = "
S$RESONANCE$DIAGRAM[FREQUENCY,ORDER,RANGE] -- generate resonance diagram with frequencies <FREQUENCY> (list of symbols) for given resonance order <ORDER> (integer) and tune ranges <RANGE> (list)
" ;
Options[S$RESONANCE$DIAGRAM] = List[
  Rule["CONSTRAINT",Function[True]],(* -- SET GENERATION CONSTRAINT (FUNCTION) *)
  Rule["FILTER",Identity],          (* -- FILTER (FUNCTION) *)
  Rule["PRIME",False]               (* -- PRIME FLAG (LOGICAL) *)
] ;
S$RESONANCE$DIAGRAM[                (* -- GENERATE RESONANCE DIAGRAM *)
  FREQUENCY_List,                   (* -- LIST SYMBOLS (LIST), {...,NU_I,...} *)
  ORDER_Integer,                    (* -- RESONANCE ORDER (INTEGER) *)
  RANGE_List,                       (* -- LIST OF TUNE RANGES (LIST), {...,{MIN_I,MAX_I},...} WITH MIN_I < MAX_J *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Catch[
  Block[
    {DIMENSION,LIST,CONSTRAINT,SET,TUNE,FILTER,MIN,MAX,DATA,EQUATION},
    (* SET DIMENSION *)
    DIMENSION = Length[FREQUENCY] ;
    (* CHECK DIMENSION *)
    If[UnsameQ[Length[RANGE],DIMENSION],Throw[$Failed]] ;
    (* GENERATE SET *)
    CONSTRAINT = OptionValue["CONSTRAINT"] ;
    SET = Array[LIST,DIMENSION] ;
    SET = ReplaceAll[SET,Quiet[NSolve[And[Apply[CONSTRAINT,SET],Equal[Total[Abs[SET]],ORDER]],SET,Integers]]] ;
    (* FILTER SET *)
    FILTER = OptionValue["FILTER"] ;
    SET = Map[FILTER,SET] ;
    (* COMPUTE MIN AND MAX VALUES *)
    LIST = Array[TUNE,DIMENSION] ;
    CONSTRAINT = Transpose[RANGE] ;
    CONSTRAINT = Transpose[List[First[CONSTRAINT],LIST,Last[CONSTRAINT]]] ;
    CONSTRAINT = Map[Apply[LessEqual],CONSTRAINT] ;
    MIN = Ceiling[Map[Function[First[Minimize[List[Dot[Slot[1],LIST],Apply[Sequence,CONSTRAINT]],LIST]]],SET]] ;
    MAX = Ceiling[Map[Function[First[Maximize[List[Dot[Slot[1],LIST],Apply[Sequence,CONSTRAINT]],LIST]]],SET]] ;
    (* GENERATE AND SIMPLIFY LINES *)
    DATA = Map[Apply[Range],Transpose[List[MIN,MAX]]] ;
    DATA = Transpose[List[Map[Curry[Dot][LIST],SET],DATA]] ;
    DATA = Map[Composition[Thread,Apply[Equal]],DATA] ;
    DATA = Composition[DeleteDuplicates,Simplify,Flatten][DATA] ;
    (* FORMAT LINES *)
    DATA = Map[
      Function[
        Block[
          {LINE,COEFFICIENT,CONSTANT,FACTOR},
          LINE = Slot[1] ;
          LINE = Apply[Subtract,LINE] ;
          COEFFICIENT = Coefficient[LINE,LIST] ;
          CONSTANT = ReplaceAll[LINE,Thread[Rule[LIST,0]]] ;
          FACTOR = ReplaceAll[FACTOR,First[Solve[Equal[Times[FACTOR,Total[Abs[COEFFICIENT]]],ORDER],FACTOR,Integers]]] ;
          If[
            OptionValue["PRIME"],
            If[
              SameQ[FACTOR,1],
              List[
                Equal[Expand[Times[FACTOR,Dot[COEFFICIENT,LIST]]],Times[Subtract[0,1],FACTOR,CONSTANT]],
                Times[FACTOR,Flatten[List[COEFFICIENT,CONSTANT]]]
              ],
              Nothing
            ],
            List[
              Equal[Expand[Times[FACTOR,Dot[COEFFICIENT,LIST]]],Times[Subtract[0,1],FACTOR,CONSTANT]],
              Times[FACTOR,Flatten[List[COEFFICIENT,CONSTANT]]]
            ]
          ]
        ]
      ],
      DATA
    ] ;
    (* RETURN *)
    DATA = ReplaceAll[DATA,Thread[Rule[LIST,FREQUENCY]]] ;
    DATA = Transpose[DATA] ;
    List[EQUATION,LIST] = DATA ;
    Association[
      Rule["EQUATION",EQUATION],
      Rule["LIST",LIST]
    ]
  ]
] /; GreaterEqual[ORDER,1] ;
(* ################################################################################################################################################################ *)
(* SELECT POINTS CLOSE TO GIVEN RESONANCE *)
(* ################################################################################################################################################################ *)
ClearAll[S$SELECT] ;
S$SELECT::usage = "
S$SELECT[LINE,ERROR,DATA,MASK] -- select points close (frequency distance <ERROR> (real)) to a given 2d resonance line <LINE> (list of integers) for given data <DATA> with frequencies at positions <MASK> (list of integers)
" ;
S$SELECT[                           (* -- SELECT POINTS CLOSE TO GIVEN RESONANCE *)
  LINE_List,                        (* -- RESONANCE LINE (LIST), {A,B,C}: A*QX+B*QY+C = 0 *)
  ERROR_Real,                       (* -- DISTANCE (REAL) *)
  DATA_List,                        (* -- DATA TO PROCESS (LIST) *)
  MASK_List                         (* -- FREQUENCY POSITIONS (LIST OF INTEGERS) *)
] := Block[
  {HOR,VER,INT},
  {HOR,VER,INT} = LINE ;
  Developer`ToPackedArray[Select[DATA,Function[LessEqual[Abs[Plus[Dot[List[HOR,VER],Part[Slot[1],MASK]],INT]],ERROR]]]]
] ;
(* ################################################################################################################################################################ *)
(* PERTURB DATA *)
(* ################################################################################################################################################################ *)
ClearAll[S$PERTURB] ;
S$PERTURB::usage = "
S$PERTURB[DUPLICATION,AMPLITUDE,DATA] -- perturb data <DATA> (list) with amplitudes <AMPLITUDE> (list) using duplication number <DUPLICATION> (integer)
" ;
S$PERTURB[                          (* -- PERTURB INITIAL CONDITIONS *)
  DUPLICATION_Integer,              (* -- DUPLICATION NUMBER (INTEGER) *)
  AMPLITUDE_?MatrixQ,               (* -- LIST OF PERTURBATION RANGES (LIST) *)
  DATA_?MatrixQ                     (* -- DATA TO PERTURB *)
] := Block[
  {SELECT,LENGTH,RANDOM},
  SELECT = DATA ;
  LENGTH = Length[SELECT] ;
  RANDOM = Map[Curry[RandomReal][List[DUPLICATION,LENGTH]],AMPLITUDE] ;
  RANDOM = Transpose[RANDOM] ;
  Flatten[Map[Composition[Transpose,Curry[Plus][Transpose[SELECT]]],RANDOM],1]
] ;
(* ################################################################################################################################################################ *)
(* FILTER DATA *)
(* ################################################################################################################################################################ *)
ClearAll[S$FILTER$DATA] ;
S$FILTER$DATA::usage = "
S$FILTER$DATA[{Q,P},{B_Q,{MIN_Q,MAX_Q}},{B_P,{MIN_P,MAX_P}},DATA] -- filter data <DATA> = {DAT_1,DAT_2,...} = {{X_1,X_2,...,X_Q,...,X_P}_1,...} based on {Q,P} entries within region {{MIN_Q,MAX_Q},{MIN_P,MAX_P}} using number of bins {B_Q,B_P}
" ;
S$FILTER$DATA[                      (* -- FILTER DATA *)
  INDEX_,                           (* -- DATA INDEXES, {Q,P} (INTEGER) *)
  HOR_,                             (* -- 1ST INDEX PARAMETERS, {B_Q,{MIN_Q,MAX_Q}} *)
  VER_,                             (* -- 2ND INDEX PARAMETERS, {B_P,{MIN_P,MAX_P}} *)
  DATA_                             (* -- DATA TO FILTER, {DAT_1,DAT_2,...} = {{X_1,X_2,...,X_Q,...,X_P,...}_1,...} *)
] := Block[
  {IND$HOR,IND$VER,BIN$HOR,MIN$HOR,MAX$HOR,DEL$HOR,BIN$VER,MAX$VER,MIN$VER,DEL$VER,KEY,LIST},
  {IND$HOR,IND$VER} = INDEX ;
  {BIN$HOR,{MIN$HOR,MAX$HOR}} = HOR ;
  DEL$HOR = Divide[Abs[Apply[Subtract,{MIN$HOR,MAX$HOR}]],BIN$HOR] ;
  {BIN$VER,{MIN$VER,MAX$VER}} = VER ;
  DEL$VER = Divide[Abs[Apply[Subtract,{MIN$VER,MAX$VER}]],BIN$VER] ;
  KEY = Part[DATA,All,INDEX] ;
  LIST = Association[Thread[Rule[KEY,DATA]]] ;
  KEY = Select[KEY,Function[LessEqual[MIN$HOR,Part[Slot[1],1],MAX$HOR]]] ;
  KEY = Select[KEY,Function[LessEqual[MIN$VER,Part[Slot[1],2],MAX$VER]]] ;
  KEY = BinLists[KEY,{MIN$HOR,MAX$HOR,DEL$HOR},{MIN$VER,MAX$VER,DEL$VER}] ;
  KEY = Map[First,DeleteCases[Flatten[KEY,1],List[]]] ;
  Map[LIST,KEY]
] ;
(* ################################################################################################################################################################ *)
(* DISCRETE HILBERT TRANSFORM *)
(* ################################################################################################################################################################ *)
ClearAll[S$DHT] ;
S$DHT::usage = "
S$DHT[SIGNAL] -- return discrete Hilbert transform S=U+I*H[U] for given real signal U=<SIGNAL> (list), instantaneous frequency F(N)=arg(S(n)*conjugate(S(n-1))/(2*pi), envelope W(n)=abs(S(n)))
" ;
S$DHT[                              (* -- DISCRETE HILBERT TRANSFORM *)
  SIGNAL_List                       (* -- INPUT SIGNAL (LIST) *)
] := Block[
  {FOURIER,LENGTH,FLAG},
  FOURIER = Rule[FourierParameters,List[1,-1]] ;
  LENGTH = Length[SIGNAL] ;
  FLAG = Boole[EvenQ[LENGTH]];
  Plus[SIGNAL,Times[I,Im[InverseFourier[Times[Fourier[SIGNAL,FOURIER],PadRight[ArrayPad[ConstantArray[Plus[1,1],Subtract[Quotient[LENGTH,Plus[1,1]],FLAG]],List[1,FLAG],1],LENGTH]],FOURIER]]]]
] /; Apply[And,Thread[Equal[Im[SIGNAL],0]]] ;
(* ################################################################################################################################################################ *)
(* POWER SPECTRA LENGTH *)
(* ################################################################################################################################################################ *)
ClearAll[S$PSL] ;
S$PSL::usage = "
S$PSL[SIGNAL] -- compute (normalized with respect to maximum value) power spectra length
" ;
Options[S$PSL] = List[
  Rule["WINDOW",1],                 (* -- APPLY COSINE WINDOW OF GIVEN ORDER (INTEGER) *)
  Rule["PROCESS",True],             (* -- SIGNAL PROCESSING FLAG (LOGICAL) *)
  Rule["LEVEL",N[Power[10,-12]]],   (* -- SPECTRA CHOP LEVEL (REAL) *)
  Rule["COMPLEX",False]             (* -- COMPLEX FLAG (LOGICAL) *)
] ;
S$PSL[                              (* -- POWER SPECTRA LENGTH *)
  SIGNAL_List,                      (* -- INPUT SIGNAL (LIST) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Block[
  {DATA,LENGTH,FOURIER,DELTA},
  DATA = If[OptionValue["PROCESS"],Standardize[SIGNAL],Subtract[SIGNAL,S$WEIGHTED$MEAN[S$COSINE$WINDOW[OptionValue["WINDOW"],Length[SIGNAL]],SIGNAL]]] ;
  LENGTH = Length[DATA] ;
  DATA = Times[DATA,S$COSINE$WINDOW[OptionValue["WINDOW"],LENGTH]] ;
  FOURIER = N[Chop[Rest[S$SPECTRA[DATA,Rule["NORMALIZED",True],Rule["COMPLEX",OptionValue["COMPLEX"]]]],OptionValue["LEVEL"]]] ;
  DELTA = N[Divide[1,LENGTH]] ;
  Total[Sqrt[Plus[Power[Differences[FOURIER],Plus[1,1]],Power[DELTA,Plus[1,1]]]]]
] ;
(* ################################################################################################################################################################ *)
(* FRACTAL DIMENSION ESTIMATION *)
(* ################################################################################################################################################################ *)
ClearAll[S$FDE] ;
S$FDE::usage = "
S$FDE[RANGE,SIGNAL] -- fractal dimension estimation for given list of sample length <RANGE> (list of integers) and signal <SIGNAL>
" ;
Options[S$FDE] = List[
  Rule["WINDOW",1],                 (* -- APPLY COSINE WINDOW OF GIVEN ORDER (INTEGER) *)
  Rule["PROCESS",True],             (* -- SIGNAL PROCESSING FLAG (LOGICAL) *)
  Rule["LEVEL",N[Power[10,-10]]],   (* -- SPECTRA CHOP LEVEL (REAL) *)
  Rule["COMPLEX",False],            (* -- COMPLEX FLAG (LOGICAL) *)
  Rule["METHOD","FIT"]              (* -- METHOD (STRING) *)
] ;
S$FDE[                              (* -- FRACTAL DIMENSION ESTIMATION *)
  RANGE_,                           (* -- LIST OF SAMPLE LENGTH, E.G. {2^10,2^11,...} *)
  SIGNAL_,                          (* -- INPUT SIGNAL (LIST) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Block[
  {VALUE,LENGTH,DATA,MODEL,PARAMETER,FDE},
  VALUE = Log10[N[Divide[1,RANGE]]] ;
  LENGTH = Log10[
    Map[
      Composition[
        Function[
          S$PSL[
            Slot[1],
            Rule["WINDOW",OptionValue["WINDOW"]],
            Rule["PROCESS",OptionValue["PROCESS"]],
            Rule["COMPLEX",OptionValue["COMPLEX"]],
            Rule["LEVEL",OptionValue["LEVEL"]]
          ]
        ],
        Function[Take[SIGNAL,List[1,Slot[1]]]]
      ],
      RANGE
    ]
  ] ;
  DATA = Transpose[List[VALUE,LENGTH]] ;
  If[
    SameQ[OptionValue["METHOD"],"FIT"],
    MODEL = LinearModelFit[DATA,List[PARAMETER],PARAMETER] ;
    FDE = Subtract[1,Last[MODEL["BestFitParameters"]]] ;
  ] ;
  FDE
] ;
(* ################################################################################################################################################################ *)
(* INDEX LEVEL (GLOBAL VARIABLE FOR SALI/GALI INDICATORS) *)
(* ################################################################################################################################################################ *)
(* Skokos C.., Manos T. (2016) *)
(* The Smaller (SALI) and the Generalized (GALI) Alignment Indices: Efficient Methods of Chaos Detection. *)
(* In: Skokos C., Gottwald G., Laskar J. (eds) Chaos Detection and Predictability. Lecture Notes in Physics, vol 915. Springer, Berlin, Heidelberg *)
(* ################################################################################################################################################################ *)
ClearAll[S$INDEX$LEVEL] ;
S$INDEX$LEVEL::usage = "S$INDEX$LEVEL -- index level for SALI/GALI indicators" ;
S$INDEX$LEVEL = N[Power[10,-14]] ;
ClearAll[S$SALI] ;
S$SALI::usage = "S$SALI[{DV1,DV2}] -- compute SALI indicator for given normalized deviation vectors <DV1> and <DV2> " ;
S$SALI = Function[Log10[Plus[S$INDEX$LEVEL,Min[Map[Norm,List[Apply[Plus,Slot[1]],Apply[Subtract,Slot[1]]]]]]]] ;
ClearAll[S$GALI] ;
S$GALI::usage = "S$GALI[{DV1,DV2,...}] -- compute GALI indicator for given normalized deviation vectors <DV1>, <DV2>, ... " ;
S$GALI = Function[Log10[Plus[S$INDEX$LEVEL,Apply[Times,SingularValueList[Transpose[Slot[1]],Rule[Tolerance,0]]]]]] ;
(* ################################################################################################################################################################ *)
(* FLI INDEX (FAST LYAPUNOV INDICATOR) *)
(* ################################################################################################################################################################ *)
(* Lega E., Guzzo M., Froeschlé C. (2016) *)
(* Theory and Applications of the Fast Lyapunov Indicator (FLI) Method. *)
(* In: Skokos C., Gottwald G., Laskar J. (eds) Chaos Detection and Predictability. Lecture Notes in Physics, vol 915. Springer, Berlin, Heidelberg *)
(* ################################################################################################################################################################ *)
ClearAll[S$FLI$FUN] ;
S$FLI$FUN::usage = "S$FLI$FUN -- FLI folding function " ;
S$FLI$FUN = Max ;
ClearAll[S$FLI$ONE] ;
S$FLI$ONE::usage = "S$FLI$ONE[VEC] -- compute FLI indicator using single deviation vector <VEC> given for all iterations, <VEC> = {<VEC>_1, <VEC>_2, ...} " ;
S$FLI$ONE = Composition[FoldList[Composition[S$FLI$FUN,Flatten,List]],Log,Map[Norm]] ;
ClearAll[S$FLI$ALL] ;
S$FLI$ALL::usage = "S$FLI$ALL[VEC] -- compute FLI indicator using several (orthonormal) deviation vectors <VEC> given for all iterations, <VEC> = {<VEC>_1, <VEC>_2, ...}, <VEC>_I = {<VEC>_I_1, <VEC>_I_2, ...} " ;
S$FLI$ALL = Composition[Map[Max],Transpose,Map[S$FLI$ONE],List] ;
ClearAll[S$FLI$LAS] ;
S$FLI$LAS::usage = "S$FLI$LAS[VEC] -- compute FLI indicator using several (orthonormal) deviation vectors <VEC> given for all iterations, <VEC> = {<VEC>_1, <VEC>_2, ...}, <VEC>_I = {<VEC>_I_1, <VEC>_I_2, ...} (returns only last value) " ;
S$FLI$LAS = Composition[Last,S$FLI$ALL] ;
ClearAll[S$FLI] ;
S$FLI::usage = "S$FLI[VEC] -- compute FLI indicator using several (orthonormal) deviation vectors <VEC> given for all iterations, <VEC> = {<VEC>_1, <VEC>_2, ...}, <VEC>_I = {<VEC>_I_1, <VEC>_I_2, ...} (returns only last value) " ;
S$FLI = S$FLI$LAS ;
(* ################################################################################################################################################################ *)
(* OPTICS *)
(* ################################################################################################################################################################ *)
ClearAll[S$WOLSKI] ;
S$WOLSKI::usage = "
S$WOLSKI[MATRIX] -- compute coupled optics in Wolski representation for given single-turn (2N)D real symplectic matrix <MATRIX>
" ;
Options[S$WOLSKI] = List[
  Rule["TOLERANCE",10.0^-15],       (* -- CHOP TOLERANCE (REAL) *)
  Rule["IDENTIFY",True],            (* -- FLAG TO IDENTIFY PLANES (LOGICAL) *)
  Rule["METHOD","TWISS"],           (* -- PLANES IDENTIFICATION METHOD ("TWISS", "GUESS" OR "TRACK") *)
  Rule["GUESS",Null]                (* -- LIST GUESS TUNES FOR PLANES (LIST OF REALS) *)
] ;
S$WOLSKI[
  MATRIX_?MatrixQ,                  (* -- SINGLE-TURN (2N)D REAL SYMPLECTIC MATRIX *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Catch[
  Block[
    {TOLERANCE, DIMENSION, SYMPLECTIC, COMPLEX, PROJECTION, ROTATION, VALUE, VECTOR, TUNE, FACTOR, EIGEN, NORMAL, WOLSKI, ORDERING, GAUGE, START, ORBIT},
    TOLERANCE = OptionValue["TOLERANCE"] ;
    DIMENSION = 1/2*Length[MATRIX] ;
    SYMPLECTIC = ArrayFlatten[Outer[Times,IdentityMatrix[DIMENSION], {{0, 1}, {-1, 0}}]] ;
    COMPLEX = 1/Sqrt[2]*ArrayFlatten[Outer[Times, IdentityMatrix[DIMENSION], {{1, +I}, {1, -I}}]] ;
    PROJECTION = Table[SparseArray[{2*{X, X} + 1, 2*{X, X} + 2} -> {1, 1}, 2*DIMENSION], {X, 0, DIMENSION - 1}] ;
    ROTATION = Table[SparseArray[{2*{X, X} + 1, 2*{X, X} + 2, {Y_, Y_}} -> {-1, -1, 1}, 2*DIMENSION], {X, 0, DIMENSION - 1}] ;
    {VALUE, VECTOR} = Chop[Eigensystem[MATRIX], TOLERANCE] ;
    {VALUE, VECTOR} = Transpose[SortBy[Transpose[{VALUE, VECTOR}], Composition[Abs, Log, First]]] ;
    VALUE = Partition[VALUE, 2] ;
    VECTOR = Partition[VECTOR, 2] ;
    TUNE = Map[Function[Abs[Im[Log[Mean[{First[Slot[1]], Conjugate[Last[Slot[1]]]}]]]]], VALUE] ;
    FACTOR = Map[Function[Sqrt[Abs[Chop[-I*Dot[First[Slot[1]], SYMPLECTIC, Last[Slot[1]]], TOLERANCE]]]], VECTOR] ;
    EIGEN = Transpose[Flatten[VECTOR/FACTOR, 1]] ;
    NORMAL = Chop[Dot[EIGEN, COMPLEX], TOLERANCE] ;
    WOLSKI = Map[Function[Dot[NORMAL, Slot[1], Transpose[NORMAL]]], PROJECTION] ;
    If[
      Not[OptionValue["IDENTIFY"]],
      Throw[Association[Rule["TUNE", TUNE], Rule["WOLSKI", WOLSKI], Rule["NORMAL", NORMAL], Rule["EIGEN", EIGEN]]] ;
    ] ;
    If[
      SameQ[OptionValue["METHOD"], "TWISS"],
      ORDERING  = Flatten[Map[Function[Position[Slot[1], Max[Slot[1]]]], Table[WOLSKI[[X, 2*Y - 1, 2*Y - 1]], {X, 1, DIMENSION}, {Y, 1, DIMENSION}]]] ;
    ] ;
    If[
      SameQ[OptionValue["METHOD"], "GUESS"],
      ORDERING = OptionValue["GUESS"] ;
      ORDERING = Flatten[Map[Function[Position[Abs[ORDERING - Slot[1]], Min[Abs[ORDERING - Slot[1]]]]], TUNE]]
    ] ;
    If[
      SameQ[OptionValue["METHOD"], "TRACK"],
      ORDERING = ConstantArray[0, DIMENSION] ;
      Do[
        {
          START = ConstantArray[0, 2*DIMENSION] ;
          START[[2*X-1]] = 1 ;
          ORBIT = NestList[Function[Dot[MATRIX, Slot[1]]], START, S$WINDOW$LENGTH] ;
          ORBIT = Part[ORBIT, All, 2*X-1] ;
          ORDERING[[X]] = 2*Pi*S$FREQUENCY[ORBIT] ;
        },
        {X, 1, DIMENSION}
      ] ;
      ORDERING = Flatten[Map[Function[Position[Abs[ORDERING - Slot[1]], Min[Abs[ORDERING - Slot[1]]]]], TUNE]] ;
    ] ;
    If[
      UnsameQ[ORDERING, Range[DIMENSION]],
      VALUE = VALUE[[ORDERING]] ;
      TUNE = Map[Function[Abs[Im[Log[Mean[{First[Slot[1]], Conjugate[Last[Slot[1]]]}]]]]], VALUE] ;
      VECTOR = VECTOR[[ORDERING]] ;
      FACTOR = FACTOR[[ORDERING]] ;
      EIGEN = Transpose[Flatten[VECTOR/FACTOR, 1]] ;
      NORMAL = Chop[Dot[EIGEN, COMPLEX], TOLERANCE] ;
      WOLSKI = Map[Function[Dot[NORMAL, Slot[1], Transpose[NORMAL]]], PROJECTION] ;
    ] ;
    GAUGE = Sign[Table[NORMAL[[2*X - 1, 2*X - 1]], {X, 0, DIMENSION - 1}]] ;
    If[
      UnsameQ[GAUGE, ConstantArray[1, DIMENSION]],
      NORMAL = Dot[NORMAL, Apply[Times,Pick[ROTATION, GAUGE, 1]]] ;
      WOLSKI = Map[Function[Dot[NORMAL, Slot[1], Transpose[NORMAL]]], PROJECTION] ;
    ] ;
    Association[Rule["TUNE", TUNE], Rule["WOLSKI", WOLSKI], Rule["NORMAL", NORMAL], Rule["EIGEN", EIGEN]]
  ]
] ;
ClearAll[S$ADVANCE] ;
S$ADVANCE::usage = "
S$ADVANCE[MATRIX12, WOLSKI1] -- advance Wolski representation <WOLSKI1> at point 1 throught given trasport matrix <MATRIX12>
" ;
S$ADVANCE[MATRIX12_, WOLSKI1_] := Map[Function[Dot[MATRIX12, Slot[1], Transpose[MATRIX12]]], WOLSKI1] ;
ClearAll[S$SHIFT] ;
S$SHIFT::usage = "
S$SHIFT[MATRIX12, NORMAL1] -- compute phase advance for given transport matrix <MATRIX12> and normalization matrix <NORMAL1> at point 1
" ;
S$SHIFT[MATRIX12_, NORMAL1_] := Block[
	{ROTATION},
	ROTATION = MATRIX12.NORMAL1 ;
	Table[ArcTan[ROTATION[[2*X+1, 2*X+1]], ROTATION[[2*X+1, 2*X+2]]], {X, 0, 1/2*Length[MATRIX12] - 1}]
] ;
ClearAll[S$LB$NORMAL] ;
S$LB$NORMAL::usage = "
S$LB$NORMAL[{B1X, B1Y, B2X, B2Y, A1X, A1Y, A2X, A2Y, U, CV1, SV1, CV2, SV2}] -- compute normalization matrix from Lebedev-Bogacz representation, where CVI = Cos[VI] and SVI = Sin[VI]
" ;
S$LB$NORMAL[{B1X_, B1Y_, B2X_, B2Y_, A1X_, A1Y_, A2X_, A2Y_, U_, CV1_, SV1_, CV2_, SV2_}] := {
  {Sqrt[B1X], 0, Sqrt[B2X]*CV2, -Sqrt[B2X]*SV2},
  {-(A1X/Sqrt[B1X]), (1-U)/Sqrt[B1X], If[B2X > 0, (-A2X*CV2 + SV2*U)/Sqrt[B2X], 0], If[B2X >0, (A2X*SV2 + CV2*U)/Sqrt[B2X], 0]},
  {Sqrt[B1Y]*CV1, -Sqrt[B1Y]*SV1, Sqrt[B2Y], 0},
  {If[B1Y > 0, (-A1Y*CV1 + SV1*U)/Sqrt[B1Y], 0], If[B1Y > 0, (A1Y*SV1 + CV1*U)/Sqrt[B1Y], 0], -(A2Y/Sqrt[B2Y]),(1-U)/Sqrt[B2Y]}
} ;
ClearAll[S$LB$CONVERT] ;
S$LB$CONVERT::usage = "
S$LB$CONVERT[{W1, W2}] -- convert Wolski representation to Lebedev-Bogacz representation
S$LB$CONVERT[{B1X, B1Y, B2X, B2Y, A1X, A1Y, A2X, A2Y, U, CV1, SV1, CV2, SV2}] -- convert Lebedev-Bogacz representation to Wolski representation, where CVI = Cos[VI] and SVI = Sin[VI]
" ;
S$LB$CONVERT[{W1_, W2_}] := Block[
  {B1X, B1Y, B2X, B2Y, A1X, A1Y, A2X, A2Y, U, CV1, SV1, CV2, SV2},
  B1X = +W1[[1,1]] ;
  B1Y = +W1[[3,3]] ;
  B2X = +W2[[1,1]] ;
  B2Y = +W2[[3,3]] ;
  A1X = -W1[[1,2]] ;
  A1Y = -W1[[3,4]] ;
  A2X = -W2[[1,2]] ;
  A2Y = -W2[[3,4]] ;
  U = 1/2*(1 + A1X^2 - A1Y^2 - B1X*W1[[2,2]] + B1Y*W1[[4,4]]) ;
  {CV1, SV1} = If[B1Y > 0, {1/Sqrt[B1X*B1Y]*W1[[1,3]], 1/U*(A1Y*CV1 + 1/Sqrt[B1X]*(Sqrt[B1Y]*W1[[1,4]]))}, {-1, 0}] ;
  {CV2, SV2} = If[B2X > 0, {1/Sqrt[B2X*B2Y]*W2[[1,3]], 1/U*(A2X*CV2 + 1/Sqrt[B2Y]*(Sqrt[B2X]*W2[[2,3]]))}, {+1, 0}] ;
  {B1X, B1Y, B2X, B2Y, A1X, A1Y, A2X, A2Y, U, CV1, SV1, CV2, SV2}
] ;
S$LB$CONVERT[{B1X_, B1Y_, B2X_, B2Y_, A1X_, A1Y_, A2X_, A2Y_, U_, CV1_, SV1_, CV2_, SV2_}] := Block[
  {NORMAL, P1, P2, W1, W2},
  NORMAL = S$LB$NORMAL[{B1X, B1Y, B2X, B2Y, A1X, A1Y, A2X, A2Y, U, CV1, SV1, CV2, SV2}] ;
  P1 = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}} ;
  P2 = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}} ;
  W1 = Dot[NORMAL, P1, Transpose[NORMAL]] ;
  W2 = Dot[NORMAL, P2, Transpose[NORMAL]] ;
  {W1, W2}
] ;
(* ################################################################################################################################################################ *)
(* COMPLETION *)
(* ################################################################################################################################################################ *)
ClearAll[S$COMPLETION] ;
S$COMPLETION = Composition[FE`Evaluate,FEPrivate`AddSpecialArgCompletion] ;
ClearAll[S$COMLETION$OPTION] ;
S$COMLETION$OPTION[] := Block[
  {},
  S$COMPLETION[Rule["S$SET$VARIABLES",List[0,0,0,0,0,0,0,Apply[Sequence,ConstantArray[Map[First,Options[S$SET$VARIABLES]],Length[Options[S$SET$VARIABLES]]]]]]] ;
  S$COMPLETION[Rule["S$SPECTRA",List[0,Apply[Sequence,ConstantArray[Map[First,Options[S$SPECTRA]],Length[Options[S$SPECTRA]]]]]]] ;
  S$COMPLETION[Rule["S$FREQUENCY",List[0,Apply[Sequence,ConstantArray[Map[First,Options[S$FREQUENCY]],Length[Options[S$FREQUENCY]]]]]]] ;
  S$COMPLETION[Rule["S$DECOMPOSITION",List[0,Apply[Sequence,ConstantArray[Map[First,Options[S$DECOMPOSITION]],Length[Options[S$DECOMPOSITION]]]]]]] ;
  S$COMPLETION[Rule["S$GENERATE$MODEL",List[0,Apply[Sequence,ConstantArray[Map[First,Options[S$GENERATE$MODEL]],Length[Options[S$GENERATE$MODEL]]]]]]] ;
  S$COMPLETION[Rule["S$FIT$MODEL",List[0,0,Apply[Sequence,ConstantArray[Map[First,Options[S$FIT$MODEL]],Length[Options[S$FIT$MODEL]]]]]]] ;
  S$COMPLETION[Rule["S$QUASIPERIODIC$SUBTRACT",List[0,0,Apply[Sequence,ConstantArray[Map[First,Options[S$QUASIPERIODIC$SUBTRACT]],Length[Options[S$QUASIPERIODIC$SUBTRACT]]]]]]] ;
  S$COMPLETION[Rule["S$QUASIPERIODIC$PEAK",List[0,0,Apply[Sequence,ConstantArray[Map[First,Options[S$QUASIPERIODIC$PEAK]],Length[Options[S$QUASIPERIODIC$PEAK]]]]]]] ;
  S$COMPLETION[Rule["S$QUASIPERIODIC$HARMONICS",List[0,0,Apply[Sequence,ConstantArray[Map[First,Options[S$QUASIPERIODIC$HARMONICS]],Length[Options[S$QUASIPERIODIC$HARMONICS]]]]]]] ;
  S$COMPLETION[Rule["S$PRONY",List[0,0,Apply[Sequence,ConstantArray[Map[First,Options[S$PRONY]],Length[Options[S$PRONY]]]]]]] ;
  S$COMPLETION[Rule["S$FMA",List[0,Apply[Sequence,ConstantArray[Map[First,Options[S$FMA]],Length[Options[S$FMA]]]]]]] ;
  S$COMPLETION[Rule["S$BINARIZE$REGION",List[0,0,0,0,0,Apply[Sequence,ConstantArray[Map[First,Options[S$BINARIZE$REGION]],Length[Options[S$BINARIZE$REGION]]]]]]] ;
  S$COMPLETION[Rule["S$COLOR$PLOT",List[0,0,0,0,0,Apply[Sequence,ConstantArray[Map[First,Options[S$COLOR$PLOT]],Length[Options[S$COLOR$PLOT]]]]]]] ;
  S$COMPLETION[Rule["S$RESONANCE$DIAGRAM",List[0,0,0,Apply[Sequence,ConstantArray[Map[First,Options[S$RESONANCE$DIAGRAM]],Length[Options[S$RESONANCE$DIAGRAM]]]]]]] ;
  S$COMPLETION[Rule["S$PSL",List[0,Apply[Sequence,ConstantArray[Map[First,Options[S$PSL]],Length[Options[S$PSL]]]]]]] ;
  S$COMPLETION[Rule["S$FDE",List[0,Apply[Sequence,ConstantArray[Map[First,Options[S$FDE]],Length[Options[S$FDE]]]]]]] ;
] ;
ClearAll[S$COMPLETION$RULE] ;
S$COMPLETION$RULE[DIRECTORY_String] := Block[
  {PATH,FILE,LIST},
  PATH = FileNameJoin[List[NotebookDirectory[]]] ;
  If[Not[FileExistsQ[DIRECTORY]],CreateDirectory[FileNameJoin[List[NotebookDirectory[],DIRECTORY]]]] ;
  PATH = FileNameJoin[List[NotebookDirectory[],DIRECTORY]] ;
  FILE = FileNameJoin[List[DIRECTORY,"S$SET$VARIABLES.m"]] ;
  LIST = List[
    Rule["VERBOSE",List["True","False"]]
  ] ;
  Export[FILE,LIST] ;
  FILE = FileNameJoin[List[PATH,"S$SPECTRA.m"]] ;
  LIST = List[
    Rule["FOURIER",List["{-1,1}","{0,1}","{1,-1}"]],
    Rule["NORMALIZED",List["True","False"]],
    Rule["COMPLEX",List["True","False"]],
    Rule["2D",List["True","False"]],
    Rule["LEVEL",List["-15.0","-10.0"]]
  ] ;
  Export[FILE,LIST] ;
  FILE = FileNameJoin[List[PATH,"S$FREQUENCY.m"]] ;
  LIST = List[
    Rule["PROCESS",List["True","False"]],
    Rule["PADDING",List["0"]],
    Rule["PEAK",List["0","1"]],
    Rule["METHOD",List["\"FFT\"","\"REFINE\"","\"PARABOLA\"","\"INTERPOLATION\"","\"FIT\""]]
  ] ;
  Export[FILE,LIST] ;
  FILE = FileNameJoin[List[PATH,"S$DECOMPOSITION.m"]] ;
  LIST = List[
    Rule["TYPE",List["\"PEAK\"","\"MAXIMUM\""]],
    Rule["CHOP",List["10.^-10"]],
    Rule["PADDING",List["0"]],
    Rule["METHOD",List["\"FFT\"","\"REFINE\"","\"PARABOLA\"","\"INTERPOLATION\"","\"FIT\""]]
  ] ;
  Export[FILE,LIST] ;
  FILE = FileNameJoin[List[PATH,"S$GENERATE$MODEL.m"]] ;
  LIST = List[
    Rule["DAMPING",List["True","False"]]
  ] ;
  Export[FILE,LIST] ;
  FILE = FileNameJoin[List[PATH,"S$FIT$MODEL.m"]] ;
  LIST = List[
    Rule["DAMPING",List["True","False"]],
    Rule["NONLINEAR",List["True","False"]],
    Rule["DECOMPOSITION",List["True","False"]],
    Rule["TYPE",List["\"PEAK\"","\"MAXIMUM\""]],
    Rule["CHOP",List["10.^-10"]],
    Rule["PADDING",List["0"]],
    Rule["METHOD",List["\"FFT\"","\"REFINE\"","\"PARABOLA\"","\"INTERPOLATION\"","\"FIT\""]]
  ] ;
  Export[FILE,LIST] ;
  FILE = FileNameJoin[List[PATH,"S$QUASIPERIODIC$SUBTRACT.m"]] ;
  LIST = List[
    Rule["DAMPING",List["True","False"]],
    Rule["NONLINEAR",List["True","False"]],
    Rule["DECOMPOSITION",List["True","False"]],
    Rule["TYPE",List["\"PEAK\"","\"MAXIMUM\""]],
    Rule["CHOP",List["10.^-10"]],
    Rule["PADDING",List["0"]],
    Rule["METHOD",List["\"FFT\"","\"REFINE\"","\"PARABOLA\"","\"INTERPOLATION\"","\"FIT\""]]
  ] ;
  Export[FILE,LIST] ;
  FILE = FileNameJoin[List[PATH,"S$QUASIPERIODIC$PEAK.m"]] ;
  LIST = List[
    Rule["DAMPING",List["True","False"]],
    Rule["NONLINEAR",List["True","False"]],
    Rule["DECOMPOSITION",List["True","False"]],
    Rule["TYPE",List["\"PEAK\"","\"MAXIMUM\""]],
    Rule["CHOP",List["10.^-10"]],
    Rule["PADDING",List["0"]],
    Rule["METHOD",List["\"FFT\"","\"REFINE\"","\"PARABOLA\"","\"INTERPOLATION\"","\"FIT\""]]
  ] ;
  Export[FILE,LIST] ;
  FILE = FileNameJoin[List[PATH,"S$QUASIPERIODIC$HARMONICS.m"]] ;
  LIST = List[
    Rule["DAMPING",List["True","False"]],
    Rule["NONLINEAR",List["True","False"]],
    Rule["DECOMPOSITION",List["True","False"]],
    Rule["TYPE",List["\"PEAK\"","\"MAXIMUM\""]],
    Rule["CHOP",List["10.^-10"]],
    Rule["PADDING",List["0"]],
    Rule["METHOD",List["\"FFT\"","\"REFINE\"","\"PARABOLA\"","\"INTERPOLATION\"","\"FIT\""]]
  ] ;
  Export[FILE,LIST] ;
  FILE = FileNameJoin[List[PATH,"S$PRONY.m"]] ;
  LIST = List[
    Rule["CHOP",List["10.^-10"]]
  ] ;
  Export[FILE,LIST] ;
  FILE = FileNameJoin[List[PATH,"S$FMA.m"]] ;
  LIST = List[
    Rule["TYPE",List["\"PEAK\"","\"MAXIMUM\""]],
    Rule["CHOP",List["10.^-10"]],
    Rule["PADDING",List["0"]],
    Rule["METHOD",List["\"FFT\"","\"REFINE\"","\"PARABOLA\"","\"INTERPOLATION\"","\"FIT\""]],
    Rule["SCAN",List["0"]],
    Rule["MODE",List["\"SINGLE\"","\"MULTIPLE\""]],
    Rule["PEAK",List["1","2","3"]]
  ] ;
  Export[FILE,LIST] ;
  FILE = FileNameJoin[List[PATH,"S$BINARIZE$REGION.m"]] ;
  LIST = List[
    Rule["FUNCTION",List["Max","Mean","Min"]],
    Rule["PARALLEL",List["True","False"]],
    Rule["LEVEL",List["0."]],
    Rule["ERROR",List["10.^-10"]]
  ] ;
  Export[FILE,LIST] ;
  FILE = FileNameJoin[List[PATH,"S$COLOR$PLOT.m"]] ;
  LIST = List[
    Rule["LEVEL",List["LightGray"]],
    Rule["SCHEME",List["\"SunsetColors\""]],
    Rule["SIZE",List["1000","2000","3000"]],
    Rule["PADDING",List["20","40"]]
  ] ;
  Export[FILE,LIST] ;
  FILE = FileNameJoin[List[PATH,"S$RESONANCE$DIAGRAM.m"]] ;
  LIST = List[
    Rule["CONSTRAINT",List["(True&)"]],
    Rule["FILTER",List["Identity"]],
    Rule["PRIME",List["True","False"]]
  ] ;
  Export[FILE,LIST] ;
  FILE = FileNameJoin[List[PATH,"S$PSL.m"]] ;
  LIST = List[
    Rule["WINDOW",List["0","1"]],
    Rule["STANDARDIZE",List["True","False"]],
    Rule["LEVEL",List["10.^-10"]],
    Rule["COMPLEX",List["True","False"]]
  ] ;
  Export[FILE,LIST] ;
  FILE = FileNameJoin[List[PATH,"S$FDE.m"]] ;
  LIST = List[
    Rule["WINDOW",List["0","1"]],
    Rule["STANDARDIZE",List["True","False"]],
    Rule["LEVEL",List["10.^-10"]],
    Rule["COMPLEX",List["True","False"]],
    Rule["METHOD",List["\"FIT\""]]
  ] ;
  Export[FILE,LIST] ;
  StringTemplate["cp `1`/*.m `2`/"][PATH,FileNameJoin[List[$InstallationDirectory,"SystemFiles","Components","AutoCompletionData","Main","OptionValues"]]]
] ;
(* ################################################################################################################################################################ *)
(* COMPILED *)
(* ################################################################################################################################################################ *)
ClearAll[S$COMPILED$FFT] ;
S$COMPILED$FFT::usage = "
S$COMPILED$FFT[LENGTH,DIRECTION,SAMPLE] -- radix-two fft of sample <SAMPLE> = {...,SR_I,SI_I,...} (list) of length 2^(N+1), <LENGTH> = 2^N (integer) and direction <DIRECTION> = +1/-1 (integer) corresponds to {1,+1}/{1,-1} fourier parameters of Fourier[] function
" ;
S$COMPILED$FFT = Compile[
  {
    {LENGTH,_Integer},
    {DIRECTION,_Integer},
    {SAMPLE,_Real,1}
  },
  Block[
    {PI,LOCAL,II,IJ,IN,IM,LIMIT,STEP,PIM,PRE,ANG,WR,WI,WPR,WPI,WCOPY,FACTOR},
    PI = Pi ;
    LOCAL = SAMPLE ;
    IN = 2*LENGTH ;
    IJ = 1 ;
    FACTOR = 2.0*PI*DIRECTION ;
    Do[
      Block[
        {},
        If[
          IJ > II,
          PRE = Compile`GetElement[LOCAL,IJ] ;
          PIM = Compile`GetElement[LOCAL,IJ+1] ;
          LOCAL[[IJ]] = Compile`GetElement[LOCAL,II] ;
          LOCAL[[IJ+1]] = Compile`GetElement[LOCAL,II+1] ;
          LOCAL[[II]] = PRE ;
          LOCAL[[II+1]] = PIM ;
        ] ;
        IM = Floor[IN/2] ;
        While[
          IM >= 2 && IJ > IM,
          IJ = IJ-IM ;
          IM = Floor[IM/2] ;
        ] ;
        IJ = IJ+IM ;
      ],
      {II,1,IN,2}
    ] ;
    LIMIT = 2 ;
    While[
      IN > LIMIT,
      STEP = 2*LIMIT ;
      ANG = FACTOR/LIMIT ;
      WPI = Sin[ANG] ;
      ANG = Sin[0.5*ANG] ;
      WPR = -2.0*ANG*ANG ;
      WR = 1.0 ;
      WI = 0.0 ;
      Do[
        Block[
          {},
          Do[
            Block[
              {},
              IJ = II + LIMIT ;
              PRE = WR*Compile`GetElement[LOCAL,IJ]-WI*Compile`GetElement[LOCAL,IJ+1] ;
              PIM = WR*Compile`GetElement[LOCAL,IJ+1]+WI*Compile`GetElement[LOCAL,IJ] ;
              LOCAL[[IJ]] = Compile`GetElement[LOCAL,II]-PRE ;
              LOCAL[[IJ+1]] = Compile`GetElement[LOCAL,II+1]-PIM ;
              LOCAL[[II]] = Compile`GetElement[LOCAL,II] + PRE ;
              LOCAL[[II+1]] = Compile`GetElement[LOCAL,II+1]+PIM ;
            ],
            {II,IM,IN,STEP}
          ] ;
          WCOPY = WR ;
          WR = WR*WPR-WI*WPI+WR ;
          WI = WI*WPR+WCOPY*WPI+WI ;
        ],
        {IM,1,LIMIT,2}
      ] ;
      LIMIT = STEP ;
    ] ;
    LOCAL
  ],
  CompilationTarget -> "C",
  RuntimeOptions -> "Speed"
] ;
ClearAll[S$COMPILED$FFRFT] ;
S$COMPILED$FFRFT::usage = "
S$COMPILED$FFRFT[LENGTH,ARGUMENT,SAMPLE] -- ffrft of sample <SAMPLE> = {...,SR_I,SI_I,...} (list) of length 2^(N+1), <LENGTH> = 2^N (integer) and argument <ARGUMENT> (real) corresponds to {1,<ARGUMENT>} fourier parameters of Fourier[] function
" ;
S$COMPILED$FFRFT = Compile[
  {
    {LENGTH,_Integer},
    {ARGUMENT,_Real},
    {SAMPLE,_Real,1}
  },
  Block[
    {PI,LOCAL,FACTOR,MUL,COS$MUL,SIN$MUL,ONE,TWO,TRE,COPY},
    PI = Pi ;
    LOCAL = SAMPLE ;
    FACTOR = 1.0/LENGTH*ARGUMENT*PI ;
    MUL = FACTOR*(Range[LENGTH]-1.0)^2 ;
    COS$MUL = Cos[MUL] ;
    SIN$MUL = Sin[MUL] ;
    ONE = TWO = TRE = Table[0.0,List[4*LENGTH]] ;
    COPY = Table[0.0,List[2*LENGTH]] ;
    ONE[[1;;2*LENGTH;;2]] = LOCAL[[1;;2*LENGTH;;2]]*COS$MUL-LOCAL[[2;;2*LENGTH;;2]]*SIN$MUL ;
    ONE[[2;;2*LENGTH;;2]] = LOCAL[[1;;2*LENGTH;;2]]*SIN$MUL+LOCAL[[2;;2*LENGTH;;2]]*COS$MUL ;
    TWO[[1;;2*LENGTH;;2]] = + COS$MUL ;
    TWO[[2;;2*LENGTH;;2]] = - SIN$MUL ;
    MUL = -FACTOR*(Range[LENGTH+1,2*LENGTH,1]-1.0-2.0*LENGTH)^2 ;
    TWO[[2*LENGTH+1;;4*LENGTH;;2]] = Cos[MUL] ;
    TWO[[2*LENGTH+2;;4*LENGTH;;2]] = Sin[MUL] ;
    ONE = S$COMPILED$FFT[2*LENGTH,+1,ONE] ;
    TWO = S$COMPILED$FFT[2*LENGTH,+1,TWO] ;
    TRE = ONE ;
    ONE[[1;;4*LENGTH;;2]] = TRE[[1;;4*LENGTH;;2]]*TWO[[1;;4*LENGTH;;2]]-TRE[[2;;4*LENGTH;;2]]*TWO[[2;;4*LENGTH;;2]] ;
    ONE[[2;;4*LENGTH;;2]] = TRE[[1;;4*LENGTH;;2]]*TWO[[2;;4*LENGTH;;2]]+TRE[[2;;4*LENGTH;;2]]*TWO[[1;;4*LENGTH;;2]] ;
    ONE = S$COMPILED$FFT[2*LENGTH,-1,ONE] ;
    LOCAL = ONE[[1;;2*LENGTH;;1]]/(2.0*LENGTH) ;
    COPY = LOCAL ;
    LOCAL[[1;;2*LENGTH;;2]] = COPY[[1;;2*LENGTH;;2]]*COS$MUL-COPY[[2;;2*LENGTH;;2]]*SIN$MUL ;
    LOCAL[[2;;2*LENGTH;;2]] = COPY[[1;;2*LENGTH;;2]]*SIN$MUL+COPY[[2;;2*LENGTH;;2]]*COS$MUL ;
    LOCAL
  ],
  CompilationTarget -> "C",
  CompilationOptions -> {"InlineExternalDefinitions" -> True},
  RuntimeOptions -> "Speed"
] ;
ClearAll[S$COMPILED$WINDOW] ;
S$COMPILED$WINDOW::usage = "
S$COMPILED$WINDOW[<LENGTH>,<ORDER>] -- generate cosine window of order <ORDER> (real) and length <LENGTH> (integer)
" ;
S$COMPILED$WINDOW = Compile[
  {
    {LENGTH,_Integer},
    {ORDER,_Real}
  },
  If[
    SameQ[N[ORDER],N[0]],
    Table[N[1],LENGTH],
    Times[
      Power[Plus[1,1],ORDER],
      Power[Gamma[Plus[1,ORDER]],Plus[1,1]],
      Power[Gamma[Plus[1,Times[Plus[1,1],ORDER]]],Subtract[0,1]],
      Power[Plus[1,Cos[Times[Plus[1,1],Pi,Subtract[Divide[Subtract[N[Range[LENGTH]],1],LENGTH],Divide[1,Plus[1,1]]]]]],ORDER]
    ]
  ],
  CompilationTarget -> "C",
  RuntimeOptions -> "Speed"
] ;
ClearAll[S$COMPILED$FREQUENCY] ;
S$COMPILED$FREQUENCY::usage = "
S$COMPILED$FREQUENCY[FLAG,LENGTH,TOTAL,WINDOW,SAMPLE] -- estimate frequency for given complex flag <FLAG> (0/1),  window total <TOTAL> (real),window data <WINDOW> (list), signal length <LENGTH> = 2^N (integer), input signal <SAMPLE> = {...,SR_I,SI_I,...} (list) of length 2^(N+1)
" ;
S$COMPILED$FREQUENCY = Compile[
  {
    {FLAG,_Integer},
    {LENGTH,_Integer},
    {TOTAL,_Real},
    {WINDOW,_Real,1},
    {SAMPLE,_Real,1}
  },
  Block[
    {PI,LOCAL,FOURIER,FST,CND,FACTOR,MUL,COS$MUL,SIN$MUL,FREQUENCY},
    PI = Pi ;
    LOCAL = SAMPLE ;
    (* REMOVE (WEIGHTED) MEAN AND APPLY WINDOW *)
    LOCAL[[1;;2*LENGTH;;2]] = (LOCAL[[1;;2*LENGTH;;2]]-Total[WINDOW*LOCAL[[1;;2*LENGTH;;2]]]/TOTAL)*WINDOW ;
    LOCAL[[2;;2*LENGTH;;2]] = (LOCAL[[2;;2*LENGTH;;2]]-Total[WINDOW*LOCAL[[2;;2*LENGTH;;2]]]/TOTAL)*WINDOW ;
    (* FFT APPROXIMATION *)
    FOURIER = S$COMPILED$FFT[LENGTH,+1,LOCAL] ;
    MUL = FOURIER[[1;;2*LENGTH;;2]]^2+FOURIER[[2;;2*LENGTH;;2]]^2 ;
    FST = Part[Position[MUL,Max[MUL[[1;;Floor[LENGTH/(2-FLAG)]]]]],1,1] ;
    (* MODULATE SAMPLE *)
    FACTOR = 2.0*PI*(FST-2.0)/LENGTH ;
    MUL = FACTOR*(Range[LENGTH]-1.0) ;
    COS$MUL = Cos[MUL] ;
    SIN$MUL = Sin[MUL] ;
    FOURIER = LOCAL ;
    LOCAL[[1;;2*LENGTH;;2]] = FOURIER[[1;;2*LENGTH;;2]]*COS$MUL - FOURIER[[2;;2*LENGTH;;2]]*SIN$MUL ;
    LOCAL[[2;;2*LENGTH;;2]] = FOURIER[[1;;2*LENGTH;;2]]*SIN$MUL + FOURIER[[2;;2*LENGTH;;2]]*COS$MUL ;
    (* FFRFT APPROXIMATION *)
    LOCAL = S$COMPILED$FFRFT[LENGTH,2.0/LENGTH,LOCAL] ;
    MUL = Log10[10.^-16+Sqrt[LOCAL[[1;;2*LENGTH;;2]]^2+LOCAL[[2;;2*LENGTH;;2]]^2]] ;
    CND = Part[Position[MUL,Max[MUL]],1,1] ;
    (* PARABOLA APPROXIMATION *)
    FREQUENCY = CND-0.5+(MUL[[CND-1]]-MUL[[CND]])/(MUL[[CND-1]]-2*MUL[[CND]]+MUL[[CND+1]]) ;
    (* RESULT *)
    (FST-2.0+2.0*(FREQUENCY-1.0)/LENGTH)/LENGTH
  ],
  CompilationTarget -> "C",
  CompilationOptions -> {
    "InlineExternalDefinitions" -> True,
    "InlineCompiledFunctions" -> True
  },
  RuntimeOptions -> "Speed",
  Parallelization -> True,
  RuntimeAttributes -> {Listable}
] ;
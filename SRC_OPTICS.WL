(* ################################################################################################################################################################ *)
(* OPTICS (QUASIPERIODIC SIGNAL ANOMALY DETECTION, LINEAR OPTICS COMPUTATION), 2020-2023 *)
(* ################################################################################################################################################################ *)
ClearSystemCache[] ;                (* -- CLEAR SYSTEM CACHE *)
Set[$HistoryLength,0] ;             (* -- NO HISTORY *)
Get["CompiledFunctionTools`"] ;     (* -- LOAD COMPILED FUNCTION PACKAGE *)
Get["CCompilerDriver`"] ;           (* -- LOAD COMPILER DRIVER PACKAGE *)
(* ################################################################################################################################################################ *)
(* MODEL (LIST) *)
(* ################################################################################################################################################################ *)
ClearAll[O$MODEL$LIST] ;
O$MODEL$LIST::usage = "
O$MODEL$LIST -- list of column names to use in O$MODEL$LOAD[] function
" ;
O$MODEL$LIST := List["TYPE", "NAME", "S", "BX", "AX", "FX", "BY", "AY", "FY"] ;
(* ################################################################################################################################################################ *)
(* MODEL (RIZE) *)
(* ################################################################################################################################################################ *)
ClearAll[O$MODEL$RIZE] ;
O$MODEL$RIZE::usage = "
O$MODEL$RIZE -- first BPM name (string)
" ;
O$MODEL$RIZE := "STP0" ;
(* ################################################################################################################################################################ *)
(* MODEL (LOAD) *)
(* ################################################################################################################################################################ *)
ClearAll[O$MODEL$LOAD] ;
O$MODEL$LOAD::usage = "
O$MODEL$LOAD[\"UNCOUPLED\", FILE] -- load linear uncoupled model from file <FILE> (string)
" ;
Options[O$MODEL$LOAD] = List[
  Rule["NAME", O$MODEL$LIST],       (* -- LIST OF COLUMN NAMES (O$MODEL$LIST GLOBAL VARIABLE OR LIST OF STRINGS)*)
  Rule["REST", True],               (* -- DROP FIRST RECORD FLAG (LOGICAL) *)
  Rule["FULL", True],               (* -- SAVE FULL DATA FLAG (LOGICAL) *)
  Rule["TYPE", "MARK"],             (* -- BPM TYPE (STRING) *)
  Rule["RISE", O$MODEL$RIZE],       (* -- FIRST BPM NAME (STRING) *)
  Rule["LIST", List[]]              (* -- LIST OF BPM NAMES TO REMOVE (LIST OF STRINGS) *)
] ;
O$MODEL$LOAD[                       (* -- LOAD (UNCOUPLED) LINEAR MODEL *)
  "UNCOUPLED",                      (* -- MODEL TYPE (STRING) *)
  FILE_String,                      (* -- FILE (STRING) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Catch[
  Block[
    {TWISS, TABLE, UNIQUE, LENGTH, QX, QY, FX, FY, FIRST},
    If[SameQ[OptionValue["NAME"], List[]], Throw[$Failed]] ;
    If[Not[AllTrue[OptionValue["NAME"], StringQ]], Throw[$Failed]] ;
    If[Not[BooleanQ[OptionValue["REST"]]], Throw[$Failed]] ;
    If[Not[BooleanQ[OptionValue["FULL"]]], Throw[$Failed]] ;
    If[Not[StringQ[OptionValue["TYPE"]]], Throw[$Failed]] ;
    If[Not[StringQ[OptionValue["RISE"]]], Throw[$Failed]] ;
    If[Not[AllTrue[OptionValue["LIST"], StringQ]], Throw[$Failed]] ;
    If[Not[FileExistsQ[FILE]], Throw[$Failed]] ;
    TWISS = N[Import[FILE, "Table"]] ;
    TWISS = If[OptionValue["REST"], Rest[TWISS], TWISS] ;
    If[UnsameQ[Length[First[TWISS]], Length[OptionValue["NAME"]]], Throw[$Failed]] ;
    TABLE = If[OptionValue["FULL"], AssociationThread[Rule[OptionValue["NAME"], Transpose[TWISS]]], Nothing] ; 
    TWISS = Cases[TWISS, List[OptionValue["TYPE"], __]] ;
    If[SameQ[TWISS, List[]], Throw[$Failed]] ;
    TWISS = Map[Rest, TWISS] ;
    If[UnsameQ[First[First[TWISS]], First[Last[TWISS]]], Throw[$Failed]] ;
    If[Not[DuplicateFreeQ[Map[First, Rest[TWISS]]]], Throw[$Failed]] ;
    TWISS = DeleteCases[TWISS, List[NAME_, __] /; MemberQ[OptionValue["LIST"], NAME]] ;
    If[SameQ[TWISS, List[]], Throw[$Failed]] ;
    UNIQUE = Length[Rest[TWISS]] ;
    TWISS = Association[Thread[Rule[Rest[OptionValue["NAME"]], Transpose[TWISS]]]] ;
    LENGTH = Last[TWISS["S"]] ;
    QX = Divide[Last[TWISS["FX"]], Times[Plus[1, 1], Pi]] ;
    QY = Divide[Last[TWISS["FY"]], Times[Plus[1, 1], Pi]] ;
    FX = Differences[TWISS["FX"]] ;
    FY = Differences[TWISS["FY"]] ;
    FIRST = OptionValue["RISE"] ;
    If[Not[MemberQ[TWISS["NAME"], FIRST]], Throw[$Failed]] ;
    FIRST = Subtract[First[Flatten[Position[TWISS["NAME"], FIRST]]], 1] ;
    TWISS = Table[
      Association[
        Rule["THIS", TWISS[["NAME", INDEX]]],
        Rule["NEXT", TWISS[["NAME", INDEX + 1]]],
        Rule["S", TWISS[["S", INDEX]]],
        Rule["BX", TWISS[["BX", INDEX]]],
        Rule["BY", TWISS[["BY", INDEX]]],
		    Rule["AX", TWISS[["AX", INDEX]]],
		    Rule["AY", TWISS[["AY", INDEX]]],
		    Rule["FX", FX[[INDEX]]],
		    Rule["FY", FY[[INDEX]]],
		    Rule["FX.DEG", Divide[Times[180, FX[[INDEX]]], Pi]],
		    Rule["FY.DEG", Divide[Times[180, FY[[INDEX]]], Pi]],
		    Rule["FX.COT", Cot[FX[[INDEX]]]],
		    Rule["FY.COT", Cot[FY[[INDEX]]]]
      ],
      {INDEX, 1, UNIQUE}
    ] ;
    TWISS = RotateLeft[TWISS, FIRST] ;
    TWISS[[All, "S"]] = Mod[TWISS[[All, "S"]] + LENGTH - TWISS[[1, "S"]], LENGTH] ;
    TWISS = Association[Thread[Rule[Range[UNIQUE], TWISS]]] ;
    Association[
      Rule["UNIQUE", UNIQUE],
      Rule["LENGTH", LENGTH],
      Rule["QX", QX],
      Rule["QY", QY],
      Rule["NX", Abs[Mod[QX, 1, -1/2]]],
      Rule["NY", Abs[Mod[QY, 1, -1/2]]],
      Rule["NAME", TWISS[[All, "THIS"]]],
      Rule["ID", Association[KeyValueMap[Composition[Reverse, Rule], TWISS[[All, "THIS"]]]]],
      Rule["DATA", TWISS],
      Rule["TABLE", TABLE]
    ]
  ]
] ;
(* ################################################################################################################################################################ *)
(* MODEL (NAME) *)
(* ################################################################################################################################################################ *)
ClearAll[O$MODEL$GET$NAME] ;
O$MODEL$GET$NAME::usage = "
O$MODEL$GET$NAME[ID, MODEL] -- get BPM name (string) for given id <ID> (integer) and model data <MODEL> (association)
" ;
O$MODEL$GET$NAME[                   (* -- GET BPM NAME FOR GIVEN ID *)
  ID_Integer,                       (* -- BPM ID (INTEGER) *)
  MODEL_Association                 (* -- MODEL (ASSOCIATION, O$MODEL$LOAD[] OUTPUT) *)
] := MODEL["NAME", ID] ;
(* ################################################################################################################################################################ *)
(* MODEL (ID) *)
(* ################################################################################################################################################################ *)
ClearAll[O$MODEL$GET$ID] ;
O$MODEL$GET$ID::usage = "
O$MODEL$GET$ID[ID, MODEL] -- get BPM id (integer) for given name <NAME> (string) and model <MODEL> (association)
" ;
O$MODEL$GET$ID[                     (* -- GET BPM ID FOR GIVEN NAME *)
  NAME_String,                      (* -- BPM NAME (STRING) *)
  MODEL_Association                 (* -- MODEL (ASSOCIATION, O$MODEL$LOAD[] OUTPUT) *)
] := MODEL["ID", NAME] ;
(* ################################################################################################################################################################ *)
(* MODEL (TWISS) *)
(* ################################################################################################################################################################ *)
ClearAll[O$MODEL$GET$TWISS] ;
O$MODEL$GET$TWISS::usage = "
O$MODEL$GET$TWISS[ID, MODEL] -- get twiss data for given BPM id <ID> (integer) and model <MODEL> (association)
O$MODEL$GET$TWISS[NAME, MODEL] -- get twiss data for given BPM name <NAME> (string) and model <MODEL> (association)
" ;
O$MODEL$GET$TWISS[                  (* -- GET TWISS BY ID *)
  ID_Integer,                       (* -- BPM ID (INTEGER) *)
  MODEL_Association                 (* -- MODEL (ASSOCIATION, O$MODEL$LOAD[] OUTPUT) *)
] := MODEL["DATA", ID] ;
O$MODEL$GET$TWISS[                  (* -- GET TWISS BY NAME *)
  NAME_String,                      (* -- BPM NAME (STRING) *)
  MODEL_Association                 (* -- MODEL (ASSOCIATION, O$MODEL$LOAD[] OUTPUT) *)
] := O$MODEL$GET$TWISS[O$MODEL$GET$ID[NAME, MODEL], MODEL] ;
(* ################################################################################################################################################################ *)
(* MODEL (TRANSPORT) *)
(* ################################################################################################################################################################ *)
ClearAll[O$MODEL$TRANSPORT] ;
O$MODEL$TRANSPORT::usage = "
O$MODEL$TRANSPORT[\"UNCOUPLED\", {{BX1, AX1}, {BX2, AX2}, FX12}, {{BY1, AY1}, {BY2, AY2}, FY12}] -- generate uncoupled 4x4 transport matrix for given twiss parameters between two positions
" ;
O$MODEL$TRANSPORT[                  (* -- GENERATE UNCOUPLED TRANSPORT MATRIX *)
  "UNCOUPLED",                      (* -- TRANSPORT TYPE (STRING) *)
  {                                 (* -- HORIZONTAL PLANE *)
    {BX1_, AX1_},                   (* -- BETA AND ALPHA AT ONE *)
    {BX2_, AX2_},                   (* -- BETA AND ALPHA AT TWO *)
    FX12_                           (* -- PHASE ADVANCE BETWEEN ONE AND TWO *)
  },
  {                                 (* -- VERTICAL PLANE *)
    {BY1_, AY1_},                   (* -- BETA AND ALPHA AT ONE *)
    {BY2_, AY2_},                   (* -- BETA AND ALPHA AT TWO *)
    FY12_                           (* -- PHASE ADVANCE BETWEEN ONE AND TWO *)
  }
] := {
  {Sqrt[BX2/BX1]*(Cos[FX12] + AX1*Sin[FX12]), Sqrt[BX1*BX2]*Sin[FX12], 0, 0},
  {-(1 + AX1*AX2)/Sqrt[BX1*BX2]*Sin[FX12] + (AX1 - AX2)/Sqrt[BX1*BX2]*Cos[FX12], Sqrt[BX1/BX2]*(Cos[FX12] - AX2*Sin[FX12]), 0, 0},
  {0, 0, Sqrt[BY2/BY1]*(Cos[FY12] + AY1*Sin[FY12]), Sqrt[BY1*BY2]*Sin[FY12]},
  {0, 0, -(1 + AY1*AY2)/Sqrt[BY1*BY2]*Sin[FY12] + (AY1 - AY2)/Sqrt[BY1*BY2]*Cos[FY12], Sqrt[BY1/BY2]*(Cos[FY12] - AY2*Sin[FY12])}
} ;
(* ################################################################################################################################################################ *)
(* MODEL (MAKE) *)
(* ################################################################################################################################################################ *)
ClearAll[O$MODEL$MAKE] ;
O$MODEL$MAKE::usage = "
O$MODEL$MAKE[\"UNCOUPLED\", MODEL] -- generate uncoupled 4x4 transport matrices and maps between BPMs, defines globals (UNIQUE, TWISS, MATRIX and MAP)
" ;
Options[O$MODEL$MAKE] = List[
  Rule["TWISS", True]               (* -- FLAG TO (RE)COMPUTE TWISS (LOGICAL) *)
] ;
O$MODEL$MAKE[                       (* -- GENERATE TRANSPORT MATRICES AND MAPS FOR GIVEN MODEL (GLOBAL: UNIQUE, TWISS, MATRIX, MAP) *)
  "UNCOUPLED",                      (* -- MODEL TYPE (STRING) *)
  MODEL_Association,                (* -- MODEL (ASSOCIATION, O$MODEL$LOAD[] OUTPUT) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Catch[
  Block[
    {FLAG, TABLE, SHIFT},
    If[SameQ[MODEL, $Failed], Throw[$Failed]] ;
    ClearAll[UNIQUE] ;
    UNIQUE = MODEL["UNIQUE"] ;
    ClearAll[TWISS] ;
    TWISS = MODEL["DATA"] ;
    ClearAll[MATRIX] ;
    MATRIX[Q_Integer, P_Integer] := (MATRIX[Q, P] = Block[
	    {ID1, ID2, BX1, AX1, BY1, AY1, BX2, AX2, BY2, AY2, FX12, FY12},
	    ID1 = Mod[Q, UNIQUE, 1] ;
	    ID2 = Mod[P, UNIQUE, 1] ;
	    BX1 = TWISS[ID1, "BX"] ;
	    AX1 = TWISS[ID1, "AX"] ;
	    BY1 = TWISS[ID1, "BY"] ;
	    AY1 = TWISS[ID1, "AY"] ;
	    BX2 = TWISS[ID2, "BX"] ;
	    AX2 = TWISS[ID2, "AX"] ;
	    BY2 = TWISS[ID2, "BY"] ;
	    AY2 = TWISS[ID2, "AY"] ;
	    FX12 = TWISS[ID1, "FX"] ;
	    FY12 = TWISS[ID1, "FY"] ;
	  Developer`ToPackedArray[N[O$MODEL$TRANSPORT["UNCOUPLED", {{BX1, AX1}, {BX2, AX2}, FX12}, {{BY1, AY1}, {BY2, AY2}, FY12}]]]
    ]) /; Q + 1 == P ;
    MATRIX[Q_Integer, P_Integer] := (MATRIX[Q, P] = Dot[MATRIX[Q + 1, P], MATRIX[Q, Q + 1]]) /; Q < P ;
    MATRIX[Q_Integer, P_Integer] := (MATRIX[Q, P] = MATRIX[Q, P + UNIQUE]) /; Q == P ;
    FLAG = MATRIX[1, 1] ;
    FLAG = Chop[Norm[Re[Eigenvalues[MatrixLog[FLAG]]]]] ;
    If[UnsameQ[FLAG, 0], Throw[$Failed]] ;
    ClearAll[MAP] ;
    MAP[Q_Integer, P_Integer][List[Q1_, P1_, Q2_, P2_]] := Dot[MATRIX[Q, P], List[Q1, P1, Q2, P2]] ;
    MAP[UNIQUE, 1][List[Q1_, P1_, Q2_, P2_]] := Dot[MATRIX[UNIQUE, UNIQUE + 1], List[Q1, P1, Q2, P2]] ;
    MAP[Q_Integer, P_Integer] := Apply[Composition, Flatten[Reverse[RotateLeft[Table[List[MAP[ID, Mod[ID + 1, UNIQUE, 1]], Sow], {ID, 1, UNIQUE}], Q - 1]]]] /; Q == P ;
    If[
      OptionValue["TWISS"],
      TABLE = Table[S$CS$CONVERT[S$WOLSKI[MATRIX[ID, ID]]["WOLSKI"]], {ID, 1, UNIQUE}] ;
      SHIFT = Map[First, S$SHIFT[S$WOLSKI[MATRIX[1, 1]]["NORMAL"], Table[MATRIX[ID, ID + 1], {ID, 1, UNIQUE}]]] ;
      TWISS = Table[
        Association[
          Rule["THIS", TWISS[[ID, "THIS"]]],
          Rule["NEXT", TWISS[[ID, "NEXT"]]],
          Rule["S", TWISS[[ID, "S"]]],
          Rule["BX", TABLE[[ID, 1]]],
          Rule["BY", TABLE[[ID, 2]]],
          Rule["AX", TABLE[[ID, 3]]],
          Rule["AY", TABLE[[ID, 4]]],
          Rule["FX", SHIFT[[ID, 1]]],
          Rule["FY", SHIFT[[ID, 2]]],
		      Rule["FX.DEG", Divide[Times[180, SHIFT[[ID, 1]]], Pi]],
		      Rule["FY.DEG", Divide[Times[180, SHIFT[[ID, 2]]], Pi]],
		      Rule["FX.COT", Cot[SHIFT[[ID, 1]]]],
		      Rule["FY.COT", Cot[SHIFT[[ID, 2]]]]
        ],
        {ID, 1, UNIQUE}
      ] ;
      TWISS = Association[Thread[Rule[Range[UNIQUE], TWISS]]] ;
    ] ;
  ]
] ;
(* ################################################################################################################################################################ *)
(* MODEL (TEST) *)
(* ################################################################################################################################################################ *)
ClearAll[O$MODEL$TEST] ;
O$MODEL$TEST::usage = "
O$MODEL$TEST[\"UNCOUPLED\", MODEL] -- generate test uncoupled model with thin normal quadrupole errors between each BPM pair
" ;
Options[O$MODEL$TEST] = List[
  Rule["TWISS", True],              (* -- FLAG TO COMPUTE TWISS (LOGICAL) *)
  Rule["SEED", 1],                  (* -- RANDOM SEED *)
  Rule["SIGMA", 0.0001],            (* -- THIN QUADRUPOLE AMPLITUDE ERROR SIGMA (REAL) *)
  Rule["TRUNCATED", True],          (* -- FLAG TO GENERATE TRUNCATED MATRIX (LOGICAL) *)
  Rule["ORDER", 1]                  (* -- TRUNCATION ORDER (INTEGER) *)
] ;
O$MODEL$TEST[                       (* -- GENERATE UNCOUPLED TEST MODEL WITH THIN QUADRUPOLE ERRORS (GLOBAL: ERROR$LIST, ERROR$MATRIX, MATRIX$ERROR, MAP$ERROR, EPSILON, ERROR, MATRIX$TRUNCATED ) *)
  "UNCOUPLED",                      (* -- MODEL TYPE (STRING) *)
  MODEL_Association,                (* -- MODEL (ASSOCIATION, O$MODEL$LOAD[] OUTPUT) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Catch[
  Block[
    {SIGMA, FLAG, TABLE, SHIFT},
    If[SameQ[MODEL, $Failed], Throw[$Failed]] ;
    If[Not[IntegerQ[UNIQUE]], Throw[$Failed]] ;
    SeedRandom[OptionValue["SEED"]] ;
    SIGMA = OptionValue["SIGMA"] ;
    ClearAll[ERROR$LIST] ;
    ERROR$LIST = If[Greater[SIGMA, 0.0], RandomVariate[NormalDistribution[0.0, SIGMA], UNIQUE], ConstantArray[0.0, UNIQUE]] ;
    ClearAll[ERROR$MATRIX] ;
    ERROR$MATRIX = Table[{{1.0, 0.0, 0.0, 0.0}, {+ERROR$LIST[[ID]], 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, -ERROR$LIST[[ID]], 1.0}}, {ID, 1, UNIQUE}] ;
    ClearAll[MATRIX$ERROR] ;
    MATRIX$ERROR[Q_Integer, P_Integer] := (MATRIX$ERROR[Q, P] = Dot[ERROR$MATRIX[[Mod[Q, UNIQUE, 1]]], MATRIX[Q, P]]) /; Q + 1 == P ;
    MATRIX$ERROR[Q_Integer, P_Integer] := (MATRIX$ERROR[Q, P] = Dot[MATRIX$ERROR[Q + 1, P], MATRIX$ERROR[Q, Q + 1]]) /; Q < P ;
    MATRIX$ERROR[Q_Integer, P_Integer] := (MATRIX$ERROR[Q, P] = MATRIX$ERROR[Q, P + UNIQUE]) /; Q == P ;
    FLAG = MATRIX$ERROR[1, 1] ;
    FLAG = Chop[Norm[Re[Eigenvalues[MatrixLog[FLAG]]]]] ;
    If[UnsameQ[FLAG, 0], Throw[$Failed]] ;
    ClearAll[MAP$ERROR] ;
    MAP$ERROR[Q_Integer, P_Integer][List[Q1_, P1_, Q2_, P2_]] := Dot[MATRIX$ERROR[Q, P], List[Q1, P1, Q2, P2]] ;
    MAP$ERROR[UNIQUE, 1][List[Q1_, P1_, Q2_, P2_]] := Dot[MATRIX$ERROR[UNIQUE, UNIQUE + 1], List[Q1, P1, Q2, P2]] ;
    MAP$ERROR[Q_Integer, P_Integer] := Apply[Composition, Flatten[Reverse[RotateLeft[Table[List[MAP$ERROR[ID, Mod[ID + 1, UNIQUE, 1]], Sow], {ID, 1, UNIQUE}], Q - 1]]]] /; Q == P ;
    If[
      OptionValue["TWISS"],
      TABLE = Table[S$CS$CONVERT[S$WOLSKI[MATRIX$ERROR[ID, ID]]["WOLSKI"]], {ID, 1, UNIQUE}] ;
      SHIFT = Map[First, S$SHIFT[S$WOLSKI[MATRIX$ERROR[1, 1]]["NORMAL"], Table[MATRIX$ERROR[ID, ID + 1], {ID, 1, UNIQUE}]]] ;
      TWISS$ERROR = Table[
        Association[
          Rule["THIS", TWISS[[ID, "THIS"]]],
          Rule["NEXT", TWISS[[ID, "NEXT"]]],
          Rule["S", TWISS[[ID, "S"]]],
          Rule["BX", TABLE[[ID, 1]]],
          Rule["BY", TABLE[[ID, 2]]],
          Rule["AX", TABLE[[ID, 3]]],
          Rule["AY", TABLE[[ID, 4]]],
          Rule["FX", SHIFT[[ID, 1]]],
          Rule["FY", SHIFT[[ID, 2]]],
		      Rule["FX.DEG", Divide[Times[180, SHIFT[[ID, 1]]], Pi]],
		      Rule["FY.DEG", Divide[Times[180, SHIFT[[ID, 2]]], Pi]],
		      Rule["FX.COT", Cot[SHIFT[[ID, 1]]]],
		      Rule["FY.COT", Cot[SHIFT[[ID, 2]]]]
        ],
        {ID, 1, UNIQUE}
      ] ;
      TWISS$ERROR = Association[Thread[Rule[Range[UNIQUE], TWISS$ERROR]]] ;
    ] ;
    If[Not[OptionValue["TRUNCATED"]], Throw[Null]] ;
    ClearAll[EPSILON] ;
    With[
      {ORDER = OptionValue["ORDER"]},
      EPSILON /: Power[EPSILON, DEGREE_] /; DEGREE > ORDER := 0
    ] ;
    ClearAll[ERROR] ;
    ERROR[AMPLITUDE_] := {{1, 0, 0, 0}, {+AMPLITUDE*EPSILON, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, -AMPLITUDE*EPSILON, 1}} ;
    ClearAll[MATRIX$TRUNCATED] ;
    MATRIX$TRUNCATED[Q_Integer, P_Integer] := (MATRIX$TRUNCATED[Q, P] = Collect[Dot[ERROR[KNOB[Mod[Q, UNIQUE, 1]]], MATRIX[Q, P]], EPSILON, Expand]) /; Q + 1 == P ;
    MATRIX$TRUNCATED[Q_Integer, P_Integer] := (MATRIX$TRUNCATED[Q, P] = Collect[Dot[MATRIX$TRUNCATED[Q + 1, P], MATRIX$TRUNCATED[Q, Q + 1]], EPSILON, Expand]) /; Q < P ;
    MATRIX$TRUNCATED[Q_Integer, P_Integer] := (MATRIX$TRUNCATED[Q, P] = Collect[MATRIX$TRUNCATED[Q, P + UNIQUE], EPSILON, Expand]) /; Q == P ;
    SeedRandom[] ;
  ]
] ;
(* ################################################################################################################################################################ *)
(* TBT (START) *)
(* ################################################################################################################################################################ *)
ClearAll[O$TBT$START] ;
O$TBT$START::usage = "
O$TBT$START -- starting turn association for BPMs (<BPM> -> <START>), <BPM> -- BPM name (string), <START> -- starting (first) turn to use (integer, first index is one, inclusive)
" ;
O$TBT$START = Association[Rule["STP0", 10], Rule["STP2", 9], Rule["STP4", 9], Rule["SRP1", 9], Rule["SRP2", 9], Rule["SRP3", 9], Rule["SRP4", 9], Rule["SRP5", 9], Rule["SRP6", 9], Rule["SRP7", 9], Rule["SRP8", 9], Rule["SRP9", 9], Rule["SIP1", 9], Rule["SIP2", 9], Rule["SRP10", 9], Rule["SRP11", 9], Rule["SRP12", 9], Rule["SRP13", 9], Rule["SRP14", 9], Rule["SRP15", 9], Rule["SRP16", 9], Rule["SRP17", 9], Rule["SEP5", 10], Rule["SEP4", 10], Rule["SEP3", 10], Rule["SEP1", 10], Rule["SEP0", 10], Rule["NEP0", 10], Rule["NEP1", 10], Rule["NEP3", 10], Rule["NEP4", 10], Rule["NEP5", 10], Rule["NRP17", 9], Rule["NRP16", 9], Rule["NRP15", 9], Rule["NRP14", 9], Rule["NRP13", 10], Rule["NRP12", 10], Rule["NRP11", 10], Rule["NRP10", 10], Rule["NIP3", 10], Rule["NIP1", 10], Rule["NRP9", 10], Rule["NRP8", 10] ,Rule["NRP7", 10], Rule["NRP6", 10], Rule["NRP5", 10], Rule["NRP4", 10], Rule["NRP3", 10], Rule["NRP2", 10], Rule["NRP1", 10], Rule["NTP4", 10], Rule["NTP2", 10], Rule["NTP0", 10]] ;
(* ################################################################################################################################################################ *)
(* TBT (LOAD) *)
(* ################################################################################################################################################################ *)
ClearAll[O$TBT$LOAD] ;
O$TBT$LOAD::usage = "
O$TBT$LOAD[MODEL, FILE] -- load TbT data archive with file name prexix <FILE> (string) (full filename should be \"<FILE>.tar.gz\"), data archive is assumed to contain files with individual BPM signals for each BPM present in given model <MODEL> (O$MODEL$MAKE[] function output), each BPM filename is of the form \"<FILE>-<NAME>.dat\", <NAME> (string) -- BPM name (lowcase), each BPM file is assumed to be in \"CSV\" format and contain three columns (X, Y and I)
" ;
Options[O$TBT$LOAD] = List[
  Rule["FIRST", O$TBT$START],       (* -- LIST OF FIRST TURNS FOR EACH BPM (ASSOCIATION) *) 
  Rule["SIZE", 1024],               (* -- MAXIMUN NUMBER OF TURNS TO KEEP (INTEGER) *)
  Rule["TEST", 512]                 (* -- SAMPLE TEST SIZE (INTEGER) *)
] ;
O$TBT$LOAD[                         (* -- LOAD TBT *)
  MODEL_Association,                (* -- MODEL (ASSOCIATION, O$MODEL$MAKE[] OUTPUT) *)
  FILE_String,                      (* -- DATA ARCHIVE PREXIF (STRING) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Catch[
  Block[
    {NAME, RIZE, SIZE, TEST, DATA},
    If[Not[FileExistsQ[StringTemplate["`1`.tar.gz"][FILE]]], Throw[$Failed]] ;
    If[UnsameQ[Run["tar --version"], 0], Throw[$Failed]] ;
    Run[StringTemplate["tar -xzvf `1`.tar.gz"][FILE]] ;
    NAME = Values[MODEL["NAME"]] ;
    RIZE = OptionValue["FIRST"] ;
    RIZE = Map[RIZE, NAME] ;
    SIZE = OptionValue["SIZE"] ;
    TEST = OptionValue["TEST"] ;
    DATA = ParallelTable[
      Block[
        {LOAD, FLAG},
        LOAD = StringTemplate["`1`-`2`.dat"][FILE, ToLowerCase[NAME[[INDEX]]]] ;
        FLAG = FileExistsQ[LOAD] ;
        LOAD = If[FLAG, N[Import[LOAD, "CSV"]], ConstantArray[0.0, List[SIZE, Plus[1, 1, 1]]]] ;
        LOAD = If[UnsameQ[Length[LOAD], 1], Take[LOAD, List[RIZE[[INDEX]], RIZE[[INDEX]] + SIZE - 1]], ConstantArray[0.0, List[SIZE, Plus[1, 1, 1]]]] ;
        Transpose[LOAD]
      ],
      {INDEX, 1, Length[NAME]}
    ] ;
    DATA = Association[Thread[Rule[List["X", "Y", "Z"], Transpose[DATA]]]] ;
    Run[StringTemplate["rm `1`*.dat"][FILE]] ;
    Association[
      Rule["SOURCE", "LOAD"],
      Rule["ARCHIVE", FILE],
      Rule["SIZE", SIZE],
      Rule["TEST", TEST],
      Rule["DATA", DATA]
    ]
  ]
] ;
(* ################################################################################################################################################################ *)
(* TBT (TEST) *)
(* ################################################################################################################################################################ *)
ClearAll[O$TBT$TEST] ;
O$TBT$TEST::usage = "
O$TBT$TEST[MODEL, INITIAL] -- generate test TbT data for given model <MODEL> (O$MODEL$MAKE[] function output) and initial condition <INITIAL> (list of reals, should match MAP)
" ;
Options[O$TBT$TEST] = List[
  Rule["SIZE", 1024],               (* -- NUMBER OF TURNS TO GENERATE (INTEGER) *)
  Rule["TEST", 256],                (* -- SAMPLE TEST SIZE (INTEGER) *)
  Rule["CURRENT", List[1.0, 0.01]], (* -- DEFAULT CURRENT MEAN AND SIGMA VALUES (REAL) *)
  Rule["FLAG", False],              (* -- FLAG TO USE ENSEMBLE OF INITIAL CONDITIONS (LOGICAL) *)
  Rule["BEAM", 10^3],               (* -- TOTAL NUMBER OF INITIAL CONDITIONS (INTEGER) *)
  Rule["SQX", 0.5*10.0^-3],         (* -- H-SIZE (REAL) *)
  Rule["SPX", 4.0*10.0^-5],         (* -- H-SIZE (REAL) *)
  Rule["SQY", 0.1*10.0^-3],         (* -- V-SIZE (REAL) *)
  Rule["SPY", 2.0*10.0^-5],         (* -- H-SIZE (REAL) *)
  Rule["NONLINEAR", False],         (* -- FLAG TO USE THIN NONLINEAR KICK (LOGICAL) *)
  Rule["KN", 0.5],                  (* -- N-AMPLITUDE (REAL) *)
  Rule["KS", 0.5],                  (* -- S-AMPLITUDE (REAL) *)
  Rule["CHROMATIC", True],          (* -- FLAG TO USE CHROMATIC MATRIX (LOGICAL) *)
  Rule["DELTA", 10.0^-4],           (* -- ENERGY SPREAD (REAL) *)
  Rule["CX", 0.5],                  (* -- H-CHROMATIC FACTOR (REAL) *)
  Rule["CY", 0.5]                   (* -- V-CHROMATIC FACTOR (REAL) *)
] ;
O$TBT$TEST[                         (* -- GENERATE TBT DATA *)
  MODEL_Association,                (* -- MODEL (ASSOCIATION, O$MODEL$MAKE[] OUTPUT) *)
  INITIAL_List,                     (* -- INITIAL CONDITION AT BPM WITH INDEX ONE (LIST OF REALS, SHOULD MATCH 'MAP' FUNCTION DEFINED BY O$MODEL$MAKE[]) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Catch[
  Block[
    {SIZE, TEST, FLAG, DATA, X, Y, Z, BEAM, NONLINEAR, KN, KS, QX, PX, QY, PY, CHROMATIC, CX, CY, AX, BX, AY, BY, DELTA},
    SIZE = OptionValue["SIZE"] ;
    TEST = OptionValue["TEST"] ;
    FLAG = OptionValue["FLAG"] ;
    If[
      Not[FLAG],
      DATA = First[Last[Reap[NestList[MAP$ERROR[1, 1], INITIAL, SIZE]]]] ;
      DATA = DATA[[All, List[1, Plus[1, 1, 1]]]] ;
      DATA = Transpose[Partition[DATA, MODEL[["UNIQUE"]]]] ;
      DATA = Transpose[Map[Transpose, DATA]] ;
      {X, Y} = DATA ;
      Z = RandomVariate[Apply[NormalDistribution, OptionValue["CURRENT"]], List[MODEL[["UNIQUE"]], OptionValue["SIZE"]]] ;
      DATA = Association[Thread[Rule[List["X", "Y", "Z"], List[X, Y, Z]]]] ;
      DATA = Association[
        Rule["ORIGIN", "TEST"],
        Rule["INITIAL", INITIAL],
        Rule["SIZE", SIZE],
        Rule["TEST", TEST],
        Rule["DATA", DATA]
      ] ;
      Throw[DATA] ;
    ] ;
    BEAM = OptionValue["BEAM"] ;
    BEAM = RandomVariate[
      MultinormalDistribution[
        INITIAL,
        DiagonalMatrix[List[OptionValue["SQX"], OptionValue["SPX"], OptionValue["SQY"], OptionValue["SPY"]]^2]
      ],
      BEAM
    ] ;
    NONLINEAR = If[
      OptionValue["NONLINEAR"],
      KN = OptionValue["KN"] ;
      KS = OptionValue["KS"] ;
      Apply[Function, List[List[QX, PX, QY, PY], List[QX, PX + 2*KS*QX*QY - KN*(QX^2 - QY^2), QY, PY + 2*KN*QX*QY + KS*(QX^2 - QY^2)]]],
      Apply[Function, List[List[QX, PX, QY, PY], List[QX, PX, QY, PY]]]
    ] ;
    CHROMATIC = If[
      OptionValue["CHROMATIC"],
      CX = OptionValue["CX"] ;
      CY = OptionValue["CY"] ;
      AX = TWISS$ERROR[[1, "AX"]] ;
      BX = TWISS$ERROR[[1, "BX"]] ;
      AY = TWISS$ERROR[[1, "AY"]] ;
      BY = TWISS$ERROR[[1, "BY"]] ;
      Apply[Function, List[List[QX, PX, QY, PY, DELTA], List[QX*(Cos[2*Pi*CX*DELTA] + AX*Sin[2*Pi*CX*DELTA]) + BX*PX*Sin[2*Pi*CX*DELTA], -1/BX*((1 + AX^2)*QX*Sin[2*Pi*CX*DELTA]) + PX*(Cos[2*Pi*CX*DELTA] - AX*Sin[2*Pi*CX*DELTA]), QY*(Cos[2*Pi*CY*DELTA] + AY*Sin[2*Pi*CY*DELTA]) + BY*PY*Sin[2*Pi*CY*DELTA], -1/BY*((1 + AY^2)*QY*Sin[2*Pi*CY*DELTA]) + PY*(Cos[2*Pi*CY*DELTA] - AY*Sin[2*Pi*CY*DELTA]), DELTA]]],
      Apply[Function, List[List[QX, PX, QY, PY, DELTA], List[QX, PX, QY, PY, DELTA]]]
    ] ;
    BEAM = If[
      OptionValue["CHROMATIC"],
      Map[Flatten, Transpose[List[BEAM, RandomVariate[NormalDistribution[N[0], OptionValue["DELTA"]], Length[BEAM]]]]],
      Map[Flatten, Transpose[List[BEAM, ConstantArray[N[0], Length[BEAM]]]]]
    ] ;
    BEAM = ParallelMap[
      Function[
        Block[
          {KEEP, Q1, P1, Q2, P2, DP, LIST},
          {Q1, P1, Q2, P2, DP} = Slot[1] ;
          KEEP = ConstantArray[N[0], SIZE] ;
          Do[
            {{Q1, P1, Q2, P2}, LIST} = Reap[MAP$ERROR[1, 1][{Q1, P1, Q2, P2}]] ;
            {Q1, P1, Q2, P2, DP} = Apply[CHROMATIC, {Q1, P1, Q2, P2, DP}] ;
            {Q1, P1, Q2, P2} = Apply[NONLINEAR, {Q1, P1, Q2, P2}] ;
            KEEP[[ID]] = First[LIST],
            {ID, 1, SIZE}
          ] ;
          KEEP = Transpose[KEEP] ;
          KEEP[[All, All, List[1, 3]]]
        ]
      ],
      BEAM
    ] ;
    X = Table[Map[Mean, Transpose[BEAM[[All, INDEX, All, +1]]]], {INDEX, 1, MODEL[["UNIQUE"]]}] ;
    Y = Table[Map[Mean, Transpose[BEAM[[All, INDEX, All, -1]]]], {INDEX, 1, MODEL[["UNIQUE"]]}] ;
    Z = RandomVariate[Apply[NormalDistribution, OptionValue["CURRENT"]], List[MODEL[["UNIQUE"]], OptionValue["SIZE"]]] ;
    DATA = Association[Thread[Rule[List["X", "Y", "Z"], List[X, Y, Z]]]] ;
    Association[
      Rule["ORIGIN", "TEST"],
      Rule["INITIAL", INITIAL],
      Rule["SIZE", SIZE],
      Rule["TEST", TEST],
      Rule["DATA", DATA]
    ]
  ]
] ;
(* ################################################################################################################################################################ *)
(* TBT (ERROR) *)
(* ################################################################################################################################################################ *)
ClearAll[O$TBT$ERROR] ;
O$TBT$ERROR::usage = "
O$TBT$ERROR[MODEL, TBT, TILT, GAIN] -- add calibration errors (<TILT> and <GAIN>) to TbT data <TBT> (association) with model <MODEL> (association)
" ;
O$TBT$ERROR[                        (* -- GENERATE TBT DATA WITH CALIBRATION ERRORS *)
  MODEL_Association,                (* -- MODEL (ASSOCIATION, O$MODEL$MAKE[] OUTPUT) *)
  TBT_Association,                  (* -- TBT DATA (ASSOCIATION) *)
  TILT_List,                        (* -- TITL LIST (LIST) *)
  GAIN_List                         (* -- GAIN LIST (LIST) *)
] := Catch[
  Block[
    {COS, SIN, MIX, X, Y, OUT},
    If[UnsameQ[Dimensions[TILT], List[MODEL[["UNIQUE"]]]], Throw[$Failed]] ;
    If[UnsameQ[Dimensions[GAIN], List[MODEL[["UNIQUE"]], Plus[1, 1]]], Throw[$Failed]] ;
    COS = Cos[TILT] ;
    SIN = Sin[TILT] ;
    MIX[{C_, S_}, {GX_, GY_}][{X_, Y_}] := Dot[{{C, S}, {-S, C}}, Dot[{{GX, 0}, {0, GY}}, {X, Y}]] ;
    X = TBT[["DATA", "X"]] ;
    Y = TBT[["DATA", "Y"]] ;
    {X, Y} = Transpose[Table[MIX[List[COS[[INDEX]], SIN[[INDEX]]], GAIN[[INDEX]]][List[X[[INDEX]], Y[[INDEX]]]], {INDEX, 1, MODEL[["UNIQUE"]]}]] ;
    OUT = TBT ;
    OUT[["DATA", "X"]] = X ;
    OUT[["DATA", "Y"]] = Y ;
    Association[
      OUT,
      Rule["TILT", TILT],
      Rule["GAIN", GAIN]
    ]
  ]
] ;
(* ################################################################################################################################################################ *)
(* TBT (NOISE) *)
(* ################################################################################################################################################################ *)
ClearAll[O$TBT$NOISE] ;
O$TBT$NOISE::usage = "
O$TBT$NOISE[MODEL, TBT, SIGMA] -- add noise errors <LIST> to TbT data <TBT> (association) with model <MODEL> (association)
" ;
Options[O$TBT$NOISE] = List[
  Rule["SEED", 1]                   (* -- RANDOM SEED *)
] ;
O$TBT$NOISE[                        (* -- GENERATE TBT DATA WITH RANDOM NOISE *)
  MODEL_Association,                (* -- MODEL (ASSOCIATION, O$MODEL$MAKE[] OUTPUT) *)
  TBT_Association,                  (* -- TBT DATA (ASSOCIATION) *)
  LIST_List,                        (* -- NOISE SIGMA LIST *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Catch[
  Block[
    {SIZE, NOISE, X, Y, OUT},
    If[UnsameQ[Dimensions[LIST], List[MODEL[["UNIQUE"]], Plus[1, 1]]], Throw[$Failed]] ;
    SIZE = TBT["SIZE"] ;
    SeedRandom[OptionValue["SEED"]] ;
    NOISE = Flatten[LIST] ;
    NOISE = RandomVariate[MultinormalDistribution[ConstantArray[0.0, Times[Plus[1, 1], MODEL[["UNIQUE"]]]], DiagonalMatrix[NOISE*NOISE]], SIZE] ;
    NOISE = Transpose[Partition[Transpose[NOISE], Plus[1, 1]]] ;
    SeedRandom[] ;
    X = TBT[["DATA", "X"]] ;
    Y = TBT[["DATA", "Y"]] ;
    {X, Y} += NOISE ;
    OUT = TBT ;
    OUT[["DATA", "X"]] = X ;
    OUT[["DATA", "Y"]] = Y ;
    Association[
      OUT,
      Rule["NOISE", LIST]
    ]
  ]
] ;
(* ################################################################################################################################################################ *)
(* ANOMALY (NOISE) *)
(* ################################################################################################################################################################ *)
ClearAll[O$ANOMALY$NOISE] ;
O$ANOMALY$NOISE::usage = "
O$ANOMALY$NOISE[SIGMA, LENGTH] -- generate (addative) noise anomaly for given noise sigma <SIGMA> (real) and signal length <LENGTH> (integer)
" ;
O$ANOMALY$NOISE[                    (* -- GENERATE NOISE ANOMALY (ADDATIVE) *)
  SIGMA_Real,                       (* -- NOISE SIGMA (REAL), ABSOLUTE VALUE *)
  LENGTH_Integer                    (* -- SIGNAL LENGTH (INTEGER) *)
] := If[
  Greater[SIGMA, 0.0],
  RandomVariate[NormalDistribution[0.0, SIGMA], LENGTH],
  ConstantArray[0.0, LENGTH]
] ;
(* ################################################################################################################################################################ *)
(* ANOMALY (SPIKE) *)
(* ################################################################################################################################################################ *)
ClearAll[O$ANOMALY$SPIKE] ;
O$ANOMALY$SPIKE::usage = "
O$ANOMALY$SPIKE[COUNT, SIGMA, LENGTH] -- generate (addative) spike anomaly for given number of spikes <COUNT> (integer), spike amplitude sigma <SIGMA> (real) and signal length <LENGTH> (integer)
" ;
O$ANOMALY$SPIKE[                    (* -- GENERATE SPIKE ANOMALY (ADDATIVE) *)
  COUNT_Integer,                    (* -- NUMBER OF SPIKES (INTEGER) *)
  SIGMA_Real,                       (* -- SPIKE SIGMA (REAL) *)
  LENGTH_Integer                    (* -- SIGNAL LENGTH (INTEGER) *)
] := Block[
  {INDEX, VALUE, SPIKE},
  INDEX = RandomSample[Range[LENGTH], COUNT] ;
  VALUE = If[Greater[SIGMA, 0.0], RandomVariate[NormalDistribution[0.0, SIGMA], COUNT], ConstantArray[0.0, COUNT]] ;
  SPIKE = ConstantArray[0.0, LENGTH] ;
  SPIKE[[INDEX]] = VALUE ;
  SPIKE
] ;
(* ################################################################################################################################################################ *)
(* ANOMALY (SHIFT) *)
(* ################################################################################################################################################################ *)
ClearAll[O$ANOMALY$SHIFT] ;
O$ANOMALY$SHIFT::usage = "
O$ANOMALY$SHIFT[START, END, VALUE, LENGTH] -- generate (addative) shift anomaly for given anomaly start <START> (integer), end <END> (integer), shift value <VALUE> (real) and signal length <LENGTH> (integer)
" ;
O$ANOMALY$SHIFT[                    (* -- GENERATE SHIFT ANOMALY (ADDATIVE) *)
  START_Integer,                    (* -- SHIFT START (INTEGER) *)
  END_Integer,                      (* -- SHIFT END (INTEGER) *)
  VALUE_Real,                       (* -- SHIFT VALUE (REAL) *)
  LENGTH_Integer                    (* -- SIGNAL LENGTH (INTEGER) *)
] := Block[
  {INDEX, SHIFT},
  INDEX = Range[START, END] ;
  SHIFT = ConstantArray[0.0, LENGTH] ;
  SHIFT[[INDEX]] = VALUE ;
  SHIFT
] ;
(* ################################################################################################################################################################ *)
(* ANOMALY (SCALE) *)
(* ################################################################################################################################################################ *)
ClearAll[O$ANOMALY$SCALE] ;
O$ANOMALY$SCALE::usage = "
O$ANOMALY$SCALE[START, END, VALUE, LENGTH] -- generate (multiplicative) scale anomaly for given anomaly start <START> (integer), end <END> (integer), scale value <VALUE> (real) and signal length <LENGTH> (integer)
" ;
O$ANOMALY$SCALE[                    (* -- GENERATE SCALE ANOMALY (MULTIPLICATIVE) *)
  START_Integer,                    (* -- SHIFT START (INTEGER) *)
  END_Integer,                      (* -- SHIFT END (INTEGER) *)
  VALUE_Real,                       (* -- SCALE VALUE (REAL) *)
  LENGTH_Integer                    (* -- SIGNAL LENGTH (INTEGER) *)
] := Block[
  {INDEX, SCALE},
  INDEX = Range[START, END] ;
  SCALE = ConstantArray[1.0, LENGTH] ;
  SCALE[[INDEX]] = 1.0+VALUE ;
  SCALE
] ;
(* ################################################################################################################################################################ *)
(* TBT (ANOMALY) *)
(* ################################################################################################################################################################ *)
ClearAll[O$TBT$ANOMALY] ;
O$TBT$ANOMALY::usage = "
O$TBT$ANOMALY[MODEL, TBT] -- generate TbT data with anomalies for given data <TBT> (association) with model <MODEL> (association)
" ;
Options[O$TBT$ANOMALY] = List[
  Rule["PROBABILITY", 0.05],        (* -- ANOMALY PROBABILITY (REAL), CHANCE FOR EACH INDEX TO BE AN ANOMALY *)
  Rule["NOISE", {1.0, 5.0}*10^-3],  (* -- NOISE ANOMALY MAX SIGMA RANGE (REAL, REAL) *)
  Rule["SPIKE", {2^4, 0.01}],       (* -- SPIKE ANOMALY MAX NUMBER OF SPIKES AND SPIKE SIGMA (INTEGER, REAL) *)
  Rule["SHIFT", {2^6, 2^8, 0.001}], (* -- SHIFT ANOMALY MAX START, MIN END AND MAX ABS SHIFT VALUE (INTEGER, INTEGER, REAL) *)
  Rule["SCALE", {2^6, 2^8, 0.1}],   (* -- SCALE ANOMALY MAX START, MIN END AND MAX SCALE VALUE (INTEGER, INTEGER, REAL) *)
  Rule["RANDOM", True]              (* -- RANDOM FLAG (LOGICAL) *)
] ;
O$TBT$ANOMALY[                      (* -- GENERATE TBT DATA WITH ANOMALIES *)
  MODEL_Association,                (* -- MODEL (ASSOCIATION, O$MODEL$MAKE[] OUTPUT) *)
  TBT_Association,                  (* -- TBT DATA (ASSOCIATION) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Catch[
  Block[
    {INDEX, LENGTH, NOISE, COUNT, SIGMA, SPIKE, START, END, VALUE, SHIFT, SCALE, OUT, RANDOM},
    INDEX = RandomVariate[UniformDistribution[], MODEL[["UNIQUE"]]] ;
    INDEX = Flatten[Position[INDEX, PROBABILITY_ /; PROBABILITY <= OptionValue["PROBABILITY"]]] ;
    If[SameQ[INDEX, List[]], Throw[Association[Rule["ANOMALY", Null], TBT]]] ;
    LENGTH = TBT["SIZE"] ;
    NOISE = OptionValue["NOISE"] ;
    NOISE = RandomVariate[UniformDistribution[NOISE]] ;
    NOISE = Table[O$ANOMALY$NOISE[NOISE, LENGTH], Times[Plus[1, 1], Length[INDEX]]] ;
    NOISE = Partition[NOISE, Plus[1, 1]] ;
    {COUNT, SIGMA} = OptionValue["SPIKE"] ;
    SPIKE = Table[O$ANOMALY$SPIKE[RandomChoice[Range[COUNT]], SIGMA, LENGTH], Times[Plus[1, 1], Length[INDEX]]] ;
    SPIKE = Partition[SPIKE, Plus[1, 1]] ;
    {START, END, VALUE} = OptionValue["SHIFT"] ;
    VALUE = Abs[VALUE] ;
    SHIFT = Table[O$ANOMALY$SHIFT[RandomChoice[Range[1, START]], RandomChoice[Range[END, LENGTH]], RandomChoice[List[-1, 1]]*RandomVariate[UniformDistribution[List[0.5*VALUE, VALUE]]], LENGTH], Times[Plus[1, 1], Length[INDEX]]] ;
    SHIFT = Partition[SHIFT, Plus[1, 1]] ;
    {START, END, VALUE} = OptionValue["SCALE"] ;
    SCALE = Table[O$ANOMALY$SCALE[RandomChoice[Range[1, START]], RandomChoice[Range[END, LENGTH]], RandomVariate[UniformDistribution[List[0.5*VALUE, VALUE]]], LENGTH], Times[Plus[1, 1], Length[INDEX]]] ;
    SCALE = Partition[SCALE, Plus[1, 1]] ;
    OUT = TBT ;
    COUNT = 0 ;
    RANDOM := If[OptionValue["RANDOM"], RandomInteger[], 1] ;
    Do[
      COUNT++;
      OUT[["DATA", "X", THIS]] = (OUT[["DATA", "X", THIS]] + NOISE[[COUNT, 1]] + RANDOM*SPIKE[[COUNT, 1]] + RANDOM*SHIFT[[COUNT, 1]])*(RANDOM*(SCALE[[COUNT, 1]] - 1) + 1) ;
      OUT[["DATA", "Y", THIS]] = (OUT[["DATA", "Y", THIS]] + NOISE[[COUNT, Plus[1, 1]]] + RANDOM*SPIKE[[COUNT, Plus[1, 1]]] + RANDOM*SHIFT[[COUNT, Plus[1, 1]]])*(RANDOM*(SCALE[[COUNT, Plus[1, 1]]] - 1) + 1),
      {THIS, INDEX}
    ] ;
    Association[
      OUT,
      Rule["ANOMALY", INDEX]
    ]
  ]
] ;
(* ################################################################################################################################################################ *)
(* TBT (MEAN) *)
(* ################################################################################################################################################################ *)
ClearAll[O$TBT$MEAN] ;
O$TBT$MEAN::usage = "
O$TBT$MEAN[TBT] -- remove mean values for x and y signals for given TbT data <TBT> (association)
" ;
O$TBT$MEAN[                         (* -- REMOVE MEAN VALUES *)
  TBT_Association                   (* -- TBT DATA (ASSOCIATION) *)
] := Block[
  {OUT, X, Y, Z, MX, MY, MZ},
  OUT = TBT ;
  X = OUT[["DATA", "X"]] ;
  Y = OUT[["DATA", "Y"]] ;
  Z = OUT[["DATA", "Z"]] ;
  MX = Map[Mean, X] ;
  MY = Map[Mean, Y] ;
  MZ = Map[Mean, Z] ;
  OUT[["DATA", "X"]] = X - MX ;
  OUT[["DATA", "Y"]] = Y - MY ;
  OUT[["DATA", "Z"]] = Z - MZ ;
  Association[
    OUT,
    Rule["MEAN", List[MX, MY, MZ]]
  ]
] ;
(* ################################################################################################################################################################ *)
(* ESTIMATE NOISE (OPTIMAL SVD TRUNCATION) *)
(* ################################################################################################################################################################ *)
ClearAll[O$ESTIMATE$NOISE] ;
O$ESTIMATE$NOISE::usage = "
O$ESTIMATE$NOISE[TBT] -- estimate noise for given TbT data <TBT> (association)
" ;
O$ESTIMATE$NOISE[                   (* -- ESTIMATE NOISE *)
  TBT_Association                   (* -- TBT DATA (ASSOCIATION) *)
] := Block[
  {X, Y, NX, NY},
  X = TBT[["DATA", "X"]] ;
  Y = TBT[["DATA", "Y"]] ;
  NX = Map[Composition[Last, S$OPTIMAL$SVD$TRUNCATION, S$MAKE$MATRIX], X] ;
  NY = Map[Composition[Last, S$OPTIMAL$SVD$TRUNCATION, S$MAKE$MATRIX], Y] ;
  Transpose[List[NX, NY]]
] ;
(* ################################################################################################################################################################ *)
(* GENERATOR (AMPLITUDE) *)
(* ################################################################################################################################################################ *)
ClearAll[O$GENERATOR$AMPLITUDE] ;
O$GENERATOR$AMPLITUDE::usage = "
O$GENERATOR$AMPLITUDE[] -- (feature generator) maximum absolute amplitude feature (acts on a single sample, recommended transformation 'Standardize')
" ;
O$GENERATOR$AMPLITUDE[              (* -- (FEATURE GENERATOR) MAXIMUM ABSOLUTE AMPLITUDE *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Composition[
  Max,
  Abs
] ;
(* ################################################################################################################################################################ *)
(* GENERATOR (SPECTRA) *)
(* ################################################################################################################################################################ *)
ClearAll[O$GENERATOR$SPECTRA] ;
O$GENERATOR$SPECTRA::usage = "
O$GENERATOR$SPECTRA[RANGE] -- (feature generator) fourier spectra mean in given range <RANGE> (min and max) of frequencies (acts on a single sample, window length should match sample length, recommended transformation 'Standardize')
" ;
O$GENERATOR$SPECTRA[                (* -- (FEATURE GENERATOR) FOURIER SPECTRA MEAN IN GIVEN RANGE OF FREQUENCIES *)
  RANGE_,                           (* -- RANGE OF FREQUENCIES, {MIN, MAX} *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Composition[
  Mean,
  Last,
  Transpose,
  Select[Composition[Between[RANGE], First]],
  Function[S$SPECTRA[Times[S$WINDOW$DATA, Slot[1]], FilterRules[List[OPTIONS], Options[S$SPECTRA]], Rule["NORMALIZED", True], Rule["2D", True]]]
] ;
(* ################################################################################################################################################################ *)
(* GENERATOR (FREQUENCY) *)
(* ################################################################################################################################################################ *)
ClearAll[O$GENERATOR$FREQUENCY] ;
O$GENERATOR$FREQUENCY::usage = "
O$GENERATOR$FREQUENCY[] -- (feature generator) sample frequency (acts on a single sample, window length should match sample length, recommended transformation 'Standardize')
" ;
O$GENERATOR$FREQUENCY[              (* -- (FEATURE GENERATOR) SAMPLE FREQUENCY *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Composition[
  Function[S$FREQUENCY[Slot[1], FilterRules[List[OPTIONS], Options[S$FREQUENCY]]]]
] ;
(* ################################################################################################################################################################ *)
(* GENERATOR (DECOMPOSITION) *)
(* ################################################################################################################################################################ *)
ClearAll[O$GENERATOR$DECOMPOSITION] ;
O$GENERATOR$DECOMPOSITION::usage = "
O$GENERATOR$DECOMPOSITION[NUMBER] -- (feature generator) decomposition error for given number of harmonics <NUMBER> (integer) (acts on a single sample, window length should match sample length, recommended transformation 'Standardize')
" ;
O$GENERATOR$DECOMPOSITION[          (* -- (FEATURE GENERATOR) DECOMPOSITION ERROR *)
  NUMBER_,                          (* -- NUMBER OF HARMONICS (INTEGER) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Composition[
  Norm,
  Function[Subtract[Slot[1], S$DECOMPOSITION[NUMBER, Slot[1], FilterRules[List[OPTIONS], Options[S$DECOMPOSITION]]]["FUNCTION"][N[Range[S$WINDOW$LENGTH]]]]]
] ;
(* ################################################################################################################################################################ *)
(* GENERATOR (SVD) *)
(* ################################################################################################################################################################ *)
ClearAll[O$GENERATOR$SVD] ;
O$GENERATOR$SVD::usage = "
O$GENERATOR$SVD[] -- (feature generator) svd variation (acts on full sample matrix, recommended transformation 'Standardize')
" ;
O$GENERATOR$SVD[                    (* -- (FEATURE GENERATOR) SVD VARIATION *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Composition[
  Map[Max],
  Abs,
  First,
  SingularValueDecomposition[Slot[1], Min[Dimensions[Slot[1]]], FilterRules[List[OPTIONS], Options[SingularValueDecomposition]]]
] ;
(* ################################################################################################################################################################ *)
(* GENERATOR (RANK) *)
(* ################################################################################################################################################################ *)
ClearAll[O$GENERATOR$RANK] ;
O$GENERATOR$RANK::usage = "
O$GENERATOR$RANK[] -- (feature generator) rank estimation (acts on a single sample, recommended transformation 'Standardize' and non-robust dispersion estimator)
" ;
O$GENERATOR$RANK[                   (* -- (FEATURE GENERATOR) RANK ESTIMATION *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Composition[
  N,
  First,
  S$OPTIMAL$SVD$TRUNCATION,
  S$MAKE$MATRIX
] ;
(* ################################################################################################################################################################ *)
(* GENERATOR (SIGMA) *)
(* ################################################################################################################################################################ *)
ClearAll[O$GENERATOR$SIGMA] ;
O$GENERATOR$SIGMA::usage = "
O$GENERATOR$SIGMA[] -- (feature generator) noise estimation (acts on a single sample, recommended transformation 'Standardize')
" ;
O$GENERATOR$SIGMA[                  (* -- (FEATURE GENERATOR) NOISE ESTIMATION *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Composition[
  Last,
  S$OPTIMAL$SVD$TRUNCATION,
  S$MAKE$MATRIX
] ;
(* ################################################################################################################################################################ *)
(* GENERATOR (HANKEL) *)
(* ################################################################################################################################################################ *)
ClearAll[O$GENERATOR$HANKEL] ;
O$GENERATOR$HANKEL::usage = "
O$GENERATOR$HANKEL[NUMBER] -- (feature generator) hankel filter noise estimation for given number of singular values <NUMBER> (integer) (acts on a single sample, recommended transformation 'Standardize')
" ;
O$GENERATOR$HANKEL[                 (* -- (FEATURE GENERATOR) HANKEL FILTER NOISE ESTIMATION *)
  NUMBER_,                          (* -- NUMBER OF SINGULAR VALUES TO USE (INTEGER) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Composition[
  Norm,
  Function[Subtract[Slot[1], S$SVD$FILTER[NUMBER, Slot[1]]]]
] ;
(* ################################################################################################################################################################ *)
(* GENERATOR (LEVEL) *)
(* ################################################################################################################################################################ *)
ClearAll[O$GENERATOR$LEVEL] ;
O$GENERATOR$LEVEL::usage = "
O$GENERATOR$LEVEL[NUMBER, DROP] -- (feature generator) hankel matrix singular values trend for total number of singular values <NUMBER> (integer) and number of singular values to drop <DROP> (integer) (acts on a single sample, recommended transformation 'Standardize')
" ;
O$GENERATOR$LEVEL[                  (* -- (FEATURE GENERATOR) HANKEL MATRIX SINGULAR VALUES TREND *)
  NUMBER_,                          (* -- NUMBER OF SINGULAR VALUES TO USE (INTEGER) *)
  DROP_,                            (* -- NUMBER OF SINGULAR VALUES TO DROP (INTEGER) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Composition[
  Median,
  Log10,
  Function[Drop[Slot[1], DROP]],
  Function[Divide[Slot[1], Max[Slot[1]]]],
  Function[S$SVD$LIST[NUMBER, Slot[1]]]
] ;
(* ################################################################################################################################################################ *)
(* GENERATOR (RPCA) *)
(* ################################################################################################################################################################ *)
ClearAll[O$GENERATOR$RPCA] ;
O$GENERATOR$RPCA::usage = "
O$GENERATOR$RPCA[FACTOR, LIMIT] -- (feature generator) rpca error matrix columns norm for given tolerance factor <FACTOR> (real) and number of iterations <LIMIT> (integer) (acts on full sample matrix, recommended transformation 'Standardize')
" ;
O$GENERATOR$RPCA[                   (* -- (FEATURE GENERATOR) RPCA ERROR MATRIX COLUMNS NORM *)
  FACTOR_,                          (* -- TOLERANCE FACTOR (REAL) *)
  LIMIT_,                           (* -- NUMBER OF ITERATIONS (INTEGER) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Composition[
  Map[Norm],
  First,
  Rest,
  Function[S$RPCA[Slot[1], FACTOR, LIMIT, FilterRules[List[OPTIONS], Options[S$RPCA]]]]
] ;
(* ################################################################################################################################################################ *)
(* GENERATOR (DMD) *)
(* ################################################################################################################################################################ *)
ClearAll[O$GENERATOR$DMD] ;
O$GENERATOR$DMD::usage = "
O$GENERATOR$DMD[SIZE, COUNT, FLAG] -- (feature generator) dmd prdiction error for sample size <SIZE> (integer), number of steps <COUNT> (integer) and dmd matrix flag <FLAG> (logical) (acts on a single sample, recommended transformation 'Standardize')
" ;
O$GENERATOR$DMD[                    (* -- (FEATURE GENERATOR) DMD PREDICTION ERROR *)
  SIZE_,                            (* -- SAMPLE SIZE TO USE (INTEGER) *)
  COUNT_,                           (* -- NUMBER OF PREDICTION STEPS (INTEGER) *)
  FLAG_,                            (* -- FLAG TO RECOMPUTE DMD MATRIX (LOGICAL) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Composition[
  Norm,
  Function[Subtract[Take[Slot[1], List[Plus[SIZE, 1], Plus[SIZE, COUNT]]], Take[S$DMD$FORWARD[Slot[1], SIZE, COUNT, FLAG], List[Plus[SIZE, 1], Plus[SIZE, COUNT]]]]]
] ;
(* ################################################################################################################################################################ *)
(* ANOMALY (CENTER) *)
(* ################################################################################################################################################################ *)
ClearAll[O$CENTER] ;
O$CENTER::usage = "
O$CENTER -- (robust) feature center estimation function used for anomaly detection, e.g. Mean, Median (default) or Fit[#, List[1], Null, Rule[NormFunction, Norm[#, 1] &]] &
" ;
O$CENTER = Median ;
(* ################################################################################################################################################################ *)
(* ANOMALY (SPREAD) *)
(* ################################################################################################################################################################ *)
ClearAll[O$SPREAD] ;
O$SPREAD::usage = "
O$SPREAD -- (robust) fearure spread estimation function used for anomaly detction, e.g. Composition[Sqrt, Variance] or Composition[Sqrt, BiweightMidvariance] (default)
" ;
O$SPREAD = Composition[Sqrt, BiweightMidvariance] ;
(* ################################################################################################################################################################ *)
(* ANOMALY (DETECT) *)
(* ################################################################################################################################################################ *)
ClearAll[O$DETECT] ;
O$DETECT::usage = "
O$DETECT[PLANE, TBT, GENERATOR] -- feature generation and anomaly detection for given TbT plane <PLANE> (string), TbT data <TBT> (association) and feature generation function <GENERATOR> (function)
" ;
Options[O$DETECT] = List[
  Rule["SHIFT", True],              (* -- FLAG TO GENERATE SHIFTED SAMPLES (LOGICAL), GENERATOR FUNCTION ACTS ON FLATTENED LIST OF SAMPLES FROM ALL SIGNALS *)
  Rule["STEP", 2^4],                (* -- SHIFT STEP SIZE (INTEGER), SAMPLES OF LENGTH 'TEST' (SEE TBT ASSOCIATION) ARE GENERATED FROM EACH SIGNAL USING GIVEN SHIFT STEP, I.E. STARTING POSITION IS VARIED *)
  Rule["PROCESS", Standardize],     (* -- DATA TRANSFORMATION FUNCTION TO APPLY TO EACH SAMPLE (FUNCTION), E.G. 'Rescale' OR 'Standardize' *)
  Rule["CENTER", O$CENTER],         (* -- (ROBUST) FEATURE CENTER ESTIMATOR (FUNCTION) *)
  Rule["SPREAD", O$SPREAD],         (* -- (ROBUST) FEATURE SPREAD ESTIMATOR (FUNCTION) *)
  Rule["THRESHOLD", 4.0],           (* -- THRESHOLD VALUE IN UNITS OF SPEAD (REAL), USE 3.0-6.0 WITHOUT SCAN AND LARGE VALUE WITH SCAN *)
  Rule["SCAN", False],              (* -- FLAG TO SCAN THRESHOLD VALUE (LOGICAL), VALUE IS SCANED STARTING FROM ESTIMATED SPREAD UPTO SPREAD TIMES GIVEN THRESHOLD VALUE *)
  Rule["SIZE", 250],                (* -- NUMBER OF SCAN STEPS (INTEGER) *)
  Rule["FRACTION", 0.75],           (* -- FRACTION CUTOFF FOR SCAN MODE (REAL) USED TO SELECT ANOMALY CANDIDATES, IGNORED IF PASS IS ONE *)
  Rule["PASS", 1],                  (* -- NUMBER OF PASSES (INTEGER), EACH PASS REMOVES CURRENT ANOMALY CANDIDATES FROM INPUT TBT DATA USING GIVEN THRESHOLD VALUE, DO NOT USE IN SCAN MODE *)
  Rule["MAP", ParallelMap]          (* -- MAP FUNCTION FOR FEATURE COMPUTATION *)
] ;
O$DETECT[                           (* -- FEATURE GENERATION AND ANOMALY DETECTION *)
  PLANE_String,                     (* -- PLANE TO USE (STRING) *)
  TBT_Association,                  (* -- TBT DATA (ASSOCIATION) *)
  GENERATOR_,                       (* -- FEATURE GENERATOR FUNCTION (FUNCTION) *)
  OPTIONS:OptionsPattern[]          (* -- OPTIONS(S) *)
] := Block[
  {SHIFT, STEP, PROCESS, CENTER, SPREAD, THRESHOLD, SCAN, SIZE, FRACTION, PASS, MAP, PICK, MATRIX, LENGTH, COUNT, TOTAL, LOOP, LIST, DATA, LINE, GATE, CASE, EXIT},
  {SHIFT, STEP, PROCESS, CENTER, SPREAD, THRESHOLD, SCAN, SIZE, FRACTION, PASS, MAP} = OptionValue[O$DETECT, Map[First, Options[O$DETECT]]] ;
  PICK = Range[Length[TBT["DATA", PLANE]]] ;
  MATRIX = TBT["DATA", PLANE] ;
  LENGTH = If[SHIFT, TBT["TEST"], TBT["SIZE"]] ;
  MATRIX = Map[Function[Partition[Slot[1], LENGTH, STEP]], MATRIX] ;
  MATRIX = Map[Map[PROCESS], MATRIX] ;
  COUNT = Length[First[MATRIX]] ;
  MATRIX = Flatten[MATRIX, 1] ;
  MATRIX = MAP[GENERATOR, MATRIX] ;
  MATRIX = Partition[MATRIX, COUNT] ;
  TOTAL = 0 ;
  Do[
    TOTAL++ ;
    LIST = Flatten[Transpose[ConstantArray[PICK, COUNT]]] ;
    DATA = Flatten[MATRIX[[PICK]], 1] ;
    LINE = CENTER[DATA] ;
    GATE = SPREAD[DATA] ;
    DATA = Abs[Subtract[DATA, LINE]] ;
    GATE = Rescale[GATE, MinMax[DATA]] ;
    DATA = Rescale[DATA, MinMax[DATA]] ;
    GATE = If[SCAN, Subdivide[GATE, Min[List[1, Times[GATE, THRESHOLD]]], Subtract[SIZE, 1]], Times[GATE, THRESHOLD]] ;
    CASE = If[
      SCAN,
      Select[Merge[Map[Function[Counts[Pick[LIST, Thread[GreaterEqual[DATA, Slot[1]]], True]]], GATE], Length], Curry[GreaterEqual][Times[FRACTION, Length[GATE]]]],
      Counts[Pick[LIST, Thread[GreaterEqual[DATA, GATE]], True]]
    ] ;
    EXIT[LOOP] = Association[
      Rule["DATA", Transpose[Partition[DATA, COUNT]]],
      Rule["CASE", CASE]
    ] ;
    PICK = DeleteCases[PICK, Apply[Alternatives, Keys[CASE]]],
    List[LOOP, 1, PASS]
  ] ;
  Merge[Map[EXIT, Range[TOTAL]], Identity]
] ;
(* ################################################################################################################################################################ *)
(* ANOMALY (NORMAL) *)
(* ################################################################################################################################################################ *)
ClearAll[O$DETECT$IS$NORMAL] ;
O$DETECT$IS$NORMAL::usage = "
O$DETECT$IS$NORMAL[DATA, LIMIT] -- test data normality for given list of candidate associations from different methods <DATA> (list) and maximum allowed empty associations <LIMIT> (integer)
" ;
O$DETECT$IS$NORMAL[                 (* -- NORMAL DATA TEST (COUNT NUMBER OF EMPTY CANDIDATES) *)
  DATA_List,                        (* -- LIST OF CANDIDATES (LIST) *)
  LIMIT_Integer                     (* -- MAXIMIM NUMBER OF ALLOWED EMPLY LISTS (INTEGER) *)
] := LessEqual[LIMIT, Count[DATA, Association[]]] ;
(* ################################################################################################################################################################ *)
(* ANOMALY (VALIDATE) *)
(* ################################################################################################################################################################ *)
ClearAll[O$DETECT$VALIDATE$CUT] ;
O$DETECT$VALIDATE$CUT::usage = "
O$DETECT$VALIDATE$CUT[DATA, LIMIT] -- validate detectors for given anomaly candidates <DATA> (list of associations) and acceptance factor <LIMIT> (real)
" ;
O$DETECT$VALIDATE$CUT[              (* -- VALIDATE DETECTORS *)
  DATA_List,                        (* -- LIST OF CANDIDATES (LIST OF ASSOCIATIONS) *)
  LIMIT_Real:0.1                    (* -- ACCEPTANCE LIMIT (REAL), FRACTION WITH RESPECT TO MAXIMUM TOTAL SCORE *)
] := Keys[Select[Merge[DATA, Total], Curry[GreaterEqual][Times[LIMIT, Max[Values[Merge[DATA, Total]]]]]]] ;
(* ################################################################################################################################################################ *)
(* ANOMALY (VALIDATE) *)
(* ################################################################################################################################################################ *)
ClearAll[O$DETECT$VALIDATE$FIT] ;
O$DETECT$VALIDATE$FIT::usage = "
O$DETECT$VALIDATE$FIT[DATA, LIMIT] -- validate detectors for given anomaly candidates <DATA> (list of associations) and acceptance factor <LIMIT> (real)
" ;
O$DETECT$VALIDATE$FIT[              (* -- VALIDATE DETECTORS *)
  DATA_List,                        (* -- LIST OF CANDIDATES (LIST OF ASSOCIATIONS) *)
  LIMIT_Real:4.0                    (* -- ACCEPTANCE LIMIT (REAL), FACTOR WITH RESPECT TO BASELINE *)
] := Block[
  {KEY, VALUE, LINE, SELECT},
  {KEY, VALUE} = Transpose[Map[Apply[List], Normal[Merge[DATA, Total]]]] ;
  LINE = Fit[N[VALUE], List[1], Null, Rule[NormFunction, Function[Norm[Slot[1], 1]]]] ;
  Keys[Select[Association[Thread[Rule[KEY, VALUE]]], Curry[GreaterEqual][Times[LIMIT, LINE]]]]
] ;
(* ################################################################################################################################################################ *)
(* ANOMALY (CONTAMINATION) *)
(* ################################################################################################################################################################ *)
ClearAll[O$DETECT$CONTAMINATION] ;
O$DETECT$CONTAMINATION::usage = "
O$DETECT$CONTAMINATION[TOTAL, DATA, LIMIT] -- estimate data contamination fraction for given total number of candidates <TOTAL> (integer), list of anomaly candidates <DATA> (list of associations), validator <VALIDATOR> (function) and corresponding threshold value <LIMIT> (real)
" ;
O$DETECT$CONTAMINATION[             (* -- ESTIMATE CONTAMINATION FRACTION *)
  TOTAL_Integer,                    (* -- TOTAL NUMBER OF CANDIDATES (INTEGER) *)
  DATA_List,                        (* -- LIST OF CANDIDATES (LIST) *)
  VALIDATOR_:O$DETECT$VALIDATE$FIT, (* -- VALIDATOR (FUNCTION) *)
  LIMIT_Real:4.0                    (* -- ACCEPTANCE LIMIT (REAL) *)
] := N[Divide[Length[VALIDATOR[DATA, LIMIT]], TOTAL]] ;
(* ################################################################################################################################################################ *)
(* ANOMALY (KDE) *)
(* ################################################################################################################################################################ *)
ClearAll[O$KDE] ;
O$KDE::usage = "
O$KDE[DATA, CONTAMINATION] -- KDE based anomaly detection for given data matrix <DATA> (matrix) and contamination factor <CONTAMINATION> (real)
" ;
O$KDE[                              (* -- KDE ANOMALY DETECTION *)
  DATA_?MatrixQ,                    (* -- DATA MATRIX *)
  CONTAMINATION_Real,               (* -- CONTAMINATION/ACCEPTANCE THRESHOLD (REAL) *)
  OPTIONS:OptionsPattern[]          (* -- OPTION(S) *)
] := Block[
  {DETECTOR, INDEX},
  DETECTOR = AnomalyDetection[
    RandomSample[DATA, Length[DATA]],
    FilterRules[List[OPTIONS], Options[AnomalyDetection]],
    Rule[PerformanceGoal, "Quality"],
    Rule[Method, "KernelDensityEstimation"],
    Rule[AcceptanceThreshold, CONTAMINATION],
    Rule[TrainingProgressReporting, None]
  ] ;
  INDEX = Map[DETECTOR, DATA] ;
  Flatten[Position[INDEX, True]]
] ;
(* ################################################################################################################################################################ *)
(* ANOMALY (LOF) *)
(* ################################################################################################################################################################ *)
ClearAll[O$LOF] ;
O$LOF::usage = "
O$LOF[DATA, COUNT, FACTOR] -- (SKLEARN) LOF anomaly detection for given data matrix <DATA> (matrix), number of neighbors <COUNT> (integer) and contamination factor <CONTAMINATION> (real)
" ;
O$LOF[                              (* -- (SKLEARN) LOF ANOMALY DETECTION *)
  DATA_?MatrixQ,                    (* -- DATA MATRIX *)
  COUNT_Integer,                    (* -- NUMBER OF NEIGHBORS TO USE (INTEGER) *)
  CONTAMINATION_Real                (* -- CONTAMINATION (REAL) *)
] := Block[
  {SESSION, FUNCTION, INDEX},
  SESSION = StartExternalSession["Python"] ;
  FUNCTION = ExternalFunction[
    SESSION,
    "def function(neighbors, contamination, matrix):
	from sklearn.neighbors import LocalOutlierFactor
	indicator = LocalOutlierFactor(n_neighbors = neighbors, contamination = contamination)
	data = indicator.fit_predict(matrix)
	return data
	"
  ] ;
  INDEX = Normal[FUNCTION[COUNT, CONTAMINATION, DATA]] ;
  DeleteObject[SESSION] ;
  Flatten[Position[INDEX, Subtract[0, 1]]]
] ;
(* ################################################################################################################################################################ *)
(* ANOMALY (IFOREST) *)
(* ################################################################################################################################################################ *)
ClearAll[O$IFOREST] ;
O$IFOREST::usage = "
O$IFOREST[DATA, COUNT, CONTAMINATION] -- (SKLEARN) Isolation Forest anomaly detection for given data matrix <DATA> (matrix), number of estimators <COUNT> (integer) and contamination factor <CONTAMINATION> (real)
" ;
O$IFOREST[                          (* -- (SKLEARN) ISOLATION FOREST ANOMALY DETECTION *)
  DATA_?MatrixQ,                    (* -- DATA MATRIX *)
  COUNT_Integer,                    (* -- NUMBER OF NEIGHBORS TO USE (INTEGER) *)
  CONTAMINATION_Real                (* -- CONTAMINATION (REAL) *)
] := Block[
  {SESSION, FUNCTION, INDEX},
  SESSION = StartExternalSession["Python"] ;
  FUNCTION = ExternalFunction[
    SESSION,
    "def function(estimators, contamination, matrix):
	from sklearn.ensemble import IsolationForest
	indicator = IsolationForest(n_estimators = estimators, contamination = contamination)
	data = indicator.fit(matrix).predict(matrix)
	return data
	"
  ] ;
  INDEX = Normal[FUNCTION[COUNT, CONTAMINATION, DATA]] ;
  DeleteObject[SESSION] ;
  Flatten[Position[INDEX, Subtract[0, 1]]]
] ;

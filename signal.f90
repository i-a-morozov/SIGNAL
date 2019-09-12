! SIGNAL,2019,I.A.MOROZOV@INP.NSK.SU
! SIMPLIFIED FORTRAN IMPLEMENTATION OF S$FREQUENCY[] FUNCTION (WITH PARABOLIC INTERPOLATION OF REFINED SPECTRA OR MAXLOC BRENT)

! MAIN FUNCTION(S):
! <OUT>=FREQUENCY(<AVE>,<WIN>,<SIG>)         -- WITH PARABOLIC INTERPOLATION
! <OUT>=FREQUENCY_BRENT(<AVE>,<WIN>,<SIG>)   -- BRENT (MIGHT BE LESS ACCURATE THEN PARABOLIC, NOT COMPATIBLE WITH OPENACC ) 
! <AVE> -- WINDOW SUM (REAL)
! <WIN> -- WINDOW (REAL) WITH LENGTH = NUM (MUST BE POWER OF TWO)
! <SIG> -- INPUT SIGNAL (REAL) WITH LENGTH = 2*NUM (MUST BE POWER OF TWO), [R_1,I_1,R_2,I_2,...,R_NUM,I_NUM]
! <OUT> -- FREQUENCY ESTIMATION (REAL)

! SEE EXAMPLES AT THE END OF FILE 

MODULE SIGNAL

    USE ISO_FORTRAN_ENV    

    IMPLICIT NONE
    
    PRIVATE
    
    INTEGER,PARAMETER  :: RK=REAL64 ! REAL KIND
    INTEGER,PARAMETER  :: NUM=1024                      ! SIGNAL LENGTH (POWER OF TWO)
    INTEGER,PARAMETER  :: FLA=1                         ! COMPLEX SIGNAL FLAG (0/1 FOR REAL/COMPLEX INPUT)
    REAL(RK),PARAMETER :: PI=3.141592653589793238460_RK ! PI
    REAL(RK),PARAMETER :: TWO_PI=2.0_RK*PI              ! 2*PI
    
    INTEGER,PARAMETER  :: BRENT_MAX=250                 ! MAX # OF ITERATIONS S FOR BRENT
    REAL(RK),PARAMETER :: BRENT_TOL=1.E-14_RK           ! BRENT TOLERANCE
    INTEGER,PARAMETER  :: INTERPOLATION_CASE=1          ! INTERPOLATION TYPE
    INTEGER,PARAMETER  :: INTERPOLATION_POINTS=10       ! # OF INTERPOLATING POINTS
    
    
    PUBLIC :: NUM
    PUBLIC :: FLA
    PUBLIC :: RK
    PUBLIC :: PI
    PUBLIC :: TWO_PI
    PUBLIC :: WINDOW_
    PUBLIC :: FREQUENCY_
    
    PUBLIC :: BRENT_MAX
    PUBLIC :: BRENT_TOL
    PUBLIC :: FREQUENCY_BRENT_
    REAL(RK),DIMENSION(2*INTERPOLATION_POINTS+1) :: X_LIST
    REAL(RK),DIMENSION(2*INTERPOLATION_POINTS+1) :: Y_LIST
    ABSTRACT INTERFACE
        PURE FUNCTION GENERIC_FUNCTION(ARG)
            USE ISO_FORTRAN_ENV
            IMPLICIT NONE
            INTEGER,PARAMETER  :: RK=REAL64
            REAL(RK), INTENT(IN) :: ARG
            REAL(RK) :: GENERIC_FUNCTION
        END FUNCTION GENERIC_FUNCTION
    END INTERFACE    
            
    CONTAINS
    
    ! DISCRETE FOURIER TRANSFORM (NRF77)
    PURE SUBROUTINE FFT_(NUM,DIR,ARR)
        !$ACC ROUTINE SEQ
        INTEGER,INTENT(IN) :: NUM
        INTEGER,INTENT(IN) :: DIR
        REAL(RK),DIMENSION(2*NUM),INTENT(INOUT) :: ARR
        INTEGER :: N,I,J,M,LIM,STE
        REAL(RK) :: PIM,PRE,ANG,WR,WI,WPR,WPI,WX,MUL
        N=2*NUM
        J=1
        MUL = REAL(DIR,RK)*TWO_PI
        DO I=1,N,2
            IF(J.GT.I)THEN
                PRE=ARR(J)
                PIM=ARR(J+1)
                ARR(J)=ARR(I)
                ARR(J+1)=ARR(I+1)
                ARR(I)=PRE
                ARR(I+1)=PIM
            ENDIF
            M=N/2
1           IF ((M.GE.2).AND.(J.GT.M)) THEN
                J=J-M
                M=M/2
                GOTO 1
            ENDIF
            J=J+M
        ENDDO
        LIM=2
2       IF (N.GT.LIM) THEN
            STE=2*LIM
            ANG=MUL/REAL(LIM,RK)
            WPI=SIN(ANG)
            ANG=SIN(0.5_RK*ANG)
            WPR=-2.0_RK*ANG*ANG
            WR=1.0_RK
            WI=0.0_RK
            DO M=1,LIM,2
                DO I=M,N,STE
                    J=I+LIM
                    PRE=WR*ARR(J)-WI*ARR(J+1)
                    PIM=WR*ARR(J+1)+WI*ARR(J)
                    ARR(J)=ARR(I)-PRE
                    ARR(J+1)=ARR(I+1)-PIM
                    ARR(I)=ARR(I)+PRE
                    ARR(I+1)=ARR(I+1)+PIM
                ENDDO
                WX=WR
                WR=WR*WPR-WI*WPI+WR
                WI=WI*WPR+WX*WPI+WI
            ENDDO
            LIM=STE
            GOTO 2
        ENDIF
    END SUBROUTINE FFT_
    
    ! DISCRETE (LINEAR) FRACTIONAL FOURIER TRANSFORM
    PURE SUBROUTINE FFRFT_(NUM,PAR,ARR)
        !$ACC ROUTINE SEQ
        INTEGER,INTENT(IN) :: NUM
        REAL(RK),INTENT(IN) :: PAR
        REAL(RK),DIMENSION(2*NUM),INTENT(INOUT) :: ARR
        INTEGER :: I
        REAL(RK) :: FAC
        REAL(RK),DIMENSION(NUM)   :: MUL,COS_MUL,SIN_MUL
        REAL(RK),DIMENSION(4*NUM) :: ONE,TWO,TRE
        REAL(RK),DIMENSION(2*NUM) :: TMP
        FAC=PAR*PI/REAL(NUM,RK)
        MUL=FAC*REAL([(I-1,I=1,NUM,1)],RK)**2
        COS_MUL=COS(MUL)
        SIN_MUL=SIN(MUL)
        ONE(1:2*NUM:2)=ARR(1:2*NUM:2)*COS_MUL-ARR(2:2*NUM:2)*SIN_MUL
        ONE(2:2*NUM:2)=ARR(1:2*NUM:2)*SIN_MUL+ARR(2:2*NUM:2)*COS_MUL
        TWO(1:2*NUM:2) = +COS_MUL
        TWO(2:2*NUM:2) = -SIN_MUL
        ONE(2*NUM+1:4*NUM:1)=0.0_RK
        MUL=-FAC*REAL([(I-1-2*NUM,I=NUM+1,2*NUM,1)],RK)**2
        TWO(2*NUM+1:4*NUM:2)=COS(MUL)
        TWO(2*NUM+2:4*NUM:2)=SIN(MUL)
        CALL FFT_(2*NUM,+1,ONE)
        CALL FFT_(2*NUM,+1,TWO)
        TRE=ONE
        ONE(1:4*NUM:2)=TRE(1:4*NUM:2)*TWO(1:4*NUM:2)-TRE(2:4*NUM:2)*TWO(2:4*NUM:2)
        ONE(2:4*NUM:2)=TRE(1:4*NUM:2)*TWO(2:4*NUM:2)+TRE(2:4*NUM:2)*TWO(1:4*NUM:2)
        CALL FFT_(2*NUM,-1,ONE)
        ARR=ONE(1:2*NUM:1)/REAL(2*NUM,RK)
        TMP=ARR
        ARR(1:2*NUM:2)=TMP(1:2*NUM:2)*COS_MUL-TMP(2:2*NUM:2)*SIN_MUL
        ARR(2:2*NUM:2)=TMP(1:2*NUM:2)*SIN_MUL+TMP(2:2*NUM:2)*COS_MUL
    END SUBROUTINE FFRFT_

    ! FACTORIAL
    INTEGER PURE FUNCTION FACTORIAL_(N) RESULT(M)
        !$ACC ROUTINE SEQ
        INTEGER,INTENT(IN) :: N
        INTEGER :: I
        M = 1
        DO I=1,N,1
            M=M*I
        END DO
    END FUNCTION FACTORIAL_
    
    ! WINDOW
    PURE SUBROUTINE WINDOW_(ORD,ARR) 
        !$ACC ROUTINE SEQ
        INTEGER,INTENT(IN):: ORD
        REAL(RK),INTENT(OUT),DIMENSION(NUM) :: ARR
        INTEGER :: I
        REAL(RK) :: MUL
        MUL=2._RK**ORD*REAL(FACTORIAL_(ORD),RK)**2/REAL(FACTORIAL_(2*ORD),RK)
        ARR=MUL*(1._RK+COS(TWO_PI*(REAL([(I-1,I=1,NUM,1)],RK)/REAL(NUM,RK)-0.5_RK)))**ORD
    END SUBROUTINE WINDOW_
    
    ! MAXLOC (PGFORTRAN 19.4)
    PURE INTEGER FUNCTION MAXLOC_(ARR)
        !$ACC ROUTINE SEQ
        REAL(RK),DIMENSION(:),INTENT(IN) :: ARR
        REAL(RK),DIMENSION(SIZE(ARR)) :: LIS
        REAL(RK) :: VAL
        INTEGER :: I
        VAL=MAXVAL(ARR)
        DO I=1,SIZE(ARR),1
            IF(ARR(I)==VAL)THEN
                MAXLOC_=I
                RETURN
            END IF
        END DO
    END FUNCTION MAXLOC_
    
    ! MINLOC (PGFORTRAN 19.4)
    PURE INTEGER FUNCTION MINLOC_(ARR)
        !$ACC ROUTINE SEQ
        REAL(RK),DIMENSION(:),INTENT(IN) :: ARR
        REAL(RK),DIMENSION(SIZE(ARR)) :: LIS
        REAL(RK) :: VAL
        INTEGER :: I
        VAL=MINVAL(ARR)
        DO I=1,SIZE(ARR),1
            IF(ARR(I)==VAL)THEN
                MINLOC_=I
                RETURN
            END IF
        END DO
    END FUNCTION MINLOC_
    
    ! FREQUENCY SEARCH (LARGEST AMPLITUDE)
    ! N.B. INTRISIC FUNCTION MAXLOC IS NOT SUPPORTED IN PGFORTRAN 19.4 FOR ACC (NOT ACC ROUTINE INFO WHEN USED IN SUBS)
    ! MAXLOC_ FUNCTION (SEE ABOVE) IS USED (PERFORMANCE IS NEAR IDENTICAL)
    ! COMMENT/UNCOMMENT FOR GFORTRAN
    PURE REAL(RK) FUNCTION FREQUENCY_(WINDOW_TOTAL,WINDOW,SIGNAL)
        !$ACC ROUTINE SEQ
        REAL(RK),INTENT(IN) :: WINDOW_TOTAL
        REAL(RK),INTENT(IN),DIMENSION(NUM) :: WINDOW
        REAL(RK),INTENT(IN),DIMENSION(2*NUM) :: SIGNAL
        REAL(RK),DIMENSION(2*NUM) :: ARR,FOU        
        INTEGER :: FST,CND
        REAL(RK) :: FAC
        REAL(RK),DIMENSION(NUM) :: MUL,COS_MUL,SIN_MUL
        INTEGER :: I
        ! REMOVE (WEIGHTED) MEAN AND APPLY WINDOW
        ARR(1:2*NUM:2)=(SIGNAL(1:2*NUM:2)-SUM(WINDOW*SIGNAL(1:2*NUM:2))/WINDOW_TOTAL)*WINDOW
        ARR(2:2*NUM:2)=(SIGNAL(2:2*NUM:2)-SUM(WINDOW*SIGNAL(2:2*NUM:2))/WINDOW_TOTAL)*WINDOW
        ! FFT APPROXIMATION
        FOU=ARR
        CALL FFT_(NUM,1,FOU)
        MUL = FOU(1:2*NUM:2)**2+FOU(2:2*NUM:2)**2
        !PGFORTRAN/ACC
        FST=MAXLOC_(MUL(1:NUM/(2-FLA):1))
        !GFORTRAN
        !FST=MAXLOC(MUL(1:NUM/(2-FLA):1),1)
        ! MODULATE SIGNAL
        FAC=TWO_PI*REAL(FST-2,RK)/REAL(NUM,RK)
        MUL=FAC*REAL([(I-1,I=1,NUM,1)],RK)
        COS_MUL=COS(MUL)
        SIN_MUL=SIN(MUL)
        FOU=ARR
        ARR(1:2*NUM:2)=FOU(1:2*NUM:2)*COS_MUL-FOU(2:2*NUM:2)*SIN_MUL
        ARR(2:2*NUM:2)=FOU(1:2*NUM:2)*SIN_MUL+FOU(2:2*NUM:2)*COS_MUL
        ! FFRFT APPROXIMATION
        CALL FFRFT_(NUM,2.0_RK/REAL(NUM,RK),ARR)
        MUL=LOG10(SQRT(ARR(1:2*NUM:2)**2+ARR(2:2*NUM:2)**2)+1.E-16_RK)
        !PGFORTRAN/ACC
        CND=MAXLOC_(MUL)
        !GFORTRAN
        !CND=MAXLOC(MUL,1)
        ! PARABOLA APPROXIMATION
        FREQUENCY_=REAL(CND,RK)-0.5_RK+(MUL(-1+CND)-MUL(CND))/(MUL(-1+CND)-2._RK*MUL(CND)+MUL(1+CND))
        ! RESULT
        FREQUENCY_=(REAL(FST,RK)-2._RK+2._RK*(FREQUENCY_-1._RK)/REAL(NUM,RK))/REAL(NUM,RK)   
    END FUNCTION FREQUENCY_
    
    ! POLINT (NRF90)
    PURE SUBROUTINE POLINT_(XA,YA,X,Y,DY)
        REAL(RK),DIMENSION(:),INTENT(IN) :: XA
        REAL(RK),DIMENSION(:),INTENT(IN) :: YA
        REAL(RK),INTENT(IN)  :: X
        REAL(RK),INTENT(OUT) :: Y,DY
        INTEGER :: N,M,NS,I
        REAL(RK), DIMENSION(SIZE(XA)) :: C,D,DEN,HO
        N=SIZE(XA)
        C=YA
        D=YA
        HO=XA-X
        NS=MINLOC_(ABS(HO))
        Y=YA(NS)
        DO M=1,N-1
            I=N-M
            DEN(1:I)=HO(1:I)-HO(1+M:N)
            DEN(1:I)=(C(2:I+1)-D(1:I))/DEN(1:I)
            D(1:I)=HO(1+M:N)*DEN(1:I)
            C(1:I)=HO(1:I)*DEN(1:I)
            IF(2*NS<I)THEN
                DY=C(NS+1)
            ELSE
                DY=D(NS)
                NS=NS-1
            END IF
            Y = Y+DY
        END DO
    END SUBROUTINE POLINT_
    
    ! RATINT (NRF90)
    PURE SUBROUTINE RATINT_(XA,YA,X,Y,DY)
        REAL(RK),DIMENSION(:),INTENT(IN) :: XA
        REAL(RK),DIMENSION(:),INTENT(IN) :: YA
        REAL(RK),INTENT(IN)  :: X
        REAL(RK),INTENT(OUT) :: Y,DY
        INTEGER :: N,M,NS,I
        REAL(RK), DIMENSION(SIZE(XA)) :: C,D,DD,H,T
        REAL(RK), PARAMETER :: EPS = 1.0E-25_RK
        N=SIZE(XA)
        H=XA-X
        NS=MINLOC_(ABS(H))
        Y=YA(NS)
        IF(X==XA(NS)) THEN
            DY=0.0_RK
            RETURN
        END IF
        C=YA
        D=YA+EPS
        NS=NS-1
        DO M=1,N-1
            I=N-M
            T(1:I)=(XA(1:I)-X)*D(1:I)/H(1+M:N)
            DD(1:I)=T(1:I)-C(2:I+1)
            DD(1:I)=(C(2:I+1)-D(1:I))/DD(1:I)
            D(1:I)=C(2:I+1)*DD(1:I)
            C(1:I)=T(1:I)*DD(1:I)
            IF(2*NS<I)THEN
                DY=C(NS+1)
            ELSE
                DY=D(NS)
                NS=NS-1
            END IF
            Y=Y+DY
        END DO
    END SUBROUTINE RATINT_
    
    ! INTERPOLATING FUNCTION
    PURE FUNCTION INTERPOLATION_(ARG)
        REAL(RK),INTENT(IN) :: ARG
        REAL(RK) :: INTERPOLATION_
        REAL(RK) :: X,Y,DY
        X = ARG
        SELECT CASE(INTERPOLATION_CASE)
            CASE(1)
                CALL POLINT_(X_LIST,Y_LIST,X,INTERPOLATION_,DY)
            CASE(2)
                CALL RATINT_(X_LIST,Y_LIST,X,INTERPOLATION_,DY)
            CASE DEFAULT
                CALL POLINT_(X_LIST,Y_LIST,X,INTERPOLATION_,DY)
        END SELECT
    END FUNCTION INTERPOLATION_
    
    ! BRENT (NRF90)
    PURE FUNCTION BRENT_(AX,BX,CX,FUN)
        REAL(RK),INTENT(IN) :: AX
        REAL(RK),INTENT(IN) :: BX
        REAL(RK),INTENT(IN) :: CX
        PROCEDURE(GENERIC_FUNCTION) :: FUN
        REAL(RK),DIMENSION(2) :: BRENT_
        INTEGER,PARAMETER :: ITMAX=BRENT_MAX
        REAL(RK), PARAMETER :: CGOLD=0.3819660112501052_RK
        REAL(RK), PARAMETER :: ZEPS=TINY(1.0_RK)
        INTEGER :: ITER
        REAL(RK) :: A,B,D,E,ETEMP,FU,FV,FW,FX,P,Q,R,TOL1,TOL2,U,V,W,X,XM,XMIN
        A=MIN(AX,CX)
        B=MAX(AX,CX)
        V=BX
        W=V
        X=V
        E=0.0_RK
        FX=FUN(X)
        FV=FX
        FW=FX
        DO ITER=1,ITMAX,1
            XM=0.5_RK*(A+B)
            TOL1=BRENT_TOL*ABS(X)+ZEPS
            TOL2=2.0_RK*TOL1
            IF(ABS(X-XM)<=(TOL2-0.5_RK*(B-A)))THEN
                XMIN=X
                BRENT_=[XMIN,FX]
                RETURN
            END IF
            IF(ABS(E)>TOL1)THEN
                R=(X-W)*(FX-FV)
                Q=(X-V)*(FX-FW)
                P=(X-V)*Q-(X-W)*R
                Q=2.0_RK*(Q-R)
                if(Q>0.0_RK) P=-P
                Q=ABS(Q)
                ETEMP=E
                E=D
                IF(ABS(P)>=ABS(0.5_RK*Q*ETEMP).OR.P<=Q*(A-X).OR.P>=Q*(B-X))then
                    E=MERGE(A-X,B-X,X>=XM)
                    D=CGOLD*E
                ELSE
                    D=P/Q
                    U=X+D
                    IF(U-A<TOL2.OR.B-U<TOL2)D=SIGN(TOL1,XM-X)
                END IF
            ELSE
                E=MERGE(A-X,B-X,X>=XM)
                D=CGOLD*E
            END IF
            U=MERGE(X+D,X+SIGN(TOL1,D),ABS(D)>=TOL1)
            FU=FUN(U)
            IF(FU<=FX)THEN
                IF(U>=X)THEN
                    A=X
                ELSE
                    B=X
                END IF
                CALL SHFT_(V,W,X,U)
                CALL SHFT_(FV,FW,FX,FU)
            ELSE
                IF(U<X)THEN
                    A=U
                ELSE
                    B=U
                END IF
                IF (FU<=FW.OR.W==X) THEN
                    V=W
                    FV=FW
                    W=U
                    FW=FU
                ELSE IF (FU<=FV.OR.V==X.OR.V==W) THEN
                    V=U
                    FV=FU
                END IF
            END IF
        END DO
        CONTAINS
        PURE SUBROUTINE SHFT_(A,B,C,D)
            REAL(RK),INTENT(OUT) :: A
            REAL(RK),INTENT(INOUT) :: B
            REAL(RK),INTENT(INOUT) :: C
            REAL(RK),INTENT(IN) :: D
            A=B
            B=C
            C=D
        END SUBROUTINE SHFT_
    END FUNCTION BRENT_
    
    ! FREQUENCY SEARCH (BRENT) (LARGEST AMPLITUDE) (NOT PURE)
    REAL(RK) FUNCTION FREQUENCY_BRENT_(WINDOW_TOTAL,WINDOW,SIGNAL)
        REAL(RK),INTENT(IN) :: WINDOW_TOTAL
        REAL(RK),INTENT(IN),DIMENSION(NUM) :: WINDOW
        REAL(RK),INTENT(IN),DIMENSION(2*NUM) :: SIGNAL
        REAL(RK),DIMENSION(2*NUM) :: ARR,FOU        
        INTEGER :: FST,CND
        REAL(RK) :: FAC
        REAL(RK),DIMENSION(NUM) :: MUL,COS_MUL,SIN_MUL
        INTEGER :: I
        REAL(RK),DIMENSION(2) :: OUTPUT
        REAL(RK) :: AX,BX,CX        
        ! REMOVE (WEIGHTED) MEAN AND APPLY WINDOW
        ARR(1:2*NUM:2)=(SIGNAL(1:2*NUM:2)-SUM(WINDOW*SIGNAL(1:2*NUM:2))/WINDOW_TOTAL)*WINDOW
        ARR(2:2*NUM:2)=(SIGNAL(2:2*NUM:2)-SUM(WINDOW*SIGNAL(2:2*NUM:2))/WINDOW_TOTAL)*WINDOW
        ! FFT APPROXIMATION
        FOU=ARR
        CALL FFT_(NUM,1,FOU)
        MUL = FOU(1:2*NUM:2)**2+FOU(2:2*NUM:2)**2
        FST=MAXLOC(MUL(1:NUM/(2-FLA):1),1)
        ! MODULATE SIGNAL
        FAC=TWO_PI*REAL(FST-2,RK)/REAL(NUM,RK)
        MUL=FAC*REAL([(I-1,I=1,NUM,1)],RK)
        COS_MUL=COS(MUL)
        SIN_MUL=SIN(MUL)
        FOU=ARR
        ARR(1:2*NUM:2)=FOU(1:2*NUM:2)*COS_MUL-FOU(2:2*NUM:2)*SIN_MUL
        ARR(2:2*NUM:2)=FOU(1:2*NUM:2)*SIN_MUL+FOU(2:2*NUM:2)*COS_MUL
        ! FFRFT APPROXIMATION
        CALL FFRFT_(NUM,2.0_RK/REAL(NUM,RK),ARR)
        MUL=LOG10(SQRT(ARR(1:2*NUM:2)**2+ARR(2:2*NUM:2)**2)+1.E-16_RK)
        CND=MAXLOC(MUL,1)
        ! BRENT
        X_LIST=[(REAL(I+CND,RK),I=-INTERPOLATION_POINTS,INTERPOLATION_POINTS,1)]
        Y_LIST=MUL(CND-INTERPOLATION_POINTS:CND+INTERPOLATION_POINTS)
        Y_LIST=1.0_RK/Y_LIST
        AX=X_LIST(INTERPOLATION_POINTS+0)
        BX=X_LIST(INTERPOLATION_POINTS+1)
        CX=X_LIST(INTERPOLATION_POINTS+2)
        OUTPUT=BRENT_(AX,BX,CX,INTERPOLATION_)
        FREQUENCY_BRENT_=OUTPUT(1)       
        ! RESULT
        FREQUENCY_BRENT_=(REAL(FST,RK)-2._RK+2._RK*(FREQUENCY_BRENT_-1._RK)/REAL(NUM,RK))/REAL(NUM,RK)    
    END FUNCTION FREQUENCY_BRENT_

END MODULE SIGNAL

! ! EXAMPLE (NUM=2**9,FLA=1)
! ! gfortran -o example -O3 -ffast-math -march=native signal.f90 example.f90
! ! example.f90:
! PROGRAM EXAMPLE
!   USE SIGNAL
!   IMPLICIT NONE
!   REAL(RK),DIMENSION(2*NUM) :: ARR
!   REAL(RK),DIMENSION(NUM)   :: WIN
!   REAL(RK)                  :: AVE
!   INTEGER                   :: I
!   REAL(RK),PARAMETER        :: FRE=0.623456789_RK
!   ARR(1:2*NUM:2) = COS(2._RK*PI*FRE*REAL([(I,I=1,NUM,1)],RK))
!   ARR(1:2*NUM:2) = ARR(1:2*NUM:2) + 0.1_RK
!   IF(FLA==0) THEN 
!     ARR(2:2*NUM:2) = 0.0_RK
!   ELSE 
!     ARR(2:2*NUM:2) = -SIN(2._RK*PI*FRE*REAL([(I,I=1,NUM,1)],RK))
!   END IF
!   CALL WINDOW_(1,WIN)
!   AVE = SUM(WIN)
!   WRITE(*,*) FRE
!   WRITE(*,*) FREQUENCY_(AVE,WIN,ARR)
!   WRITE(*,*) FREQUENCY_BRENT_(AVE,WIN,ARR)
! END PROGRAM EXAMPLE

! ! EXAMPLE (NUM=2**9,FLA=1) (OPENMP)
! ! gfortran -o example -O3 -ffast-math -march=native -fopenmp signal.f90 example.f90
! ! example.f90:
! PROGRAM EXAMPLE
!   USE SIGNAL
!   IMPLICIT NONE
!   REAL(RK),DIMENSION(2*NUM)         :: ARR
!   REAL(RK),DIMENSION(NUM)           :: WIN
!   REAL(RK)                          :: AVE
!   INTEGER                           :: I
!   REAL(RK),PARAMETER                :: FRE=0.623456789_RK
!   REAL(RK),DIMENSION(10000,2*NUM)   :: DAT
!   REAL(RK),DIMENSION(10000)         :: OUT
!   ARR(1:2*NUM:2) = COS(2._RK*PI*FRE*REAL([(I,I=1,NUM,1)],RK))
!   ARR(1:2*NUM:2) = ARR(1:2*NUM:2) + 0.1_RK
!   IF(FLA==0) THEN 
!     ARR(2:2*NUM:2) = 0.0_RK
!   ELSE 
!     ARR(2:2*NUM:2) = -SIN(2._RK*PI*FRE*REAL([(I,I=1,NUM,1)],RK))
!   END IF
!   CALL WINDOW_(1,WIN)
!   AVE = SUM(WIN)
!   DO I=1,10000,1
!     DAT(I,:) = ARR
!   END DO
!   !$OMP PARALLEL DO
!   DO I=1,10000,1
!     OUT(I) = FREQUENCY_(AVE,WIN,DAT(I,:))
!   END DO
!   !$OMP END PARALLEL DO
!   WRITE(*,*) OUT(1)-FRE
! END PROGRAM EXAMPLE

! ! EXAMPLE (NUM=2**9,FLA=1) (OPENACC) (COMMENT ALL BRENT RELATED PARTS PART TO COMPILE)
! ! pgfortran -o example -fast -ta=nvidia signal.f90 example.f90
! ! example.f90:
! PROGRAM EXAMPLE
!   USE SIGNAL
!   IMPLICIT NONE
!   REAL(RK),DIMENSION(2*NUM)         :: ARR
!   REAL(RK),DIMENSION(NUM)           :: WIN
!   REAL(RK)                          :: AVE
!   INTEGER                           :: I
!   REAL(RK),PARAMETER                :: FRE=0.623456789_RK
!   REAL(RK),DIMENSION(10000,2*NUM)   :: DAT
!   REAL(RK),DIMENSION(10000)         :: OUT
!   ARR(1:2*NUM:2) = COS(2._RK*PI*FRE*REAL([(I,I=1,NUM,1)],RK))
!   ARR(1:2*NUM:2) = ARR(1:2*NUM:2) + 0.1_RK
!   IF(FLA==0) THEN 
!     ARR(2:2*NUM:2) = 0.0_RK
!   ELSE 
!     ARR(2:2*NUM:2) = -SIN(2._RK*PI*FRE*REAL([(I,I=1,NUM,1)],RK))
!   END IF
!   CALL WINDOW_(1,WIN)
!   AVE = SUM(WIN)
!   DO I=1,10000,1
!     DAT(I,:) = ARR
!   END DO
!   !$ACC DATA COPYIN(DAT(:,:),WIN(:),AVE) COPYOUT(OUT(:))
!   !$ACC PARALLEL LOOP
!   DO I=1,10000,1
!     OUT(I) = FREQUENCY_(AVE,WIN,DAT(I,:))
!   END DO
!   !$ACC END PARALLEL LOOP
!   !$ACC END DATA
!   WRITE(*,*) OUT(1)-FRE
! END PROGRAM EXAMPLE

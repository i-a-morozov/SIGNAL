! SIGNAL,2019,I.A.MOROZOV@INP.NSK.SU
! SIMPLIFIED FORTRAN IMPLEMENTATION OF S$FREQUENCY[] FUNCTION (WITH PARABOLIC INTERPOLATION OF REFINED SPECTRA OR BRENT)

! MAIN FUNCTION(S):
! <OUT>=FREQUENCY(<AVE>,<WIN>,<SIG>)         -- WITH PARABOLIC INTERPOLATION
! <OUT>=FREQUENCY_BRENT(<AVE>,<WIN>,<SIG>)   -- INTERPOLATION + BRENT (MIGHT BE LESS ACCURATE THEN PARABOLIC) 
! <AVE> -- WINDOW SUM (COMPLEX)
! <WIN> -- WINDOW (COMPLEX) WITH POWER OF TWO LENGTH, NO LENGTH CHECK IS PERFORMED
! <SIG> -- INPUT SIGNAL (COMPLEX) WITH POWER OF TWO LENGTH, NO LENGTH CHECK IS PERFORMED
! <OUT> -- FREQUENCY ESTIMATION (REAL)

! ! EXAMPLE (NUM=512,FLA=1,BRENT_MAX=250,BRENT_TOL=1.E-14_RK,INTERPOLATION_CASE=1,INTERPOLATION_POINTS=10)
! PROGRAM EXAMPLE
!   USE SIGNAL
!   
!   IMPLICIT NONE
!   
!   COMPLEX(RK),DIMENSION(NUM) :: ARR
!   COMPLEX(RK),DIMENSION(NUM) :: WIN
!   COMPLEX(RK)                :: AVE
!   INTEGER                    :: I
!   REAL(RK),PARAMETER         :: FRE=0.123456789_RK
!   
!   SELECT CASE(FLA)
!     CASE(0)
!       DO I=1,NUM,1
!         ARR(I)=SIN(2._RK*PI*FRE*REAL(I,RK))
!       END DO
!     CASE(1)
!       DO I=1,NUM,1
!         ARR(I)=CMPLX(SIN(2._RK*PI*FRE*REAL(I,RK)),COS(2._RK*PI*0.123456789_RK*REAL(I,RK)),RK)
!       END DO
!     CASE DEFAULT
!       DO I=1,NUM,1
!         ARR(I)=SIN(2._RK*PI*FRE*REAL(I,RK))
!       END DO
!   END SELECT
!   ARR=ARR+0.1_RK
! 
!   CALL WINDOW(0,WIN)
!   AVE = SUM(WIN)
!   WRITE(*,*) ABS(FRE-FREQUENCY(AVE,WIN,ARR))
!   WRITE(*,*) ABS(FRE-FREQUENCY_BRENT(AVE,WIN,ARR))
!   
!   CALL WINDOW(2,WIN)
!   AVE = SUM(WIN)
!   WRITE(*,*) ABS(FRE-FREQUENCY(AVE,WIN,ARR))
!   WRITE(*,*) ABS(FRE-FREQUENCY_BRENT(AVE,WIN,ARR))
!   
! END PROGRAM EXAMPLE

MODULE SIGNAL

    IMPLICIT NONE
    
    PRIVATE
    
    INTEGER,PARAMETER  :: RK=SELECTED_REAL_KIND(15,307) ! REAL KIND
    INTEGER,PARAMETER  :: NUM=2048                      ! SIGNAL LENGTH (POWER OF TWO)
    INTEGER,PARAMETER  :: FLA=1                         ! COMPLEX SIGNAL FLAG (FLA=0 FOR REAL INPUT SIGNAL AND FLA=1 FOR COMPLEX INPUT SIGNAL)
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
    PUBLIC :: BRENT_MAX
    PUBLIC :: BRENT_TOL
    PUBLIC :: WINDOW
    PUBLIC :: FREQUENCY
    PUBLIC :: FREQUENCY_BRENT
    
    REAL(RK),DIMENSION(2*INTERPOLATION_POINTS+1) :: X_LIST
    REAL(RK),DIMENSION(2*INTERPOLATION_POINTS+1) :: Y_LIST
    
    ABSTRACT INTERFACE
        PURE FUNCTION GENERIC_FUNCTION(ARG)
            INTEGER,PARAMETER  :: RK=SELECTED_REAL_KIND(15,307)
            REAL(RK), INTENT(IN) :: ARG
            REAL(RK) :: GENERIC_FUNCTION
        END FUNCTION GENERIC_FUNCTION
    END INTERFACE
        
    CONTAINS
    
    ! DISCRETE FOURIER TRANSFORM (NRF77)
    PURE SUBROUTINE FFT(NUM,DIR,ARR)
        INTEGER,INTENT(IN) :: NUM
        INTEGER,INTENT(IN) :: DIR
        COMPLEX(RK),DIMENSION(NUM),INTENT(INOUT) :: ARR
        INTEGER :: N,I,J,M,LIM,STE
        REAL(RK) :: PIM,PRE,ANG,WR,WI,WPR,WPI,WX,MUL
        REAL(RK),DIMENSION(2*NUM) :: DAT
        DAT(1:2*NUM:2)=REAL(ARR,RK)
        DAT(2:2*NUM:2)=REAL(ARR*CMPLX(0.0_RK,-1.0_RK,RK),RK)
        N=2*NUM
        J=1
        MUL = REAL(DIR,RK)*TWO_PI
        DO I=1,N,2
            IF(J.GT.I)THEN
                PRE=DAT(J)
                PIM=DAT(J+1)
                DAT(J)=DAT(I)
                DAT(J+1)=DAT(I+1)
                DAT(I)=PRE
                DAT(I+1)=PIM
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
                    PRE=WR*DAT(J)-WI*DAT(J+1)
                    PIM=WR*DAT(J+1)+WI*DAT(J)
                    DAT(J)=DAT(I)-PRE
                    DAT(J+1)=DAT(I+1)-PIM
                    DAT(I)=DAT(I)+PRE
                    DAT(I+1)=DAT(I+1)+PIM
                ENDDO
                WX=WR
                WR=WR*WPR-WI*WPI+WR
                WI=WI*WPR+WX*WPI+WI
            ENDDO
            LIM=STE
            GOTO 2
        ENDIF
        ARR=CMPLX(DAT(1:2*NUM:2),DAT(2:2*NUM:2),RK)
    END SUBROUTINE FFT
    
    ! DISCRETE (LINEAR) FRACTIONAL FOURIER TRANSFORM
    PURE SUBROUTINE FFRFT(NUM,PAR,ARR)
        INTEGER,INTENT(IN) :: NUM
        REAL(RK),INTENT(IN) :: PAR
        COMPLEX(RK),DIMENSION(NUM),INTENT(INOUT) :: ARR
        COMPLEX(RK),DIMENSION(2*NUM) :: ONE,TWO
        INTEGER :: I
        REAL(RK) :: FAC,MUL
        FAC=PAR*PI/REAL(NUM,RK)
        DO I=1,NUM,1
            MUL=FAC*REAL(I-1)**2
            ONE(I)=ARR(I)*EXP(CMPLX(0.0_RK,MUL,RK))
            TWO(I)=EXP(CMPLX(0.0_RK,-MUL,RK))
        END DO
        ONE(NUM+1:2*NUM:1)=0._RK 
        DO I=NUM+1,2*NUM,1
            TWO(I)=EXP(CMPLX(0.0_RK,-FAC*REAL(I-1-2*NUM,RK)**2,RK))
        END DO
        CALL FFT(2*NUM,+1,ONE)
        CALL FFT(2*NUM,+1,TWO)
        ONE=ONE*TWO
        CALL FFT(2*NUM,-1,ONE)
        ARR=ONE(1:NUM:1)/REAL(2*NUM,RK)
        DO I=1,NUM,1
            ARR(I)=ARR(I)*EXP(CMPLX(0.0_RK,FAC*REAL(I-1)**2,RK))
        END DO
    END SUBROUTINE FFRFT
    
    ! FACTORIAL
    INTEGER PURE RECURSIVE FUNCTION FACTORIAL(N) RESULT(M)
        INTEGER,INTENT(IN) :: N
        IF(N==0) THEN
            M=1
        ELSE
            M=N*FACTORIAL(N-1)
        END IF
    END FUNCTION FACTORIAL
    
    ! WINDOW
    PURE SUBROUTINE WINDOW(ORD,ARR) 
        INTEGER,INTENT(IN):: ORD
        COMPLEX(RK),INTENT(OUT),DIMENSION(NUM) :: ARR
        INTEGER :: I
        REAL(RK) :: MUL
        MUL=2._RK**ORD*REAL(FACTORIAL(ORD),RK)**2/REAL(FACTORIAL(2*ORD),RK)
        DO I=1,NUM,1
            ARR(I)=MUL*(1._RK+COS(TWO_PI*(REAL(I-1,RK)/REAL(NUM,RK)-0.5_RK)))**ORD
        END DO
    END SUBROUTINE WINDOW
    
    ! FREQUENCY SEARCH (LARGEST AMPLITUDE)
    PURE REAL(RK) FUNCTION FREQUENCY(WINDOW_TOTAL,WINDOW,SIGNAL)
        COMPLEX(RK),INTENT(IN) :: WINDOW_TOTAL
        COMPLEX(RK),INTENT(IN),DIMENSION(NUM) :: WINDOW
        COMPLEX(RK),INTENT(IN),DIMENSION(NUM) :: SIGNAL
        COMPLEX(RK),DIMENSION(NUM) :: ARR
        COMPLEX(RK),DIMENSION(NUM) :: FOU
        REAL(RK),DIMENSION(NUM) :: SPE
        REAL(RK) :: FAC
        INTEGER :: FST,CND
        INTEGER :: I
        ! REMOVE MEAN & APPLY WINDOW
        ARR=(SIGNAL-SUM(WINDOW*SIGNAL)/WINDOW_TOTAL)*WINDOW
        ! FFT
        FOU=ARR
        CALL FFT(NUM,1,FOU)
        FST=MAXLOC(ABS(FOU(1:NUM/(2-FLA):1)),1)
        FAC=TWO_PI*REAL(FST-2,RK)/REAL(NUM,RK)
        ! MODULATE SIGNAL
        DO I=1,NUM,1
            ARR(I)=ARR(I)*EXP(CMPLX(0.0_RK,FAC*REAL(I-1,RK),RK))
        END DO
        ! FFRFT
        CALL FFRFT(NUM,2.0_RK/REAL(NUM,RK),ARR)
        SPE=LOG10(ABS(ARR)+1.E-16_RK)
        CND=MAXLOC(SPE,1)
        ! PARABOLA
        FREQUENCY=REAL(CND,RK)-0.5_RK+(SPE(-1+CND)-SPE(CND))/(SPE(-1+CND)-2._RK*SPE(CND)+SPE(1+CND))
        ! RESULT
        FREQUENCY=(REAL(FST,RK)-2._RK+2._RK*(FREQUENCY-1._RK)/REAL(NUM,RK))/REAL(NUM,RK)    
    END FUNCTION FREQUENCY
    
    ! POLINT (NRF90)
    PURE SUBROUTINE POLINT(XA,YA,X,Y,DY)
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
        NS=MINLOC(ABS(HO),1)
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
    END SUBROUTINE POLINT
    
    ! RATINT (NRF90)
    PURE SUBROUTINE RATINT(XA,YA,X,Y,DY)
        REAL(RK),DIMENSION(:),INTENT(IN) :: XA
        REAL(RK),DIMENSION(:),INTENT(IN) :: YA
        REAL(RK),INTENT(IN)  :: X
        REAL(RK),INTENT(OUT) :: Y,DY
        INTEGER :: N,M,NS,I
        REAL(RK), DIMENSION(SIZE(XA)) :: C,D,DD,H,T
        REAL(RK), PARAMETER :: EPS = 1.0E-25_RK
        N=SIZE(XA)
        H=XA-X
        NS=MINLOC(ABS(H),1)
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
    END SUBROUTINE RATINT
    
    ! INTERPOLATING FUNCTION
    PURE FUNCTION INTERPOLATION(ARG)
        REAL(RK),INTENT(IN) :: ARG
        REAL(RK) :: INTERPOLATION
        REAL(RK) :: X,Y,DY
        X = ARG
        SELECT CASE(INTERPOLATION_CASE)
            CASE(1)
                CALL POLINT(X_LIST,Y_LIST,X,Y,DY)
            CASE(2)
                CALL RATINT(X_LIST,Y_LIST,X,Y,DY)
            CASE DEFAULT
                CALL POLINT(X_LIST,Y_LIST,X,Y,DY)
        END SELECT
        INTERPOLATION = Y
    END FUNCTION INTERPOLATION
  
    ! BRENT (NRF90)
    PURE FUNCTION BRENT(AX,BX,CX,FUN)
        REAL(RK),INTENT(IN) :: AX
        REAL(RK),INTENT(IN) :: BX
        REAL(RK),INTENT(IN) :: CX
        PROCEDURE(GENERIC_FUNCTION) :: FUN
        REAL(RK),DIMENSION(2) :: BRENT
        INTEGER,PARAMETER :: ITMAX=BRENT_MAX
        REAL(RK), PARAMETER :: CGOLD=0.3819660112501052_RK
        REAL(RK), PARAMETER :: ZEPS=1.0E-3_RK*EPSILON(AX)
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
                BRENT=[XMIN,FX]
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
                CALL SHFT(V,W,X,U)
                CALL SHFT(FV,FW,FX,FU)
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
        
        PURE SUBROUTINE SHFT(A,B,C,D)
            REAL(RK),INTENT(OUT) :: A
            REAL(RK),INTENT(INOUT) :: B
            REAL(RK),INTENT(INOUT) :: C
            REAL(RK),INTENT(IN) :: D
            A=B
            B=C
            C=D
        END SUBROUTINE SHFT
        
    END FUNCTION BRENT
    
    ! FREQUENCY SEARCH (BRENT) (LARGEST AMPLITUDE) (NOT PURE)
    REAL(RK) FUNCTION FREQUENCY_BRENT(WINDOW_TOTAL,WINDOW,SIGNAL)
        COMPLEX(RK),INTENT(IN) :: WINDOW_TOTAL
        COMPLEX(RK),INTENT(IN),DIMENSION(NUM) :: WINDOW
        COMPLEX(RK),INTENT(IN),DIMENSION(NUM) :: SIGNAL
        COMPLEX(RK),DIMENSION(NUM) :: ARR
        COMPLEX(RK),DIMENSION(NUM) :: FOU
        REAL(RK),DIMENSION(NUM) :: SPE
        REAL(RK) :: FAC
        INTEGER :: FST,CND
        INTEGER :: I
        REAL(RK),DIMENSION(2) :: OUTPUT
        REAL(RK) :: AX,BX,CX
        ! REMOVE MEAN & APPLY WINDOW
        ARR=(SIGNAL-SUM(WINDOW*SIGNAL)/WINDOW_TOTAL)*WINDOW
        ! FFT
        FOU=ARR
        CALL FFT(NUM,1,FOU)
        FST=MAXLOC(ABS(FOU(1:NUM/(2-FLA):1)),1)
        FAC=TWO_PI*REAL(FST-2,RK)/REAL(NUM,RK)
        ! MODULATE SIGNAL
        DO I=1,NUM,1
            ARR(I)=ARR(I)*EXP(CMPLX(0.0_RK,FAC*REAL(I-1,RK),RK))
        END DO
        ! FFRFT
        CALL FFRFT(NUM,2.0_RK/REAL(NUM,RK),ARR)
        SPE=LOG10(ABS(ARR)+1.E-16_RK)
        CND=MAXLOC(SPE,1)
        ! BRENT
        X_LIST=[(REAL(I+CND,RK),I=-INTERPOLATION_POINTS,INTERPOLATION_POINTS,1)]
        Y_LIST=SPE(CND-INTERPOLATION_POINTS:CND+INTERPOLATION_POINTS)
        Y_LIST=1.0_RK/Y_LIST
        AX = X_LIST(INTERPOLATION_POINTS+0)
        BX = X_LIST(INTERPOLATION_POINTS+1)
        CX = X_LIST(INTERPOLATION_POINTS+2)
        OUTPUT = BRENT(AX,BX,CX,INTERPOLATION)
        FREQUENCY_BRENT = OUTPUT(1)
        ! RESULT
        FREQUENCY_BRENT=(REAL(FST,RK)-2._RK+2._RK*(FREQUENCY_BRENT-1._RK)/REAL(NUM,RK))/REAL(NUM,RK)    
    END FUNCTION FREQUENCY_BRENT
    
END MODULE SIGNAL

!     ABSTRACT INTERFACE
!         PURE SUBROUTINE TRANSFORMATION(STATE)
!             INTEGER,PARAMETER  :: RK=SELECTED_REAL_KIND(15,307)
!             REAL(RK),DIMENSION(:), INTENT(INOUT) :: STATE
!         END SUBROUTINE TRANSFORMATION
!     END INTERFACE
!     
!     ! FMA (LARGEST AMPLITUDE)
!     PURE SUBROUTINE FMA(AVE,WIN,MAP,INI,PAR,CAN,ANS)
!         COMPLEX(RK),INTENT(IN) :: AVE
!         COMPLEX(RK),INTENT(IN),DIMENSION(NUM) :: WIN
!         PROCEDURE(TRANSFORMATION) :: MAP
!         REAL(RK),DIMENSION(:),INTENT(INOUT) :: INI
!         INTEGER,INTENT(IN) :: PAR
!         INTEGER,INTENT(IN) :: CAN
!         REAL(RK),DIMENSION(:),INTENT(OUT) :: ANS
!         COMPLEX(RK),ALLOCATABLE,DIMENSION(:,:) :: ORB
!         COMPLEX(RK),ALLOCATABLE,DIMENSION(:) :: SIG
!         REAL(RK),ALLOCATABLE,DIMENSION(:) :: FRE
!         REAL(RK),ALLOCATABLE,DIMENSION(:) :: MEA, DEV
!         REAL(RK) :: MUL
!         INTEGER :: I
!         ALLOCATE(ORB(SIZE(INI),2*NUM))
!         ALLOCATE(SIG(NUM))
!         ALLOCATE(FRE(CAN))
!         ALLOCATE(MEA(CAN/2))
!         ALLOCATE(DEV(CAN/2))
!         MUL=REAL(FLA,RK)
!         ANS=0.0_RK
!         ANS(1:SIZE(INI):1) = INI
!         DO I=1,2*NUM,1
!             CALL MAP(INI)
!             ORB(:,I) = INI
!         END DO
!         ANS(1)=SUM(ORB(1,:))
!         IF(ANS(1)==REAL(2*NUM,RK))THEN        
!             DO I=1,CAN,2
!                 SIG=ORB(1+PAR+I,1:NUM:1)+MUL*CMPLX(0.0_RK,1.0_RK)*ORB(1+PAR+I+1,1:NUM:1)
!                 FRE(I)=FREQUENCY(AVE,WIN,SIG)
!                 SIG=ORB(1+PAR+I,NUM+1:2*NUM:1)+MUL*CMPLX(0.0_RK,1.0_RK)*ORB(1+PAR+I+1,NUM+1:2*NUM:1)
!                 FRE(I+1)=FREQUENCY(AVE,WIN,SIG)
!             END DO
!             MEA = 0.5_RK*(FRE(1:CAN:2)+FRE(2:CAN:2))
!             DEV = FRE(1:CAN:2)-FRE(2:CAN:2)
!             ANS = [ANS(1:SIZE(INI):1),MEA,LOG10(SUM(DEV**2)+1.E-16_RK)]
!         END IF        
!     END SUBROUTINE FMA

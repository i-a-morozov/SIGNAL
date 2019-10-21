! SIGNAL,2019,I.A.MOROZOV@INP.NSK.SU
! SIMPLIFIED FORTRAN IMPLEMENTATION OF S$FREQUENCY[] FUNCTION (PARABOLIC INTERPOLATION OF REFINED SPECTRA) AND S$NAFF[] FUNCTION

! MAIN FUNCTION(S):

! <OUT> = FREQUENCY_(<AVE>,<WIN>,<SIG>)
! <OUT> - FREQUENCY ESTIMATION (REAL)
! <AVE> - WINDOW SUM (REAL)
! <WIN> - WINDOW (REAL) WITH LENGTH = NUM (MUST BE POWER OF TWO)
! <SIG> - INPUT SIGNAL (REAL) WITH LENGTH = 2*NUM (MUST BE POWER OF TWO)

! DECOMPOSITION_(<AVE>,<WIN>,<SIG>,<NUM>,<FRE>,<COS>,<SIN>)
! <AVE> - WINDOW SUM (REAL)
! <WIN> - WINDOW (REAL) WITH LENGTH = NUM (MUST BE POWER OF TWO)
! <SIG> - INPUT SIGNAL (REAL) WITH LENGTH = 2*NUM (MUST BE POWER OF TWO)
! <NUM> - NUMBER OF HARMONICS (INTEGER)
! <FRE> - LIST OF ESTIMATED FREQUENCIES (REAL ARRAY WITH DIMENSION <NUM>)
! <COS> - LIST OF ESTIMATED COS AMPLITUDES (REAL ARRAY WITH DIMENSION <NUM>)
! <SIN> - LIST OF ESTIMATED SIN AMPLITUDES (REAL ARRAY WITH DIMENSION <NUM>)

! SEE EXAMPLES AT THE END OF FILE 

MODULE SIGNAL

    USE ISO_FORTRAN_ENV    

    IMPLICIT NONE
    
    PRIVATE
    
    INTEGER,PARAMETER  :: RK=REAL64                     ! REAL KIND
    INTEGER,PARAMETER  :: NUM=512                       ! SIGNAL LENGTH (POWER OF TWO)
    INTEGER,PARAMETER  :: FLA=1                         ! COMPLEX SIGNAL FLAG (0/1 FOR REAL/COMPLEX INPUT)
    
    REAL(RK),PARAMETER :: PI=3.141592653589793238460_RK ! PI
    REAL(RK),PARAMETER :: TWO_PI=2.0_RK*PI              ! 2*PI
    
    PUBLIC :: NUM
    PUBLIC :: FLA
    PUBLIC :: RK
    PUBLIC :: PI
    PUBLIC :: TWO_PI
    PUBLIC :: WINDOW_
    PUBLIC :: FREQUENCY_
    PUBLIC :: DECOMPOSITION_
                
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
    
    PURE SUBROUTINE DECOMPOSITION_(AVE,WIN,SIG,ITR,FRE,COS_AMP,SIN_AMP)
        !$ACC ROUTINE SEQ
        REAL(RK),INTENT(IN) :: AVE
        REAL(RK),DIMENSION(NUM),INTENT(IN) :: WIN
        REAL(RK),DIMENSION(2*NUM),INTENT(IN) :: SIG
        INTEGER,INTENT(IN) :: ITR
        REAL(RK),DIMENSION(ITR),INTENT(OUT) :: FRE
        REAL(RK),DIMENSION(ITR),INTENT(OUT) :: COS_AMP
        REAL(RK),DIMENSION(ITR),INTENT(OUT) :: SIN_AMP
        INTEGER :: I
        REAL(RK),DIMENSION(2*NUM) :: ORB
        REAL(RK),DIMENSION(NUM) :: RAN
        REAL(RK),DIMENSION(2*NUM) :: DEL
        ORB=SIG-AVE
        RAN=REAL([(I,I=1,NUM,1)],RK)
        DO I=1,ITR,1
            FRE(I)=FREQUENCY_(AVE,WIN,ORB)
            DEL(1:2*NUM:2)=+COS(TWO_PI*FRE(I)*RAN)
            DEL(2:2*NUM:2)=-SIN(TWO_PI*FRE(I)*RAN)
            COS_AMP(I)=SUM(SIG(1:2*NUM:2)*DEL(1:2*NUM:2)*WIN)+SUM(SIG(2:2*NUM:2)*DEL(2:2*NUM:2)*WIN)
            SIN_AMP(I)=SUM(SIG(2:2*NUM:2)*DEL(1:2*NUM:2)*WIN)-SUM(SIG(1:2*NUM:2)*DEL(2:2*NUM:2)*WIN)
            COS_AMP(I)=COS_AMP(I)/REAL(NUM,RK)
            SIN_AMP(I)=SIN_AMP(I)/REAL(NUM,RK)
            ORB(1:2*NUM:2)=ORB(1:2*NUM:2)-COS_AMP(I)*DEL(1:2*NUM:2)+SIN_AMP(I)*DEL(2:2*NUM:2)
            ORB(2:2*NUM:2)=ORB(2:2*NUM:2)-COS_AMP(I)*DEL(2:2*NUM:2)-SIN_AMP(I)*DEL(1:2*NUM:2)
        END DO
        IF(FLA==0)THEN
            COS_AMP = 2.0_RK*COS_AMP
            SIN_AMP = 2.0_RK*SIN_AMP
        END IF
    END SUBROUTINE DECOMPOSITION_
    
END MODULE SIGNAL

! ! EXAMPLE-01 (NUM=2**9,FLA=0/1) REAL/COMPLEX SIGNAL FREQUENCY
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
! END PROGRAM EXAMPLE

! ! EXAMPLE-02 (NUM=2**9,FLA=0/1) REAL/COMPLEX SIGNAL FREQUENCY (OPENMP) 
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

! ! EXAMPLE-03 (NUM=2**9,FLA=0/1) REAL/COMPLEX SIGNAL FREQUENCY (OPENACC) 
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

! ! EXAMPLE-04 (NUM=2**12,FLA=0) REAL SIGNAL DECOMPOSITION
! ! gfortran -o example -O3 -ffast-math -march=native signal.f90 example.f90
! ! example.f90:
! PROGRAM EXAMPLE
!     USE SIGNAL
!     IMPLICIT NONE
!     REAL(RK),PARAMETER :: F1=1._RK*0.123456789_RK
!     REAL(RK),PARAMETER :: F2=2._RK*0.123456789_RK
!     REAL(RK),PARAMETER :: F3=3._RK*0.123456789_RK
!     REAL(RK),PARAMETER :: F4=4._RK*0.123456789_RK
!     REAL(RK),PARAMETER :: A1=1.0_RK
!     REAL(RK),PARAMETER :: B1=1.E-1_RK
!     REAL(RK),PARAMETER :: A2=1.E-1_RK
!     REAL(RK),PARAMETER :: B2=5.E-2_RK
!     REAL(RK),PARAMETER :: A3=1.E-3_RK
!     REAL(RK),PARAMETER :: B3=0.0_RK
!     REAL(RK),PARAMETER :: A4=0.0_RK
!     REAL(RK),PARAMETER :: B4=1.E-4_RK
!     REAL(RK),DIMENSION(2*NUM) :: SIG
!     REAL(RK),DIMENSION(NUM)   :: WIN
!     REAL(RK)                  :: AVE
!     INTEGER                   :: I
!     REAL(RK),DIMENSION(4)     :: FRE, COS_AMP, SIN_AMP
!     SIG(1:2*NUM:2)=REAL([(I,I=1,NUM,1)],RK)
!     SIG(1:2*NUM:2)=A1*COS(TWO_PI*F1*SIG(1:2*NUM:2))+B1*SIN(TWO_PI*F1*SIG(1:2*NUM:2)) + &
!                    A2*COS(TWO_PI*F2*SIG(1:2*NUM:2))+B2*SIN(TWO_PI*F2*SIG(1:2*NUM:2)) + &
!                    A3*COS(TWO_PI*F3*SIG(1:2*NUM:2))+B3*SIN(TWO_PI*F3*SIG(1:2*NUM:2)) + &
!                    A4*COS(TWO_PI*F4*SIG(1:2*NUM:2))+B4*SIN(TWO_PI*F4*SIG(1:2*NUM:2))
!     SIG(2:2*NUM:2)=0.0_RK
!     CALL WINDOW_(4,WIN)
!     AVE = SUM(WIN)
!     FRE = 0.0_RK
!     COS_AMP=0.0_RK
!     SIN_AMP=0.0_RK
!     CALL DECOMPOSITION_(AVE,WIN,SIG,4,FRE,COS_AMP,SIN_AMP)
!     DO I=1,4,1
!         WRITE(*,*) FRE(I),COS_AMP(I),SIN_AMP(I)
!     END DO
! END PROGRAM EXAMPLE


! ! EXAMPLE-05 (NUM=2**12,FLA=1) COMPLEX SIGNAL DECOMPOSITION
! ! gfortran -o example -O3 -ffast-math -march=native signal.f90 example.f90
! ! example.f90:
! PROGRAM EXAMPLE
!     USE SIGNAL
!     IMPLICIT NONE
!     REAL(RK),PARAMETER :: F1=1._RK*0.123456789_RK
!     REAL(RK),PARAMETER :: F2=2._RK*0.123456789_RK
!     REAL(RK),PARAMETER :: F3=3._RK*0.123456789_RK
!     REAL(RK),PARAMETER :: F4=4._RK*0.123456789_RK
!     REAL(RK),PARAMETER :: A1=1.0_RK
!     REAL(RK),PARAMETER :: B1=1.E-1_RK
!     REAL(RK),PARAMETER :: A2=1.E-1_RK
!     REAL(RK),PARAMETER :: B2=5.E-2_RK
!     REAL(RK),PARAMETER :: A3=1.E-3_RK
!     REAL(RK),PARAMETER :: B3=0.0_RK
!     REAL(RK),PARAMETER :: A4=0.0_RK
!     REAL(RK),PARAMETER :: B4=1.E-4_RK
!     REAL(RK),DIMENSION(2*NUM) :: SIG
!     REAL(RK),DIMENSION(NUM)   :: WIN
!     REAL(RK)                  :: AVE
!     INTEGER                   :: I
!     REAL(RK),DIMENSION(4)     :: FRE, COS_AMP, SIN_AMP    
!     SIG(1:2*NUM:2)=REAL([(I,I=1,NUM,1)],RK)
!     SIG(1:2*NUM:2)=A1*COS(TWO_PI*F1*SIG(1:2*NUM:2))+B1*SIN(TWO_PI*F1*SIG(1:2*NUM:2)) + &
!                    A2*COS(TWO_PI*F2*SIG(1:2*NUM:2))+B2*SIN(TWO_PI*F2*SIG(1:2*NUM:2)) + &
!                    A3*COS(TWO_PI*F3*SIG(1:2*NUM:2))+B3*SIN(TWO_PI*F3*SIG(1:2*NUM:2)) + &
!                    A4*COS(TWO_PI*F4*SIG(1:2*NUM:2))+B4*SIN(TWO_PI*F4*SIG(1:2*NUM:2))
!     
!     SIG(2:2*NUM:2)=REAL([(I,I=1,NUM,1)],RK)
!     SIG(2:2*NUM:2)=B1*COS(TWO_PI*F1*SIG(2:2*NUM:2))-A1*SIN(TWO_PI*F1*SIG(2:2*NUM:2)) + &
!                    B2*COS(TWO_PI*F2*SIG(2:2*NUM:2))-A2*SIN(TWO_PI*F2*SIG(2:2*NUM:2)) + &
!                    B3*COS(TWO_PI*F3*SIG(2:2*NUM:2))-A3*SIN(TWO_PI*F3*SIG(2:2*NUM:2)) + &
!                    B4*COS(TWO_PI*F4*SIG(2:2*NUM:2))-A4*SIN(TWO_PI*F4*SIG(2:2*NUM:2))                
!     CALL WINDOW_(4,WIN)
!     AVE = SUM(WIN)
!     FRE = 0.0_RK
!     COS_AMP=0.0_RK
!     SIN_AMP=0.0_RK
!     CALL DECOMPOSITION_(AVE,WIN,SIG,4,FRE,COS_AMP,SIN_AMP)
!     DO I=1,4,1
!         WRITE(*,*) FRE(I),COS_AMP(I),SIN_AMP(I)
!     END DO
! END PROGRAM EXAMPLE

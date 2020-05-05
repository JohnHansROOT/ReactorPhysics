C THIS PROGRAM SOLVE FOR THE TRANSIENT, TWO-DIMENSIONAL CENTER TEMPERATURE IN A SQUARE SLAB
C AUTHOR: WENBIN HAN

C VARIABLE NAME   UNIT    MEANING
C     L           M       WIDTH OF THE SQUEARE SLAB
C     Q           W/M^3   VOLUMETRIC HEAT SOURCE
C     T_W         ℃      TEMPERATURE OF THE RIGHT WALL
C     T_INF       ℃      TEMPERATURE OF THE LEFT FLUID
C     FLUX        W/M^2   HEAT FLUX OF THE UPPER WALL
C     LAMBDA      W/MK    THERMAL CONDUCTIVITY OF THE SQUARE
C     ALPHA       W/M^2K  CONVECTION COEFFICIENT
C     N           1       NUMBER OF MESHES       
C     H           M       MESH SIZE
C     ERR_LIMIT   1       ERROR LIMIT FOR CONVERGENCE
C     OMEGE       1       RELAXATION FACTOR OF SOR METHOD
C     T           ℃      TEMPERATURE MATRIX
C     RHO         KG/M^3  MATERIAL DENSITY
C     CP          J/(KG K) SPECIFIC HEAT
C     THETA               IMPLICIT FACTOR
C     FO                  FOURIER NUMBER


C BASIC PARAMETERS MODULE
      MODULE BASIC_PARAMETERS
          REAL*8 L,Q,T_W,T_INF,FLUX,LAMBDA,ALPHA,RHO,CP,FO,DTIME,THETA,H
          REAL*8 QH2L,HRL,HTRL,FHL,R_X
          INTEGER N
      END MODULE
      
C PARAMETERS FOR STEP CALCULATIONS
      MODULE STEP_CAL_PARAMETERS
          USE BASIC_PARAMETERS
          REAL*8,ALLOCATABLE::A(:),B(:),C(:),D(:),X(:),ERR(:)
          REAL*8 ERR_LIMIT,ERR_MAX,OMEGA
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
CCCC             MAIN PROGRAM                 CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM MAIN
      USE BASIC_PARAMETERS
      USE STEP_CAL_PARAMETERS
      REAL*8 T_INIT,T_CENTER,T_CENTER_NEW,DT_CENTER,T_SS,SS_LIMIT
      REAL*8,ALLOCATABLE::T(:,:),T_NEW(:,:)
      INTEGER*8::K=0
C BASIC PHYSICAL PARAMETERS
      L = 0.4
      Q = 2.E4
      T_W = 40.
      T_INF = 30.
      FLUX = 1200.
      LAMBDA = 2.5
      ALPHA = 20.
      !RHO = 2000.
      RHO = 1000.
      CP = 400.
      CP = 800.
C NUMBER OF MESHES!
      N = 80
      H = L / N
      ALLOCATE(A(N))
      ALLOCATE(B(N))
      ALLOCATE(C(N))
      ALLOCATE(D(N))
      ALLOCATE(X(N))
      ALLOCATE(ERR(N))
      ALLOCATE(T(N,N))
      ALLOCATE(T_NEW(N,N))
C SOME INTERMIDIATE CONSTANTS
      R_X = H / (2. * LAMBDA) + 1. / ALPHA 
      QH2L = Q * H ** 2 / LAMBDA
      HRL = H / (R_X * LAMBDA)
      HTRL = H * T_INF / (R_X * LAMBDA)
      FHL = FLUX * H / LAMBDA
      
C RELAXATION FACTOR, CHANGEBLE!
      OMEGA = 1.6
C IMPLICIT FACTOR, CHANGEBLE!
      THETA = 1.0
C FOURIER NUMBER, TIME STEP
      FO = 0.1
C ERROR LIMIT FOR STEP CALCULATIONS
      ERR_LIMIT = 1.D-7
      
C STEADY STATE LIMIT
      !SS_LIMIT = 1.D-6
      !DT_CENTER = 1.0
      
C CENTER TEMPERATURE DIFFERENCE LIMIT
      DT_LIMIT = 0.999
      DT_CENTER = 0.
      
C TIME SIZE
      DTIME = FO * RHO * CP * H ** 2 / LAMBDA
C STEADY STATE TEMPERETURE
      T_SS = 300.3717
      !T_SS = 300.37
C INITIANL TEMPERATURE
      T_INIT = 20.d0

      T(:,:)=T_INIT
      T_NEW(:,:)=T_INIT
      T_CENTER = (T(N/2,N/2)+T(N/2+1,N/2+1)+T(N/2,N/2+1)+T(N/2+1,N/2))/4


      OPEN(UNIT = 10, FILE = "T_CENTER.TXT")
      WRITE(10,200) "Steps","Time","Center temperature"
      !WRITE(*,200) "Steps","Time","Center temperature"
      WRITE(10,300) K, K * DTIME, T_CENTER
      !WRITE(*,300) K, K * DTIME, T_CENTER

C START ITERATIONS
      !DO WHILE(DT_CENTER .GT. SS_LIMIT)
      DO WHILE(DT_CENTER .LT. DT_LIMIT)
          CALL STEP_CALCULATE(T,T_NEW)
          T_CENTER_NEW = (T_NEW(N/2,N/2)+T_NEW(N/2+1,N/2+1)+T_NEW(N/2,
     1    N/2+1)+T_NEW(N/2+1,N/2))/4
          K = K + 1
          !DT_CENTER = ABS(T_CENTER_NEW - T_CENTER)
          DT_CENTER = (T_CENTER_NEW-T_INIT)/(T_SS-T_INIT)
          T=T_NEW
          T_CENTER = T_CENTER_NEW
          WRITE(10,300) K, K * DTIME, T_CENTER
          !WRITE(*,300) K, K * DTIME, T_CENTER
      ENDDO
      !WRITE(*,300) K, K * DTIME, T_CENTER
200   FORMAT(A10,A10,A20)
300   FORMAT(I12,F12.2,F12.6)
      END
      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC      TIME STEP CALCULATION      CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE STEP_CALCULATE(T,T_NEW)
      USE BASIC_PARAMETERS
      USE STEP_CAL_PARAMETERS
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 T(N,N),T_NEW(N,N)
      INTEGER::J
      
      ERR_MAX = 1.
          DO WHILE(ERR_MAX .GT. ERR_LIMIT)
              J = 1
              CALL COEFFICIENT_ASSIGN(T,T_NEW,J)
              CALL TDMA(A,B,C,D,N,X)
              ERR(J) = MAXVAL(ABS((X-T_NEW(:,J))/X))
              T_NEW(:,J)=OMEGA * X + (1- OMEGA) * T_NEW(:,J)
              DO J = 2, N - 1
                  CALL COEFFICIENT_ASSIGN(T,T_NEW,J)
                  CALL TDMA(A,B,C,D,N,X)
                  ERR(J) = MAXVAL(ABS((X-T_NEW(:,J))/X))
                  T_NEW(:,J) = OMEGA * X + (1 - OMEGA) * T_NEW(:,J)
              ENDDO
              J=N
              CALL COEFFICIENT_ASSIGN(T,T_NEW,J)
              CALL TDMA(A,B,C,D,N,X)
              ERR(J) = MAXVAL(ABS((X-T_NEW(:,J))/X))
              T_NEW(:,J) = OMEGA * X + (1 - OMEGA) * T_NEW(:,J)
              ERR_MAX = MAXVAL(ERR(:))
          ENDDO
          RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC     ASSIGN THE VETORS OF A B C D         CCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE COEFFICIENT_ASSIGN(T,T_NEW,J)
      USE BASIC_PARAMETERS
      USE STEP_CAL_PARAMETERS
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 T(N,N),T_NEW(N,N)
      INTEGER::I,J
      
      IF(J .EQ. 1) THEN
          I = 1
          A(I)=1/FO+THETA*(2+HRL)
          B(I)=THETA
          C(I)=0.
          D(I)=THETA*T_NEW(I,J+1)+(1/FO-(1-THETA)*(2+HRL))*T(I,J)+
     1    (1-THETA)*(T(I+1,J)+T(I,J+1))+HTRL+QH2L
          DO I = 2, N-1
              A(I)=1/FO+3*THETA
              B(I)=THETA
              C(I)=THETA
              D(I)=THETA*T_NEW(I,J+1)+(1/FO-3*(1-THETA))*T(I,J)+
     1        (1-THETA)*(T(I+1,J)+T(I-1,J)+T(I,J+1))+QH2L
          ENDDO
          I = N
          A(I)=1/FO+4*THETA
          B(I)=0.
          C(I)=THETA
          D(I)=THETA*T_NEW(I,J+1)+(1/FO-4*(1-THETA))*T(I,J)+(1-THETA)*
     1    (T(I-1,J)+T(I,J+1))+2*T_W+QH2L
      ELSE IF(J .EQ. N) THEN
          I = 1
          A(I)=1/FO+THETA*(2+HRL)
          B(I)=THETA
          C(I)=0.
          D(I)=THETA*T_NEW(I,J-1)+(1/FO-(1-THETA)*(2+HRL))*T(I,J)+
     1    (1-THETA)*(T(I+1,J)+T(I,J-1))+HTRL+FHL+QH2L
          DO I = 2, N-1
              A(I)=1/FO+3*THETA
              B(I)=THETA
              C(I)=THETA
              D(I)=THETA*T_NEW(I,J-1)+(1/FO-3*(1-THETA))*T(I,J)+
     1        (1-THETA)*(T(I+1,J)+T(I-1,J)+T(I,J-1))+QH2L+FHL
          ENDDO
          I = N
          A(I)=1/FO+4*THETA
          B(I)=0.
          C(I)=THETA
          D(I)=THETA*T_NEW(I,J-1)+(1/FO-4*(1-THETA))*T(I,J)+(1-THETA)*
     1    (T(I-1,J)+T(I,J-1))+2*T_W+FHL+QH2L
      ELSE
          I = 1
          A(I)=1/FO+THETA*(3+HRL)
          B(I)=THETA
          C(I)=0.
          D(I)=THETA*(T_NEW(I,J-1)+T_NEW(I,J+1))+(1/FO-(1-THETA)*(3+HRL)
     1    )*T(I,J)+(1-THETA)*(T(I+1,J)+T(I,J-1)+T(I,J+1))+HTRL+QH2L
          DO I = 2, N-1
              A(I)=1/FO+4*THETA
              B(I)=THETA
              C(I)=THETA
              D(I)=THETA*(T_NEW(I,J-1)+T_NEW(I,J+1))+(1/FO-4*(1-THETA))
     1        *T(I,J)+(1-THETA)*(T(I+1,J)+T(I-1,J)+T(I,J-1)+T(I,J+1))
     1        +QH2L
          ENDDO
          I = N
          A(I)=1/FO+5*THETA
          B(I)=0.
          C(I)=THETA
          D(I)=THETA*(T_NEW(I,J-1)+T_NEW(I,J+1))+(1/FO-5*(1-THETA))*
     1    T(I,J)+(1-THETA)*(T(I-1,J)+T(I,J-1)+T(I,J+1))+2*T_W+QH2L
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC       TRI-DIAGONAL SOLUTION ALGORITHM            CCCC
CCCC FOR THE SOLUTION OF: -C*T(I-1)+A*T(I)-B*T(I+1)=D CCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TDMA(A,B,C,D,N,X)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 A(N),B(N),C(N),D(N),X(N),B0(N),D0(N)
      B0(1) = B(1) / A(1)
      D0(1) = D(1) / A(1)
      DO I = 2 , N
          B0(I) = B(I) / (A(I) - C(I) * B0(I-1))
          D0(I) = (D(I) + C(I) * D0(I-1)) / (A(I) - C(I) * B0(I-1))
      ENDDO
      X(N) = D0(N)
      DO I = N - 1 , 1 , -1
          X(I) = B0(I) * X(I+1) + D0(I)
      ENDDO
      RETURN
      END
C THIS PROGRAM SOLVE FOR THE STEADY-STATE, TWO-DIMENSIONAL TEMPERATURE DISTRIBUTION IN A SQUARE SLAB
C AUTHOR: WENBIN HAN

C VARIABLE NAME   UNIT    MEANING
C     L           M       WIDTH OF THE SQUEARE SLAB
C     Q           W/M^3   VOLUMETRIC HEAT SOURCE
C     T_W         ¡æ      TEMPERATURE OF THE RIGHT WALL
C     T_INF       ¡æ      TEMPERATURE OF THE LEFT FLUID
C     FLUX        W/M^2   HEAT FLUX OF THE UPPER WALL
C     LAMBDA      W/MK    THERMAL CONDUCTIVITY OF THE SQUARE
C     ALPHA       W/M^2K  CONVECTION COEFFICIENT
C     N           1       NUMBER OF MESHES       
C     H           M       MESH SIZE
C     ERR_LIMIT   1       ERROR LIMIT FOR CONVERGENCE
C     OMEGE       1       RELAXATION FACTOR OF SOR METHOD
C     ITERATIONS  1       NUMBER OF ITERATIONS
C     T           ¡æ      TEMPERATURE MATRIX
C     TOTAL_HEAT_TRANSFER W/M TOTAL HEAT TRANSFER FROM EACH WALL, POSITIVE IN

      PROGRAM MAIN
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 L, Q, T_W, T_INF, FLUX, LAMBDA, ALPHA, H, ERR_LIMIT, OMEGA
      REAL*8 TOTAL_HEAT_TRANSFER(4),LIMIT(10)
      REAL*8,ALLOCATABLE::T(:,:),A(:),B(:),C(:),D(:),X(:),ERR(:)
      INTEGER*8 ITERATIONS
C SOME PARAMETERS
      L = 0.4
      Q = 2.E4
      T_W = 40.
      T_INF = 30.
      FLUX = 1200.
      LAMBDA = 2.5
      ALPHA = 20.
C NUMBER OF MESHES, CHANGEBLE!
      N = 80
      H = L / N
C MATRIX PREPARATION
      ALLOCATE(T(N,N))
      ALLOCATE(A(N))
      ALLOCATE(B(N))
      ALLOCATE(C(N))
      ALLOCATE(D(N))
      ALLOCATE(X(N))
      ALLOCATE(ERR(N))
C SOME INTERMIDIATE CONSTANTS
      R_X = H / (2. * LAMBDA) + 1. / ALPHA 
      QH2L = Q * H ** 2 / LAMBDA
      HRL = H / (R_X * LAMBDA)
      HTRL = H * T_INF / (R_X * LAMBDA)
      FHL = FLUX * H / LAMBDA
C INITIAL TEMPERATURE
      T_INIT = 0.d0

C MAXIMUM RELATIVE ERROR, INITIAL SET TO 1.

C RELAXATION FACTOR
      OMEGA = 1.9398

C INITIAL TEMPERATURE DISTRIBUTION
      T(:,:)=T_INIT
      DO I = 1, 10
          LIMIT(I) = 1./ (10.**I)
      ENDDO
      
C START ITERATIONS
      OPEN(UNIT=10,FILE='ERR_LIMITS.TXT')
      WRITE(*,"(A9,A18)") "ERR LIMIT", "CENTER TEMPERATURE"
      DO I = 1, 10, 1
          ERR_LIMIT = LIMIT(I)
          ITERATIONS = 0
          ERR_MAX = 1.
          T(:,:)=T_INIT
          DO WHILE(ERR_MAX .GT. ERR_LIMIT)
              CALL COEFFICIENT_ASSIGN(2.+HRL,3.d0,4.d0,1.d0,1.d0,0.d0
     1            ,0.d0 ,1.d0,1.d0,HTRL+QH2L,QH2L,2*T_W+QH2L,A,B,C,D,N)
              CALL TDMA(A,B,C,T(:,2)+D,N,X)
              ERR(1) = MAXVAL(ABS((X-T(:,1))/X))
              T(:,1) = OMEGA * X + (1 - OMEGA) * T(:,1)
              CALL COEFFICIENT_ASSIGN(3.+HRL,4.d0,5.d0,1.d0,1.d0,0.d0
     1            ,0.d0,1.d0,1.0d0,HTRL+QH2L,QH2L,2*T_W+QH2L,A,B,C,D,N)
              DO J = 2, N - 1 
                  CALL TDMA(A,B,C,T(:,J+1)+T(:,J-1)+D,N,X)
                  ERR(J) = MAXVAL(ABS((X-T(:,J))/X))
                  T(:,J) = OMEGA * X + (1 - OMEGA) * T(:,J)
              ENDDO
              CALL COEFFICIENT_ASSIGN(2.+HRL,3.d0,4.d0,1.d0,1.d0,0.d0
     1            ,0.d0,1.d0,1.d0,HTRL+FHL+QH2L,FHL+QH2L,2*T_W+FHL+QH2L
     1            ,A,B,C,D,N)
              CALL TDMA(A,B,C,T(:,N-1)+D,N,X)
              ERR(N) = MAXVAL(ABS((X-T(:,N))/X))
              T(:,N) = OMEGA * X + (1 - OMEGA) * T(:,N)
              ERR_MAX = MAXVAL(ERR(:))
              ITERATIONS = ITERATIONS + 1
          ENDDO
          WRITE(*,"(D12.6,F12.6)") ERR_LIMIT, (T(N/2,N/2)+T(N/2+1,N/2
     1    +1)+T(N/2,N/2+1)+T(N/2+1,N/2))/4
      ENDDO
      END
      

      
C ASSIGN THE VETORS OF A B C D
      SUBROUTINE COEFFICIENT_ASSIGN(A1,AX,AN,B1,BX,BN,C1,CX,CN,D1,DX,DN
     &,A,B,C,D,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(N),B(N),C(N),D(N)
      A(1) = A1
      B(1) = B1
      C(1) = C1
      D(1) = D1
      DO I = 2, N - 1
          A(I) = AX
          B(I) = BX
          C(I) = CX
          D(I) = DX
      ENDDO
      A(N) = AN
      B(N) = BN
      C(N) = CN
      D(N) = DN
      RETURN
      END
      
      
C TRI-DIAGONAL SOLUTION ALGORITHM
C FOR THE SOLUTION OF: -C*T(I-1)+A*T(I)-B*T(I+1)=D
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
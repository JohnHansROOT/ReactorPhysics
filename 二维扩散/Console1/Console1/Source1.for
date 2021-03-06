CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCC         PARAMETERS MODULE               CCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCC GROUP AND CROSS SECTION PARAMETERS
      MODULE PARAMETERS0
      REAL*8 CHI(2,3),NUSIGMA_F(2,3),SIGMA_A(2,3),SIGMA_S_1_TO_2(2,3)
      REAL*8 D(2,3),SIGMA_R(2,3)
      DATA CHI /1.0,0.0,1.0,0.0,0.0,0.0/
      DATA NUSIGMA_F /0.0085,0.1851,0.006,0.150,0.0,0.0/
      DATA SIGMA_A /0.0121,0.121,0.010,0.100,0.0004,0.020/
      DATA SIGMA_S_1_TO_2 /0.,0.0241,0.,0.016,0.,0.0493/
      DATA D /1.267,0.354,1.280,0.400,1.130,0.166/
      END MODULE
CCCCCC SOME COMMON PARAMTERS
      MODULE PARAMETERS1
      REAL*8 DELTAX,DELTAY,INNER_lIMIT,OMEGA,KEFF
      REAL*8,ALLOCATABLE::PHI1(:,:),PHI2(:,:)
      REAL*8,ALLOCATABLE::SOURCE(:,:),SOURCE_NEW(:,:)
      INTEGER*8 M,N
      END MODULE
CCCCCC TDMA COEFFICIENTS
      MODULE PARAMETERS2
      REAL*8,ALLOCATABLE::A(:),B(:),C(:),H(:),X(:)
      END MODULE
CCCCCC GEOMETRAY PARAMTERS
      MODULE PARAMETERS3
      REAL*8 XMAX,X_REF_IN,YMAX,Y_REF_UP_IN,Y_JUNC,Y_REF_DOWN_IN
      END MODULE
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC      MAIN PROGRAM         CCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM MAIN
      USE PARAMETERS0
      USE PARAMETERS1
      USE PARAMETERS2
      USE PARAMETERS3
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 OUTTER_ITERATIONS
      REAL*8 OUTTER_LIMIT,KEFF_NEW,ERR_KEFF
CCCCCC GEOMETRY PARAMETERS
      XMAX=65.; X_REF_IN=50.; YMAX=120.; Y_REF_UP_IN=105.; Y_JUNC=55.;
      Y_REF_DOWN_IN=15.;
CCCCCC MESH GRID
      M = 130; N = 240;
      DELTAX = XMAX / M; DELTAY = YMAX / N
CCCCCC GET SIGMA_R
      SIGMA_R(1,:)=SIGMA_A(1,:)+SIGMA_S_1_TO_2(2,:)
      SIGMA_R(2,:)=SIGMA_A(2,:)

      ALLOCATE(PHI1(M,N))
      ALLOCATE(PHI2(M,N))
      ALLOCATE(SOURCE(M,N))
      ALLOCATE(SOURCE_NEW(M,N))
      ALLOCATE(A(M))
      ALLOCATE(B(M))
      ALLOCATE(C(M))
      ALLOCATE(H(M))
      ALLOCATE(X(M))
CCCCCC RELAXATION FACTOR
      OMEGA = 1.9
CCCCCC INNER AND OUTTER LIMIT
      INNER_LIMIT=1E-5; OUTTER_LIMIT=1.E-10;
CCCCCC INITIAL KEFF AND FLUX
      KEFF=1.0
      PHI1(:,:)=1.D12
      PHI2(:,:)=1.D12
CCCCCC INITIAL SOURCE
      CALL UPDATE_SOURCE()
      SOURCE(:,:)=SOURCE_NEW(:,:)
CCCCCC NUMBER OF OUTTER ITERATIONS 
      OUTTER_ITERATIONS = 0
      ERR_KEFF=1.
      DO WHILE(ERR_KEFF .GT. OUTTER_LIMIT)
CCCCCC CALCULATE THE FIRST GROUP
          CALL CAL_GROUP(1,PHI1)
CCCCCC CALCULATE THE SECOND GROUP
          CALL CAL_GROUP(2,PHI2)
CCCCCC UPDATE THE SOURCE
          CALL UPDATE_SOURCE()
CCCCCC CALCULATE THE NEW KEFF
          KEFF_NEW=SUM(SOURCE_NEW(:,:))/SUM(SOURCE(:,:))*KEFF
CCCCCC CALCULATE THE KEFF ERROR
          ERR_KEFF=ABS(1.-KEFF/KEFF_NEW)
          KEFF=KEFF_NEW
          OUTTER_ITERATIONS = OUTTER_ITERATIONS + 1
          SOURCE(:,:)=SOURCE_NEW(:,:)
          WRITE(*,*) OUTTER_ITERATIONS,KEFF
      ENDDO
CCCCCC PRITE OUT THE FLUX
      OPEN(UNIT=10,FILE='PHI1OUT.TXT')
      OPEN(UNIT=20,FILE='PHI2OUT.TXT')
      OPEN(UNIT=30,FILE='PHIOUT.TXT')
      DO J = N, 1, -1
              WRITE(10,100),(PHI1(I,J),I=1,M)
              WRITE(20,100),(PHI2(I,J),I=1,M)
              WRITE(30,100),(PHI1(I,J)+PHI2(I,J),I=1,M)
      ENDDO
100   FORMAT(<M>(1X,E14.6))
      END PROGRAM
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCC       CALCULATE FOR EACH GROUP        CCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CAL_GROUP(ID_GROUP,PHI)
      USE PARAMETERS1
      USE PARAMETERS2
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PHI(M,N),PHI_NEW(M,N),ERR_PHI
      INTEGER*8 J
      ERR_PHI=1.0
      PHI_NEW(:,:)=PHI(:,:)
      DO WHILE(ERR_PHI .GT. INNER_LIMIT)
          DO J=1,N
              CALL COEFFICIENT_ASSIGN(J,ID_GROUP,PHI_NEW)
              CALL TDMA(A,B,C,H,M,X)
              PHI_NEW(:,J)=OMEGA *X + (1-OMEGA)*PHI_NEW(:,J)
          ENDDO
          ERR_PHI=SUM(ABS(PHI_NEW(:,:)-PHI(:,:)))/SUM(ABS(PHI_NEW(:,:)))
          PHI(:,:)=PHI_NEW(:,:)
      ENDDO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC          UPDATE SOURCE       CCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE UPDATE_SOURCE()
      USE PARAMETERS1
      USE PARAMETERS0
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*8 I,J
      DO I=1,M
          DO J=1,N
              CALL FIND_REGION(I,J,ID_REGION)
              SOURCE_NEW(I,J)=NUSIGMA_F(1,ID_REGION)*PHI1(I,J)+
     &        NUSIGMA_F(2,ID_REGION)*PHI2(I,J)
          ENDDO
      ENDDO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC      ASSIGN THE COEFFICIENTS     CCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE COEFFICIENT_ASSIGN(J,ID_GROUP,PHI)
      USE PARAMETERS0
      USE PARAMETERS1
      USE PARAMETERS2
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PHI(M,N),DE,DN,DS,DW,E,F,G
      INTEGER*8::I,J
      
      IF(J .EQ. 1) THEN
          I = 1
          CALL FIND_REGION(I,J,ID_REGION)
          CALL FIND_REGION(I+1,J,ID_REGION_E)
          CALL FIND_REGION(I,J+1,ID_REGION_N)
          DE=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_E))/2.
          DN=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_N))/2.
          DS=D(ID_GROUP,ID_REGION)
          B(I)=DE*DELTAY/DELTAX
          C(I)=0.
          E=DN*DELTAX/DELTAY
          F=0.
          A(I)=B(I)+C(I)+E+F+2*DS*DELTAX/(4*DS+DELTAY)+
     &    SIGMA_R(ID_GROUP,ID_REGION)*DELTAX*DELTAY
          G=(SIGMA_S_1_TO_2(ID_GROUP,ID_REGION)*PHI1(I,J)+
     &    SOURCE(I,J)*CHI(ID_GROUP,ID_REGION)/KEFF)  *DELTAX*DELTAY
          H(I)=E*PHI(I,J+1)+G
          DO I = 2, M-1
              CALL FIND_REGION(I,J,ID_REGION)
              CALL FIND_REGION(I+1,J,ID_REGION_E)
              CALL FIND_REGION(I-1,J,ID_REGION_W)
              CALL FIND_REGION(I,J+1,ID_REGION_N)
              DE=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_E))/2.
              DN=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_N))/2.
              DW=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_W))/2.
              DS=D(ID_GROUP,ID_REGION)
              B(I)=DE*DELTAY/DELTAX
              C(I)=DW*DELTAY/DELTAX
              E=DN*DELTAX/DELTAY
              F=0.
              A(I)=B(I)+C(I)+E+F+2*DS*DELTAX/(4*DS+DELTAY)+
     &        SIGMA_R(ID_GROUP,ID_REGION)*DELTAX*DELTAY
              G=(SIGMA_S_1_TO_2(ID_GROUP,ID_REGION)*PHI1(I,J)
     &         +SOURCE(I,J)*CHI(ID_GROUP,ID_REGION)/KEFF)*DELTAX*DELTAY
              H(I)=E*PHI(I,J+1)+G
          ENDDO
          I = M
          CALL FIND_REGION(I,J,ID_REGION)
          CALL FIND_REGION(I-1,J,ID_REGION_W)
          CALL FIND_REGION(I,J+1,ID_REGION_N)
          DE=D(ID_GROUP,ID_REGION)
          DN=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_N))/2.
          DW=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_W))/2.
          DS=D(ID_GROUP,ID_REGION)
          B(I)=0.
          C(I)=DW*DELTAY/DELTAX
          E=DN*DELTAX/DELTAY
          F=0.
          A(I)=B(I)+C(I)+E+F+2*DS*DELTAX/(4*DS+DELTAY)+2*DE*DELTAY/
     &    (4*DE+DELTAX)+SIGMA_R(ID_GROUP,ID_REGION)*DELTAX*DELTAY
          G=(SIGMA_S_1_TO_2(ID_GROUP,ID_REGION)*PHI1(I,J)+
     &    SOURCE(I,J)*CHI(ID_GROUP,ID_REGION)/KEFF)*DELTAX*DELTAY
          H(I)=E*PHI(I,J+1)+G
      ELSE IF(J .EQ. N) THEN
          I = 1
          CALL FIND_REGION(I,J,ID_REGION)
          CALL FIND_REGION(I+1,J,ID_REGION_E)
          CALL FIND_REGION(I,J-1,ID_REGION_S)
          DE=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_E))/2.
          DS=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_S))/2.
          DN=D(ID_GROUP,ID_REGION)
          B(I)=DE*DELTAY/DELTAX
          C(I)=0.
          F=DS*DELTAX/DELTAY
          E=0.
          A(I)=B(I)+C(I)+E+F+2*DN*DELTAX/(4*DN+DELTAY)+
     &    SIGMA_R(ID_GROUP,ID_REGION)*DELTAX*DELTAY
          G=(SIGMA_S_1_TO_2(ID_GROUP,ID_REGION)*PHI1(I,J)+
     &    SOURCE(I,J)*CHI(ID_GROUP,ID_REGION)/KEFF)*DELTAX*DELTAY
          H(I)=F*PHI(I,J-1)+G
          DO I = 2, M-1
              CALL FIND_REGION(I,J,ID_REGION)
              CALL FIND_REGION(I+1,J,ID_REGION_E)
              CALL FIND_REGION(I-1,J,ID_REGION_W)
              CALL FIND_REGION(I,J-1,ID_REGION_S)
              DE=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_E))/2.
              DS=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_S))/2.
              DW=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_W))/2.
              DN=D(ID_GROUP,ID_REGION)
              B(I)=DE*DELTAY/DELTAX
              C(I)=DW*DELTAY/DELTAX
              F=DS*DELTAX/DELTAY
              E=0.
              A(I)=B(I)+C(I)+E+F+2*DN*DELTAX/(4*DN+DELTAY)+
     &        SIGMA_R(ID_GROUP,ID_REGION)*DELTAX*DELTAY
              G=(SIGMA_S_1_TO_2(ID_GROUP,ID_REGION)*PHI1(I,J)
     &         +SOURCE(I,J)*CHI(ID_GROUP,ID_REGION)/KEFF)*DELTAX*DELTAY
              H(I)=F*PHI(I,J-1)+G
          ENDDO
          I = M
          CALL FIND_REGION(I,J,ID_REGION)
          CALL FIND_REGION(I-1,J,ID_REGION_W)
          CALL FIND_REGION(I,J-1,ID_REGION_S)
          DE=D(ID_GROUP,ID_REGION)
          DS=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_S))/2.
          DW=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_W))/2.
          DN=D(ID_GROUP,ID_REGION)
          B(I)=0.
          C(I)=DW*DELTAY/DELTAX
          F=DS*DELTAX/DELTAY
          E=0.
          A(I)=B(I)+C(I)+E+F+2*DN*DELTAX/(4*DN+DELTAY)+2*DE*DELTAY/
     &    (4*DE+DELTAX)+SIGMA_R(ID_GROUP,ID_REGION)*DELTAX*DELTAY
          G=(SIGMA_S_1_TO_2(ID_GROUP,ID_REGION)*PHI1(I,J)+
     &    SOURCE(I,J)*CHI(ID_GROUP,ID_REGION)/KEFF)*DELTAX*DELTAY
          H(I)=F*PHI(I,J-1)+G
      ELSE
          I = 1
          CALL FIND_REGION(I,J,ID_REGION)
          CALL FIND_REGION(I+1,J,ID_REGION_E)
          CALL FIND_REGION(I,J+1,ID_REGION_N)
          CALL FIND_REGION(I,J-1,ID_REGION_S)
          DE=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_E))/2.
          DS=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_S))/2.
          DN=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_N))/2.
          B(I)=DE*DELTAY/DELTAX
          C(I)=0.
          F=DS*DELTAX/DELTAY
          E=DN*DELTAX/DELTAY
          A(I)=B(I)+C(I)+E+F+SIGMA_R(ID_GROUP,ID_REGION)*DELTAX*DELTAY
          G=(SIGMA_S_1_TO_2(ID_GROUP,ID_REGION)*PHI1(I,J)+
     &    SOURCE(I,J)*CHI(ID_GROUP,ID_REGION)/KEFF)*DELTAX*DELTAY
          H(I)=E*PHI(I,J+1)+F*PHI(I,J-1)+G
          DO I = 2, M-1
              CALL FIND_REGION(I,J,ID_REGION)
              CALL FIND_REGION(I+1,J,ID_REGION_E)
              CALL FIND_REGION(I-1,J,ID_REGION_W)
              CALL FIND_REGION(I,J+1,ID_REGION_N)
              CALL FIND_REGION(I,J-1,ID_REGION_S)
              DE=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_E))/2.
              DS=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_S))/2.
              DW=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_W))/2.
              DN=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_N))/2.
              B(I)=DE*DELTAY/DELTAX
              C(I)=DW*DELTAY/DELTAX
              F=DS*DELTAX/DELTAY
              E=DN*DELTAX/DELTAY
              A(I)=B(I)+C(I)+E+F+SIGMA_R(ID_GROUP,ID_REGION)
     &        *DELTAX*DELTAY
              G=(SIGMA_S_1_TO_2(ID_GROUP,ID_REGION)*PHI1(I,J)
     &        +SOURCE(I,J)*CHI(ID_GROUP,ID_REGION)/KEFF)*DELTAX*DELTAY
              H(I)=E*PHI(I,J+1)+F*PHI(I,J-1)+G
          ENDDO
          I = M
          CALL FIND_REGION(I,J,ID_REGION)
          CALL FIND_REGION(I-1,J,ID_REGION_W)
          CALL FIND_REGION(I,J+1,ID_REGION_N)
          CALL FIND_REGION(I,J-1,ID_REGION_S)
          DE=D(ID_GROUP,ID_REGION)
          DS=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_S))/2.
          DW=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_W))/2.
          DN=(D(ID_GROUP,ID_REGION)+D(ID_GROUP,ID_REGION_N))/2.
          B(I)=0.
          C(I)=DW*DELTAY/DELTAX
          F=DS*DELTAX/DELTAY
          E=DN*DELTAX/DELTAY
          A(I)=B(I)+C(I)+E+F+2*DE*DELTAY/(4*DE+DELTAX)
     &    +SIGMA_R(ID_GROUP,ID_REGION)*DELTAX*DELTAY
          G=(SIGMA_S_1_TO_2(ID_GROUP,ID_REGION)*PHI1(I,J)+
     &    SOURCE(I,J)*CHI(ID_GROUP,ID_REGION)/KEFF)*DELTAX*DELTAY
          H(I)=E*PHI(I,J+1)+F*PHI(I,J-1)+G
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC     FIND THE REGION             CCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FIND_REGION(I,J,ID_REGION)
      USE PARAMETERS1
      USE PARAMETERS3
      REAL*8 X,Y
      INTEGER*8 I,J
      INTEGER ID_REGION
      X=(I-0.5)*DELTAX
      Y=(J-0.5)*DELTAY
      IF((X .LT. X_REF_IN) .AND. (Y .GT. Y_REF_DOWN_IN) .AND. 
     &    (Y .LE. Y_JUNC)) THEN
          ID_REGION=1
      ELSE IF((X .LT. X_REF_IN) .AND. (Y .GT. Y_JUNC) .AND. 
     &    (Y .LT. Y_REF_UP_IN)) THEN
          ID_REGION=2
      ELSE
          ID_REGION=3
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC       TRI-DIAGONAL SOLUTION ALGORITHM            CCCC
CCCC FOR THE SOLUTION OF: -C*T(I-1)+A*T(I)-B*T(I+1)=H CCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TDMA(A,B,C,H,M,X)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 M
      REAL*8 A(M),B(M),C(M),H(M),X(M),B0(M),H0(M)
      B0(1) = B(1) / A(1)
      H0(1) = H(1) / A(1)
      DO I = 2 , M
          B0(I) = B(I) / (A(I) - C(I) * B0(I-1))
          H0(I) = (H(I) + C(I) * H0(I-1)) / (A(I) - C(I) * B0(I-1))
      ENDDO
      X(M) = H0(M)
      DO I = M - 1 , 1 , -1
          X(I) = B0(I) * X(I+1) + H0(I)
      ENDDO
      RETURN
      END
C ������������S4�������㵥��һά����ƽ�����
C AUTHOR: ���ı�
C DATE: 2019/11/10

      
C �������Բ���
      MODULE PARAMETERS
      REAL*8 THICKNESS
      REAL*8 SIGMA_T,SIGMA_S,NU_SIGMA_F,Q
      REAL*8 OMEGA(4),MU(4)
      DATA OMEGA/0.3478548451,0.6521451549,0.6521451549,0.3478548451/
      DATA MU/-0.8611363116,-0.3399810436,0.3399810436,0.8611363116/
      INTEGER*4 N,K
      END MODULE

C ������
      PROGRAM MAIN
      USE PARAMETERS
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8,ALLOCATABLE::FLUX(:,:),FLUX_OLD(:,:),FLUX_ALL(:)
      REAL*8,ALLOCATABLE::L_OPT(:),FLUX_ALL_OLD(:),SOURCE(:)
      REAL*8 DELTA_X,ERR_LIMIT,ERR_EFF,ERR_FLUX,KEFF,KEFF_OLD,TOTAL_FLUX
      LOGICAL::FLAG=.FALSE.
      
      SIGMA_T = 0.050
      SIGMA_S = 0.030
      NU_SIGMA_F = 0.0225
      Q = 0.
      THICKNESS = 66.0053
      !THICKNESS = 1e6
      N = 4
      K = 1000
      
      ALLOCATE(FLUX(N,K+1))
      ALLOCATE(FLUX_OLD(N,K+1))
      ALLOCATE(FLUX_ALL(K))
      ALLOCATE(FLUX_ALL_OLD(K))
      ALLOCATE(SOURCE(K))
      ALLOCATE(L_OPT(N))
      
      DELTA_X = THICKNESS / K
C ��ѧ����
      L_OPT(:) = SIGMA_T * DELTA_X /MU(:)
C ��ʼԴ�ֲ�
      SOURCE(:) = 1.
C ��ʼͨ���ܶ�
      FLUX_ALL_OLD(:) = 1.
C ��ʼ��Ч��ֵϵ��
      KEFF_OLD=1.
      KEFF=KEFF_OLD
C ��ʼ���������
      ERR_EFF = 1.
      ERR_FLUX = 1.
      ERR_LIMIT = 1.E-12
      
C ������ʼ
      DO WHILE((ERR_EFF .GT. ERR_LIMIT) .AND. (ERR_FLUX .GT. ERR_LIMIT))
C ��ձ߽�����
          FLUX(1:(N/2),K+1) = 0.
C \mu_m < 0����
          DO M = 1 , N / 2
              DO J = K , 1, -1 
                  FLUX(M,J) = (2. + L_OPT(M)) / (2. - L_OPT(M)) * 
     &            FLUX(M,J + 1) - 2 * L_OPT(M) / (2 - L_OPT(M)) * 
     &            SOURCE(J) / SIGMA_T
              ENDDO
          ENDDO
C ����߽�����
          DO M = 1, N / 2
              FLUX(N - M + 1,1) = FLUX(M , 1)
          ENDDO
C \mu_m > 0����
          DO M = N / 2 + 1, N
              DO J = 1 , K, 1
                  FLUX(M,J + 1) = (2. - L_OPT(M)) / (2. + L_OPT(M)) * 
     &            FLUX(M,J) + 2 * L_OPT(M) / (2 + L_OPT(M)) * SOURCE(J)
     &            / SIGMA_T
              ENDDO
          ENDDO
C ����ͨ���ܶ�
          DO J = 1, K , 1
              FLUX_ALL(J) = 0.5 * SUM(OMEGA(:) * (FLUX(:,J) + 
     &        FLUX(:,J+1)))
          ENDDO
C flag���״ε�����־����ʼͨ���ֲ�������keff�����
          IF(FLAG) THEN
              KEFF = SUM(NU_SIGMA_F*FLUX_ALL(:)) / 
     &        SUM(NU_SIGMA_F*FLUX_ALL_OLD(:)) * KEFF_OLD
              ERR_EFF = ABS((KEFF - KEFF_OLD)/KEFF_OLD)
              KEFF_OLD = KEFF
              ERR_FLUX = MAXVAL(ABS(FLUX_ALL(:)-FLUX_ALL_OLD(:))
     &        /FLUX_ALL_OLD(:))
          ENDIF
          FLAG = .TRUE.
          FLUX_OLD(:,:) = FLUX(:,:)
          FLUX_ALL_OLD(:) = FLUX_ALL(:)
c ����Դ��
          CALL SOURCE_CAL(FLUX,SOURCE,KEFF)
      ENDDO
      
C ���������ļ�����Ļ
C �����ͨ���ܶ�
      OPEN(UNIT=10,FILE='ANGULAR_FLUX.TXT')
      WRITE(10,400) "x","flux_1","flux_2","flux_3","flux_4"
      DO J = 1, K , 1
          WRITE(10,300) (J-1)*DELTA_X, FLUX(:,J)
      ENDDO
C ���ͨ���ܶ�
      OPEN(UNIT=20,FILE='FLUX_ALL.TXT')
      WRITE(20,100) "x","flux"
      DO J = 1, K , 1
          WRITE(20,200) (J-1)*DELTA_X, FLUX_ALL(J)/FLUX_ALL(1)
      ENDDO
C �����Ч������ֵϵ��
      WRITE(*,*) "keff=",KEFF
      
100   FORMAT(2(A10))
200   FORMAT(2(F12.6)) 
300   FORMAT(<N+1>(F12.6))
400   FORMAT(<N+1>(A10))
      END PROGRAM

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCC ����Դ���      CCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SOURCE_CAL(FLUX,S,KEFF)
      USE PARAMETERS
      REAL*8 FLUX(N,K+1),S(K),KEFF
      INTEGER J
      DO J = 1, K, 1
          S(J)=(SIGMA_S+NU_SIGMA_F/KEFF)/4. * SUM(OMEGA(:) * (FLUX(:,J)+
     &    FLUX(:,J+1)))+ Q
      ENDDO
      RETURN
      END SUBROUTINE
      SUBROUTINE FUNC_RK(NEQ,X,Y,DELTAY,H)
C ODE TO BE INTEGRATED BY FOURTH-ORDER RK METHOD
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 :: Y(NEQ),DELTAY(10)
      DELTAY(1) = Y(2) * H
      DELTAY(2) = Y(3) * H
      DELTAY(3) = Y(4) * H
      DELTAY(4) = Y(5) * H
      DELTAY(5) = 120. * H
      RETURN
      END

      SUBROUTINE RUNGE(NEQ,X,Y,H)
C FOURTH-ORDER RUNGE KUTTA INTEGRATION ROUTINE
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y(NEQ),DELTAY1(10),DELTAY2(10),DELTAY3(10),
     1DELTAY4(10),Z(10)
      
C FIRST STEP
      DO I = 1, NEQ
          Z(I) = Y(I)
      ENDDO
      CALL FUNC_RK(NEQ,X,Z,DELTAY1,H)
C SECOND STEP
      DO I = 1, NEQ
          Z(I) = Y(I) + DELTAY1(I) / 2.
      ENDDO
      X = X + H / 2.
      CALL FUNC_RK(NEQ,X,Z,DELTAY2,H)
C THIRD STEP
      DO I = 1, NEQ
          Z(I) = Y(I) + DELTAY2(I) / 2.
      ENDDO
      CALL FUNC_RK(NEQ,X,Z,DELTAY3,H)
C FOURTH STEP
      DO I = 1, NEQ
          Z(I) = Y(I) + DELTAY3(I)
      ENDDO
      X = X + H / 2.
      CALL FUNC_RK(NEQ,X,Z,DELTAY4,H)
C SUMMARY
      DO I = 1, NEQ
          Y(I) = Y(I) + (DELTAY1(I) + 2. * DELTAY2(I) + 2. * DELTAY3(I)+
     1    DELTAY4(I)) / 6.
      ENDDO
      RETURN
      END



C TEST RK ROUTINE
      PROGRAM TEST
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 Y(5),DY(5)
      
      WRITE(*,*) "ENTER NMAX:"
      READ(*,*) NMAX
      XMAX = 1.
      H = XMAX / NMAX
      NEQ = 5
C RK INTEGRATION LOOP
      Y(1) = 3.
      Y(2) = 1.
      Y(3) = 2.
      y(4) = 6.
      Y(5) = 24.
      X = 0.
      ERROR = 0.
      
      DO I = 1, NMAX
          CALL RUNGE(NEQ,X,Y,H)
C COMPARE RESULTS TO ANALYTUCAL EXPRESSION
          Y_ANAL = X ** 5 + X ** 4 + X ** 3 + X ** 2 + X + 3
          ERROR = ERROR + (Y(1) - Y_ANAL) ** 2
          WRITE(*,*) X,Y(1),Y(2),Y(3),Y(4),Y(5)
      ENDDO
      WRITE(*,*) "RUNGE KUTTA TOTAL ERROR:",SQRT(ERROR)
      STOP
      END
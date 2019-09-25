C SHOOTING ROUTINE - FIND ROOTS OF A POLYNOMIAL
      PROGRAM FIND_ROOTS
      IMPLICIT REAL*8(A-H,O-Z)
      A1 = 0
      CALL FUNC(A1,B1)
      WRITE(*,*) A1,B1
      A2 = -5
      CALL FUNC(A2,B2)
      WRITE(*,*) A2, B2
      B3 = 1.
      DO WHILE(ABS(B3) .GT. 1.D-6)
          !A3 = A1 - B1 * (A2-A1)/(B2-B1)
          A3 = (A1 + A2) / 2
          CALL FUNC(A3,B3)
          WRITE(*,*) A3,B3
          IF(B3*B1 .GT. 0) THEN
              B1 = B3
              A1 = A3
          ELSE
              B2 = B3
              A2 = A3
          ENDIF
      ENDDO
      STOP
      END

      SUBROUTINE FUNC(X,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      Y = X ** 2 + 2 * X - 3
      RETURN
      END
      
      
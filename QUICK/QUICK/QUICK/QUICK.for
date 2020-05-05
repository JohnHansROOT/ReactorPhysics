      PROGRAM MAIN
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 T_W,T_E,T_WW,P
      REAL*8 T_P,T_LE,T_LW
      REAL*8 T_P_OLD,T_LE_OLD,T_LW_OLD
      REAL*8 ERR,ERR_LIMIT
      T_W = 10.
      T_E = 6.
      T_WW = 12.1
      P = 2.5
      
      T_P_OLD = 0.
      T_LE_OLD = 0.
      T_LW_OLD = 0.
      T_LE = T_LE_OLD
      T_LW = T_LW_OLD
      
      ERR_LIMIT = 1.E-6
      ERR = 1.0
      OPEN(UNIT = 10, FILE = "OUT.txt")
      WRITE(10,200) "T_w","T_P","T_e"
      DO WHILE(ERR .GT. ERR_LIMIT)
          T_P = (T_E+(1+P)*T_W+P*(T_LW-T_LE)+P*(T_P_OLD-T_W))/(2+P)
          WRITE(10,100) T_LW,T_P,T_LE
          CALL QUICK(T_W,T_P,T_E,T_LE)
          CALL QUICK(T_WW,T_W,T_P,T_LW)
          ERR = MAX(ABS(T_P-T_P_OLD)/T_P,ABS(T_LE-T_LE_OLD)/T_LE,
     &    ABS(T_LW-T_LW_OLD)/T_LW)
          T_P_OLD = T_P
          T_LW_OLD = T_LW
          T_LE_OLD = T_LE
100       FORMAT(3(F12.4))
200       FORMAT(3(A12))
      ENDDO
      
      END PROGRAM MAIN
      
      SUBROUTINE QUICK(T1,T2,T3,T)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 T1,T2,T3,T
      T = -1./8. * T1 + 3./4. * T2 + 3./8. * T3
      RETURN
      END SUBROUTINE QUICK
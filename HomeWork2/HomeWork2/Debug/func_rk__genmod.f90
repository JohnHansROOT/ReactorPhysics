        !COMPILER-GENERATED INTERFACE MODULE: Tue Sep 24 22:28:44 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FUNC_RK__genmod
          INTERFACE 
            SUBROUTINE FUNC_RK(NEQ,X,Y,DELTAY,H)
              INTEGER(KIND=4) :: NEQ
              REAL(KIND=8) :: X
              REAL(KIND=8) :: Y(NEQ)
              REAL(KIND=8) :: DELTAY(10)
              REAL(KIND=8) :: H
            END SUBROUTINE FUNC_RK
          END INTERFACE 
        END MODULE FUNC_RK__genmod

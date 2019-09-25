        !COMPILER-GENERATED INTERFACE MODULE: Tue Sep 24 22:28:44 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RUNGE__genmod
          INTERFACE 
            SUBROUTINE RUNGE(NEQ,X,Y,H)
              INTEGER(KIND=4) :: NEQ
              REAL(KIND=8) :: X
              REAL(KIND=8) :: Y(NEQ)
              REAL(KIND=8) :: H
            END SUBROUTINE RUNGE
          END INTERFACE 
        END MODULE RUNGE__genmod

      SUBROUTINE SWTN
!     APEX1501
!     THIS SUBPROGRAM CALCULATES SOIL WATER TENSION AS A FUNCTION OF
!     WATER CONTENT IN RELATION TO WILTING POINT & FIELD CAPACITY.
      USE PARM 
      DO J=1,NBSL(ISA)
          ISL=LID(J,ISA)
          IF(Z(ISL,ISA)>=.15)EXIT
      END DO
      XX=LOG10(S15(ISL,ISA))
      WTN=10.**(3.1761-1.6576*((LOG10(SWST(ISL,ISA))-XX)/(LOG10(FC(ISL,&
      ISA))-XX)))
      RETURN
      END
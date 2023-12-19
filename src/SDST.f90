      SUBROUTINE SDST(X,DG,DG1,X1,X2,I,ISA,MSL)
!     APEX1501
!     THIS SUBPROGRAM IS USED TO ESTIMATE MISSING SOIL DATA.
      DIMENSION X(MSL,ISA)
      IF(X(I,ISA)>0.)RETURN
      IF(I>1)THEN
          X(I,ISA)=X(I-1,ISA)*DG*EXP(-X2*DG)/DG1
      ELSE
          X(1,ISA)=X1
      END IF
      RETURN
      END
      SUBROUTINE HGASP(A,PT,Q1,RX)
!     APEX1501
!     THIS SUBPROGRAM SOLVES THE GREEN & AMPT INFILTRATION EQ ASSUMING
!     F1 IS INCREMENTED BY TOTAL RAIN DURING DT.
      USE PARM
      F1=PT-QVOL(IDO)
      ZI=SATK(ISA)*(SCN/F1+1.)
      IF(RX<=ZI)THEN
          Q1=0.
      ELSE
          Q1=A*(RX-ZI)/RX
      END IF
      QVOL(IDO)=QVOL(IDO)+Q1
!     WRITE(KW(1),2)ISA,IYR,MO,KDA,SATK(ISA),RX,ZI,F1,Q1,PT,QVOL(IDO)
!   2 FORMAT(1X,4I4,10F10.3)
      RETURN
      END
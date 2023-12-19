      SUBROUTINE HWTBL
!     APEX1501
!     THIS SUBPROGRAM SIMULATES WATER TABLE DYNAMICS AS A FUNCTION
!     OF RAIN AND EVAP.
      USE PARM
!     RTO=100.*GWST(ISA)/GWMX(ISA)
!     WTBL(ISA)=WTMX(ISA)+(WTMN(ISA)-WTMX(ISA))*RTO/(RTO+EXP(SCRP(26,1)-&
!     SCRP(26,2)*RTO))
      RTO=(SMRF(ISA)-SMEO(ISA))/SMEO(ISA)
      IF(RTO>0.)THEN
          XX=WTMN(ISA)
          X1=1.
      ELSE
          XX=WTMX(ISA)
          X1=PRMT(87)*(REAL(IDA)/REAL(ND))**PRMT(89)
      END IF
      X2=MIN(PRMT(88),ABS(RTO)*X1)
      IF(IMF==0)WTBL(ISA)=(WTBL(ISA)-(WTBL(ISA)-XX)*X2)

      SUM=0.
      TOT=0.
      IF(WTBL(ISA)<=Z(LID(NBSL(ISA),ISA),ISA))THEN
          XX=0.
          NN=0
          DO K=1,NBSL(ISA)
              ISL=LID(K,ISA)
              SUM=SUM+SWST(ISL,ISA)
              IF(WTBL(ISA)<=Z(ISL,ISA))THEN
                  IF(NN>0)THEN
                      SWST(ISL,ISA)=PO(ISL,ISA)
                  ELSE
                      NN=1
                      SWST(ISL,ISA)=(SWST(ISL,ISA)*(WTBL(ISA)-XX)+PO(ISL,ISA)*&
                      (Z(ISL,ISA)-WTBL(ISA)))/(Z(ISL,ISA)-XX)
                  END IF
              END IF
              TOT=TOT+SWST(ISL,ISA)
              XX=Z(ISL,ISA)
          END DO
      END IF
      XX=TOT-SUM
      QIN(MO,ISA)=QIN(MO,ISA)+XX
      XX=10.*XX*WSA(ISA)
      SMM(19,MO,ISA)=SMM(19,MO,ISA)+XX
      VAR(19,ISA)=XX
      RETURN
      END
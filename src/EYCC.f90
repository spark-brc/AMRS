      SUBROUTINE EYCC(ISX)
!     APEX1501
!     THIS SUBPROGRAM ESTIMATES THE USLE C FACTOR BASED ON PLANT POP &
!     BIOMASS & RESIDUE COVER
      USE PARM
      NN=NCP(IRO(ISX),ISX)
      SUM=0.
      TOT=0.
      DO I=1,NBSL(ISX)
          L=LID(I,ISX)
          IF(Z(L,ISX)>PMX(ISX))GO TO 3
          IF(I>1)TOT=TOT+RSD(L,ISX)
          DO K=1,NN
              SUM=SUM+RWT(L,JE(K,ISX),ISX)
          END DO
      END DO
      GO TO 4
    3 KK=LID(I-1,ISX)
      RTO=(PMX(ISX)-Z(KK,ISX))/(Z(L,ISX)-Z(KK,ISX))
      TOT=TOT+RTO*RSD(L,ISX)
      DO K=1,NN
          SUM=SUM+RTO*RWT(L,JE(K,ISX),ISX)
      END DO
    4 SUM=SUM/PMX(ISX)
      TOT=TOT/(PMX(ISX)-.01)
      FRUF=MIN(1.,EXP(-.026*(RRUF(ISX)-6.1)))
      IF(CVRS(ISX)<15.)THEN
          FRSD=EXP(-PRMT(46)*CVRS(ISX))
      ELSE
          FRSD=.0001
      END IF    
      X1=MAX(FGC(ISX),FGSL(ISX))
      FBIO=1.-X1*EXP(-PRMT(47)*CHMX(ISX))
      IF(IRP(ISX)>0)THEN
          FPPL=.9*(1.-CVP(ISX))+.1
      ELSE
          FPPL=1.
      END IF
      CVF(ISX)=FRSD*FBIO*FRUF*FPPL
      CVF(ISX)=CVF(ISX)*EXP(-.03*ROK(LID(1,ISX),ISX))
      RETURN
      END
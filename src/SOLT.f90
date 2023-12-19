      SUBROUTINE SOLT
!     APEX1501
!     THIS SUBPROGRAM ESTIMATES DAILY AVEAGE TEMPERATURE AT THE CENTER
!     OF EACH SOIL LAYER.
      USE PARM 
      DATA ITMP/0/
      DATA XLAG/.8/
      XLG1=1.-XLAG
      F=ABD(ISA)/(ABD(ISA)+686.*EXP(-5.63*ABD(ISA)))
      DP=1.+2.5*F
      WW=.356-.144*ABD(ISA)
      B=LOG(.5/DP)
      WC=.001*SW(ISA)/(WW*Z(LID(NBSL(ISA),ISA),ISA))
      F=MAX(.2,EXP(B*((1.-WC)/(1.+WC))**2))
      DD=F*DP
      IF(ITMP==0)THEN
          X4=.5*AMPX
          XZ=.5*(TMX(IRF(ISA))-TMN(IRF(ISA)))*ST0(ISA)/10.
          X2=TX
          X3=(1.-BCV(ISA))*X2+BCV(ISA)*STMP(LID(2,ISA),ISA)
          TG=.5*(X2+X3)
          DST0(ISA)=AVT+X4*COS((IDA-200)/PIT)
          X1=TX-DST0(ISA)
          DO J=1,NBSL(ISA)
              ISL=LID(J,ISA)
              X7=Z(ISL,ISA)/DD
              X5=EXP(-X7)
              STMP(ISL,ISA)=(AVT+X4*COS((IDA-200)/PIT-X7)*X5)+X1*X5
          END DO
      ELSE        
          X2=TX+.5*(TMX(IRF(ISA))-TMN(IRF(ISA)))*ST0(ISA)/30.
          X3=(1.-BCV(ISA))*X2+BCV(ISA)*STMP(LID(2,ISA),ISA)
          DST0(ISA)=.5*(X2+X3)
          ZZ=2.*DD
          XX=0.
          X1=AVT-DST0(ISA)
          DO J=1,NBSL(ISA)
              ISL=LID(J,ISA)
              ZD=(XX+Z(ISL,ISA))/ZZ
              F=ZD/(ZD+EXP(-.8669-2.0775*ZD))
              STMP(ISL,ISA)=XLAG*STMP(ISL,ISA)+XLG1*(F*X1+DST0(ISA))
              XX=Z(ISL,ISA)
          END DO
      END IF    
      RETURN
      END
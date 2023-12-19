      SUBROUTINE TBURN
!     APEX1501
!     THIS SUBPROGRAM BURNS ALL STANDING AND FLAT CROP RESIDUE.
      USE PARM
      DIMENSION JGO(MNC)
      !write(*,*) iyr,ida
      LD1=LID(1,ISA)
      X2=1.-PRMT(56)
      WSAX=WSA(ISA)
      ADD=0.
      SUM=0.
      IF(IGO(ISA)>0)THEN
          JGO=0
          JBG(ISA)=JBG(ISA)+1
          IF(JBG(ISA)>IGO(ISA))JBG(ISA)=1
          I=JBG(ISA)
          DO J=1,LC
              IF(KGO(J,ISA)>0)THEN
                  JGO(I)=KGO(J,ISA)
                  I=I+1
                  IF(I>IGO(ISA))I=1
              END IF    
          END DO
      ENDIF
      
      I1=0
      DO IN2=1,IGO(ISA)
         I1=I1+1
         JJK=JGO(I1)
          !IF(IDC(J)/=NDC(7).AND.IDC(J)/=NDC(8).AND.IDC(J)/=NDC(10))THEN
              RTO=MIN(.99,STL(JJK,ISA)/(DM(JJK,ISA)+1.E-10))
              X1=PRMT(56)*STL(JJK,ISA) 
              DM(JJK,ISA)=DM(JJK,ISA)-X1
              STL(JJK,ISA)=STL(JJK,ISA)-X1
              X3=PRMT(56)*UN1(JJK,ISA)*RTO
              UN1(JJK,ISA)=UN1(JJK,ISA)-X3
              SLAI(JJK,ISA)=MAX(.1,SLAI(JJK,ISA)*X2)
              CPHT(JJK,ISA)=MAX(.01,CPHT(JJK,ISA)*X2)
              HU(JJK,ISA)=MAX(.1,HU(JJK,ISA)*X2)
          !END IF
          X5=STD(JJK,ISA)*PRMT(56)
          STD(JJK,ISA)=STD(JJK,ISA)-X5
          ADD=ADD+420.*X5
          STDL(JJK,ISA)=STDL(JJK,ISA)*X2
          X1=STDN(JJK,ISA)*PRMT(56)
          STDN(JJK,ISA)=MAX(1.E-5,STDN(JJK,ISA)-X1)
          SUM=SUM+X1+X3
          X1=STDP(JJK,ISA)*PRMT(56)
          FOP(LD1,ISA)=FOP(LD1,ISA)+X1
          STDP(JJK,ISA)=STDP(JJK,ISA)-X1
      END DO
      FOP(LD1,ISA)=FOP(LD1,ISA)+STDOP(ISA)
      RSDM(LD1,ISA)=RSDM(LD1,ISA)*X2
      SMM(11,MO,ISA)=SMM(11,MO,ISA)+SWLT(ISA)
      WLS(LD1,ISA)=WLS(LD1,ISA)*X2
      WLM(LD1,ISA)=WLM(LD1,ISA)*X2
      X1=PRMT(56)*WLSN(LD1,ISA)
      WLSN(LD1,ISA)=WLSN(LD1,ISA)-X1
      X3=PRMT(56)*WLMN(LD1,ISA)
      WLMN(LD1,ISA)=WLMN(LD1,ISA)-X3
      SUM=SUM+X1+X3+STDON(ISA)
      WLSL(LD1,ISA)=WLSL(LD1,ISA)*X2
      X1=PRMT(56)*WLSC(LD1,ISA)
      WLSC(LD1,ISA)=WLSC(LD1,ISA)-X1
      X3=PRMT(56)*WLMC(LD1,ISA)
      WLMC(LD1,ISA)=WLMC(LD1,ISA)-X3
      WLSLC(LD1,ISA)=WLSLC(LD1,ISA)*X2
      WLSLNC(LD1,ISA)=WLSC(LD1,ISA)-WLSLC(LD1,ISA)
      ADD=ADD+X1+X3
      ADD=ADD+WSAX
      SMM(140,MO,ISA)=SMM(140,MO,ISA)+ADD
      VAR(140,ISA)=ADD
      SUM=SUM*WSAX
      SMM(82,MO,ISA)=SMM(82,MO,ISA)+SUM
      VAR(82,ISA)=SUM
      RSD(LD1,ISA)=.001*(WLS(LD1,ISA)+WLM(LD1,ISA))
      SWLT(ISA)=0.
      STDO(ISA)=0.
      STDOK(ISA)=0.
      STDON(ISA)=0.
      STDOP(ISA)=0.
      RETURN
      END
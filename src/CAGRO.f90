      SUBROUTINE CAGRO
!     APEX1501
!     THIS SUBPROGRAM CALCULATES THE DAILY INCREASE IN PLANT BIOMASS,
!     ROOT WEIGHT, AND YIELD BY ADJUSTING THE POTENTIAL VALUES WITH THE
!     ACTIVE STRESS CONSTRAINT.
      USE PARM
      !AD1=0.
      !DO K=1,NBSL(ISA)
          !ISL=LID(K,ISA)
          !AD1=AD1+ST(ISL,ISA)
      !END DO
      LD1=LID(1,ISA)
      X1=REG(JJK,ISA)*AJWA*SHRL
      RWL=RW(JJK,ISA)
      RGD=DDM(JJK)*X1
      X2=100.*HUI(JJK,ISA)
      AJH0=AJH1(JJK,ISA)
      IF(SCLM(3)>0.)X2=MIN(X2,SCLM(3))
      AJH1(JJK,ISA)=HI(JJK)*X2/(X2+EXP(SCRP(3,1)-SCRP(3,2)*X2))
      DHI=AJH1(JJK,ISA)-AJH0
      X3=TOPC(JJK)-TX
      IF(X3<0..AND.HUI(JJK,ISA)>.7)THEN
          XX=EXP(PRMT(30)*X3/TOPC(JJK))
          DHI=DHI*XX
      END IF    
      AJHI(JJK,ISA)=AJHI(JJK,ISA)+DHI
      X3=DM(JJK,ISA)-DDM(JJK)
      X4=MAX(1.E-5,X3+RGD)
      DM(JJK,ISA)=X4
      DM1(JJK,ISA)=DM1(JJK,ISA)+RGD
      X1=HUI(JJK,ISA)
      RF=MIN(.99,MAX(.2,RWPC(1,JJK)*(1.-X1)+RWPC(2,JJK)*X1))
      RW(JJK,ISA)=PRMT(97)*RW(JJK,ISA)+(1.-PRMT(97))*RF*DM(JJK,ISA)
	   DRW=RW(JJK,ISA)-RWL
      STL(JJK,ISA)=DM(JJK,ISA)-RW(JJK,ISA)
      IF(IDC(JJK)==NDC(7).OR.IDC(JJK)==NDC(8).OR.IDC(JJK)==NDC(10))THEN
          X1=MO*MO
	       FALF=STL(JJK,ISA)*X1*XMTU(JJK)
          SMM(109,MO,ISA)=SMM(109,MO,ISA)+FALF
          SLAI(JJK,ISA)=MAX(.05,SLAI(JJK,ISA)-FALF/STL(JJK,ISA)*SLAI(JJK,ISA))
          DM(JJK,ISA)=DM(JJK,ISA)-FALF
          STL(JJK,ISA)=STL(JJK,ISA)-FALF
          X11=FALF*1000.
          X5=CNY(JJK)*X11
          X6=CPY(JJK)*X11
          X10=CKY(JJK)*X11
          IF(FALF>1.E-5)CALL NCNSTD(FALF,X5,LD1)
          FOP(LD1,ISA)=FOP(LD1,ISA)+X6
          UN1(JJK,ISA)=UN1(JJK,ISA)-X5
          UP1(JJK,ISA)=UP1(JJK,ISA)-X6
          UK1(JJK,ISA)=UK1(JJK,ISA)-X10
          SOLK(LD1,ISA)=SOLK(LD1,ISA)+X10
      END IF
      IF(IDC(JJK)==NDC(3).OR.IDC(JJK)==NDC(6))THEN
          X7=.01*(HUI(JJK,ISA)+.01)**10*STL(JJK,ISA)
          STL(JJK,ISA)=STL(JJK,ISA)-X7
          DM(JJK,ISA)=DM(JJK,ISA)-X7F
          STD(JJK,ISA)=STD(JJK,ISA)+X7
          STDL(JJK,ISA)=STDL(JJK,ISA)+CLG(ISA)*X7
          X8=X7*BN(3,JJK)
          XUN=UN1(JJK,ISA)
          IF(XUN-X8<.01)X8=XUN-.01
          X9=X7*BP(3,JJK)
          XUP=UP1(JJK,ISA)
          IF(XUP-X9<.01)X9=XUP-.01
          STDN(JJK,ISA)=STDN(JJK,ISA)+X8
          STDP(JJK,ISA)=STDP(JJK,ISA)+X9
          UN1(JJK,ISA)=XUN-X8
          UP1(JJK,ISA)=XUP-X9
          IF(HUI(JJK,ISA)>.6.AND.STL(JJK,ISA)<.1)HU(JJK,ISA)=0.
      END IF
      SUM=0.
      DO J=1,NBSL(ISA)
          ISL=LID(J,ISA)
          RTO=WNMU(ISL,ISA)/(WNMU(ISL,ISA)+WNO3(ISL,ISA)+1.E-5)
          UU=UN(ISL)*RTO
          UN(ISL)=UN(ISL)-UU
          IF(WNO3(ISL,ISA)<UN(ISL))UN(ISL)=WNO3(ISL,ISA)
          IF(WNMU(ISL,ISA)<UU)UU=WNMU(ISL,ISA)
          WNO3(ISL,ISA)=WNO3(ISL,ISA)-UN(ISL)
          WNMU(ISL,ISA)=WNMU(ISL,ISA)-UU
          IF(DRW<0.)THEN      
              UTO=RWT(ISL,JJK,ISA)/RWL
              X1=-DRW*UTO
              X2=MIN(UN1(JJK,ISA),1000.*BN(3,JJK)*X1)
              CALL NCNSTD(X1,X2,ISL)
              UN1(JJK,ISA)=UN1(JJK,ISA)-X2
          ELSE
              UTO=UW(ISL)/(AEP(JJK)+1.E-20)
          END IF
          SWST(ISL,ISA)=MAX(1.E-10,SWST(ISL,ISA)-UW(ISL))
          X1=WPML(ISL,ISA)+WPMU(ISL,ISA)
          IF(X1>UP(ISL))THEN
              X2=WPML(ISL,ISA)/X1
              WPML(ISL,ISA)=WPML(ISL,ISA)-UP(ISL)*X2
              WPMU(ISL,ISA)=WPMU(ISL,ISA)-UP(ISL)*(1.-X2)
          ELSE
              UP(ISL)=X1
              WPML(ISL,ISA)=0.
              WPMU(ISL,ISA)=0.
          END IF
          X1=SOLK(ISL,ISA)+WKMU(ISL,ISA)
          IF(X1>UP(ISL))THEN
              X2=SOLK(ISL,ISA)/X1
              SOLK(ISL,ISA)=SOLK(ISL,ISA)-UK(ISL)*X2
              WKMU(ISL,ISA)=WKMU(ISL,ISA)-UK(ISL)*(1.-X2)
          ELSE
              UP(ISL)=X1
              SOLK(ISL,ISA)=0.
              WKMU(ISL,ISA)=0.
          END IF
          RWT(ISL,JJK,ISA)=RWT(ISL,JJK,ISA)+DRW*UTO
          SUM=SUM+RWT(ISL,JJK,ISA)
      END DO
      RW(JJK,ISA)=SUM
      !AD2=0.
      !AD3=0.
      !DO K=1,NBSL(ISA)
          !ISL=LID(K,ISA)
          !AD2=AD2+ST(ISL,ISA)
          !AD3=AD3+UW(ISL)
      !END DO
      !DF=AD1-AD2-AD3
      !IF(ABS(DF)>.001)WRITE(KW(1),1)IY,MO,KDA,AD1,AD2,AD3,DF
    !1 FORMAT(5X,'CAGRO',3I4,10E13.5)  
      RETURN
      END
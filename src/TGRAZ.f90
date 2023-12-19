      SUBROUTINE TGRAZ(JRT)
!     APEX1501
!     THIS SUBPROGRAM SIMULATES ANIMAL GRAZING.
      USE PARM 
	  JRT=0
	  LGZ=1
      II=IHC(JT1)
      LD1=LID(1,ISA)
      NN=NCP(IRO(ISA),ISA)
	  YY=0.
      YLN=0.
      YLP=0.
      IOW=IDON(ISA)
	  KOMP(KT(ISA),ISA)=0
	  IF(IHRD<2)THEN
	      X1=RST0(ISA)
	      IHD=IHDM(ISA)
	      GCOW(IHD,ISA)=WSA(ISA)/X1
	      IGZ(ISA)=1
	  ELSE
          DO IHD=1,NHRD(IOW)
              IF(IGZX(IHD,IOW)==ISA)EXIT
          END DO
          IF(IHD>NHRD(IOW))RETURN
          GCOW(IHD,ISA)=NCOW(IHD,IOW)
      END IF
      IF(IHRD==1)IGZ(ISA)=1          
      NGZ(IHD,ISA)=IHD
      IZ=IFED(IHD,IOW)
      IF(IZ==ISA)THEN
          IF(NCOW(IHD,IOW)>0)GCOW(IHD,ISA)=NCOW(IHD,IOW)*FFED(IHD,IOW)
      ELSE
          GCOW(IHD,ISA)=NCOW(IHD,IOW)*(1.-FFED(IHD,IOW))
	  END IF
      ! DEMAND  
	  DMD=GCOW(IHD,ISA)*GZRT(IHD,IOW)/(WSA(ISA)*HE(JT1))
      GZLX=GZLM(IHD,ISA)
	  IF(ORHI(JT1)>0.)GZLX=ORHI(JT1)
      GZLX=GZLX*GZRF(ISA)
      IF(AGPM(ISA)<GZLX.OR.GCOW(IHD,ISA)<1.E-5)THEN
          GCOW(IHD,ISA)=0.
          NGZ(IHD,ISA)=0
          IGZ(ISA)=0
          GZRF(ISA)=PRMT(48)
          IF(IHAY>0)THEN
              SPFD(ISA)=DMD*HE(JT1)
              DMD=0.
          ELSE
              RETURN
          END IF    
      ELSE
          GZRF(ISA)=1.
          SPFD(ISA)=0.
	  END IF
	  ! N CONC WEIGHTED FORAGE SUPPLY
      AD2=0.
      DO K=1,NN
          JJ=JE(K,ISA)
          IF(IDC(JJ)==NDC(7).OR.IDC(JJ)==NDC(8).OR.IDC(JJ)==NDC(10).OR.&
          IDC(JJ)==NDC(11))CYCLE
          IF(JJ==MNC)CYCLE
          IF(JP(JJ,ISA)<=0)THEN
              JP(JJ,ISA)=1
              NCR(JJ,ISA)=NCR(JJ,ISA)+1
          END IF
          CNLV(JJ)=MAX(BN(3,JJ),MIN(BN(1,JJ),UN1(JJ,ISA)/((DM(JJ,ISA)+1.E-5)*1000.)))
          GRLV(JJ)=STL(JJ,ISA)*CNLV(JJ)
          GRDD(JJ)=STD(JJ,ISA)*BN(3,JJ)
          AD2=AD2+GRLV(JJ)+GRDD(JJ)
      END DO
      IF(AD2>.01)THEN
          RTO=.001*DMD/AD2
      ELSE
          RTO=0.
      END IF    
      ! CORRECT N WT FORAGE SUPPLY WITH DEMAND
      DO K=1,NN
          IF(IGO(ISA)>0)THEN
              JJ=JE(K,ISA)
              IF(JJ==MNC)CYCLE
          ELSE
              JJ=1
          END IF    
          IF(IDC(JJ)==NDC(7).OR.IDC(JJ)==NDC(8).OR.IDC(JJ)==NDC(10).OR.&
          IDC(JJ)==NDC(11))CYCLE
          GRLV(JJ)=GRLV(JJ)*RTO
          GRDD(JJ)=GRDD(JJ)*RTO
      END DO 
      ! CALCULATE GRAZING YIELD FOR EACH FORAGE           
      DO K=1,NN
          IF(IGO(ISA)>0)THEN
              JJ=JE(K,ISA)
              IF(JJ==MNC)CYCLE
          ELSE
              JJ=1
          END IF    
          IF(IDC(JJ)==NDC(7).OR.IDC(JJ)==NDC(8).OR.IDC(JJ)==NDC(10).OR.&
          IDC(JJ)==NDC(11))CYCLE
          HUF(JJ,ISA)=MAX(HUF(JJ,ISA),HU(JJ,ISA))
          DMF(JJ,ISA)=DM1(JJ,ISA)
          TRA(JJ,ISA)=SRA(JJ,ISA)+TRA(JJ,ISA)
          IF(RD(JJ,ISA)>RDF(JJ,ISA))RDF(JJ,ISA)=RD(JJ,ISA)
          AJHI(JJ,ISA)=0.
          ! HARVEST INDEX FOR STANDING LIVE
          X1=MIN(GRLV(JJ)/(STL(JJ,ISA)+1.E-5),.9)
          ZZ=MAX(.01,1.-X1)
          ZZ2=ZZ**5
	      CPHT(JJ,ISA)=MAX(.001,CPHT(JJ,ISA)*ZZ2)
          HU(JJ,ISA)=MAX(.1,HU(JJ,ISA)*ZZ)
          SLAI(JJ,ISA)=MAX(.0001,SLAI(JJ,ISA)*ZZ2)
          STD(JJ,ISA)=MAX(.01,STD(JJ,ISA)*ZZ)
          STDO(ISA)=MAX(.01,STDO(ISA)*ZZ)
          ! GRAZING YIELD FOR STANDING LIVE 
          YLD(JJ)=GRLV(JJ)*HE(JT1)
          ! GRAZING YIELD FOR STANDING DEAD
          YLSD=GRDD(JJ)*HE(JT1)
          Y4=GRDD(JJ)*BN(3,JJ)
          Y5=GRDD(JJ)*BP(3,JJ)
          STDN(JJ,ISA)=STDN(JJ,ISA)-Y4
          STDP(JJ,ISA)=STDP(JJ,ISA)-Y5
          X4=MIN(GRLV(JJ)*CNLV(JJ),UN1(JJ,ISA))
          X3=UP1(JJ,ISA)/((DM(JJ,ISA)+.01)*1000.)
          X5=MIN(GRLV(JJ)*X3,UP1(JJ,ISA))
          Z2=YLSD*BN(3,JJ)
          Z3=YLSD*BP(3,JJ)
          YLN=MIN(.9*(UN1(JJ,ISA)+STDN(JJ,ISA)),YLD(JJ)*CNLV(JJ))
          YLP=MIN(.9*(UP1(JJ,ISA)+STDP(JJ,ISA)),YLD(JJ)*X3)
          X11=GRLV(JJ)-YLD(JJ)+GRDD(JJ)-YLSD
          X10=X4-YLN+Y4-Z2
	      JJK=JJ
          CALL NCNSTD(X11,X10,LD1)
          FOP(LD1,ISA)=MAX(.01,FOP(LD1,ISA)+X5-YLP+Y5-Z3)
          YY=YLD(JJ)+YLSD
          YLD2(JJ,ISA)=YLD2(JJ,ISA)+YY
          JD(ISA)=JJ
          SRA(JJ,ISA)=0.
          UN1(JJ,ISA)=UN1(JJ,ISA)-X4
          UP1(JJ,ISA)=UP1(JJ,ISA)-X5
          DM(JJ,ISA)=DM(JJ,ISA)-GRLV(JJ)
          IF(DM(JJ,ISA)<RW(JJ,ISA))RW(JJ,ISA)=DM(JJ,ISA)
          STL(JJ,ISA)=DM(JJ,ISA)-RW(JJ,ISA)
          YLN=YLN+Z2
          YLP=YLP+Z3
          YLNF(JJ,ISA)=YLNF(JJ,ISA)+YLN
          YLPF(JJ,ISA)=YLPF(JJ,ISA)+YLP
          TYN(ISA)=TYN(ISA)+YLN
          TYP(ISA)=TYP(ISA)+YLP
          Y1=1000.*YLD(JJ)
          Y2=1000.*YLSD
          IF(NOP>0.OR.NBSA(ISA)==ISAP)WRITE(KW(1),29)ISA,NBSA(ISA),IYR,&
          MO,KDA,IOW,IHD,TIL(JT1),CPNM(JD(ISA)),Y1,Y2,AGPM(ISA),STL(JJ,ISA),&
          STD(JJ,ISA),CNLV(JJ),BN(3,JJ),XHSM(ISA)
          IF(KFL(18)>0)THEN
              IF(K==1)THEN
                  X1=SPFD(ISA)
              ELSE
                  X1=0.
              END IF    
              WRITE(KW(18),3)ISA,NBSA(ISA),IYR,MO,KDA,IY,IOW,IHD&
              ,TIL(JT1),CPNM(JD(ISA)),Y1,Y2,X1,AGPM(ISA),STL(JJ,ISA),&
              STD(JJ,ISA),CNLV(JJ),BN(3,JJ)
          END IF    
      END DO 
      SMM(141,MO,ISA)=SMM(141,MO,ISA)+1.         
      RETURN
    3 FORMAT(1X,2I8,1X,6I4,1X,A8,1X,A4,6F10.2,6F10.4)
   29 FORMAT(1X,2I8,1X,3I4,2X,'IDON=',I4,2X,'HRD#=',I3,2X,A8,2X,A4,2X&
      ,'YLD=',F7.2,'kg/ha',2X,'YSD=',F7.2,'kg/ha',2X,'AGPM=',F7.2,'t/ha',&
      2X,'STL=',F7.2,'t/ha',2X,'STD=',F7.2,'t/ha',2X,'CNLV=',F7.4,'g/g',&
      2X,'CNDD=',F7.4,'g/g',2X,'HUSC=',F5.2)
      END
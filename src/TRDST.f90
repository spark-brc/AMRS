      SUBROUTINE TRDST
!     APEX1501
!     THIS SUBPROGRAM PLACES 10% OF CROP RESIDUE AND ITS N AND P CONTENT
!     IN THE SURFACE LAYER STORAGE WHEN THE CROP IS HARVESTED, CONVERTS
!     ROOT WEIGHT TO RESIDUE, AND ZEROS CROP GROWTH ACCUMULATORS.
      USE PARM 
      STD(JJK,ISA)=STD(JJK,ISA)+STL(JJK,ISA)
      X1=STL(JJK,ISA)+RW(JJK,ISA)
      XX=UN1(JJK,ISA)/X1
      X3=UP1(JJK,ISA)/X1
      STDN(JJK,ISA)=STDN(JJK,ISA)+XX*STL(JJK,ISA)
      STDP(JJK,ISA)=STDP(JJK,ISA)+X3*STL(JJK,ISA)
      STDL(JJK,ISA)=STDL(JJK,ISA)+CLG(ISA)*STL(JJK,ISA)
      DO KK=1,LRD(ISA)
          K=LID(KK,ISA)
          X1=RWT(K,JJK,ISA)
          X2=X1*XX
          CALL NCNSTD(X1,X2,K)
          FOP(K,ISA)=FOP(K,ISA)+X1*X3
          RWT(K,JJK,ISA)=0.
      END DO
      DM(JJK,ISA)=0.
      STL(JJK,ISA)=0.
      UN1(JJK,ISA)=0.
      UP1(JJK,ISA)=0.
      UK1(JJK,ISA)=0.
      RW(JJK,ISA)=0.
      RD(JJK,ISA)=0.
      CPHT(JJK,ISA)=0.
      AJH1(JJK,ISA)=0.
      X1=TRSD(ISA)+STD(JJK,ISA)+STDO(ISA)
      SRSD(ISA)=SRSD(ISA)+X1
      VIRT(ISA)=0.
      WS(ISA)=1.
      IGO(ISA)=MAX(0,IGO(ISA)-1)
      KGO(JJK,ISA)=0
      JE(JJK,ISA)=MNC
      HU(JJK,ISA)=0.
      HUI(JJK,ISA)=0.
      HSM(ISA)=0.
      SLAI(JJK,ISA)=0.
      WLV(JJK,ISA)=0.
      WCHT(JJK,ISA)=0.
      IYH(JJK,ISA)=0
      NII(ISA)=IRI(ISA)
      CSTF(JJK,ISA)=COST(ISA)
      COST(ISA)=0.
      IHU(JJK,ISA)=IHU(JJK,ISA)+1
      IF(IHU(JJK,ISA)>NHU(JJK,ISA))IHU(JJK,ISA)=1
      CAW(JJK,ISA)=AWC(JJK,ISA)
      ETG(JJK,ISA)=ACET(JJK,ISA)+ETG(JJK,ISA)
      PSTS(ISA)=0.
      IPST(ISA)=0
      NGD(JJK,ISA)=0
      FGC(ISA)=0.
      FGSL(ISA)=0.
      AJHI(JJK,ISA)=0.
      RETURN
      END
      SUBROUTINE NLGB(WTI,WTO,SOFU,WTB,WTE,ISA,KW,NBSA,MSO)
!     APEX1501
!     THIS SUBPROGRAM CALCULATES THE LAGOON MANURE BALANCE.
      DIMENSION NBSA(ISA),KW(MSO)
      DF=WTB+WTI-WTO-SOFU-WTE
	  PER=200.*DF/(WTB+WTE)
!	  IF(ABS(PER)<1.)RETURN
	  WRITE(KW(1),2)ISA,NBSA(ISA)
	  WRITE(KW(1),1)PER,DF,WTB,WTI,WTO,SOFU,WTE     
      RETURN
    1 FORMAT(5X,'PER =',E13.6,2X,'DF  =',E13.6,2X,'WTB =',E13.6,2X,&
      'WTI =',E13.6,2X,'WTO =',E13.6,2X,'OFMU=',E13.6/5X,'WTE =',E13.6)
    2 FORMAT(/T10,'LAGOON MANURE BALANCE',2X,'SA#= ',I8,1X,'ID= ',I8)
      END
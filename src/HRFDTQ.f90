      SUBROUTINE HRFDTQ
!     APEX1501
!     THIS SUBPROGRAM COMPUTES RAINFALL EXCESS USING THE GREEN & AMPT EQ.
      USE PARM
      X1=RFDT(1)
	  DO I=2,NRF
	      Q1=0.
	   	  X1=RFDT(I)-X1
		  RX=X1/DTHY
		  IF(RX>REP)REP=RX
		  IF(RFDT(I)>0.)CALL HGAWY(X1,RFDT(I),Q1,RX)
		  QGA(I)=Q1
		  X1=RFDT(I)
	  END DO
	  RETURN
      END
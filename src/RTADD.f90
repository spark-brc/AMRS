      SUBROUTINE RTADD
!     APEX1501
!     THIS SUBPROGRAM ADDS SUBAREA OUTPUTS TO ROUTED OUTPUTS TO
!     DETERMINE TOTAL OUTPUT FROM A ROUTING REACH.
      USE PARM
      IDO=IDOT(ICMD)
      IDN1=IDN1T(ICMD)
      IDN2=IDN2T(ICMD)
      Z1=RWSA(IDN1)
      Z2=RWSA(IDN2)
      !IDOR(IDN2)=IDO
      Z3=Z1+Z2
      RWSA(IDO)=Z3
      RQRB(IDO)=0.
      QVOL(IDO)=QVOL(IDN1)+QVOL(IDN2)
      QMM=.1*QVOL(IDO)/Z3
      VARH(2,IDO)=QVOL(IDO)/86400.
      SMMH(2,MO,IDO)=SMMH(2,MO,IDO)+VARH(2,IDO)
      RSSF(IDO)=RSSF(IDN1)+RSSF(IDN2)
      QRF(IDO)=QRF(IDN1)+QRF(IDN2)
      QDR(IDO)=QDR(IDN1)+QDR(IDN2)
      CPVH(IDO)=(CPVH(IDN1)*Z1+CPVH(IDN2)*Z2)/Z3
      QN(IDO)=QN(IDN1)+QN(IDN2)
      QC(IDO)=QC(IDN1)+QC(IDN2)
      QP(IDO)=QP(IDN1)+QP(IDN2)
      QPU(IDO)=(QPU(IDN1)*Z1+QPU(IDN2)*Z2)/Z3
      YSD(NDRV,IDO)=YSD(NDRV,IDN1)+YSD(NDRV,IDN2)
      YC(IDO)=YC(IDN1)+YC(IDN2)
      YN(IDO)=YN(IDN1)+YN(IDN2)
      YP(IDO)=YP(IDN1)+YP(IDN2)
      YPM(IDO)=(YPM(IDN1)*Z1+YPM(IDN2)*Z2)/Z3
      YMNU(IDO)=YMNU(IDN1)+YMNU(IDN2)
      YCOU(IDO)=YCOU(IDN1)+YCOU(IDN2)
      YNOU(IDO)=YNOU(IDN1)+YNOU(IDN2)
      YPOU(IDO)=YPOU(IDN1)+YPOU(IDN2)
      
      !rtb salt - add salt ion mass (kg) from the two converging subareas
      if(ISALT>0) then
      do m=1,mion
        surfqsalt(IDO,m) = surfqsalt(IDN1,m) + surfqsalt(IDN2,m)
        gwsalt(IDO,m) = gwsalt(IDN1,m) + gwsalt(IDN2,m)
        latsalt(IDO,m) = latsalt(IDN1,m) + latsalt(IDN2,m)
        surfsalt(IDO,m) = surfsalt(IDN1,m) + surfsalt(IDN2,m)
        qrfsalt(IDO,m) = qrfsalt(IDN1,m) + qrfsalt(IDN2,m)
      enddo
      endif
      
      QDRP(IDO)=(QDRP(IDN1)*Z1+QDRP(IDN2)*Z2)/Z3
      IF(ICDT(ICMD-2)==2)THEN
          II=NISA(IDNB(IDN1))
          IJ=IDOA(II)
          IDOR(II)=IDO
          SST(IDO)=SST(IJ)
      ELSE
          SST(IDO)=SST(IDN1)+SST(IDN2)
      END IF
      VARH(15,IDO)=SST(IDO)
      SMMH(15,MO,IDO)=SMMH(15,MO,IDO)+SST(IDO)
      TSFN(IDO)=TSFN(IDN1)+TSFN(IDN2)
      RSFN(IDO)=RSFN(IDN1)+RSFN(IDN2)
      QDRN(IDO)=QDRN(IDN1)+QDRN(IDN2)
      QRFN(IDO)=QRFN(IDN1)+QRFN(IDN2)
      SYC(IDO)=SYC(IDO)+YC(IDO)
      VYC(IDO)=YC(IDO)
      SQC(IDO)=SQC(IDO)+QC(IDO)
      VQC(IDO)=QC(IDO)
      SMQN(IDO)=QN(IDO)+RSFN(IDO)+QRFN(IDO)+QDRN(IDO) !total nitrogen from both subareas
      SQN(IDO)=SQN(IDO)+SMQN(IDO)
      
      !rtb salt
      !add salt ion mass (kg) from both subareas
      if(ISALT>0) then
      do m=1,mion
        SMQS(IDO,m) = surfqsalt(IDO,m) + gwsalt(IDO,m) + latsalt(IDO,m) + surfsalt(IDO,m) + qrfsalt(IDO,m)
        if(SMQS(IDO,m).lt.0) then
          SMQS(IDO,m) = 0. !this can happen if stream seepage (gwsw) removes all salt from the channel
        endif
        SMQS_total(IDO,m) = SMQS_total(IDO,m) + SMQS(IDO,m)
      enddo
      endif
      
      VAQN(IDO)=SMQN(IDO)
      SQP(IDO)=SQP(IDO)+QP(IDO)
      VAQP(IDO)=QP(IDO)
      SYN(IDO)=SYN(IDO)+YN(IDO)
      VAYN(IDO)=YN(IDO)
      SYP(IDO)=SYP(IDO)+YP(IDO)
      VAYP(IDO)=YP(IDO)
      WYLD(IDO)=QVOL(IDO)+RSSF(IDO)+QRF(IDO)+QDR(IDO)+CPVH(IDO)
      if(WYLD(IDO).lt.0) WYLD(IDO) = 0. !rtb modflow
      VARH(34,IDO)=WYLD(IDO)
      SMMH(34,MO,IDO)=SMMH(34,MO,IDO)+WYLD(IDO)
      VARH(35,IDO)=WYLD(IDO)/86400.
	  SMMH(35,MO,IDO)=SMMH(35,MO,IDO)+VARH(35,IDO)
      VSSN(IDO)=TSFN(IDO)
      SSN(IDO)=SSN(IDO)+TSFN(IDO)
      SP2=0.
      SP4=0.
      DO K=1,NDP
          QPST(K,IDO)=(QPST(K,IDN1)*Z1+QPST(K,IDN2)*Z2)/Z3
          TSPS(K,IDO)=(TSPS(K,IDN1)*Z1+TSPS(K,IDN2)*Z2)/Z3
          RSPS(K,IDO)=(RSPS(K,IDN1)*Z1+RSPS(K,IDN2)*Z2)/Z3
          YPST(K,IDO)=(YPST(K,IDN1)*Z1+YPST(K,IDN2)*Z2)/Z3
          VARP(2,K,IDO)=QPST(K,IDO)
          SMMP(2,K,MO,IDO)=SMMP(2,K,MO,IDO)+QPST(K,IDO)
          VARP(4,K,IDO)=TSPS(K,IDO)
          SMMP(4,K,MO,IDO)=SMMP(4,K,MO,IDO)+TSPS(K,IDO)
          VARP(5,K,IDO)=YPST(K,IDO)
          SMMP(5,K,MO,IDO)=SMMP(5,K,MO,IDO)+YPST(K,IDO)
          VARP(11,K,IDO)=RSPS(K,IDO)
          SMMP(11,K,MO,IDO)=SMMP(11,K,MO,IDO)+RSPS(K,IDO)
          SP2=SP2+QPST(K,IDO)
	      SP4=SP4+YPST(K,IDO)
      END DO
      VARH(6,IDO)=YSD(NDRV,IDO)
	  SMMH(6,MO,IDO)=SMMH(6,MO,IDO)+VARH(6,IDO)
	  VARH(9,IDO)=YN(IDO)+YNOU(IDO)
	  SMMH(9,MO,IDO)=SMMH(9,MO,IDO)+VARH(9,IDO)
	  VARH(11,IDO)=YP(IDO)+YPOU(IDO)
	  SMMH(11,MO,IDO)=SMMH(11,MO,IDO)+VARH(11,IDO)
	  VARH(13,IDO)=QN(IDO)    !+RSFN(IDO)+QRFN(IDO)+QDRN(IDO)
	  SMMH(13,MO,IDO)=SMMH(13,MO,IDO)+VARH(13,IDO)
	  VARH(19,IDO)=QP(IDO)
	  SMMH(19,MO,IDO)=SMMH(19,MO,IDO)+VARH(19,IDO)
	  VARH(27,IDO)=Z3*SP2
	  SMMH(27,MO,IDO)=SMMH(27,MO,IDO)+VARH(27,IDO)
      VARH(29,IDO)=Z3*SP4
	  SMMH(29,MO,IDO)=SMMH(29,MO,IDO)+VARH(29,IDO)
	  IF(WYLD(IDO)>0.)THEN
		  IF(QVOL(IDO)>0.)THEN
              TC(IDO)=(TC(IDN1)*QVOL(IDN1)+TC(IDN2)*QVOL(IDN2))/&
              (QVOL(IDN1)+QVOL(IDN2))
              IF(TC(IDO)>0.)THEN
	              TCAV(IDO)=TCAV(IDO)+TC(IDO)
			      ALTC=1.-EXP(-TC(IDO)*PRFF)
			      RQRB(IDO)=ALTC*QMM/TC(IDO)
			      PRAV(IDO)=PRAV(IDO)+RQRB(IDO)
			      IF(RQRB(IDO)>PRB(IDO))PRB(IDO)=RQRB(IDO)
			      NQRB(IDO)=NQRB(IDO)+1
			  END IF
	      ELSE
		      TC(IDO)=MAX(TC(IDN1),TC(IDN2))
			  PRFF=.042
			  ALTC=1.-EXP(-TC(IDO)*PRFF) 
	      END IF
	      !IF(KFL(9)>0.AND.ABS(Z3-RWSA(NCMD))<1.E-5)THEN
			  !X1=MAX(RQRB(IDO),X1/24.)
              !WRITE(KW(9),1202)IY,IYR,MO,KDA,QVOL(IDO),SST(IDO),&
              !QRF(IDO),RSSF(IDO),WYLD(IDO),X1,TC(IDO),ALTC
		  !END IF
	  END IF
      X1=YSD(NDRV,IDN1)*Z1
      X2=YSD(NDRV,IDN2)*Z2
      X3=X1+X2
      TOT=0.
      DO I=1,NSZ
          PCT(I,IDO)=MAX(.01,(PCT(I,IDN1)*X1+PCT(I,IDN2)*X2)/(X3+1.E-5))
          TOT=TOT+PCT(I,IDO)
      END DO
      PSZM(IDO)=0.
      DO I=1,NSZ
          PCT(I,IDO)=PCT(I,IDO)/(TOT+1.E-10)
          PSZM(IDO)=PSZM(IDO)+PSZY(I)*PCT(I,IDO)
      END DO
      IF(KFL(NOFL)>0.AND.ICMD==ICMO(IOF))THEN
	      IF(IY==1.AND.IDA==1)WRITE(KW(NOFL),4082)Z3
          XTP(1)=QVOL(IDO)+RSSF(IDO)+QRF(IDO)+SST(IDO)+QDR(IDO)
	      SMSO(IOF)=SMSO(IOF)+XTP(1)
          XTP(2)=Z3*YSD(NDRV,IDO)
          XTP(3)=YN(IDO)
          XTP(4)=YP(IDO)
          XTP(5)=QN(IDO)+RSFN(IDO)+QRFN(IDO)+QDRN(IDO)
          XTP(6)=QP(IDO)+QPU(IDO)
          WRITE(KW(NOFL),902)IDA,IYR,(XTP(I),I=1,6),(PSZ(I),PCT(I,IDO),&
          I=1,NSZ)
	      IOF=IOF+1
	      NOFL=NOFL+1
	  END IF
      IF(IHY==0)RETURN
      NHY(IDO)=MAX(NHY(IDN1),NHY(IDN2))
      IF(NHY(IDO)==0)RETURN
      IWH=0
      IF(QMM>QTH)IWH=1
      IF(KFL(26)>0.AND.IWH>0)THEN
          WRITE(KW(26),11)IY,MO,KDA,IDN1,IDN2,IDO
          WRITE(KW(26),9)
      END IF
      IF(IHY>0)THEN
          T1=0.
          QPK=0.
          SUM=0.
          DO K=1,NPD
              QHY(K,IDO,IHX(1))=QHY(K,IDN1,IHX(1))+QHY(K,IDN2,IHX(1))
              SUM=SUM+QHY(K,IDO,IHX(1))
              IF(QHY(K,IDO,IHX(1))>QPK)THEN
                  QPK=QHY(K,IDO,IHX(1))
                  TPK=T1
              END IF
              IF(KFL(26)>0)WRITE(KW(26),6)T1,QHY(K,IDN1,IHX(1)),QHY(K,IDN2,&
              IHX(1)),QHY(K,IDO,IHX(1))
              T1=T1+DTHY
          END DO
          X2=Z3/(DTHY*360.)
          RQRB(IDO)=QPK/X2
          SUM=(SUM-.5*(QHY(1,IDO,IHX(1))+QHY(NPD,IDO,IHX(1))))/X2
          HYDV(IDO)=SUM
          TC(IDO)=TPK
          !NQRB(IDO)=NQRB(IDO)+1
          !PRAV(IDO)=PRAV(IDO)+RQRB(IDO)
          !PRSD(IDO)=PRSD(IDO)+RQRB(IDO)*RQRB(IDO)
          !IF(RQRB(IDO)>PRB(IDO))PRB(IDO)=RQRB(IDO)
          !X2=RQRB(IDO)/(SUM+1.E-10)
          !IF(X2>QRQB(IDO))QRQB(IDO)=X2
          !QRBQ(IDO)=QRBQ(IDO)+X2
          TCAV(IDO)=TCAV(IDO)+TC(IDO)
          IF(TC(IDO)<TCMX(IDO))THEN
              IF(TC(IDO)<TCMN(IDO))TCMN(IDO)=TC(IDO)
          ELSE
              TCMX(IDO)=TC(IDO)
          END IF
          IF(KFL(26)>0.AND.IWH>0)WRITE(KW(26),13)Z1,Z2,Z3,HYDV(IDN1),&
          HYDV(IDN2),SUM,QPK,TPK
          I=24
          DO J=2,MHX
              T1=T1-DTHY
              DO K=1,NPD
                  I=I+1
                  QHY(K,IDO,IHX(J))=QHY(K,IDN1,IHX(J))+QHY(K,IDN2,IHX(J))
                  IF(KFL(26)>0)WRITE(KW(26),6)T1,QHY(K,IDN1,IHX(J)),QHY(K,IDN2,&
                  IHX(J)),QHY(K,IDO,IHX(J))
                  T1=T1+DTHY
              END DO
              IF(I>NHY(IDO).AND.QHY(NPD,IDO,IHX(J))<.1)EXIT
	      END DO
	  END IF
      RETURN
    6 FORMAT(9X,F8.3,4F10.3)
    9 FORMAT(//13X,'Th',5X,'QHY1m3/s',2X,'QHY2m3/s',2X,'QHYOm3/s')
   11 FORMAT(//T10,'ADD HYD'/T10,'DATE(Y-M-D)=',3I4,2X,'IDN1= ',I8,2X&
      ,'IDN2= ',I8,2X,'IDO= ',I8)
   13 FORMAT(T10,'WSA(IDN1)= ',F12.3,' ha',2X,'WSA(IDN2)= ',F12.3,' ha',&
      2X,'WSA(IDO)= ',F12.3,' ha'/T10,'QID1= ',F8.3,' mm',2X,'QID2= ',&
      F8.3,' mm',2X, 'QO= ',F8.3,' mm'/T10,'PEAK RATE= ',F8.3,' m3/s',2X,&
      'TP= ',F7.2,' h')
   14 FORMAT(T10,'ADD SED= ',F7.3,' t/ha')  
   12 FORMAT(5X,A2,3I8,3I4,5F10.2)
   17 FORMAT(20X,10E13.5)
  902 FORMAT(1X,I4,1X,I4,1X,20(1X,E16.6))
 1202 FORMAT(1X,4I4,15F10.2)  
 4082 FORMAT(///10X,'WATERSHED AREA = ',F10.2,' HA'/1X,'JDA   YR     TOT &
      FLOWm3        YSDt             YNkg             YPkg          TOT &
      NO3kg          QPkg',11X,3('PSZum',11X,'FRACTION',8X))  
      END
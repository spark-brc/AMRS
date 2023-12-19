      SUBROUTINE ARESET 
!     APEX1501
!     THIS SUBPROGRAM RESETS VARIABLES TO INITIAL VALUES FOR MULTI RUNS.
      USE PARM
      IHRL=0
      NWDA=0
      LW=1
!     IYR=IYR0
      NWPD=0
      SHYD=0.
      SQVL=0.
      SRMX=0.
      THRL=0.
      YLK=0.
      YLN=0.
      YLP=0.
      DO J=1,NRN0
          IX(J)=IX0(J)
      END DO
      XIM=0.
      V3=AUNIF(IDG(3))
      V1=AUNIF(IDG(2))
      JCN=0
      JCN0=0
      JCN1=0
      CST1=0.
      VALF1=0.
      PRSD=0.
      QRQB=0.
      QRBQ=0.
      TCPA=0.
      TCPY=0.
      TLMF=0.
      TYK=0.
      TYN=0.
      TYP=0.
      CYAV=0.
      CYSD=0.
      CYMX=0.
      FPF=0.
      SSFI=0.
      TR=0.
      TSN=0.
      TSY=0.
      TVGF=0.
      TYW=0.
      TQ=0.
      TCN=0.
      TYON=0.
      TYTP=0.
      TQP=0.
      TQN=0.
      TXMX=0.
      TXMN=0.
      TAMX=0.
      TSR=0.
      TCVF=0.
      CX=1.E-10
      TEI=0.
      SET=0.
      TET=0.
      ASW=0.
      SRD=0.
      QIN=0.
      TRHT=0.
      NYLN=0
      STV=0.
      IF(NDP>0)THEN
          VQ=0.
          VY=0.
          SQB=0.
          SYB=0.
          !SPQ=0.
          !SPY=0.
          APQ=0.
          APY=0.
          AQB=0.
          AYB=0.
          PVY=0.
          PVQ=0.
          SMAP=0.
          SMYP=0.
          SMMP=0.
      END IF
      NCR=0
      ACET=0.
      TYL1=0.
      TYL2=0.
      TYLK=0.
      TYLN=0.
      TYLP=0.
      TETG=0.
      TCAW=0.
      TDM=0.
      TRA=0.
      THU=0.
      TRD=0.
      TSFC=0.
      STDA=0.
      DO ISA=1,MSA
          DO I=1,NBSL(ISA)
              J=LID(I,ISA)
              SOL(1,J,ISA)=WPMA(J,ISA)
              SOL(2,J,ISA)=WON(J,ISA)
              SOL(4,J,ISA)=WPO(J,ISA)
              SOL(5,J,ISA)=WPMS(J,ISA)
              SOL(7,J,ISA)=WOC(J,ISA)
          END DO
      END DO
      IHT=0
      CNDS=0.
      NQRB=0
      PRAV=0.
      PRB=0.
      DRAV=0.
      ERAV=0.
      DPMT=0.
      TCMX=0.
      TCMN=100.
      TCAV=0.
      IF(NDP/=0)THEN
          SMYRP=0.
          SMMRP=0.
          SMRP=0.
      END IF
      TSFN=0.
      TSFK=0.
      QRFN=0.
      SRCH=0.
      SM=0.
      SMY=0.
      SMM=0.
      SMYR=0.
      SMMR=0.
      SMR=0.
      TSPS=0.
      VARP=0.
      RETURN
      END
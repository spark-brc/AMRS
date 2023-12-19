      SUBROUTINE ADEALLOCATE                                                      
!     THIS SUBROUTINE DEDEALLOCATES ARRAY SIZES                                          
      USE PARM                                                                       
      DEALLOCATE (CPNM,FPSO,FTNM,PSTN,TITOP,TITSO,HEDH,HED)                                                            
	  DEALLOCATE (KW,ICDT,IDN1T,IDN2T,IDNB,IDOT,IDRO,NHY,NQRB,NTX,NISA,&
	  ICUS,IHC,MASA,NBE,NBT,IDC,KDC,NTP,KDF,KFL)                                                     
	  DEALLOCATE (ISAL,ISAS,NFED,NHRD,NSAL,NSAO,NSAS)                                                               
	  DEALLOCATE (IIR,KIR,NIR)                                               
      DEALLOCATE (IAC,IAMF,IAPL,IAUF,IAUI,IAUL,IBSA,IDFH,IDNF,IDOA,IDON,&
      IDOR,IDR,IDRL,IDS,IEXT,IFA,IFD,IFLS,IGO,IGZ,IHDM,ILQF,IMW,IPMP,&
      IPSO,IPST,IPTS,IRF,IRI,IRO,IRP,IRR,IRRS,IRRW,ISAO,ISCP)
      DEALLOCATE (ISG,ISPF,IWTH,JBG,JCN,JCN0,JCN1,JD,JSA,KC,KP1,KT,KTF,&
      KTMX,KTT,LM,LRD,LUN,LUNS,MXSR,NBCF,NBCT,NBFF,NBFT,NBSA,NBSL,NBW,&
      NDFA,NII,NMW,NPSF,NRO,NVCN,NWDA)
      DEALLOCATE (IPSF,IDG,IX,IX0,KPSN,JPC,KPC,NPC,IHX) 
	  DEALLOCATE (IHU,IYH,JE,JP,JPL,KGO,NCR,NGD,NHU,NYLN)             
	  DEALLOCATE (NCP,NFRT,NPST,NTL)                    
      DEALLOCATE (IDFA,IDFD,IDMU,IGZO,IGZX,IHBS,IYHO,LGIR,NCOW,NGZA,&
      NHBS,NYHO)                  
      DEALLOCATE (IDSL,IDSS,IDOW,IHT,KOMP,IFED,NGZ,LID,LORG,IFLO,IHRL,&
      IDFT,IDF0)                          
	  DEALLOCATE (ITL,JH,KDT,LFT,LT,LYR,LPC,LY,IHDT,NGIX,NBSX)                                 
      DEALLOCATE (OSAA,OWSA,PKRZ,RNMN,UK,UN,UP,UW)
      DEALLOCATE (FCST,FK,FN,FNMA,FNMN,FNO,FOC,FP,FPO,FSLT)                                                 
	  DEALLOCATE (PCST,PHLF,PHLS,PKOC,PLCH,PSOL,PWOF,SSPS)                                                 
      DEALLOCATE (COOP,COTL,DKH,DKI,EFM,EMX,FPOP,FRCP,FULU,HE,HMO,ORHI,&
      RHT,RIN,RR,STIR,TIL,TLD)                                      
      DEALLOCATE (AEP,ALT,CAF,CKY,CNLV,CNY,CSTS,CPY,DDM,DLAI,DMLA,&
      DMLX,EP,EXTC,FLT,FTO,GMHU,GRDD,GRLV,GSI,HI,HMX,PHUX,PLAX,POPX,&
      PRYF,PRYG,PST,RBMD,RDMX,RLAD,SDW,TBSC,TCPA,TCPY,TOPC,VPD2,VPTH,&
      WA,WAI,WAVP,WAX,WCY,WSYF,WXYF,XDLAI,XMTU,YLD,YLX)                    
      DEALLOCATE (ABD,AFLG,AGPM,ALGI,ALQ,ARMN,ARMX,ARSD,BA1,BA2,BCOF,&
      BCV,BFFL,BFSN,BFT,BGWS,BIG,BIR,BR1,BR2,BRSV,BSALA,BSNO,BSW,BTC,&
      BTCX,BTCZ,BTK,BTN,BTNX,BTNZ,BTP,BTPX,BTPZ,BV1,BV2,BVIR,CFNP,CHL,&
      CHMX,CHN,CHS,CHXA,CHXP,CLG,CN0,CN2,CNSX)
      DEALLOCATE (COST,COWW,CPMX,CST1,CV,CVF,CVP,CVRS,CYAV,CYMX,CYSD,&
      DALG,DDLG,DEPC,DHT,DKIN,DKHL,DRT,DST0,DWOC,EFI,EK,EM10,EVRS,EVRT,&
      FBAR,FBM,FCMN,FCMP,FDSF,FFC,FFPQ,FGC,FGSL,FHP,FIRG,FPF,FPSC,FSFN,&
      FSFP,GMA,GRDL,GWE0,GWEL,GWSN,GWST,GWMX,GZRF,HCLD,HCLN,HLMN,HR0,&
      HSM,OCPD)
      DEALLOCATE (OMAP,ORSD,PAW,PCOF,PDAW,PDPL,PDPL0,PDPLC,PDPLX,PDSKC,&
      PDSW,PEC,PMX,PM10,PRSD,PSTF,PSTM,PSTS,QCAP,QRBQ,QRQB,RCBW,RCF,RCHC,&
      RCHD,RCHK,RCHL,RCHN,RCHS,RCHX,RCSS,RCTW,REPI,RFPK,RFPL,RFPS,RFPW,&
      RFPX,RFQN,RFTT,RFV,RFV0,RHD,RHTT)
      DEALLOCATE (RINT,RLF,RMXS,ROSP,RRUF,RSAE,RSAP,RSBD,RSDP,RSEE,RSEP,&
      RSF,RSHC,RSK,RSLK,RSOC,RSON,RSOP,RSO3,RSRR,RSSA,RSSP,RST0,RSV,RSVB,&
      RSVE,RSVF,RSVP,RSPK,RSYB,RSYF,RSYN,RSYS,RVE0,RVP0,RZ,RZSW,&
      SALA,SALB,SAMA,SATK,SCI,SCNX,SDVR,SLF,SLT0,SLTX,SMAS,SMEO,SMFN,&
      SMFU,SMKS,SMLA,SMMU,SMNS,SMNU,SMPL,SMPQ,SMPS,SMPY,SMRF,SMSS,SMST,&
      SMTS,SMWS,SMX,SMY1,SMY2,SNO,SOLQ,SPFD,SPLG,SRAD,SRSD,SSFI,SSIN,SSW)                     
      DEALLOCATE (ST0,STDO,STDOK,STDON,STDOP,STKR,STLT,STP,SW,&
      SWLT,SYN,S3,TAGP,TCC,TCS,TFLG,THK,TILG,TKR,TLMF,TMN,TMX,TNOR,&
      TOC,TPSF,TRSD,TSLA,TSMY,TSNO,TVGF,TYK,TYN,TYP,U10,UB1,UOB,UPSX,&
      URBF,USL,VAC,VALF1,VAP,VCHA,VCHB,VFPA,VFPB,VIMX,VIRT,VLG,VLGB,&
      VLGI,VLGM,VLGN,VPU,VRSE,VSK,VSLT,WDRM)                                                     
      DEALLOCATE (WK,WS,WSA,WSX,WTMB,WTBL,WTMN,WTMU,WTMX,XCT,XHSM,XIDK,&
      XIDS,XMAP,XNS,XRFI,YCT,YLC,YLS,YTN,YTX,ZBMC,ZBMN,ZCO,ZCOB,ZEK,ZFK,&
      ZFOP,ZHPC,ZHPN,ZHSC,ZHSN,ZLM,ZLMC,ZLMN,ZLS,ZLSC,ZLSL,ZLSLC,ZLSLNC,&
      ZLSN,ZNH3,ZNO3,ZNMU,ZNOA,ZNOS,ZNOU,ZOC,ZON,ZPMA,ZPML,ZPMS,ZPMU,&
      ZPO,ZPOU,ZSK,ZSLT,ZTP)                                               
      DEALLOCATE (CPVH,DPMT,DRAV,ERAV,HYDV,RCTC,PRAV,PRB,PSZM,QC,QDR)
      DEALLOCATE (QDRN,QDRP,QN,QP,QPR,QPU,QRF,QRFN,QRFP,QRP,QURB,QVOL,&
      RQRB,RSFN)
      DEALLOCATE (RSSF,RWSA,SHYD,SMIO,SMQN,SMYN,SQC,SQN,SQP,SQVL,SSN,SST,&
      STY,SYC)
      DEALLOCATE (SYP,TC,TCAV,TCMN,TCMX,TSFK,TSFN,VAQN,VAQP,VAYN,VAYP,&
      VQC,VSSN,VYC,WYLD)
      DEALLOCATE (YC,YCOU,YCWN,YMNU,YN,YNOU,YNWN,YP,YPM,YPOU,YPWN,YW)
      DEALLOCATE (PSO3,PSON,PSOP,PSOQ,PSOY,PQPS,PSSP,PYPS)
      DEALLOCATE (QGA,RFDT,VARW)
      DEALLOCATE (EO5,RF5,SCFS,XMS,ASW,CX,QIN,SET,SRD,SRMX,TAMX,TCN,TCVF,&
      TEI,TET,THRL,TQ,TR,TRHT,TSN,TSR,TSY,TXMX,TXMN,TQN,TQP,TQPU,TYON,&
      TYTP,TYW,CNSC)
      DEALLOCATE (ALS,BD,BDD,BDM,BDP,BPT,CAC,CBN,CDG,CEC,CLA,CNDS,&
      CNRT,CPRH,CPRV,DHN,ECND,EQKE,EQKS,EXCK,FE26,FIXK,FOP,HCL,PH,PO,&
      PSP,ROK,RSD,RSDM,SAN,SATC,SEV,SIL,SMB,SOLK,SULF,SUT)
      DEALLOCATE (SWST,STFR,STMP,VNO3,WBMN,WCMU,WCOU,WHPC,WHPN,WHSC,&
      WHSN,WKMU,WLM,WLMC,WLMN,WLS,WLSC,WLSL,WLSLC,WLSLNC,WLSN,&
      WNMU,WNOU,WOC,WON,WPMA,WPML,WPMS,WPMU,WPO,WPOU,WSLT,WT,Z)
      DEALLOCATE (ACET,AJH1,AJHI,AWC,CAW,CNT,CPHT,CSTF,DM,DMF,DM1,ETG,&
      FRTK,FRTN,FRTP,HU,HUF,HUI,PPL0,RD,RDF,REG,RW,SLAI,SLA0,SRA,STD,&
      STDK,STDL,STDN)               
      DEALLOCATE (STDP,STL,SWH,SWP,TCAW,TDM,TETG,TFTK,TFTN,TFTP,THU,TRA,&
      TRD,TVIR,TYL1,TYL2,TYLK,TYLN,TYLP,UK1,UNA,UN1,UP1,VIR,WCHT,WLV,&
      XLAI,XDLA0,YLD1,YLD2,YLKF,YLNF,YLPF)
      DEALLOCATE(ACO2C,AFP,AN2OC,AO2C,CGCO2,CGN2O,CGO2,CLCO2,CLN2O,CLO2)
      DEALLOCATE(DCO2GEN,DN2G,DN2OG,DO2CONS,DPRC,DPRN,DPRO,DRWX,EAR,FC,&
      HKPC,HKPN)
      DEALLOCATE(HKPO,RSPC,RWTZ,S15,SADST,SMEA,SMES,SOT,SSFCO2,SSFN2O,&
      SSFO2)
      DEALLOCATE(TPOR,VCO2,VGA,VGN,VN2O,VO2,VFC,VWC,VWP,WBMC,WCO2G,&
      WCO2L,WN2O,WN2OG,WN2OL,WNO2,WNH3,WNO3,WO2G,WO2L,XN2O,ZC)    
      DEALLOCATE (YHY,SMH,SMYH,VARH,QPST,RSPS,TSPS,YPST,SRCH,CPFH,QSF,&
      SSF,SM,SMUA,SMY,SMYUA,VAR,VARUA,CTSA,VQ,VY,GWPS,PFOL,ANA,FNMX,&
      GCOW,GZLM,DUMP,FFED,GZRT,VURN,PPX,YSD,PCT,PCTH,SQB,SYB,FNP,SMYRP,&
      BK,BN,BP,BLG,BWN,DLAP,FRST,PPCF,PPLP,RWPC,STX,WAC2)    
      DEALLOCATE (SMMH,PVQ,PVY,PSTE,PSTR,PSSF,PSTZ,CND,SMM,SMMUA,STV,&
      SMS,SOL,FIRX,QIR,RSTK,TIR,VIRR,WFA,WMUCH,PHU,POP,PPLA,VARC,SMAP,&
      SMYP,VARP,SMMRP,SMRP)
      !SPQ(5,20,MHY),SPQC(5,20,MHY),SPY(5,20,MHY)
      DEALLOCATE (TSFC,SOIL,HUSC,RWT,STDA,SFCP,SFMO,XZP,QHY)
      DEALLOCATE (SMMP,SMMC)
      ! DEALLOCATE (APQ(5,20,100,MHY),APQC(5,20,100,MHY),APY(5,20,100,MHY)&              
      ! ,AQB(5,20,100,MHY),AYB(5,20,100,MHY),SMMP(20,MPS,12,MSA),SMMC(15,&             
      ! MNC,12,MSA))
      ! deallocate paddy arrays Jaehak 2016                                                                                        
      DEALLOCATE(SUB_D50,PADDY_STO,PADDY_OUT,PADDY_HWEIR,PADDY_HMIN)
      DEALLOCATE(HWEIR,LWEIR,PADDY_LWEIR,HMIN_IRR)
      deallocate(xgo2,xgco2,xgn2o)    ! Variables added by Luca Doro
      RETURN                                                                         
      END                                                                            
                                                                                     
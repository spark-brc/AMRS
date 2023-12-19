      MODULE PARM
      !IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER(80)::FCROP,FFERT,FHERD,FMLRN,FOPSC,FPARM,FPEST,&
      FPRNT,FSITE,FSOIL,FSUBA,FTILL,FTR55,FWIND,FWLST,FWPM,SAFILE,&
      SITEFILE
      CHARACTER(80)::FWTH(1000)
	  CHARACTER(4)::TITLE(60),SID(11),HDRP(5),HEDP(13),HEDR(20),HEDC(17)&
      ,HEDS(20)
	  CHARACTER(2)::CMDX(5)
	  CHARACTER(80),DIMENSION(:),ALLOCATABLE::FPSO
	  CHARACTER(20),DIMENSION(:),ALLOCATABLE::TITOP,TITSO
      CHARACTER(16),DIMENSION(:),ALLOCATABLE::PSTN
      CHARACTER(8),DIMENSION(:),ALLOCATABLE::FTNM,TIL
      CHARACTER(4),DIMENSION(:),ALLOCATABLE::CPNM,HED,HEDH
      INTEGER::IBD,IBDT,ICDP,ICMD,ICO2,ICP,IDA,IDAY,IDIR,IDN1,IDN2,IDNT,&
      IDO,IERT,IET,IGC,IGN,IGSD,IHAY,IHD,IHRD,IHV,IHUS,IHY,III,IKAT,IMON,&
      INFL,INP,IOF,IOW,IOX,IPAT,IPC,IPCD,IPD,IPF,IPL,IPRK,IPY,IPYI,IREM,&
      IRFT,IRGX,IRH,IRUN,ISA,ISAP,ISCN,ISLF,ISL,ISTA,ISW,IT1,IT2,IT3,&
      ITYP,IUN,IWI,IWTB,IY,IYER,IYR,IYR0,IYX,JD0,JDA,JDE,JDHU,JJK,JT1,JT2
      INTEGER::KDA,KF,KI,KP,LBP,KND,LC,LGRZ,LGZ,LND,LNS,LPD,LPYR,LW,MASP,&
      MBS,MCA12,MFT,MHD,MHP,MHX,MHY,MIR,ML1,MNC,MNT,MNUL,MO,MO1,MOW,MPO,&
      MPS,MRO,MSA,MSC,MSCP,MSL,MXT,MXW,NAQ,NBCL,NBCX,NBDT,NBFX,NBMX,&
      NBON,NBYR,NCMD,NCMO,ND,NDF,NDP,NDRV,NDT,NDVSS,NDWT,NEV,NGN,NGN0,&
      NJC,NKA,NKD,NKS,NKY,NOFL,NOP,NPD,NPRC,NPSO,NRF,NRN0,NRNG,NSH,NSM,&
      MSO,NSNN,NSTP,NSX,NSZ,NYD,NT0,NT1,NTV,NUPC,NWP,NWTH
      INTEGER::KDT1(800),KDC1(1500),KDF1(1200),KDP1(1000),KA(100),&
      NXP(90),KR(60),KDT2(50),KD(40),KY(40),KDP(50),NHC(27),KS(20),&
      NWDR(16),NC(13),IAD(10),ICMO(10),NWPD(10),NDC(11),JX(7),KGN(5),&
      JC(4)     
	  INTEGER,DIMENSION(:),ALLOCATABLE::IAC,IAMF,IAPL,IAUF,IAUI,IAUL,&
      IBSA,ICUS,IDFH,IDG,IDNB,IDNF,IDOA,IDON,IDOR,IDR,IDRL,IDRO,IDS,IEXT,&
      IFA,IFD,IFLS,IGO,IGZ,IHDM,IHX,ILQF,IMW,IPMP,IPSF,IPSO,IPST,IPTS,&
      IRF,IRI,IRO,IRP,IRR,IRRS,IRRW,ISAL,ISAO,ISAS,ISCP,ISPF,IX,IX0,JBG,&
      JCN,JCN0,JCN1,JD,JSA,KC,KFL,KP1,KPSN,KT,KTF,KTMX,KTT,KW,LM,LRD,LUN,&
      LUNS,NBSL,NBW,NDFA,NFED,NHRD,NII,NISA,NPSF,NRO,NSAL,NSAO,NSAS,NTP,&
      NVCN,NBE,NBT,IHC,ICDT,IDC,IDN1T,IDN2T,IDOT,IIR,ISG,IWTH,JPC,KDC,&
      KDF,KIR,KPC,MASA,MXSR,NBCF,NBCT,NBFF,NBFT,NBSA,NHY,NIR,NMW,NPC,&
      NQRB,NTX,NWDA,KRST	
      INTEGER,DIMENSION(:,:),ALLOCATABLE::IDF0,IDFA,IDFD,IDFT,IDMU,IDOW,&
      IDSL,IDSS,IFED,IFLO,IGZO,IGZX,IHBS,IHRL,IHT,IHU,IYH,IYHO,JE,JP,JPL,&
      KGO,KOMP,LGIR,LID,LORG,NCOW,NCP,NCR,NFRT,NGD,NGZ,NGZA,NHBS,NHU,&
      NPST,NTL,NYHO,NYLN
      INTEGER,DIMENSION(:,:,:),ALLOCATABLE::IHDT,ITL,JH,KDT,LFT,LPC,LT,&
      LY,LYR,NBSX,NGIX
      REAL::ACW,ADEO,ADHY,ADRF,AEPT,AET,AFN,AGP,AHSM,AJWA,ALG,ALMN,ALTC,&
      AL5,AMPX,ANG,AMSR,AVT,BETA,BRNG,BS1,BS2,BXCT,&
      BYCT,CBVT,CLF,CLT,CMM,CMS,CN,CNO3I,COIR,COL,CON,COP,CO2,COS1,CPH0,&
      CPV0,CPVV,CQNI,CRLNC,CRUNC,CSFX,CSLT,D150,DAYL,DD,DEMR,DN2O,DRSW
      REAL::DRTO,DTG,DTHY,DUR,DURG,DVOL,DZ10,DZDN,EI,ELEV,EO,ERTN,ERTO,&
      ERTP,ES,EXNN,EXPK,FL,FRSD,FULP,FW,GWBX,GX,HGX,HMN,HRLT,HR1,PB,PI2,&
      PIT,PMOEO,PMORF,PRFF,PSTX,QAPY,QPQ,QRB,QSFN,QSFP,QTH,RAMX,REP,RFNC,&
      RFRA,RGIN,RGRF,RHCF,RHM,RHP,RHS,RM,RMNR,RMX0
      REAL::RNO3,RRF,RTP,RUNT,RWO,SAET,SALF,SAT,SBAL,SCN,SCN2,SDEG,&
      SDEP,SDN,SDRF,SEP,SGMN,SHRL,SK,SML,SMNR,SMP,SN,SN2,SN2O,SNA15,&
      SNA30,SNIT,SNMN,SNN15,SNN30,SNOF,SNPKT,SOFU,SP,SPRC,SPRF,SRPM,&
      SSFK,SSFN,SSST,STND,SUK,SUN,SUP,SVOL
      REAL::SW15,SW30,SX,SYMU,SYW,TA,TDEG,TDEP,TEVP,TEV1,THW,TLA,&
      TMAF,TMAP,TMPD,TPRK,TPSP,TREV,TRFR,TRSP,TSAE,TSFS,TX,TXXM,UK2,UKM,&
      UNM,UNR,UN2,UPM,UPR,UP2,USTRT,USTT,USTW,UX,UXP,VGF,VMU,VPD,&
      V1,V3,V56,V57,WAGE,WB,WBMX,WCF,WCYD,WDN,WFX,WIM,WIP,WKA,WKMNH3,&
      WKMNO2,WKMNO3,WMP,WNCMAX,WNCMIN,WTN,XCS,XDA,XDA1,XET,XK1,XK2,XKN1,&
      XKN3,XKN5,XKP1,XKP2,XSA,XSL,YCS,YERO,YEW,YEWN,YLAT,YLAZ,YLN,YLP,&
      YSL,YWKS,ZF,ZMIX,ZQT      
      REAL::XTP(100),XYP(100),PRMT(110),SCLM(38),SCMN(30),SRNG(30),&
      SMYR(21),VARS(20),RNCF(12),TAV(12),TMNF(12),TMXF(12),UAV0(12),&
      UAVM(12),AWXP(10),RFSG(10),SMSO(9),OPV(7),SMSW(6),PSZ(5),PSZX(5),&
      PSZY(5),WX(3),XAV(3),XDV(3),XIM(3),XRG(3)
	  REAL::SMMR(21,12),SMR(21,12),DIR(12,16),OBMN(10,12),&
	  OBMNX(10,12),OBMX(10,12),OBSL(10,12),PCF(10,12),RH(10,12),SDTMN(10,12),&
	  SDTMX(10,12),WFT(10,12),WI(10,12),SCRP(35,2)
	  REAL::CQRB(8,17,4),RST(3,10,12),PRW(2,10,12)
	  REAL,DIMENSION(:),ALLOCATABLE::ABD,AEP,AFLG,AGPM,ALGI,ALQ,ALT,ARMN,&
      ARMX,ARSD,BA1,BA2,BCOF,BCV,BFFL,BFSN,BFT,BGWS,BIG,BIR,BR1,BR2,&
      BRSV,BSALA,BSNO,BSW,BTC,BTCX,BTK,BTCZ,BTN,BTNX,BTNZ,BTP,BTPX,BTPZ,&
      BV1,BV2,BVIR,CAF,CAP,CFNP,CHL,CHMX,CHN,CHS,CHXA,CHXP,CKY,CLG,CN0,&
      CN2,CNLV,CNSX,CNY,COOP,CSTS,COST,COTL,COWW,CPMX,CPVH,CPY,CST1,CV,&
      CVF,CVP,CVRS,CYAV,CYMX,CYSD,DALG,DDLG,DDM,DEPC,DHT,DKH,DKI,DKHL,&
      DKIN,DLAI,DMLA,DMLX,DPMT,DRAV,DRT,DST0,DWOC,EFI,EFM,EK,EM10,EMX,&
      EP,ERAV,EVRS,EVRT,EXTC,FBAR,FBM,FCMN,FCMP,FCST,FDSF,FDSR,FFC,FFPQ,&
      FGC,FGSL,FHP,FIRG,FK,FLT,FN,FNMA,FNMN,FNO,FOC,FP,FPO,FPF,FPOP,&
      FPSC,FRCP,FSLT,FSFN,FSFP,FTO
      REAL,DIMENSION(:),ALLOCATABLE::FULU,GMA,GMHU,GRDD,GRDL,GRLV,GSI,&
      GWE0,GWEL,GWMX,GWSN,GWST,GZRF,HCLD,HCLN,HE,HI,HLMN,HMO,HMX,HR0,&
      HSM,HYDV,OCPD,OMAP,ORHI,ORSD,OSAA,OWSA,PAW,PCOF,PCST,PDAW,PDPL,&
      PDPL0,PDPLC,PDPLX,PDSKC,PDSW,PEC,PHLF,PHLS,PHUX,PKOC,PKRZ,PLCH,&
      PM10,PMX,PLAX,POPX,PQPS,PRAV,PRB,PRSD,PRYF,PRYG,PSO3,PSOL,PSON,&
      PSOP,PSOQ,PSOY,PSSP,PST,PSTF,PSTM,PSTS,PSZM,PWOF,PYPS,QC,QCAP,QDR,&
      QDRN,QDRP,QGA,QN,QP,QPR,QPU,QRBQ,QRF,QRFN,QRFP,QRP,QRQB,QURB,QVOL,&
      RBMD,RCBW,RCF,RCHC,RCHD,RCHK,RCHL,RCHN,RCHS,RCHX,RCSS,RCTC,RCTW,&
      RDMX,REPI,RFDT,RFPK,RFPL,RFPS,RFPW,RFPX,RFQN,RFTT,RFV,RFV0,RHD,RHT,&
      RHTT,RIN,RINT,RLAD,RLF,RMXS,ROSP,RQRB,RR,RRUF,RSAE,RSAP,RSBD,RSDP,&
      RSEE,RSEP,RSF,RSFN,RSHC,RSK,RSLK,RSO3
      REAL,DIMENSION(:),ALLOCATABLE::RSOC,RSON,RSOP,RSRR,RSSA,RSSF,&
      RSSP,RST0,RSV,RSVB,RSVE,RSVF,RSVP,RSPK,RSYB,RSYF,RSYN,RSYS,RVE0,&
      RVP0,RWSA,RZ,RZSW,S3,SALA,SALB,SAMA,SATK,SCI,SCNX,SDVR,SDW,SHYD,&
      SLF,SLT0,SLTX,SMAS,SMEO,SMFN,SMFU,SMIO,SMKS,SMLA,SMMU,SMNS,SMNU,&
      SMPL,SMPQ,SMPS,SMPY,SMQN,SMRF,SMSS,SMST,SMTS,SMWS,SMX,SMY1,SMY2,&
      SMYN,SNO,SOLQ,SPFD,SPLG,SQC,SQN,SQP,SQVL,SRAD,SRSD,SSFI,SSIN,SSN,&
      SSPS,SST,SSW,ST0,STDO,STDOK,STDON,STDOP,STIR,STKR,STLT,STP,STY,SW,&
      SWLT,SYC,SYN,SYP,TAGP,TBSC,TC,TCAV,TCC,TCMN,TCMX,TCPA,TCPY,TCS,&
      TFLG,THK,TILG,TKR,TLD,TLMF,TMN,TMX,TNOR,TOC,TOPC,TPSF,TRSD,TSFK,&
      TSFN,TSLA,TSMQ,TSMY,TSNO,TVGF,TYK,TYN,TYP,U10,UB1,UK,UN,UOB,UP,&
      UPSX,URBF,USL
      REAL,DIMENSION(:),ALLOCATABLE::UW,VAC,VALF1,VAP,VAQN,VAQP,VARW,&
      VAYN,VAYP,VCHA,VCHB,VFPA,VFPB,VIMX,VIRT,VLG,VLGB,VLGI,VLGM,VLGN,&
      VPD2,VPTH,VPU,VQC,VYC,VRSE,VSK,VSLT,VSSN,WA,WAI,WAVP,WAX,WCY,WDRM,&
      WK,WS,WSA,WSX,WSYF,WTBL,WTMB,WTMN,WTMU,WTMX,WXYF,WYLD,XCT,XDLAI,&
      XHSM,XIDK,XIDS,XMAP,XMTU,XNS,XRFI,YC,YCOU,YCT,YCWN,YLC,YLD,YLS,&
      YLX,YMNU,YN,YNOU,YNWN,YP,YPM,YPOU,YPWN,YTN,YTX,YW,ZBMC,ZBMN,ZCO,&
      ZCOB,ZEK,ZFK,ZFOP,ZHPC,ZHPN,ZHSC,ZHSN,ZLM,ZLMC,ZLMN,ZLS,ZLSC,ZLSL,&
      ZLSLC,ZLSLNC,ZLSN,ZNH3,ZNO3,ZNMU,ZNOA,ZNOS,ZNOU,ZOC,ZON,ZPMA,ZPML,&
      ZPMS,ZPMU,ZPO,ZPOU,ZSK,ZSLT,ZTP        
      REAL,DIMENSION(:,:), ALLOCATABLE::ACO2C,AFP,AN2OC,AO2C,CDG,CGCO2,&
      CGN2O,CGO2,CLCO2,CLN2O,CLO2,DCO2GEN,DN2G,DN2OG,DO2CONS,DPRC,DPRN,&
      DPRO,EAR,FC,HKPC,HKPN,HKPO,S15,SADST,SMEA,SMES,SOT,TPOR,VFC,&
      VO2,VWC,VWP,WBMC,WCO2G,WCO2L,WN2O,WN2OG,WN2OL,WNO2,WNO3,WO2G,WO2L,&
      XN2O,ZC,xgo2,xgco2,xgn2o    ! xgo2,xgco2,xgn2o added by Luca Doro - 20190827
      REAL,DIMENSION(:,:),ALLOCATABLE::ACET,AJH1,AJHI,ALS,ANA,ASW,AWC,&
      CX,BD,BDD,BDM,BDP,BK,BLG,BN,BP,BPT,BWN,CAC,CAW,CBN,CEC,CLA,&
      CNDS,CNRT,CNSC,CNT,CPFH,CPHT,CPRH,CPRV,CTSA,CSTF,DHN,DLAP,DM,DMF,&
      DM1,DRWX,DUMP,ECND,EO5,EQKE,EQKS,ETG,EXCK,FE26,FFED,FIXK,FNMX,FNP,&
      FOP,FRST,FRTK,FRTN,FRTP,GCOW,GWPS,GZLM,GZRT,HCL,HU,HUF,HUI,PCT,&
      PCTH,PFOL,PH,PO,PPCF,PPL0,PPLP,PPX,PSP,QIN,QPST,QSF,RD,RDF,REG,&
      RF5,RNMN,ROK,RSD,RSDM,RSPC,RSPS,RTN,RW,RWPC,RWTZ,SAN,SATC,SCFS,&
      SET,SEV,SIL,SLAI,SLA0,SM,SMB,SMUA,SMH
      REAL,DIMENSION(:,:),ALLOCATABLE::SMYH,SMY,SMYUA,SMYRP,SOLK,SQB,&
      SRA,SRCH,SRD,SRMX,SSF,SSFCO2,SSFN2O,SSFO2,SWST,STD,STDK,STDL,STDN,&
      STDP,STFR,STL,STMP,STX,SULF,SUT,SWH,SWP,SYB,TAMX,TCAW,TCN,TCVF,&
      TDM,TEI,TET,TETG,TFTK,TFTN,TFTP,THRL,THU,TQ,TR,TRA,TRD,TRHT,TSN,&
      TSR,TSY,TXMX,TXMN,TQN,TQP,TQPU,TSPS,TVIR,TYL1,TYL2,TYLK,TYLN,TYLP,&
      TYON,TYTP,TYW,UK1,UNA,UN1,UP1,VAR,VARUA,VARH,VCO2,VGA,VGN,VIR,&
      VN2O,VNO3,VQ,VURN,VY,WAC2,WBMN,WCHT,WCMU,WCOU,WHPC,WHPN,WHSC,WHSN,&
      WKMU,WNH3,WLM,WLMC,WLMN,WLS,WLSC,WLSL,WLSLC,WLSLNC,WLSN,WLV,WNMU,&
      WNOU,WOC,WON,WPMA,WPML,WPMS,WPMU,WPO,WPOU,WSLT,WT,XDLA0,XLAI,XMS,&
      XSP,YHY,YLD1,YLD2,YLKF,YLNF,YLPF,YPST,YSD,Z
      REAL,DIMENSION(:,:,:),ALLOCATABLE::CND,FIRX,HUSC,PHU,POP,PPLA,PSTE,&
      PSTR,PSSF,PSTZ,PVQ,PVY,QHY,QIR,RSTK,RWT,SFCP,SFMO,SMAP,SMM,SMMUA,&
      SMMH,SMMRP,SMRP,SMS,SMYP,SOIL,SOL,STDA,STV,TIR,TSFC,VARC,VARP,&
      VIRR,WFA,WMUCH,XZP
      !SPQ,SPQC,SPY,
	  REAL,DIMENSION(:,:,:,:),ALLOCATABLE::APQ,APQC,APY,AQB,AYB,SMMP,SMMC
      !PADDY MODEL
      REAL,DIMENSION(:),ALLOCATABLE::SUB_D50,PADDY_HWEIR,PADDY_LWEIR,PADDY_HMIN
      REAL,DIMENSION(:,:),ALLOCATABLE::PADDY_STO,PADDY_OUT
      REAL,DIMENSION(:,:,:),ALLOCATABLE::HWEIR,LWEIR,HMIN_IRR
      REAL::PADDY_SED_NORM,PADDY_IRRMM,LAI_INIT
      !MODFLOW
      INTEGER:: IMF=0
      INTEGER:: idaf
      !SALINITY !rtb salt -----------------------------------------------------------------------------------------------------------------
      INTEGER:: ISALT=0
      integer:: salt_point=0
      integer :: mion=8 !number of salt ions in salinity module
      character (len=80) :: salt_species(8)
      !soil
      real, dimension (:,:,:), allocatable :: wsalt !salt ion mass in soil layers (kg/ha)
      real, dimension (:,:,:), allocatable :: soil_csalt !salt ion concentration in soil water (g/m3)
      real, dimension (:,:), allocatable :: vsalt !salt ion mass in percolating water (kg/ha)
      
      !surface runoff
      real, dimension (:,:), allocatable :: surfqsalt !salt ion mass in surface runoff (kg/ha)
      real, dimension (:,:), allocatable :: surfqsalt_sub

      !water erosion
      real, dimension (:,:), allocatable :: surfsalt !salt ion mass in soil erosion from water (kg/ha)
      real, dimension (:,:), allocatable :: surfsalt_sub
      
      !irrigation water
      real, dimension (:,:), allocatable :: irrsalt !salt ion mass in irrigation water (kg/ha)
      real, dimension (:,:), allocatable :: irrig_csalt !salt ion concentration in irrigation water (g/m3)

      !lateral flow
      real, dimension (:,:), allocatable :: latsalt !salt ion mass in lateral return flow (kg/ha)
      real, dimension (:,:), allocatable :: latsalt_sub
      
      !quick return flow
      real, dimension (:,:), allocatable :: qrfsalt !salt ion mass in quick return flow (kg/ha)
      real, dimension (:,:), allocatable :: qrfsalt_sub
      
      !groundwater return flow
      real, dimension (:,:), allocatable :: gwsalt !salt ion mass in groundwater return flow (kg/ha)
      real, dimension (:,:), allocatable :: gwswsalt_sub,swgwsalt_sub
      
      !point source
      real, dimension (:,:), allocatable :: pts_sub !salt ion mass in point source
      
      !river water routing
      real, dimension (:,:), allocatable :: SMQS,SMQS_total !total salt ion mass in the subarea channel

      !APEX-RT3D linkage
      real, dimension (:,:),   allocatable :: sub_salt_exchange !salt ion loadings via gw-sw exchange (kg)
      real, dimension (:,:),   allocatable :: rchrg_salt !salt ion leaching to water table (kg/ha)
      real, dimension (:,:),   allocatable :: rchrg_csalt !salt ion concentration in leaching water to water table (g/m3)
      real, dimension (:,:),   allocatable :: chan_csalt !salt ion concentration in subarea channel (g/m3)
      
      !equilibrium chemistry
      real, dimension (:,:,:,:),   allocatable :: salt_solid_aqu !salt mineral fractions for each grid cell
      real, dimension (:,:,:),   allocatable :: salt_solid_soil !salt mineral fractions for each subarea soil layer
      real, dimension (:,:,:),   allocatable :: init_salt_fraction !initial salt mineral fractions for each subarea soil layer
      real, dimension (:,:,:),   allocatable :: init_salt_conc
      double precision, dimension (:), allocatable :: Sul_Conc,Cal_Conc,&
		                                                  Mg_Conc,Sod_Conc,&
                                                      Pot_Conc,Cl_Conc,&
                                                     Car_Conc,BiCar_Conc
      integer c11,c22,salt_c3,salt_c4,c5
      double precision salt_K1,salt_K2,salt_K3,salt_K4,salt_K5
      double precision Ksp11,Ksp21,Ksp31,Ksp41,Ksp51
      double precision Ksp12,Ksp22,Ksp32,Ksp42,Ksp52
      double precision upion1,upion2,upion3,upion4,&
                       upion5,upion6,upion7,upion8
      double precision Sol_CaCO3(1000),Sol_MgCO3(1000),Sol_CaSO4(1000),&
                       Sol_MgSO4(1000),Sol_NaCl(1000)
      double precision, dimension (:), allocatable :: LAMDA
      
      real, dimension (:,:),   allocatable :: saltdissolve_soil !change in salt mass due to salt mineral precipitation-dissolution (kg)
      real, dimension (:),   allocatable :: salt_soil !total salt in soil (kg)
      
      !point sources
      real, dimension (:,:,:,:),   allocatable :: salt_ptloads
      
      !time-variable soil hydrologic group
      logical :: hsg_vary
      integer :: num_soil_type
      real, dimension (:,:,:),   allocatable :: hsg_values
      real, dimension (:,:,:),   allocatable :: cnsc1_vary,cnsc2_vary
      
      END 

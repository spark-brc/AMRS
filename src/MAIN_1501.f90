!     PROGRAM APEX1501
!     THIS MODEL IS CALLED APEX(AGRICULTURAL POLICY/ENVIRONMENTAL
!     EXTENDER)
!     IT IS A COMPREHENSIVE AGRICULTURAL MANAGEMENT MODEL THAT IS USEFUL
!     IN SOLVING A VARIETY OF PROBLEMS INVOLVING SUSTAINABLE
!     AGRICULTURE, WATER QUALITY, WATER SUPPLY, AND GLOBAL CLIMATE
!     & CO2 CHANGE.  THE MODEL SIMULATES MOST COMMONLY USED
!     MANAGEMENT PRACTICES LIKE TILLAGE, IRRIGATION, FERTILIZATION,
!     PESTICIDE APPLICATION, LIMING, AND FURROW DIKING.
!     THE WATERSHED MAY BE DIVIDED INTO SUBAREAS(SUBWATERSHEDS,
!     FIELDS, LANDSCAPE POSITIONS, ETC). RUNOFF, SEDIMENT, AND AG
!     CHEMICALS ARE ROUTED AS SURFACE, SUBSURFACE, AND CHANNEL FLOW.
!     THE MAIN PROGRAM INITIALIZES VARIABLES, ACCUMULATES
!     MONTHLY AND ANNUAL VALUES, AND COMPUTES AVERAGE VALUES FOR THE
!     SIMULATION PERIOD.
	  PROGRAM MAIN
	   USE PARM
      use amrt_parm !APEX-MODFLOW
	  INTEGER::  YY0,MM0,DD0 
	   CHARACTER(1)::ASG
      CHARACTER(2)::RFPT
      CHARACTER(4)::HEDW
      CHARACTER(7)::AUNT
	   CHARACTER(8)::ASOL,RCMD
      CHARACTER(10)::FMT
      CHARACTER(80)::ASOT,ASTN,FPSOD,OPSCFILE,FRFDT,RFTFILE,SOILFILE,&
      TITWP,WINDFILE,WPMFILE,filename
      CHARACTER(80)::ADFL,ADUM,AWTH,TITWN,TMPSTR,STR
      REAL::TMPFL
      integer :: j,yrc,dayc,year_index,year_days !rtb hsg
      DIMENSION HEDW(16),TITWP(10),ASG(4)
	   DIMENSION IDOS(1000),NX(100),NY(20),KFL0(51),IWPM(10)
      DIMENSION SAV(1000),SMWD(3000),&!XGCO2(200),XGN2O(30),XGO2(30),&    ! Commented, added in allocated_parms, deallocate, modparms - Luca Doro 20190827
      TWB(11),YTP(16),FWXP(10),COCH(6),RCMD(5),RFPT(4)
      DIMENSION RMO(10,12),SCRX(35,2)
	   DATA RFPT/' 1','1A',' 2',' 3'/,RCMD/'SUBAREA ','ROUTE   ',&
      'ADD     ','ROUTE RS','ROUTE PD'/,ASG/'A','B','C','D'/,HEDW/'   N'&
      ,' NNE','  NE',' ENE','   E',' ESE','  SE',' SSE','   S',' SSW',&
      '  SW',' WSW','   W',' WNW','  NW',' NNW'/
     
      ! Version description
      ! STR='APEX1501-MODFLOW for the White River Watershed. rev.05112020 '
      WRITE(*,*) ''
      STR='  APEX-MODFLOW-RT3D-Salt (AMRS) rev.24-003  '
      WRITE(*,*) STR
      WRITE(*,*) ''
     
      CALL AHEAD
      ADUM='APEXRUN.DAT'
      CALL OPENV(KR(11),ADUM,0)
      ADUM='APEXFILE.DAT'
      CALL OPENV(KR(12),ADUM,0)
      ADUM='APEXCONT.DAT'
      CALL OPENV(KR(20),ADUM,0)
      ADUM='APEXDIM.DAT'
      CALL OPENV(KR(26),ADUM,0)
      ! READ DIMENSIONS
      ! 1 MPS = MAX # PESTICIDES
      ! 2 MRO = MAX # YRS CROP ROTATION
      ! 3 MNT = MAX # TILLAGE OPERATIONS
      ! 4 MNC = MAX # CROPS USED
      ! 5 MHD = MAX # ANIMAL HERDS
      ! 6 MBS = MAX # BUY/SELL LIVESTOCK TRANSACTIONS
      ! 7 MFT = MAX # FERTILIZER
      ! 8 MPO = MAX # POINT SOURCES
      ! 9 MHP = MAX# HYDROGRAPH POINTS
      !10 MHX = MAX# DAYS FOR STORM HYDROGRAPH BASE
      !11 MSA = MAX# SUBAREAS
      !12 MIR = MAX# IRRIGATION APPLICATIONS
      READ(KR(26),*)MPS,MRO,MNT,MNC,MHD,MBS,MFT,MPO,MHP,MHX,MSA,MIR
      !READ FILE NAMES
      READ(KR(12),509)FSITE,FSUBA,FWPM,FWIND,FCROP,FTILL,FPEST,FFERT,&
      FSOIL,FOPSC,FTR55,FPARM,FMLRN,FPRNT,FHERD,FWLST,FPSOD,FRFDT
      !  WRITE(*,508)'FSITE',FSITE,'FSUBA',FSUBA,'FWPM',FWPM,'FWIND',FWIND,&
      !'FCROP',FCROP,'FTILL',FTILL,'FPEST',FPEST,'FFERT',FFERT,'FSOIL',&
      !FSOIL,'FOPSC',FOPSC,'FTR55',FTR55,'FPARM',FPARM,'FMLRN',FMLRN,&
      !'FPRNT',FPRNT,'FHERD',FHERD,'FWLST',FWLST,'FPSOD',FPSOD,'FRFDT',FRFDT

      !  1  NBY0 = NUMBER OF YEARS OF SIMULATION
      !  2  IYR0 = BEGINNING YEAR OF SIMULATION
      !  3  IMO  = MONTH SIMULATION BEGINS
      !  4  IDA  = DAY OF MONTH SIMULATION BEGINS
      !  5  IPD  = N0 FOR ANNUAL WATERSHED OUTPUT
      !          = N1 FOR ANNUAL OUTPUT                       | N YEAR INTERVAL
      !          = N2 FOR ANNUAL WITH SOIL TABLE              | N=0 SAME AS
      !          = N3 FOR MONTHLY                             | N=1 EXCEPT
      !          = N4 FOR MONTHLY WITH SOIL TABLE             | N=0 PRINTS
      !          = N5 FOR MONTHLY WITH SOIL TABLE AT HARVEST  | OPERATIONS
      !          = N6 FOR N DAY INTERVAL
      !          = N7 FOR SOIL TABLE ONLY N DAY INTERVAL
      !          = N8 FOR SOIL TABLE ONLY DURING GROWING SEASON N DAY INTERVAL
      !          = N9 FOR N DAY INTERVAL DURING GROWING SEASON
      !  6  NGN0 = ID NUMBER OF WEATHER VARIABLES INPUT.  RAIN=1,  TEMP=2,
      !            RAD=3,  WIND SPEED=4,  REL HUM=5.  IF ANY VARIABLES ARE INP
      !            RAIN MUST BE INCLUDED.  THUS, IT IS NOT NECESSARY TO SPECIF
      !            ID=1 UNLESS RAIN IS THE ONLY INPUT VARIABLE.
      !            NGN=-1 ALL VARIABLES GENERATED(SAME VALUES ALL SUBAREAS).
      !            NGN=0  ALL VARIABLES GENERATED(SPATIALLY DISTRIBUTED).  
      !            EXAMPLES
      !            NGN=1 INPUTS RAIN.
      !            NGN=23 INPUTS RAIN, TEMP, AND RAD.
      !            NGN=2345 INPUTS ALL 5 VARIABLES.
      !  7  IGN  = NUMBER TIMES RANDOM NUMBER GEN CYCLES BEFORE
      !            SIMULATION STARTS
      !  8  IGS0 = 0 FOR NORMAL OPERATION OF WEATHER MODEL.
      !          = N NO YRS INPUT WEATHER BEFORE REWINDING (USED FOR REAL TIME
      !            SIMULATION).
      !  9  LPYR = 0 IF LEAP YEAR IS CONSIDERED
      !          = 1 IF LEAP YEAR IS IGNORED
      ! 10  IET  = PET METHOD CODE
      !          = 1 FOR PENMAN-MONTEITH
      !          = 2 FOR PENMAN
      !          = 3 FOR PRIESTLEY-TAYLOR
      !          = 4 FOR HARGREAVES
      !          = 5 FOR BAIER-ROBERTSON
      ! 11  ISCN = 0 FOR STOCHASTIC CURVE NUMBER ESTIMATOR.
      !          > 0 FOR RIGID CURVE NUMBER ESTIMATOR.
      ! 12  ITYP = 0 FOR MODIFIED RATIONAL EQ STOCHASTIC PEAK RATE ESTIMATE.
      !          < 0 FOR MODIFIED RATIONAL EQ RIGID PEAK RATE ESTIMATE.
      !          > 0 FOR SCS TR55 PEAK RATE ESTIMATE.
      !          = 1 FOR TYPE 1 RAINFALL PATTERN
      !          = 2     TYPE 1A
      !          = 3     TYPE 2
      !          = 4     TYPE 3
      !          = 5 FOR SCS UNIT HYD PEAK RATE ESTIMATE
      ! 13  ISTA = 0 FOR NORMAL EROSION OF SOIL PROFILE
      !          = 1 FOR STATIC SOIL PROFILE
      ! 14  IHUS = 0 FOR NORMAL OPERATION
      !          = 1 FOR AUTOMATIC HEAT UNIT SCHEDULE(PHU MUST BE INPUT AT
      !              PLANTING)
      ! 15  NVCN0= 0 VARIABLE DAILY CN NONLINEAR CN/SW WITH DEPTH SOIL WATER
      !              WEIGHTING
      !          = 1 VARIABLE DAILY CN NONLINEAR CN/SW NO DEPTH WEIGHTING
      !          = 2 VARIABLE DAILY CN LINEAR CN/SW NO DEPTH WEIGHTING
      !          = 3 NON-VARYING CN--CN2 USED FOR ALL STORMS
      !          = 4 VARIABLE DAILY CN SMI(SOIL MOISTURE INDEX)
      ! 16  INFL = 0 FOR CN ESTIMATE OF Q
      !          = 1 FOR GREEN & AMPT ESTIMATE OF Q, RF EXP DST, PEAK RF RATE
      !              SIMULATED.
      !          = 2 FOR G&A Q, RF EXP DST, PEAK RF INPUT
      !          = 3 FOR G&A Q, RF UNIFORMLY DST, PEAK RF INPUT
      !          = 4 FOR G&A Q, RF INPUT AT TIME INTERVAL DTHY
      ! 17  MASP = 1 FOR PESTICIDE APPLIED IN g/ha
      !          = 1000 FOR PESTICIDE APPLIED IN kg/ha
      ! 18  IERT = 0 FOR EPIC ENRICHMENT RATIO METHOD
      !          = 1 FOR GLEAMS ENRICHMENT RATIO METHOD
      ! 19  LBP  = 0 FOR SOL P RUNOFF ESTIMATE USING GLEAMS PESTICIDE EQ
      !          = 1 FOR LANGMUIR EQ
      ! 20  NUPC = N AND P PLANT UPTAKE CONCENTRATION CODE
      !          = 0 FOR SMITH CURVE
      !          > 0 FOR S CURVE
      ! 21  MNUL = MANURE APPLICATION CODE
      !          = 0 FOR AUTO APPLICATION TO SUBAREA WITH MINIMUM LAB P CONC
      !          = 1 FOR VARIABLE P RATE LIMITS ON ANNUAL APPLICATION BASED ON
      !            JAN 1 LAB P CONC.
      !          = 2 FOR VARIABLE N RATE LIMITS ON ANNUAL APPLICATION BASED ON
      !            JAN 1 LAB P CONC.
      !          = 3 SAME AS 1 EXCEPT APPLICATIONS OCCUR ON ONE SUBAREA AT A 
      !            TIME UNTIL LAB P CONC REACHES 200 ppm. THEN ANOTHER SUBAREA
      !            IS USED, ETC.
      ! 22  LPD  = DAY OF YEAR TO TRIGGER LAGOON PUMPING DISREGARDING NORMAL
      !            PUMPING TRIGGER--USUALLY BEFORE WINTER OR HIGH RAINFALL 
      !            SEASON.
      !          = 0 DOES NOT TRIGGER EXTRA PUMPING
      ! 23  MSCP = INTERVAL FOR SCRAPING SOLID MANURE FROM FEEDING AREA (d)
      ! 24  ISLF = 0 FOR RUSLE SLOPE LENGTH/STEEPNESS FACTOR
      !          > 0 FOR MUSLE SLOPE LENGTH/STEEPNESS FACTOR
      ! 25  NAQ  > 0 FOR AIR QUALITY ANALYSIS
      !          = 0 NO AIR QUALITY ANALYSIS
      ! 26  IHY  = 0 NO FLOOD ROUTING
      !          = 1 VSC FLOOD ROUTING
      !          = 2 SVS FLOOD ROUTING
      !          = 3 MUSKINGUM-CUNGE VC
      !          = 4 MUSKINGUM-CUNGE M_CVC4 
      ! 27  ICO2 = 0 FOR CONSTANT ATMOSPHERIC CO2
      !          = 1 FOR DYNAMIC ATMOSPHERIC CO2
      !          = 2 FOR INPUTTING CO2
      ! 28  ISW  = 0 FIELD CAP/WILTING PT EST RAWLS METHOD DYNAMIC.
      !          = 1 FIELD CAP/WILTING PT INP RAWLS METHOD DYNAMIC.
      !          = 2 FIELD CAP/WILTING PT EST RAWLS METHOD STATIC.
      !          = 3 FIELD CAP/WILTING PT INP STATIC.
      !          = 4 FIELD CAP/WILTING PT NEAREST NEIGHBOR DYNAMIC
      !          = 5 FIELD CAP/WILTING PT NEAREST NEIGHBOR STATIC
      !          = 6 FIELD CAP/WILTING PT BEHRMAN-NORFLEET-WILLIAMS (BNW) DYNAMIC                         
      !          = 7 FIELD CAP/WILTING PT BEHRMAN-NORFLEET-WILLIAMS (BNW) STATIC
      ! 29  IGMX = # TIMES GENERATOR SEEDS ARE INITIALIZED FOR A SITE.
      ! 30  IDIR = 0 FOR READING DATA FROM WORKING DIRECTORY
      !          > 0 FOR READING FROM WEATDATA DIRECTORY
      ! 31  IMW0 = MIN INTERVAL BETWEEN AUTO MOW 
      ! 32  IOX  = 0 FOR ORIGINAL EPIC OXYGEN/DEPTH FUNCTION
      !          > 0 FOR ARMEN KEMANIAN CARBON/CLAY FUNCTION
      ! 33  IDNT = 1 FOR ORIGINAL EPIC DENITRIFICATION SUBPROGRAM
      !          = 2 FOR ARMEN KEMANIAN DENITRIFICATION SUBPROGRAM
      !          = 3 FOR CESAR IZAURRALDE DENITRIFICATION SUBPROGRAM (ORIGINAL DW)
      !          = 4 FOR CESAR IZAURRALDE DENITRIFICATION SUBPROGRAM (NEW DW)     
      ! 34  IAZM = 0 FOR USING INPUT LATITUDES FOR SUBAREAS
      !          > 0 FOR COMPUTING EQUIVALENT LATITUDE BASED ON AZIMUTH 
      !            ORIENTATION OF LAND SLOPE.
      ! 35  IPAT = 0 TURNS OFF AUTO P APPLICATION
      !          > 0 FOR AUTO P APPLICATION
      ! 36  IHRD = 0 FOR LEVEL 0(MANUAL) GRAZING MODE(NO HERD FILE REQUIRED)
      !          = 1 FOR LEVEL 1(HYBRID) GRAZING MODE(HERD FILE REQUIRED) 
      !          = 2 FOR LEVEL 2(AUTOMATIC) GRAZING MODE(HERD FILE REQUIRED)  
      ! 37  IWTB = DURATION OF ANTEDEDENT PERIOD FOR RAINFALL AND PET 
      !            ACCUMULATION TO DRIVE WATER TABLE.
      ! 38  IKAT = 0 TURNS OFF AUTO K APPLICATION
      !          > 0 FOR AUTO K APPLICATION
      ! 39  NSTP = REAL TIME DAY OF YEAR
      ! 40  IPRK = 0 FOR HPERC 
      !          > 0 FOR HPERC1 (4MM SLUG FLOW)
      ! 41  ICP  = 0 FOR NCNMI_PHOENIX
      !          > 0 FOR NCNMI_CENTURY
      ! 42  NTV  = 0 FOR ORIGINAL APEX NITVOL EQS
      !          > 0 FOR IZAURRALDE REVISED NITVOL EQS
      ! 43  IREM = 0 SSK FROM REMX
      !          > 0 SSK FROM USLE
      ! 44  IHAY = 0 FOR NO HAY FEEDING
      !          > 0 FOR HAY FEEDING WHEN FORAGE IS < GZLM
      ! 45  ISAP = DOUBLE VARIABLE (ISA/IPCD) TO PRINT DAILY BALANCES USING SR NCONT
      !            IPCD=1 FOR C
      !                =2 FOR H2O
      !                =3 FOR P
      !                =4 FOR N
      !            NBSA=SELECTED SA ID #   
      ! 46   IMF = 0 NO MODFLOW SIMULATION
      !           >0 TURN ON MODFLOW SIMULATION
      ! 47   ISALT = 0 NO SALT ION SIMULATION !rtb salt
      !             > 0 TURN ON SALT ION SIMULATION
      !     LINE 1/2
      READ(KR(20),*)NBY0,IYR0,IMO,IDA0,IPD,NGN0,IGN,IGSD,LPYR,IET,&
      ISCN,ITYP,ISTA,IHUS,NVCN0,INFL,MASP,IERT,LBP,NUPC,MNUL,LPD,MSCP,&
      ISLF,NAQ,IHY,ICO2,ISW,IGMX,IDIR,IMW0,IOX,IDNT,IAZM,IPAT,IHRD,IWTB,&
      IKAT,NSTP,IPRK,ICP,NTV,IREM,IHAY,ISAP,IMF,ISALT !rtb salt
      
      IF(IMF>0)THEN
         !rtb MODFLOW: prepare APEX-MODFLOW linkage
         call amrt_read_link           
         print *
         print *, 'Preparing APEX-MODFLOW linkage...'
         call MapRiver_Grid
         call MapSubArea_Grid
         call MapGrid_SubArea  
         print *         
      ENDIF
      
      IF(IWTB==0)IWTB=15
      CALL AISPL(ISAP,JJ)
      IPCD=ISAP
      ISAP=JJ
      IF(ISW==4.OR.ISW==5)THEN
          ADUM='SOIL35K.DAT'
          CALL OPENV(KR(28),ADUM,IDIR)
          READ(KR(28),140)XAV,XDV,XRG,BRNG,NSX
          ALLOCATE(XSP(NSX,5))
		  NSNN=.655*NSX**.493
		  EXNN=.767*NSX**.049
	  	  DO I=1,NSX
			  READ(KR(28),140)(XSP(I,J),J=1,5)
          END DO
		  CLOSE(KR(28))
      END IF
      CALL OPENV(KR(2),FPARM,0)
      CALL OPENV(KR(17),FPRNT,0)
      CALL OPENV(KR(6),FMLRN,0)
      CALL OPENV(KR(18),FWPM,IDIR)
      CALL OPENV(KR(19),FWIND,IDIR)
      CALL OPENV(KR(4),FCROP,IDIR)
      CALL OPENV(KR(3),FTILL,IDIR)
      CALL OPENV(KR(8),FPEST,IDIR)
      CALL OPENV(KR(9),FFERT,IDIR)
      CALL OPENV(KR(10),FTR55,IDIR)      
      CALL OPENV(KR(13),FSOIL,IDIR)
      CALL OPENV(KR(15),FOPSC,IDIR)      
	  IF(IHRD>0)CALL OPENV(KR(21),FHERD,IDIR) 
      CALL OPENV(KR(23),FSITE,IDIR)
      CALL OPENV(KR(24),FSUBA,IDIR)
	  CALL OPENV(KR(25),FWLST,IDIR)
	  CALL OPENV(KR(27),FPSOD,IDIR)
	  IRUN=0
	  IDAY=0
	  IMON=0
	  IT1=0
	  IT2=0
	  IT3=0
	  IYER=0
	  NMC=0
      ! SCRP = S CURVE SHAPE PARAMETERS(CONSTANTS EXCEPT FOR
      !        EXPERIMENTAL PURPOSES)
      ! LINE 1/35
      READ(KR(2),2510)((SCRP(I,J),J=1,2),I=1,35)
      ! MISCELLANEOUS PARAMETERS (CONSTANTS EXCEPT FOR EXPERIMENTAL PURPO
      ! LINE 36/41
      READ(KR(2),3120)PRMT
      IF(PRMT(96)==0.)PRMT(96)=1.
      ! READ ECONOMIC DATA
      ! 1  COIR = COST OF IRRIGATION WATER ($/mm)
      ! 2  COL  = COST OF LIME ($/t)
      ! 3  FULP = COST OF FUEL ($/l)
      ! 4  WAGE = LABOR COST ($/h)
      ! LINE 42
      READ(KR(2),3120)COIR,COL,FULP,WAGE
      ! LINE 43
      READ(KR(2),3120)XKN5,XKN3,XKN1,CBVT                              
      CLOSE(KR(2))
      ! CQRB  = COEFS OF 7TH DEG POLY IN TR55 QRB EST
      ! LINE 1/68
      READ(KR(10),396)CQRB
      ! COCH = COEFS FOR CHANNEL GEOMETRY X=COCH(N)*WSA**COCH(N+1)
      ! X=WIDTH FOR COCH(1) & COCH(2)
      ! X=DEPTH FOR COCH(3) & COCH(4)
      ! X=LENGTH FOR COCH(5) & COCH(6)
      ! LINE 69
      READ(KR(10),3120)COCH
      IF(COCH(1)<1.E-10)COCH(1)=.0814
      IF(COCH(2)<1.E-10)COCH(2)=.6
      IF(COCH(3)<1.E-10)COCH(3)=.0208
      IF(COCH(4)<1.E-10)COCH(4)=.4
      IF(COCH(5)<1.E-10)COCH(5)=.0803
      IF(COCH(6)<1.E-10)COCH(6)=.6
      CLOSE(KR(10))
      !     KA   = OUTPUT VARIABLE ID NO(ACCUMULATED AND AVERAGE VALUES)
      !     LINE 1/5
      READ(KR(17),3100)(KDC1(I),I=1,NKA)
      DO I=1,NKA
          IF(KDC1(I)<=0)EXIT
          KA(I)=KDC1(I)
      END DO
      NKA=I-1
      ! JC = OUTPUT VARIABLE ID NO(CONCENTRATION VARIABLES)
      ! LINE 6
      READ(KR(17),3100)(KDC1(I),I=1,NJC)
      DO I=1,NJC
          IF(KDC1(I)<=0)EXIT
          JC(I)=KDC1(I)
      END DO
      NJC=I-1
      ! KS = OUTPUT VARIABLE ID NO(STATE VARIABLES)
      ! LINE 7
      READ(KR(17),3100)(KDC1(I),I=1,NKS)
      DO I=1,NKS
          IF(KDC1(I)<=0)EXIT
          KS(I)=KDC1(I)
      END DO
      NKS=I-1
      ! KD = DAILY OUTPUT VARIABLE ID NO
      ! LINE 8/9
      READ(KR(17),3100)(KDC1(I),I=1,NKD)
      DO I=1,NKD
          IF(KDC1(I)<=0)EXIT
          KD(I)=KDC1(I)
      END DO
      NKD=I-1
      ! KY = ANNUAL OUTPUT--VARIABLE ID NOS(ACCUMULATED AND AVERAGE
      ! LINE 10/11
      READ(KR(17),3100)(KDC1(I),I=1,NKY)
      DO I=1,NKY
          IF(KDC1(I)<=0)EXIT
          KY(I)=KDC1(I)
      END DO
      NKY=I-1
      ! SELECT OUTPUT FILES--KFL=0(NO OUTPUT); KFL>0(GIVES OUTPUT FOR
      ! SELECTED FILE)
      ! 1  OUT=STANDARD OUTPUT FILE
      ! 2  MAN=SPECIAL MANURE MANAGEMENT SUMMARY FILE
      ! 3  SUS=SUBAREA SUMMARY FILE
      ! 4  ASA=ANNUAL SUBAREA FILE
      ! 5  SWT=WATERSHED OUTPUT TO SWAT
      ! 6  DPS=DAILY SUBAREA PESTICIDE FILE
      ! 7  MSA=MONTHLY SUBAREA FILE
      ! 8  AWP=ANNUAL CEAP FILE
      ! 9  DHY=DAILY SUBAREA HYDROLOGY FILE
      !10  WSS=WATERSHED SUMMARY FILE
      !11  SAD=DAILY SUBAREA FILE
      !12  HYC=CONTINUOUS HYDROGRAPHS(DTHY) AT WATERSHED OUTLET
      !13  DRS=DAILY RESERVOIR FILE
      !14  APEXBUF.OUT SPECIAL FILE FOR BUFFER STRIPS
      !15  MWS=MONTHLY WATERSHED FILE
      !16  DWS=DAILY WATERSHED OUTLET FILE
      !17  AWS=ANNUAL WATERSHED OUTLET FILE
      !18  DGZ=DAILY GRAZING
      !19  DUX=DAILY MANURE APPLICATION
      !20  DDD=DAILY DUST DISTRIBUTION
      !21  ACN=ANNUAL SOIL ORGANIC C & N TABLE
      !22  DCN=DAILY SOIL ORGANIC C & N TABLE
      !23  SCX=SUMMARY SOIL ORGANIC C & N TABLE
      !24  ACY=ANNUAL SUBAREA CROP YIELD
      !25  EFR=RUNOFF EVENT FLOOD ROUTING
      !26  EHY=RUNOFF EVENT HYDROGRAPHS
      !27  APS=ANNUAL SUBAREA/WATERSHED PESTICIDE
      !28  MSW=MONTHLY OUTPUT TO SWAT.
      !29  DPW=DAILY WATERSHED PESTICIDE FILE
      !30  SPS=PESTICIDE SUBAREA SUMMARY
      !31  ACO=ANNUAL COST 
      !32  SWN=SPECIAL WATERSHED SUMMARY FOR NRCS FARM PLANNING
      !33  NOT USED
      !34  SAO=SPECIAL SUBAREA FILE FOR GIS
      !35  RCH=SPECIAL REACH FILE FOR GIS
      !36  ERX=ERROR FILE
      !37  DMR=DAILY WATERSHED NUTRIENT & SEDIMENT CONC NRCS MRBI
      !38  STR=SUMMARY OF SUBARES & WATERSHED FOR STAR
      !39  MRH=MONTHLY REACH FILE ANNUAL GIS REACH FILE FOR SELECTED COMMAND #'S ICMO(FROM .SIT)
      !40  MGZ=MONTHLY GRAZING FILE
      !41  DNC=DAILY NITROGEN/CARBON CESAR IZAURRALDE 
      !42  DHS=DAILY HYDROLOGY/SOIL                                                                                                     
      !43  MSX=MONTHLY SOIL TABLE
      !44  DGN=DAILY GENERAL OUTPUT (VAR AFTER COMMAND LOOP IN BSIM)
      !45  DPD=DAILY PADDY OUTPUT
      !46  ASL=ANNUAL SOIL TABLE
      !47  MS5=MONTHLY SOIL PROPERTIES 0-0.05m
      !48  AS5=ANNUAL SOIL PROPERTIES 0-0.05m 
      !49  RUN1501.SUM
      !    MSO
      !    +1/ SOT=SUBAREA FINAL SOIL TABLE FOR USE IN OTHER RUNS
      !    MSA
      !     LINE 12/14
      READ(KR(17),3100)KFL0
      CLOSE(KR(17))
!  1  RFN0 = AVE CONC OF N IN RAINFALL(ppm)
!  2  CO20 = CO2 CONCENTRATION IN ATMOSPHERE(ppm)
!  3  CQN0 = CONC OF N IN IRRIGATION WATER(ppm)
!  4  PSTX = PEST DAMAGE SCALING FACTOR (0.-10.)--0. SHUTS OFF PEST
!            DAMAGE FUNCTION. PEST DAMAGE FUNCTION CAN BE REGULATED FROM
!            VERY MILD(0.05-0.1) TO VERY SEVERE(1.-10.)
!  5  YWI  = NO Y RECORD MAX .5H RAIN(0 IF WI IS NOT
!            INPUT)
!  6  BTA  = COEF(0-1)GOVERNING WET-DRY PROBABILITIES GIVEN DAYS
!            OF RAIN(0 IF UNKNOWN OR IF W|D PROBS ARE
!            INPUT)
!  7  EXPK = PARAMETER USED TO MODIFY EXPONENTIAL RAINFALL AMOUNT
!            DISTRIBUTION(0 IF UNKNOWN OR IF ST DEV & SK CF ARE
!            INPUT)
!  8  QG   = 2 YEAR FREQ 24-H RAINFALL (mm)--ESTIMATES REACH CH GEOMETRY
!            IF UNKNOWN(0 IF CH GEO IS INPUT)
!  9  QCF  = EXPONENT IN WATERSHED AREA FLOW RATE EQ
! 10  CHS0 = AVE UPLAND SLOPE (m/m) IN WATERSHED
! 11  BWD  = CHANNEL BOTTOM WIDTH/DEPTH(QG>0)
! 12  FCW  = FLOODPLAIN WIDTH/CHANNEL WIDTH
! 13  FPS0 = FLOODPLAIN SAT COND ADJUSTMENT FACTOR(.0001_10.)
! 14  GWS0 = MAXIMUM GROUNDWATER STORAGE(mm)
! 15  RFT0 = GROUNDWATER RESIDENCE TIME(d)
! 16  RFP0 = RETURN FLOW / (RETURN FLOW + DEEP PERCOLATION)
! 17  SAT0 = SATURARTED CONDUCTIVITY ADJUSTMENT FACTOR(.1_10.)
! 18  FL0  = FIELD LENGTH(km)(0 IF UNKNOWN)
! 19  FW0  = FIELD WIDTH(km)(0 IF UNKNOWN)
! 20  ANG0 = CLOCKWISE ANGLE OF FIELD LENGTH FROM NORTH(deg)(0 IF
!            UNKNOWN)
! 21  UXP  = POWER PARAMETER OF MODIFIED EXP DIST OF WIND SPEED(0
!            IF UNKNOWN)
! 22  DIAM = SOIL PARTICLE DIAMETER(UM)(0 IF UNKNOWN)
! 23  ACW  = WIND EROSION CONTROL FACTOR
!            = 0.0 NO WIND EROSION
!            = 1.0 FOR NORMAL SIMULATION
!            > 1.0 ACCELERATES WIND EROSION(CONDENSES TIME)
! 24  GZL0 = GRAZING LIMIT--MINIMUM PLANT MATERIAL(t/ha)
! 25  RTN0 = NUMBER YEARS OF CULTIVATION AT START OF SIMULATION
! 26  BXCT = LINEAR COEF OF CHANGE IN RAINFALL FROM E TO W(PI/P0/KM)
! 27  BYCT = LINEAR COEF OF CHANGE IN RAINFALL FROM S TO N(PI/P0/KM)
! 28  DTHY = TIME INTERVAL FOR FLOOD ROUTING(h)
! 29  QTH  = ROUTING THRESHOLD(mm)--VSC ROUTING USED ON QVOL>QTH
! 30  STND = VSC ROUTING USED WHEN REACH STORAGE >STND
! 31  DRV  = SPECIFIES WATER EROSION DRIVING EQ.
!            (1=RUSL2;  2=USLE;  3=MUSS;  4=MUSL;   5=MUST;  6=REMX)
! 32  PCO0 = FRACTION OF SUBAREAS CONTROLLED BY PONDS.
! 33  RCC0 = REACH CHANNEL USLE C FACTOR
! 34  CSLT = SALT CONC IN IRRIGATION WATER (ppm)
! 35  CPV0 = FRACTION INFLOW PARTITIONED TO VERTICLE CRACK OR PIPE FLOW
! 36  CPH0 = FRACTION INFLOW PARTITIONED TO HORIZONTAL CRACK OR PIPE FLOW 
! 37  DZDN = LAYER THICKNESS FOR DIFFERENTIAL EQ SOLN TO GAS DIFF EQS(m)
! 38  DTG  = TIME INTERVAL FOR GAS DIFF EQS (h)
! 39  QPQ  = RATIO VOLUME BEFORE TP TO TOTAL VOLUME OF UNIT HYD
!     LINES 3/6
      READ(KR(20),*)RFN0,CO20,CQN0,PSTX,YWI,BTA,EXPK,QG,QCF,CHS0,&
      BWD,FCW,FPS0,GWS0,RFT0,RFP0,SAT0,FL0,FW0,ANG0,UXP,DIAM,ACW,GZL0,&
      RTN0,BXCT,BYCT,DTHY,QTH,STND,DRV,PCO0,RCC0,CSLT,CPV0,CPH0,DZDN,&
      DTG,QPQ
      IF(DTHY<1.E-10)DTHY=1.
      IF(IHY>0)NPD=24./DTHY+1.
      IF(BTA<1.E-10)BTA=.75
      DO I=1,NSZ
          PSZX(I)=.411*PSZ(I)*PSZ(I)
          PSZY(I)=SQRT(PSZ(I))
      END DO
      ZMIX=PRMT(31)
!  1  ASTN = RUN #
!  2  ISIT = SITE #
!  3  IWPM = WEATHER STATION # FROM KR(18) 
!  4  IWND = WIND STATION # FROM KR(19) 
!  5  ISUB = SUBAREA # FROM KR(24)
!  6  ISOL = 0 FOR NORMAL RUN
!          > 0 FOR RUNS USING .SOT FILES
!  7  IRFT = WITHIN STORM RAINFALL STATION # FROM KR(29) 

      !rtb MODFLOW: initialize MODFLOW linkage
      IF(IMF>0) call amrt_init_mf
      DO ! RUN LOOP
          READ(KR(11),*,IOSTAT=NFL)ASTN,ISIT,IWPM(1),IWND,ISUB,ISOL,IRFT
          IF(NFL/=0)EXIT
          IF(ISIT==0)EXIT
          ASTN=ADJUSTR(ASTN)
          I2=-1
          DO WHILE(I2/=ISUB) 
	          READ(KR(24),*,IOSTAT=NFL)I2,SAFILE
              IF(NFL/=0)THEN
                  WRITE(*,*)'SUBAREA NO = ',ISUB,' NOT IN SUBAREA LIST FILE'
                  STOP
              END IF
          END DO
          REWIND KR(24)
          I2=-1
          DO WHILE(I2/=ISIT) 
              READ(KR(23),*,IOSTAT=NFL)I2,SITEFILE
              IF(NFL/=0)THEN
                  WRITE(*,*)'SITE NO = ',ISIT,' NOT IN SITE LIST FILE'
                  STOP
              END IF
          END DO
          REWIND KR(23)
          CALL OPENV(KR(1),SITEFILE,IDIR)
	      CALL OPENV(KR(5),SAFILE,IDIR)
          WRITE(*,2143)IRUN,SAFILE
          MOW=1
          NBMX=0
          I=0     
          DO 
              I=I+1 
              READ(KR(5),'(2I8)',IOSTAT=NFL)N1,ISA
              IF(NFL/=0.OR.N1==0)EXIT
              READ(KR(5),*)INPS,IOPS,IOW,I1,I2,I3,I4,I5,I6
              IF(IOW>MOW)MOW=IOW
              IF(N1>NBMX)NBMX=N1
              DO J=1,10
                  READ(KR(5),'()')
              END DO
          END DO
          MSA=I-1
          XSA=MSA
          REWIND KR(5)
          I=1
          DO
              READ(KR(27),'(A80)',IOSTAT=NFL)TITWN
              IF(NFL/=0.OR.TITWN=='')EXIT
              I=I+1
          END DO
          MPO=I-1
          REWIND KR(27)
          CALL ALLOCATE_PARMS
          IF(IRUN==0)NBMX0=NBMX
          HED=(/" TMX"," TMN","SRAD","PRCP","SNOF","SNOM","WSPD","RHUM",&
          " VPD"," PET","  ET","  EP","   Q","  CN"," SSF"," PRK"," QDR",&
          "IRGA"," QIN","TLGE","TLGW","TLGQ","TLGF","  EI","   C","USLE",&
          "MUSL","REMX","MUSS","MUST","RUS2"," WK1","RHTT","RRUF","RGRF",&
          "YWND","  YN","  QN","SSFN","PRKN"," GMN","  DN","NFIX"," NMN",&
          "NITR","AVOL","QDRN","  YP","  QP"," MNP","PRKP","  ER"," FNO",&
          "FNMN","FNMA"," FPO"," FPL","LIME"," TMP","SW10","LGMI","LGMO",&
          " EPP","RSQI","RSQO","RSEV","RSLK","RSYI","RSYO","RSYD","DPRK",&
          "RSSF","RSDC","RSPC","PRKC","  QC","  YC","RSDA"," QFP","RSFN",&
          " MAP","BUNL"," QRF","QRFN","RFIC","RSBK","CPVH","YMNU","SNOU",&
          "SPOU","DNMO","DPMO","DEMR","P10D","SSFI","DPKN","CPVV"," FPF",&
          " FOC"," RFV","SCOU","DEPC","DECR","PSOQ","PSON","PSOP","    ",&
          " QPU","FALF","IRDL"," QRP"," YRP","YNRP","YPRP","QNRP","QPRP",&
          "WYLD"," YPM"," YPO","  SW","PSOY","PQPS","PYPS","    ","  QI",&
          "QARS","RFRA"," DN2","SLTI","SLTQ","SLTS","SLTF","SLTV","YNWN",&
          "YPWN","YCWN","PSO3","PSSP","YWKS","CBUR","GRZD","QRFP","QDRP",&
          "YTHS","YWTH","RSIR","EVRT"," FSK","YEFK"," QSK"," SSK"," VSK",&
          "SKOU","DN2O","ETM3","WLIR","GWTR","SMQN"/)
          HEDH=(/"  QI","  QO","  ET"," FPF","  YI","  YO","  CY","YONI",&
          "YONO","YOPI","YOPO","NO3I","NO3O","NH4I","NH4O","NO2I","NO2O",&
          " QPI"," QPO","ALGI","ALGO","BODI","BODO","DO2I","DO2O","QPSI",&
          "QPSO","YPSI","YPSO","RPST","VPST","DPST","WYLI","WYLO"," RFV"/)
          MASA=(/0,0,0,2,2,2,0,0,0,2, & !10  ! array for unit conversion of outputs 
                 0,0,2,0,2,2,2,2,2,2, & !20  !    fixed by Jaehak 2020:  143,62, 84, 88
                 2,2,2,0,0,1,1,1,1,1, & !30
                 1,0,0,0,0,1,1,1,1,1, & !40
                 0,1,1,0,1,1,1,1,1,1, & !50
                 1,0,1,1,1,1,1,0,0,0, & !60
                 0,0,0,2,2,2,2,1,1,1, & !70
                 2,2,1,1,1,1,1,0,0,1, & !80
                 1,0,2,1,0,0,0,1,0,0, & !90
                 0,0,0,0,0,1,0,0,0,0, & !100
                 0,0,0,0,0,0,0,0,0,2, & !110
                 2,0,0,0,0,0,2,1,1,0, & !120
                 0,0,0,0,0,0,0,1,0,0, & !130
                 0,0,0,1,1,1,0,0,0,0, & !140
                 0,0,1,0,0,0,2,0,0,0, & !150
                 0,0,0,1,2,0,0,0/)      !158

          DO I=1,MHX
	          IHX(I)=I
          END DO		 
          J=101
          J1=2*MSA+MSO       	 
          DO I=1,J1
	          KW(I)=J
     	      J=J+1
          END DO
          KFL=KFL0
          IF(KFL(14)>0)OPEN(KW(14),FILE='APEXBUF.OUT')
          CALL OPENF(ASTN)
	      NBYR=NBY0
	      NGN=NGN0
          DO I=1,60
              TITLE(I)='    '
          END DO
          IF(KFL(1)>0)WRITE(KW(1),508)'FSITE',FSITE,'FSUBA',FSUBA,'FWPM',&
          FWPM,'FWIND',FWIND,'FCROP',FCROP,'FTILL',FTILL,'FPEST',FPEST,&
          'FFERT',FFERT,'FSOIL',FSOIL,'FOPSC',FOPSC,'FTR55',FTR55,&
          'FPARM',FPARM,'FMLRN',FMLRN,'FPRNT',FPRNT,'FHERD',FHERD,&
          'FWLST',FWLST,'FPSOD',FPSOD,'FRFDT',FRFDT
          KGN=0
	      CALL AINLZ
          CALL AINIX
	      DO J=1,NRN0
              IDG(J)=J
          END DO
          IF(KFL(49)<0)THEN
              OPEN(KW(49),FILE='RUN1501.SUM')
          ELSE
              OPEN(KW(49),FILE='RUN1501.SUM',ACCESS='APPEND')
          END IF
          IRUN=IRUN+1
          IF(IRUN==1)THEN
              CALL ATIMER(0)
              CALL APAGE(0)
              IF(KFL(1)>0)WRITE(KW(1),2570)
              IF(KFL(49)<0)THEN
                  WRITE(KW(49),'(//10X,A)')'APEX1501 RUN1501.SUM'
                  WRITE(KW(49),104)HED(4),HED(13),HED(16),HED(24),&
                  HED(25),HED(26),HED(27),HED(28),HED(29),HED(30)
              END IF    
              DO I=1,2
                  DO J=1,35
                      SCRX(J,I)=SCRP(J,I)
                  END DO
              END DO        
              DO I=1,34
                  IF(SCRP(I,1)<1.E-10)CYCLE
                  X1=ASPLT(SCRP(I,1))
                  X2=ASPLT(SCRP(I,2))
                  IF(I==13)THEN
                      SRNG(I)=SCRP(I,2)-SCRP(I,1)
                      SCMN(I)=SCRP(I,1)
                      SCRP(I,1)=.01
                      SCRP(I,2)=.99
                  END IF    
                  CALL ASCRV(SCRP(I,1),SCRP(I,2),X1,X2,SCLM,I,KW(1))
              END DO
              DO J=1,35
                  IF(KFL(1)>0)WRITE(KW(1),2540)J,PRMT(J),(SCRX(J,I),I=1,2),(SCRP(J,I),I=1,2)
              END DO
              DO J=31,SIZE(PRMT)
                  IF(KFL(1)>0)WRITE(KW(1),2540)J,PRMT(J)
              END DO
              LGRZ=PRMT(60)
          END IF
          RHS=PRMT(100)
          RHP=PRMT(101)
          CRLNC=PRMT(102) 
          CRUNC=PRMT(103)
          WKA=PRMT(104)
          WNCMIN=PRMT(105)
          WNCMAX=PRMT(106)
          VMU=PRMT(107)
          WKMNH3=PRMT(108)
          WKMNO2=PRMT(109)
          WKMNO3=PRMT(110)                                                                        
    !     TITLE = PROBLEM DESCRIPTION( 3 LINES)
    !     LINES 1/3
          READ(KR(1),3090)TITLE
    !  1  YLAT = LATITUDE(deg)
    !  2  XLOG = LONGITUDE(deg)
    !  3  ELEV = ELEVATION OF WATERSHED (m)
    !  4  APM  = PEAK RATE - EI ADJUSTMENT FACTOR (BLANK IF UNKNOWN)
    !  5  CO2X = CO2 CONCENTRATION IN ATMOSPHERE (ppm)--NON ZERO VALUE
    !            OVERRIDES CO2 INPUT IN APEXCONT.DAT
    !  6  CQNX = CONC OF NO3 IN IRRIGATION WATER (ppm)--NON ZERO VALUE
    !            OVERRIDES CNO30 INPUT IN APEXCONT.DAT 
    !  7  RFNX = AVE CONC OF N IN RAINFALL (ppm)
    !  8  UPR  = MANURE APPL RATE TO SUPPLY P UPTAKE RATE (kg/ha/y)
    !  9  UNR  = MANURE APPL RATE TO SUPPLY N UPTAKE RATE (kg/ha/y)
    ! 10  FIR0 = FACTOR TO ADJUST AUTO IRRIGATION VOLUME (FIRG*FC)
    !     LINE 4
          READ(KR(1),3120)YLAT,XLOG,ELEV,APM,CO2X,CQNX,RFNX,UPR,UNR,FIR0
    ! 1/5 X1/X5= DUMMY  
    !  6  GWSP = GROUND WATER SLOPE (m/m)
    ! 7/8 X7/X8= DUMMY
    !  9  BCHL = SWAT BASIN CHANNEL LENGTH(km) 
    ! 10  BCHS = SWAT BASIN CHANNEL SLOPE(m/m)
    !     LINE 5
          READ(KR(1),3120)X1,X2,X3,X4,X5,GWSP,X7,X8,BCHL,BCHS
    !     LINE 6
          READ(KR(1),'()')
	      CO2=CO20
	      IF(CO2X>0.)CO2=CO2X
	      CQNI=CQN0
	      IF(CQNX>0.)CQNI=CQNX
          RFNC=RFN0
          IF(RFNX>0.)RFNC=RFNX
	      YLTR=YLAT/CLT
	      SIN1=SIN(YLTR)
          COS1=COS(YLTR)
          TAN1=TAN(YLTR)
          IF(DTG<1.E-10)DTG=1.
          NBDT=24./DTG+.9 
          CH=.4349*ABS(TAN1)
          IF(CH>=1.)THEN
              H=0.
          ELSE
              H=ACOS(CH)
          END IF
          CH=-CH
          IF(CH>=1.)THEN
              H1=0.
          ELSE
              IF(CH<=-1.)THEN
                  H1=3.1416
              ELSE
                  H1=ACOS(CH)
              END IF
          END IF  
          HLMX=7.72*H1
	      DO ISA=1,MSA
	          YTN(ISA)=TAN1
              YLS(ISA)=SIN1
              YLC(ISA)=COS1
              HLMN(ISA)=7.72*H
              DHL=.1*(HLMX-HLMN(ISA))
              JT1=1
              IF(YLAT<0.)JT1=6
              JDHU=400
              IF(HLMN(ISA)<=11.)CALL ADAJ(NC,JDHU,JT1,15,1)
              WDRM(ISA)=PRMT(6)+HLMN(ISA)
              HLMN(ISA)=HLMN(ISA)+MAX(DHL,PRMT(6))      
          END DO    
          IF(IWND>0)THEN
	          I2=-1
	          DO WHILE(I2/=IWND) 
                  READ(KR(19),*,IOSTAT=NFL)I2,WINDFILE
                  IF(NFL/=0)THEN
                      WRITE(*,*)'WIND NO = ',IWND,' NOT IN WIND LIST FILE'
                      STOP
                  END IF              
              END DO
          ELSE
              D0=1.E20
              DO 
                  READ(KR(19),*,IOSTAT=NFL)I2,OPSCFILE,Y,X
                  IF(NFL/=0)EXIT
                  RY=Y/CLT
	              XX=SIN1*SIN(RY)+COS1*COS(RY)*COS((X-XLOG)/CLT)
                  D=6378.8*ACOS(XX)
                  IF(D>D0)CYCLE
                  D0=D
                  WINDFILE=OPSCFILE
              END DO
          END IF
          REWIND KR(19)
          CALL OPENV(KR(22),WINDFILE,IDIR)
	      TITWN=WINDFILE
          READ(KR(22),'()')
          READ(KR(22),'()')
          IYR=IYR0
          IDA=IDA0
          IBDT=IDA+100*IMO+10000*IYR
          CALL AISPL(IPD,INP)
          IF(INP>0)THEN
              NOP=0
          ELSE
              NOP=1
              INP=1
          END IF
          IF(IPD<=5)IPYI=INP
          CALL APAGE(0)
          IF(KFL(1)>0)THEN
              WRITE(KW(1),3540)
              WRITE(KW(1),'(/1X,A)')'-----VARIABLE NAMES & UNIT CONVERSION CODE'
              WRITE(KW(1),2410)
              WRITE(KW(1),'(20(1X,A4,I2))')(HED(I),MASA(I),I=1,NSM)
          END IF    
          NDRV=DRV
          IF(YWI<1.E-10)YWI=10.
          IF(GWS0<1.E-5)GWS0=50.
          IF(RFT0<1.E-5)RFT0=10.
          IF(RFP0<1.E-5)RFP0=.5
          CALL APAGE(0)
          IF(KFL(1)>0)THEN
              WRITE(KW(1),3540)
              WRITE(KW(1),3110)NBYR,IYR,IMO,IDA
              IF(LPYR>0)THEN
                  WRITE(KW(1),'(T10,A)')'LEAP YEAR IGNORED'
              ELSE
                  WRITE(KW(1),'(T10,A)')'LEAP YEAR CONSIDERED'
              END IF
              WRITE(KW(1),3730)YLAT,XLOG,ELEV
          END IF    
          SELECT CASE(NDRV)
	          CASE(1)
                  NDVSS=31
	          CASE(2)
                  NDVSS=26
	          CASE(3)
	              NDVSS=29
	          CASE(4)
                  NDVSS=27
	          CASE(5)
                  NDVSS=30
	          CASE(6)
	              NDVSS=28
	      END SELECT
          IF(KFL(1)>0)THEN    
              WRITE(KW(1),'(T10,A,A4)')'WATER EROSION FACTORS--DRIVING EQ = ',&
              HED(NDVSS)
              WRITE(KW(1),'(T10,A)')'DAILY RUNOFF ESTIMATION'
              SELECT CASE(INFL)
              CASE(0)
                  WRITE(KW(1),'(T15,A)')'NRCS CURVE NUMBER EQ'
              CASE(1)
                  WRITE(KW(1),'(T15,A/T15,A)')'GREEN & AMPT EQ','RF EXP DST&
                  --PEAK RF RATE SIM'
              CASE(2)
                  WRITE(KW(1),'(T15,A/T15,A)')'GREEN & AMPT EQ','RF EXP DST&
                  --PEAK RF RATE INPUT'
              CASE(3)
                  WRITE(KW(1),'(T15,A/T15,A)')'GREEN & AMPT EQ','RF UNIF DST--&
                  PEAK RF RATE INP'
              CASE(4)
                  WRITE(KW(1),'(T15,A/T15,A,F6.3,A)')'GREEN & AMPT EQ','RF INPUT AT &
                  INTERVAL = ',DTHY,' h'
              END SELECT
              IF(ISCN==0)THEN
                  WRITE(KW(1),'(T15,A)')'DAILY CN--STOCHASTIC'
              ELSE
                  WRITE(KW(1),'(T15,A)')'DAILY CN--DETERMINISTIC'
              END IF
              SELECT CASE(ITYP)
                  CASE(1,2,3,4)
                      WRITE(KW(1),'(T10,A,A2)')'PEAK RATE EST WITH TR55--RF TYPE =',&
                      RFPT(ITYP)
                  CASE(-1)
                      WRITE(KW(1),'(T10,A)')'PEAK RATE MOD RATIONAL EQ - RIGID'
                  CASE(0)
                      WRITE(KW(1),'(T10,A)')'PEAK RATE MOD RATIONAL EQ - &
                      STOCHASTIC'
                  CASE(5)    
                      WRITE(KW(1),'(T10,A,F6.3)')'PEAK RATE SCS UNIT &
                      HYD EQ QPQ = ',QPQ
              END SELECT
          END IF    
          IF(QG>0.)THEN
              X2=(QG-8.96)**2/(QG+35.86)
              X1=.1732
              TCX=.0663*X1**.77/CHS0**.385
              QG0=.1*X2/TCX
              IF(KFL(1)>0)WRITE(KW(1),401)QG,X2,QCF,X1,CHS0,TCX,QG0,BWD,&
              FCW,FPS0,GWS0,RFT0,RFP0,PCO0,RCC0
          END IF
          SELECT CASE(IHY)
              CASE(1)
                  IF(KFL(1)>0)WRITE(KW(1),'(T15,A)')'VSC METHOD'
                  IF(KFL(12)>0)WRITE(KW(12),'(T15,A)')'VSC METHOD'
                  IF(KFL(26)>0)WRITE(KW(26),'(T15,A)')'VSC METHOD'
              CASE(2)
                  IF(KFL(1)>0)WRITE(KW(1),'(T15,A)')'SVS METHOD'
                  IF(KFL(12)>0)WRITE(KW(12),'(T15,A)')'SVS METHOD'
                  IF(KFL(26)>0)WRITE(KW(26),'(T15,A)')'SVS METHOD'
              CASE(3)    
                  IF(KFL(1)>0)WRITE(KW(1),'(T15,A)')'M_CVC METHOD'
                  IF(KFL(12)>0)WRITE(KW(12),'(T15,A)')'M_CVC METHOD'
                  IF(KFL(26)>0)WRITE(KW(26),'(T15,A)')'M_CVC METHOD'
              CASE(4)
                  IF(KFL(1)>0)WRITE(KW(1),'(T15,A)')'M_CVC4 METHOD'    
                  IF(KFL(12)>0)WRITE(KW(12),'(T15,A)')'M_CVC4 METHOD'    
                  IF(KFL(26)>0)WRITE(KW(26),'(T15,A)')'M_CVC4 METHOD'    
              CASE DEFAULT
                  IF(KFL(1)>0)WRITE(KW(1),'(T10,A)')'RUNOFF VOLUME ROUTING'
           END SELECT
           IF(KFL(1)>0)THEN    
              IF(IPRK==0)THEN
                  WRITE(KW(1),'(T10,A)')'PERCOLATION--HPERC'
              ELSE
                  WRITE(KW(1),'(T10,A)')'PERCOLATION--HPERC1 (4mm SLUG FLOW)'
              END IF
              IF(IREM==0)THEN
                  WRITE(KW(1),'(T10,A)')'SSK FROM REMX'
              ELSE
                  WRITE(KW(1),'(T10,A)')'SSK FROM USLE'
              END IF        
              IF(IERT>0)THEN
                  WRITE(KW(1),'(T10,A)')'GLEAMS ENRICHMENT RATIO'
              ELSE
                  WRITE(KW(1),'(T10,A)')'EPIC ENRICHMENT RATIO'
              END IF
              IF(LBP==0)THEN   
                  WRITE(KW(1),'(T10,A)')'GLEAMS PESTICIDE EQ SOL P RUNOFF'
              ELSE
                  WRITE(KW(1),'(T10,A)')'LANGMUIR EQ SOL P RUNOFF'    
	          END IF
              IF(NUPC>0)THEN
                  WRITE(KW(1),'(T10,A)')'N & P UPTAKE CONC S CURVE'
              ELSE
                  WRITE(KW(1),'(T10,A)')'N & P UPTAKE CONC SMITH CURVE'
              END IF
              IF(ICP==0)THEN
                  WRITE(KW(1),'(T10,A)')'PHOENIX N-C MODEL'
              ELSE    
                  WRITE(KW(1),'(T10,A)')'CENTURY N-C MODEL'
              END IF    
              IF(IOX==0)THEN
                  WRITE(KW(1),'(T10,A)')'EPIC O2=F(Z)'
              ELSE
                  WRITE(KW(1),'(T10,A)')'O2=F(C/CLA)'
              END IF
              IF(NTV==0)THEN
                  WRITE(KW(1),'(T10,A)')'ORIGINAL APEX NITVOL EQS'
              ELSE                                          
                  WRITE(KW(1),'(T10,A)')'IZAURRALDE REVISED NITVOL EQS'
              END IF    
              SELECT CASE(IDNT)
                  CASE(1)	          
                      WRITE(KW(1),'(T10,A)')'EPIC DNIT'
                  CASE(2)              
	                  WRITE(KW(1),'(T10,A)')'KEMANIAN DNIT'                                                                                              	                                                                                     
                  CASE(3)
                      WRITE(KW(1),'(T10,A,F5.2,A)')'IZAURRALDE DNIT DZDN=',DZDN,' m'
                      WRITE(KW(1),'(T15,3(A,E16.6)A)')'XKN1=',XKN1,' XKN3=',XKN3,&
                      ' XKN5=',XKN5,' ORIGINAL DW'                                                          
                      IF(JRRS==0.)THEN
                          WRITE(KW(1),'(T10,A)')'ROOT RESPIRATION CALCULATED IN NDNITCI'
                      ELSE                                                              
                          WRITE(KW(1),'(T10,A)')'ROOT RESPIRATION NOT CONSIDERED IN NDNITCI'                       
                      END IF
                  CASE(4)
                      WRITE(KW(1),'(T10,A,F5.2,A)')'IZAURRALDE DNIT DZDN=',DZDN,' m'
                      WRITE(KW(1),'(T15,3(A,E16.6)A)')'XKN1=',XKN1,' XKN3=',XKN3,&
                      ' XKN5=',XKN5,' NEW DW'
                      IF(JRRS==0.)THEN
                          WRITE(KW(1),'(T10,A)')'ROOT RESPIRATION CALCULATED IN NDNITCI'
                      ELSE                                                              
                          WRITE(KW(1),'(T10,A)')'ROOT RESPIRATION NOT CONSIDERED IN NDNITCI'                       
                      END IF
                  CASE DEFAULT
                      WRITE(KW(1),'(T10,A)')'KEMANIAN DNIT'
                      !IDNT=2                                                                                              	                                                                                     
              END SELECT
              IF(ISLF>0)THEN
                  WRITE(KW(1),'(T10,A)')'MUSLE SLOPE/LENGTH FACTOR'
              ELSE
                  WRITE(KW(1),'(T10,A)')'RUSLE SLOPE/LENGTH FACTOR'
              END IF
              IF(PRMT(81)>0)THEN
                  WRITE(KW(1),'(T10,A)')'DYNAMIC TECHNOLOGY'
              ELSE
                  WRITE(KW(1),'(T10,A)')'STATIC TECHNOLOGY'
              END IF
              WRITE(KW(1),412)LPD,MSCP
              WRITE(KW(1),2450)APM,RFNC,CQNI,UPR,UNR,MNUL
	          WRITE(KW(1),330)COL,COIR,FULP,WAGE
              WRITE(KW(1),12)COCH
          END IF    
          RFNC=RFNC/1000.
          CQNI=CQNI*.01
          IF(KFL(1)>0)THEN
              CALL APAGE(0)
              WRITE(KW(1),'(//1X,A/)')'____________________WEATHER DATA_________&
              ______________'
              WRITE(KW(1),3750)CO2
              SELECT CASE(ICO2)
	              CASE(1)
                      WRITE(KW(1),'(T10,A)')'DYNAMIC ATMOSPHERIC CO2'
                  CASE(2)          
                      WRITE(KW(1),'(T10,A)')'ATMOSPHERIC CO2 INPUT'
                  CASE DEFAULT
                      WRITE(KW(1),'(T10,A)')'STATIC ATMOSPHERIC CO2'
              END SELECT
              WRITE(KW(1),2520)YWI
          END IF    
          XY2=.5/YWI
          SUM=0.
          IWPX=IWPM(1)
          ! WEATHER PARM # (FROM WPM1US.DAT) FOR SPATIAL WEATHER GENERATION
          ! LINE 7
          READ(KR(1),811)(IWPM(J),J=1,10)
          DO J=1,10
              IF(IWPM(J)==0)EXIT
          END DO
          NWP=J-1
          ! FRACTION OF WSA REPRESENTED BY IWPM. FOR SPATIAL WEATHER GENERATION
          ! LINE 8
          READ(KR(1),3120)(FWXP(J),J=1,NWP)
          IF(NWP>1)THEN
              DO J=1,NWP
                  SUM=SUM+FWXP(J)
              END DO
              XX=0.
              DO J=1,NWP
                  FWXP(J)=FWXP(J)/SUM
                  AWXP(J)=FWXP(J)+XX
                  XX=AWXP(J)
              END DO
          ELSE
              FWXP(1)=1.
          END IF
          IF(IWPM(1)==0)IWPM(1)=IWPX
          K=35
          NWP=MAX(1,NWP)
          DO IWI=1,NWP
              IF(IWI>1.OR.IWPM(1)>0)THEN
                  I2=-1 
                  DO WHILE(I2/=IWPM(IWI))
                      READ(KR(18),*,IOSTAT=NFL)I2,WPMFILE
                      IF(NFL/=0)THEN
                          WRITE(*,*)'WPM1 NO = ',IWPM(IWI),' NOT IN MO WEATHER LIST FILE'
                          STOP
                      END IF
                  END DO
              ELSE
                  D0=1.E20
                  DO 
	                  READ(KR(18),*,IOSTAT=NFL)I2,OPSCFILE,Y,X
	                  IF(NFL/=0.OR.I2==0)EXIT
	                  RY=Y/CLT
	                  XX=SIN1*SIN(RY)+COS1*COS(RY)*COS((X-XLOG)/CLT)
                      D=6378.8*ACOS(XX)
                      IF(D>D0)CYCLE
                      D0=D
                      WPMFILE=OPSCFILE
                  END DO
              END IF          
              REWIND KR(18)
              CALL OPENV(KR(K),WPMFILE,IDIR)
	          TITWP(IWI)=WPMFILE
              !     LINES 1/2
              READ(KR(K),'()')
              READ(KR(K),'()')
              !  3  OBMX   = AV MO MAX TEMP (c)
              !  4  OBMN   = AV MO MIN TEMP (c)
              !  5  SDTMX  = MO STANDARD DEV MAX TEMP (c)OR EXTREME MAXIMUM TEMP (c)
              !              IF STANDARD DEV IS NOT AVAILABLE (BLANK IF TEMP IS INPUT
              !              DAILY)
              !  6  SDTMN  = MO STANDARD DEV MIN TEMP (c)OR EXTREME MIN TEMP (c)
              !              IF STANDARD DEV IS NOT AVAILABLE (BLANK IF TEMP IS INPUT
              !              DAILY)
              !  7  RMO    = AV MO PRECIPITATION (mm)
              !  8  RST(2) = MONTHLY ST DEV OF DAILY RAINFALL (mm)(BLANK IF UNKNOWN
              !              OR IF DAILY PRECIPITATION IS INPUT)
              !  9  RST(3) = MONTHLY SKEW COEF OF DAILY RAINFALL (BLANK IF UNKNOWN OR
              !              DAILY PRECIPITATION IS INPUT)
              ! 10  PRW(1) = MONTHLY PROBABILITY OF WET DAY AFTER DRY DAY (BLANK IF
              !              UNKNOWN OR IF DAILY PRECIPITATION IS INPUT)
              ! 11  PRW(2) = MONTHLY PROBABILITY OF WET DAY AFTER WET DAY (BLANK IF
              !              UNKNOWN OR IF DAILY PRECIPITATION IS INPUT)
              ! 12  UAVM   = AV NO DAYS OF PRECIPITATION/MO (BLANK IF PRECIP IS
              !              GENERATED AND IF PRW 1&2 ARE INPUT)
              ! 13  WI     = 3 OPTIONS--(1)MO MAX .5 H RAIN FOR PERIOD = YWI (mm)
              !                         (2)ALPHA (MEAN .5 H RAIN/MEAN STORM
              !                             AMOUNT)
              !                         (3)BLANK IF UNKNOWN
              ! 14  OBSL   = AV MO SOL RAD (MJ/m2 OR LY)(BLANK IF UNKNOWN)
              ! 15  RH     = 3 OPTIONS--(1)AV MO RELATIVE HUMIDITY (FRACTION)
              !                         (2)AV MO DEW POINT TEMP deg C
              !                         (3)BLANK IF UNKNOWN
              !              USED IN PENMAN OR PENMAN-MONTEITH EQS
              ! 16  UAV0   = AV MO WIND SPEED(m/s)

              !Check the format of the wp1 first and read monthly values. Jaehak 2019
              READ(KR(K),'(A10)')TMPSTR
              Read (TMPSTR,'(F10.0)',iostat=io) TMPFL
              IF( io/= 0)THEN
                  FMT='(12F6.0)'
              ELSE
                  FMT='(12F10.0)'
              ENDIF
              BACKSPACE(KR(K))
              
              READ(KR(K),FMT)(OBMX(IWI,I),I=1,12)
              READ(KR(K),FMT)(OBMN(IWI,I),I=1,12)
              READ(KR(K),FMT)(SDTMX(IWI,I),I=1,12)
              READ(KR(K),FMT)(SDTMN(IWI,I),I=1,12)
              READ(KR(K),FMT)(RMO(IWI,I),I=1,12)
              READ(KR(K),FMT)(RST(2,IWI,I),I=1,12)
              READ(KR(K),FMT)(RST(3,IWI,I),I=1,12)
              READ(KR(K),FMT)(PRW(1,IWI,I),I=1,12)
              READ(KR(K),FMT)(PRW(2,IWI,I),I=1,12)
              READ(KR(K),FMT)(UAVM(I),I=1,12)
              READ(KR(K),FMT)(WI(IWI,I),I=1,12)
              READ(KR(K),FMT)(OBSL(IWI,I),I=1,12)
              READ(KR(K),FMT)(RH(IWI,I),I=1,12)
              READ(KR(K),FMT)(UAV0(I),I=1,12)
              REWIND KR(K)
              K=K+1
              IF(IWI==1)THEN
	              TAV(1)=.25*(OBMX(IWI,12)+OBMN(IWI,12)+OBMX(IWI,1)+OBMN(IWI,1))
	              BAV=TAV(1)
                  SVA=BAV   
                  JT1=1
                  IF(OBMN(IWI,1)>OBMN(IWI,12))JT1=12
                  TMN(1)=OBMN(IWI,JT1)
                  DO I=2,12
                      I1=I-1
	                  TAV(I)=.25*(OBMX(IWI,I)+OBMN(IWI,I)+OBMX(IWI,I1)+OBMN(IWI,I1))
	                  IF(TAV(I)>BAV)THEN
                          BAV=TAV(I)                   
                      ELSE
                          IF(TAV(I)<SVA)SVA=TAV(I)
                      END IF
                      IF(OBMN(IWI,I)>TMN(1))CYCLE
                      JT1=I
                      TMN(1)=OBMN(IWI,I)
                  END DO
              END IF
              RHCF=1.
              AMPX=BAV-SVA
              DO I=1,12
	              IF(RST(2,IWI,I)<1.E-5.OR.RST(3,IWI,I)<1.E-5)EXIT
	          END DO
	          IF(I>12)THEN
                  ICDP=0
              ELSE
                  ICDP=1
                  SUM=0.
                  DO I=1,10000
                      XX=AUNIF(IDG(3))
                      SUM=SUM+(-LOG(XX))**EXPK
                  END DO
                  REXP=10100./SUM
              END IF
              ISA=1
              DO I=1,12
                  I1=I+1
                  XM=NC(I1)-NC(I)
                  JDA=(NC(I1)+NC(I))*.5
                  CALL WHLRMX(JDA)
                  SRMX(I,1)=RAMX
                  THRL(I,1)=HRLT
                  XYP(I)=0.
                  XX=SDTMX(IWI,I)-SDTMN(IWI,I)
                  IF(XX>10.)THEN
                      SDTMX(IWI,I)=(SDTMX(IWI,I)-OBMX(IWI,I))*.25
                      SDTMN(IWI,I)=(OBMN(IWI,I)-SDTMN(IWI,I))*.25
                  END IF
                  IF(PRW(1,IWI,I)>0.)THEN
                      UAVM(I)=XM*PRW(1,IWI,I)/(1.-PRW(2,IWI,I)+PRW(1,IWI,I))
                  ELSE
                      PRW(1,IWI,I)=BTA*(UAVM(I)+.0001)/XM
                      PRW(2,IWI,I)=1.-BTA+PRW(1,IWI,I)
                  END IF
                  RST(1,IWI,I)=RMO(IWI,I)/(UAVM(I)+.01)
                  IF(OBSL(IWI,I)<=0.)THEN                                                         
	                  X1=MAX(.8,.21*SQRT(OBMX(IWI,I)-OBMN(IWI,I)))                                                                                    
	                  OBSL(IWI,I)=X1*RAMX                                                                
	              END IF 
                  V3=AUNIF(IDG(3))
                  IF(ICDP>0)THEN
                      RST(1,IWI,I)=RST(1,IWI,I)*REXP
                      PCF(IWI,I)=1.
                      CYCLE
                  END IF
                  SUM=0.
                  R6=RST(3,IWI,I)/6.
                  DO J=1,1000
                      V4=AUNIF(IDG(3))
                      XX=ADSTN(V3,V4)
                      V3=V4
                      R1=WRAIN(R6,XX,RST,IWI,I)
                      SUM=SUM+R1
                  END DO
                  PCF(IWI,I)=1010.*RST(1,IWI,I)/SUM
              END DO          
              XYP(1)=OBMX(IWI,1)
              BIG(1)=OBSL(IWI,1)
              UPLM=RH(IWI,1)
              BLM=RH(IWI,1)
              RFMX=RMO(IWI,1)
              EXTM=WI(IWI,1)
              RHCF=1. 
              DO I=2,12
                  TX=.5*(OBMX(IWI,I)+OBMN(IWI,I))
                  RTO=RH(IWI,I)/TX
                  IF(RTO>1.)RHCF=100.
                  IF(OBSL(IWI,I)>BIG(1))BIG(1)=OBSL(IWI,I)
                  IF(RMO(IWI,I)>RFMX)RFMX=RMO(IWI,I)
                  IF(RH(IWI,I)>UPLM)THEN
                      UPLM=RH(IWI,I)
                  ELSE
                      IF(RH(IWI,I)<BLM)BLM=RH(IWI,I)
                  END IF
                  IF(WI(IWI,I)>EXTM)EXTM=WI(IWI,I)
                  XYP(1)=XYP(1)+OBMX(IWI,I)
              END DO
              RUNT=1.
              IF(BIG(1)>100.)RUNT=.04184
              XYP(1)=XYP(1)/12.
              X3=.3725/(XYP(1)+20.)
              DO I=1,12
                  XM=NC(I+1)-NC(I)
                  WFT(IWI,I)=UAVM(I)/XM
                  XYP(2)=XYP(2)+OBMN(IWI,I)
                  XYP(3)=XYP(3)+RMO(IWI,I)
                  XYP(4)=XYP(4)+UAVM(I)
                  OBSL(IWI,I)=RUNT*OBSL(IWI,I)
                  XYP(5)=XYP(5)+OBSL(IWI,I)
                  X1=MAX(RMO(IWI,I),12.7)
                  TX=.5*(OBMX(IWI,I)+OBMN(IWI,I))
                  IF(UPLM>1.)THEN
                      IF(BLM<0..OR.RHCF<1.1)THEN
                          RH(IWI,I)=ASVP(RH(IWI,I)+273.)/ASVP(TX+273.)
                      ELSE
                          RH(IWI,I)=RH(IWI,I)/RHCF
                      END IF
                  ELSE
                      IF(RH(IWI,I)<1.E-10)THEN
                          XX=OBMX(IWI,I)-OBMN(IWI,I)
                          RH(IWI,I)=.9-.8*XX/(XX+EXP(5.122-.1269*XX))
                      END IF
                  END IF              
                  X2=MAX(TX,-1.7)
                  XYP(6)=XYP(6)+((X1/25.4)/(1.8*X2+22.))**1.111
                  X1=RMO(IWI,I)/(UAVM(I)+1.E-10)
                  IF(EXTM>1.)THEN
                      F=XY2/(UAVM(I)+.01)
                      XTP(I)=WI(IWI,I)
                      WI(IWI,I)=-XTP(I)/LOG(F)
                      WI(IWI,I)=APM*WI(IWI,I)/(X1+1.)
                      IF(WI(IWI,I)<.1)WI(IWI,I)=.1
                      IF(WI(IWI,I)>.95)WI(IWI,I)=.95
                  ELSE
                      IF(EXTM<1.E-10)WI(IWI,I)=APM*X3*(OBMX(IWI,I)+20.)
                      XTP(I)=5.3*X1*WI(IWI,I)
                  END IF
              END DO          
              XYP(2)=XYP(2)/12.
              XYP(5)=XYP(5)/12.
              XYP(6)=115.*XYP(6)
              IF(KFL(1)>0)THEN
                  WRITE(KW(1),3050)
                  WRITE(KW(1),602)TITWP(IWI),FWXP(IWI)
                  WRITE(KW(1),3300)HED(1),(OBMX(IWI,I),I=1,12),XYP(1),HED(1)
                  WRITE(KW(1),3300)HED(2),(OBMN(IWI,I),I=1,12),XYP(2),HED(2)
                  WRITE(KW(1),3010)'SDMX',(SDTMX(IWI,I),I=1,12),'SDMX'
                  WRITE(KW(1),3010)'SDMN',(SDTMN(IWI,I),I=1,12),'SDMN'
                  WRITE(KW(1),2550)HED(4),(RMO(IWI,I),I=1,12),XYP(3),HED(4)
                  WRITE(KW(1),3440)'SDRF',(RST(2,IWI,I),I=1,12),'SDRF'
                  WRITE(KW(1),3010)'SKRF',(RST(3,IWI,I),I=1,12),'SKRF'
                  WRITE(KW(1),3460)'PW/D',(PRW(1,IWI,I),I=1,12),'PW/D'
                  WRITE(KW(1),3460)'PW/W',(PRW(2,IWI,I),I=1,12),'PW/W'
                  WRITE(KW(1),3300)'DAYP',(UAVM(I),I=1,12),XYP(4),'DAYP'
                  WRITE(KW(1),3440)'P5MX',(XTP(I),I=1,12),'P5MX'
                  WRITE(KW(1),2550)HED(3),(OBSL(IWI,I),I=1,12),XYP(5),HED(3)
                  WRITE(KW(1),3440)'RAMX',(SRMX(I,1),I=1,12),'RAMX'
                  WRITE(KW(1),3010)'HRLT',(THRL(I,1),I=1,12),'HRLT'
                  WRITE(KW(1),3010)'RHUM',(RH(IWI,I),I=1,12),'RHUM'
                  WRITE(KW(1),3010)'ALPH',(WI(IWI,I),I=1,12),'ALPH'
                  WRITE(KW(1),3010)' PCF',(PCF(IWI,I),I=1,12),' PCF'
              END IF    
          END DO
          KND=K-1
          IWI=NWP
          IF(NGN<=0)THEN
              IF(KFL(1)>0)THEN
                  WRITE(KW(1),'(/T10,A)')'__________RAIN, TEMP, RAD, WIND &
                  SPEED, & REL HUM ARE GENERATED__________'
                  IF(NGN<0)THEN
                      WRITE(KW(1),'(T10,A)')'RAINFALL IS SAME FOR ALL SUBAREAS'
                  ELSE
                      WRITE(KW(1),'(T10,A)')'RAINFALL IS SPATIALLY DISTRIBUTED'
                  END IF
              END IF    
          ELSE
              L=1
              N1=NGN
              DO J=4,1,-1
                  CALL AISPL(N1,N2)
                  IF(N1==0)EXIT
                  KGN(N1)=1
                  SELECT CASE(N1)
	                  CASE(1)
	                      N1=N2
	                      CYCLE
	                  CASE(2)
	                      KDT2(L)=59
	                  CASE(3)
	                      KDT2(L)=3
	                  CASE(4)
	                      KDT2(L)=7
	                  CASE(5)
	                      KDT2(L)=8
	                  CASE DEFAULT
	              END SELECT
                  L=L+1
                  N1=N2
              END DO
              IF(KFL(1)>0)THEN
                  SELECT CASE(L)
                      CASE(1)
                          WRITE(KW(1),'(/T10,A)')'__________RAIN IS INPUT__________'
	                  CASE(2)
	                      WRITE(KW(1),'(/T10,A,1X,A4,A)')'__________RAIN,',HED&
	                      (KDT2(1)),' ARE INPUT__________'     
	                  CASE(3)
	                      WRITE(KW(1),3650)(HED(KDT2(J)),J=1,2)
	                  CASE(4)
	                      WRITE(KW(1),3660)(HED(KDT2(J)),J=1,3)
	                  CASE(5)
	                      WRITE(KW(1),3670)(HED(KDT2(J)),J=1,4)
                  END SELECT
              END IF        
	      END IF
	      IF(INFL==4)THEN
	          CALL OPENV(KR(29),FRFDT,IDIR)
	          I2=-1
              DO WHILE(I2/=IRFT) 
	              READ(KR(29),*,IOSTAT=NFL)I2,RFTFILE
                  IF(NFL/=0)THEN
                      WRITE(*,*)'RAINDT NO = ',IRFT,' NOT IN RFTFILE LIST FILE'
                      STOP
                  END IF
              END DO
              REWIND KR(29)
              CALL OPENV(KR(30),RFTFILE,IDIR)
              IF(KFL(1)>0)WRITE(KW(1),'(/T10,A,A80)')'WITHIN STORM RAINFALL READ FROM ',RFTFILE
	      END IF
          ! ICMO=COMMAND NUMBERS OF WATERSHED OUTLETS
          ! LINE 9  
          READ(KR(1),811)ICMO
          DO I=1,10
              IF(ICMO(I)==0)EXIT
          END DO
          NCMO=I-1	      
          IF(KFL(1)>0)WRITE(KW(1),'(/T10,A,F6.3)')'WET-DRY PROB COEF = ',BTA
          XTP(2)=XTP(2)/12.
          AAP=XYP(3)
          RFAD=AAP/365.25
          AVT=(XYP(1)+XYP(2))*.5
	      CLF=RFAD/AVT
          X1=MAX(1.,AVT)
	      CLF0=SQRT(AAP/X1)
          IF(UXP<1.E-10)UXP=.5
          IF(DIAM<1.E-10)DIAM=500.
          USTRT=.0161*SQRT(DIAM)
          USTT=USTRT*USTRT
          ! UAVM = AV MO WIND SPEED(m/s)(REQUIRED TO SIMULATE WIND
          ! EROSION--ACW>0 LINE 23  AND POTENTIAL ET IF PENMAN OR
          ! PENMAN-MONTEITH EQS ARE USED--LINE 4)
          ! LINE 7
          READ(KR(22),3191)UAVM
          AWV=0.
          WB=0.
          XTP=0.
          DO I=1,12
              IF(UAVM(I)<1.E-5)UAVM(I)=UAV0(I)
              AWV=AWV+UAVM(I)
              DO J=1,100
                  RN2=AUNIF(IDG(5))
                  WV=UAVM(I)*(-LOG(RN2))**UXP
                  IF(WV>6.)THEN
                      EV=193.*EXP(1.103*(WV-30.)/(WV+1.))
                      XTP(I)=XTP(I)+EV
                  END IF
              END DO
              WB=WB+XTP(I)
          END DO
          AWV=AWV/12.
          WCF=(3.86*AWV**3/XYP(6)**2)**.3484
          IF(KFL(1)>0)WRITE(KW(1),3300)HED(7),UAVM,AWV,HED(7)
          IF(IET==4.AND.PRMT(34)>0.)THEN
              HGX=PRMT(34)
          ELSE
              IF(WCF<.5)THEN
                  HGXW=.5
              ELSE
                  HGXW=MIN(.6,WCF)
              END IF    
              IF(CLF0>7.)THEN
                  HGX0=.5
              ELSE
                  HGX0=MIN(.6,1.2-.1*CLF0)
              END IF 
              HGX=.5*(HGXW+HGX0)       
          END IF                        
          ! DIR  = AV MO FRACTION OF WIND FROM 16 DIRECTIONS(BLANK UNLESS
          ! WIND EROSION IS SIMULATED--ACW>0 LINE 23).
          DO J=1,16
              ! LINES 25/40
              READ(KR(22),3191)(DIR(I,J),I=1,12)
              IF(DIR(1,J)>0.)CYCLE
              DO I=1,12
                  DIR(I,J)=1.
              END DO
          END DO        
          REWIND KR(22)
          IF(KFL(1)>0)THEN
              WRITE(KW(1),'(/T10,A)')'WIND DIRECTION DISTRIBUTION'
              WRITE(KW(1),'(T20,A,A80)')'WIND = ',TITWN
              WRITE(KW(1),3280)NBYR
          END IF    
          XTP=0.
          DO I=1,12
              XTP(1)=XTP(1)+DIR(I,1)
              IF(UAVM(1)>0.)THEN
                  CALL AEXINT(UXP,SUM)
                  UAVM(I)=UAVM(I)/SUM
              END IF
              DO J=2,16
                  XTP(J)=XTP(J)+DIR(I,J)
                  DIR(I,J)=DIR(I,J)+DIR(I,J-1)
              END DO
              DO J=1,16
                  DIR(I,J)=DIR(I,J)/DIR(I,16)
              END DO
          END DO        
          SUM=0.
          DO J=1,16
              SUM=SUM+XTP(J)
          END DO
          DO J=1,16
              XTP(J)=XTP(J)/SUM
              IF(KFL(1)>0)WRITE(KW(1),3300)HEDW(J),(DIR(I,J),I=1,12),XTP(J),HEDW(J)
          END DO
          IF(KFL(1)>0)THEN
              WRITE(KW(1),'(T10,A,F7.3)')'WIND EROS CLIMATIC FACTOR = ',WCF                                                                              
              WRITE(KW(1),'(T10,A,F7.2)')'CLIMATIC FACTOR = ',CLF0
              SELECT CASE(IET)
                  CASE(1)
                      WRITE(KW(1),'(/T10,A)')'___________PENMAN-MONTEITH  EQ &
                      USED TO EST POT ET___________'
                  CASE(2)
                      WRITE(KW(1),'(/T10,A)')'___________PENMAN EQ USED TO EST &
                      POT ET___________'
                  CASE(3)
                      WRITE(KW(1),'(/T10,A)')'___________PRIESTLEY-TAYLOR EQ &
                      USED TO EST POT ET___________'
                  CASE(4)
                      WRITE(KW(1),'(/T10,A,F8.4,A,F5.2)')'___________HARGREAVES EQ USED TO &
                      EST POT ET ___________ PARM(23)=',PRMT(23),' PARM(34)=',HGX
                  CASE(5)
                      WRITE(KW(1),'(/T10,A)')'__________ BAIER-ROBERTSON EQ &
                      USED TO EST POT ET __________'
                  CASE DEFAULT
                      WRITE(KW(1),'(/T10,A,F8.4,A,F5.2)')'___________HARGREAVES EQ USED TO &
                      EST POT ET ___________ PARM(23)=',PRMT(23),' PARM(34)=',HGX
              END SELECT
          END IF        
          CALL ALPYR(IYR,NYD,LPYR)
	      CALL ADAJ(NC,IBD,IMO,IDA,NYD)
	      JDA=IBD
          MO=1
          CALL AXMON
          MO1=MO
          TX=(OBMX(1,MO)+OBMN(1,MO))/2.
          ST0=OBSL(1,MO)
          DST0=TX
          DO I=1,100
              KDC1(I)=0
              NX(I)=I
          END DO
     2140 LC=0
          NDT=0
          NDP=0
          NDF=0
          ISV=0
          X1=IWTB*RFAD
          DO I=1,MSA
              TITSO(I)=''
              TITOP(I)=''
              SMRF(I)=X1
              SMEO(I)=X1
              DO J=1,MSL
                  LID(J,I)=J
                  LORG(J,I)=J
              END DO
              LID(ML1,I)=ML1
          END DO
          CALL AINLZ
	      CALL AINIX
	      DO J=1,NRN0
              IDG(J)=J
          END DO
          IYX=IYR0-1880
          ! RANDOM NUMBER GENERATOR ID NUMBERS
          ! IDG = 1 DETERMINES WET AND DRY DAYS
          !     = 2 RELATES WEATHER VARIABLES TO RAIN
          !     = 3 RAINFALL AMOUNT
          !     = 4 RAINFALL ENERGY(EI)- PEAK RUNOFF RATE(QRB)
          !     = 5 WIND SPEED
          !     = 6 WIND DIRECTION
          !     = 7 RELATIVE HUMIDITY
          !     = 8 RUNOFF CURVE NUMBER
          !     = 9 WITHIN DAY WIND SPEED DIST
          !     =10 X COORDINATE OF SPATIAL RAINFALL GENERATOR
          !     =11 Y COORDINATE OF SPATIAL RAINFALL GENERATOR
          !     =12 MULTIPLE WEATHER GEN PARM SELECTOR
          !     =13 PADDY SEDIMENT CONC
          IF(IGN>0)THEN
              DO KK=1,IGN
                  DO J=1,NRN0
                      XX=AUNIF(NRNG)
                      IX(J)=IX(NRNG)
                  END DO
              END DO        
              CALL AISHFL
          END IF
          DO J=1,NRN0
              IX0(J)=IX(J)
          END DO
          IF(KFL(1)>0)THEN
              WRITE(KW(1),3070)IGN,(IX(IDG(I)),I=1,NRN0),(IDG(I),I=1,NRN0)
	          WRITE(KW(1),714)
          END IF    
  	      V3=AUNIF(IDG(3))
          V1=AUNIF(IDG(2))
          IF(IHRD>0)THEN 
	          IOW1=1
	          DO 
                  IHD=IHD+1
                  ! INITIAL HERD DATA
                  ! IDON = OWNER ID#
                  ! NCOW = NUMBER OF COWS PER HERD BY OWNER(HD)
                  ! IDMU = MANURE ID # (FERTCOM.DAT)
                  ! FFED = MINIMUM FRACTION OF DAY HERD IN FEEDING AREA.
                  ! GZRT = GRAZING RATE FOR EACH HERD (kg/hd/d)
                  ! DUMP = DAILY UNIT MANURE PROD(kg/hd/d)
                  ! VURN = VOLUME OF URINE (l/hd/d)
                  READ(KR(21),2320,IOSTAT=NFL)IOW,(XTP(J),J=1,6)
                  IF(NFL/=0)EXIT
	              IF(IOW==0.OR.IOW>MOW)EXIT
	              IF(IOW1/=IOW)THEN
	                  IOW1=IOW
                      IHD=1
	              END IF
                  NHRD(IOW)=NHRD(IOW)+1
	              NCOW(IHD,IOW)=XTP(1)
	              IDMU(IHD,IOW)=XTP(2)
	              FFED(IHD,IOW)=XTP(3)
	              GZRT(IHD,IOW)=XTP(4)
	              DUMP(IHD,IOW)=XTP(5)
	              VURN(IHD,IOW)=XTP(6)
                  NMC=NMC+NCOW(IHD,IOW)
                  IF(KFL(1)>0)WRITE(KW(1),646)IOW,IHD,NCOW(IHD,IOW),IDMU&
                  (IHD,IOW),(XTP(J),J=3,6)
              END DO
          END IF
          IF(IHRD<2)THEN    
              IHDM=1
              NHRD=1          
          END IF
          IF(KFL(1)>0)WRITE(KW(1),101)NMC
          IHBS=0
	      NYHO=0
          !IF(NMC>0)THEN
	          !IF(KFL(1)>0)WRITE(KW(1),713)
  	          !IYHO=1 
              ! HERD MANAGEMENT DATA
              ! READ OPERATION SCHEDULE
              ! 1    IOW  = OWNER NUMBER
              ! 2    IHD  = HERD NUMBER
              ! 3    I1   = YR OF BUY/SELL
              ! 4    I2   = MO OF BUY/SELL
              ! 5    I3   = DAY OF BUY/SELL
              ! 6    I4   = NUMBER OF ANIMALS IN HERD AFTER BUY/SELL
              !DO
                  !READ(KR(21),3100)IOW,IHD,I1,I2,I3,I4
	              !IF(IOW==0)EXIT
	              !IHBS(IHD,IOW)=IHBS(IHD,IOW)+1
	              !NHBS(IHD,IOW)=IHBS(IHD,IOW)
	              !NBSX(IHBS(IHD,IOW),IHD,IOW)=I4
       	          !IF(I1>IYHO(IHD,IOW))IYHO(IHD,IOW)=I1
	              !NYHO(IHD,IOW)=IYHO(IHD,IOW)
	              !IHDT(IHBS(IHD,IOW),IHD,IOW)=I1*10000+I2*100+I3
	              !IF(KFL(1)>0)WRITE(KW(1),848)I1,I2,I3,IOW,IHD,I4
	          !END DO
	      !END IF
  	      IYHO=1
          REWIND KR(21)
          IF(KFL(36)>0)THEN
              WRITE(KW(36),3381)IYER,IMON,IDAY,IT1,IT2,IT3
              WRITE(KW(36),'(T10,A,I4,1X,A,1X,A80)')'RUN # ',IRUN,'NAME',ASTN
          END IF
          K1=1
          SUMA=0.
          TOT=0.
	      NFED=0
          XCU=0.
          YCU=0.
          XCS=1.E10
          YCS=1.E10
          ELMN=1.E10
          XZP=0.
          TWMB=0.
	      NBON=0
	      TBTN=0.
	      NCP=1
        
        
        !rtb hsg
        !read in temporally (daily) variable hydrologic soil group (HSG) for each soil type, if 
        !the file (soil_hsg_variable) is present in the model folder
        inquire (file='soil_hsg_variable',exist=hsg_vary)
        if(hsg_vary) then
          open(8855,file='soil_hsg_variable')
          read(8855,*)
          read(8855,*) num_soil_type
          allocate(hsg_values(num_soil_type,NBYR,366))
          allocate(cnsc1_vary(msa,NBYR,366))
          allocate(cnsc2_vary(msa,NBYR,366))
          hsg_values = 0.
          cnsc1_vary = 0.
          cnsc2_vary = 0.
          read(8855,*)
          year_index = IYR0
          do yrc=1,NBYR
            !number of days in the year
            if(mod(year_index,4) == 0) then
              year_days = 366
            else
              year_days = 365
            endif
            do dayc=1,year_days
              read(8855,*) day,(hsg_values(j,yrc,dayc),j=1,num_soil_type)
            enddo
            year_index = year_index + 1
          enddo
        endif
        
        
          DO ISA=1,MSA  !SUBAREA LOOP
              KRST(ISA)=ISA+1000
              IBSA(ISA)=ISA
              JSA(ISA)=ISA
              XX=0.
              SUM=0.
              BIG(ISA)=.15
              RZ(ISA)=3.
              PAW(ISA)=0.
              RZSW(ISA)=0.
              PMX(ISA)=PRMT(43)
              ! NBSA=SUBAREA 
              ! LINE 1
              READ(KR(5),'(I8)')NBSA(ISA)
              IF(NBSA(ISA)==0)EXIT
              IF(NBSA(ISA)>NBMX0)THEN
                  NBMX0=NBSA(ISA)
                  DEALLOCATE(NISA,FPSO)
                  ALLOCATE(NISA(NBMX),FPSO(NBMX))
              END IF
              ! SUBAREA DATA
              !  1  INPS = SOIL # FROM TABLE KR(13)
              !  2  IOPS = OP SCHED # FROM TABLE KR(15)
              !  3  IOW  = OWNER ID #
              !  4  II   = 0 FOR NON FEEDING AREAS
              !          = HERD # FOR FEEDING AREAS
              !  5  IAPL = 0 FOR NON MANURE APPL AREAS
              !          = - FEED AREA ID # FOR LIQUID MANURE APPL AREAS
              !          =   FEED AREA ID # FOR SOLID MANURE APPL AREAS
              !  6  NOT USED
              !  7  NVCN = 0 VARIABLE DAILY CN NONLINEAR CN/SW WITH DEPTH SOIL WATER
              !              WEIGHTING
              !          = 1 VARIABLE DAILY CN NONLINEAR CN/SW NO DEPTH WEIGHTING
              !          = 2 VARIABLE DAILY CN LINEAR CN/SW NO DEPTH WEIGHTING
              !          = 3 NON-VARYING CN--CN2 USED FOR ALL STORMS
              !          = 4 VARIABLE DAILY CN SMI(SOIL MOISTURE INDEX)
              !  8  IWTH = INPUT DAILY WEATHER STATION NUMBER
              !  9  IPTS = POINT SOURCE NUMBER
              ! 10  ISAO = 0 FOR NORMAL RES PRINCIPAL SPILLWAY RELEASE
              !          = ID OF SUBAREA RECEIVING OUTFLOW FROM BURRIED PIPE OUTLET
              ! 11  LUNS = LAND USE NUMBER FROM NRCS LAND USE-HYDROLOGIC SOIL GROUP
              !            TABLE.  OVERRIDES LUN IN .OPC FILE.  
              ! 12  IMW  = MIN INTERVAL BETWEEN AUTO MOW
              !     LINE 2
              READ(KR(5),*)INPS,IOPS,IOW,II,IAPL(ISA),I1,NVCN(ISA),IWTH(ISA),&
              IPTS(ISA),ISAO(ISA),LUNS(ISA),IMW(ISA)
              IF(IMW(ISA)==0)IMW(ISA)=IMW0 
              ! INITIAL CONDITIONS
              !  1  SNO  = WATER CONTENT OF SNOW COVER(mm)
              !  2  STDO = STANDING DEAD CROP RESIDUE(t/ha)
              !  3  YCT  = LATITUDE OF SUBAREA CENTROID
              !  4  XCT  = LONGITUDE OF SUBAREA CENTROID
              !  5  AZM  = AZIMUTH ORIENTATION OF LAND SLOPE (DEGREES CLOCKWISE FROM NORTH)
              !  6  SAEL = SUBAREA ELEVATION(m)
              !  7  FL   = FIELD LENGTH(km)(0 IF UNKNOWN)
              !  8  FW   = FIELD WIDTH(km)(0 IF UNKNOWN)
              !  9  ANGL = CLOCKWISE ANGLE OF FIELD LENGTH FROM NORTH(deg)(0 IF
              !            UNKNOWN)
              !     LINE 3
              READ(KR(5),*)SNO(ISA),STDO(ISA),YCT(ISA),XCT(ISA),AZM,SAEL,FL,FW,ANGL
              XCU=MAX(ABS(XCT(ISA)),XCU)
              XCS=MIN(ABS(XCT(ISA)),XCS)
              YCU=MAX(ABS(YCT(ISA)),YCU)
              YCS=MIN(ABS(YCT(ISA)),YCS)
              IF(SAEL<ELMN)THEN
                  ELMN=SAEL
                  KSA=ISA
              END IF    
              IF(FL<1.E-10)FL=FL0
              IF(FW<1.E-10)FW=FW0
              IF(ANGL<1.E-10)ANGL=ANG0
              IF(FL<1.E-10)FL=.632
              IF(FW<1.E-10)FW=.316
              ANG=ANGL/CLT
              IF(ABS(XCT(ISA))>0..OR.ABS(YCT(ISA))>0.)THEN
                  YLAZ=YCT(ISA)
              ELSE
                  YLAZ=YLAT
              END IF
              IF(IOW>NBON)NBON=IOW
              IDON(ISA)=IOW
              NSAO(IOW)=NSAO(IOW)+1
              IDOW(NSAO(IOW),IOW)=ISA
              IDFH(ISA)=II
              IF(II>0)THEN
                  IFED(II,IOW)=ISA
                  NFED(IOW)=NFED(IOW)+1
                  IDFD(NFED(IOW),IOW)=NBSA(ISA)
                  IDFA(NFED(IOW),IOW)=ISA
              END IF 
              IF(NVCN(ISA)==0)NVCN(ISA)=NVCN0
              WCHS=0.
              ! CATCHMENT CHARACTERISTICS
              !  1  WSA  = DRAINAGE AREA(ha)
              !  2  CHL  = CHANNEL LENGTH(km)(0 IF UNKNOWN)
              !  3  CHD  = CHANNEL DEPTH(m)(0 IF UNKNOWN)
              !  4  CHS  = CHANNEL SLOPE(m/m)(0 IF UNKNOWN)
              !  5  CHN  = MANNINGS N FOR CHANNEL(0 IF UNKNOWN)
              !  6  STP  = AVE UPLAND SLOPE(m/m)
              !  7  SPLG = AVE UPLAND SLOPE LENGTH(m)
              !  8  UPN  = MANNINGS N FOR UPLAND(0 IF UNKNOWN)
              !  9  FFPQ = FRACTION FLOODPLAIN FLOW--PARTITIONS FLOW THRU FILTER STRIPS
              ! 10  URBF = URBAN FRACTION OF SUBAREA
              !     LINE 4
              READ(KR(5),*)WSA(ISA),CHL(ISA),CHD,CHS(ISA),CHN(ISA),STP(ISA),&
              SPLG(ISA),UPN,FFPQ(ISA),URBF(ISA)
              IF(IAZM>0)THEN
                  X1=ASIN(STP(ISA))
                  YLAZ=YLAZ/CLT
                  AZR=AZM/CLT
                  YLAZ=CLT*ASIN(STP(ISA)*COS(AZR)*COS(YLAZ)+COS(X1)*SIN(YLAZ))
                  IF(KFL(1)>0)THEN
                      WRITE(KW(1),'(T10,A,F8.3)')'AZIMUTH ORIENTATION OF LAND SLOPE',&
                      '(DEGREES CLOCKWISE FROM NORTH)',AZM
                      WRITE(KW(1),'(T10,A,F8.3)')'EQUIVALENT LATITUDE = ',YLAZ
                  END IF    
              END IF
              XX=YLAZ/CLT
              YLS(ISA)=SIN(XX)
              YLC(ISA)=COS(XX)
              YTN(ISA)=TAN(XX)
              PB=101.3-ELEV*(.01152-5.44E-7*ELEV)
              GMA(ISA)=6.595E-4*PB
              IRF(ISA)=ISA
              IF(NGN>0)CALL WDLYSTA	
              AHSM=CAHU(1,365,0.,1)
              Z1=ABS(WSA(ISA))
              GRDL(ISA)=SQRT(Z1*1.E4)
              TOT=TOT+Z1
              X3=1.+Z1
              IF(STP(ISA)<1.E-10)STP(ISA)=.001
              ISCH=1
              IF(WSA(ISA)<0.)ISCH=0
              IF(CHD<1.E-10)CHD=COCH(3)*Z1**COCH(4)
              IF(CHL(ISA)<1.E-10)CHL(ISA)=COCH(5)*Z1**COCH(6)
              IF(CHN(ISA)<1.E-10)CHN(ISA)=.05
              IF(UPN<1.E-10)UPN=.15
              X8=SQRT(STP(ISA))
              X1=100.*STP(ISA)
              IF(SCLM(18)>0.)X1=MIN(X1,SCLM(18))
              CNSX(ISA)=1.1-.3*X1/(X1+EXP(SCRP(18,1)-SCRP(18,2)*X1))
              IF(PRMT(92)>0.)THEN
                  SCNX(ISA)=PRMT(92)
              ELSE
                  IF(SCLM(28)>0.)X1=MIN(X1,SCLM(28))
                  SCNX(ISA)=MIN(1.7,2.-X1/(X1+EXP(SCRP(28,1)-SCRP(28,2)*X1)))
                  SCNX(ISA)=MAX(1.,SCNX(ISA))
              END IF    
              IF(CHS(ISA)<1.E-10)CHS(ISA)=STP(ISA)*X3**(-.3)
              UPSX(ISA)=SPLG(ISA)/22.127
              XM=.3*STP(ISA)/(STP(ISA)+EXP(-1.47-61.09*STP(ISA)))+.2
              SLF(ISA)=UPSX(ISA)**XM*(STP(ISA)*(65.41*STP(ISA)+4.56)+.065)
              X1=3.*STP(ISA)**.8+.56
              BETA=STP(ISA)/(.0896*X1)
              RXM=BETA/(1.+BETA)
              RLF(ISA)=UPSX(ISA)**RXM
              IF(SPLG(ISA)>4.57)THEN
                  IF(STP(ISA)>.09)THEN
                      RSF(ISA)=16.8*STP(ISA)-.5 
                  ELSE
                      RSF(ISA)=10.8*STP(ISA)+.03
                  END IF
              ELSE    
                  RSF(ISA)=X1
              END IF
              ACET(1,ISA)=Z1
              ACET(2,ISA)=CHL(ISA)
              ACET(3,ISA)=CHD
              ACET(4,ISA)=CHS(ISA)
              ACET(5,ISA)=CHN(ISA)
              ACET(6,ISA)=STP(ISA)
              ACET(7,ISA)=SPLG(ISA)
              ACET(8,ISA)=UPN
              X1=MIN(SPLG(ISA),SQRT(10000.*Z1))
              IF(ITYP>0)THEN
                  IF(CHL(ISA)>.1)THEN
                      SFL=50.
                  ELSE
                      IF(CHL(ISA)>.05)THEN
                          SFL=100.*(CHL(ISA)-.05)
                      ELSE
                          SFL=0.
                      END IF
                  END IF
                  TSF=SFL/MIN(2160.,17712.*X8)
                  X1=MAX(CHL(ISA)-(SPLG(ISA)+SFL)*.001,0.)
                  TCC(ISA)=X1/(3.6*CHD**.66667*SQRT(CHS(ISA))/CHN(ISA))
                  TCS(ISA)=.0913*(SPLG(ISA)*UPN)**.8/STP(ISA)**.4
                  TCC(ISA)=TCC(ISA)+TSF
                  TC(ISA)=TCC(ISA)+TCS(ISA)
              ELSE
                  TCS(ISA)=.0216*(SPLG(ISA)*UPN)**.75/STP(ISA)**.375
                  TCC(ISA)=1.75*CHL(ISA)*CHN(ISA)**.75/(Z1**.125*CHS(ISA)**.375)
                  X4=MIN(SPLG(ISA)/360.,TCS(ISA))
                  TC(ISA)=X4+TCC(ISA)
              END IF
              ! CHANNEL GEOMETRY OF ROUTING REACH THRU SUBAREA
              !  1  RCHL = CHANNEL LENGTH OF ROUTING REACH(km)
              !  2  RCHD = CHANNEL DEPTH(m)(0 IF UNKNOWN)
              !  3  RCBW = BOTTOM WIDTH OF CHANNEL(m)(0 IF UNKNOWN)
              !  4  RCTW = TOP WIDTH OF CHANNEL(m)(0 IF UNKNOWN)
              !  5  RCHS = CHANNEL SLOPE(m/m)(0 IF UNKNOWN)
              !  6  RCHN = MANNINGS N VALUE OF CHANNEL(0 IF UNKNOWN)
              !  7  RCHC = USLE C FOR CHANNEL
              !  8  RCHK = USLE K FOR CHANNEL
              !  9  RFPW = FLOODPLAIN WIDTH(m)(0 IF UNKNOWN)
              ! 10  RFPL = FLOODPLAIN LENGTH(km)(0 IF UNKNOWN)
              ! 11  SAT1 = SATURARTED CONDUCTIVITY(GREEN & AMPT) ADJUSTMENT FACTOR(.01_10.)
              ! 12  FPS1 = FLOODPLAIN SATURARTED CONDUCTIVITY ADJUSTMENT FACTOR(.0001_10.)
              SAT1=0.
              FPS1=0.
              ! LINE 5
              READ(KR(5),*)RCHL(ISA),RCHD(ISA),RCBW(ISA),RCTW(ISA),RCHS(ISA),&
              RCHN(ISA),RCHC(ISA),RCHK(ISA),RFPW(ISA),RFPL(ISA),SAT1,FPS1
              IF(RCHC(ISA)<1.E-10)RCHC(ISA)=RCC0
              ! RESERVOIR DATA
              !  1  RSEE = ELEV AT EMERGENCY SPILLWAY ELEV(m)
              !  2  RSAE = SURFACE AREA AT EMERGENCY SPILLWAY ELEV(ha)
              !  3  RVE0 = VOLUME AT EMERGENCY SPILLWAY ELEV(mm)
              !  4  RSEP = ELEV AT PRINCIPAL SPILLWAY ELEV(m)
              !  5  RSAP = SURFACE AREA AT PRINCIPAL SPILLWAY ELEV(ha)
              !  6  RVP0 = VOLUME AT PRINCIPAL SPILLWAY ELEV(mm)
              !  7  RSV  = INITIAL VOLUME(mm)
              !  8  RSRR = TIME TO RELEASE FLOOD STORAGE(d)
              !  9  RSYS = INITIAL SEDIMENT CONCENTRATION(ppm)
              ! 10  RSYN = NORMAL SEDIMENT CONC(ppm)
              ! 11  RSHC = BOTTOM HYDRAULIC CONDUCTIVITY(mm/h)
              ! 12  RSDP = TIME REQUIRED TO RETURN TO NORMAL SED CONC AFTER RUNOFF
              !            EVENT(d)
              ! 13  RSBD = BULK DENSITY OF SEDIMENT IN RESERVOIR(t/m^3)
              ! 14  PCOF = FRACTION OF SUBAREA CONTROLLED BY PONDS
              ! 15  BCOF = FRACTION OF SUBAREA CONTROLLED BY BUFFERS
              ! 16  BFFL = BUFFER FLOW LENGTH (m)
              ! 17  WTMN = MIN DEPTH TO WATER TABLE(m)(0 IF UNKNOWN)
              ! 18  WTMX = MAX DEPTH TO WATER TABLE(m)(0 IF UNKNOWN)
              ! 19  WTBL = INITIAL WATER TABLE HEIGHT(m)(0 IF UNKNOWN)
              ! 20  GWST = GROUNDWATER STORAGE (mm)
              ! 21  GWMX = MAXIMUM GROUNDWATER STORAGE (mm)
              ! 22  RFTT = GROUNDWATER RESIDENCE TIME(d)(0 IF UNKNOWN)
              ! 23  RFPK = RETURN FLOW / (RETURN FLOW + DEEP PERCOLATION)
              ! LINE 6/7
              READ(KR(5),*)RSEE(ISA),RSAE(ISA),RVE0(ISA),RSEP(ISA),&
              RSAP(ISA),RVP0(ISA),RSV(ISA),RSRR(ISA),RSYS(ISA),RSYN(ISA),&
              RSHC(ISA),RSDP(ISA),RSBD(ISA),PCOF(ISA),BCOF(ISA),&
              BFFL(ISA),WTMN(ISA),WTMX(ISA),WTBL(ISA),GWST(ISA),&
              GWMX(ISA),RFTT(ISA),RFPK(ISA)
              SVX=0.
              IF(RCHL(ISA)<1.E-10)RCHL(ISA)=CHL(ISA)
              IF(RFPL(ISA)<1.E-10)RFPL(ISA)=.95*RCHL(ISA)
              IF(RSBD(ISA)<1.E-10)RSBD(ISA)=.8
              IF(ABS(RCHL(ISA)-CHL(ISA))<1.E-10)THEN
                  ! SUBAREA	
	              IEXT(ISA)=1
                  ISV=ISV+1
                  SAV(ISV)=SUMA
                  SUMA=Z1
                  ICDT(K1)=1
                  IDOT(K1)=K1
                  IDN1T(K1)=ISA
                  IDN2T(K1)=0
                  IDOS(ISV)=K1-1
                  TC(K1)=TC(ISA)
              ELSE
                  ! ROUTE      
                  SUMA=SUMA+Z1
                  ICDT(K1)=2
                  IDOT(K1)=K1
                  IDN1T(K1)=IDOT(K1-1)
                  IDN2T(K1)=ISA
                  TC(K1)=TC(IDN1T(K1))
                  K1=K1+1
                  ! SUBAREA      
                  ICDT(K1)=1
                  IDOT(K1)=K1
                  IDN1T(K1)=ISA
                  IDN2T(K1)=0
                  TC(K1)=TC(ISA)
                  K1=K1+1
                  ! ADD      
                  ICDT(K1)=3
                  IDOT(K1)=K1
                  IDN1T(K1)=IDOT(K1-1)
                  IDN2T(K1)=IDOT(K1-2)
                  TC(K1)=TC(IDN1T(K1))+TC(IDN2T(K1))
              END IF
              K1=K1+1
              IF(PCOF(ISA)<1.E-20)PCOF(ISA)=PCO0 
              X1=ABS(RSAE(ISA))
              ! RESERVOIR
              IF(X1>0.)THEN
                  ICDT(K1)=4
                  IDOT(K1)=K1
                  IDN1T(K1)=IDOT(K1-1)
	              IDN2T(K1)=ISA
                  K1=K1+1
                  ! jjeong wrote
				  IF(RSEE(ISA)==-1)THEN !Special arrangement for reading daily discharge rates from .dat
                    write (filename, "(I3,A4)") NBSA(ISA), ".DAT"
                    OPEN(ISA+100001,FILE=filename)
                    READ(ISA+100001,*) YY0,MM0,DD0
                    DO 
                        IF(YY0==IYR0.AND.MM0==IMO.AND.DD0==IDA0) THEN
                            BACKSPACE(ISA+100001)
                            EXIT
                        ENDIF
                        READ(ISA+100001,*) YY0,MM0,DD0
                    END DO
                    
                        
                !ELSEIF(RSEE(ISA)==-2)THEN !Read monthly average discharge from .dat
                !ELSEIF(RSEE(ISA)==-3)THEN !Read yearly average discharge from .dat
                !ELSEIF(RSEE(ISA)==-4)THEN !Read average annual discahrge from .dat
                ENDIF  
                
              ELSE
                  IF(PCOF(ISA)>0.)THEN
                      ! GENERIC PONDS
                      ICDT(K1)=5
                      IDOT(K1)=K1
                      IDN1T(K1)=IDOT(K1-1)
	                  IDN2T(K1)=ISA
                      K1=K1+1
	                  RVP0(ISA)=50.
                      RVE0(ISA)=55.
	                  RSV(ISA)=RVP0(ISA)
	                  XX=1000.*Z1*PCOF(ISA)
	                  IF(XX<260.)THEN
	                      RSAP(ISA)=3.E-5*XX
                      ELSE
	                      XX=XX-256.
	                      TW=.5*(16.+SQRT(XX))
	                      RSAP(ISA)=.0001*TW*TW
	                  END IF
                      RSAE(ISA)=1.01*RSAP(ISA)
	                  RSYS(ISA)=300.
	                  RSYN(ISA)=100.
	                  RSHC(ISA)=.00001
	                  RSDP(ISA)=5.
	                  RSBD(ISA)=.8
	                  RSRR(ISA)=20.
                  END IF    
              END IF    
              IF(ISCH==0)THEN
                  ! ADD      
                  WSA(ISA)=Z1
                  SUMA=SUMA+SAV(ISV)
                  ICDT(K1)=3
                  IDOT(K1)=K1
                  IDN1T(K1)=IDOT(K1-1)
                  IDN2T(K1)=IDOS(ISV)
                  TC(K1)=TC(IDN1T(K1))+TC(IDN2T(K1))
                  SVX=SAV(ISV)
                  K1=K1+1
                  ISV=ISV-1
              END IF    
              RWSA(ISA)=SUMA-SVX
              X3=RWSA(ISA)
              X1=X3
              IF(RCHN(ISA)<1.E-10)RCHN(ISA)=.05
              IF(RCHL(ISA)<1.E-10)RCHL(ISA)=.001*GRDL(ISA)
              GRDL(ISA)=1.414*GRDL(ISA)
              X3=X3+1.
              IF(RCHS(ISA)<1.E-10)RCHS(ISA)=STP(ISA)*X3**(-.3)
              WCHS=WCHS+RCHS(ISA)*Z1
              RFPS(ISA)=RCHS(ISA)*RCHL(ISA)/RFPL(ISA)
              IF(QG>0.)THEN
                  QX=.002778*X1*QG0*(1./RWSA(ISA))**QCF
                  X1=BWD+2.
                  X1=(QX*RCHN(ISA)/(X1*(X1/(BWD+4.472))**.6667*SQRT(RCHS(ISA))))**&
                  .375
                  IF(RCBW(ISA)<1.E-10)RCBW(ISA)=BWD*X1
                  IF(RCTW(ISA)<1.E-10)RCTW(ISA)=RCBW(ISA)+4.*X1
                  IF(RCHD(ISA)<1.E-10)RCHD(ISA)=X1
                  IF(RFPW(ISA)<1.E-10)RFPW(ISA)=FCW*RCTW(ISA)
              ELSE
                  IF(RCTW(ISA)<1.E-10)RCTW(ISA)=COCH(1)*X3**COCH(2)
                  IF(RFPW(ISA)<1.E-10)RFPW(ISA)=10.*RCTW(ISA)
                  IF(RCHD(ISA)<1.E-10)RCHD(ISA)=COCH(3)*X3**COCH(4)
                  IF(RCBW(ISA)<1.E-10)RCBW(ISA)=MAX(RCTW(ISA)-4.*RCHD(ISA),.1*&
                  RCTW(ISA))
              END IF
              IF(RCBW(ISA)>RCTW(ISA))RCBW(ISA)=.75*RCTW(ISA)
              RFPX(ISA)=SQRT(RFPS(ISA))*RFPW(ISA)/UPN
              BFSN(ISA)=X8/UPN
              ZCH=RCHD(ISA)
              CBW=RCBW(ISA)
              RCSS(ISA)=.5*(RCTW(ISA)-CBW)/ZCH
              X3=SQRT(RCSS(ISA)*RCSS(ISA)+1.)
              RCHX(ISA)=SQRT(RCHS(ISA))/RCHN(ISA)
              CHXA(ISA)=ZCH*(CBW+ZCH*RCSS(ISA))
              CHXP(ISA)=CBW+2.*ZCH*X3
              QCAP(ISA)=CHXA(ISA)**1.66667*RCHX(ISA)/CHXP(ISA)**.66667
              WCHS=WCHS/TOT
              IF(IHY==3.AND.IEXT(ISA)==0)THEN
                  !XL=.001*SQRT(3.E4*RWSA(ISA))
                  !TCX=.0663*XL**.77/RCHS(ISA)**.385
                  !QPX=.9*RWSA(ISA)/TCX
                  !QPX=13.*SQRT(RWSA(ISA))
                  QPX=3591.
                  SSS=SQRT(RCSS(ISA)*RCSS(ISA)+1.)
                  !XCTW=RCTW(ISA)
                  !XFPW=RFPW(ISA)
                  CALL HQDAV(AO2,CBW,QPX,SSS,ZCH,ZI2,RCTW(ISA),RFPW(ISA),ISA)
                  !DZRT=ZI2/NPRC
                  DZRT=.3333333
                  CALL RATECVM_CM(CBW,DZRT,SSS,ZCH)
              END IF
              IF(KFL(28)>0.OR.KFL(5)>0)THEN
                  IF(BCHL>0..AND.BCHS>0.)THEN
                      WCHL=.3*SQRT(TOT/3.)
                      TCB=.0663*BCHL**.77/BCHS**.385
                      TCW=.0663*WCHL**.77/WCHS**.385
                      DRSW=MIN(.95,(TCW/TCB)**PRMT(37))
                      BS2=LOG10(DRSW)/2.699
                      BS1=1./.1**BS2
                  ELSE
                      BS1=0.
                      BS2=1.
                      DRSW=1.
                  END IF
              END IF
              IF(KFL(1)>0)WRITE(KW(1),'(T10,A,F8.4)')'DELIVERY RATIO TO SWAT = ',DRSW
              IF(ISOL==0)THEN
                  I2=-1
                  DO WHILE(I2/=INPS)
                      READ(KR(13),*,IOSTAT=NFL)I2,SOILFILE
                      IF(NFL/=0)THEN
                          WRITE(*,*)'SOIL NO = ',INPS,' NOT IN SOIL LIST FILE &
                          SAID = ',NBSA(ISA)
                          STOP
                      END IF
                  END DO
                  CALL OPENV(KR(14),SOILFILE,IDIR)
                  REWIND KR(13)
              ELSE
                  WRITE(ASOL,'(I8.8)')NBSA(ISA)
                  SOILFILE=ASOL//".SOT"
                  CALL OPENV(KR(14),SOILFILE,0)
              END IF
              TITSO(ISA)=SOILFILE
              ! LINE 1
              READ(KR(14),'()')
              ! SOIL PROPERTIES
              !  1  SALB = SOIL ALBEDO
              !  2  HSG  = HYDROLOGIC SOIL GROUP--1.=A; 2.=B; 3.=C; 4.=D
              !  3  FFC  = FRACTION OF FIELD CAP FOR INITAL WATER STORAGE(BLANK IF
              !            UNKNOWN)
              ! 4-10 DUM
              ! 11  TSLA = MAXIMUM NUMBER OF SOIL LAYERS(3-30)
              ! 12  XIDS = 0. FOR CALCAREOUS SOILS AND NON CALCAREOUS
              !                  WITHOUT WEATHERING INFORMATION
              !          = 1. FOR NON CACO3 SLIGHTLY WEATHERED
              !          = 2. NON CACO3 MODERATELY WEATHERED
              !          = 3. NON CACO3 HIGHLY WEATHERED
              !          = 4. INPUT PSP
              
                ! 13  RTN1 = NUMBER YEARS OF CULTIVATION AT START OF SIMULATION. BLANK                                                  
              !            DEFAULTS TO RTN0.                                                                              
              ! 14  XIDK = 1 FOR KAOLINITIC SOIL GROUP                                                                    
              !          = 2 FOR MIXED SOIL GROUP                                                                         
              !          = 3 FOR SMECTITIC SOIL GROUP                                                                     
              ! 15  ZQT  = MINIMUM THICKNESS OF MAXIMUM LAYER(m)(SPLITTING                                                
              !            STOPS WHEN ZQT IS REACHED)                                                                     
              ! 16  ZF   = MINIMUM PROFILE THICKNESS(m)--STOPS SIMULATION.                                                
              ! 17  ZTK  = MINIMUM LAYER THICKNESS FOR BEGINNING SIMULATION LAYER                                         
              !            SPLITTING--MODEL SPLITS FIRST LAYER WITH THICKNESS GREATER                                     
              !            THAN ZTK(m); IF NONE EXIST THE THICKEST LAYER IS SPLIT.                                        
              ! 18  FBM  = FRACTION OF ORG C IN BIOMASS POOL(.03-.05)                                                     
              ! 19  FHP  = FRACTION OF HUMUS IN PASSIVE POOL(.3-.7)                                                       
              ! 20  XCC  = CODE WRITTEN AUTOMATICALLY BY .SOT (NOT USER INPUT)                                            
              !     LINE 2/3                                                                                              
              READ(KR(14),3120)SALB(ISA),HSG,FFC(ISA),DUM,DUM,DUM,DUM,&
              DUM,DUM,DUM,TSLA(ISA),XIDS(ISA),RTN1,XIDK(ISA),ZQT,ZF,ZTK,&
              FBM(ISA),FHP(ISA),XCC                                                             
              NCC=XCC                                                                                                     
              IF(HSG<1.E-10.OR.HSG>4.)THEN                                                                                
                  WRITE(*,*)'HYDROLOGIC SOIL GROUP NO MISSING SOIL = ',SOILFILE                                           
                  STOP                                                                                                    
              END IF                                                                                                      
              GWSN(ISA)=0.                                                                                                
              IF(WTMX(ISA)<1.E-10)THEN                                                                                    
                  WTMN(ISA)=50.                                                                                           
                  WTMX(ISA)=100.                                                                                          
                  WTBL(ISA)=75.                                                                                           
              END IF                                                                                                      
              IDS(ISA)=XIDS(ISA)+1.1                                                                                      
              IF(RTN1<1.E-5)RTN1=RTN0                                                                                     
              IF(FBM(ISA)<1.E-10)FBM(ISA)=.04                                                                             
              IF(FHP(ISA)<1.E-10)FHP(ISA)=.7-.3*EXP(-.0277*RTN1)                                                          
              IF(GWST(ISA)<1.E-10)GWST(ISA)=GWS0                                                                          
              IF(GWMX(ISA)<1.E-10)GWMX(ISA)=GWS0                                                                          
              IF(RFPK(ISA)<1.E-10)RFPK(ISA)=RFP0                                                                          
              IF(FFC(ISA)<1.E-10)FFC(ISA)=AAP/(AAP+EXP(9.043-.002135*AAP))                                                
              IDSK=MAX(XIDK(ISA),1.)                                                                                      
              MXLA=TSLA(ISA)                                                                                              
              IF(ZQT<1.E-10)ZQT=.1                                                                                        
              IF(ZF<1.E-10)ZF=.1                                                                                          
              IF(ZTK<1.E-10)ZTK=.15                                                                                       
              IF(RFTT(ISA)<1.E-10)RFTT(ISA)=RFT0                                                                          
              ZNH3(ISA)=0.                                                                                                
              ZFOP(ISA)=0.                                                                                                
              ZLS(ISA)=0.                                                                                                 
              ZLM(ISA)=0.                                                                                                 
              ZLSL(ISA)=0.                                                                                                
              ZLSC(ISA)=0.                                                                                                
              ZLMC(ISA)=0.                                                                                                
              ZLSLC(ISA)=0.                                                                                               
              ZLSLNC(ISA)=0.                                                                                              
              ZBMC(ISA)=0.                                                                                                
              ZHSC(ISA)=0.                                                                                                
              ZHPC(ISA)=0.                                                                                                
              ZLSN(ISA)=0.                                                                                                
              ZLMN(ISA)=0.                                                                                                
              ZBMN(ISA)=0.                                                                                                
              ZHSN(ISA)=0.                                                                                                
              ZHPN(ISA)=0.                                                                                                
              !     READ SOIL DATA                                                                                        
              !     THE SOIL IS DIVIDED VERTICALLY INTO LAYERS(MAX OF 30 LAYERS                                           
              !     OF USER SPECIFIED THICKNESS).  DATA ARE INPUT 10 LAYERS AT A TIME.                                    
              !     THUS 10 VALUES OF THE FOLLOWING DATA ARE INPUT ON SPECIFIED LINES.                                    
              !     LINES 9/29 CONTAIN DATA FOR LAYERS 1/10                                                               
              !  4  Z    = DEPTH TO BOTTOM OF LAYERS(m)                                                                   
              !  5  BD   = BULK DENSITY(t/m3)                                                                             
              !  6  UW   = SOIL WATER CONTENT AT WILTING POINT(1500 KPA)(m/m)                                             
              !            (BLANK IF UNKNOWN)                                                                             
              !  7  FC   = WATER CONTENT AT FIELD CAPACITY(33KPA)(m/m)                                                    
              !            (BLANK IF UNKNOWN)                                                                             
              !  8  SAN  = % SAND                                                                                         
              !  9  SIL  = % SILT                                                                                         
              ! 10  WN   = INITIAL ORGANIC N CONC(g/t)       (BLANK IF UNKNOWN)                                           
              ! 11  PH   = SOIL PH                                                                                        
              ! 12  SMB  = SUM OF BASES(cmol/kg)              (BLANK IF UNKNOWN)                                          
              ! 13  WOC  = ORGANIC CARBON CONC(%)                                                                         
              ! 14  CAC  = CALCIUM CARBONATE(%)                                                                           
              ! 15  CEC  = CATION EXCHANGE CAPACITY(cmol/kg)(BLANK IF UNKNOWN)                                            
              ! 16  ROK  = COARSE FRAGMENTS(% VOL)              (BLANK IF UNKNOWN)                                        
              ! 17  CNDS = INITIAL SOL N CONC(g/t)            (BLANK IF UNKNOWN)                                          
              ! 18  SSF  = INITIAL SOL P CONC(g/t)       (BLANK IF UNKNOWN)                                               
              ! 19  RSD  = CROP RESIDUE(t/ha)                (BLANK IF UNKNOWN)                                           
              ! 20  BDD  = BULK DENSITY(OVEN DRY)(t/m3)   (BLANK IF UNKNOWN)                                              
              ! 21  PSP  = P SORPTION RATIO                   (BLANK IF UNKNOWN)                                          
              ! 22  SATC = SATURATED CONDUCTIVITY(mm/h)     (BLANK IF UNKNOWN)                                            
              ! 23  HCL  = LATERAL HYDRAULIC CONDUCTIVITY(mm/h)                                                           
              ! 24  WPO  = INITIAL ORGANIC P CONC(g/t)      (BLANK IF UNKNOWN)                                            
              ! 25  DHN  = EXCHANGEABLE K CONC (g/t)                                                                      
              ! 26  ECND = ELECTRICAL COND (mmho/cm)                                                                      
              ! 27  STFR = FRACTION OF STORAGE INTERACTING WITH NO3 LEACHING                                              
              !                                               (BLANK IF UNKNOWN)                                          
              ! 28  SWST = INITIAL SOIL WATER STORAGE (m/m)                                                               
              ! 29  CPRV = FRACTION INFLOW PARTITIONED TO VERTICLE CRACK OR PIPE FLOW                                     
              ! 30  CPRH = FRACTION INFLOW PARTITIONED TO HORIZONTAL CRACK OR PIPE                                        
              !            FLOW                                                                                           
              ! 31  WLS  = STRUCTURAL LITTER(kg/ha)           (BLANK IF UNKNOWN)                                          
              ! 32  WLM  = METABOLIC LITTER(kg/ha)            (BLANK IF UNKNOWN)                                          
              ! 33  WLSL = LIGNIN CONTENT OF STRUCTURAL LITTER(kg/ha)(B I U)                                              
              ! 34  WLSC = CARBON CONTENT OF STRUCTURAL LITTER(kg/ha)(B I U)                                              
              ! 35  WLMC = C CONTENT OF METABOLIC LITTER(kg/ha)(B I U)                                                    
              ! 36  WLSLC= C CONTENT OF LIGNIN OF STRUCTURAL LITTER(kg/ha)(B I U)                                         
              ! 37  WLSLNC=N CONTENT OF LIGNIN OF STRUCTURAL LITTER(kg/ha)(BIU)                                           
              ! 38  WBMC = C CONTENT OF BIOMASS(kg/ha)(BIU)                                                               
              ! 39  WHSC = C CONTENT OF SLOW HUMUS(kg/ha)(BIU)                                                            
              ! 40  WHPC = C CONTENT OF PASSIVE HUMUS(kg/ha)(BIU)                                                         
              ! 41  WLSN = N CONTENT OF STRUCTURAL LITTER(kg/ha)(BIU)                                                     
              ! 42  WLMN = N CONTENT OF METABOLIC LITTER(kg/ha)(BIU)                                                      
              ! 43  WBMN = N CONTENT OF BIOMASS(kg/ha)(BIU)                                                               
              ! 44  WHSN = N CONTENT OF SLOW HUMUS(kg/ha)(BIU)                                                            
              ! 45  WHPN = N CONTENT OF PASSIVE HUMUS(kg/ha)(BIU)                                                         
              ! 46  FE26 = IRON CONTENT(%)                                                                                
              ! 47  SULF = SULFUR CONTENT(%)                                                                              
              ! 48  ASHZ = SOIL HORIZON(A,B,C)                                                                            
              ! 49  CGO2 = O2 CONC IN GAS PHASE (g/m3 OF SOIL AIR)                                                        
              ! 50  CGCO2= CO2 CONC IN GAS PHASE (g/m3 OF SOIL AIR)                                                       
              ! 51  CGN2O= N2O CONC IN GAS PHASE (g/m3 OF SOIL AIR)                                                       
              !     LINES 4/47
              READ(KR(14),3111)(Z(I,ISA),I=1,MSL)                                                                         
              READ(KR(14),3111)(BD(I,ISA),I=1,MSL)                                                                        
              READ(KR(14),3111)(UW(I),I=1,MSL)                                                                            
              READ(KR(14),3111)(FC(I,ISA),I=1,MSL)                                                                        
              READ(KR(14),3111)(SAN(I,ISA),I=1,MSL)                                                                       
              READ(KR(14),3111)(SIL(I,ISA),I=1,MSL)                                                                       
              READ(KR(14),3111)(WON(I,ISA),I=1,MSL)                                                                       
              READ(KR(14),3111)(PH(I,ISA),I=1,MSL)                                                                        
              READ(KR(14),3111)(SMB(I,ISA),I=1,MSL)                                                                       
              READ(KR(14),3111)(WOC(I,ISA),I=1,MSL)                                                                       
              READ(KR(14),3111)(CAC(I,ISA),I=1,MSL)                                                                       
              READ(KR(14),3111)(CEC(I,ISA),I=1,MSL)                                                                       
              READ(KR(14),3111)(ROK(I,ISA),I=1,MSL)                                                                       
              READ(KR(14),3111)(CNDS(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(SSF(I,ISA),I=1,MSL)                                                                       
              READ(KR(14),3111)(RSD(I,ISA),I=1,MSL)                                                                       
              READ(KR(14),3111)(BDD(I,ISA),I=1,MSL)                                                                       
              READ(KR(14),3111)(PSP(I,ISA),I=1,MSL)                                                                       
              READ(KR(14),3111)(SATC(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(HCL(I,ISA),I=1,MSL)                                                                       
              READ(KR(14),3111)(WPO(I,ISA),I=1,MSL)                                                                       
              READ(KR(14),3111)(DHN(I,ISA),I=1,MSL)                                                                       
              READ(KR(14),3111)(ECND(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(STFR(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(SWST(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(CPRV(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(CPRH(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(WLS(I,ISA),I=1,MSL)                                                                       
              READ(KR(14),3111)(WLM(I,ISA),I=1,MSL)                                                                       
              READ(KR(14),3111)(WLSL(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(WLSC(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(WLMC(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(WLSLC(I,ISA),I=1,MSL)                                                                     
              READ(KR(14),3111)(WLSLNC(I,ISA),I=1,MSL)                                                                    
              READ(KR(14),3111)(WBMC(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(WHSC(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(WHPC(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(WLSN(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(WLMN(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(WBMN(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(WHSN(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(WHPN(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(FE26(I,ISA),I=1,MSL)                                                                      
              READ(KR(14),3111)(SULF(I,ISA),I=1,MSL)                                                                      
              ! LINE 48                                                                                                   
              READ(KR(14),'(15A8)')ASHZ                                                                                   
              ! LINES 49/51                                                                                               
              IF(IDNT>2)THEN                                                                                              
                  READ(KR(14),3111)(CGO2(I,ISA),I=1,MSL)                                                                  
                  READ(KR(14),3111)(CGCO2(I,ISA),I=1,MSL)                                                                 
                  READ(KR(14),3111)(CGN2O(I,ISA),I=1,MSL)                                                                 
              ELSE                                                                                                        
                  CGO2=0.                                                                                                 
                  CGCO2=0.                                                                                                
                  CGN2O=0.                                                                                                
              END IF                                                                                                      
              L=1                                                                                                         
              LZ=0                                                                                                        
              XX=0.                                                                                                       
              XCB=.2                                                                                                      
              CLOSE(KR(14))                                                                                               
              DO I=1,MSL !SOIL LOOP                                                                                             
                  IF(Z(I,ISA)<1.E-10)EXIT                                                                           
                  IF(SAN(I,ISA)<1.E-5)THEN                                                                                
                      WRITE(KW(36),*)'SUBAREA ',ISA,' SOIL LAYER ',I,' SAN=0.0'                                           
                      IF(I==1)THEN                                                                                        
                          SAN(I,ISA)=.33                                                                                  
                      ELSE                                                                                                
                          SAN(I,ISA)=SAN(I-1,ISA)                                                                         
                      END IF                                                                                              
                  END IF                                                                                                  
                  IF(SIL(I,ISA)<1.E-5)THEN                                                                                
                      WRITE(KW(36),*)'SUBAREA ',ISA,' SOIL LAYER ',I,' SIL=0.0'                                           
                      IF(I==1)THEN                                                                                        
                          SIL(I,ISA)=.33                                                                                  
                      ELSE                                                                                                
                          SIL(I,ISA)=SIL(I-1,ISA)                                                                         
                      END IF                                                                                              
                  END IF                                                                                                  
                  CLA(I,ISA)=100.-SAN(I,ISA)-SIL(I,ISA)
                  !Calculate median particle size for the top soil layer, paddy modeling Jaehak Jeong 2014
                  IF(I==1) THEN
                      SUB_D50(ISA)=EXP(.0041*CLA(I,ISA)+.0271*SIL(I,ISA)+.057*SAN(I,ISA)) 
                  END IF
                  IF(BD(I,ISA)<1.E-5)THEN                                                                                 
                      WRITE(KW(36),*)'SUBAREA ',ISA,' SOIL LAYER ',I,' BD=0.0'                                            
                      BD(I,ISA)=1.25+.005*SAN(I,ISA)                                                                      
                  END IF                                                                                                  
                  X1=BD(I,ISA)                                                                                            
                  DG=1000.*(Z(I,ISA)-XX)                                                                                  
                  CALL SBDSC(BD(I,ISA),PRMT(2),F,I,1)                                                                     
                  CALL SDST(RSD,DG,DG1,1.,.01,I,ISA,MSL)                                                                 
                  CALL SDST(SSF,DG,DG,20.,.001,I,ISA,MSL)                                                                 
                  CALL SDST(CNDS,DG,DG,10.,.001,I,ISA,MSL)                                                                
                  CALL SDST(DHN,DG,DG,100.,.001,I,ISA,MSL)                                                                
                  TRSD(ISA)=TRSD(ISA)+RSD(I,ISA)                                                                          
                  ZD=.25*(XX+Z(I,ISA))                                                                                    
                  F=ZD/(ZD+EXP(-.8669-2.0775*ZD))                                                                         
                  STMP(I,ISA)=F*(AVT-TX)+TX                                                                               
                  IF(WOC(I,ISA)<=0.)WOC(I,ISA)=XCB*EXP(-.001*DG)                                                          
                  XCB=WOC(I,ISA)                                                                                          
                  XZ=WOC(I,ISA)*.0172                                                                                     
                  ZZ=1.-XZ                                                                                                
                  WT(I,ISA)=BD(I,ISA)*DG*10.                                                                              
                  DG1=DG                                                                                                  
                  WT1=WT(I,ISA)/1000.                                                                                     
                  X1=10.*WOC(I,ISA)*WT(I,ISA)                                                                             
                  WOC(I,ISA)=X1                                                                                           
                  IF(WON(I,ISA)>0.)THEN                                                                                   
                      WON(I,ISA)=WT1*WON(I,ISA)                                                                           
                      KK=0                                                                                                
                  ELSE                                                                                                    
                      WON(I,ISA)=.1*WOC(I,ISA)                                                                            
                      KK=1                                                                                                
                  END IF                                                                                                  
                  IF(NCC==0)THEN                                                                                          
                      WBM=FBM(ISA)*X1                                                                                     
                      WBMC(I,ISA)=WBM                                                                                     
                      IF(KK>0)THEN                                                                                        
                          RTO=.1                                                                                          
                      ELSE                                                                                                
                          RTO=WON(I,ISA)/WOC(I,ISA)                                                                       
                      END IF                                                                                              
                      WBMN(I,ISA)=RTO*WBMC(I,ISA)                                                                         
                      WHP=FHP(ISA)*(X1-WBM)                                                                               
                      WHS=X1-WBM-WHP                                                                                      
                      WHSC(I,ISA)=WHS                                                                                     
                      WHSN(I,ISA)=RTO*WHSC(I,ISA)                                                                         
                      WHPC(I,ISA)=WHP                                                                                     
                      WHPN(I,ISA)=RTO*WHPC(I,ISA)                                                                         
                      X1=RSD(I,ISA)                                                                                       
                      IF(I==1)X1=X1+STDO(ISA)                                                                             
                      WLM(I,ISA)=500.*X1                                                                                  
                      WLS(I,ISA)=WLM(I,ISA)                                                                               
                      WLSL(I,ISA)=.8*WLS(I,ISA)                                                                           
                      WLMC(I,ISA)=.42*WLM(I,ISA)                                                                          
                      WLMN(I,ISA)=.1*WLMC(I,ISA)                                                                          
                      WLSC(I,ISA)=.42*WLS(I,ISA)                                                                          
                      WLSLC(I,ISA)=.8*WLSC(I,ISA)                                                                         
                      WLSLNC(I,ISA)=.2*WLSC(I,ISA)                                                                        
                      WLSN(I,ISA)=WLSC(I,ISA)/150.                                                                        
                      WOC(I,ISA)=WOC(I,ISA)+WLSC(I,ISA)+WLMC(I,ISA)                                                       
                      WON(I,ISA)=WON(I,ISA)+WLSN(I,ISA)+WLMN(I,ISA)                                                       
                  END IF                                                                                                  
                  FOP(I,ISA)=RSD(I,ISA)*1.1                                                                               
                  WPMA(I,ISA)=0.                                                                                          
                  IF(WPO(I,ISA)>0.)THEN                                                                                   
                      WPO(I,ISA)=WT1*WPO(I,ISA)                                                                           
                  ELSE                                                                                                    
                      WPO(I,ISA)=.125*WON(I,ISA)                                                                          
                  END IF                                                                                                  
                  PO(I,ISA)=1.-BD(I,ISA)/2.65                                                                             
                  X2=.1*WOC(I,ISA)/WT(I,ISA)                                                                              
                  SELECT CASE(ISW)                                                                                      
                  CASE(0,2)                                                                                               
                      CALL SWRTNR(CLA(I,ISA),SAN(I,ISA),X2,UW(I),FC(I,ISA))                                               
                  CASE(1,3)                                                                                               
                      IF(UW(I)<1.E-10.OR.FC(I,ISA)<1.E-10)CALL SWRTN_BNW(CLA(I,ISA),&                                       
                      SAN(I,ISA),X2,BD(I,ISA),UW(I),FC(I,ISA))                                                            
                  CASE(4,5)                                                                                               
                      CALL SWNN(CLA(I,ISA),SAN(I,ISA),X2,UW(I),FC(I,ISA))                                                 
                  CASE(6,7)                                                                                               
                      CALL SWRTN_BNW(CLA(I,ISA),SAN(I,ISA),X2,BD(I,ISA),UW(I),FC(I,ISA))                                    
                  CASE DEFAULT                                                                                            
                      CALL SWRTN_BNW(CLA(I,ISA),SAN(I,ISA),X2,BD(I,ISA),UW(I),FC(I,ISA))                                    
                      ISW=6                                                                                               
                  END SELECT 
                  IF(ROK(I,ISA)>99.)ROK(I,ISA)=90.                                                                        
                  XY=1.-ROK(I,ISA)*.01                                                                                    
                  UW(I)=UW(I)*XY                                                                                          
                  XY=XY*DG                                                                                                
                  FC(I,ISA)=FC(I,ISA)*XY                                                                                  
                  S15(I,ISA)=UW(I)*DG                                                                                     
                  PO(I,ISA)=PO(I,ISA)*XY                                                                                  
                  CALL SPOFC(I)
                  VGA(I,ISA)=EXP(-4.3003-.0097*CLA(I,ISA)+.0138*SAN(I,ISA)-.1706*X2)
                  VGA(I,ISA)=max(VGA(I,ISA),0.001)
                  VGN0=1.
                  VGA0=VGA(I,ISA)
                  UW0=S15(I,ISA)
                  PO0=PO(I,ISA)
                  FC0=FC(I,ISA)
                  DO IT=1,10
                      FU=UW0+(PO0-UW0)/(1.+(VGA0*336.27)**VGN0)**(1.-1./VGN0)-FC0
                      IF(ABS(FU)<1.E-5)EXIT
                      IF(IT==1)THEN
                          DF=.01
                      ELSE
                          if(abs(VGN0-VG1)>1e-5)then
                             DFDN=(FU-FU1)/(VGN0-VG1)
                             if(abs(DFDN)>1e-5)then
                                DF=FU/DFDN
                             else
                                DF=.01
                             endif
                             
                          else
                             DF=.01
                          endif
                      END IF
                      VG1=VGN0
                      FU1=FU
                      VGN0=VGN0-DF
                      IF(VGN0<1.2)VGN0=1.2
                      IF(ABS(DF)<1.E-5)EXIT
                  END DO    
                  VGN(I,ISA)=VGN0
                  IF(NCC==0)THEN
                      IF(SWST(I,ISA)<1.E-10)SWST(I,ISA)=FFC(ISA)
                      SWST(I,ISA)=SWST(I,ISA)*(FC(I,ISA)-S15(I,ISA))+S15(I,ISA)
                  END IF    
                  SEV(1,ISA)=SEV(1,ISA)+XY                                                                                
                  SEV(3,ISA)=SEV(3,ISA)+WT(I,ISA)                                                                         
                  XX=Z(I,ISA)                                                                                             
                  IF(HCL(I,ISA)<1.E-10)HCL(I,ISA)=STP(ISA)*SATC(I,ISA)                                                    
                  IF(CEC(I,ISA)>0.)THEN                                                                                   
                      IF(CAC(I,ISA)>0.)SMB(I,ISA)=CEC(I,ISA)                                                              
                      IF(SMB(I,ISA)>CEC(I,ISA))SMB(I,ISA)=CEC(I,ISA)                                                      
                      BSA=SMB(I,ISA)*100./(CEC(I,ISA)+1.E-20)                                                             
                      IF(PH(I,ISA)>5.6)THEN                                                                               
                          ALS(I,ISA)=0.                                                                                   
                      ELSE                                                                                                
                          X1=.1*WOC(I,ISA)/WT(I,ISA)                                                                      
                          ALS(I,ISA)=154.2-1.017*BSA-3.173*X1-14.23*PH(I,ISA)                                             
                          IF(ALS(I,ISA)<0.)THEN                                                                           
                              ALS(I,ISA)=0.                                                                               
                          ELSE                                                                                            
                              IF(ALS(I,ISA)>95.)ALS(I,ISA)=95.                                                            
                          END IF                                                                                          
                      END IF                                                                                              
                  ELSE                                                                                                    
                      CEC(I,ISA)=PH(I,ISA)                                                                                
                      SMB(I,ISA)=CEC(I,ISA)                                                                               
                      ALS(I,ISA)=0.                                                                                       
                  END IF                                                                                                  
                  SELECT CASE(IDS(ISA))                                                                                   
                  CASE(1)                                                                                                 
                      IF(CAC(I,ISA)>0.)THEN                                                                               
                          PSP(I,ISA)=.58-.0061*CAC(I,ISA)                                                                 
                          BPT(I,ISA)=7.6E-4                                                                               
                      ELSE                                                                                                
                          PSP(I,ISA)=.5                                                                                   
                          BPT(I,ISA)=EXP(-1.77*PSP(I,ISA)-7.05)                                                           
                      END IF                                                                                              
                  CASE(2)                                                                                                 
                      PSP(I,ISA)=.02+.0104*SSF(I,ISA)                                                                     
                      BPT(I,ISA)=EXP(-1.77*PSP(I,ISA)-7.05)                                                               
                  CASE(3)                                                                                                 
                      PSP(I,ISA)=.0054*BSA+.116*PH(I,ISA)-.73                                                             
                      BPT(I,ISA)=EXP(-1.77*PSP(I,ISA)-7.05)                                                               
                  CASE(4)                                                                                                 
                      PSP(I,ISA)=.46-.0916*LOG(CLA(I,ISA))                                                                
                      BPT(I,ISA)=EXP(-1.77*PSP(I,ISA)-7.05)                                                               
                  CASE(5)                                                                                                 
                      IF(CAC(I,ISA)>0.)THEN                                                                               
                          BPT(I,ISA)=7.6E-4                                                                               
                      ELSE                                                                                                
                          BPT(I,ISA)=EXP(-1.77*PSP(I,ISA)-7.05)                                                           
                      END IF                                                                                              
                  END SELECT                                                                                              
                  IF(PSP(I,ISA)<.05)PSP(I,ISA)=.05                                                                        
                  IF(PSP(I,ISA)>.75)PSP(I,ISA)=.75                                                                        
                  BPT(I,ISA)=EXP(-1.77*PSP(I,ISA)-7.05)                                                                   
                  SELECT CASE(IDSK)                                                                                       
                  CASE(1)                                                                                                 
                      SOLK(I,ISA)=MAX(.05*DHN(I,ISA),.052*DHN(I,ISA)-.12)                                                 
                      FIXK(I,ISA)=374.+236.*CLA(I,ISA)                                                                    
                  CASE(2)                                                                                                 
                      SOLK(I,ISA)=.026*DHN(I,ISA)+.5                                                                      
                      FIXK(I,ISA)=1781.+316.*CLA(I,ISA)                                                                   
                  CASE(3)                                                                                                 
                      SOLK(I,ISA)=.026*DHN(I,ISA)+.5                                                                      
                      FIXK(I,ISA)=1781.+316.*CLA(I,ISA)                                                                   
                  END SELECT                                                                                              
                  EQKS(I,ISA)=SOLK(I,ISA)/DHN(I,ISA)                                                                      
                  EQKE(I,ISA)=DHN(I,ISA)/FIXK(I,ISA)                                                                      
                  EXCK(I,ISA)=DHN(I,ISA)*WT1                                                                              
                  SOLK(I,ISA)=SOLK(I,ISA)*WT1                                                                             
                  FIXK(I,ISA)=FIXK(I,ISA)*WT1                                                                             
                  ZSK(ISA)=ZSK(ISA)+SOLK(I,ISA)                                                                           
                  ZEK(ISA)=ZEK(ISA)+EXCK(I,ISA)                                                                           
                  ZFK(ISA)=ZFK(ISA)+FIXK(I,ISA)                                                                           
                  WPMA(I,ISA)=SSF(I,ISA)*(1.-PSP(I,ISA))/PSP(I,ISA)                                                       
                  WPMS(I,ISA)=4.*WPMA(I,ISA)                                                                              
                  WSLT(I,ISA)=6.4*ECND(I,ISA)*SWST(I,ISA)                                                                 
                  ZSLT(ISA)=ZSLT(ISA)+WSLT(I,ISA)                                                                         
                  WPMA(I,ISA)=WPMA(I,ISA)*WT1                                                                             
                  WPMS(I,ISA)=WPMS(I,ISA)*WT1                                                                             
                  WPML(I,ISA)=SSF(I,ISA)*WT1                                                                              
                  ZPML(ISA)=ZPML(ISA)+WPML(I,ISA)                                                                         
                  ZPMA(ISA)=ZPMA(ISA)+WPMA(I,ISA)                                                                         
                  ZFOP(ISA)=ZFOP(ISA)+FOP(I,ISA)                                                                          
                  ZPMU(ISA)=ZPMU(ISA)+WPMU(I,ISA)                                                                         
                  ZPOU(ISA)=ZPOU(ISA)+WPOU(I,ISA)                                                                         
                  ZLSC(ISA)=ZLSC(ISA)+WLSC(I,ISA)                                                                         
                  ZLMC(ISA)=ZLMC(ISA)+WLMC(I,ISA)                                                                         
                  ZLSLC(ISA)=ZLSLC(ISA)+WLSLC(I,ISA)                                                                      
                  ZLSLNC(ISA)=ZLSLNC(ISA)+WLSLNC(I,ISA)                                                                   
                  ZBMC(ISA)=ZBMC(ISA)+WBMC(I,ISA)                                                                         
                  ZHSC(ISA)=ZHSC(ISA)+WHSC(I,ISA)                                                                         
                  ZHPC(ISA)=ZHPC(ISA)+WHPC(I,ISA)                                                                         
                  ZLSN(ISA)=ZLSN(ISA)+WLSN(I,ISA)                                                                         
                  ZLMN(ISA)=ZLMN(ISA)+WLMN(I,ISA)                                                                         
                  ZBMN(ISA)=ZBMN(ISA)+WBMN(I,ISA)                                                                         
                  ZHSN(ISA)=ZHSN(ISA)+WHSN(I,ISA)                                                                         
                  ZHPN(ISA)=ZHPN(ISA)+WHPN(I,ISA)                                                                         
                  IF(Z(I,ISA)<=PMX(ISA))THEN                                                                              
                      SUM=SUM+WT(I,ISA)                                                                                   
                      PDPL(ISA)=PDPL(ISA)+WPML(I,ISA)+WPMU(I,ISA)                                                         
                      OCPD(ISA)=OCPD(ISA)+WOC(I,ISA)                                                                      
                      PDSKC(ISA)=PDSKC(ISA)+SOLK(I,ISA)                                                                   
                      PDSW(ISA)=PDSW(ISA)+SWST(I,ISA)-S15(I,ISA)                                                          
                      PDAW(ISA)=PDAW(ISA)+FC(I,ISA)-S15(I,ISA)                                                            
                      L=I                                                                                                 
                  END IF                                                                                                  
                  WNO3(I,ISA)=CNDS(I,ISA)*WT1                                                                             
                  ZNO3(ISA)=ZNO3(ISA)+WNO3(I,ISA)                                                                         
                  IF(Z(I,ISA)<=RZ(ISA))THEN                                                                               
                      RZSW(ISA)=RZSW(ISA)+SWST(I,ISA)-S15(I,ISA)                                                          
                      PAW(ISA)=PAW(ISA)+FC(I,ISA)-S15(I,ISA)                                                              
                      LZ=I                                                                                                
                  END IF                                                                                                  
                  IF(BDD(I,ISA)<1.E-10)BDD(I,ISA)=BD(I,ISA)                                                               
                  BDP(I,ISA)=BD(I,ISA)                                                                                    
                  BDD(I,ISA)=BDD(I,ISA)/BD(I,ISA)                                                                         
                  ZPMS(ISA)=ZPMS(ISA)+WPMS(I,ISA)                                                                         
                  ZPO(ISA)=ZPO(ISA)+WPO(I,ISA)                                                                            
              END DO  !SOIL LOOP
              IF(I>MSL)THEN
                  NBSL(ISA)=MSL
              ELSE    
                  L1=LZ+1                                                                                               
                  NBSL(ISA)=I-1                                                                                               
                  IF(L1/=I)THEN                                                                                               
                      ZZ=RZ(ISA)-Z(LZ,ISA)                                                                                    
                      RTO=ZZ/(Z(L1,ISA)-Z(LZ,ISA))                                                                            
                      RZSW(ISA)=RZSW(ISA)+(SWST(L1,ISA)-S15(L1,ISA))*RTO                                                      
                      PAW(ISA)=PAW(ISA)+RTO*(FC(L1,ISA)-S15(L1,ISA))
                  END IF    
              END IF                                                                                                      
              LRD(ISA)=NBSL(ISA)                                                                                    
              IF(MXLA<NBSL(ISA))MXLA=NBSL(ISA)                                                                            
              IF(Z(1,ISA)>.01)THEN                                                                                        
                  NBSL(ISA)=NBSL(ISA)+1                                                                                   
                  DO I=NBSL(ISA),2,-1                                                                                     
                      LID(I,ISA)=LID(I-1,ISA)                                                                             
                  END DO                                                                                                  
                  LID(1,ISA)=NBSL(ISA)                                                                                    
                  LORG(NBSL(ISA),ISA)=1                                                                                   
                  RTO=.01/Z(1,ISA)                                                                                        
                  CALL SPLA(1,1,NBSL(ISA),RTO,1)                                                                          
                  Z(NBSL(ISA),ISA)=.01                                                                                    
              ELSE                                                                                                        
                  IF(Z(1,ISA)<.01)Z(1,ISA)=.01                                                                            
              END IF                                                                                                      
              IF(L/=NBSL(ISA))THEN                                                                                        
                  L1=LID(L+1,ISA)                                                                                         
                  X1=0.                                                                                                   
                  IF(L>0)X1=Z(LID(L,ISA),ISA)                                                                             
                  RTO=(PMX(ISA)-X1)/(Z(L1,ISA)-X1)                                                                        
                  SUM=SUM+WT(L1,ISA)*RTO                                                                                  
                  OCPD(ISA)=OCPD(ISA)+RTO*WOC(L1,ISA)                                                                     
                  PDPLX(ISA)=PDPL(ISA)+(WPML(L1,ISA)+WPMU(L1,ISA))*RTO                                                    
                  PDPL(ISA)=1000.*PDPLX(ISA)/SUM                                                                          
                  PDPL0(ISA)=PDPL(ISA)                                                                                    
                  PDPLC(ISA)=PDPL(ISA)                                                                                    
                  PDSKC(ISA)=PDSKC(ISA)+RTO*SOLK(L1,ISA)                                                                  
                  PDSKC(ISA)=1000.*PDSKC(ISA)/SUM                                                                         
                  PDSW(ISA)=PDSW(ISA)+RTO*(SWST(L1,ISA)-S15(L1,ISA))                                                      
                  PDAW(ISA)=PDAW(ISA)+RTO*(FC(L1,ISA)-S15(L1,ISA))                                                        
              ELSE                                                                                                        
                  PMX(ISA)=Z(NBSL(ISA),ISA)                                                                               
              END IF                                                                                                      
              OCPD(ISA)=.1*OCPD(ISA)/SUM                                                                                  
              SUM=1.E-4*SUM/PMX(ISA)                                                                                      
              ABD(ISA)=1.E-4*SEV(3,ISA)/XX
              DO WHILE(NBSL(ISA)<MXLA)
                  L1=LID(1,ISA)                                                                                               
                  ZMX=0.                                                                                                      
                  MXZ=2                                                                                                       
                  DO I=2,NBSL(ISA)                                                                                            
                      ISL=LID(I,ISA)                                                                                          
                      ZZ=Z(ISL,ISA)-Z(L1,ISA)                                                                                 
                      IF(ZZ>=ZTK)THEN                                                                                         
                          MXZ=I                                                                                               
                          EXIT
                      END IF                                                                                                  
                      IF(ZZ>ZMX+.01)THEN                                                                                      
                          ZMX=ZZ                                                                                              
                          MXZ=I                                                                                               
                      END IF                                                                                                  
                      L1=ISL                                                                                                  
                  END DO
                  IF(I>NBSL(ISA))THEN
                      ISL=LID(MXZ,ISA)                                                                                            
                      L1=LID(MXZ-1,ISA)  
                  END IF
                  NBSL(ISA)=NBSL(ISA)+1                                                                                 
                  CALL SPLA(ISL,L1,NBSL(ISA),.5,0)                                                                            
                  DO K=NBSL(ISA),MXZ+1,-1                                                                                     
                      LID(K,ISA)=LID(K-1,ISA)                                                                                 
                  END DO                                                                                                      
                  LID(MXZ,ISA)=NBSL(ISA)                                                                                      
                  LORG(NBSL(ISA),ISA)=LORG(ISL,ISA)                                                                           
              END DO
              DO J=1,NBSL(ISA)                                                                                      
                  L=LID(J,ISA)                                                                                            
                  XGO2(L,isa)=CGO2(L,ISA)    ! Original statement XGO2(L)=CGO2(L,ISA) - Luca Doro
                  XGCO2(L,isa)=CGCO2(L,ISA)    ! Original statement XGCO2(L)=CGCO2(L,ISA) - Luca Doro
                  XGN2O(L,isa)=CGN2O(L,ISA)    ! Original statement XGN2O(L)=CGN2O(L,ISA) - Luca Doro
                  IF(NCC>0)THEN                                                                                           
                      WOC(L,ISA)=WLSC(L,ISA)+WLMC(L,ISA)+WBMC(L,ISA)+WHSC(L,ISA)+WHPC(L,ISA)                              
                      WON(L,ISA)=WLSN(L,ISA)+WLMN(L,ISA)+WBMN(L,ISA)+WHSN(L,ISA)+WHPN(L,ISA)                              
                  ELSE                                                                                                    
                      WLSC(L,ISA)=.42*WLS(L,ISA)                                                                          
                      WLMC(L,ISA)=.42*WLM(L,ISA)                                                                          
                      WLSLC(L,ISA)=.42*WLSL(L,ISA)                                                                        
                      WLSLNC(L,ISA)=WLSC(L,ISA)-WLSLC(L,ISA)                                                              
                  END IF                                                                                                  
                  SOL(1,L,ISA)=WHSC(L,ISA)                                                                                
                  SOL(2,L,ISA)=WHPC(L,ISA)                                                                                
                  SOL(3,L,ISA)=WLSC(L,ISA)                                                                                
                  SOL(4,L,ISA)=WLMC(L,ISA)                                                                                
                  SOL(5,L,ISA)=WBMC(L,ISA)                                                                                
                  SOL(6,L,ISA)=WOC(L,ISA)                                                                                 
                  SOL(7,L,ISA)=WHSN(L,ISA)                                                                                
                  SOL(8,L,ISA)=WHPN(L,ISA)                                                                                
                  SOL(9,L,ISA)=WLSN(L,ISA)                                                                                
                  SOL(10,L,ISA)=WLMN(L,ISA)                                                                               
                  SOL(11,L,ISA)=WBMN(L,ISA)                                                                               
                  SOL(12,L,ISA)=WON(L,ISA)                                                                                
                  SOL(13,L,ISA)=WPMA(L,ISA)                                                                               
                  SOL(14,L,ISA)=WPMS(L,ISA)                                                                               
                  SOL(15,L,ISA)=WPO(L,ISA)                                                                                
                  SOL(16,L,ISA)=EXCK(L,ISA)                                                                               
                  SOL(17,L,ISA)=FIXK(L,ISA)                                                                               
                  SOL(18,L,ISA)=SWST(L,ISA)                                                                               
                  SOL(19,L,ISA)=WLS(L,ISA)                                                                                
                  SOL(20,L,ISA)=WLM(L,ISA)                                                                                
                  SOL(21,L,ISA)=WLSL(L,ISA)                                                                               
                  SOL(22,L,ISA)=WLSLC(L,ISA)                                                                              
                  SOL(23,L,ISA)=WLSLNC(L,ISA)                                                                             
                  XZP(1,L,ISA)=WHSC(L,ISA)                                                                                
                  XZP(2,L,ISA)=WHPC(L,ISA)                                                                                
                  XZP(3,L,ISA)=WLSC(L,ISA)                                                                                
                  XZP(4,L,ISA)=WLMC(L,ISA)                                                                                
                  XZP(5,L,ISA)=WBMC(L,ISA)                                                                                
                  XZP(6,L,ISA)=WOC(L,ISA)                                                                                 
                  XZP(7,L,ISA)=WHSN(L,ISA)                                                                                
                  XZP(8,L,ISA)=WHPN(L,ISA)                                                                                
                  XZP(9,L,ISA)=WLSN(L,ISA)                                                                                
                  XZP(10,L,ISA)=WLMN(L,ISA)                                                                               
                  XZP(11,L,ISA)=WBMN(L,ISA)                                                                               
                  XZP(12,L,ISA)=WON(L,ISA)                                                                                
                  XZP(13,L,ISA)=WOC(L,ISA)/WON(L,ISA)                                                                     
              END DO                                                                                                      
              ZON(ISA)=ZLSN(ISA)+ZLMN(ISA)+ZBMN(ISA)+ZHSN(ISA)+ZHPN(ISA)                                                  
              ZOC(ISA)=ZLSC(ISA)+ZLMC(ISA)+ZBMC(ISA)+ZHSC(ISA)+ZHPC(ISA)                                                  
              BTC(ISA)=ZOC(ISA)*WSA(ISA)                                                                                           
              BTCX(ISA)=BTC(ISA)                                                                                          
              BTCZ(ISA)=.001*BTC(ISA)                                                                                          
              SLT0(ISA)=ZSLT(ISA)                                                                                         
              SLTX(ISA)=ZSLT(ISA)                                                                                         
              XZP(1,11,ISA)=ZHSC(ISA)                                                                                     
              XZP(2,11,ISA)=ZHPC(ISA)                                                                                     
              XZP(3,11,ISA)=ZLSC(ISA)                                                                                     
              XZP(4,11,ISA)=ZLMC(ISA)                                                                                     
              XZP(5,11,ISA)=ZBMC(ISA)                                                                                     
              XZP(6,11,ISA)=ZOC(ISA)                                                                                      
              XZP(7,11,ISA)=ZHSN(ISA)                                                                                     
              XZP(8,11,ISA)=ZHPN(ISA)                                                                                     
              XZP(9,11,ISA)=ZLSN(ISA)                                                                                     
              XZP(10,11,ISA)=ZLMN(ISA)                                                                                    
              XZP(11,11,ISA)=ZBMN(ISA)                                                                                    
              XZP(12,11,ISA)=ZON(ISA)                                                                                     
              XZP(13,11,ISA)=ZOC(ISA)/ZON(ISA)                                                                            
              IF(IDNT>2)THEN                                                                                              
                  NBCL=Z(LID(NBSL(ISA),ISA),ISA)/DZDN+.999                                                                
                  IF(NBCL>30)THEN                                                                                         
                      NBCL=30                                                                                             
                      DZDN=Z(LID(NBSL(ISA),ISA),ISA)/30.                                                                  
                  ELSE                                                                                                    
                      IF(NBCL<NBSL(ISA))THEN                                                                              
                          NBCL=NBSL(ISA)                                                                                  
                          X1=NBCL                                                                                         
                          DZDN=Z(LID(NBSL(ISA),ISA),ISA)/X1                                                               
                      END IF                                                                                              
                  END IF                                                                                                  
                  DZ10=10.*DZDN                                                                                           
                  TOT=0.                                                                                                  
                  DO I=1,NBCL                                                                                             
                      TOT=TOT+DZDN                                                                                        
                      ZC(I,ISA)=TOT                                                                                       
                  END DO                                                                                                  
                  CALL AINTRIC(XGO2,CGO2,NBSL(ISA),NBCL)                                                                  
                  CALL AINTRIC(XGCO2,CGCO2,NBSL(ISA),NBCL)                                                                
                  CALL AINTRIC(XGN2O,CGN2O,NBSL(ISA),NBCL)                                                                
                  IUN=NBCL-1                                                                                              
              END IF                                                                                                      
              IF(LBP>0)THEN                                                                                               
                  X1=CLA(LID(2,ISA),ISA)                                                                                  
                  CPMX(ISA)=1000.*X1/(X1+EXP(3.519-.027*X1))                                                              
              END IF                                                                                                      
              CALL SPRNT                                                                                                  
              IF(KFL(1)>0)THEN                                                                                            
                  WRITE(KW(1),2950)                                                                                       
                  WRITE(KW(1),2960)FL,FW,ANG0,UXP,ACW                                                                     
                  WRITE(KW(1),'(/T10,A,F8.0,A)')'ANNUAL HEAT UNITS = ',AHSM,' C'                                          
                  CALL APAGE(1)                                                                                           
                  WRITE(KW(1),'(//1X,A/)')'____________________SOIL DATA____________&                                     
                  ________'                                                                                               
                  WRITE(KW(1),3000)SALB(ISA),MXLA,ZQT,ZF,ZTK,FBM(ISA),FHP(ISA),XIDS&                                      
                  (ISA),OCPD(ISA),RTN1                                                                                    
                  WRITE(KW(1),3001)IWTB,WTMN(ISA),WTMX(ISA),WTBL(ISA)                                                     
                  SELECT CASE(ISW)                                                                                         
                  CASE(0)                                                                                                 
                      WRITE(KW(1),'(T10,A)')'FIELD CAP/WILTING PT EST RAWLS &                                             
                      METHOD DYNAMIC'                                                                                     
                  CASE(1)                                                                                                 
                      WRITE(KW(1),'(T10,A)')'FIELD CAP/WILTING PT INP RAWLS &                                             
                      METHOD DYNAMIC'                                                                                     
                  CASE(2)                                                                                                 
                      WRITE(KW(1),'(T10,A)')'FIELD CAP/WILTING PT EST RAWLS &                                             
                      METHOD STATIC'                                                                                      
                  CASE(3)                                                                                                 
                      WRITE(KW(1),'(T10,A)')'FIELD CAP/WILTING PT INP STATIC '                                            
                  CASE(4)                                                                                                 
                      WRITE(KW(1),'(T10,A)')'FIELD CAP/WILTING PT EST NEAREST &                                           
                      NEIGHBOR METHOD DYNAMIC'                                                                            
                  CASE(5)                                                                                                 
                      WRITE(KW(1),'(T10,A)')'FIELD CAP/WILTING PT EST NEAREST &                                           
                      NEIGHBOR METHOD STATIC'	                                                                            
                  CASE(6)                                                                                                 
                      WRITE(KW(1),'(T10,A)')'FIELD CAP/WILTING PT EST BNW &                                               
                      METHOD DYNAMIC'	                                                                                    
                  CASE(7)                                                                                                 
                      WRITE(KW(1),'(T10,A)')'FIELD CAP/WILTING PT EST BNW &                                               
                      METHOD STATIC'	                                                                                     
                  END SELECT                                                                                              
              END IF                                                                                                      
              IF(KFL(1)>0.AND.ISTA>0)WRITE(KW(1),'(T10,A)')'STATIC SOIL PROFILE'                                                       
              IF(SAT1>0.)THEN                                                                                             
                  SATZ=SAT1                                                                                               
              ELSE                                                                                                        
                  IF(SAT0>0.)THEN                                                                                         
                      SATZ=SAT0                                                                                           
                  ELSE                                                                                                    
                      SATZ=1.                                                                                             
                  END IF                                                                                                  
              END IF                                                                                                      
              SATK(ISA)=SATC(LID(2,ISA),ISA)*SATZ                                                                         
              IF(FPS1>0.)THEN                                                                                             
                  FPSZ=FPS1                                                                                               
              ELSE                                                                                                        
                  IF(FPS0>0.)THEN                                                                                         
                      FPSZ=FPS0                                                                                           
                  ELSE                                                                                                    
                      FPSZ=1.                                                                                             
                  END IF                                                                                                  
              END IF                                                                                                      
              FPSC(ISA)=SATC(LID(2,ISA),ISA)*FPSZ                                                                         
              SW(ISA)=UW(2)*XX*1000.                                                                                      
              ZCOB(ISA)=ZCO(ISA)                                                                                          
              I2=-1                                                                                                       
              DO WHILE(I2/=IOPS)                                                                                          
                  READ(KR(15),*,IOSTAT=NFL)I2,OPSCFILE                                                                    
                  IF(NFL/=0)THEN                                                                                          
                      WRITE(*,*)'OPS NO = ',IOPS,' NOT IN OPSC LIST FILE &                                                
                      SAID = ',NBSA(ISA)                                                                                  
                      STOP                                                                                                
                  END IF                                                                                                  
              END DO                                                                                                      
              REWIND KR(15)                                                                                               
              CALL OPENV(KR(16),OPSCFILE,IDIR)                                                                            
              TITOP(ISA)=OPSCFILE                                                                                         
              ! LINE 1                                                                                                    
              READ(KR(16),'()')                                                                                           
              ! 1 LUN  = LAND USE NUMBER FROM NRCS LAND USE-HYDROLOGIC SOIL GROUP                                         
              !          TABLE                                                                                            
              ! 2 AUTO APPLICATION EQUIPMENT FROM TILLCOM.DAT                                                             
              ! 3 IAUI = AUTO IRRIGATION                                                                                  
              ! 4 IAUF = AUTO N FERT                                                                                      
              ! 5 IAMF = AUTO MANURE FROM GRAZING ANIMALS                                                                 
              ! 6 ISPF = AUTO MANURE APPLICATION FROM FEEDYARD STOCKPILE                                                  
              ! 7 ILQF = AUTO LIQUID MANURE APPLICATION FROM LAGOONS                                                      
              ! 8 IAUL = AUTO LIME APPLICATION                                                                            
              !     LINE 2                                                                                                
              READ(KR(16),3100)LUN(ISA),IAUI(ISA),IAUF(ISA),IAMF(ISA),ISPF(ISA),&                                         
              ILQF(ISA),IAUL(ISA)                                                                                         
              IF(IAUI(ISA)==0)IAUI(ISA)=500                                                                               
              IF(IAUF(ISA)==0)IAUF(ISA)=261                                                                               
              IF(IAMF(ISA)==0)IAMF(ISA)=268                                                                               
              IF(ISPF(ISA)==0)ISPF(ISA)=266                                                                               
              IF(ILQF(ISA)==0)ILQF(ISA)=265                                                                               
              IF(IAUL(ISA)==0)IAUL(ISA)=267                                                                               
              IF(LUNS(ISA)/=0)THEN                                                                                        
                  I1=LUNS(ISA)                                                                                            
                  LUNS(ISA)=LUNS(ISA)-LUN(ISA)                                                                            
                  LUN(ISA)=I1                                                                                             
              END IF                                                                                                      
              !     READ MANAGEMENT INFORMATION                                                                           
              !  1  IRR  = N0 FOR DRYLAND AREAS          | N = 0 APPLIES MINIMUM OF                                       
              !          = N1 FROM SPRINKLER IRRIGATION  | VOLUME INPUT, ARMX, & FC-SW                                    
              !          = N2 FOR FURROW IRRIGATION      | N = 1 APPLIES INPUT VOLUME                                     
              !          = N3 FOR IRR WITH FERT ADDED    | OR ARMX                                                        
              !          = N4 FOR IRRIGATION FROM LAGOON |                                                                
              !          = N5 FOR DRIP IRR               |                                                                
              !  2  IRI  = N DAY APPLICATION INTERVAL FOR AUTOMATIC IRRIGATION                                            
              !  3  IFA  = MIN FERT APPL INTERVAL(0 FOR USER SPECIFIED)                                               
              !  4  LM   = 0 APPLIES LIME                                                                                 
              !          = 1 DOES NOT APPLY LIME                                                                          
              !  5  IFD  = 0 WITHOUT FURROW DIKES                                                                         
              !            1 WITH FURROW DIKES                                                                            
              !  6  IDR  = 0 NO DRAINAGE                                                                                  
              !          = DEPTH OF DRAINAGE SYSTEM(mm)                                                                   
              !     IDF0 = FERT #                                                                                         
              !  7         1 FOR FERTIGATION FROM LAGOON                                                                  
              !  8         2 FOR AUTOMATIC SOLID MANURE APPL FROM FEEDING AREA STOCK                                      
              !              PILE.                                                                                        
              !  9         3 AUTO COMERCIAL P FERT APPL (DEFAULTS TO ELEM P)                                              
              !  10        4 FOR AUTOMATIC COMERCIAL FERT APPL(DEFALTS TO ELEM N)                                         
              !  11        5 FOR AUTOMATIC SOLID MANURE APPLICATION.                                                      
              !  12        6 AUTO COMERCIAL K FERT APPL (DEFAULTS TO ELEM K)                                              
              !  13 IRRS = ID OF SA SUPPLYING IRRIGATION WATER FROM A RESERVOIR                                           
              !            0 NO RESERVOIR SUPPLY OR NO IRRIGATION
              !  14 IRRW = ID OF SA SUPPLYING IRRIGATION WATER FROM A WELL
              !            0 NO WELL SUPPLY
              !     LINE 8                                                                                                
              READ(KR(5),*)IRR(ISA),IRI(ISA),IFA(ISA),LM(ISA),IFD(ISA),IDR&                                               
              (ISA),(IDF0(I,ISA),I=1,6),IRRS(ISA),IRRW(ISA)                                                                         
              IF(IDR(ISA)>0.AND.HSG>2)HSG=2                                                                               
              ISG(ISA)=HSG                                                                                                
              CALL HSGCN                                                                                                  
              CALL HCNSLP(CN2(ISA),X3) 
              
              !rtb hsg - if requested, calculate daily values of CNSC
              if(hsg_vary) then
                year_index = IYR0
                do yrc=1,NBYR
                  !number of days in the year
                  if(mod(year_index,4) == 0) then
                    year_days = 366
                  else
                    year_days = 365
                  endif
                  do dayc=1,year_days
                    X3_const = X3
                    ISG(ISA) = hsg_values(INPS,yrc,dayc) !INPS = soil type for the subarea   
                    CALL HSGCN !get CN2 value
                    call HCNSLP(CN2(ISA),X3_const) !get CNSC values
                    cnsc1_vary(ISA,yrc,dayc) = CNSC(1,ISA)
                    cnsc2_vary(ISA,yrc,dayc) = CNSC(2,ISA)
                  enddo
                  year_index = year_index + 1  
							  enddo
              endif !rtb hsg
                            
              SCI(ISA)=MAX(3.,SMX(ISA)*(1.-RZSW(ISA)/PAW(ISA)))                                                           
              IF(KFL(1)>0)THEN                                                                                            
                  SELECT CASE(NVCN(ISA))                                                                                
                  CASE(0)                                                                                                 
                      WRITE(KW(1),'(T10,A)')'VARIABLE CN NONLINEAR DEPTH/SW &                                             
                      WEIGHTING'                                                                                          
                  CASE(1)                                                                                                 
                      WRITE(KW(1),'(T10,A)')'VARIABLE CN NONLINEAR NO DEPTH/SW &                                          
                      WEIGHTING'                                                                                          
                  CASE(2)                                                                                                 
                      WRITE(KW(1),'(T10,A)')'VARIABLE CN LINEAR NO DEPTH/SW &                                             
                      WEIGHTING'                                                                                          
                  CASE(3)                                                                                                 
                      WRITE(KW(1),'(T10,A)')'CONSTANT CN'                                                                 
                  CASE(4)                                                                                                 
                      WRITE(KW(1),'(T10,A)')'VARIABLE CN SOIL MOIST INDEX NO &                                            
                      DEPTH/SW WEIGHTING'                                                                                 
                  END SELECT                                                                                              
                  WRITE(KW(1),3970)ASG(ISG(ISA)),LUN(ISA),CN2(ISA),X3,&                                                   
                  SCNX(ISA),SCRP(35,1),SCRP(35,2),CNSC(1,ISA),CNSC(2,ISA)                                                 
              END IF                                                                                                      
              CN2(ISA)=X3                                                                                                 
              CN0(ISA)=X3                                                                                                 
              IF(KFL(1)>0)THEN                                                                                            
                  CALL APAGE(1)                                                                                           
                  WRITE(KW(1),'(//1X,A/)')'____________________SOIL PHYSICAL DATA___&                                     
                  _________________'                                                                                      
                  CALL SOLIOP                                                                                             
                  CALL APAGE(1)                                                                                           
                  WRITE(KW(1),'(//1X,A/)')'____________________SOIL CHEMICAL DATA___&                                     
                  _________________'                                                                                      
                  CALL SOLIOC                                                                                             
              END IF                                                                                                      
              !  1  BIR  = IRRIGATION TRIGGER--3 OPTIONS                                                                  
              !            1. PLANT WATER STRESS FACTOR (0-1)                                                             
              !            2. SOIL WATER TENSION IN TOP 200 MM(> 1 KPA)                                                   
              !            3. PLANT AVAILABLE WATER DEFICIT IN ROOT ZONE (-mm)                                            
              !  2  EFI  = RUNOFF VOL / VOL IRR WATER APPLIED(0 IF IRR=0)                                             
              !  3  VIMX = MAXIMUM ANNUAL IRRIGATION VOLUME ALLOWED FOR EACH CROP (mm)                                    
              !  4  ARMN = MINIMUM SINGLE APPLICATION VOLUME ALLOWED (mm)                                                 
              !  5  ARMX = MAXIMUM SINGLE APPLICATION VOLUME ALLOWED (mm)                                                 
              !  6  BFT  = AUTO FERTILIZER TRIGGER--2 OPTIONS                                                             
              !            1. PLANT N STRESS FACTOR (0-1)                                                                 
              !            2. SOIL N CONC IN ROOT ZONE (g/t)                                                              
              !  7  FNP4 = AUTO FERT FIXED APPLICATION RATE (kg/ha)                                                       
              !  8  FMX  = MAXIMUM ANNUAL N FERTILIZER APPLICATION FOR A CROP (kg/ha)                                     
              !  9  DRT  = TIME REQUIRED FOR DRAINAGE SYSTEM TO REDUCE PLANT STRESS(d)                                    
              !            (0 IF DRAINAGE NOT USED)                                                                   
              ! 10  FDSF = FURROW DIKE SAFETY FACTOR(0-1.)                                                                
              ! 11  PEC  = CONSERVATION PRACTICE FACTOR(=0.0 ELIMINATES WATER EROSION)                                    
              ! 12  DALG = FRACTION OF SUBAREA CONTROLLED BY LAGOON.                                                      
              ! 13  VLGN = LAGOON VOLUME RATIO--NORMAL / MAXIMUM                                                          
              ! 14  COWW = LAGOON INPUT FROM WASH WATER (m3/hd/d)                                                         
              ! 15  DDLG = TIME TO REDUCE LAGOON STORAGE FROM MAX TO NORM (d)                                             
              ! 16  SOLQ = RATIO LIQUID/TOTAL MANURE PRODUCED.                                                            
              ! 17  SFLG = SAFETY FACTOR FOR LAGOON DESIGN (VLG=VLG0/(1.-SFLG)                                            
              ! 18  FNP2 = FEEDING AREA STOCK PILE AUTO SOLID MANURE APPL RATE (kg/ha)                                    
              ! 19  FNP5 = AUTOMATIC MANURE APPLICATION RATE (kg/ha)                                                      
              ! 20  FIRG = FACTOR TO ADJUST AUTO IRRIGATION VOLUME (FIRG*FC)                                              
              !     LINE 9/10                                                                                             
              READ(KR(5),*)BIR(ISA),EFI(ISA),VIMX(ISA),ARMN(ISA),ARMX(ISA),&                                              
              BFT(ISA),FNP(4,ISA),FMX,DRT(ISA),FDSF(ISA),PEC(ISA),DALG(ISA),VLGN&                                         
              (ISA),COWW(ISA),DDLG(ISA),SOLQ(ISA),SFLG,FNP(2,ISA),FNP(5,ISA),&                                            
              FIRG(ISA)                                                                                                   
              IF(FIRG(ISA)<1.E-10)THEN                                                                                    
                  IF(FIR0<1.E-10)THEN                                                                                     
                      FIRG(ISA)=1.                                                                                        
                  ELSE                                                                                                    
                      FIRG(ISA)=FIR0                                                                                      
                  END IF                                                                                                  
              END IF                                                                                                      
              NZ=MAX(1,NHRD(IOW))                                                                                         
              !     NY   = 0 FOR NON GRAZING AREA                                                                         
              !          = HERD NUMBERS FOR GRAZING AREA                                                                  
              !     LINE 11                                                                                               
              READ(KR(5),*)(NY(J),J=1,NZ)                                                                                 
              !     XTP  = GRAZING LIMIT FOR EACH HERD--MINIMUM PLANT MATERIAL(t/ha)                                      
              !     LINE 12                                                                                               
              READ(KR(5),*)(XTP(J),J=1,NZ)                                                                                
              DO J=1,NZ                                                                                               
                  J1=NY(J)                                                                                            
                  IF(J1==0)CYCLE                                                                                      
                  NGZ(J1,ISA)=J1                                                                                      
                  GZLM(J1,ISA)=XTP(J)                                                                                 
              END DO                                                                                                  
              IHDM(ISA)=NY(1)                                                                                         
              IF(DRT(ISA)<1.E-5)DRT(ISA)=1.                                                                               
              IF(KFL(1)>0)THEN                                                                                            
                  CALL APAGE(1)                                                                                           
                  WRITE(KW(1),'(//1X,A/)')'____________________MANAGEMENT DATA&                                           
                  ____________________'                                                                                   
                  WRITE(KW(1),'(T10,A,I4)')'OWNER ID=',IOW                                                                
                  IF(II>0)WRITE(KW(1),'(T10,A/T15,A,I3)')'FEEDING AREA','HERD ID=',II                                     
              END IF                                                                                                      
              DO IHD=1,NZ                                                                                                 
                  IF(NGZ(IHD,ISA)>0)EXIT                                                                                  
              END DO                                                                                                      
              IF(KFL(1)>0.AND.IHD<=NZ)THEN                                                                                             
                  WRITE(KW(1),'(T10,A/T15,A)')'GRAZING AREA','GRAZING MODE'                                               
                  SELECT CASE(IHRD)                                                                                     
                  CASE(0)                                                                                                 
                      WRITE(KW(1),'(T15,A)')'LEVEL 0 MANUAL_NO HERD FILE REQUIRED'                                        
                  CASE(1)                                                                                                 
                      WRITE(KW(1),'(T15,A)')'LEVEL 1 HYBRID_HERD FILE REQUIRED'                                           
                  CASE(2)                                                                                                 
                      WRITE(KW(1),'(T15,A)')'LEVEL 2 AUTOMATIC_HERD FILE REQUIRED'                                        
                  END SELECT
                  IF(IHAY>0)THEN    
                      WRITE(KW(1),'(T15,A)')'HAY FEED AS NEEDED'
                  ELSE
                      WRITE(KW(1),'(T15,A)')'NO HAY FEED--GRAZERS REMOVED WHEN AGPM<GZLM'
                  END IF    
              END IF                                                                                                      
              DO IHD=1,NZ                                                                                                 
                  IF(NGZ(IHD,ISA)==0)CYCLE                                                                                
                  IF(GZLM(IHD,ISA)<1.E-5)GZLM(IHD,ISA)=GZL0                                                               
                  X1=24.*FFED(IHD,IOW)                                                                                    
                  IF(KFL(1)>0)WRITE(KW(1),3536)IHD,NCOW(IHD,IOW),X1,GZLM(IHD,ISA)                                         
              END DO                                                                                                      
              IFA(ISA)=MAX(IFA(ISA),1)                                                                                    
              CALL AISPL(IRR(ISA),IAC(ISA))                                                                               
              IF(ARMX(ISA)<1.E-10)ARMX(ISA)=1000.                                                                         
              IF(FMX<1.E-10)FMX=200.                                                                                      
              IF(KFL(1)>0)THEN                                                                                            
                  IF(IAPL(ISA)<0)WRITE(KW(1),'(T10,A)')'LIQUID MANURE APPL AREA'                                          
                  IF(IAPL(ISA)>0)WRITE(KW(1),'(T10,A)')'SOLID MANURE APPL AREA'                                           
                  IF(FFPQ(ISA)>0)WRITE(KW(1),'(T10,A)')'FILTER STRIP'                                                     
              END IF                                                                                                      
              IF(IRR(ISA)==0)THEN                                                                                         
                  IF(DALG(ISA)>0.)THEN                                                                                    
                      X3=COWW(ISA)*NCOW(IDFH(ISA),IOW)                                                                    
                      RFM0=2.*RFMX                                                                                        
                      X1=RFM0-5.64                                                                                        
                      QLG=X1*X1/(RFM0+22.6)                                                                               
                      DALG(ISA)=DALG(ISA)*WSA(ISA)                                                                        
                      X1=10.*DALG(ISA)                                                                                    
                      QWW=30.*X3/X1                                                                                       
                      VLGM(ISA)=(QLG+QWW)/(1.-VLGN(ISA))                                                                  
                      VLGN(ISA)=VLGN(ISA)*VLGM(ISA)                                                                       
                      VLGM(ISA)=VLGM(ISA)/(1.-SFLG)                                                                       
                      IF(KFL(1)>0)WRITE(KW(1),39)DALG(ISA),VLGN(ISA),VLGM(ISA),DDLG(ISA)&                                 
                      ,COWW(ISA),SFLG                                                                                     
                      COWW(ISA)=X3                                                                                        
                      VLGN(ISA)=X1*VLGN(ISA)                                                                              
                      VLGB(ISA)=VLGN(ISA)                                                                                 
                      VLGM(ISA)=X1*VLGM(ISA)                                                                              
                      VLG(ISA)=VLGN(ISA)                                                                                  
                      CFNP(ISA)=0.                                                                                        
                      WTMU(ISA)=CFNP(ISA)*VLG(ISA)                                                                        
                      WTMB(ISA)=WTMU(ISA)                                                                                 
                      TWMB=TWMB+WTMU(ISA)                                                                                 
                      VLGI(ISA)=(VLGM(ISA)-VLGN(ISA))/(DDLG(ISA)+1.E-5)                                                   
                  ELSE                                                                                                    
                      IF(KFL(1)>0)WRITE(KW(1),'(T10,A)')'DRYLAND AGRICULTURE'                                             
                  END IF                                                                                                  
              ELSE                                                                                                        
                  IF(VIMX(ISA)<1.E-10)VIMX(ISA)=2000.                                                                     
                  IF(KFL(1)>0)THEN                                                                                        
                      IF(BIR(ISA)<0.)THEN                                                                                 
                          WRITE(KW(1),'(T10,A)')'AUTOMATIC IRRIGATION'                                                    
                          WRITE(KW(1),4080)BIR(ISA),IRI(ISA)                                                              
                      ELSE                                                                                                
                          IF(BIR(ISA)>0.)THEN                                                                             
                              WRITE(KW(1),'(T10,A)')'AUTOMATIC IRRIGATION'                                                
                              IF(BIR(ISA)<1.)THEN                                                                         
                                  WRITE(KW(1),3680)BIR(ISA),IRI(ISA)                                                      
                              ELSE                                                                                        
                                  WRITE(KW(1),3630)BIR(ISA),IRI(ISA)                                                      
                              END IF                                                                                      
                          ELSE                                                                                            
                              WRITE(KW(1),'(T10,A)')'USER SPECIFIED IRRIGATION'                                           
                          END IF                                                                                          
                      END IF                                                                                              
                      WRITE(KW(1),3620)VIMX(ISA),ARMN(ISA),ARMX(ISA)                                                      
                      IF(IAC(ISA)==0)THEN                                                                                 
                          WRITE(KW(1),'(T15,A)')'FLEXIBLE APPL VOLUMES'                                                   
                      ELSE                                                                                                
                          WRITE(KW(1),'(T15,A)')'RIGID APPL VOLUMES'                                                      
                      END IF                                                                                              
                      SELECT CASE(IRR(ISA))                                                                               
                      CASE(1)                                                                                             
                          WRITE(KW(1),'(T15,A)')'SPRINKLER IRRIGATION'                                                    
                          WRITE(KW(1),3140)EFI(ISA)                                                                       
                      CASE(2)                                                                                             
                          WRITE(KW(1),'(T15,A)')'FURROW IRRIGATION'                                                       
                          WRITE(KW(1),3140)EFI(ISA)                                                                       
                      CASE(3)                                                                                             
                          WRITE(KW(1),'(T15,A)')'IRRIGATION WITH FERT ADDED'                                              
                          WRITE(KW(1),3140)EFI(ISA)                                                                       
                      CASE(4)                                                                                             
                          WRITE(KW(1),'(T15,A)')'IRRIGATION FROM LAGOON'                                                  
                      CASE(5)                                                                                             
                          WRITE(KW(1),'(T15,A)')'DRIP IRRIGATION'	                                                        
                      END SELECT                                                                                          
                      WRITE(KW(1),'(T15,A,F6.3)')'FRACTION FC TO CALCULATE AUTO IRR VOL =',&                              
                      FIRG(ISA)                                                                                           
                  END IF                                                                                                  
              END IF                                                                                                      
              IF(KFL(1)>0)THEN	                                                                                           
                  IF(BFT(ISA)>0.)THEN                                                                                     
                      WRITE(KW(1),'(T10,A/T15,A,I3,A)')'AUTO SCHEDULED FERT','MIN APPL &                                  
                      INTERVAL = ',IFA(ISA),' d'                                                                          
                      IF(BFT(ISA)>1.)THEN                                                                                 
                          WRITE(KW(1),'(T15,A,F4.0,A)')'SOIL SOL N CONC TRIGGER = ',BFT(ISA)&                             
                          ,' g/t'                                                                                         
                      ELSE                                                                                                
                          WRITE(KW(1),'(T15,A,F5.2)')'PLANT STRESS TRIGGER = ',BFT(ISA)                                   
                      END IF                                                                                              
                      WRITE(KW(1),'(T15,A,F8.1,A)')'FIXED RATE, SURFACE APPLIED = ',FNP&                                  
                      (4,ISA),' kg/ha/APP'                                                                                
                  ELSE                                                                                                    
                      WRITE(KW(1),'(T10,A)')'USER SPECIFIED FERT'                                                         
                  END IF                                                                                                  
                  IF(IDF0(5,ISA)>0.AND.FNP(5,ISA)>0.)WRITE(KW(1),'(T15,A)')'AUTO MANURE APPL'                             
                  WRITE(KW(1),'(T10,A,F8.0,A)')'MAX N FERT/CROP = ',FMX,' kg/ha'                                          
                  IF(IPAT>0)THEN                                                                                          
                      WRITE(KW(1),'(T10,A)')'AUTO P FERT'                                                                 
                  ELSE                                                                                                    
                      WRITE(KW(1),'(T10,A)')'MANUAL P FERT'                                                               
                  END IF                                                                                                  
                  IF(IKAT>0)THEN                                                                                          
                      WRITE(KW(1),'(T10,A)')'AUTO K FERT'                                                                 
                  ELSE                                                                                                    
                      WRITE(KW(1),'(T10,A)')'MANUAL K FERT'                                                               
                  END IF                                                                                                  
              END IF                                                                                                      
              NII(ISA)=IRI(ISA)                                                                                           
              IF(IFD(ISA)>0)THEN                                                                                          
                  IF(FDSF(ISA)<1.E-10)FDSF(ISA)=.9                                                                        
                  IF(KFL(1)>0)WRITE(KW(1),2560)FDSF(ISA)                                                                  
                  FDSF(ISA)=FDSF(ISA)*1000.                                                                               
              END IF                                                                                                      
              IF(KFL(1)>0)THEN                                                                                            
                  IF(LM(ISA)==0)THEN                                                                                      
                      WRITE(KW(1),'(T10,A)')'LIME APPLIED AS NEEDED'                                                      
                  ELSE                                                                                                    
                      WRITE(KW(1),'(T10,A)')'NO LIME APPLICATIONS'                                                        
                  END IF                                                                                                  
                  WRITE(KW(1),'(T10,A,F10.3)')'USLE P FACTOR = ',PEC(ISA)                                                 
              END IF                                                                                                      
              IF(IDR(ISA)>0)THEN                                                                                          
                  X1=IDR(ISA)                                                                                             
                  X1=.001*IDR(ISA)                                                                                        
                  DO I=1,NBSL(ISA)                                                                                        
                      L=LID(I,ISA)                                                                                        
                      IF(Z(L,ISA)>X1)EXIT                                                                                 
                  END DO                                                                                                  
                  IDR(ISA)=L                                                                                              
                  HCLN(ISA)=HCL(L,ISA)                                                                                    
                  HCL(L,ISA)=MAX(PRMT(83)*SATC(L,ISA),(PO(L,ISA)-S15(L,ISA))/(24.*DRT(ISA)))                              
                  HCLD(ISA)=HCL(L,ISA)                                                                                    
                  IF(KFL(1)>0)WRITE(KW(1),3360)X1,DRT(ISA),HCL(L,ISA)                                                     
              END IF                                                                                                      
              IF(KFL(1)>0.AND.IRRS(ISA)>0)WRITE(KW(1),'(T10,A,I8)')'IRR SUPPLIED FROM RES IN SA',IRRS(ISA)                             
              IF(KFL(1)>0.AND.IRRW(ISA)>0)WRITE(KW(1),'(T10,A,I8)')'IRR SUPPLIED FROM WELL IN SA',IRRW(ISA)                             
              HU(1,ISA)=0.                                                                                                
              IBGN=1                                                                                                      
              JJK=1                                                                                                       
              IPL=0                                                                                                       
              I=1                                                                                                         
              IY1=1                                                                                                       
              IF(KFL(1)>0)WRITE(KW(1),'(/1X,A)')'-----OPERATION SCHEDULE'
              DO !OPSC LOOP
                  NCRP=IGO(ISA)                                                                                         
                  K=1                                                                                                         
                  IF(KFL(1)>0)THEN
                      WRITE(KW(1),'(/T10,A,I2)')'YR',I                                                                            
                      WRITE(KW(1),316)                                                                                            
                  END IF    
                  I1=I-1                                                                                                
                  FNMX(I,ISA)=FMX                                                                                             
                  KI=0                                                                                                        
                  KF=0                                                                                                        
                  KP=0                                                                                                        
                  HU0=0.                                                                                                      
                  IF(JDHU<366)HU(JJK,ISA)=0.                                                                                  
                  IF(IGO(ISA)>0)THEN                                                                                          
                      DO J=1,MNC                                                                                              
                          IF(NHU(J,ISA)==0)CYCLE                                                                              
                          LY(I,K,ISA)=NHU(J,ISA)                                                                              
                          K=K+1                                                                                               
                      END DO                                                                                                  
                  END IF                                                                                                      
                  J=0                                                                                                         
                  DO                                                                                                          
                      J=J+1                                                                                                   
                      ! READ OPERATION SCHEDULE                                                                               
                      !  1  JX(1)= YR OF OPERATION                                                                            
                      !  2  JX(2)= MO OF OPERATION                                                                            
                      !  3  JX(3)= DAY OF OPERATION                                                                           
                      !  4  JX(4)= EQUIPMENT ID NO                                                                            
                      !  5  JX(5)= TRACTOR ID NO                                                                              
                      !  6  JX(6)= CROP ID NO                                                                                 
                      !  7  JX(7)= XMTU--TIME FROM PLANTING TO MATURITY (Y)(FOR TREE                                          
                      !            CROPS AT PLANTING ONLY).                                                                   
                      !          = LYR--TIME FROM PLANTING TO HARVEST (Y)(HARVEST ONLY)                                       
                      !          = PESTICIDE ID NO (FOR PESTICIDE APPLICATION ONLY)                                           
                      !          = FERTILIZER ID NO (FOR FERTILIZER APPLICATION ONLY)                                         
                      !  8  OPV1 = POTENTIAL HEAT UNITS FOR PLANTING (0 IF UNKNOWN)                                       
                      !          = APPLICATION VOLUME (mm)FOR IRRIGATION                                                      
                      !          = FERTILIZER APPLICATION RATE (kg/ha) = 0 FOR VARIABLE RATE                                  
                      !          = PESTICIDE APPLICATION RATE (kg/ha)                                                         
                      !          = STOCKING RATE FOR GRAZE (ha/hd)                                                            
                      !  9  OPV2 = 2 CONDITION SCS RUNOFF CURVE NUMBER (OPTIONAL)                                             
                      !          = PEST CONTROL FACTOR FOR PEST APPLICATION (FRACTION OF PESTS                                
                      !            CONTROLLED)                                                                                
                      ! 10  OPV3 = PLANT WATER STRESS FACTOR(0-1); SOIL WATER TENSION(>1KPA);                                 
                      !            OR PLANT AVAILABLE WATER DEFICIT IN ROOT ZONE(-MM)TO                                       
                      !            TRIGGER AUTO IRR. (0. OR BLANK DOES NOT CHANGE TRIGGER)                                    
                      ! 11  OPV4 = RUNOFF VOL/VOL IRRIGATION WATER APPLIED                                                    
                      ! 12  OPV5 = PLANT POPULATION (PLANTS/M**2 OR PLANTS/ha IF P/m2<1.)                                     
                      !            (FOR PLANTING ONLY)                                                                        
                      !          = FACTOR TO ADJUST AUTO IRRIGATION VOLUME (FIRG*FC)                                          
                      ! 13  OPV6 = MAX ANNUAL N FERTILIZER APPLIED TO A CROP (0. OR BLANK                                     
                      !            DOES NOT CHANGE FMX; > 0 SETS NEW FMX)(FOR PLANTING ONLY)                                  
                      ! 14  OPV7 = TIME OF OPERATION AS FRACTION OF GROWING SEASON (ENTER                                     
                      !            EARLIEST POSSIBLE MO & DAY -- JX(2) & JX(3))                                               
                      !     LINE 3/L                                                                                          
                      IF(I==1.OR.J>1)READ(KR(16),2470,IOSTAT=NFL)JX,(OPV(L),L=1,7)                                            
                      IF(NFL/=0)EXIT                                                                                          
                      IF(JX(1)/=IY1)EXIT                                                                                      
                      CALL TILTBL                                                                                             
                      LT(I,J,ISA)=NDT                                                                                         
                      JH(I,J,ISA)=JX(6)                                                                                       
                      IJ=LT(I,J,ISA)                                                                                          
                      CALL ADAJ(NC,ITL(I,J,ISA),JX(2),JX(3),NYD)                                                              
                      X4=TLD(IJ)*1000.                                                                                        
                      I3=IHC(IJ)                                                                                              
                      IF(IBGN<ITL(I,J,ISA))THEN
                          IF(IGO(ISA)>0)HU(JJK,ISA)=HU(JJK,ISA)+CAHU(IBGN,ITL(I,J,ISA),&                                    
                          BASE,0)/(PHU(JJK,IHU(JJK,ISA),ISA)+1.)                                                                  
                          HU0=HU0+CAHU(IBGN,ITL(I,J,ISA),0.,1)/AHSM                                                               
                          IBGN=ITL(I,J,ISA)                                                                                       
                      ELSE
                          IF(IBGN/=ITL(I,J,ISA))THEN                                                                        
                              IF(IGO(ISA)>0)THEN
                                  HU(JJK,ISA)=HU(JJK,ISA)+CAHU(IBGN,365,BASE,0)/&                                           
                                  (PHU(JJK,IHU(JJK,ISA),ISA)+1.)                                                                          
                                  IBGN=1                                                                                                  
                                  HU(JJK,ISA)=HU(JJK,ISA)+CAHU(IBGN,ITL(I,J,ISA),&
                                  BASE,0)/(PHU(JJK,IHU(JJK,ISA),ISA)+1.)
                                  HU0=HU0+CAHU(IBGN,ITL(I,J,ISA),0.,1)/AHSM
                                  IBGN=ITL(I,J,ISA)
                              END IF
                          END IF
                      END IF    
                      IF(OPV(7)>0..OR.IHUS==0)THEN                                                                
                          HUSC(I,J,ISA)=OPV(7)                                                                              
                      ELSE                       
                          IF(IGO(ISA)==0)THEN 
                              HUSC(I,J,ISA)=HU0                                                                                 
                          ELSE                            
                              IF(IDC(JJK)==NDC(1).OR.IDC(JJK)==NDC(2).OR.&
                              IDC(JJK)==NDC(4).OR.IDC(JJK)==NDC(5).OR.&
                              IDC(JJK)==NDC(9))HUSC(I,J,ISA)=HU(JJK,ISA)                                                                         
                          END IF
                      END IF
                      CALL INIFP(I3,I,J,JRT)                                                                            
                      ! PRINTOUT OPERATION SCHEDULE                                                                       
                      IF(KFL(1)>0)WRITE(KW(1),317)I,JX(2),JX(3),TIL(IJ),I3,COTL(IJ),&                                         
                      EMX(IJ),RR(IJ),X4,FRCP(IJ),RHT(IJ),RIN(IJ),DKH(IJ),DKI(IJ),HE&                                          
                      (IJ),ORHI(IJ),CN2(ISA),BIR(ISA),EFI(ISA),HUSC(I,J,ISA)                                                  
                      IF(TLD(IJ)>BIG(ISA))BIG(ISA)=TLD(IJ)                                                                    
                      IF(I3==NHC(5).OR.I3==NHC(6))THEN                                                                        
                          NCRP=NCRP+1                                                                                         
                          IGO(ISA)=IGO(ISA)+1                                                                                 
                          X3=OPV(5)                                                                                           
                          CALL CPTBL                                                                                          
                          IF(OPV(6)>0.)FMX=OPV(6)                                                                             
                          FNMX(I,ISA)=FMX                                                                                     
                          BASE=TBSC(JJK)                                                                                      
                          IPL=ITL(I,J,ISA)+365*I1                                                                             
                          LY(I,K,ISA)=JJK                                                                                     
                          NHU(JJK,ISA)=JJK                                                                                    
                          K=K+1                                                                                               
                          IF(KFL(1)>0)WRITE(KW(1),2480)CPNM(JJK)                                                              
                          CYCLE                                                                                               
                      END IF                                                                                                  
                      IF(I3/=NHC(1).AND.I3/=NHC(2).AND.I3/=NHC(3))CYCLE                                                       
                      IF(KDC1(JX(6))==0)CYCLE                                                                                 
                      JJK=KDC1(JX(6))                                                                                         
                      IF(I3==NHC(1))THEN                                                                                      
                          IHV=ITL(I,J,ISA)+365*I1                                                                             
                          NHU(JJK,ISA)=0                                                                                      
                          IGO(ISA)=MAX(0,IGO(ISA)-1)                                                                          
                      END IF                                                                                                  
                      HU(JJK,ISA)=0.                                                                                          
                      IF(IDC(JJK)==NDC(7).OR.IDC(JJK)==NDC(8).OR.IDC(JJK)==NDC(10))LYR&                                       
                      (I,J,ISA)=MAX(1,JX(7))                                                                                  
                      IF(KFL(1)>0)WRITE(KW(1),2480)CPNM(JJK)                                                                  
                  END DO                                                                                                      
                  ITL(I,J,ISA)=367                                                                                            
                  HUSC(I,J,ISA)=10.                                                                                           
                  NCP(I,ISA)=MAX(NCRP,1)
                  NTP=0
                  N2=0
                  DO J3=1,NCP(I,ISA)
                      J4=LY(I,J3,ISA)
                      IF(J4==0)CYCLE
                      IF(NTP(J4)>0)N2=N2+1      
                      NTP(J4)=1
                  END DO
                  NCP(I,ISA)=NCP(I,ISA)-N2
                  NTL(I,ISA)=J-1                                                                                              
                  NPST(I,ISA)=KP                                                                                              
                  NFRT(I,ISA)=KF                                                                                              
                  LT(I,J,ISA)=1                                                                                               
                  JH(I,J,ISA)=0                                                                                               
                  CND(I,J,ISA)=CN2(ISA)                                                                                       
                  QIR(I,J,ISA )=EFI(ISA)                                                                                      
                  TIR(I,J,ISA)=BIR(ISA)                                                                                       
                  FIRX(I,J,ISA)=FIRG(ISA)                                                                                     
                  RSTK(I,J,ISA)=RST0(ISA)                                                                                     
                  IF(KFL(1)>0)THEN                                                                                            
                      IF(KI>0)THEN                                                                                            
                          WRITE(KW(1),2590)                                                                                   
                          DO L=1,KI                                                                                           
                              N1=KIR(L)                                                                                       
                              ! PRINTOUT IRRIGATION SCHEDULE                                                              
                              WRITE(KW(1),2420)I,NIR(L),IIR(L),VIRR(I,N1,ISA),HUSC(I,N1,ISA)                                  
                          END DO                                                                                              
                      END IF                                                                                                  
                      IF(KF>0)THEN                                                                                            
                          WRITE(KW(1),328)                                                                                    
                          MO=1                                                                                                
                          JJ=367                                                                                              
                          KK=0                                                                                                
                          DO L=1,NTL(I,ISA)                                                                                   
                              J2=LT(I,L,ISA)                                                                                  
                              IF(IHC(J2)/=NHC(9))CYCLE                                                                        
                              X1=MAX(0.,COTL(J2))                                                                             
                              JDA=ITL(I,L,ISA)                                                                                
                              IF(NYD==0)JDA=JDA+1                                                                             
                              IF(JDA==JJ.AND.NBT(J2)==0.AND.NBE(J2)==KK)X1=0.                                                 
                              JJ=JDA                                                                                          
                              KK=NBE(J2)                                                                                      
                              CALL AXMON                                                                                      
                              CALL AICL                                                                                       
                              M=LFT(I,L,ISA)                                                                                  
                              XZ=FCST(M)*WFA(I,L,ISA)                                                                         
                              XY=X1+XZ                                                                                        
                              X4=TLD(J2)*1000.                                                                                
                              ! PRINTOUT FERTILIZER SCHEDULE                                                              
                              WRITE(KW(1),329)I,MO,KDA,FTNM(M),KDF(M),NBE(J2),NBT(J2),XY,WFA&                                 
                              (I,L,ISA),X4,FN(M),FNMA(M),FNO(M),FP(M),FPO(M),HUSC(I,L,ISA)                                    
                          END DO                                                                                              
                      END IF                                                                                                  
                      IF(KP>0)THEN                                                                                            
                          IF(MASP==1)THEN                                                                                     
                              AUNT='(g/ha)'                                                                                   
                          ELSE                                                                                                
                              AUNT='(kg/ha)'                                                                                  
                          END IF                                                                                              
                          WRITE(KW(1),3850)AUNT                                                                               
                          DO L=1,KP                                                                                           
                              N1=NPC(L)                                                                                       
                              M=LPC(I,L,ISA)                                                                                  
                              ! PRINTOUT PESTICIDE SCHEDULE                                                               
                              WRITE(KW(1),3860)I,KPC(L),JPC(L),PSTN(M),PSTR(I,L,ISA),PSTE&                                    
                              (I,L,ISA),PCST(M),HUSC(I,N1,ISA)                                                                
                          END DO                                                                                              
                      END IF                                                                                                  
                  END IF                                                                                                      
                  IF(NFL==0.AND.JX(1)>0)THEN                                                                                  
                      I=JX(1)                                                                                                 
                      IY1=I                                                                                                   
                  ELSE
                      EXIT
                  END IF    
              END DO !OPSC LOOP                                                                                                      
              NRO(ISA)=IY1                                                                                                
              JX(4)=IAUI(ISA)                                                                                             
              JX(5)=0                                                                                                     
              CALL TILTBL                                                                                                 
              IAUI(ISA)=NDT                                                                                               
              IF(KFL(1)>0)WRITE(KW(1),677)TIL(NDT),TLD(NDT)                                                               
              JX(4)=IAUF(ISA)                                                                                             
              JX(5)=12                                                                                                    
              CALL TILTBL                                                                                                 
              IAUF(ISA)=NDT                                                                                               
              JX(4)=IAMF(ISA)                                                                                             
              JX(5)=0                                                                                                     
              CALL TILTBL                                                                                                 
              IAMF(ISA)=NDT                                                                                               
              JX(4)=ISPF(ISA)                                                                                             
              JX(5)=12                                                                                                    
              CALL TILTBL                                                                                                 
              ISPF(ISA)=NDT                                                                                               
              JX(4)=ILQF(ISA)                                                                                             
              JX(5)=12                                                                                                    
              CALL TILTBL                                                                                                 
              ILQF(ISA)=NDT                                                                                               
              JX(4)=IAUL(ISA)                                                                                             
              JX(5)=0                                                                                                     
              CALL TILTBL                                                                                                 
              IAUL(ISA)=NDT                                                                                               
              K=1                                                                                                         
              IF(IDF0(1,ISA)==0)IDF0(1,ISA)=69                                                                            
              IF(IDF0(2,ISA)==0)IDF0(2,ISA)=69                                                                            
              IF(IDF0(3,ISA)==0)IDF0(3,ISA)=53                                                                            
              IF(IDF0(4,ISA)==0)IDF0(4,ISA)=52                                                                            
              IF(IDF0(5,ISA)==0.AND.FNP(5,ISA)>0.)IDF0(5,ISA)=69                                                          
              IF(IDF0(6,ISA)==0)IDF0(6,ISA)=54                                                                            
              DO K2=1,6                                                                                                   
                  IF(IDF0(K2,ISA)>0)THEN                                                                                  
                      JX(7)=IDF0(K2,ISA)                                                                                  
                      CALL NFTBL(K)                                                                                       
                      IDF0(K2,ISA)=K                                                                                      
                  END IF                                                                                                  
              END DO                                                                                                      
              IF(KFL(1)>0)THEN                                                                                            
                  WRITE(KW(1),716)TIL(IAUF(ISA)),TLD(IAUF(ISA)),FTNM(IDF0(4,ISA))                                         
                  WRITE(KW(1),689)TIL(IAMF(ISA)),TLD(IAMF(ISA)),FTNM(IDF0(2,ISA))                                         
                  WRITE(KW(1),690)TIL(ISPF(ISA)),TLD(ISPF(ISA)),FTNM(IDF0(2,ISA))                                         
                  WRITE(KW(1),691)TIL(ILQF(ISA)),TLD(ILQF(ISA)),FTNM(IDF0(1,ISA))                                         
              END IF                                                                                                      
              DO I=1,NRO(ISA)                                                                                             
                  I1=I-1                                                                                                  
                  IF(I1==0)I1=NRO(ISA)                                                                                    
                  ANA(I,ISA)=0.                                                                                           
                  IF(LY(I,1,ISA)>0)CYCLE                                                                                  
                  IF(LY(I1,NCP(I1,ISA),ISA)>0)THEN                                                                        
                      LY(I,1,ISA)=LY(I1,NCP(I1,ISA),ISA)                                                                  
                  ELSE                                                                                                    
                      I2=I1-1                                                                                             
                      LY(I,1,ISA)=LY(I2,NCP(I2,ISA),ISA)                                                                  
                  END IF                                                                                                  
              END DO                                                                                                      
              IF(IGO(ISA)>0)THEN                                                                                          
                  DO J=1,LC                                                                                               
                      IF(NHU(J,ISA)==0)CYCLE                                                                              
                      DO I=1,NCP(1,ISA)                                                                                   
                          IF(LY(1,I,ISA)==J)EXIT                                                                          
                      END DO                                                                                              
                      IF(I<=NCP(1,ISA))CYCLE                                                                              
                      NCP(1,ISA)=NCP(1,ISA)+1                                                                             
                      L=NCP(1,ISA)                                                                                        
                      L1=1000                                                                                             
                      DO WHILE(L1>1)                                                                                      
                          L1=L-1                                                                                          
                          LY(1,L,ISA)=LY(1,L1,ISA)                                                                        
                          L=L1                                                                                            
                      END DO                                                                                              
                      LY(1,1,ISA)=NHU(J,ISA)                                                                              
                  END DO                                                                                                  
              END IF                                                                                                      
              I=NRO(ISA)                                                                                                  
              JD(ISA)=MAX(1,LY(I,NCP(I,ISA),ISA))                                                                         
              IF(RZ(ISA)>XX)RZ(ISA)=XX                                                                                    
              IF(BIG(ISA)>XX)BIG(ISA)=XX                                                                                  
              TNOR(ISA)=0.                                                                                                
              STDN(JD(ISA),ISA)=8.29*STDO(ISA)                                                                            
              STDL(JD(ISA),ISA)=.1*STDO(ISA)                                                                              
              BTN(ISA)=WSA(ISA)*(ZNO3(ISA)+ZON(ISA)+STDN(JD(ISA),ISA))                                                               
              BTNX(ISA)=BTN(ISA)                                                                                          
              BTNZ(ISA)=BTN(ISA)                                                                                          
              TBTN=TBTN+BTN(ISA)                                                                                 
              STDP(JD(ISA),ISA)=1.04*STDO(ISA)                                                                            
              BTP(ISA)=WSA(ISA)*(ZPML(ISA)+ZPMA(ISA)+ZPMS(ISA)+ZPO(ISA)+&
              ZFOP(ISA)+STDP(JD(ISA),ISA))                                                                                               
              BTPX(ISA)=BTP(ISA)                                                                                          
              BTPZ(ISA)=BTP(ISA)                                                                                          
              STDK(JD(ISA),ISA)=8.29*STDO(ISA)                                                                            
              BTK(ISA)=ZSK(ISA)+ZEK(ISA)+ZFK(ISA)+STDK(JD(ISA),ISA)                                                       
              KK=NTL(1,ISA)                                                                                               
              KP1(ISA)=1                                                                                                  
              DO KT2=1,KK                                                                                                 
                  IF(ITL(1,KT2,ISA)>=IBD)EXIT
                  II=IHC(LT(1,KT2,ISA))                                                                                   
                  IF(II==NHC(7))KP1(ISA)=KP1(ISA)+1                                                                       
              END DO
              IF(KT2>KK)KT2=KK
              KT(ISA)=KT2                                                                                           
              JT2=LT(1,KT(ISA),ISA)                                                                                       
              UB1(ISA)=RZ(ISA)*5.                                                                                         
              UOB(ISA)=1.-EXP(-UB1(ISA))                                                                                  
              IGO(ISA)=0                                                                                                  
              JJK=LY(1,1,ISA)                                                                                             
              AWC(JJK,ISA)=0.                                                                                             
              IRO(ISA)=0                                                                                                  
              JP(JJK,ISA)=0                                                                                               
              KC(ISA)=0                                                                                                   
              MO=MO1                                                                                                      
              CLOSE(KR(16))                                                                                               
              DO K=1,6                                                                                                    
                  IDFT(K,ISA)=IDF0(K,ISA)                                                                                 
              END DO                                                                                                         
          END DO   ! SUBAREA LOOP
          XSL=333.6*(XCU-XCS)*COS1
          YSL=333.6*(YCU-YCS)
          XCS=111.2*(2.*XCS-XCU)*COS1
          YCS=111.2*(2.*YCS-YCU)
          SEV=0.
          !DO I=1,NWTH
              !L=I+KND
	          !CALL WREAD(L,2)
	      !END DO
  	      IF(NGN<0.OR.INFL==4)IRF=1
          TWMB=.001*TWMB
          DO IOW=1,NBON
              OSAA(IOW)=0.
              NSAS(IOW)=0
              NSAL(IOW)=0
              OWSA(IOW)=0.
              DO ISA=1,MSA
                  IF(IDON(ISA)/=IOW)CYCLE
                  IF(IAPL(ISA)==0)CYCLE
                  OWSA(IOW)=OWSA(IOW)+WSA(ISA)
                  IF(IAPL(ISA)<0)THEN
                      I1=-1
                      NSAL(IOW)=NSAL(IOW)+1
                      IDSL(NSAL(IOW),IOW)=ISA
                  ELSE
                      OSAA(IOW)=OSAA(IOW)+WSA(ISA)
                      NSAS(IOW)=NSAS(IOW)+1
                      IDSS(NSAS(IOW),IOW)=ISA      
                      I1=1
                  END IF
                  IF(NFED(IOW)==0)CYCLE
                  DO J=1,NFED(IOW)
                      IF(IDFD(J,IOW)==I1*IAPL(ISA))EXIT
                  END DO
                  IAPL(ISA)=I1*IDFA(J,IOW)
              END DO
              ISAS(IOW)=IDSS(1,IOW)
          END DO
          NY=0
          DO IOW=1,NBON
              NGIX(1,1,IOW)=1
              IF(NHRD(IOW)==0)CYCLE
              DO IHD=1,NHRD(IOW)
                  JX(7)=IDMU(IHD,IOW)
                  CALL NFTBL(K)
                  IDMU(IHD,IOW)=K
                  IGZO(IHD,IOW)=0
                  DO ISA=1,MSA
                      IF(IDON(ISA)/=IOW)CYCLE
                      IF(IFED(IHD,IOW)==ISA)THEN
                          GCOW(IHD,ISA)=NCOW(IHD,IOW)*FFED(IHD,IOW)
                      ELSE
                          IF(NGZ(IHD,ISA)==0)CYCLE
                          GCOW(IHD,ISA)=NCOW(IHD,IOW)*(1.-FFED(IHD,IOW))
                          IGZO(IHD,IOW)=IGZO(IHD,IOW)+1
                          NGIX(IGZO(IHD,IOW),IHD,IOW)=ISA
                          IGZ(ISA)=0
                      END IF
                  END DO
                  ISA=MSA
                  NGZA(IHD,IOW)=IGZO(IHD,IOW)
                  I2=NGZA(IHD,IOW)
                  IF(I2==0)CYCLE
                  IGZX(IHD,IOW)=NGIX(I2,IHD,IOW)
                  IF(IHD>1)THEN
                      DO J=1,IHD
                          J1=IGZO(IHD,IOW)
                          DO WHILE(IGZX(IHD,IOW)==NY(J))
                              J1=J1-1      
                              IGZX(IHD,IOW)=NGIX(J1,IHD,IOW)
                          END DO
                          IGZO(IHD,IOW)=J1
                      END DO
                  END IF
                  NY(IHD)=IGZX(IHD,IOW)
              END DO
              IHD=NHRD(IOW)
          END DO
          IF(NDP>0)THEN
	          IF(KFL(1)/=0)THEN
                  WRITE(KW(1),'(//1X,A/)')'____________________PESTICIDE &
                  DATA____________________'
                  WRITE(KW(1),3900)
              END IF
              DO I=1,NDP
                  ! PRINTOUT PESTICIDE DATA
                  IF(KFL(1)>0)WRITE(KW(1),3880)PSTN(I),PSOL(I),PHLS(I),PHLF(I)&
                  ,PWOF(I),PKOC(I),PCST(I)
                  PSOL(I)=PSOL(I)*10.
                  PHLS(I)=1.-EXP(-.693/PHLS(I))
                  PHLF(I)=1.-EXP(-.693/PHLF(I))
              END DO
          END IF
          CALL CRPIO
          DO I=2,MSO-1
              IF(KFL(I)==0.OR.(I>33.AND.I<37.OR.I==49))CYCLE
              WRITE(KW(I),3381)IYER,IMON,IDAY,IT1,IT2,IT3
              WRITE(KW(I),'(T10,A,I4)')'RUN # ',IRUN
              WRITE(KW(I),'(10X,20A4)')TITLE
              WRITE(KW(I),5031)SITEFILE
              WRITE(KW(I),5031)SAFILE
          END DO
          IF(KFL(45)>0)THEN
              WRITE(KW(45),'(3A6,30A10)') 'YEAR','DAY','SUB','PRECIP',&
              'IRR','PKRZ','ET','Q','LAI','SEDYLD','QN','QP','PAD_QVOL',&
              'PAD_SED','PAD_MinN','PAD_MinP'
              WRITE(KW(45),'(3A6,30A10)') ' ',' ',' ','(mm)','(mm)',&
              '(mm)','(mm)','(mm)',' ','(t/ha)','(kg/ha)','(kg/ha)',&
              '(mm)','(t/ha)','(kg/ha)','(kg/ha)'
          END IF    
          CALL ASORTI(NBSA,IBSA,MSA)
          NPSO=1
          L=2
          MXW=1
          LND=KND+NWTH
          DO ISA=1,MSA
              IF(IPTS(ISA)>0)THEN
              ! DAILY POINT SOURCE FILE NAME      
                  JJ=-1
                  DO WHILE(JJ/=IPTS(ISA))
                      READ(KR(27),*,IOSTAT=NFL)JJ,FPSO(JJ)
                      IF(NFL/=0)THEN
                          WRITE(*,*)'POINT SOURCE NO = ',IPTS(ISA),' NOT IN PT SO &
                          LIST FILE     SAID = ',NBSA(ISA)
	                      STOP
	                  END IF
	              END DO
	              REWIND KR(27)
                  IF(NPSO>1)THEN
                      DO L=1,MXW
                          IF(NBW(L)==IPTS(ISA))EXIT
                      END DO
                  END IF
                  IF(L>MXW)THEN
                      MXW=NPSO
                      L=MXW+LND
                      NBW(MXW)=IPTS(ISA)
                      IPSO(ISA)=NPSO
	                  NPSO=NPSO+1
	                  CALL OPENV(KR(L),FPSO(JJ),IDIR)
	                  DO J=1,6
			              READ(KR(L),'()')
		              END DO
	                  CALL WREAD(L,3)
	              ELSE
	                  IPSO(ISA)=L  
	              END IF
              END IF	
              LRG=0
        !     IF(KFL(21)>0)CALL SOCIOA(1)
              IF(KFL(22)>0)CALL SOCIOD(1)
              DO J=1,LC
                  STD(J,ISA)=0.
                  NHU(J,ISA)=IHU(J,ISA)
                  IF(NHU(J,ISA)>LRG)LRG=NHU(J,ISA)
                  IF(RDMX(J)>RZ(ISA))RZ(ISA)=RDMX(J)
                  IHU(J,ISA)=1
                  UNA(J,ISA)=625.*BN(3,J)*WA(J)
              END DO
              YLTI=YCT(ISA)/CLT
              SINI=SIN(YLTI)
              COSI=COS(YLTI)
              DO J=ISA+1,MSA
                  RY=YCT(J)/CLT
                  XX=SINI*SIN(RY)+COSI*COS(RY)*COS((XCT(J)-XCT(ISA))/CLT)
                  XX=MIN(1.,XX)
                  SADST(ISA,J)=6378.8*ACOS(XX)
                  !SADST(ISA,J)=ATRI(.25,.5,1.,9)
                  SADST(J,ISA)=SADST(ISA,J)
              END DO
          END DO    
          NPSO=NPSO-1
          SMMC=0.
          IF(KFL(1)>0)THEN
              WRITE(KW(1),2122)'POP',(POPX(J),J=1,LC)
              WRITE(KW(1),2122)'MXLA',(PLAX(J),J=1,LC)
              WRITE(KW(1),2123)'PHU',(PHUX(J),J=1,LC)
              CALL APAGE(0)
              WRITE(KW(1),'(//1X,A/)')'________________SUBAREA HYDROLOGIC DATA__&
              _______________________'
              WRITE(KW(1),2111)
              SMWD=0.
              DO ISA=1,MSA
                  AWTH=' '
                  IF(NGN>0)THEN
                      AWTH=FWTH(IRF(ISA))
                      SMWD(IRF(ISA))=SMWD(IRF(ISA))+WSA(ISA)
                  END IF
                  IF(ISLF==0)THEN
	                  X1=RLF(ISA)*RSF(ISA)
	              ELSE
	                  X1=SLF(ISA)
	              END IF
                  WRITE(KW(1),3121)ISA,NBSA(ISA),ACET(1,ISA),XCT(ISA),&
                  YCT(ISA),(ACET(I,ISA),I=2,8),X1,TC(ISA),RFPK(ISA),&
                  GWST(ISA),GWMX(ISA),RFTT(ISA),SNO(ISA),STDO(ISA),&
                  URBF(ISA),BCOF(ISA),BFFL(ISA),CN2(ISA),TITSO(ISA),&
                  TITOP(ISA),NVCN(ISA),SATK(ISA),AWTH
	          END DO
	          WRITE(KW(1),'(T5,A,3E20.10)')'SUM WSA = ',TOT,SUMA,RWSA(MSA)
	          IF(NGN>0)THEN
                  DO I=1,NWTH
                      SMWD(I)=SMWD(I)/TOT
                      WRITE(KW(1),2124)FWTH(I),SMWD(I)
                  END DO
              END IF
              CALL APAGE(0)
              WRITE(KW(1),'(//1X,A/)')'___________________ROUTING REACH DATA____&
              _______________________'
              WRITE(KW(1),2126)
              DO ISA=1,MSA
                  ! PRINTOUT REACH DATA
                  IF(ABS(ACET(2,ISA)-RCHL(ISA))>1.E-5)WRITE(KW(1),3122)ISA,NBSA&
                  (ISA),WSA(ISA),RWSA(ISA),RCHL(ISA),RCHD(ISA),RCHS(ISA),RCBW(ISA)&
                  ,RCTW(ISA),RCHN(ISA),RCHC(ISA),RCHK(ISA),RFPW(ISA),RFPL(ISA),&
                  RFPS(ISA)
                  RCHK(ISA)=RCHK(ISA)*RCHC(ISA)
              END DO
              CALL APAGE(0)
              WRITE(KW(1),'(//1X,A/)')'___________________RESERVOIR DATA________&
              ___________________'
              WRITE(KW(1),2128)
          END IF    
	      SRYB=0.
	      SRVA=0.
	      SPDA=0.
          BSWB=0.
          BSNOB=0.
          BRSVB=0.
          BGWSB=0.
          BSWLB=0.
	      JDA=IBD-1
          IF(JDA<=0)JDA=365
          DO ISA=1,MSA !INITIAL WATER CONTENT LOOP
              X1=RWSA(ISA)
              IF(IEXT(ISA)>0)X1=WSA(ISA)
              RFTT(ISA)=1.-EXP(-1./RFTT(ISA))
              IF(ISA==KSA)THEN
                  GWE0(ISA)=ELMN-50.
                  GWEL(ISA)=GWE0(ISA)
              ELSE    
                  GWE0(ISA)=ELMN-50.+GWSP*SADST(ISA,KSA)
                  GWEL(ISA)=GWE0(ISA)
              END IF    
              WSAX1=10.*WSA(ISA)
              GWST(ISA)=GWST(ISA)*WSAX1
              GWMX(ISA)=GWMX(ISA)*WSAX1
              BSW(ISA)=WSAX1*SW(ISA)
              BSWB=BSWB+BSW(ISA)
              BSNO(ISA)=WSAX1*SNO(ISA)
              BSNOB=BSNOB+BSNO(ISA)
              BGWS(ISA)=GWST(ISA)
              BGWSB=BGWSB+BGWS(ISA)
              SAET=0.
              CALL WHLRMX(JDA)
              HR0(ISA)=HRLT
              SALA(ISA)=WSA(ISA)
              IF(RSAE(ISA)>0.)THEN
                  IF(PCOF(ISA)>0..AND.PCOF(ISA)<1.)THEN
                      X1=WSA(ISA)*PCOF(ISA)
	                  SPDA=SPDA+X1
	              ELSE
                      SRVA=SRVA+X1
	              END IF
                  ! PRINTOUT RESERVOIR DATA
                  IF(KFL(1)>0)WRITE(KW(1),3123)ISA,NBSA(ISA),RWSA(ISA),RSAE(ISA),RVE0(ISA),RSAP&
                  (ISA),RVP0(ISA),RSRR(ISA),RSV(ISA),RSYS(ISA),RSYN(ISA),RSHC(ISA),&
                  RSDP(ISA)
                  A10=10.*X1
                  ! SET INITIAL RESERVOIR PARAMETERS
                  RSYS(ISA)=RSYS(ISA)*1.E-6
                  RSYN(ISA)=RSYN(ISA)*1.E-6
	              RVE0(ISA)=RVE0(ISA)*A10
                  RVP0(ISA)=RVP0(ISA)*A10
	              RSVP(ISA)=RVP0(ISA)
	              RSVE(ISA)=RVE0(ISA)
                  RSV(ISA)=RSV(ISA)*A10
                  RSVB(ISA)=RSV(ISA)
            !     SWST(LID(1,ISA),ISA)=SWST(LID(1,ISA),ISA)+.1*RSV(ISA)/WSA(ISA)
                  RSRR(ISA)=(RSVE(ISA)-RSVP(ISA))/RSRR(ISA)
                  RSDP(ISA)=EXP(-4.605/RSDP(ISA))
                  RSO3(ISA)=0.
                  RSSP(ISA)=.0001*RSV(ISA)
                  X1=LOG10(RSAE(ISA))-LOG10(RSAP(ISA))
	              X2=LOG10(RSVE(ISA))-LOG10(RSVP(ISA))
                  !BV2(ISA)=X2/X1 
                  !BV1(ISA)=RSVE(ISA)/RSAE(ISA)**BV2(ISA)
	              !BA2(ISA)=(LOG10(RSEE(ISA))-LOG10(RSEP(ISA)))/X1
                  !BA1(ISA)=RSEE(ISA)/RSAE(ISA)**BA2(ISA)
                  BR2(ISA)=X1/X2
                  BR1(ISA)=RSAE(ISA)/RSVE(ISA)**BR2(ISA)
                  RSSA(ISA)=BR1(ISA)*RSV(ISA)**BR2(ISA)
                  SALA(ISA)=WSA(ISA)-RSSA(ISA)
                  X1=WSAX1*BSW(ISA)
	              IF(KFL(13)>0)WRITE(KW(13),4094)ISA,X1
                  RSYS(ISA)=RSYS(ISA)*RSV(ISA)
                  SRYB=SRYB+RSYS(ISA)
                  RSON(ISA)=1.5*RSYS(ISA)
                  RSOP(ISA)=RSYS(ISA)
	              RSOC(ISA)=15.*RSYS(ISA)
                  RSYB(ISA)=RSYS(ISA)
              ELSE
                  RSV(ISA)=0.
              END IF    
              BRSV(ISA)=RSV(ISA)
              BRSVB=BRSVB+RSV(ISA)
              BSALA(ISA)=SALA(ISA)
          END DO  !INITIAL WATER CONTENT LOOP
          SRYX=SRYB
          RTO=SRVA/RWSA(MSA)
	      IF(KFL(1)>0)WRITE(KW(1),'(/T10,A,F8.4)')'FRACTION OF WATERSHED CONTROLLED BY &
	      RESERVOIRS = ',RTO
	      NCMD=K1-1
	      RWSA(NCMD)=SUMA 
          IF(KFL(1)>0)THEN
              CALL APAGE(0)
              WRITE(KW(1),'(//1X,A/)')'____________________ROUTING SCHEME&
              ______________'
              WRITE(KW(1),16)
          END IF
          ! PRINTOUT ROUTING COMMAND TABLE
          DO I=1,NCMD
              IF(ICDT(I)>2)THEN
                  IF(KFL(1)>0)WRITE(KW(1),17)RCMD(ICDT(I)),ICDT(I),IDOT(I),&
                  IDN1T(I),IDN2T(I)
              ELSE
                  IDO=IDOT(I)
                  IF(ICDT(I)==2)THEN
                      ISA=IDN2T(I)
                      I1=NBSA(ISA)
                      IF(KFL(1)>0)WRITE(KW(1),17)RCMD(ICDT(I)),ICDT(I),IDOT(I),&
                      IDN1T(I),IDN2T(I),I1
                  ELSE
                      ISA=IDN1T(I)
                      I1=NBSA(ISA)
                      IF(KFL(1)>0)WRITE(KW(1),2151)RCMD(ICDT(I)),ICDT(I),IDOT(I),&
                      IDN1T(I),I1,IDN2T(I)
                  END IF
                  IDOA(ISA)=IDO              
              END IF     
          END DO
          DO    ! OUTPUT LOOP
              IF(KFL(2)>0)WRITE(KW(2),2284)
              IF(KFL(3)>0)WRITE(KW(3),4085)(HED(KY(J1)),J1=1,NKY),HED(NDVSS)
              IF(KFL(4)>0)WRITE(KW(4),4081)(HED(KY(J1)),J1=1,NKY)
              IF(KFL(5)>0)THEN
                  WRITE(KW(5),935)SUMA,(PSTN(K),K=1,NDP)
                  WRITE(KW(5),937)                                        
              END IF
              IF(KFL(37)>0)THEN
                  WRITE(KW(37),935)SUMA
                  WRITE(KW(37),939)                                        
              END IF
              IF(KFL(8)>0)WRITE(KW(8),1001)
              IF(KFL(9)>0)WRITE(KW(9),2112)
              IF(KFL(43)>0)THEN
                  WRITE(KW(43),'(//1X,A/)')'____________________SOIL DATA___&
                  _________________'
                  WRITE(KW(43),900)
              END IF 
              IF(KFL(46)>0)THEN
                  WRITE(KW(46),'(//1X,A/)')'____________________SOIL DATA___&
                  _________________'
                  WRITE(KW(46),900)
              END IF    
              IF(KFL(47)>0)THEN
                  WRITE(KW(47),'(//1X,A/)')'____________________SOIL DATA___&
                  _________________'
                  WRITE(KW(47),900)
              END IF    
              IF(KFL(48)>0)THEN
                  WRITE(KW(48),'(//1X,A/)')'____________________SOIL DATA___&
                  _________________'
                  WRITE(KW(48),900)
              END IF    
              IF(KFL(11)>0)WRITE(KW(11),4082)HEDC,(HED(KD(J1)),J1=1,NKD),&
              (HEDS(KS(J1)),J1=1,NKS)
              !IF(KFL(11)>0)WRITE(KW(11),4082)(HED(KD(J1)),J1=1,NKD)
              !(HEDS(KS(J1)),J1=1,NKS),HEDC
              !IF(KFL(11)>0)WRITE(KW(11),4082)'GWST','GWSN',HED(16),HED(71),&
              !HED(72),HED(40),HED(96),HED(80)
              IF(KFL(44)>0)WRITE(KW(44),3990)(HED(KD(J1)),J1=1,NKD)
              IF(KFL(13)>0)WRITE(KW(13),4084)
              IF(KFL(14)>0)WRITE(KW(14),2313)
              IF(KFL(15)>0)WRITE(KW(15),603)
              IF(KFL(16)>0)WRITE(KW(16),4083)' RFV',(HED(KD(J1)),J1=1,NKD)
	          IF(KFL(29)>0)WRITE(KW(29),2072)
              IF(KFL(17)>0)WRITE(KW(17),3284)HED(4),(HEDR(J),HEDR(J+10),J=1,10)
              IF(KFL(18)>0)WRITE(KW(18),3285)
              IF(KFL(19)>0)WRITE(KW(19),3286)
              IF(KFL(20)>0)WRITE(KW(20),3287)
              IF(KFL(24)>0)WRITE(KW(24),4091)
              IF(KFL(25)>0)WRITE(KW(25),3288)
	          IF(KFL(27)>0)WRITE(KW(27),729)
	          IF(KFL(28)>0)THEN
                  WRITE(KW(28),934)SUMA,(PSTN(K),K=1,NDP)                                        
                  WRITE(KW(28),938)
              END IF
	          IF(KFL(30)>0)WRITE(KW(30),'(/1X,A)')'-----AVE ANNUAL VALUES(g/ha)'
	          IF(KFL(31)>0)WRITE(KW(31),583)
	          IF(KFL(32)>0)WRITE(KW(32),584)
	          IF(KFL(38)>0)WRITE(KW(38),585)
	          IF(IPD>0.AND.KFL(34)>0)THEN
	              WRITE(KW(34),2073)IYER,IMON,IDAY,IT1,IT2,IT3
                  WRITE(KW(34),854)ASTN
                  WRITE(KW(34),856)SITEFILE
                  WRITE(KW(34),856)SAFILE
                  WRITE(KW(34),2142)HED(4),'mm  ',HED(5),'mm  ',HED(6),'mm  ',HED&
                  (18),'mm  ',HED(10),'mm  ',HED(11),'mm  ',HEDS(8),'mm  ',HED(16),&
                  'mm  ',HED(71),'mm  ','QSUR','mm  ',HED(15),'mm  ',HED(72),'mm  ',&
                  HED(117),'mm  ',HED(14),'-   ',HED(1),'c   ',HED(2),'c   ',HED(59)&
                  ,'c   ',HED(3),'mjm2',HED(27),'t/ha',HED(31),'t/ha',HED(53),'kgha'&
                  ,HED(54),'kgha',HED(55),'kgha',HED(56),'kgha',HED(57),'kgha',HED&
                  (43),'kgha',HED(42),'kgha',HED(37),'kgha',HED(119),'kgha',HED(38),&
                  'kgha',HED(49),'kgha',HED(118),'kgha',HED(39),'kgha',HED(80),&
                  'kgha',(HEDC(1),'-   ',HEDC(2),'-   ',HEDC(3),'m   ',HEDC(4),&
                  't/ha',HEDC(5),'t/ha',HEDC(6),'t/ha',HEDC(7),'m   ',HEDC(8),&
                  't/ha',HEDC(9),'t/ha',HEDC(10),'d   ',HEDC(11),'d   ',HEDC(12),&
                  'd   ',HEDC(13),'d   ',HEDC(14),'d   ',HEDC(15),'d   ',HEDC(16),&
                  'd   ',HEDC(17),'d   ',J=1,5)
	          END IF
	          DO I1=35,39,4
                  IF(KFL(I1)>0)THEN
                      WRITE(KW(I1),2073)IYER,IMON,IDAY,IT1,IT2,IT3
                      WRITE(KW(I1),854)ASTN
                      WRITE(KW(I1),856)SITEFILE
                      WRITE(KW(I1),856)SAFILE
                      WRITE(KW(I1),2162)HEDH(1),'m3/s',HEDH(2),'m3/s',HEDH(33),'m3/s',&
                      HEDH(34),'m3/s',HEDH(3),'m3/s',HEDH(4),'m3/s',HEDH(5),'t   ',&
                      HEDH(6),'t   ',HEDH(7),'ppm ',HEDH(8),'kg  ',HEDH(9),'kg  ',&
                      HEDH(10),'kg  ',HEDH(11),'kg  ',HEDH(12),'kg  ',HEDH(13),'kg  ',&
                      HEDH(14),'kg  ',HEDH(15),'kg  ',HEDH(16),'kg  ',HEDH(17),'kg  ',&
                      HEDH(18),'kg  ',HEDH(19),'kg  ',HEDH(20),'-   ',HEDH(21),'-   ',&
                      HEDH(22),'-   ',HEDH(23),'-   ',HEDH(24),'-   ',HEDH(25),'-   ',&
                      HEDH(26),'g   ',HEDH(27),'g   ',HEDH(28),'g   ',HEDH(29),'g   ',&
                      HEDH(30),'g   ',HEDH(31),'g   ',HEDH(32),'g   '
                  END IF
              END DO
              IF(IPD>0.AND.KFL(40)>0)THEN
	              WRITE(KW(40),2073)IYER,IMON,IDAY,IT1,IT2,IT3
                  WRITE(KW(40),854)ASTN
                  WRITE(KW(40),856)SITEFILE
                  WRITE(KW(40),856)SAFILE
                  WRITE(KW(40),2129)HED(4),'mm  ',HED(13),'mm  ',HED(117),'mm  ',&
                  HED(16),'mm  ','RSTKha/au ',('CPNM-   ',HEDC(1),'-   ',HEDC(2),&
                  '-   ',HEDC(5),'t/ha',HEDC(6),'t/ha',HEDC(8),'t/ha',HEDC(10),&
                  'd   ',HEDC(11),'d   ',HEDC(12),'d   ',HEDC(13),'d   ',HEDC(14),&
                  'd   ',HEDC(15),'d   ','CNSLg/g','CNSDg/g','YLDkg/ha','GRZDd',J=1,5)
              END IF
              IF(KFL(7)>0)WRITE(KW(7),3282)
              IF(KFL(6)>0)WRITE(KW(6),3283)HEDP
	          IHBS=1
              
              !rtb salt
              if(ISALT>0) call salt_read
            
              ! BEGIN ANNUAL SIMULATION LOOP
              CALL BSIM
              ADHY=360.*ADHY*DTHY/SUMA
              IF(KFL(12)>0)WRITE(KW(12),'(25X,A,F10.2,A)')'OUTFLOW = ',ADHY,' mm'
              SM1=0.
              SM2=0.
              SM3=0.
              SM4=0.
              SM6=0.
              SM7=0.
              SM8=0.
              SM9=0.
	          SM11=0.
              SM18=0.
              SM19=0.
              SM40=0.
              SM42=0.
              SM43=0.
              SM46=0.
              SM53=0.
              SM54=0.
              SM55=0.
              SM66=0.
              SM68=0.
              SM69=0.
              SM71=0.
              SM81=0.
              SM82=0.
              SM95=0.
              SM96=0.
              SM104=0.
              SM105=0.
              SM106=0.
              SM110=0.
              SM134=0.
              SM146=0.
              SM147=0.
              SM155=0.
              SM156=0.
              IY=IY-1
              XYR=IY
              RR=0.
	          KSO=MSO
	          WRITE(KW(28),4993)SMSW
	          SMOC=0.
	          SRYF=0.
	          SFTN=0.
	          SRQN=0.
	          TRMN=0.
	          TRON=0.
	          SMWS=0.
	          SMNS=0.
	          SMPS=0.
	          SMKS=0.
	          SMTS=0.
	          SMAS=0.
	          SMSS=0.
	          SMPL=0.
	          SMPQ=0.
	          SMPY=0.
	          SMY1=0.
	          SMY2=0.
              FSNOB=0.
              FSWB=0.
              FRSVB=0.
              FGWSB=0.
              FSWLB=0.
              DO ISA=1,MSA !SUBAREA LOOP OUTPUT
                  IDX=IDOA(ISA)
                  WSAX=WSA(ISA)
                  WSAX1=WSAX*10.
                  IF(KFL(1)>0)CALL APAGE(1)
                  CALL SPRNT
                  DO I=1,12
                      X1=MAX(1,IHRL(I,ISA))
                      THRL(I,ISA)=THRL(I,ISA)/X1
                      SRMX(I,ISA)=SRMX(I,ISA)/X1
                  END DO
                  SW(ISA)=1000.*Z(LID(NBSL(ISA),ISA),ISA)*UW(2)
                  IF(KFL(MSO+1)>0)THEN
                      I2=NBSL(ISA)
	                  KSO=KSO+1
	                  XCC=1.
                      X1=0.
	                  X2=ISG(ISA)
	                  WRITE(ASOL,'(I8.8)')NBSA(ISA)
	                  ASOT=ASOL//".SOT"
	                  OPEN(KW(KSO),FILE=ASOT)
	                  TITSO(ISA)=ASOT
                      WRITE(KW(KSO),526)ASOT,IYER,IMON,IDAY
                      WRITE(KW(KSO),3114)SALB(ISA),X2,FFC(ISA),TSLA(ISA),&
                      XIDS(ISA),X1,XIDK(ISA),ZQT,ZF,ZTK,FBM(ISA),&
                      FHP(ISA),XCC
                      WRITE(KW(KSO),3113)(Z(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3113)(BDP(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3113)(SSF(LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3113)(SOIL(9,LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3113)(SAN(LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3113)(SIL(LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3111)(SOIL(6,LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3113)(PH(LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3113)(SMB(LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3113)(SOIL(7,LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3113)(CAC(LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3113)(CEC(LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3113)(ROK(LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3113)(SOIL(5,LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3111)(SOIL(1,LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3113)(RSD(LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3113)(SOIL(13,LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3113)(PSP(LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3113)(SATC(LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3113)(HCL(LID(I,ISA),ISA),I=1,I2) 
                      WRITE(KW(KSO),3111)(SOIL(4,LID(I,ISA),ISA),I=1,I2)
	                  WRITE(KW(KSO),3113)(SOIL(14,LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3113)(ECND(LID(I,ISA),ISA),I=1,I2)
	                  WRITE(KW(KSO),3113)(STFR(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3113)(SOIL(12,LID(I,ISA),ISA),I=1,I2)
	                  WRITE(KW(KSO),3111)(CPRV(LID(I,ISA),ISA),I=1,I2)
	                  WRITE(KW(KSO),3111)(CPRH(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3111)(WLS(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3111)(WLM(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3111)(WLSL(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3111)(WLSC(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3111)(WLMC(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3111)(WLSLC(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3111)(WLSLNC(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3111)(WBMC(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3111)(WHSC(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3111)(WHPC(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3113)(WLSN(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3113)(WLMN(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3113)(WBMN(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3111)(WHSN(LID(I,ISA),ISA),I=1,I2)
                      WRITE(KW(KSO),3111)(WHPN(LID(I,ISA),ISA),I=1,I2)
	                  WRITE(KW(KSO),3113)(SOIL(16,LID(I,ISA),ISA),I=1,I2)
	                  WRITE(KW(KSO),3113)(SOIL(17,LID(I,ISA),ISA),I=1,I2)
	              END IF	      
                  IF(KFL(1)>0)THEN
                      WRITE(KW(1),'(//1X,A/)')'____________________FINAL SOIL &
                      PHYSICAL DATA____________________'
                      CALL SOLIOP
                      CALL APAGE(1)
                      WRITE(KW(1),'(//1X,A/)')'____________________FINAL SOIL &
                      CHEMICAL DATA____________________'
                      CALL SOLIOC
                      WRITE(KW(1),'(/T10,A,F7.1,A)')'ERODED SOIL THICKNESS = ',&
                      THK(ISA),' mm'
                      WRITE(KW(1),'(/T10,A,F7.2,A)')'FINAL WATER CONTENT OF &
                      SNOW = ',SNO(ISA),' mm'
                      SRSD(ISA)=SRSD(ISA)/XYR
                      WRITE(KW(1),'(/T10,A,F6.1,A)')'AVE ANNUAL CROP RESIDUE = ',&
                      SRSD(ISA),' T'
                  END IF
                  FSNO=SNO(ISA)*WSAX1
                  FSNOB=FSNOB+FSNO
                  FSW=SW(ISA)*WSAX1
                  FSWB=FSWB+FSW
                  FGWS=GWST(ISA)
                  FGWSB=FGWSB+FGWS
                  FSWL=SWLT(ISA)*WSAX1
                  FSWLB=FSWLB+FSWL
                  IPR=2
                  IF(RSAE(ISA)>0.)THEN
                     X5=SM(65,ISA)
                     IPR=1
                  ELSE
                      IF(IEXT(ISA)>0)THEN
                          X5=SM(117,ISA)
                          X15=SM(15,ISA)
                          X17=0.
                          X2=0.
                      ELSE    
                          X5=SMH(34,IDOR(ISA))
                          X2=SMIO(IDNF(ISA))
                          X15=SMH(15,IDOR(ISA))
                          X17=SMH(15,IDNF(ISA))
                      END IF    
                  END IF    
                  IF(PCOF(ISA)>0.)IPR=0
	              IF(IPR<2)THEN
	                  CALL RESPQB(SM(64,ISA),SM(66,ISA),SM(65,ISA),SM(67,ISA),&
                      SM(111,ISA),RSVF(ISA),RSVB(ISA),SM(146,ISA),ISA,KFL,&
                      KW,NBSA,MSO,IPR)
                      CALL RESYB(SM(68,ISA),SM(69,ISA),SM(112,ISA),SM(70,ISA),&
                      RSYF(ISA),RSYB(ISA),ISA,KFL,KW,NBSA,MSO,IPR)
     	              SRYF=SRYF+RSYF(ISA)
     	              TRMN=TRMN+RSO3(ISA)
     	              TRON=TRON+RSON(ISA)
                  END IF  
                  ! PRINTOUT SA WATER BALANCE
                  IF(KFL(1)>0)CALL HSWBL(SM(4,ISA),X2,X5,X17,X15,SM(155,ISA),&
                  SM(147,ISA),SM(66,ISA),SM(71,ISA),SM(18,ISA),SM(110,ISA),&
                  SM(19,ISA),SM(67,ISA),SM(146,ISA),SM(156,ISA),SM(157,ISA),&
                  BSNO(ISA),FSNO,BSW(ISA),FSW,BRSV(ISA),RSV(ISA),&
                  BGWS(ISA),FGWS,FSWL,ISA,KW,NBSA,MSO)
	              SRYF=SRYF+RSYF(ISA) 
                  SFNT=ZNO3(ISA)+ZON(ISA)+STDN(JD(ISA),ISA)
                  IF(DALG(ISA)>0.)THEN
                      IF(KW(1)>0)CALL HLGB(SM(22,ISA),SM(20,ISA),SM(23,ISA),SM(52,ISA),VLG&
                      (ISA),VLGB(ISA),SM(21,ISA),ISA,KW,NBSA,MSO)
                      II=IFED(IDFH(ISA),IDON(ISA))
                      !CALL NLGB(SM(61,ISA),SM(62,II),SOFU,WTMB(ISA),WTMU(ISA),&
                      !ISA,KW,NBSA,MSO)
                      SM6=SM6+SM(62,II)
                  END IF
                  SM7=SM7+FCMP(ISA)*WSAX
                  SM8=SM8+FCMN(ISA)*WSAX
                  SM9=SM9+SM(4,ISA)
                  SM42=SM42+SM(42,ISA)
                  SM46=SM46+SM(46,ISA)
                  SM43=SM43+SM(43,ISA)
                  SM53=SM53+SM(53,ISA)
                  SM54=SM54+SM(54,ISA)
                  SM55=SM55+SM(55,ISA)
                  SM82=SM82+SM(82,ISA)
                  SM96=SM96+SM(96,ISA)
                  SM105=SM105+SM(105,ISA)
                  RFQN(ISA)=SM(4,ISA)*RFNC
                  SRQN=SRQN+RFQN(ISA)
                  SM134=SM134+SM(134,ISA)
                  AD1=0.
                  AD2=0.
                  AD3=0.
                  ADD=0.
                  SUM=0.
                  TOT=0.
                  DO J=1,LC
                      AD1=AD1+STDP(J,ISA)
                      AD2=AD2+STDK(J,ISA)
                      ADD=ADD+STDN(J,ISA)
                      TOT=TOT+UP1(J,ISA)
                      SUM=SUM+UN1(J,ISA)
                      AD3=AD3+UK1(J,ISA)
                  END DO
                  FTN=WSAX*(ZNO3(ISA)+ZNH3(ISA)+ZON(ISA)+ADD+STDON(ISA)+SUM+ZNOU(ISA))+GWSN(ISA) 
                  SFTN=SFTN+FTN
                  ! SA N BALANCE
                  IF(KFL(1)>0)THEN
                      IF(IEXT(ISA)>0)THEN
                          X5=SM(37,ISA)
                          X4=SM(158,ISA)
                          X7=SM(39,ISA)
                          X2=0.
                          X3=0.
                          X6=0.
                          X7=0.
                      ELSE    
                          X5=SYN(IDOR(ISA))
                          X4=SQN(IDOR(ISA))
                          X2=SYN(IDNF(ISA))
                          X3=SQN(IDNF(ISA))
                          X6=SSN(IDNF(ISA))
                          X7=SSN(IDOR(ISA))
                      END IF
                      TYN(ISA)=TYN(ISA)*WSAX
                      CALL NBL(BTN(ISA),X2,X5,X3,X4,X6,X7,SM(96,ISA),RFQN(ISA),&
                      SM(134,ISA),SM(42,ISA),SM(53,ISA),TYN(ISA),SM(46,ISA),&
                      SM(54,ISA),SM(55,ISA),SM(43,ISA),FTN,GWSN(ISA),SM(82,ISA),&
                      SM(89,ISA),SM(105,ISA),ISA,1,KW,NBSA,MSO,JRT)
                      !(BTN,QNI,QNO,YNI,YNO,DPKN,RN,YNWN,DN,TFO,YLN,VOL,FNMN,FNMA,&
                      !FX,FTN,BURN,SCOU,PSON,ISA,KBL,KFL,KW,NBSA,MSO,JRT)
	                  IF(JRT>0)WRITE(KW(1),4088)ZNO3(ISA),ZNH3(ISA),&
                      ZON(ISA),ADD,STDON(ISA),SUM,ZNMU(ISA),ZNOU(ISA)
                      ! SA C BALANCE
                      IF(IEXT(ISA)>0)THEN
                          X5=SM(77,ISA)
                          X4=SM(76,ISA)
                          X2=0.
                          X3=0.
                      ELSE    
                          X5=SYC(IDOR(ISA))
                          X4=SQC(IDOR(ISA))
                          X2=SYC(IDNF(ISA))
                          X3=SQC(IDNF(ISA))
                      END IF 
                      FTC=ZOC(ISA)*WSAX
                      CALL NCBL(BTC(ISA),X2,X5,X3,X4,SM(136,ISA),SM(75,ISA),&
                      SM(74,ISA),SM(73,ISA),SM(99,ISA),SM(140,ISA),&
                      SM(101,ISA),FTC,ISA,KW,NBSA,MSO,JRT)
	                  IF(JRT>0)WRITE(KW(1),637)ZLSC(ISA),ZLMC(ISA),&
                      ZBMC(ISA),ZHSC(ISA),ZHPC(ISA)
                      FTP=WSAX*(ZPML(ISA)+ZPMU(ISA)+ZPMS(ISA)+ZPMA(ISA)+&
                      ZPO(ISA)+ZFOP(ISA)+ZPOU(ISA)+AD1+STDOP(ISA)+TOT)
                      ! SA P BALANCE
                      IF(IEXT(ISA)>0)THEN
                          X5=SM(48,ISA)
                          X4=SM(49,ISA)
                          X2=0.
                          X3=0.
                      ELSE    
                          X5=SYP(IDOR(ISA))
                          X4=SQP(IDOR(ISA))
                          X2=SYP(IDNF(ISA))
                          X3=SQP(IDNF(ISA))
                      END IF
                      TYP(ISA)=TYP(ISA)*WSAX
                      CALL NBL(BTP(ISA),X2,X5,X3,X4,0.,0.,SM(51,ISA),0.,&
                      SM(135,ISA),0.,SM(56,ISA),TYP(ISA),0.,SM(57,ISA),&
                      0.,0.,FTP,0.,0.,SM(90,ISA),SM(106,ISA),ISA,2,&
                      KW,NBSA,MSO,JRT)
                      IF(JRT>0)WRITE(KW(1),4089)ZPML(ISA),ZPMA(ISA),&
                      ZPMS(ISA),ZPO(ISA),ZFOP(ISA),AD1,STDOP(ISA),TOT,ZPMU(ISA),&
                      ZPOU(ISA)
                      FTK=ZSK(ISA)+ZEK(ISA)+ZFK(ISA)+AD2+STDOK(ISA)+AD3
                      !CALL NBL(BTK(ISA),0.,SM(149,ISA),SM(150,ISA),SM(151,ISA),0.,&
                      !SM(152,ISA),0.,0.,TYK(ISA),0.,SM(148,ISA),0.,0.,FTK,0.,0.,0.,&
                      !0.,0.,0.,0.,ISA,3,KFL,KW,NBSA,MSO,JRT)                                                 
                      IF(JRT>0)WRITE(KW(1),639)ZSK(ISA),ZEK(ISA),&
                      ZFK(ISA),AD2,STDOK(ISA),AD3 
                      CALL SLTB(SM(129,ISA),SM(132,ISA),SM(133,ISA),SM(131,ISA),SM&
                      (130,ISA),SLT0(ISA),ZSLT(ISA),ISA,NBSA,KW,MSO)
                  END IF    
                  XQRB=NQRB(IDX)
                  XCN=JCN(ISA)
                  XX=XQRB+.01
                  PRAX=PRAV(IDX)/XX
                  PRSD(ISA)=PRSD(ISA)-PRAV(IDX)*PRAX
                  IF(PRSD(ISA)>0.)PRSD(ISA)=SQRT(PRSD(ISA)/XX)
                  CYSD(ISA)=CYSD(ISA)-CYAV(ISA)*CYAV(ISA)/XX
                  CYAV(ISA)=CYAV(ISA)/XX
                  IF(CYSD(ISA)>0.)CYSD(ISA)=SQRT(CYSD(ISA)/XX)
                  QRBQ(ISA)=QRBQ(ISA)/XX
                  X5=ERAV(IDX)/XX
                  TCAX=TCAV(IDX)/XX
                  ! DETERMINE AVE ANNUAL VALUES
                  SMY(13,ISA)=0.
                  DO I=1,7
                      VAR(I,ISA)=0.
                  END DO
                  YTX(ISA)=0.
                  YWSA=XYR*WSAX
                  YWSA1=10.*YWSA
                  DO I=1,12
                      RR(I)=RR(I)+TR(I,ISA)/XYR
                      TR(I,ISA)=TR(I,ISA)/YWSA1
                      TSN(I,ISA)=TSN(I,ISA)/XYR
                      TSY(I,ISA)=TSY(I,ISA)/YWSA
                      TYW(I,ISA)=TYW(I,ISA)/XYR
                      TQ(I,ISA)=TQ(I,ISA)/(XYR*WSAX1)
	                  TCN(I,ISA)=TCN(I,ISA)/XYR
                      TYON(I,ISA)=TYON(I,ISA)/YWSA
                      TYTP(I,ISA)=TYTP(I,ISA)/YWSA
                      TQP(I,ISA)=TQP(I,ISA)/YWSA
	                  TQPU(I,ISA)=TQPU(I,ISA)/XYR
                      TQN(I,ISA)=TQN(I,ISA)/YWSA
                      VAR(4,ISA)=VAR(4,ISA)+TYTP(I,ISA)
                      VAR(6,ISA)=VAR(6,ISA)+TQP(I,ISA)
	                  VAR(8,ISA)=VAR(8,ISA)+TQP(I,ISA)
                      TXMX(I,ISA)=TXMX(I,ISA)/XYR
                      TXMN(I,ISA)=TXMN(I,ISA)/XYR
                      TSR(I,ISA)=TSR(I,ISA)/XYR
                      TCVF(I,ISA)=TCVF(I,ISA)/(TEI(I,ISA)+1.E-20)
                      CX(I,ISA)=CX(I,ISA)/XYR
                      VAR(3,ISA)=VAR(3,ISA)+CX(I,ISA)
                      TEI(I,ISA)=TEI(I,ISA)/XYR
                      SMY(I,ISA)=SRD(I,ISA)/XYR
                      YTX(ISA)=YTX(ISA)+SRD(I,ISA)
                      SET(I,ISA)=SET(I,ISA)/XYR
                      TET(I,ISA)=TET(I,ISA)/XYR
                      ASW(I,ISA)=ASW(I,ISA)/XYR
                      QIN(I,ISA)=QIN(I,ISA)/XYR
                      TRHT(I,ISA)=TRHT(I,ISA)/XYR
                      TAMX(I,ISA)=TAMX(I,ISA)/XYR
                      VAR(7,ISA)=VAR(7,ISA)+QIN(I,ISA)
                      VAR(2,ISA)=VAR(2,ISA)+TAMX(I,ISA)
                      VAR(1,ISA)=VAR(1,ISA)+ASW(I,ISA)
                      SMY(13,ISA)=SMY(13,ISA)+SMY(I,ISA)
                  END DO
                  SM(14,ISA)=SM(14,ISA)/(XCN+1.E-10)
                  SMUA(14,ISA)=SM(14,ISA)
                  SM(25,ISA)=SM(25,ISA)/(SM(24,ISA)+1.E-20)
                  SMUA(25,ISA)=SM(25,ISA)
                  VAR(1,ISA)=VAR(1,ISA)/12.
                  DO K=1,13 ! SM(/XYR
                      SM(K,ISA)=SM(K,ISA)/XYR
                      SMUA(K,ISA)=SMUA(K,ISA)/XYR
                  END DO
                  DO K=15,24
                      SM(K,ISA)=SM(K,ISA)/XYR
                      SMUA(K,ISA)=SMUA(K,ISA)/XYR
                  END DO
                  DO K=26,NSM
                      SM(K,ISA)=SM(K,ISA)/XYR
                      SMUA(K,ISA)=SMUA(K,ISA)/XYR
                  END DO
                  IF(NDP>0)THEN
                      IF(KFL(1)>0)THEN
	                      CALL APAGE(1)
                          WRITE(KW(1),'(//1X,A/)')'______________PESTICIDE SUMMARY &
                          TABLE________________'
                          WRITE(KW(1),460)
                      END IF
                      CALL PSTSUM(1.,XYR,NX,ISA)
                      DO K=1,NDP
                          SMPQ(ISA)=SMPQ(ISA)+SMAP(2,K,IDX)
                          SMPY(ISA)=SMPY(ISA)+SMAP(5,K,IDX)
                          SMPL(ISA)=SMPL(ISA)+GWPS(K,ISA)
                      END DO
                  END IF          
                  X1=CST1(ISA)/XYR
                  X2=VALF1(ISA)/XYR
	              IF(KFL(1)>0)THEN
                      CALL APAGE(1)
	                  WRITE(KW(1),'(//1X,A/)')'____________________SUMMARY TABLE&
	                  ____________________'
                      WRITE(KW(1),3320)PRB(IDX),PRAX,PRSD(ISA),QRQB(ISA),QRBQ(ISA),&
                      NQRB(IDX)
                      WRITE(KW(1),417)TCAX,TCMN(IDX),TCMX(IDX)
                      WRITE(KW(1),448)CYAV(ISA),CYMX(ISA),CYSD(ISA)
                      AD1=0.
                      DO K=1,10
                          AD1=AD1+CNDS(K,ISA)
                      END DO
                      DO K=1,10
                          CNDS(K,ISA)=CNDS(K,ISA)/AD1
                      END DO
                      WRITE(KW(1),2460)(CNDS(K,ISA),K=1,10)
                  END IF
                  IF(LUN(ISA)==35)THEN
                      DWOC(ISA)=0.
                  ELSE    
                      DO J=1,6
                          XTP(J)=0.
                          DO I=1,NBSL(ISA)
                              ISL=LID(I,ISA)
                              SMS(J,ISL,ISA)=SMS(J,ISL,ISA)/(SMS(11,ISL,ISA)&
                              +1.E-5)
                              XTP(J)=XTP(J)+SMS(J,ISL,ISA)
                          END DO
                      END DO
                      DO J=7,10
                          XTP(J)=0.
                          DO I=1,NBSL(ISA)
                              ISL=LID(I,ISA)
                              SMS(J,ISL,ISA)=SMS(J,ISL,ISA)/XYR
                              XTP(J)=XTP(J)+SMS(J,ISL,ISA)
                          END DO
                      END DO
                      DO I=1,NBSL(ISA)
                          ISL=LID(I,ISA)
                          XYP(ISL)=WOC(ISL,ISA)-XZP(6,ISL,ISA)
                          YTP(ISL)=WON(ISL,ISA)-XZP(12,ISL,ISA)
                      END DO
                      DWOC(ISA)=ZOC(ISA)-XZP(6,11,ISA)
                      YTP(11)=ZON(ISA)-XZP(12,11,ISA)
                      SMOC=SMOC+.001*WSA(ISA)*DWOC(ISA)
                      XNS(ISA)=NBSL(ISA)
                      DO J=1,6
                          XTP(J)=XTP(J)/XNS(ISA)
                      END DO
                      IF(KFL(23)>0)THEN
                          WRITE(KW(23),635)ISA,NBSA(ISA),(SID(LORG(LID(J,ISA),&
                          ISA)),J=1,10),SID(11)
                          WRITE(KW(23),649)'   Z',(Z(LID(I,ISA),ISA),I=1,10),Z(LID(10,ISA),ISA),' Z   '
                          WRITE(KW(23),649)' SWF',(SMS(1,LID(I,ISA),ISA),I=1,10),XTP(1),' SWF '
                          WRITE(KW(23),649)'TEMP',(SMS(2,LID(I,ISA),ISA),I=1,10),XTP(2),' TEMP'
                          WRITE(KW(23),649)'SWTF',(SMS(3,LID(I,ISA),ISA),I=1,10),XTP(3),' SWTF'
                          WRITE(KW(23),649)'TLEF',(SMS(4,LID(I,ISA),ISA),I=1,10),XTP(4),' TLEF'
                          WRITE(KW(23),649)'SPDM',(SMS(5,LID(I,ISA),ISA),I=1,10),XTP(5),' SPDM'
                          WRITE(KW(23),576)'RSDC',(SMS(7,LID(I,ISA),ISA),I=1,10),XTP(7),' RSDC'
                          WRITE(KW(23),576)'RSPC',(SMS(8,LID(I,ISA),ISA),I=1,10),XTP(8),' RSPC'
                          WRITE(KW(23),576)'RNMN',(SMS(9,LID(I,ISA),ISA),I=1,10),XTP(9),' RNMN'
                          WRITE(KW(23),576)'DNO3',(SMS(10,LID(I,ISA),ISA),I=1,10),XTP(10),' DNO3'
                          WRITE(KW(23),576)'HSC0',(XZP(1,LID(I,ISA),ISA),I=1,10),XZP(1,11,ISA),' HSC0'
                          WRITE(KW(23),576)'HSCF',(WHSC(LID(I,ISA),ISA),I=1,10),ZHSC(ISA),' HSCF'
                          WRITE(KW(23),576)'HPC0',(XZP(2,LID(I,ISA),ISA),I=1,10),XZP(2,11,ISA),' HPC0'
                          WRITE(KW(23),576)'HPCF',(WHPC(LID(I,ISA),ISA),I=1,10),ZHPC(ISA),' HPCF'
                          WRITE(KW(23),576)'LSC0',(XZP(3,LID(I,ISA),ISA),I=1,10),XZP(3,11,ISA),' LSC0'
                          WRITE(KW(23),576)'LSCF',(WLSC(LID(I,ISA),ISA),I=1,10),ZLSC(ISA),' LSCF'
                          WRITE(KW(23),576)'LMC0',(XZP(4,LID(I,ISA),ISA),I=1,10),XZP(4,11,ISA),' LMC0'
                          WRITE(KW(23),576)'LMCF',(WLMC(LID(I,ISA),ISA),I=1,10),ZLMC(ISA),' LMCF'
                          WRITE(KW(23),576)'BMC0',(XZP(5,LID(I,ISA),ISA),I=1,10),XZP(5,11,ISA),' BMC0'
                          WRITE(KW(23),576)'BMCF',(WBMC(LID(I,ISA),ISA),I=1,10),ZBMC(ISA),' BMCF'
                          WRITE(KW(23),576)'WOC0',(XZP(6,LID(I,ISA),ISA),I=1,10),XZP(6,11,ISA),' WOC0'
                          WRITE(KW(23),576)'WOCF',(WOC(LID(I,ISA),ISA),I=1,10),ZOC(ISA),' WOCF'
                          WRITE(KW(23),576)'DWOC',(XYP(LID(I,ISA)),I=1,10),DWOC(ISA),' DWOC'
                          WRITE(KW(23),576)'HSN0',(XZP(7,LID(I,ISA),ISA),I=1,10),XZP(7,11,ISA),' HSN0'
                          WRITE(KW(23),576)'HSNF',(WHSN(LID(I,ISA),ISA),I=1,10),ZHSN(ISA),' HSNF'
                          WRITE(KW(23),576)'HPN0',(XZP(8,LID(I,ISA),ISA),I=1,10),XZP(8,11,ISA),' HPN0'
                          WRITE(KW(23),576)'HPNF',(WHPN(LID(I,ISA),ISA),I=1,10),ZHPN(ISA),' HPNF'
                          WRITE(KW(23),576)'LSN0',(XZP(9,LID(I,ISA),ISA),I=1,10),XZP(9,11,ISA),' LSN0'
                          WRITE(KW(23),576)'LSNF',(WLSN(LID(I,ISA),ISA),I=1,10),ZLSN(ISA),' LSNF'
                          WRITE(KW(23),576)'LMN0',(XZP(10,LID(I,ISA),ISA),I=1,10),XZP(10,11,ISA),' LMN0'
                          WRITE(KW(23),576)'LMNF',(WLMN(LID(I,ISA),ISA),I=1,10),ZLMN(ISA),' LMNF'
                          WRITE(KW(23),576)'BMN0',(XZP(11,LID(I,ISA),ISA),I=1,10),XZP(11,11,ISA),' BMN0'
                          WRITE(KW(23),576)'BMNF',(WBMN(LID(I,ISA),ISA),I=1,10),ZBMN(ISA),' BMNF'
                          WRITE(KW(23),576)'WON0',(XZP(12,LID(I,ISA),ISA),I=1,10),XZP(12,11,ISA),' WON0'
                          WRITE(KW(23),576)'WONF',(WON(LID(I,ISA),ISA),I=1,10),ZON(ISA),' WONF'
                          WRITE(KW(23),576)'DWON',(YTP(LID(I,ISA)),I=1,10),YTP(11),' DWON'
                          WRITE(KW(23),649)'C/N0',(XZP(13,LID(I,ISA),ISA),I=1,10),XZP(13,11,ISA),' C/N0'
                          DO I=1,100
                              XTP(I)=0.
                          END DO
                          DO I=1,NBSL(ISA)
                              ISL=LID(I,ISA)
                              XTP(ISL)=WOC(ISL,ISA)/WON(ISL,ISA)
                          END DO
                          XTP(11)=ZOC(ISA)/ZON(ISA)
                          WRITE(KW(23),649)'C/NF',(XTP(LID(I,ISA)),I=1,10),XTP(11),' C/NF'
                      END IF
                  END IF
                  PPX(1,ISA)=SMUA(38,ISA)
                  PPX(2,ISA)=SMUA(39,ISA)
                  PPX(3,ISA)=SMUA(40,ISA)
                  PPX(4,ISA)=SMUA(49,ISA)
                  IF(KFL(1)>0)THEN
                      WRITE(KW(1),3050)
                      ! PRINTOUT SUMMARY MONTHLY
                      WRITE(KW(1),3280)IY
                      WRITE(KW(1),3300)HED(1),(TXMX(I,ISA),I=1,12),SM(1,ISA),HED(1)
                      WRITE(KW(1),3300)HED(2),(TXMN(I,ISA),I=1,12),SM(2,ISA),HED(2)
                      WRITE(KW(1),2550)HED(4),(TR(I,ISA),I=1,12),SMUA(4,ISA),HED(4)
                      WRITE(KW(1),3300)'DAYP',(SMY(I,ISA),I=1,13),'DAYP'
                      WRITE(KW(1),2550)HED(16),(TSN(I,ISA),I=1,12),SMUA(16,ISA),HED(16)
                      WRITE(KW(1),2550)HED(13),(TQ(I,ISA),I=1,12),SMUA(13,ISA),HED(13)
	                  WRITE(KW(1),2550)HED(14),(TCN(I,ISA),I=1,12),SM(14,ISA),HED(14)
                      WRITE(KW(1),2550)'DAYQ',(CX(I,ISA),I=1,12),VAR(3,ISA),'DAYQ'
	                  WRITE(KW(1),3300)'RZSW',(ASW(I,ISA),I=1,12),VAR(1,ISA),'RZSW'
                      WRITE(KW(1),3200)HED(24),(TEI(I,ISA),I=1,12),SM(24,ISA),HED(24)
                      WRITE(KW(1),3301)HED(25),(TCVF(I,ISA),I=1,12),SM(25,ISA),HED(25)
                      WRITE(KW(1),3300)HED(NDVSS),(TSY(I,ISA),I=1,12),SMUA(NDVSS,ISA),HED(NDVSS)
                      WRITE(KW(1),3300)HED(37),(TYON(I,ISA),I=1,12),SMUA(37,ISA),HED(37)
                      WRITE(KW(1),3300)HED(48),(TYTP(I,ISA),I=1,12),SMUA(48,ISA),HED(48)
                      WRITE(KW(1),3300)HED(38),(TQN(I,ISA),I=1,12),SMUA(38,ISA),HED(38)
                      WRITE(KW(1),3300)HED(49),(TQP(I,ISA),I=1,12),SMUA(49,ISA),HED(49)
                      WRITE(KW(1),3300)HED(11),(TET(I,ISA),I=1,12),SM(11,ISA),HED(11)
                !     WRITE(KW(1),2550)'DAYW',(TAMX(I,ISA),I=1,12),VAR(2,ISA),'DAYW'
                      WRITE(KW(1),2550)HED(33),(TRHT(I,ISA),I=1,12),SM(33,ISA),HED(33)
                      WRITE(KW(1),3300)HED(36),(TYW(I,ISA),I=1,12),SM(36,ISA),HED(36)
                      WRITE(KW(1),3300)HED(19),(QIN(I,ISA),I=1,12),VAR(7,ISA),HED(19)
                      WRITE(KW(1),3200)HED(10),(SET(I,ISA),I=1,12),SMUA(10,ISA),HED(10)
                      WRITE(KW(1),3200)HED(3),(TSR(I,ISA),I=1,12),SM(3,ISA),HED(3)
                      WRITE(KW(1),3)'RAMX',(SRMX(I,ISA),I=1,12),'RAMX'
                      WRITE(KW(1),3010)'HRLT',(THRL(I,ISA),I=1,12),'HRLT'
                      WRITE(KW(1),'(/1X,A)')'-----AVE ANNUAL VALUES'
                      WRITE(KW(1),3030)IY,(HED(KA(K)),SMUA(KA(K),ISA),K=5,NKA),'COST',X1,&
                      'RTRN',X2
                      WRITE(KW(1),3770)(HED(JC(K)),PPX(K,ISA),K=1,NJC)
                  END IF
                  X4=SMUA(53,ISA)+SMUA(54,ISA)+SMUA(55,ISA)
                  X1=SMUA(56,ISA)+SMUA(57,ISA)
                  FSFN(ISA)=FSFN(ISA)/(X4*XYR+1.E-5)
                  FSFP(ISA)=FSFP(ISA)/(X1*XYR+1.E-5)
                  SM81=SM81+SM(81,ISA)
                  SM(81,ISA)=.001*SM(81,ISA)
                  LD1=LID(1,ISA)
                  X2=1000.*WPML(LD1,ISA)/WT(LD1,ISA)
                  X3=1000.*SM(49,ISA)/(SM(13,ISA)+.1)
                  XYRD=XYR*365.25
                  STKR(ISA)=STKR(ISA)/XYRD
                  PSTM(ISA)=PSTM(ISA)/XYR
                  K1=0
	              DO K=1,LC
                      !IF(NCR(K,ISA)==0)CYCLE
                !     K=LY(IRO(ISA),J,ISA)
                      !XX=MIN(NCR(K,ISA),IY)+1.E-10
                      XX=IY
                      TETG(K,ISA)=1000.*TYL1(K,ISA)/(TETG(K,ISA)+1.E-10)
                      SMY1(ISA)=SMY1(ISA)+TYL1(K,ISA)
                      SMY2(ISA)=SMY2(ISA)+TYL2(K,ISA)
                      TYL1(K,ISA)=TYL1(K,ISA)/XX
                      TYL2(K,ISA)=TYL2(K,ISA)/XX
                      TYLN(K,ISA)=TYLN(K,ISA)/XX
                      TYLP(K,ISA)=TYLP(K,ISA)/XX
                      TYLK(K,ISA)=TYLK(K,ISA)/XX
                      TDM(K,ISA)=TDM(K,ISA)/XX
	                  THU(K,ISA)=THU(K,ISA)/XX
	                  XX=IY
                      TRA(K,ISA)=TRA(K,ISA)/XX
                      TRD(K,ISA)=TRD(K,ISA)/XX
                      TCAW(K,ISA)=TCAW(K,ISA)/XX
                      TVIR(K,ISA)=TVIR(K,ISA)/XX
                      TFTN(K,ISA)=TFTN(K,ISA)/XX
                      TFTP(K,ISA)=TFTP(K,ISA)/XX
                      TFTK(K,ISA)=TFTK(K,ISA)/XX
                      SMWS(ISA)=SMWS(ISA)+TSFC(1,K,ISA)
                      SMNS(ISA)=SMNS(ISA)+TSFC(2,K,ISA)
                      SMPS(ISA)=SMPS(ISA)+TSFC(3,K,ISA)
                      SMKS(ISA)=SMKS(ISA)+TSFC(4,K,ISA)
                      SMTS(ISA)=SMTS(ISA)+TSFC(5,K,ISA)
                      SMAS(ISA)=SMAS(ISA)+TSFC(6,K,ISA)
                      SMSS(ISA)=SMSS(ISA)+TSFC(7,K,ISA)
                      DO L=1,7
                          TSFC(L,K,ISA)=TSFC(L,K,ISA)/XX
                      END DO
                      DO L=1,3
                          STDA(L,K,ISA)=STDA(L,K,ISA)/XX
                      END DO
                      IF(PSTM(ISA)<1.E-10)THEN
                          X6=1.
                      ELSE
                          X6=PSTM(ISA)
                      END IF        
                      ! PRINTOUT CROP SUMMARY
                      IF(KFL(1)>0)WRITE(KW(1),326)CPNM(K),TYL1(K,ISA),TYL2(K,ISA),&
                      TDM(K,ISA),TYLN(K,ISA),TYLP(K,ISA),TYLK(K,ISA),TFTN(K,ISA),&
                      TFTP(K,ISA),TFTK(K,ISA),TVIR(K,ISA),TCAW(K,ISA),TETG(K,ISA),&
                      TRA(K,ISA),THU(K,ISA),X6
                      IF(KFL(1)>0)WRITE(KW(1),577)(TSFC(L,K,ISA),L=1,7),(STDA(L,K,ISA),L=1,3)
                      IF(K1>0)THEN
                          IF(KFL(2)>0.AND.SM(81,ISA)>0.)WRITE(KW(2),2282)ISA,NBSA(ISA),&
                          IDON(ISA),CPNM(K),TYL1(K,ISA),TYL2(K,ISA),TYLN(K,ISA),TYLP(K,ISA)
                      ELSE
                          IF(KFL(2)>0.AND.SM(81,ISA)>0.)THEN
                              SM80=SM(80,ISA)
                              WRITE(KW(2),2282)ISA,NBSA(ISA),IDON(ISA),&
                              CPNM(K),TYL1(K,ISA),TYL2(K,ISA),TYLN(K,ISA),&
                              TYLP(K,ISA),STKR(ISA),WSAX,SMUA(13,ISA),&
                              SM(29,ISA),SMUA(49,ISA),SMUA(48,ISA),SMUA(38,ISA),&
                              SMUA(39,ISA),SMUA(80,ISA),SMUA(40,ISA),SMUA(37,ISA),X1,X4,&
                              SMUA(81,ISA),PDPL0(ISA),PDPLC(ISA),X3
                          END IF    
                          K1=1
                      END IF
                      SM2=SM2+TYLP(K,ISA)*WSAX
                      SM3=SM3+TYLN(K,ISA)*WSAX
                      IF(KFL(49)/=0)THEN
                          IF(K==1)THEN
                              ADFL=ADJUSTL(ASTN)
                              WRITE(KW(49),105)IYER,IMON,IDAY,ISA,&
                              NBSA(ISA),NQRB(IDX),ADFL,WSA(ISA),PRB(IDX),PRAX,&
                              PRSD(ISA),SMUA(4,ISA),SMUA(13,ISA),SMUA(16,ISA),&
                              SMUA(24,ISA),SMUA(25,ISA),SMUA(26,ISA),&
                              SMUA(27,ISA),SMUA(28,ISA),SMUA(29,ISA),&
                              SMUA(30,ISA)
                          END IF
                      END IF    
                  END DO
                  X2=NWDA(ISA)
                  YW0=100.*WCF*TVGF(ISA)/(X2+1.E-5)
                  ! PRINTOUT SUBAREA SUMMARY
                  IF(KFL(3)>0)WRITE(KW(3),498)ISA,NBSA(ISA),WSA(ISA),CN2(ISA),YW0,&
                  OCPD(ISA),FSFN(ISA),FSFP(ISA),PRB(IDX),PRAX,TCAX,CYAV(ISA),&
                  CYMX(ISA),(SMUA(KY(K),ISA),K=1,NKY),SMUA(NDVSS,ISA)
                  SM1=SM1+X1*WSAX
                  SM4=SM4+X4*WSAX
                  SM40=SM40+SM(40,ISA)
                  SM71=SM71+SM(71,ISA)
	              SM104=SM104+SM(104,ISA)
                  SM146=SM146+SM(146,ISA)
                  SM147=SM147+SM(147,ISA)
                  SM155=SM155+SM(155,ISA)
                  SM156=SM156+SM(156,ISA)
                  SM18=SM18+SM(18,ISA)
	              SM19=SM19+SM(19,ISA)
                  SM110=SM110+SM(110,ISA)
              END DO    !SUBAREA LOOP OUTPUT
              IF(KFL(1)>0)WRITE(KW(1),461)
              IF(NDP>0)THEN
                  IF(KFL(1)>0)WRITE(KW(1),460)
                  IF(MSA>1)CALL PSTSUM(SUMA,XYR,NX,NCMD)
              END IF
              XYA=XYR*SUMA
              AVRF=SM9/XYA
              DO I=1,12
                  RR(I)=.1*RR(I)/SUMA
                  RR(13)=RR(13)+RR(I)
              END DO
              X1=YTX(1)-1.
              X2=XYR*SM(100,1)
              SX2=SDRF-X2*X2/YTX(1)
              SDX=SQRT(SX2/(X1+1.E-10))
              X3=X2/YTX(1)
              SMH34=SMH(34,NCMD)
              IF(KFL(35)>0)THEN
	              PCTH=0.
	              PCT=0.
	              NTX=0
                  DO IDO=1,NCMD
                      SMH(7,IDO)=1.E5*SMH(6,IDO)/(SMH(2,IDO)*RWSA(IDO)+1.E-10)
		              DO K=1,NSH
		                  SMH(K,IDO)=SMH(K,IDO)/XYR
		                  IF(K<5.OR.K>32)THEN 
		                      X1=366-NYD
		                      SMH(K,IDO)=SMH(K,IDO)/X1
		                  END IF
		              END DO
                  END DO
	          END IF
	          AD5=0.
	          AD1=0.
	          AD2=0.
              DO ISA=1,MSA
                  WSAX=WSA(ISA)
                  WSAX1=10.*WSAX
                  YWSA=XYR*WSAX
                  YWSA1=10.*YWSA
                  IF(SM(141,ISA)>0.)THEN
                      AD1=AD1+SM(141,ISA)*WSAX
                      AD2=AD2+WSAX
                  END IF
                  AD5=AD5+TYN(ISA)
                  X2=XYR*SM(4,ISA)/WSAX1
                  SX2=SDVR(ISA)-X2*X2/(YTX(ISA)+1.E-10)
                  SDX=SQRT(SX2/X1)
                  X3=X2/YTX(ISA)
                  IF(KFL(34)>0)THEN
                      DO I=1,MSA
                          I1=NBSA(IBSA(I))
                          I2=NISA(I1)	
                          II=IDOA(I2)
                          X1=.001*ZOC(I2)
                          X71=SM(71,I2)/WSA(I2)
                          WRITE(KW(34),859)I1,NBYR,IGC,WSA(I2),SMUA(4,I2),SMUA(5,I2),&
                          SMUA(6,I2),SMUA(18,I2),SMUA(10,I2),SMUA(11,I2),RZSW(I2),&
                          SMUA(16,I2),X71,SMUA(13,I2),SMUA(15,I2),SMUA(72,I2),&
                          SMUA(117,I2),SMUA(14,I2),SMUA(1,I2),SMUA(2,I2),SMUA(59,I2),&
                          SMUA(3,I2),SMUA(27,I2),SMUA(31,I2),SMUA(53,I2),SMUA(54,I2),&
                          SMUA(55,I2),SMUA(56,I2),SMUA(57,I2),SMUA(43,I2),SMUA(42,I2),&
                          SMUA(37,I2),SMUA(119,I2),SMUA(38,I2),SMUA(49,I2),SMUA(118,I2),&
                          SMUA(39,I2),SMUA(80,I2),X1,(TYL1(LY(IRO(I2),J,I2),I2),&
                          TYL2(LY(IRO(I2),J,I2),I2),(TSFC(L,J,I2),L=1,7),CPNM&
                          (LY(IRO(I2),J,I2)),J=1,NCP(IRO(I2),I2))
                      END DO
                  END IF
                  IF(KFL(40)>0)THEN
                      DO J=1,NCP(IRO(ISA),ISA)
                          YTP(J)=1000.*TYL2(LY(IRO(ISA),J,ISA),ISA)
                      END DO	
                      WRITE(KW(40),857)NBSA(ISA),SMUA(4,ISA),SMUA(13,ISA),SMUA(117,ISA),&
                      SMUA(16,ISA),(CPNM(LY(IRO(ISA),J,ISA)),(TSFC(L,J,ISA),L=1,7),&
                      YTP(J),SMUA(141,ISA),J=1,NCP(IRO(ISA),ISA))
                  END IF
  	              IF(KFL(35)==0.OR.MSA==1)CYCLE
                  I1=NBSA(IBSA(ISA))
	              I2=NISA(I1)
	              IF(IEXT(I2)>0)THEN
	                  II=IDOA(I2)
	                  X2=WSA(I2)
	                  I3=II
	              ELSE
	                  II=IDOA(I2)-1
	                  X2=RWSA(II)+WSA(I2)
	                  I3=II+2
	              END IF
                  IF(NTX(II)>0)CYCLE
	              WRITE(KW(35),5005)I1,II,NBYR,X2,SMH(1,II),SMH(2,I3),SMH(33,II),&
	              SMH(35,I3),(SMH(K,II),K=3,NSH-3)
                  !(PCTH(J,II),J=1,NSZ),(PCT(J,II),J=1,NSZ)
 	              NTX(II)=1
              END DO
	          AD1=AD1/(AD2+1.E-10)
              IF(KFL(39)>0)THEN
                  DO I=1,NCMO
	                  II=ICMO(I)
	                  WRITE(KW(39),472)II,IYR,MO,RWSA(II),SMH(1,II),SMH(2,II),&
	                  SMH(33,II),SMH(34,II),(SMH(K,II),K=3,NSH-2)
	              END DO    
	          END IF    
              DO IWI=1,NWP
                  X1=NWPD(IWI)
                  RTO=X1/XYRD
              END DO
              IF(NDP>0.AND.(KFL(1)>0.OR.KFL(38)>0))THEN
                  CALL APAGE(0)
                  WRITE(KW(1),18)
                  WRITE(KW(1),3281)NBYR
                  YTP=0.
                  DO K=1,NDP
                      WRITE(KW(1),3060)PSTN(K)
	                  SUM=0.
                      DO I=1,12
                          SMRP(5,K,I)=SMRP(5,K,I)/XYR
                          SUM=SUM+SMRP(5,K,I)
                      END DO
                      WRITE(KW(1),3302)HEDP(1),(SMRP(5,K,I),I=1,12),SUM,HEDP(1)
                      DO J=1,2
                          J1=J+2
                          XTP(J)=0.
                          XTP(J1)=0.
                          DO I=1,12
                              SMRP(J,K,I)=SMRP(J,K,I)/XYR
                              SMRP(J1,K,I)=SMRP(J1,K,I)/XYR
                              XTP(J1)=XTP(J1)+SMRP(J1,K,I)
                              XTP(J)=XTP(J)+SMRP(J,K,I)
                          END DO
                          YTP(J1)=YTP(J1)+XTP(J1)
                          YTP(J)=YTP(J)+XTP(J)
                          ! PRINTOUT REACH PESTICIDE SUMMARY MONTHLY
                          WRITE(KW(1),3302)HDRP(J),(SMRP(J,K,I),I=1,12),XTP(J),HDRP(J)
                          WRITE(KW(1),3302)HDRP(J1),(SMRP(J1,K,I),I=1,12),XTP(J1),HDRP(J1)
                      END DO
                  END DO        
                  CALL APAGE(0)
              END IF          
              SUM=0.
              DO I=1,LC
                  TCPY(I)=TCPY(I)/(TCPA(I)+1.E-10)
                  TCPA(I)=TCPA(I)/XYR
                  SUM=SUM+TCPA(I)
              END DO
              IF(KFL(1)>0)WRITE(KW(1),5007)
              IF(KFL(10)>0)WRITE(KW(10),5007)
              DO I=1,LC
                  X1=TCPA(I)/(SUM+1.E-10)
                  IF(KFL(1)>0)WRITE(KW(1),5008)CPNM(I),TCPA(I),X1,TCPY(I)
                  IF(KFL(10)>0)WRITE(KW(10),5008)CPNM(I),TCPA(I),X1,TCPY(I)
              END DO
              IF(KFL(1)>0)THEN
                  WRITE(KW(1),'(/T10,A,F10.0/T10,A,F10.2,A)')'AREA WEIGHTED GRAZING DAYS = ',&
                  AD1,'GRAZING AREA = ',AD2,' ha'
                  IF(NGN==0)THEN
                      NAD=0
                      DO IWI=1,NWP
                          NAD=NAD+IAD(IWI)
                      END DO
                      IF(NAD>0)THEN
                          AD1=0.
                          DO IWI=1,NWP
                              XTP(IWI)=REAL(IAD(IWI))/REAL(NAD)
                              RFSG(IWI)=RFSG(IWI)/XYR
                              AD1=AD1+RFSG(IWI)
                          END DO
                          WRITE(KW(1),'(/T10,A,10(1X,I4,F10.4))')'FRACTION WPM1 USED IN SPATIAL &
                          GENERATOR = ',(IWI,XTP(IWI),IWI=1,NWP)
                          WRITE(KW(1),'(/T10,A,10(1X,I4,F10.0))')'PCP GENERATED BY &
                          WPM1 = ',(IWI,RFSG(IWI),IWI=1,NWP)
                          WRITE(KW(1),'(T10,A,F10.0)')' TOT = ',AD1
                      END IF
                  END IF    
              END IF
              IF(KFL(10)>0)THEN
                  WRITE(KW(10),'(/T42,A)')'AVE ANNUAL DATA'
                  WRITE(KW(10),18)
                  WRITE(KW(10),3281)NBYR
                  WRITE(KW(10),2550)HED(4),(RR(I),I=1,12),RR(13),HED(4)
              END IF
              IF(KFL(1)>0)THEN
                  IF(NAQ>0)THEN
                      SUM=0.
                      DO I=1,16
                          XTP(I)=NWDR(I)
                          SUM=SUM+XTP(I)
                      END DO
                      DO I=1,16
                          XTP(I)=XTP(I)/SUM
                      END DO
                      WRITE(KW(1),5012)(XTP(I),I=1,16)
                      WRITE(KW(1),5016)
                      TOT=0.
                      ADD=0.
                      DO I=1,MSA
                          TOT=TOT+PM10(I)
                          ADD=ADD+EM10(I)
                      END DO
                      CALL ASORT1(PM10,NX,MSA)
                      SUM=0.
                      DO I=1,MSA
                          I1=NX(I)
                          X1=PM10(I1)/TOT
                          SUM=SUM+X1
                          X2=PM10(I1)/WSA(I1)
                          WRITE(KW(1),5017)I1,NBSA(I1),X2,X1,SUM
                      END DO
                      WRITE(KW(1),5018)TOT,ADD
                  END IF
                  WRITE(KW(1),18)
                  WRITE(KW(1),3280)NBYR
                  WRITE(KW(1),2550)HED(4),(RR(I),I=1,12),RR(13),HED(4)
              END IF
              DO J=1,10
                  J1=J+10
                  XTP(J)=0.
                  XTP(J1)=0.
                  DO I=1,12
                      IF(J<3)THEN
                          X1=.1
                      ELSE
                          X1=1.
                      END IF    
                      SMR(J,I)=X1*SMR(J,I)/XYA                                                 
                      SMR(J1,I)=X1*SMR(J1,I)/XYA
                      XTP(J1)=XTP(J1)+SMR(J1,I)
                      XTP(J)=XTP(J)+SMR(J,I)
                  END DO
                  ! PRINTOUT REACH SUMMARY MONTHLY
                  IF(KFL(1)>0)THEN
                      WRITE(KW(1),3300)HEDR(J),(SMR(J,I),I=1,12),XTP(J),HEDR(J)
                      WRITE(KW(1),3300)HEDR(J1),(SMR(J1,I),I=1,12),XTP(J1),HEDR(J1)
                  END IF    
                  IF(KFL(10)==0)CYCLE
                  WRITE(KW(10),3301)HEDR(J),(SMR(J,I),I=1,12),XTP(J),HEDR(J)
                  WRITE(KW(10),3301)HEDR(J1),(SMR(J1,I),I=1,12),XTP(J1),HEDR(J1)
              END DO
              XTP(21)=0.
              DO I=1,12
	              SMR(21,I)=SMR(21,I)/XYA
                  XTP(21)=XTP(21)+SMR(21,I)
              END DO
              IF(KFL(32)>0)THEN
                  X1=XTP(21)-XTP(16)
	              SMOC=SMOC/SUMA
                  WRITE(KW(32),2157)XTP(11),XTP(13),XTP(16),X1,XTP(14),SMOC
              END IF
              IF(KFL(1)>0)THEN
                  CALL APAGE(0)
                  WRITE(KW(1),'(//1X,A/)')'______________ROUTING REACH SUMMARY &
                  TABLE_____________'
                  WRITE(KW(1),'(T24,A)')'INDIVIDUAL STORM DATA'
                  WRITE(KW(1),35)
                  K=0
                  DO J=1,NCMD
                      IF(NQRB(J)==0)CYCLE
                      XX=NQRB(J)
                      PRAV(J)=PRAV(J)/XX
                      TCAV(J)=TCAV(J)/XX
                      DRAV(J)=DRAV(J)/XX
                      ERAV(J)=ERAV(J)/XX
                      J1=NQRB(J)/NBYR
                      SR7=.1*SRCH(7,J)/RWSA(J)
                      X1=SR7/XX
                      IF(ICDT(J)/=1)THEN
                          ! PRINTOUT REACH STORM SUMMARY
                          WRITE(KW(1),24)CMDX(ICDT(J)),J,RWSA(J),J1,TCAV(J),DRAV(J),&
                          ERAV(J),X1,PRAV(J),PRB(J)
                      ELSE
                          K=K+1
                          X2=MAX(NBCT(K),NBCF(K))
                          X3=MAX(NBFT(K),NBFF(K))
                          VCHA(K)=VCHA(K)/(X2+1.E-10)
                          VFPA(K)=VFPA(K)/(X3+1.E-10)
                          WRITE(KW(1),4087)CMDX(ICDT(J)),J,K,NBSA(K),RWSA(J),J1,TCAV(J),&
                          DRAV(J),ERAV(J),X1,PRAV(J),PRB(J),NBCF(K),VCHA(K),VCHB(K)&
                          ,NBFF(K),VFPA(K),VFPB(K)
                      END IF
                  END DO
                  WRITE(KW(1),'(///T42,A)')'AVE ANNUAL DATA'
                  WRITE(KW(1),36)
              END IF    
              IF(KFL(10)>0)WRITE(KW(10),36)
              K=0
              DO J=1,NCMD  ! COMMAND OUTPUT LOOP
                  DO I=7,25
                      SRCH(I,J)=SRCH(I,J)/XYR
                  END DO
                  SR7=.1*SRCH(7,J)/RWSA(J)
                  SR8=SRCH(8,J)/RWSA(J)
                  SR9=SRCH(9,J)/RWSA(J)
                  SR10=SRCH(10,J)/RWSA(J)
                  SR11=SRCH(11,J)/RWSA(J)
                  SR12=SRCH(12,J)/RWSA(J)
                  SR13=SRCH(13,J)/RWSA(J)
                  SR14=SRCH(14,J)/RWSA(J)
                  SR17=SRCH(17,J)/RWSA(J)
                  SR18=SRCH(18,J)/RWSA(J)
                  SR21=SRCH(21,J)/RWSA(J)
                  IF(ICDT(J)==1)THEN
                      K=K+1
                      I1=NISA(NBSA(K))
                      IF(KFL(1)>0)WRITE(KW(1),2302)CMDX(ICDT(J)),J,K,NBSA(K),RWSA(J),&
                      SR7,SRCH(15,J),SRCH(19,J),SRCH(20,J),SRCH(16,J),SRCH(23,J),&
                      SR8,SR13,SR14,SRCH(24,J),SRCH(25,J),&
                      SR9,SR11,SR10,SR12,SR17,SR21,SRCH(22,J),SR18
                      IF(KFL(10)>0)WRITE(KW(10),2302)CMDX(ICDT(J)),J,K,NBSA(K),&
                      RWSA(J),SR7,SRCH(15,J),SRCH(19,J),SRCH(20,J),SRCH(16,&
                      J),SRCH(23,J),SR8,SR13,SR14,SRCH(24,J),&
                      SRCH(25,J),SR9,SR11,SR10,SR12,SR17,SR21,SRCH(22,J),SR18
                      CYCLE
                  END IF            
                  IF(ICDT(J)==4)THEN
                      X1=.1/RWSA(J)
                      X2=X1*RSVP(I1)
                      X3=X1*RSVE(I1)
                      J1=I1
                      J2=NBSA(K)
                  ELSE
                      X2=0.
                      X3=0.
                      J1=0
                      J2=0        
                  END IF
                  ! PRINTOUT REACH SUMMARY
                  IF(KFL(1)>0)WRITE(KW(1),2302)CMDX(ICDT(J)),J,J1,J2,RWSA(J),SR7,&
                  SRCH(15,J),SRCH(19,J),SRCH(20,J),SRCH(16,J),SRCH(23,J),&
                  SR8,SR13,SR14,SRCH(24,J),SRCH(25,J),SR9,&
                  SR11,SR10,SR12,SR17,SR21,SRCH(22,J),SR18,X2,X3
                  IF(KFL(10)>0)WRITE(KW(10),2302)CMDX(ICDT(J)),J,J1,J2,RWSA(J),&
                  SR7,SRCH(15,J),SRCH(19,J),SRCH(20,J),SRCH(16,J),&
                  SRCH(23,J),SR8,SR13,SR14,SRCH(24,J),SRCH(25,J),&
                  SR9,SR11,SR10,SR12,SR17,SR21,SRCH(22,J),SR18,X2,X3
              END DO ! COMMAND OUTPUT LOOP
              SR7=.1*SRCH(7,NCMD)/SUMA
              IF(KFL(38)>0)THEN
                  ADRF=0.
                  ADY2=0.
                  ADR2=0.
                  ADYW=0.
                  ADST=0.
                  ADDC=0.
                  ADWS=0.
                  ADKS=0.
                  ADNS=0.
                  ADPS=0.
                  ADTS=0.
                  ADAS=0.
                  ADSS=0.
                  ADLN=0.
                  ADLP=0.
                  ADPL=0.
                  ADY1=0.
                  ADRC=0.
                  ADFU=0.
                  ADFN=0.
                  ADFP=0.
                  ARFP=0. 
                  ADRP=0.
                  AYWP=0.
                  ADPQ=0.
                  ADPY=0.
                  NTX=0
                  DO ISA=1,MSA
                      WSAX=WSA(ISA)
                      ADRF=ADRF+SM(4,ISA)
                      ADR2=ADR2+SM(31,ISA)
                      ADYW=ADYW+SM(36,ISA)
                      ADLN=ADLN+SM(40,ISA)
                      ADLP=ADLP+SM(51,ISA)
                      ADRC=ADRC+SM(25,ISA)*WSAX
                      ARFP=ARFP+SM(142,ISA)
                      ADRP=ADRP+SM(143,ISA)
                      AYWP=AYWP+SM(135,ISA)
                      SMWS(ISA)=SMWS(ISA)/XYR
                      ADWS=ADWS+SMWS(ISA)*WSAX
                      SMNS(ISA)=SMNS(ISA)/XYR
                      ADNS=ADNS+SMNS(ISA)*WSAX
                      SMPS(ISA)=SMPS(ISA)/XYR
                      ADPS=ADPS+SMPS(ISA)*WSAX
                      SMKS(ISA)=SMKS(ISA)/XYR
                      ADKS=ADKS+SMKS(ISA)*WSAX
                      SMTS(ISA)=SMTS(ISA)/XYR
                      ADTS=ADTS+SMTS(ISA)*WSAX
                      SMAS(ISA)=SMAS(ISA)/XYR
                      ADAS=ADAS+SMAS(ISA)*WSAX
                      SMSS(ISA)=SMSS(ISA)/XYR
                      ADSS=ADSS+SMSS(ISA)*WSAX
                      SMY1(ISA)=SMY1(ISA)/XYR
                      ADY1=ADY1+SMY1(ISA)*WSAX
                      SMY2(ISA)=SMY2(ISA)/XYR
                      ADY2=ADY2+SMY2(ISA)*WSAX
                      SMPL(ISA)=SMPL(ISA)/XYR
                      ADPL=ADPL+SMPL(ISA)*WSAX
                      SMPQ(ISA)=SMPQ(ISA)/XYR
                      ADPQ=ADPQ+SMPQ(ISA)*WSAX
                      SMPY(ISA)=SMPY(ISA)/XYR
                      ADPY=ADPY+SMPY(ISA)*WSAX
                      SMFU(ISA)=SMFU(ISA)/XYR
                      ADFU=ADFU+SMFU(ISA)*WSAX
                      SMST(ISA)=SMST(ISA)/XYR
                      ADST=ADST+SMST(ISA)*WSAX
                      DWOC(ISA)=.001*DWOC(ISA)
                      ADDC=ADDC+DWOC(ISA)*WSAX
                      X4=SM(53,ISA)+SM(54,ISA)+SM(55,ISA)
                      ADFN=ADFN+X4
                      X1=SM(56,ISA)+SM(57,ISA)
                      ADFP=ADFP+X1
                      I1=NBSA(IBSA(ISA))
                      I2=NISA(I1)
                      IF(IEXT(I2)>0)THEN
                          II=IDOA(I2)
                          X2=WSA(I2)
                          I3=II
                      ELSE
                          II=IDOA(I2)-1
                          X2=RWSA(II)+WSA(I2)
                          I3=II+2
                      END IF
                      IF(NTX(II)>0)CYCLE
                      CM3MM=3155760./X2
	                  SMH(2,I3)=SMH(2,I3)*CM3MM
 	                  SMH(34,I3)=SMH(34,I3)*CM3MM
                      SMH(35,I3)=SMH(35,I3)*CM3MM
                      SMH(6,I3)=SMH(6,I3)/X2
                      SMH(13,I3)=SMH(13,I3)/X2
                      SMH(9,I3)=SMH(9,I3)/X2
                      SMH(19,I3)=SMH(19,I3)/X2
                      SMH(11,I3)=SMH(11,I3)/X2
                      SMH(27,I3)=SMH(27,I3)/X2
                      SMH(29,I3)=SMH(29,I3)/X2
                      SM(135,ISA)=SM(135,ISA)/X2
                      X1=X1/X2
                      X4=X4/X2
                      WRITE(KW(38),586)ISA,NBSA(ISA),SMUA(4,ISA),SMUA(13,ISA),SMUA(117,ISA),&
                      SMUA(31,ISA),SMUA(NDVSS,ISA),SMUA(36,ISA),SMY1(ISA),SMY2(ISA),SMWS(ISA),&
                      SMNS(ISA),SMPS(ISA),SMKS(ISA),SMTS(ISA),SMAS(ISA),SMSS(ISA),X4,X1,&
                      SMUA(18,ISA),SMST(ISA),SMFU(ISA),DWOC(ISA),SMUA(38,ISA),SMUA(39,ISA),&
                      SMUA(84,ISA),SMUA(47,ISA),SMUA(80,ISA),SMUA(96,ISA),SMUA(37,ISA),&
                      SMUA(134,ISA),SMUA(46,ISA),SMUA(42,ISA),SMUA(43,ISA),SMUA(49,ISA),&
                      SMUA(142,ISA),SMUA(143,ISA),SMUA(51,ISA),SMUA(48,ISA),SMUA(135,ISA),&
                      SMPQ(ISA),SMPL(ISA),SMPY(ISA),SMH(2,I3),SMH(35,I3),SMH(6,I3),&
                      SMH(13,I3),SMH(9,I3),SMH(19,I3),SMH(11,I3),SMH(27,I3),SMH(29,I3)
                      NTX(II)=1
                  END DO
                  IF(MSA>1)THEN
                      ADRF=.1*ADRF/SUMA
                      ADY2=ADY2/SUMA
                      ADY1=ADY1/SUMA
                      ADR2=ADR2/SUMA
                      ADYW=ADYW/SUMA
                      ADST=ADST/SUMA
                      ADDC=ADDC/SUMA
                      ADWS=ADWS/SUMA
                      ADNS=ADNS/SUMA
                      ADPS=ADPS/SUMA
                      ADKS=ADKS/SUMA
                      ADTS=ADTS/SUMA
                      ADAS=ADAS/SUMA
                      ADSS=ADSS/SUMA
                      ADLN=ADLN/SUMA
                      ADLP=ADLP/SUMA
                      ADPL=ADPL/SUMA
                      ADRC=ADRC/SUMA
                      ADFU=ADFU/SUMA
                      ADFN=ADFN/SUMA
                      ADFP=ADFP/SUMA
                      ARFP=ARFP/SUMA
                      ADRP=ADRP/SUMA
                      AYWP=AYWP/SUMA
                      ADIR=SM18/XYA
                      ADDP=SM96/SUMA
                      ADWN=SM134/XYA
                      ADVO=SM46/XYA
                      ADDN=SM42/XYA
                      ADFX=SM43/XYA
                      ADPQ=ADPQ/XYA
                      SR8=SRCH(8,NCMD)/SUMA
                      SR10=SRCH(10,NCMD)/SUMA
                      SR12=SRCH(12,NCMD)/SUMA
                      SR9=SRCH(9,NCMD)/SUMA
                      SR18=SRCH(18,NCMD)/SUMA
                      SR22=SRCH(22,NCMD)/SUMA
                      SR21=SRCH(21,NCMD)/SUMA
                      SR17=SRCH(17,NCMD)/SUMA
                      SR11=SRCH(11,NCMD)/SUMA
                      WRITE(KW(38),587)ADRF,SR7,XTP(12),ADR2,&
                      SR8,ADYW,ADY1,ADY2,ADWS,ADNS,ADPS,ADKS,&
                      ADTS,ADAS,ADSS,ADFN,ADFP,ADIR,ADST,ADFU,ADDC,SR11,&
                      SR17,SR21,SR22,SR18,ADDP,SR9,ADWN,ADVO,ADDN,ADFX,&
                      SR12,ARFP,ADRP,ADLP,SR10,AYWP,ADPQ,ADPL,ADPY
                  END IF    
              END IF
              IF(KFL(14)>0)THEN
                  WRITE(KW(14),373)RWSA(1),RCHL(1),RCHS(1),RFPL(1),(SRCH(I,1)&
                  ,I=7,12)
                  WRITE(KW(14),373)SUMA,RCHL(MSA),RCHS(MSA),RFPL(MSA)&
                  ,(SRCH(I,NCMD),I=7,12)
              END IF
              X11=10.*SUMA
              SM1=SM1/SUMA
              SM2=SM2/SUMA
              SM3=SM3/SUMA
              SM4=SM4/SUMA
              SM40=SM40/SUMA
              SM81=.001*SM81/SUMA
              X2=100.*SRCH(12,NCMD)/(SRCH(7,NCMD)*SUMA+.1)
              TMAF=.001*TMAF
              TMPD=.001*TMPD
              ! PRINTOUT MANURE FILE
              XX=X11*XYR
              XX1=SUMA*XYR
              AD4=0.
              TSMU=0.
              TSMN=0.                                                                        
              TLMN=0.
              IF(TMAF>0..OR.TMPD>0.)THEN
                  IF(KFL(2)>0)THEN
                      SR9=SRCH(9,NCMD)/SUMA
                      SR8=SRCH(8,NCMD)/SUMA
                      SR10=SRCH(10,NCMD)/SUMA
                      SR11=SRCH(11,NCMD)/SUMA
                      SR12=SRCH(12,NCMD)/SUMA
                      SR17=SRCH(17,NCMD)/SUMA
                      SR18=SRCH(18,NCMD)/SUMA
                      WRITE(KW(2),2303)'TOTAL',SM3,SM2,NMC,SUMA,SR7,SR8,&
                      SR12,SR10,SR11,SR17,SR18,SM40,SR9,SM1,SM4,SM81,X2
                      WRITE(KW(2),'(//T6,A)')'OW ID #  AREA(ha)  MAP(t/ha/y)'
                  END IF    
                  DO IOW=1,NBON
                      IF(OWSA(IOW)<1.E-10)CYCLE
                      OMAP(IOW)=.001*OMAP(IOW)/(XYR*OWSA(IOW))
                      IF(KFL(2)>0)WRITE(KW(2),2308)IOW,OWSA(IOW),OMAP(IOW)
                      TSMU=TSMU+SMNU(IOW)
                      IF(IDFA(1,IOW)>0)THEN
                          ISA=IDFA(1,IOW)
                          KF=IDFT(2,ISA)
                          TSMN=TSMN+SMNU(IOW)*(FN(KF)+FNO(KF))
                          DO J=1,NBSL(ISA)
                              ISL=LID(J,ISA)
                              AD4=AD4+RSDM(ISL,ISA)*(FN(KF)+FNO(KF))*WSA(ISA)
                          END DO
                          KF=IDFT(1,ISA)
                          TLMN=TLMN+WTMU(ISA)*(FN(KF)+FNO(KF))
                      END IF
                  END DO
                  TSMN=1000.*TSMN
                  TOT=0.
                  AD1=0.
                  TWMP=0.
                  AD2=0. 
                  AD3=0.
                  DO ISA=1,MSA
                      SUM=0.
                      DO J=1,NBSL(ISA)
                          ISL=LID(J,ISA)
                          SUM=SUM+RSDM(ISL,ISA)
                      END DO
                      X2=SUM*WSA(ISA)
                      TOT=TOT+X2
                      AD1=AD1+SM(103,ISA)*WSA(ISA)
                      TWMP=TWMP+WTMU(ISA)
                      IF(IDFH(ISA)>0)THEN
                          AD2=AD2+X2
                          SMMU(ISA)=SMMU(ISA)*WSA(ISA)
                          AD3=AD3+SMMU(ISA)
                      END IF
                  END DO
                  AD1=AD1*XYR
                  X1=SUMA*XYR*SRCH(24,NCMD)
                  TWMP=.001*TWMP
                  SM6=.001*SM6
                  SM7=.001*SM7
                  SM8=.001*SM8
                  DDF1=TMPD+TWMB-TWMP-TSMU-TOT-X1-AD1
                  PER=100.*DDF1/(TMPD+1.E-10)
                  IF(KFL(2)>0)WRITE(KW(2),2304)PER,DDF1,TMPD,TWMB,TWMP,TSMU,TOT,X1,AD1
                  DF=TMPD-TMAF+TWMB-TWMP-TSMU-SYMU-SOFU-AD3-AD2
                  PER=100.*DF/(TMPD+1.E-5)
                  IF(KFL(2)>0)THEN
                      WRITE(KW(2),102)PER,DF,TMPD,TMAF,TWMB,TWMP,TSMU,SYMU,SOFU,AD3,AD2
                      WRITE(KW(2),103)SM6,SM7,SM8
                  END IF    
              END IF
              AD4=1000.*AD4
              TWB(1)=XYR*SRCH(9,NCMD)
              TWB(2)=XYA*XTP(21)
              DF=TBTN-SM134+SRQN-SM42+SM43-SM46+SM53+SM54+SM55-SM82-SM96+SM105-TWB(1)&
              -TWB(2)-AD5-SFTN-TLMN-TSMN-TRMN-TRON-AD4
              PER=100.*DF/SFTN
              IF(KFL(1)>0)THEN
                  WRITE(KW(1),'(//T10,A/)')'TOTAL N BALANCE (kg)'
                  WRITE(KW(1),2314)PER,DF,TBTN,SRQN,SM42,SM43,SM46,SM53,SM54,SM55,SM82,&
                  SM96,SM105,SM134,(TWB(I),I=1,2),AD5,SFTN,TLMN,TSMN,TRMN,TRON,AD4
                  IF(IHY>0)WRITE(KW(1),2299)SQVL(NCMD),SHYD(NCMD)
              END IF    
              TWB(1)=SM9
              TWB(2)=SMH34
              XY1=10.*XYR
              TWB(3)=XYR*SM71
              TWB(4)=XYR*SM155
              TWB(5)=XYR*SM18
              TWB(6)=XYR*SM19
	          TWB(7)=XYR*SM104
              TWB(8)=SM147*XYR
              TWB(9)=SM146*XYR
              TWB(10)=SM156*XYR
              TWB(11)=XYR*SM110
              DF=BSWB+BGWSB+BRSVB+BSNOB+BSWLB+TWB(1)-TWB(2)-TWB(3)-TWB(4)&
              +TWB(5)+TWB(6)+TWB(7)-TWB(8)-TWB(9)-TWB(10)-TWB(11)-FSWB-&
              FGWSB-FRSVB-FSNOB-FSWLB
              PER=100.*DF/FSWB
              IF(KFL(1)>0)THEN
                  WRITE(KW(1),'(//T10,A/)')'TOTAL WATER BALANCE (m3)'
                  WRITE(KW(1),2312)PER,DF,BSWB,BGWSB,BRSVB,BSNOB,BSWLB,TWB,&
                  FSWB,FGWSB,FRSVB,FSNOB,FSWLB
              END IF    
              TWB(1)=XYA*XTP(3) 
              TWB(2)=XYA*XTP(13)
	          DF=SRYB+SDEG-SDEP+TWB(1)-TWB(2)-SRYF
              PER=100.*DF/(TWB(2)+1.E-10)
              IF(KFL(1)>0)THEN
	              WRITE(KW(1),'(//T10,A/)')'TOTAL SEDIMENT BALANCE'
                  WRITE(KW(1),2315)PER,DF,SRYB,TWB(1),TWB(2),SDEP,SDEG,SRYF
              END IF    
              ! 1 NX(1) = NUMBER OF Y FOR SECOND THRU LAST SIMULATION.
              ! 2 NX(2) = 0 FOR NORMAL EROSION OF SOIL PROFILE
              !         = 1 FOR STATIC SOIL PROFILE
              ! 3 NX(3) = N0 FOR ANNUAL WATERSHED OUTPUT
              !         = N1 FOR ANNUAL OUTPUT                       | N YEAR INTERVAL
              !         = N2 FOR ANNUAL WITH SOIL TABLE              | N=0 SAME AS
              !         = N3 FOR MONTHLY                             | N=1 EXCEPT
              !         = N4 FOR MONTHLY WITH SOIL TABLE             | N=0 PRINTS
              !         = N5 FOR MONTHLY WITH SOIL TABLE AT HARVEST  | OPERATIONS
              !         = N6 FOR N DAY INTERVAL
              !         = N7 FOR SOIL TABLE ONLY N DAY INTERVAL
              !         = N8 FOR SOIL TABLE ONLY DURING GROWING SEASON N DAY INTERVAL
              !         = N9 FOR N DAY INTERVAL DURING GROWING SEASON
              ! 4 NX(4) = ID NUMBER OF WEATHER VARIABLES INPUT.  RAIN=1,  TEMP=2,
              !           RAD=3,  WIND SPEED=4,  REL HUM=5.  IF ANY VARIABLES ARE INP
              !           RAIN MUST BE INCLUDED.  THUS, IT IS NOT NECESSARY TO SPECIF
              !           ID=1 UNLESS RAIN IS THE ONLY INPUT VARIABLE.
              !           LEAVE BLANK IF ALL VARIABLES ARE GENERATED.  EXAMPLES
              !           NGN=1 INPUTS RAIN.
              !           NGN=23 INPUTS RAIN, TEMP, AND RAD.
              !           NGN=2345 INPUTS ALL 5 VARIABLES.
              ! 5 NX(5) = TURNS ON .SAO FILE.
              ! 6 NX(6) = TURNS ON .RCH FILE.
              NN=NGN
              READ(KR(6),3100)(NX(I),I=1,7)
              IF(NX(1)==0)EXIT
              NBYR=NX(1)
              ISTA=NX(2)
              IPD=NX(3)
              NGN=NX(4)
	          KFL(34)=NX(5)
	          IF(KFL(34)>0)OPEN(KW(34),FILE=ASTN//".SAO")
	          KFL(35)=NX(6)
	          IF(KFL(35)>0)THEN
                  OPEN(KW(35),FILE=ASTN//".RCH")
	              SMMH=0.
              END IF
	          CALL AISPL(IPD,INP)
              IF(INP>0)THEN
                  NOP=0
              ELSE
                  NOP=1
                  INP=1
              END IF
              IF(IPD<=5)IPYI=INP
              IF(NGN>0)THEN
                  NDWT=1
                  DO ISA=1,MSA
	                  CALL WDLYSTA
	              END DO
              END IF
              IPY=1
              DO K=1,100
                  NX(K)=K
              END DO
              CALL ARESET
          END DO    ! OUTPUT LOOP
	      REWIND KR(6)
          IF(NGN0>0)THEN
              M=NDWT-1
              DO L=1,M
                  REWIND KRST(L)
              END DO
          END IF
          IGN=IGN+100
          IF(IGN<IGMX*100)THEN
              REWIND KR(1)
              REWIND KR(5)
              GO TO 2140
          ELSE
              IGN=0
          END IF
          CLOSE(KR(1))
          CLOSE(KR(5))
          DO I=2,SIZE(KW)
              CLOSE(KW(I))
          END DO
          CALL ADEALLOCATE
      END DO  ! RUN LOOP
      DO I=1,SIZE(KR)
          CLOSE(KR(I))
      END DO
      ALLOCATE (KW(2*MSA+MSO))
      KW(1)=101
      CALL ATIMER(1)
      CLOSE(KW(1))

      !rtb MODFLOW: deallocate arrays
      if(IMF>0) then
        call amrt_close
        call mf_close
        call rt_close
      endif
      STOP
    3 FORMAT(1X,A4,12F9.0,11X,A4)
   12 FORMAT(T10,'CHANNEL GEOMETRY PARMS'/T15,'WIDTH  = ',2F8.4/T15,&
      'DEPTH  = ',2F8.4/T15,'LENGTH = ',2F8.4)
   16 FORMAT(T23,'--------ID NUMBERS---------'/T10,'COMMAND',6X,'CMD',&
      5X,'OUT',5X,'IN1',5X,'SA#',5X,'IN2',5X,'SA#')
   17 FORMAT(T10,A8,3I8,8X,2I8)
   18 FORMAT(//1X,'______________WATERSHED SUMMARY TABLE________________&
      _'/T10,'AVE ANNUAL SUM OF SUBAREA OUTFLOWS/TOTAL WATERSHED OUTFLOW'/)
   19 FORMAT(11F10.0)   
   24 FORMAT(8X,A2,I8,19X,F12.2,I4,F8.2,F8.3,F8.2,3F8.1)
   33 FORMAT(8X,A2,I8,19X,F12.2,5F8.1,20F8.2)
   35 FORMAT(T89,'_____QRB______'/T16,'OUT',12X,'SA',11X,'WSA',3X,&
      '# FLOWS',2X,'TC',5X,'DR',7X,'ER',6X,'Q',6X,'MEAN',4X,'MAX',5X,&
      '# CH',3X,'VCHA',4X,'VCHB',5X,'# FP',3X,'VFPA',4X,'VFPB'/T8,'CMD',&
      5X,'ID#',9X,'#',7X,'ID',6X,'(ha)',4X,'/y',3X,'(h)',T82,'(mm)',3X,&
      '|---(mm/h)---|',3X,'FLOWS',2X,'|----(m/s)---|',3X,'FLOWS',2X,&
      '|----(m/s)---|')
   36 FORMAT(T16,'OUT',12X,'SA',12X,'WSA',6X,'Q',6X,'SSF',5X,'QRF',5X,&
      'QDR',5X,'RTF',4X,'PIPE',6X,'Y',6X,'DEP',5X,'DEG',4X,'YMNU',6X,&
      'YC',6X,'YN',6X,'QN',6X,'YP',6X,'QP',4X,'SSFN',4X,'QRFN',4X,'QDRN',&
      4X,'RTFN'/T8,'CMD',5X,'ID#',9X,'#',7X,'ID',7X,'(ha)',4X,'|--------&
      ------------(mm)-------------------|',2X,'|----------------(t/ha)-&
      --------------|',1X,'|------------------------(kg/ha)-------------&
      ----------------|')
   39 FORMAT(T10,'SUBAREA CONTAINS LAGOON'/T15,'DRAINAGE AREA = ',F7.2,&
      ' ha'/T15,'NORMAL VOL = ',E12.4,' mm'/T15,'MAX VOL = ',E12.4,' mm'&
      /T15,'DRAW DOWN TIME = ',F6.2,' d'/T15,'WASH WATER = ',F6.3,&
      ' m3/hd/d'/T15,'DESIGN SAFETY FACTOR =',F6.3)
  101 FORMAT(13X,I10)
  102 FORMAT(/5X,'MANURE APPLICATION BALANCE (T)'/5X,'PER =',E13.6,2X,&
      'DF  =',E13.6,2X,'TWMP=',E13.6,2X,'TMAP=',E13.6,2X,'TWMB=',E13.6/&
      5X,'TWMF=',E13.6,2X,'TSMU=',E13.6,2X,'YMFA=',E13.6,2X,'LGOF=',E13.6,&
      2X,'MNFA=',E13.6/5X,'RSFA=',E13.6)
  103 FORMAT(/5X,'LIQUID MANURE APPLIED (T)=',E13.6//5X,'COMMERCIAL FERT &
      APPLIED (T)'/10X,'P=',E13.6/10X,'N=',E13.6)
  104 FORMAT(//8X,'Y',3X,'M',3X,'D',1X,'ISA',1X,'ID#',1X,'NBPQ',12X,'RUN NAME',&
      26X,'WSAha',6X,'QPMX',6X,'QPMN',6X,'QPSD',20(6X,A4))
  105 FORMAT(5X,5I4,I5,1X,A40,F10.3,3F10.2,F10.0,2F10.1,F10.0,10F10.3)      
  140 FORMAT(10F10.0,I10)                     
  316 FORMAT(/T10,'TILLAGE OPERATIONS'/6X,'DATE',T21,'OP',3X,'COST',6X,'MX',&
      6X,'RR',6X,'DP',6X,'FR',5X,'RHT',5X,'RIN',5X,'DKH',5X,'DKI',6X,'HV',&
      6X,'HV',10X,'NRCS',5X,'IRR',6X,'Q/'/6X,'Y M D',3X,'NAME',2X,'CD',&
      2X,'($/ha)',5X,'EF',5X,'(mm)',4X,'(mm)',4X,'COMP',4X,'(mm)',5X,'(m)'&
      4X,'(mm)',5X,'(m)',5X,'EF',5X,'IDX',4X,'CROP',3X,'CN',4X,'TRGR',4X,&
      'VIRR',4X,'HUSC')
  317 FORMAT(5X,3I2,1X,A8,I2,2F8.2,2F8.0,F8.3,F8.0,F8.2,F8.0,2F8.2,F8.3&
      ,5X,F8.1,2F8.2,F8.3)
  326 FORMAT(/2X,A4,1X,'YLD=',F5.1,'/',F5.1,2X,'BIOM=',F5.1,'t/ha',2X,&              
      'YLN=',F5.0,2X,'YLP=',F5.0,2X,'YLK=',F5.0,2X,'FN=',F5.0,2X,'FP=',&
      F5.0,2X,'FK=',F5.0,'kg/ha'/T7,'IRGA=',F5.0,2X,'CAW=',F6.0,'mm',&
      2X,'WUEF=',F7.2,'kg/mm',2X,'RAD=',F7.0,'Mj/m2',2X,'HU=',F7.0,2X,&
      'PSTF=',F5.2)                           
  328 FORMAT(/T10,'FERTILIZER APPLIED'/T11,'DATE',T31,'FT',4X,'EQ',4X,&               
      'TR',4X,'COST',6X,'WT',4X,'DPTH',2X,'|------------FRACTION OF WT--&            
      ----------|'/T9,'Y   M   D',2X,'NAME',7X,'NO',4X,'NO',4X,'NO',4X,&                
      '$/ha',3X,'kg/ha',6X,'mm'6X,'MN',5X,'NH3',6X,'ON',6X,'MP',6X,'OP',&             
      4X,'HUSC')                                                                 
  329 FORMAT(5X,3I4,1X,A8,3I6,F8.2,2F8.0,6F8.3)                                           
  330 FORMAT(T10,'COSTS'/T15,'LIME = ',F7.2,' $/t'/T15,'IRR WATER = ',&
      F7.2,' $/mm'/T15,'FUEL PRICE = ',F7.2,' $/l'/T15,'WAGE PRICE =',&
      F7.2,' $/h')
  373 FORMAT(10E13.5)
  396 FORMAT(8F10.0)
  401 FORMAT(T10,'AUTO HYDRAULICS DESIGN'/T15,'2 Y FREQ 24 H RAINFALL ='&
      ,F6.0,' mm'/T15,'RUNOFF VOL = ',F5.1,' mm'/T15,'CH CAP-WSA EXP = '&
      ,F6.3,/T15,'BASE CH LENGTH = ',F6.2,' km'/T15,'BASE CH SLOPE = ',&
      F8.5,' m/m'/T15,'BASE TC = ',F6.3,' h'/T15,'BASE PEAK FLOW = ',&
      F6.2,' mm/h'/T15,'CH BW/D = ',F4.1/T15,'FP WIDTH/CH WIDTH = ',F6.1&
      /T10,'FLOODPLAIN SAT COND = ',F7.4,' mm/h'/T10,'MAX GROUNDWATER ST&
      ORAGE = ',F5.0,' mm'/T10,'GROUND WATER RESIDENCE TIME = ',F5.0,&
      ' d'/T10,'RF/(RF+DPRK) = ',F6.3/T10,'FRACTION SUBAREA ABOVE PONDS = '&
      ,F6.3/T10,'REACH CH C FACTOR = ',F6.3)
  412 FORMAT(T10,'LAGOON PUMP DATE = ',I4/T10,'SOLID MANURE SCRAPE &
      INTERVAL = ',I4,' d')
  417 FORMAT(/1X,'-----TIME OF CONCENTRATION STATS(h)'/T10,'MEAN = ',&
      F6.2,5X,'MIN = ',F6.2,5X,'MAX = ',F6.2)
  448 FORMAT(/1X,'-----SEDIMENT OF CONCENTRATION STATS(g/m3)'/T10,'MEAN&
      = ',F10.0,5X,'MAX = ',F10.0,5X,'STDV = ',F10.0)
  460 FORMAT(/1X,'-----FREQUENCY & DURATION OF PESTICIDE LOSS(g/ha)'/T35&
      ,'DURATION(d)'/20X,'1',12X,'4',11X,'21',11X,'60',11X,'90')
  461 FORMAT(///1X,'______________WATERSHED SUMMARY TABLE_____________')     
  470 FORMAT(5X,A4,10E12.4)
  472 FORMAT('REACH',I5,I9,I6,50E12.4)
  498 FORMAT(4X,I8,1X,I8,51F10.2)
  508 FORMAT(1X,A5,4X,A80)      
  509 FORMAT(10X,A80)
  526 FORMAT(T10,A8,2X,3I4)
  576 FORMAT(1X,A4,11F10.0,A5)
  577 FORMAT(T7,'STRESS (BIOM) WATER=',F5.1,2X,'N=',F5.1,2X,'P=',F5.1,2X&            
      ,'K=',F5.1,2X,'TEMP=',F5.1,2X,'AIR=',F5.1,2X,'SALT=',F5.1,5X,&                 
      '(ROOT) BD=',F5.1,2X,'ALSAT=',F5.1,2X,'TEMP=',F5.1,'D')                        
  583 FORMAT(T74,'COTL',6X,'COOP',6X,'MTCO',6X,'MASS',6X,'FUEL'/6X,'SA#'&
      ,6X,'ID',4X,'Y M D',5X,'OP',14X,'CROP',2X,'MT#  HC',2X,'EQ  TR',2X&
      ,'|----------($/ha)-----------|',2X,'(kg/ha)',4X,'(l/ha)')
  584 FORMAT(6X,'QS',9X,'Y',8X,'QN',6X,'SSQN',8X,'YN',6X,'DWOC'/5X,&
      '(mm)',6X,'(t/ha)',2X,'|----------(kg/ha)----------|',3X,'(t/ha)')
  585 FORMAT(6X,'SA#',4X,'SAID',6x,'RFmm'7X,'Qmm',4X,'WYLDmm',2X,'RUS2t/ha',&
      5X,'Yt/ha',2X,'YWNDt/ha',2X,'YLDGt/ha',2X,'YLDFt/ha',7X,'WSd',7X,'NSd',&
      7X,'PSd',7X,'KSd',7X,'TSd',7X,'ASd',7X,'SSd',3X,'FNkg/ha',3x,&
      'FPkg/ha',4X,'IRGAmm',6X,'STIR',2X,'FULUl/ha',2X,'DWOCt/ha',3X,&
      'QNkg/ha',1X,'SSFNkg/ha',1X,'QRFNkg/ha',1X,'QDRNkg/ha',1X,&
      'RTFNkg/ha',1X,'DPKNkg/ha',3X,'YNkg/ha',1X,'YNWNkg/ha',1X,&
      'NVOLkg/ha',1X,'DNITkg/ha',1X,'NFIXkg/ha',3X,'QPkg/ha',1X,&
      'SSFPkg/ha',1X,'QDRPkg/ha',1X,'PRKPkg/ha',3X,'YPkg/ha',1X,&
      'YPWNkg/ha',2X,'QPSTg/ha',2X,'LPSTg/ha',2X,'YPSTg/ha',7X,'Qmm',&
      4X,'WYLDmm',5X,'Yt/ha',2X,'TQNkg/ha',3X,'YNkg/ha',3X,'QPkg/ha',&
      3X,'YPkg/ha',2X,'QPSTg/ha',2X,'YPSTg/ha')
  586 FORMAT(1X,2I8,100F10.3)      
  587 FORMAT(8X,'WATERSHED',100F10.3)
  602 FORMAT(19X,'WPM1 = ',A20,2X,'FRACT WSA = ',F6.3/T11,'JAN',6X,'FEB'&
      ,6X,'MAR',6X,'APR',6X,'MAY',6X,'JUN',6X,'JUL',6X,'AUG',6X,'SEP',6X&
      ,'OCT',6X,'NOV',6X,'DEC',6X,' YR')
  603 FORMAT(T11,'JAN',9X,'FEB',9X,'MAR',9X,'APR',9X,'MAY',9X,'JUN',9X,&
      'JUL',9X,'AUG',9X,'SEP',9X,'OCT',9X,'NOV',9X,'DEC',9X,' YR')
  635 FORMAT(/10X,'SA#= ',I8,1X,'ID= ',I8/T35,'SOIL LAYER NO'/2X,11(6X,&
      A4))
  637 FORMAT(5X,'ELSC=',E13.6,2X,'ELMC=',E13.6,2X,'EBMC=',E13.6,2X,&
      'EHSC=',E13.6,2X,'EHPC=',E13.6)
  639 FORMAT(5X,'ESK =',E13.6,2X,'EEK =',E13.6,2X,'EFK =',E13.6,2X,&
      'ESDK=',E13.6/5X,'ESOK=',E13.6,2X,'EUK1=',E13.6)                    
  646 FORMAT(5X,2I4,2I10,4F10.3)
  649 FORMAT(1X,A4,11F10.3,A5)
  677 FORMAT(/T10,'AUTO IRR EQUIP  = ',A8,2X,'DEPTH = ',F6.3,' m')
  689 FORMAT(T10,'ANIMAL FERT DUMP = ',A8,2X,'DEPTH = ',F6.3,' m',2X,&
      'FERT = ',A8)
  690 FORMAT(T10,'AUTO STOCK PILE MANURE = ',A8,2X,'DEPTH = ',F6.3,' m',&
      2X,'FERT = ',A8)
  691 FORMAT(T10,'AUTO LIQUID MANURE = ',A8,2X,'DEPTH = ',F6.3,' m',&
      2X,'FERT = ',A8/)              
  713 FORMAT(/1X,'-----LIVESTOCK BUY SELL DATA'/T27,'HERD'/T12,'OWNER',&
      1X,'HERD',5X,'SIZE'/6X,'Y M D',3X,'ID  ID',5X,'(HEAD)')
  714 FORMAT(///1X,'____________________LIVESTOCK MANAGEMENT DATA_______&
      _____________'//T20,'HERD',13X,'FRACTION',3X,'GRAZE',4X,'MANURE',&
      5X,'URINE'/4X,'OWNER',1X,'HERD',5X,'SIZE',5X,'MANURE',3X,'IN FEED'&
      ,4X,'RATE',6X,'PROD',6X,'PROD',/7X,'ID  ID',5X,'(HEAD)',7X,'ID',5X&
      ,'AREA',T46,'(kg/hd/d)',1X,'(kg/hd/d)',1X,'(l/hd/d)')
  716 FORMAT(T10,'AUTO FERT EQUIP = ',A8,2X,'DEPTH = ',F6.3,' m',2X,&
      'FERT = ',A8)
  729 FORMAT(3X,'SA# SAID   YR  YR#',5X,'Q',5X,'SSF',5X,'PRK',4X,'QDRN',&
      7X,'Y',5X,'YOC',5X,'PSTN',14X,'PAPL',9X,'PSRO',9X,'PLCH',9X,'PSSF'&
      ,9X,'PSED',9X,'PDGF',9X,'PDGS',9X,'PDRN',9X,'PRSF'/23X,'|---------&
      ---(mm)------------|  (t/ha)',2X,'(t/ha)',17X,'|-----------------&
      --------------------------------(g/ha)----------------------------&
      --------------------|')          
  811 FORMAT(10I8)
  848 FORMAT(5X,3I2,1X,2I4,I10)
  854 FORMAT(/10X,A80)
  856 FORMAT(7X,A20,2X,A4)
  857 FORMAT(I9,1X,'AVEAN',4X,4F10.2,10X,5(A10,50X,7F10.2,20X,2F10.2))
  859 FORMAT(1X,'AVEAN',I6,I9,I5,1X,E9.3,6F10.2,F10.0,13F10.2,5F10.1,&
      9F10.2,F10.1,60X,6(2F10.2,90X,7F10.2,A10))
  900 FORMAT(T49,'_______WATER CONTENT________   ___SAT COND__    __BULK &
      DEN___   ___________TEXTURE__________  CROP',12X,'SUM',13X,'AL',&
      12X,'P SORP',3X,'_______MIN P_______',5X,'ORG',5X,'NO3',5X,'ORG',&
      5X,'ORG'/8X,'SUBAREA',T31,'LAY',4X,'DEPTH',4X,'POR',5X,'FC',6X,&
      'WP',4X,'CURNT',4X,'VERT',4X,'HORZ',4X,'33KPA  OV DRY',4X,'SAND',&
      4X,'SILT',4X,'CLAY',4X,'ROCK',4X,'RSD'5X,'PH',6X,'BASE',4X,'CEC',&
      4X,'SAT',5X,'CACO3',4X,'RTO',4X,'LAB',5X,'ACT',5X,'STB',6X,'P',7X,&
      'N',7X,'N',7X,'C'/7X,'#',6X,'ID',3X,'YR',1X,'YR#',2X,'MO',2X,'NO',&
      5X,'(m)',4X,'------------(m/m)-----------',3X,'---(mm/h)----',5X,&
      '--(t/m3)---',5X,'------------(%)-------------',2X,'(t/ha)',11X,&
      '-(cmol/kg)--',4X,'-----(%)-----',11X,&
      '---------------------(g/t)------------------',4X,'(%)')           
  934 FORMAT(10X,'WATERSHED AREA = ',F10.2,' ha',T84,20(A16,6X))
  935 FORMAT(10X,'WATERSHED AREA = ',F10.2,' ha',T126,20(A16,18X))
  937 FORMAT(2X,'JDA   YR',T19,'WYLDmm',11X,'Yt/ha',11X,'YNkg/ha',10X,&
      'YPkg/ha',10X,'QNkg/ha',10X,'QPkg/ha',1X,20(9X,'QPSTg/ha',9X,&
      'YPSTg/ha'))
  938 FORMAT(3X,'YR',3X,'MO',6X,'WYLDmm',6X,'Yt/ha',4X,'YNkg/ha',4X,&
      'YPkg/ha',4X,'QNkg/ha',4X,'QPkg/ha',20(3X,'QPSTg/ha',3X,'YPSTg/ha'))
  939 FORMAT(2X,'JDA   YR',T19,'WYLDmm',11X,'Yppm',12X,'YNppm',12X,&
      'YPppm',12X,'QNppm',12X,'QPppm')    
1001 FORMAT(20X,'1',11X,'2',11X,'3',11X,'4',11X,'5',11X,'6',11X,'7',&
      11X,'8',11X,'9',10X,'10',10X,'11',10X,'12',10X,'13',10X,'14',10X,&
      '15',10X,'16',10X,'17',10X,'18',10X,'19',10X,'20',10X,'21',10X,'22',&
      10X,'23',10X,'24',10X,'25',10X,'26',10X,'27',10X,'28',10X,'29',10X,&
      '30',10X,'31',10X,'32',10X,'33',10X,'34',10X,'35',10X,'36',10X,'37',&
      10X,'38',10X,'39',10X,'40',10X,'41',10X,'42',10X,'43',10X,'44',10X,&
      '45',10X,'46',10X,'47',10X,'48',10X,'49',10X,'50',10X,'51',10X,'52',&
      10X,'53',10X,'54'/7X,'YR',6X,'PRCPmm',8X,'ETmm',9X,'Qmm',7X,'SSFmm',&
      6X,'RSSFmm',7X,'QRFmm',7X,'QDRmm',7X,'PRKmm',6X,'IRGAmm',&
      8X,'WYmm',7X,'Yt/ha',4X,'YWNDt/ha',5X,'QNkg/ha',3X,'SSFNkg/ha',3X,&
      'QRFNkg/ha',3X,'RSFNkg/ha',5X,'YNkg/ha',3X,'YNWNkg/ha',3X,&
      'QDRNkg/ha',3X,'PRKNkg/ha',5X,'DNkg/ha',3X,'AVOLkg/ha',3X,&
      'NFIXkg/ha',4X,'FNOkg/ha',3X,'FNMNkg/ha',3X,'FNMAkg/ha',5X,&
      'QPkg/ha',5X,'YPkg/ha',3X,'YPWNkg/ha',3X,'PRKPkg/ha',4X,'FPOkg/ha',&
      4X,'FPLkg/ha',5X,'QCkg/ha',5X,'YCkg/ha',3X,'YCWNkg/ha',&
      4X,'RFNkg/ha',4X,'YLNkg/ha',4X,'YLPkg/ha',4X,'BTNkg/ha',4X,'BTPkg/ha',&
      4X,'FTNkg/ha',4X,'FTPkg/ha',4X,'BTCkg/ha',4X,'FTCkg/ha',3X,&
      'BPDPkg/ha',3X,'FPDPkg/ha',3X,'BSLTkg/ha',3X,'FSLTkg/ha',3X,&
      'BTC1kg/ha',3X,'FTC1kg/ha',2X,'RUS2A1t/ha',7X,'YTHSd',7X,'YWTHd',3X,&
      'QDRPkg/ha')
 2072 FORMAT(T31,'RFV',7X,'Q',5X,'SSF',5X,'YSD',10(23X,'PSRO',8X,'PSSF',&
      8X,'PSED',4X)/6X'SA#',6X,'ID',4X,'Y M D',2X,'|----------(mm)------&
      --|',1X,'(t/ha) ',15(3X,'PSTN',11X,'|--------------(g/ha)---------&
      -----|',1X))     
 2073 FORMAT(/T5,'APEX1501 v20181201',2X,3I4,2X,2(I2,':'),I2)    ! Added date to the version - Luca Doro
 2111 FORMAT(/T52,'_____________CHANNEL_____________   ________UPLAND___&
      ______',T172,'INITIAL COND'/14X,'SA',9X,'WSA',6X,'XCT',6X,'YCT',6X&
      ,'L',8X,'d',8X,'S',8X,'N',8X,'S',8X,'SL',7X,'N',8X,'SL',7X,'TC',&
      16X,'GWST',5X,'GWMX',3X,'GWRESTM',6X,'SNO',6X,'STD',T207,'BFFL',&
      T270,'SATK'/11X,'#',7X,'ID',4X,'(ha)',5X,'(km)',5X,'(km)',5X,&
      '(km)',5X,'(m)',5X,'(m/m)',13X,'(m/m)',5X,'(m)',T116,'FACT',5X,&
      '(h)',4X,'RF/DPRK',5X,'(mm)',5X,'(mm)',6X,'(d)',5X,'(mm)',5X,'(t)',&
      5X,'URBF',5X,'BCOF',7X,'(m)',6X,'CN2',T225,'SOIL',T244,'OPSC',&
      T263,'NVCN',3X,'(mm/h)',T282,'FWTH')
 2112 FORMAT(5X,'ISA',4X,'NBSA',4X,'Y   M   D',8X,'CN',7X,'SCI',5X,'RFVmm',&
      4X,'STMP2c',5X,'SMLmm',7X,'Qmm',5X,'SSFmm',5X,'QRFmm',4X,'RSSFmm',&
      4X,'WYLDmm',3X,'QRBmm/h',7X,'TCh',6X,'DURh',6X,'ALTC',7X,'AL5',3X,&
      'REPmm/h',4X,'RZSWmm',4X,'GWSTmm')           
 2122 FORMAT(T7,A4,100F10.2)
 2123 FORMAT(T7,A4,100F10.0)
 2124 FORMAT(T10,A80,F8.4)      
 2126 FORMAT(/T43,'________________________________CHANNEL______________&
      ________________    ______FLOODPLAIN_______'/14X,'SA',9X,'WSA',6X,&
      'RWSA',7X,'L',8X,'d',7X,'S',8X,'BW',7X,'TW',8X,'N',7X,'C',9X,'K',&
      8X,'W',8X,'L',7X,'S'/11X,'#',7X,'ID ',3X,'(ha)',5X,'(ha)',5X,&
      '(km)',6X,'(m)',4X,'(m/m)',6X,'(m)',6X,'(m)',33X,'(m)',5X,'(km)',&
      4X,'(m/m)')
 2128 FORMAT(/T35,'_EMERG SPLWAY_     _____PRINC SPLWAY_____'/T36,'SURF'&
      ,T54,'SURF',T71,'FLOOD',4X,'_____INIT_____',4X,'NORM',5X,'HYD',6X,&
      'SETTL'/17X,'SA',7X,'WSA',6X,'AREA',6X,'VOL',5X,'AREA',6X,'VOL',3X,&
      'ST REL',6X,'VOL',4X,'SED C    SED C',3X,'COND',7X,'TIME'/10X,'#',8X&
      ,'ID',5X,'(ha)',5X,'(ha)',5X,'(mm)',5X,'(ha)',5X,'(mm)',5X,&
      '(d)',6X,'(mm)',4X,'____(g/m3)____',2X,'(mm/h)',7X,'(d)')
 2129 FORMAT(//6X,'SAID',2X,'YR',2X,'MO',4X,4(1X,2A4,1X),A10,&
      5(1X,A8,1X,11(1X,2A4,1X),4(A8,2X)))
 2142 FORMAT(//6X,'SAID',5X,'GIS',2X,'TIME',2X,'WSAha',4X,34(1X,2A4,1X),&
      ' WOCt/ha  ','PCTI200um',1X,'PCTI10um',2X,'PCTI2um',3X,'PCTO200um',1X,&
      'PCTO10um',2X,'PCTO2um',2X,5(2X,'YLDGt/ha',2X,'YLDFt/ha',17&
      (2X,2A4),5X,'CPNM-'))
 2143 FORMAT(5X,'RUN # = ',I4,2X,'SUBAREA FILE = ',A80)
 2151 FORMAT(T10,A8,5I8)
 2157 FORMAT(1X,6F10.3)
 2162 FORMAT(//6X,'RCID',6X,'GIS',2X,'TIME',3X,'WSAha',5X,35(1X,2A4,3X),&
      1X,'PCTI200um',3X,'PCTI10um',4X,'PCTI2um',5X,'PCTO200um',3X,&
      'PCTO10um',4X,'PCTO2um')
 2232 FORMAT(//T10,'SA#=',I8,1X,'ID=',I8)
 2281 FORMAT(10X,3F10.2)
 2282 FORMAT(1X,2I8,I4,1X,A4,2F5.1,3F5.0,F10.3,F6.0,8F8.3,2F6.0,F6.1,&
      2F6.0,F6.0)
 2284 FORMAT(/10X,'AVE ANNUAL VALUES'/11X,'SA',6X,'OWN',5X,'YLD1',1X,&
      'YLD2',1X,'YLN',2X,'YLP',2X,'COW',5X,'WSA',5X,'Q',7X,'Y',6X,'QP',&
      6X,'YP',6X,'QN',4X,'SSFN',4X,'RSFN',4X,'PRKN',6X,'YN',6X,'FP',4X,&
      'FN',3X,'MAP',3X,'AP0',3X,'APF',3X,'CSP'/8X,'#',6X,'ID',3X,&
      '# CROP',1X,'|--t/ha-| |-kg/ha-|',2X,'HD',6X,'HA',6X,'mm',4X,&
      't/ha',3X,'|-----------------------------kg/ha--------------------&
      ----------|',2X,'t/ha',1X,'|---g/t---|',2X,'g/m3')
 2299 FORMAT(//T10,'FLOOD ROUTING SUMMARY WATERSHED OUTLET'/T10,&
      'SUM QVOL = ',E16.6,' mm',2X,'SUM HVOL = ',E16.6,' mm')
 2302 FORMAT(8X,A2,I8,2X,I8,1X,I8,F12.2,6F8.1,20F8.2)
 2303 FORMAT(4X,A5,27X,2F5.0,I5,30E16.6)
!2303 FORMAT(4X,A5,27X,2F5.0,I5,F10.3,F6.0,8F8.3,2F6.0,F6.1,12X,3F6.2)      
 2304 FORMAT(//5X,'MANURE BALANCE (T)'/5X,'PER =',E13.6,2X,'DF  =',E13.6,&
      2X,'TMPD=',E13.6,2X,'TWMB=',E13.6,2X,'TWMF=',E13.6/5X,'TSMU=',E13.6,&
      2X,'RSDM=',E13.6,2X,'YMNU=',E13.6,2X,'MNMU=',E13.6)            
 2308 FORMAT(7X,I4,2F10.2)
 2309 FORMAT(1X,4F10.4,6F10.3)
 2312 FORMAT(5X,'PER =',D13.6,2X,'DF  =',D13.6,2X,'BSW =',D13.6,2X,&
      'BGWS=',D13.6,2X,'BRSV=',D13.6,2X,'BSNO=',D13.6/5X,'BSWL=',D13.6,&
      2X,'PCP =',D13.6,2X,'WYLD=',D13.6,2X,'DPRK=',D13.6,2X,'ET  =',&
      D13.6,2X,'IRG =',D13.6/5X,'QIN =',D13.6,2X,'PSOQ=',D13.6,2X,&
      'EVRT=',D13.6,2X,'RSIR=',E13.6,2X,'WLIR=',E13.6,2X,'IRDL=',D13.6/&
      5X,'FSW =',D13.6,2X,'FGWS=',D13.6,2X,'FRSV=',D13.6,2X,'FSNO=',&
      D13.6,2X,'FSWL=',D13.6)
 2313 FORMAT(3X,'WSA(ha)   CHL(km)  CHS(m/m)   FPL(km)     Q(mm)',5X,&
      'Y(t/ha) YN(kg/ha) YP(kg/ha) QN(kg/ha) QP(kg/ha)')
 2314 FORMAT(5X,'PER =',D13.6,2X,'DF  =',D13.6,2X,'BSN =',D13.6,2X,&
      'PCP =',D13.6,2X,'DN  =',D13.6,2X,'NFIX=',D13.6/5X,'VOL =',D13.6,2X,&
      'FNO =',D13.6,2X,'FNMN=',D13.6,2X,'FNMA=',D13.6,2X,'BURN=',D13.6,2X,&
      'DPKN=',D13.6/5X,'PSON=',D13.6,2X,'YNWN=',D13.6,2X,'YN  =',D13.6,2X,&
      'WYLN=',D13.6,2X,'YLN =',D13.6,2X,'FSN =',D13.6/5X,'TLMN=',D13.6,2X,&
      'TSMN=',D13.6,2X,'TRMN=',D13.6,2X,'TRON=',D13.6,2X,'RSDN=',D13.6)                
 2315 FORMAT(5X,'PER =',D13.6,2X,'DF  =',D13.6,2X,'BRSY=',D13.6,2X,&
      'YS  =',D13.6,2X,'YW  =',D13.6,2X,'DEP =',D13.6/5X,'DEG =',D13.6,&
      2X,'RSYF=',D13.6/)
 2320 FORMAT(I4,6F8.0)
 2410 FORMAT(4X,'1',6X,'2',6X,'3',6X,'4',6X,'5',6X,'6',6X,'7',6X,'8',&
      6X,'9',5X,'10',5X,'11',5X,'12',5X,'13',5X,'14',5X,'15',5X,'16',5X,&
      '17',5X,'18',5X,'19',5X,'20')
 2420 FORMAT(5X,3I2,F8.1,F8.2)
 2450 FORMAT(T10,'PEAK RATE-EI ADJ FACTOR = ',F6.3/T10,'AVE N CONC IN RA&
      INFALL = ',F5.2,' ppm'/T10,'CONC OF N IN IRRIGATION WATER = ',&
      F6.1,' ppm'/T10,'MANURE APP FOR P UPTAKE RATE = ',F6.0,' kg/ha/y'/&
      T10,'MANURE APP FOR N UPTAKE RATE = ',F6.0,' kg/ha/y'/T10,&
      'MNUL = ',I3)
 2460 FORMAT(/1X,'-----CURVE NUMBER DISTRIBUTION'/T10,'>95=',F5.2,3X,'>&
      90=',F5.2,3X,'>85=',F5.2,3X,'>80=',F5.2,3X,'>75=',F5.2,3X,'>70=',&
      F5.2,3X,'>65=',F5.2,3X,'>60=',F5.2,3X,'>55=',F5.2,3X,'<55=',F5.2)
 2470 FORMAT(3I3,4I5,7F8.0)
 2480 FORMAT('+',T114,A4)
 2510 FORMAT(2F8.0)
 2520 FORMAT(T10,'PERIOD OF RECORD FOR P5MX =',F4.0,' Y')
 2540 FORMAT(I5,E13.5,2F10.3,2E13.5)
 2550 FORMAT(1X,A4,13F9.1,2X,A4)
 2560 FORMAT(T10,'FURROW DIKE SYSTEM SAFETY FACTOR = ',F5.2)
 2570 FORMAT(/1X,'-----MISCELLANEOUS PARAMETERS'//T10,'PARM',9X,'SCRP1I'&
      ,4X,'SCRP2I',4X,'SCRP1C',7X,'SCRP2C')
 2590 FORMAT(/T10,'IRRIGATION WATER APPLIED'/6X,'Y M D',2X,'VOL(mm)',&
      3X,'HUSC')
 2950 FORMAT(//1X,'-----WIND EROSION DATA')
 2960 FORMAT(T10,'FIELD LENGTH = ',F6.2,' km'/T10,'FIELD WIDTH = ',F6.2,&
      ' km'/T10,'FIELD ANGLE = ',F4.0,' deg'/T10,'WIND SPEED MOD EXP P&
      OWER PARM = ',F5.2/T10,'ACCELERATED WIND EROSION FACTOR = ',F7.3)
 3000 FORMAT(T10,'SOIL ALBEDO = ',F5.2/T10,'MAX NUMBER SOIL LAYERS = ',&
      I3/T10,'MIN THICKNESS FOR LAYER SPLITTING = ',F5.2,' m'/T10,&
      'MIN PROFILE THICKNESS--STOPS SIMULATION = ',F5.2,' m'/T10,&
      'SPLITTING PRIORITY THICKNESS = ',F5.2,' m'/T10,'BIOMASS/ORG C = '&
      ,F6.3/T10,'PASSIVE HUMUS/TOTAL HUMUS = ',F6.3/T10,'WEATHERING CODE&
       = ',F4.0/T10,'ORG C IN TOP .2 M = ',F6.3,' %'/T10,'CULT HISTORY =&
       ',F5.0,' Y')
 3001 FORMAT(T10,'WATER TABLE'/T15,'ANTECEDENT PERIOD FOR RF & PET &
      ACCOUNTING = ',I4,' d'/T15,'DEPTH(m)'/T15,'MIN= ',F6.2/T15,'MAX = ',&
      F6.2,/T15,'INITIAL = ',F6.2)
 3010 FORMAT(1X,A4,12F9.2,11X,A4)
 3030 FORMAT(//I5,6(2X,A4,E12.4)/(5X,6(2X,A4,E12.4)))
 3050 FORMAT(/1X,'-----AVE MO VALUES')
 3060 FORMAT(8X,A8,8F8.2/(10F8.3))
 3070 FORMAT(/1X,'-----GENERATOR SEEDS AFTER',I3,' CYCLES'/(5X,13I12))
 3090 FORMAT(20A4)          
 3100 FORMAT(20I4)            
 3110 FORMAT(T10,'SIMULATION DURATION = ',I4,' Y'/T10,'BEGINNING DATE =&
      ',I4,2I2)
 3111 FORMAT(20F8.0)
 3112 FORMAT(10F12.0)       
 3113  FORMAT(20F8.2)
 3114 FORMAT(10F8.2)
 3120 FORMAT(10F8.0)
 3121 FORMAT(4X,I8,1X,I8,5F9.2,F9.4,F9.3,F9.4,F9.1,4F9.3,2F9.0,2F9.1,&
      F9.2,2F9.3,F9.2,F9.1,1X,A20,1X,A20,1X,I4,F9.3,1X,A80)
 3122 FORMAT(4X,I8,1X,I8,F9.2,F9.0,2F9.2,F9.4,2F9.2,F9.3,F9.4,F9.3,F9.1,&
      F9.2,F9.4)
 3123 FORMAT(4X,I8,1X,I8,2F9.1,F9.0,F9.1,F9.0,F9.2,3F9.0,F9.3,F9.0)
 3140 FORMAT(T15,'RUNOFF RATIO = ',F6.3)
 3191 FORMAT(12F6.0)      
 3200 FORMAT(1X,A4,13F9.0,2X,A4)
 3280 FORMAT(35X,'SIM DUR=',I4,' YR'/T11,'JAN',6X,'FEB',6X,'MAR',6X,&
      'APR',6X,'MAY',6X,'JUN',6X,'JUL',6X,'AUG',6X,'SEP',6X,'OCT',6X,&
      'NOV',6X,'DEC',6X,' YR')
 3281 FORMAT(35X,'SIM DUR=',I4,' YR'/T11,'JAN',9X,'FEB',9X,'MAR',9X,&
      'APR',9X,'MAY',9X,'JUN',9X,'JUL',9X,'AUG',9X,'SEP',9X,'OCT',9X,&
      'NOV',9X,'DEC',9X,' YR')
 3282 FORMAT(9X,'SUBAREA'/8X,'#',6X,'ID',3X,'YR',1X,'YR#',T35,'JAN',10X,&
      'FEB',10X,'MAR',10X,'APR',10X,'MAY',10X,'JUN',10X,'JUL',10X,'AUG',&
      10X,'SEP',10X,'OCT',10X,'NOV',10X,'DEC',10X,'YR')
 3283 FORMAT(T34,'RFV',9X,'Q',7X,'SSF',7X,'YSD',17X,12(6X,A4)&
      /6X,'SA#',6X,'ID',3X,'Y',3X,'M',2X,'D',4X,'|------------(mm)------&
      ------|',2X,'(t/ha)',4X,'PSTN',8X,'|---------------------------&
      ----------------------------(g/ha)--------------------------------&
      --------------------------|')
 3284 FORMAT(3X,'YR ',21(4X,A4,2X))
 3285 FORMAT(/6X,'SA#',6X,'ID',4X,'Y   M   D  Y# ON# HD#',3X,'OPER',3X,&
      'CROP',2X,'YLDkg/ha',2X,'YSDkg/ha',2X,'HAYkg/ha',2X,'AGPMt/ha',3X,&
      'STLt/ha',3X,'STDt/ha',3X,'CNLVg/g',3X,'CNDDg/g')      
 3286 FORMAT(T52,'AP RATE',3X,'MN',5X,'NH3',6X,'ON',6X,'MP',6X,'OP'/6X,&
      'SA#',6X,'ID',4X,'Y M D  Y# ON# HD#',4X,'FERT',5X,'|--------------&
      ----(kg/ha)--------------------|')
 3287 FORMAT(/5X,'ORDER',1X,'|--SA--|',3X,'DP10'/8X,'#',4X,'#',2X,'ID',&
      3X,'(kg/ha)',5X,'FRACT',5X,'ACCUM')
 3288 FORMAT(4X,'CMD',5X,'IDO',5X,'ID#',3X,'Y   M   D',3X,'QPm3/s',7X,&
      'TPh',7X,'Qmm',5X,'SMQmm',5X,'SMHmm')
 3300 FORMAT(1X,A4,13F9.2,2X,A4)
 3301 FORMAT(1X,A4,13F9.3,2X,A4)
 3302 FORMAT(1X,A4,13E12.4,2X,A4)
 3320 FORMAT(/1X,'-----PEAK FLOW RATE STATS(mm/h)',T70,'UNIT PEAKS(1/h)&
      '/T10,'MAX = ',F9.2,5X,'MEAN = ',F6.2,5X,'ST DV = ',F6.2,5X,'MAX =&
       ',F10.4,5X,'MEAN = ',F9.4,5X,'NO PKS = ',I6)
 3360 FORMAT(T10,'TILE DRAIN DEPTH = ',F5.2,' m'/T10,'DRAIN TIME = ',F5.&
      2,' d'/T10,'DRAINAGE RATE = ',F6.1,' mm/h')
 3370 FORMAT(/T10,'FERTILIZER APPLIED'/T7,'DATE   ID',T32,'WT',3X,'DPTH'&
      ,4X,'|-----------FRACTION OF WT------------|'/6X,'Y M D',2X,'NO',&
      3X,'NAME',6X,'(kg/ha)',2X,'(m)',5X,'MIN N',3X,'NH3',5X,'ORG N',3X,&
      'MIN P',3X,'ORG P',3X,'HUSC')
 3380 FORMAT(3X,I4,2I2,I4,2X,A8,1X,F8.0,7F8.3)
 3381 FORMAT('1'/T5,'APEX1501 v20181201',2X,I4,2I2,2X,2(I2,':'),I2)    ! Added date to the version - Luca Doro
 3390 FORMAT(T10,'COSTS'/T15,'N FERT = ',F5.2,' $/KG'/T15,'P FERT = ',&
      F5.2,' $/KG'/T15,'LIME = ',F5.2,' $/t',/T15,'IRR WATER = ',F5.2,& 
      ' $/mm')
 3440 FORMAT(T2,A4,12F9.1,11X,A4)
 3460 FORMAT(T2,A4,12F9.3,11X,A4)
 3536 FORMAT(T15,'HERD ID=',I3,2X,'NUM=',I15,' HD',2X,'FEED AREA TIME=',&
      F6.2,' H/D',2X,'GRAZE LIMIT=',F5.2,' t/ha')
 3540 FORMAT(//1X,'____________________GENERAL INFORMATION______________&
      _______'/)
 3620 FORMAT(T15,'MAX ANNUAL VOL APPL TO A CROP = ',F6.0,' mm'/T15,'MIN &
      SINGLE APPL VOL = ',F6.0,' mm'/T15,'MAX SINGLE APPL VOL = ',F6.0,&
      ' mm')
 3630 FORMAT(T15,'SOIL-WATER TENSION TRIGGER = ',F5.0,' KPA'/T15,'MIN AP&
      PL INTERVAL = ',I3,' d')
 3650 FORMAT(/T10,'__________RAIN,',2(1X,A4,','),' ARE INPUT__________')
 3660 FORMAT(/T10,'__________RAIN,',3(1X,A4,','),' ARE INPUT__________')
 3670 FORMAT(/T10,'__________RAIN,',4(1X,A4,','),' ARE INPUT__________')     
 3680 FORMAT(T15,'PLANT WATER STRESS TRIGGER = ',F5.2/T15,'MIN APPL INT&
      ERVAL = ',I3,' d')
 3730 FORMAT(T10,'LATITUDE = ',F7.2,' deg'/T10,'LONGITUDE = ',F7.2,&
      ' deg'/T10,'ELEVATION = ',F7.1,' m')
 3750 FORMAT(T10,'CO2 CONC IN ATMOSPHERE = ',F6.0,' ppm')
 3770 FORMAT(6X,9(1X,A4,1X,E12.4))
 3850 FORMAT(/T10,'PESTICIDES APPLIED'/T26,'RATE',4X,'KILL',4X,'COST'/&
      6X,'Y M D',2X,'NAME',6X,A7,3X,'EFF',4X,'($/KG)',3X,'HUSC')
 3860 FORMAT(5X,3I2,2X,A8,4F8.2)
 3880 FORMAT(/10X,A16,3F10.0,F10.2,F10.0,F10.2)                          
 3900 FORMAT(T43,'HALF LIFE(DAYS) WASH OFF',T82,'COST'/T13,'NAME',10X,&
      'SOLUBILITY',6X,'SOIL',3X,'FOLIAGE',2X,'FRACTION',6X,'KOC',5X,'($/&
      kg)')
 3910 FORMAT(/11X,10(2X,A8,2X))
 3920 FORMAT(5X,A4,10F10.0)
 3950 FORMAT(5X,A4,10F10.4)
 3970 FORMAT(T10,'HYDROLOGIC SOIL GROUP = ',A1/T10,'LAND USE NUMBER = ',I3&
      /T10,'RUNOFF CN2 = ',F4.1/T10,'SLP ADJ CN2 = ',F4.1/T10,&
      'PARM(92)_ADJ = ',F4.1/T10,'CN SCRV SCRP(35)= ',2F6.0/T10,&
      'CN CNSC(4)= ',2E13.5)
 3990 FORMAT(11X,'SA'/8X,'#',6X,'ID',4X,'Y   M   D ',13X,'GWST',8X,&
      'GWEL',100(8X,A4))      
 4080 FORMAT(T15,'PAW DEFICIT TRIGGER = ',F5.0,' mm'/T15,'MIN APPL INTE'&
      ,'RVAL = ',I3,' d')
 4081  FORMAT(6X,'SA#',6X,'ID',3X,'YR YR#',60(8X,A4))
 4082 FORMAT(11X,'SA'/8X,'#',6X,'ID',4X,'Y   M   D',6X,'CPNM',100(6X,A4))
 4083 FORMAT(4X,'Y   M   D',100(8X,A4))
 4084 FORMAT(7X,'SA',T31,'RFRA',7X,'QSA',8X,'QI',8X,'EV',8X,'SP',8X,&
      'QO',7X,'RSV',6X,'RSVP',6X,'RSVE',8X,'YI',8X,'YO',7X,'DEP',6X,'RSSA',&
      T220,'SW BY LAYER'/4X,'#',2X,'ID',4X,'Y   M   D',3X,'|----------------&
      --------------------------(M3)------------------------------------------|',&
      1X,'|-----------(t/ha)----------|',4X,'(ha)')
 4085 FORMAT(9X,'SA#',7X,'ID',5X,'WSAha',7X,'CN2',7X,'YW0',2X,'OCPDt/ha',&
      6X,'FSFN',6X,'FSFP',3X,'PRBmm/h',2X,'PRAVmm/h',5X,'TCMNh',3X,&
      'CYAVppm',3X,'CYMXppm',40(5X,A4,1X))
 4087 FORMAT(8X,A2,I8,2X,I8,1X,I8,F12.2,I4,F8.2,F8.3,F8.2,3F8.1,I8,2F8.3&
      ,I8,2F8.3)
 4088 FORMAT(5X,'ENMN=',E13.6,2X,'ENMA=',E13.6,2X,'EON =',E13.6,2X,&
      'ESDN=',E13.6,2X,'ESON=',E13.6,2X,'EUNM=',E13.6/5X,'ENMU=',&
      E13.6,2X,'ENOU=',E13.6)
 4089 FORMAT(5X,'EPML=',E13.6,2X,'EPMA=',E13.6,2X,'EPMS=',E13.6,2X,&
      'EPO =',E13.6,2X,'EFOP=',E13.6,2X,'ESDP=',E13.6/5X,'ESOP=',E13.6,&
      2X,'EUPM=',E13.6,2X,'EPMU=',E13.6,2X,'EPOU=',E13.6)
 4091 FORMAT(6X,'SA#',6X,'ID',3X,'YR',2X,'YR#',1X,'CPNM',4X,'YLDG',6X,&
      'YLDF',6X,'BIOM',8X,'WS',8X,'NS',8X,'PS',8X,'KS',8X,'TS',8X,'AS',&
      8X,'SS',6X,'ZNO3',7X,'ZQP',6X,'AP15',7X,'ZOC',6X,'OCPD',6X,'RSDP',&
      6X,'ARSD',6X,'IRGA',8X,'FN',8X,'FP',6X,'FNMN',6X,'FNMA',7X,'FNO',&
      7X,'FPL',7X,'FPO',6X,'YTHS',6X,'YWTH')
 4094 FORMAT(5X,'ISA = ',I4,2X,'BEGINNING RZSW = ',F10.1,' M3')
 4993 FORMAT(12X,6F11.0)
 5005 FORMAT('AVEAN',I6,I9,I6,50E12.4)	    
 5007 FORMAT(/T10,'LAND USE SUMMARY'/T25,'AREA'/10X,'CROP',5X,'(ha)',4X,&
      'FRACTION',2X,'YLDt/ha')
 5008 FORMAT(10X,A4,F10.1,2F10.3)
 5012 FORMAT(/T10,'WIND DIRECTION SUMMARY'/T15,'N= ',F6.3,2X,'NNE= ',&
      F6.3,2X,'NE= ',F6.3,2X,'ENE= ',F6.3,2X,'E= ',F6.3,2X,'ESE= ',F6.3,&
      2X,'SE= ',F6.3,2X,'SSE= ',F6.3,2X,'S= ',F6.3,2X,'SSW= ',F6.3,2X,&
      'SW= ',F6.3,2X,'WSW= ',F6.3,2X,'W= ',F6.3,2X,'WNW= ',F6.3,2X,&
      'NW= ',F6.3,2X,'NNW= ',F6.3)
 5016 FORMAT(/T10,'DUST DISTRIBUTION SUMMARY'/T15,'SA(#  ID)',2X,&
      'PM10(kg/ha)   FRACTION     ACCUM')
 5017 FORMAT(15X,2I8,3E13.5)
 5018 FORMAT(11X,'TOTAL DP10=',E13.5,' KG',2X,'EM10=',E13.5,' KG')
 5031 FORMAT(10X,A80)
 5034 FORMAT(1X,'RUN NAME  ',' SITE FILE  ','SOIL FILE   ',&
      'MANAGEMENT  ','RUN#IRO# WS#',' YR  RT#',38(4X,A4),4X,'OCPD',4X,&
      ' TOC',4X,'ITOC',4X,'APBC',4X,' TAP',4X,'TNO3',2X,'CROP',4X,&
      ' STD',4X,' STL',4X,' LAI',4X,'YLDG',4X,'YLDF',4X,'BIOM',4X,&
      'YLDN',4X,'YLDP',2X,'CROP',4X,' STD',4X,' STL',4X,' LAI',&
      4X,'YLDG',4X,'YLDF',4X,'BIOM',4X,'YLDN',4X,'YLDP')
      END
      
      
      
      
      
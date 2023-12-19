      SUBROUTINE EYSED(JRT)
!     APEX1501
!     THIS SUBPROGRAM PREDICTS DAILY SOIL LOSS CAUSED BY WATER EROSION
!     AND ESTIMATES THE NUTRIENT ENRICHMENT RATIO.
      USE PARM
      
      !rtb salt
      integer m,dum
      real    sub_area_m2,water_volume,salt_mass_kg,&
              tds_conc_soil,ion_conc_soil(8),&
              sed_conc_runoff,tds_conc_runoff,&
              ion_fract,ion_conc_runoff,ion_mass_runoff, &
              ca_meq,mg_meq,na_meq,x1,x2,x3,x4,x5,x6,x7
      
      JRT=0
      LD1=LID(1,ISA)
      WSAX=WSA(ISA)
      CALL EYCC(ISA)
      IF(ISLF>0)THEN
          XX=CVF(ISA)*USL(ISA)
      ELSE
          XX=CVF(ISA)*RSLK(ISA)
      END IF
      QMM=.1*QVOL(IDO)/WSAX
      YLM=QMM
      IF(RFV(IRF(ISA))>12.7)THEN
          ! USLE
          YSD(2,IDO)=MIN(YLM,EI*XX*1.292)*WSAX
          SMM(25,MO,ISA)=SMM(25,MO,ISA)+CVF(ISA)*EI
          SMM(24,MO,ISA)=SMM(24,MO,ISA)+EI
          VAR(24,ISA)=EI
          X1=EI*CVF(ISA)
          RXM=BETA/(1.+BETA)
	      RLFX=UPSX(ISA)**RXM
	      YI=MIN(YLM,.5*X1*RSK(ISA))
          YSSK=1.5*X1*RLFX*RSK(ISA)
	      SUM=PSZX(1)*SAN(LD1,ISA)
	      SUM=SUM+PSZX(2)*SIL(LD1,ISA)
	      SUM=SUM+PSZX(3)*CLA(LD1,ISA)
	      SUM=.01*PRMT(65)*SUM/(QPR(IDO)+1.E-5)
	      T2=PRMT(66)*QPR(IDO)*STP(ISA)
          ! RULSE2
	      IF(T2>YI)THEN
	          YSD(1,IDO)=YSSK
	      ELSE
	          YSD(1,IDO)=MAX(0.,YI+SUM*(T2-YI))
	      END IF
	      YSD(1,IDO)=MIN(YLM,YSD(1,IDO))*WSAX
          IF(QMM>0.)CALL ERHEM
      END IF
      IF(QMM<1.)THEN
          JRT=1
          RETURN
      END IF
      REPC=REP*(QMM/(RFV(IRF(ISA))+.1))**.1
      QQ=QMM*RQRB(IDO)
      X3=1.+WSAX
      ! MUSS
      YSD(3,IDO)=MIN(YLM,.79*X3**.009*QQ**.65*XX)*WSAX
      X1=SQRT(QQ)
      ! MUST
      YSD(5,IDO)=MIN(YLM,PRMT(33)*X1*XX)*WSAX
      CX(MO,ISA)=CX(MO,ISA)+1.
      ! MUSL
      YSD(4,IDO)=MIN(YLM,1.586*X3**.12*QQ**.56*XX)*WSAX
      LD1=LID(1,ISA)      
      IF(RSDM(LD1,ISA)>0.)THEN
          X2=PRMT(71)*AGPM(ISA)
          IF(X2<10.)THEN
              X2=EXP(-X2)
          ELSE
              X2=.00005
          END IF
          YMNU(IDO)=MIN(.9*RSDM(LD1,ISA),PRMT(62)*X1*X2*SLF(ISA)*&
          PEC(ISA)*RSDM(LD1,ISA)**PRMT(68))
          YMNU(IDO)=MAX(0.,YMNU(IDO))*WSAX
          IF(IDFH(ISA)>0)SYMU=SYMU+YMNU(IDO)
      END IF
      IF(URBF(ISA)>0.)THEN
          DO I=1,6
              YSD(I,IDO)=YSD(I,IDO)*(1.-URBF(ISA))+URBF(ISA)*QURB(IDO)*.00003
          END DO
      END IF
      IF(YSD(NDRV,IDO)>PRMT(93))SMM(144,MO,ISA)=SMM(144,MO,ISA)+1.
      DRTO=MAX(.001,MIN(.99,(RQRB(IDO)/(REPC+1.E-5))**.56))
      DRAV(IDO)=DRAV(IDO)+DRTO
      B1=LOG(DRTO)/4.47
      SUM=0.
      DO I=1,NSZ
          PCTH(I,IDO)=PCT(I,ISA)
          X1=MAX(-10.,B1*PSZ(I))
          PCT(I,IDO)=MAX(.0001,PCT(I,ISA)*EXP(X1))
          SUM=SUM+PCT(I,IDO)
      END DO
      PSZM(IDO)=0.
      DO I=1,NSZ
          PCT(I,IDO)=PCT(I,IDO)/(SUM+1.E-10)
          PSZM(IDO)=PSZM(IDO)+PCT(I,IDO)*PSZY(I)
      END DO
      CY=YSD(NDRV,IDO)/QVOL(IDO) !sediment concentration (mg/L)
      IF(IERT==0)THEN
          B2=LOG10(DRTO)/2.699
          B1=1./.1**B2
          ERTO=B1*(CY+1.E-4)**B2
          ERTP=ERTO
      ELSE
          ERTO=PRMT(54)/(CY+1.E-4)**PRMT(55)
          ERTP=PRMT(57)/(CY+1.E-4)**PRMT(58)
      END IF
      
      
      !rtb salt
      if(ISALT>0) then
      !calculate TDS in the top soil layer
      sub_area_m2 = WSA(ISA) * 10000. !area of the subarea in m2
      water_volume = (SWST(LD1,ISA)/1000.) * sub_area_m2 !m * m2 = m3
      tds_conc_soil = 0.
      do m=1,mion
        salt_mass_kg = wsalt(LD1,ISA,m) * WSA(ISA) !kg/ha * ha = kg of salt    
        ion_conc_soil(m) = (salt_mass_kg * 1000.) / water_volume !g/m3 = mg/L
        tds_conc_soil = tds_conc_soil + ion_conc_soil(m)
      enddo
      
      !calculate TDS concentration (mg/L) in the runoff water (using polynomial regression TDS model)
      ca_meq = ion_conc_soil(2) * 0.049903
      mg_meq = ion_conc_soil(3) * 0.082288
      na_meq = ion_conc_soil(4) * 0.043498
      x1 = QMM !mm of runoff
      x2 = CY / 1000. !g/m3 --> g/L
      x3 = ROK(LD1,ISA) !rock fraction
      if((mg_meq+ca_meq).gt.0) then
        x4 = na_meq / (sqrt(0.5*(mg_meq + ca_meq))) !sar
      else
        x4 = 0.
			endif
      x5 = ion_conc_soil(2)+ion_conc_soil(3)+ion_conc_soil(4)+ion_conc_soil(5) !cec in g/m3 = ca + mg + k + na
      x6 = PH(LD1,ISA) !pH
      x7 = ECND(LD1,ISA) !electrical conductivity
      tds_conc_runoff = 160256 - &
											 (1.246 * x1) + &
                       (0.0901 * (x2**2)) - &
                       (0.00031 * (x2**3)) + &
                       (1.296 * x3) - &
                       (17.44 * x4) + &
                       (0.1557 * (x4**2)) + &
                       (1.287 * x5) - &
                       (60865 * x6) + &
                       (7659 * (x6**2)) - &
                       (319 * (x6**3)) + &
                       (171 * x7) - &
                       (14.03 * (x7**2))
      
      if(tds_conc_runoff.lt.0) then
        tds_conc_runoff = 0.
      endif
                       
      !Weltz (2020) linear relationship
      !sed_conc_runoff = CY / 1000. / 1000. !convert g/m3 --> kg/L
      !tds_conc_runoff = (5548 * sed_conc_runoff) + 70.7 !g/m3
              
      !estimate concentration of each salt ion in the runoff (based on fractions in the soil layer)
      do m=1,mion
        if(tds_conc_soil.gt.0) then
          ion_fract = ion_conc_soil(m) / tds_conc_soil
				else
          ion_fract = 0.
        endif
        ion_conc_runoff = tds_conc_runoff * ion_fract !g/m3
        ion_mass_runoff = ion_conc_runoff * QVOL(IDO) !g/m3 * m3 = g
        ion_mass_runoff = (ion_mass_runoff / 1000.) / WSA(ISA) !g --> kg/ha
        !check to see if mass is available in the first layer (kg/ha)
        if(wsalt(LD1,ISA,m).lt.0) then
          dum = 10
        endif
        if(ion_mass_runoff.gt.wsalt(LD1,ISA,m)) then
          ion_mass_runoff = wsalt(LD1,ISA,m)
          wsalt(LD1,ISA,m) = 0.
        else
          wsalt(LD1,ISA,m) = wsalt(LD1,ISA,m) - ion_mass_runoff  
        endif
        surfsalt(IDO,m) = ion_mass_runoff !kg/ha
        if(isnan(surfsalt(IDO,m))) then
          dum = 10
        endif
      enddo
      endif
      
      RETURN
      END
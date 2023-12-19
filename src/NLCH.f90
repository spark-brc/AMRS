      SUBROUTINE NLCH(L1)
!     APEX1501
!     THIS SUBPROGRAM ESTIMATES DAILY SOL N LEACHING BY PERCOLATION AND
!     LATERAL SUBSURFACE FLOW FOR ALL LAYERS EXCEPT THE SURFACE LAYER.
      USE PARM
      WNO3(ISL,ISA)=WNO3(ISL,ISA)+VNO3(L1,ISA)
      WSLT(ISL,ISA)=WSLT(ISL,ISA)+VSLT(ISA)
      
      !rtb salt
      if(ISALT>0) then
      do m=1,mion
        wsalt(ISL,ISA,m) = wsalt(ISL,ISA,m) + vsalt(ISA,m) !kg/ha
        if(wsalt(ISL,ISA,m).lt.0) then
          vsalt(ISA,m) = wsalt(ISL,ISA,m) * (-1)
          wsalt(ISL,ISA,m) = 0.
        endif
      enddo
      endif
      
      WPML(ISL,ISA)=WPML(ISL,ISA)+VAP(ISA)
      WPMU(ISL,ISA)=WPMU(ISL,ISA)+VPU(ISA)
      SSFN=0.
      SSFK=0.
      VNO3(ISL,ISA)=0.
      VAP(ISA)=0.
      VPU(ISA)=0.
      VSLT(ISA)=0.
      
      !rtb salt
      if(ISALT>0) then
      do m=1,mion
		    vsalt(ISA,m) = 0.
      enddo
      endif
      
      SSST=0.
      QSFN=0.
      QSFP=0.
      VH=QSF(ISL,ISA)+CPFH(ISL,ISA)
      V=PKRZ(ISL)+SSF(ISL,ISA)+VH
      IF(V>0.)THEN 
          X2=PO(ISL,ISA)*PRMT(4)*(1.-.01*ROK(ISL,ISA))**2   
          VP=V/X2
          IF(VP>5.)THEN
              X1=.99
          ELSE
              X1=1.-EXP(-VP)
          END IF
          IF(WNO3(ISL,ISA)>0.)THEN
              VQN=WNO3(ISL,ISA)*X1
              VQNU=WNMU(ISL,ISA)*X1
              VV=(VQN+VQNU)/(V+1.E-5)
              IF(WNO3(ISL,ISA)<VQN)VQN=WNO3(ISL,ISA)
              WNO3(ISL,ISA)=WNO3(ISL,ISA)-VQN
              IF(WNMU(ISL,ISA)<VQNU)VQNU=WNMU(ISL,ISA)
              WNMU(ISL,ISA)=WNMU(ISL,ISA)-VQNU
              VNO3(ISL,ISA)=VV*PKRZ(ISL)
              SSFN=VV*SSF(ISL,ISA)
              QSFN=VV*VH
          END IF
          IF(SOLK(ISL,ISA)>0.)THEN
              X3=SOLK(ISL,ISA)*X1
              SOLK(ISL,ISA)=SOLK(ISL,ISA)-X3
              CK1=X3/(V+1.E-5)
              VSK(ISA)=CK1*PKRZ(ISL)
              SSFK=CK1*(SSF(ISL,ISA)+VH)
          END IF
          !IF(WSLT(ISL,ISA)>0.)THEN
          !    X3=WSLT(ISL,ISA)*X1
          !    WSLT(ISL,ISA)=MAX(1.E-5,WSLT(ISL,ISA)-X3)
          !    CS1=X3/(V+1.E-5)
          !    VSLT(ISA)=CS1*PKRZ(ISL)
          !    SSST=CS1*(SSF(ISL,ISA)+VH)
          !END IF 
          
          
					!rtb salt
          if(ISALT>0) then
          do m=1,mion
            if(wsalt(ISL,ISA,m).gt.0) then
              X3 = wsalt(ISL,ISA,m) * X1
              wsalt(ISL,ISA,m) = MAX(1.E-5,wsalt(ISL,ISA,m)-X3)
              CS1 = X3 / (V + 1.E-5)
              
              !salt in percolating water (kg/ha)
              vsalt(ISA,m) = CS1 * PKRZ(ISL)
              
              !salt in lateral subsurface flow (kg/ha)
              latsalt(IDO,m) = CS1*(SSF(ISL,ISA)+VH)  
              
              !salt mass in quick return flow (kg/ha)
              CS = soil_csalt(ISL,ISA,m) !g/m3
              sub_area_m2 = WSA(ISA) * 10000. !m2
              VH_volume = sub_area_m2 * (VH/1000.) !m3
              qrfsalt(IDO,m) = ((CS*VH_volume)/1000.) / WSA(ISA) !kg/ha
              !qrfsalt(IDO,m) = VV * VH !VH = quick return flow (mm)
              !VQS = wsalt(ISL,ISA,m) * X1
              !VV = VQS / (V + 1.E-5)
              
            endif
          enddo
          endif
          
          IF(WPML(ISL,ISA)>0..OR.WPMU(ISL,ISA)>0.)THEN
              X1=CLA(ISL,ISA)
              IF(SCLM(29)>0.)X1=MIN(X1,SCLM(29))
              F=X1/(X1+EXP(SCRP(29,1)-SCRP(29,2)*X1))
              DK=F*PRMT(96)
              XX=MIN(.75,V/(WT(ISL,ISA)*DK))
              IF(LBP==0)THEN
	              ! GLEAMS LINEAR EQ
	              X3=XX*(WPML(ISL,ISA)+WPMU(ISL,ISA)) !manure portion is added, jaehak 2020
              ELSE
                  ! LANGMUIR EQ SOLUTION
                  X2=(WPML(ISL,ISA)+WPMU(ISL,ISA)) !manure portion is added, jaehak 2020
                  QQ=MIN(.1*V,5.)
                  AD1=0.
                  AD2=0.
                  QT=0.
                  IND=0
                  !WRITE(KW(1),'(T10,A,3I4,A,F8.3)')'Y-M-D',IY,MO,KDA,' QVOL=',V
                  DO 
                      AD1=AD1+QQ
                      IF(AD1>V)THEN
                          AD1=AD1-QQ
                          QQ=V-AD1
                          AD1=V
                          IND=1
                      END IF    
                      CS=1000.*X2/WT(ISL,ISA)
                      X1=MAX(1.,CPMX(ISA)-CS)
	                  CL=10.*CS/(PRMT(96)*X1)    
	                  QT=.01*CL*QQ
	                  AD2=AD2+QT
	                  !WRITE(KW(1),2)WPML(ISL,ISA),CS,CL,QT,AD1,AD2
	                  IF(IND>0)EXIT
	                  X2=X2-QT
                  END DO    
                  X3=AD2
              END IF    
              IF(X3>WPML(ISL,ISA)+WPMU(ISL,ISA)) X3=WPML(ISL,ISA)+WPMU(ISL,ISA)
              VAP(ISA)=X3*PKRZ(ISL)/V
              QSFP=X3*QSF(ISL,ISA)/V
              WPML(ISL,ISA)=WPML(ISL,ISA)-QSFP*WPML(ISL,ISA)/(WPML(ISL,ISA)+WPMU(ISL,ISA))
              WPMU(ISL,ISA)=WPMU(ISL,ISA)-QSFP*WPMU(ISL,ISA)/(WPML(ISL,ISA)+WPMU(ISL,ISA))
          END IF
      END IF
      RETURN
      END
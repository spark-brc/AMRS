      SUBROUTINE NYNIT
!     APEX1501
!     THIS SUBPROGRAM ESTIMATES DAILY SOL N LEACHING BY PERCOLATION AND
!     LATERAL SUBSURFACE FLOW AND SOL N LOSS IN RUNOFF FOR THE SURFACE
!     LAYER
      USE PARM
      LD1=LID(1,ISA)
      WSAX=WSA(ISA)
      QMM=.1*QVOL(IDO)/WSAX
      WNO3(LD1,ISA)=WNO3(LD1,ISA)+VNO3(LD1,ISA)
      WSLT(LD1,ISA)=WSLT(LD1,ISA)+VSLT(ISA)
      
      !rtb salt
      if(ISALT>0) then
      do m=1,mion
        wsalt(LD1,ISA,m) = wsalt(LD1,ISA,m) + vsalt(ISA,m) !kg/ha
      enddo
      endif
      
      VH=QSF(LD1,ISA)+CPFH(LD1,ISA)
      V=PKRZ(LD1)+QMM+SSF(LD1,ISA)+VH
      VNO3(LD1,ISA)=0.
      SSFN=0.
      SSFK=0.
      QSFN=0.
      VSLT(ISA)=0.
      
      !rtb salt
      if(ISALT>0) then
      do m=1,mion
        vsalt(ISA,m) = 0.
      enddo
      endif
      
      SSST=0.
      VSK(ISA)=0.
      IF(V>0.)THEN
          X2=V/(PRMT(4)*PO(LD1,ISA))
          IF(X2>10.)THEN
		      X3=.99
	      ELSE
		      X3=1.-EXP(-X2)
	      END IF
	      X4=PKRZ(LD1)+PRMT(14)*(QMM+SSF(LD1,ISA)+VH)
          IF(WNO3(LD1,ISA)>0.)THEN
              VQN=WNO3(LD1,ISA)*X3
              VQNU=WNMU(LD1,ISA)*X3
              CO=(VQN+VQNU)/X4
			  CS=PRMT(14)*CO
			  QN(IDO)=QN(IDO)+WSAX*(QMM*CS*(1.-URBF(ISA)))+URBF(ISA)*QURB&
			  (IDO)*.0005
			  VNO3(LD1,ISA)=PKRZ(LD1)*CO
			  SSFN=CS*SSF(LD1,ISA)
			  QSFN=CS*VH
			  WNO3(LD1,ISA)=WNO3(LD1,ISA)-VQN
			  WNMU(LD1,ISA)=MAX(0.,WNMU(LD1,ISA)-VQNU)
          END IF
		  IF(SOLK(LD1,ISA)>0.)THEN
              X1=SOLK(LD1,ISA)*X3
              SOLK(LD1,ISA)=SOLK(LD1,ISA)-X1
              COK=X1/X4
              CSK=PRMT(14)*COK
              VSK(ISA)=COK*PKRZ(LD1)
              SSFK=CSK*SSF(LD1,ISA)
              QSK=CSK*QMM
              SMM(150,MO,ISA)=SMM(150,MO,ISA)+QSK
              VAR(150,ISA)=QSK
          END IF
          !IF(WSLT(LD1,ISA)>0.)THEN
			    !  X5=(WSLT(LD1,ISA)+1.E-5)*X3
			    !  COSL=X5/X4
			    !  CSSL=PRMT(14)*COSL
			    !  SSST=CSSL*SSF(LD1,ISA)
			    !  VSLT(ISA)=COSL*PKRZ(LD1)
			    !  QSLT=CSSL*QMM
			    !  SMM(130,MO,ISA)=SMM(130,MO,ISA)+QSLT
			    !  VAR(130,ISA)=QSLT
			    !  WSLT(LD1,ISA)=WSLT(LD1,ISA)-X5
		      !END IF
          
          
          !rtb salt
          if(ISALT>0) then
          do m=1,mion
            if(wsalt(LD1,ISA,m).gt.0) then
              X5 = (wsalt(LD1,ISA,m) + 1.E-5) * X3
              COSL = X5 / X4
              CSSL = PRMT(14) * COSL
              
              !salt mass in lateral subsurface flow (kg/ha) 
              latsalt(IDO,m) = (CSSL*SSF(LD1,ISA)) 
              
              !salt mass in percolating water (kg/ha)
              vsalt(ISA,m) = COSL * PKRZ(LD1) 

              !salt mass in surface runoff (kg/ha)
              CS = soil_csalt(LD1,ISA,m) !g/m3
              sub_area_m2 = WSA(ISA) * 10000. !m2
              QM_volume = sub_area_m2 * (QMM/1000.) !m3
              surfqsalt(IDO,m) = ((CS*QM_volume)/1000.) / WSA(ISA) !kg/ha
              !surfqsalt(IDO,m) = CSSL * QMM 
              write(8642,500) IYR,IDA,isa,m,QMM,CS
              
              !salt mass in quick return flow (kg/ha)
              CS = soil_csalt(LD1,ISA,m) !g/m3
              sub_area_m2 = WSA(ISA) * 10000. !m2
              VH_volume = sub_area_m2 * (VH/1000.) !m3
              qrfsalt(IDO,m) = ((CS*VH_volume)/1000.) / WSA(ISA) !kg/ha
              !VQS = wsalt(LD1,ISA,m) * X3
              !CO = VQS / X4
              !CS = PRMT(14) * CO
              
              !update salt mass in soil layer (kg/ha)
              wsalt(LD1,ISA,m) = wsalt(LD1,ISA,m) - X5 
              if(wsalt(LD1,ISA,m).lt.0) then
                vsalt(ISA,m) = wsalt(LD1,ISA,m) * (-1)
                wsalt(LD1,ISA,m) = 0.
              endif
              
            endif
          enddo
          endif
          
      END IF	
      
      
500   format(i8,i8,i8,i8,f12.4,f12.4)

      RETURN
      END
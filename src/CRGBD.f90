      SUBROUTINE CRGBD(RGS)
!     APEX1501
!     THIS SUBPROGRAM CALCULATES ROOT GROWTH STRESSES CAUSED BY
!     TEMPERATURE,  ALUMINUM TOXICITY, AND SOIL STRENGTH AND DETERMINES
!     THE ACTIVE CONSTRAINT ON ROOT GROWTH(THE MINIMUM STRESS FACTOR).
      USE PARM
      RGS=1.
      IF(PRMT(2)>1.99)RETURN
      II=3
      XX=2.*STMP(ISL,ISA)/(TOPC(JJK)+TBSC(JJK))
      IF(XX<1.E-10)THEN
          RGS=0.
      ELSE  
          IF(XX<1.)RGS=XX**PRMT(41)
          A0=10.+(ALT(JJK)-1.)*20.
          IF(ALS(ISL,ISA)>A0)THEN
              F=(100.-ALS(ISL,ISA))/(100.-A0)
              CALL CFRG(2,II,F,RGS,.1,JRT)
              IF(JRT==0)THEN
                  CALL SBDSC(BDP(ISL,ISA),PRMT(2),F,ISL,3)
                  XX=ROK(ISL,ISA)
                  IF(SCLM(1)>0.)XX=MIN(XX,SCLM(1))
                  F=F*(1.-XX/(XX+EXP(SCRP(1,1)-SCRP(1,2)*XX)))
                  CALL CFRG(1,II,F,RGS,.1,JRT)
              END IF
          END IF    
      END IF
      STDA(II,JJK,ISA)=STDA(II,JJK,ISA)+(1.-RGS)/NBSL(ISA)
      RETURN
      END
      SUBROUTINE NKUP
!     APEX1501
!     THIS SUBPROGRAM CALCULATES THE DAILY K DEMAND FOR OPTIMAL PLANT
!     GROWTH.
      USE PARM
      IF(NUPC/=0)THEN
          X1=HUI(JJK,ISA)
          IF(SCLM(35)>0.)X1=MIN(HUI(JJK,ISA),SCLM(35))
          CPT=(BK(4,JJK)-BK(3,JJK))*(1.-X1/(X1+EXP(BK(1,JJK)-BK(2,JJK)*X1)))+BK(3,JJK)
      ELSE
          CPT=BK(2,JJK)+BK(1,JJK)*EXP(-BK(4,JJK)*HUI(JJK,ISA))
      END IF
      UK2=CPT*DM(JJK,ISA)*1000.
      IF(UK2<UK1(JJK,ISA))UK2=UK1(JJK,ISA)
      UKM=MIN(4000.*BK(3,JJK)*DDM(JJK),UK2-UK1(JJK,ISA))
      UKM=MAX(0.,UKM)
      RETURN
      END
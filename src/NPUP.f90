      SUBROUTINE NPUP
!     APEX1501
!     THIS SUBPROGRAM CALCULATES THE DAILY P DEMAND FOR OPTIMAL PLANT
!     GROWTH.
      USE PARM
      IF(NUPC/=0)THEN
          X1=HUI(JJK,ISA)
          IF(SCLM(34)>0.)X1=MIN(HUI(JJK,ISA),SCLM(34))
          CPT=(BP(4,JJK)-BP(3,JJK))*(1.-X1/(X1+EXP(BP(1,JJK)-BP(2,JJK)*X1)))+BP(3,JJK)
      ELSE
          CPT=BP(2,JJK)+BP(1,JJK)*EXP(-BP(4,JJK)*HUI(JJK,ISA))
      END IF
      UP2=CPT*DM(JJK,ISA)*1000.
      IF(UP2<UP1(JJK,ISA))UP2=UP1(JJK,ISA)
      UPM=MIN(4000.*BP(3,JJK)*DDM(JJK),UP2-UP1(JJK,ISA))
      UPM=MAX(0.,UPM)
      RETURN
      END
      SUBROUTINE NAJN(UU,AN,DMD,SUX,AJF,JNT)
!     APEX1501
!     THIS SUBPROGRAM COMPUTES ACTUAL N PLANT UPTAKE FROM EACH
!     LAYER(UPTAKE = MINIMUM OF PLANT DEMAND AND SOIL SUPPLY).
      USE PARM
      DIMENSION UU(MSL),AN(MSC,MSA)
      SUM=0.
      RT=DMD/(SUX+1.E-20)
      ! N DEMAND < ESTIMATED SUPPLY
      IF(JNT==0.AND.RT<=1.)THEN
          DO J=1,LRD(ISA)
              ISL=LID(J,ISA)
              ! REDUCE N UPTAKE TO = DEMAND 
              UU(ISL)=UU(ISL)*RT
              SUM=SUM+UU(ISL)
          END DO
      ELSE
          ! N & P DEMAND > ESTIMATED SUPPLY
          RT=AJF*(DMD-SUX)
          RT1=RT
          ! ADJUST UPTAKE TO = DEMAND IF N & P ARE AVAILABLE
          DO J=1,LRD(ISA)
              ISL=LID(J,ISA)
              XX=UU(ISL)+RT
              IF(XX<AN(ISL,ISA))THEN
                  UU(ISL)=XX
                  SUX=SUX+RT1
                  RETURN
              ELSE    
                  IF(AN(ISL,ISA)>0.)THEN
                      RT=RT-AN(ISL,ISA)+UU(ISL)
                      UU(ISL)=AN(ISL,ISA)
                      SUM=SUM+UU(ISL)
                  ELSE
                      UU(ISL)=0.
                  END IF 
              END IF
          END DO
      END IF    
      SUX=SUM
      RETURN
      END
      SUBROUTINE ADAJ(NC,JDT,M,IDA,NYD)
!     APEX1501
!     THIS SUBPROGRAM COMPUTES THE DAY OF THE YEAR, GIVEN THE MONTH AND
!     THE DAY OF THE MONTH.
      DIMENSION NC(13)
      JDT=NC(M)+IDA
      !IF(M>2)JDT=JDT-NYD    ! Original statement commented out by Luca Doro
      if(m>2)jdt=jdt-1   ! New statement added by Luca Doro. We always subtract one because we have 29 days in Feb in vector NC
      RETURN
      END
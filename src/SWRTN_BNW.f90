      SUBROUTINE SWRTN_BNW(CL,SA,OC,BD,WP,FC)
      ! APEX1501
!     THIS SUBPROGRAM USES THE BEHRMAN-NORFLEET-WILLIAMS (BNW) METHOD
!     FOR ESTIMATING SOIL WATER CONTENT AT 33 AND 1500 kpa.
      SI=100.-CL-SA
      WP=.04285+.0001041*SI+.003958*CL+.00001555*CL*SI-.005606*LOG10(OC)
      PO=1.-BD/2.65
      X1=1.72*OC
      AD1=SA+SI+CL+X1
      BDO=(1.6*SA+1.3*SI+1.1*CL+.224*X1)/AD1
      PAO=(.05*SA+.26*SI+.08*CL+.9*X1)/AD1
      RTO=BDO/BD
      PAW=RTO*PAO
      FC=WP+PAW
      RETURN
      END
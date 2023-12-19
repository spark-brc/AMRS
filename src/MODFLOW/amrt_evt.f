      subroutine amrt_evt

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine converts the APEX remaining ET (PET - Actual ET) to rates for MODFLOW's EVT package

        !Import variables
        use PARM, only: MSA !APEX
        use GLOBAL, only: NCOL,NROW,BOTM !MODFLOW
        use amrt_parm !APEX-MODFLOW
        use GWFEVTMODULE, ONLY:EVTR,EXDP,SURF
        implicit none
        
        !Define local variables
        integer i,j

        !MODFLOW variables
        !EVTR: maximum ET rate
        !EXDP: extinction depth (depth at which ET ceases)
        !SURF: surface from which extinction depth is measured

        !map APEX etremain values to MODFLOW grid cells     
        

        !Extinction depth is set in the *.evt MODFLOW input file
        !EXDP

        !store ground surface elevation in SURF
        do i=1,NROW
          do j=1,NCOL
            SURF(j,i) = BOTM(j,i,0)  
          enddo
        enddo

      end subroutine amrt_evt

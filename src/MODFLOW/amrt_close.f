      subroutine amrt_close

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine deallocates arrays for APEX-MODFLOW
      
!     Import variables
      use GWFRIVMODULE, only: RIVAUX, RIVR !MODFLOW
      use amrt_parm !APEX-MODFLOW
      implicit none

!     Deallocate amrt variables
      deallocate(g2s_map)
      deallocate(s2g_map)
      deallocate(grid2riv_id)
      deallocate(grid2riv_len)
      deallocate(riv_nsubs)
      
      deallocate(etremain)
      deallocate(sepbtm)
      deallocate(rchrg)
      deallocate(rchrg_n)
      deallocate(rchrg_p)
      deallocate(percn)
      deallocate(percp)
      
      return
      end
      subroutine amrt_grid2sa2D (mfVar2,amrtVar)

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine converts variables from MODFLOW cells to APEX subareas

!     Import variables
      use PARM      !APEX
      use GLOBAL, only: NCOL,NROW !MODFLOW
      use amrt_parm !APEX-MODFLOW
      implicit none

!     Define local variables
      real mfVar2(NCOL,NROW)
      real amrtVar(MSA)
      integer i, j, row, column
      real perc_area
      
!     Initialize variables
      amrtVar = 0.
      
      !loop through all of the APEX subareas
      do i=1,MSA 
        !loop through all of the cells that contribute to the subarea
        do j=1,size(g2s_map(i)%cell_row)
          row = g2s_map(i)%cell_row(j)
          column = g2s_map(i)%cell_col(j)
          perc_area = g2s_map(i)%cell_perc(j)
          if(column.ne.0 .and. row.ne.0)then
            amrtVar(i) = amrtVar(i) + mfVar2(column,row)*perc_area
          endif
        enddo
      enddo
      
      return
      end
      subroutine amrt_grid2sa3D (mfVar3,amrtVar,location)   !Ali

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine converts 3D MODFLOW cell variables to APEX subareas

!     Import variables
      use PARM       !APEX
      use GLOBAL, only: NCOL,NROW,NLAY !MODFLOW
      use amrt_parm  !APEX-MODFLOW
      implicit none
      
!     Define local variables
      double precision mfVar3(NCOL,NROW,NLAY)
      real amrtVar(MSA)    
      real location(NCOL,NROW)
      integer i, j, row, column, lay
      real perc_area
      
!     Initialize variables
      amrtVar = 0.
      
!     For each subarea, loop through the intersecting MODFLOW cells 
      do i=1,MSA        
        do j=1,size(g2s_map(i)%cell_row)
          row = g2s_map(i)%cell_row(j)
          column = g2s_map(i)%cell_col(j)
          perc_area = g2s_map(i)%cell_perc(j)
          if(column.ne.0 .and. row.ne.0)then
            lay = location(column,row)
            amrtVar(i) = amrtVar(i) + (mfVar3(column,row,lay)*perc_area)
          endif
        enddo
      enddo
      
      return
      end
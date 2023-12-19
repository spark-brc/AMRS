      subroutine amrt_conversion2mf

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine converts the APEX variables to MODFLOW
      
!     Import variables
      use parm !APEX
      use GLOBAL, only:LENUNI,ITMUNI !MODFLOW
      use amrt_parm  !APEX-MODFLOW
      implicit none
      
!     Define local variables
      integer mf_lengthUnit,mf_timeUnit,i,n,j,k,subareaID,
     &        cell_row,cell_col,dum,sub
      real    rchrg1
      mf_lengthUnit = LENUNI + 10
      mf_timeUnit = ITMUNI
      
      !calculate recharge to the water table, based on the percolation from the soil profile
      !(using groundwater delay terms)
      do i=1,MSA
        sub = NBSA(i) !actual subarea
        rchrg1 = rchrg(sub)
        rchrg1 = rchrg1 * 1000. !recharge was converted to m; now change back to mm
        rchrg(sub) = 0.
        rchrg(sub) = ((1.-gw_delaye(sub))*sepbtm(i)) + 
     &               (gw_delaye(sub)*rchrg1)
        if(sepbtm(i).gt.0) then
          dum = 10
        endif
        if (rchrg(sub) < 1.e-6) rchrg(sub) = 0.
      enddo

      amrt_RECH_APEX_tot_mo = amrt_RECH_APEX_tot_mo + rchrg !total recharge for the current month
      amrt_RECH_APEX_tot_yr = amrt_RECH_APEX_tot_yr + rchrg !total recharge for the current year (rtb avg)
      
      !write out recharge values to file
      if(out_APEX_recharge.eq.1 .and.
     &   day_total.eq.outapexmf(apexmf_out_ctr)) then
        write(30001,*)
        write(30001,*) 'Day:',day_total
        write(30001,*) 'Recharge (mm), Soil Perc. (mm) for each subarea'
        do i=1,MSA            
          write(30001,*) NBSA(i),rchrg(NBSA(i)),sepbtm(i)
        enddo
        write(30001,*)
      endif

!     Convert APEX subarea variables to correct MODFLOW units
      call units(rchrg, 14, mf_lengthUnit, 1, MSA, LPYR)   ! to convert units (mm/day to LENUNI/day)
      call units(rchrg, mf_timeUnit, 4, 1, MSA, LPYR)      ! to convert time units (LENUNI/days to LENUNI/ITMUNI)
      call units(etremain, 14, mf_lengthUnit, 1, MSA, LPYR) ! to convert length units (mm/day to LENUNI/day)
      call units(etremain, mf_timeUnit, 4, 1, MSA, LPYR)     ! to convert time units (LENUNI/days to LENUNI/ITMUNI)

      return
      end
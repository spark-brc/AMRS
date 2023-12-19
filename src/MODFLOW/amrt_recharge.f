      subroutine amrt_recharge
!!
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine converts the APEX recharge variable to the MODFLOW recharge array
   
        !Import variables
        use PARM !APEX
        use GLOBAL, only: NCOL,NROW !MODFLOW
        use amrt_parm  !APEX-MODFLOW
        use GWFRCHMODULE, only: RECH
        implicit none
        
        !Define local variables
        integer i,j,n,ctr,sub
        real    volume,area_m2

        !convert to MODFLOW Recharge array
        RECH = 0.
        ctr = 1
        do i=1,nrow
          do j=1,ncol
            volume = 0.
            do n=1,size(s2g_map(ctr)%subarea_id)
              sub = s2g_map(ctr)%subarea_id(n)
              if(sub.gt.0) then
                area_m2 = s2g_map(ctr)%subarea_area(n)
                volume = volume + (rchrg(sub)*area_m2) !m * m2 = m3 of recharge  
              endif
            enddo
            RECH(j,i) = volume !store in MODFLOW array
            ctr = ctr + 1
          enddo
        enddo

        amrt_RECH_tot_mo = amrt_RECH_tot_mo + RECH !total recharge for the current month
        amrt_RECH_tot_yr = amrt_RECH_tot_yr + RECH !total recharge for the current year

        !write out Recharge values
        if(out_MF_recharge.eq.1 .and.
     &    day_total.eq.outapexmf(apexmf_out_ctr)) then
          write(30002,*)
          write(30002,*) 'Day:',day_total
          write(30002,*) 'daily recharge values provided to MODFLOW'
          do i=1,nrow
            write(30002,100) (RECH(j,i),j=1,ncol)
          enddo
        endif

  100 format(1000(e17.10))

      end subroutine amrt_recharge
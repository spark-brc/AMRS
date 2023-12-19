      subroutine amrt_conversion2rt3d

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine converts the necessary APEX variables to RT3D variables
      
      !Import variables
      use parm, only:LPYR,MSA,WSA,ISALT,mion,rchrg_salt,vsalt, !APEX
     &							 wsalt,rchrg_csalt,chan_csalt,NBSA,NISA,IBSA,IEXT,IDOA
      use GLOBAL, only:LENUNI,ITMUNI,NCOL,NROW !MODFLOW
      use GWFRIVMODULE, only:RIVR,MXRIVR !MODFLOW
      use GWFRCHMODULE, only: RECH !MODFLOW
      use mf_rt_link, only: rt_criv,crch !MODFLOW-RT3D linkage
      use amrt_parm !APEX-MODFLOW
      implicit none
      
      !Define local variables
      integer mf_lengthUnit,mf_timeUnit,dum,command
      integer i,j,n,ctr,col,row,lay,mf_gridID,sub
      real rchrgn1,rchrgp1,rchrgsalt1
      real mass_no3,mass_p,mass_so4,mass_ca,mass_mg,mass_na,mass_k,
     &     mass_cl,mass_co3,mass_hco3
      real area_ha,sub_area_m2,area_m2,solute_mass
      real conc_recharge(ncol,nrow),rchrg_solute_conc(MSA)
      real rivlen,rivNitraten,rivPdsolv
      real riv_so4,riv_ca,riv_mg,riv_na,riv_k,riv_cl,riv_co3,riv_hco3

      mf_lengthUnit = LENUNI + 10
      mf_timeUnit = ITMUNI
      

      if(day_total.eq.100) then
        dum = 10
      endif
      
      
      !write out recharge concentration values to file
      write(30010,*) 'NO3 mass (kg/ha) for each subarea for day',
     &                day_total
      write(30014,*) 'P mass (kg/ha) for each subarea for day',
     &                day_total
      if(ISALT>0) then
        write(30018,*) 'Salt mass (kg/ha) for each subarea for day',
     &                day_total
        write(30018,*) 'SO4, Ca, Mg, Na, K, Cl, CO3, HCO3'
      endif

      !compute concentration (g/m3) of solutes in leaching water, for each subarea (for current day)
      do i=1,MSA
        
        sub = NBSA(i) !actual subarea
        area_ha = WSA(i) !subarea area in hectares
        
        !compute volume of recharge water for the subarea
        sub_area_m2 = area_ha * 10000. !m2
        rch_volume(sub) = rchrg(sub) * sub_area_m2 !m depth * m2 area = m3 of recharge
        
        !compute nitrate mass (g) in water reaching the water table
        rchrgn1 = 0.
        rchrgn1 = rchrg_n(sub)
        if (rchrgn1 < 1.e-6) rchrgn1 = 0.0
        rchrg_n(sub) = 0.
        rchrg_n(sub) = ((1.- gw_delaye(sub)) * percn(i)) + 
     &               (gw_delaye(sub) * rchrgn1)
        write(30010,101) sub,area_ha,rchrg_n(sub)
        
        !compute phosphorus mass (g) in water reaching the water table
        rchrgp1 = 0.
        rchrgp1 = rchrg_p(sub)
        if (rchrgp1 < 1.e-6) rchrgp1 = 0.0
        rchrg_p(sub) = 0.
        rchrg_p(sub) = ((1.- gw_delaye(sub)) * percp(i)) + 
     &               (gw_delaye(sub) * rchrgp1)
        write(30014,101) sub,area_ha,rchrg_p(sub)
        
        !compute salt mass (g) in water reaching the water table
        if(ISALT>0) then
        do j=1,mion
          rchrgsalt1 = 0.
          rchrgsalt1 = rchrg_salt(sub,j)
          if (rchrgsalt1 < 1.e-6) rchrgsalt1 = 0.0
          rchrg_salt(sub,j) = 0.
          rchrg_salt(sub,j) = ((1.- gw_delaye(sub)) * vsalt(i,j)) + 
     &                      (gw_delaye(sub) * rchrgsalt1)
        enddo
        write(30018,101) sub,area_ha,(rchrg_salt(sub,j),j=1,mion)
        endif
        
      enddo !go to next subarea
      write(30010,*)
      write(30014,*)
      if(ISALT>0) write(30018,*)


      !convert subarea mass to cell concentration for RT3D arrays
      CRCH = 0.
      ctr = 1
      do i=1,nrow
        do j=1,ncol
          if(day_total.eq.639) then
            if(i.eq.8 .and. j.eq.56) then
              dum = 10
            endif
          endif
        
          if(RECH(j,i).gt.0) then
          mass_no3 = 0.
          mass_p = 0.
          if(ISALT>0) then
            mass_so4 = 0.
            mass_ca = 0.
            mass_mg = 0.
            mass_na = 0.
            mass_k = 0.
            mass_cl = 0.
            mass_co3 = 0.
            mass_hco3 = 0.
          endif
          do n=1,size(s2g_map(ctr)%subarea_id)
            sub = s2g_map(ctr)%subarea_id(n)
            area_ha = IBSA(sub) !total area (ha) of subarea
            area_m2 = s2g_map(ctr)%subarea_area(n) !area (m2) of subarea in the grid cell
            !no3
            solute_mass = rchrg_n(sub) / 10000. * area_m2 * 1000. !kg/ha / 10000 * m2 * 1000 = g
            mass_no3 = mass_no3 + solute_mass !g of no3
            !p
            solute_mass = rchrg_p(sub) / 10000. * area_m2 * 1000.
            mass_p = mass_p + solute_mass !g of p
            if(ISALT>0) then
            !so4
            solute_mass = rchrg_salt(sub,1) / 10000. * area_m2 * 1000.
            mass_so4 = mass_so4 + solute_mass !g of so4
            !ca
            solute_mass = rchrg_salt(sub,2) / 10000. * area_m2 * 1000.
            mass_ca = mass_ca + solute_mass !g of ca
            !mg
            solute_mass = rchrg_salt(sub,3) / 10000. * area_m2 * 1000.
            mass_mg = mass_mg + solute_mass !g of mg
            !na
            solute_mass = rchrg_salt(sub,4) / 10000. * area_m2 * 1000.
            mass_na = mass_na + solute_mass !g of na
            !k
            solute_mass = rchrg_salt(sub,5) / 10000. * area_m2 * 1000.
            mass_k = mass_k + solute_mass !g of k
            !cl
            solute_mass = rchrg_salt(sub,6) / 10000. * area_m2 * 1000.
            mass_cl = mass_cl + solute_mass !g of cl
            !co3
            solute_mass = rchrg_salt(sub,7) / 10000. * area_m2 * 1000.
            mass_co3 = mass_co3 + solute_mass !g of co3
            !hco3
            solute_mass = rchrg_salt(sub,8) / 10000. * area_m2 * 1000.
            mass_hco3 = mass_hco3 + solute_mass !g of hco3
            endif
          enddo
          CRCH(j,i,1) = mass_no3 / RECH(j,i) !g/m3 of no3
          CRCH(j,i,2) = mass_p / RECH(j,i) !g/m3 of p
          if(ISALT>0) then
            CRCH(j,i,3) = mass_so4 / RECH(j,i) !g/m3 of so4
            CRCH(j,i,4) = mass_ca / RECH(j,i) !g/m3 of ca
            CRCH(j,i,5) = mass_mg / RECH(j,i) !g/m3 of mg
            CRCH(j,i,6) = mass_na / RECH(j,i) !g/m3 of na
            CRCH(j,i,7) = mass_k / RECH(j,i) !g/m3 of k
            CRCH(j,i,8) = mass_cl / RECH(j,i) !g/m3 of cl
            CRCH(j,i,9) = mass_co3 / RECH(j,i) !g/m3 of co3
            CRCH(j,i,10) = mass_hco3 / RECH(j,i) !g/m3 of hco3
          endif
          endif
          ctr = ctr + 1
        enddo
      enddo
      
      
      !write out cell-by-cell values to a file (for testing/plotting)
      if(day_total.eq.outapexmf(apexmf_out_ctr)) then
      write(30009,*) 'NO3-N perc. (mg/L) for each cell for current day'
      do i=1,nrow
        write(30009,100) (CRCH(j,i,1),j=1,ncol)
      enddo
      write(30009,*)
      write(30013,*) 'P perc. (mg/L) for each cell for current day'
      do i=1,nrow
        write(30013,100) (CRCH(j,i,2),j=1,ncol)
      enddo
      write(30013,*)
      
      if(ISALT>0) then
        write(30017,*) 'Salt perc. (mg/L) for each cell for current day'
        !SO4
        write(30017,*) 'SO4 conc. (mg/L) for each cell'
        do i=1,nrow
          write(30017,100) (CRCH(j,i,3),j=1,ncol)
        enddo
        !Ca
        write(30017,*) 'Ca conc. (mg/L) for each cell'
        do i=1,nrow
          write(30017,100) (CRCH(j,i,4),j=1,ncol)
        enddo
        !Mg
        write(30017,*) 'Mg conc. (mg/L) for each cell'
        do i=1,nrow
          write(30017,100) (CRCH(j,i,5),j=1,ncol)
        enddo
        !Na
        write(30017,*) 'Na conc. (mg/L) for each cell'
        do i=1,nrow
          write(30017,100) (CRCH(j,i,6),j=1,ncol)
        enddo
        !K
        write(30017,*) 'K conc. (mg/L) for each cell'
        do i=1,nrow
          write(30017,100) (CRCH(j,i,7),j=1,ncol)
        enddo
        !Cl
        write(30017,*) 'Cl conc. (mg/L) for each cell'
        do i=1,nrow
          write(30017,100) (CRCH(j,i,8),j=1,ncol)
        enddo
        !CO3
        write(30017,*) 'CO3 conc. (mg/L) for each cell'
        do i=1,nrow
          write(30017,100) (CRCH(j,i,9),j=1,ncol)
        enddo
        !HCO3
        write(30017,*) 'HCO3 conc. (mg/L) for each cell'
        do i=1,nrow
          write(30017,100) (CRCH(j,i,10),j=1,ncol)
        enddo
        write(30017,*)
      endif
      endif
      

      !loop through the MODFLOW River Cells that are intersected with APEX's subareas (as listed in the apexmf_river2grid.txt file)
      !for each cell, map the NO3 and P river water concentration to the RT3D array
      do i=1,nrivcells_subs
        
        !reset averaged river properties
        rivlen = 0.
        rivNitraten = 0.
        rivPdsolv = 0.

        !loop through the APEX rivers (one for each subarea)
        do j=1,riv_nsubs(i)
          
          !get the subbasin ID
          sub = grid2riv_id(i,j)
          
          !Get the river's properties to be based on a weighted average with: 
          !weights = subbasin's river segment length / total river length in grid cell (rivlen)
          rivlen = rivlen + grid2riv_len(i,j)
          rivNitraten = rivNitraten + 
     &                     (no3_chan(sub) * grid2riv_len(i,j)) !get no3 concentration from APEX reach
          rivPdsolv = rivPdsolv + 
     &                     (p_chan(sub) * grid2riv_len(i,j)) !get P concentration from APEX reach
        enddo
        
        !prevent divide by zero problems
        if(rivlen.eq.0) rivlen = 1. 
        
        !take weighted average of river nitrate concentration
        rivNitraten = rivNitraten / rivlen !mg/L
        rivPdsolv = rivPdsolv / rivlen !mg/L

        !add the concentration values to the correct MODFLOW river cell
        do j=1,MXRIVR
          row = RIVR(2,j)
          col = RIVR(3,j)
          mf_gridID = ((row-1)*NCOL) + col !grid ID  
          if(mf_gridID.eq.grid2riv_cellID(i)) then
            rt_criv(j,1) = rivNitraten !note that nitrate is hard coded to the 1st entry for RT3D
            rt_criv(j,2) = rivPdsolv !note that phosphorus is hard coded to the 2nd entry for RT3D
          endif
        enddo

      enddo !go to the next cell

      !do the same for salt ion concentrations
      if(ISALT>0) then
      do i=1,nrivcells_subs
        
        !reset averaged river properties
        rivlen = 0.
        riv_so4 = 0.
        riv_ca = 0.
        riv_mg = 0.
        riv_na = 0.
        riv_k = 0.
        riv_cl = 0.
        riv_co3 = 0.
        riv_hco3 = 0.

        !loop through the APEX rivers (one for each subarea)
        do j=1,riv_nsubs(i)
          
          !get the subbasin ID
          sub = grid2riv_id(i,j)
          
          !Get the river's properties to be based on a weighted average with: 
          !weights = subbasin's river segment length / total river length in grid cell (rivlen)
          rivlen = rivlen + grid2riv_len(i,j)
          riv_so4 = riv_so4 + 
     &                   (chan_csalt(sub,1) * grid2riv_len(i,j)) !get so4 concentration from APEX reach
          riv_ca = riv_ca + 
     &                   (chan_csalt(sub,2) * grid2riv_len(i,j)) !get ca concentration from APEX reach
          riv_mg = riv_mg + 
     &                   (chan_csalt(sub,3) * grid2riv_len(i,j)) !get mg concentration from APEX reach
          riv_na = riv_na + 
     &                   (chan_csalt(sub,4) * grid2riv_len(i,j)) !get na concentration from APEX reach
          riv_k = riv_k + 
     &                   (chan_csalt(sub,5) * grid2riv_len(i,j)) !get k concentration from APEX reach
			riv_cl = riv_cl + 
     &                   (chan_csalt(sub,6) * grid2riv_len(i,j)) !get cl concentration from APEX reach
          riv_co3 = riv_co3 + 
     &                   (chan_csalt(sub,7) * grid2riv_len(i,j)) !get co3 concentration from APEX reach
          riv_hco3 = riv_hco3 + 
     &                   (chan_csalt(sub,8) * grid2riv_len(i,j)) !get hco3 concentration from APEX reach
        enddo
        
        !prevent divide by zero problems
        if(rivlen.eq.0) rivlen = 1. 
        
        !take weighted average of river salt ion concentration
        riv_so4 = riv_so4 / rivlen !mg/L
        riv_ca = riv_ca / rivlen !mg/L
        riv_mg = riv_mg / rivlen !mg/L
        riv_na = riv_na / rivlen !mg/L
        riv_k = riv_k / rivlen !mg/L
        riv_cl = riv_cl / rivlen !mg/L
        riv_co3 = riv_co3 / rivlen !mg/L
        riv_hco3 = riv_hco3 / rivlen !mg/L

        !add the concentration values to the correct MODFLOW river cell
        do j=1,MXRIVR
          row = RIVR(2,j)
          col = RIVR(3,j)
          mf_gridID = ((row-1)*NCOL) + col !grid ID  
          if(mf_gridID.eq.grid2riv_cellID(i)) then
            rt_criv(j,3) = riv_so4 !note that these salt ions are hard coded to certain entries in the array
            rt_criv(j,4) = riv_ca
            rt_criv(j,5) = riv_mg
            rt_criv(j,6) = riv_na
            rt_criv(j,7) = riv_k
            rt_criv(j,8) = riv_cl
            rt_criv(j,9) = riv_co3
            rt_criv(j,10) = riv_hco3
          endif
        enddo

      enddo !go to the next cell
      endif
      

100   format (1000(f12.4))
101   format (i6,100(f12.4))

      return
      end
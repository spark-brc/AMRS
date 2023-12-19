      subroutine salt_channels !rtb salt
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    calculate and write out the salt ion concentrations in each subarea channel
             
      use parm
      
      implicit none

      integer i,m,I1,I2,II,I3,dum,outlet_sub
      real    vol_runoff,Q_runoff,vol_total,Q_total,salt_mass_kg,area_ha
      real    outlet_csalt(8),outlet_load(8),load(8)
      real    outlet_flow,outlet_area
      

      !loop through the subareas - write out flow rate, loadings, and concentrations for each subarea
      chan_csalt = 0.
      outlet_flow = 0.
      outlet_csalt = 0.
      outlet_load = 0.
      do i=1,MSA !loop through the subareas                                                                         
        I1 = NBSA(IBSA(i))
	  I2 = NISA(I1)
        area_ha = WSA(I2)
	  if(IBSA(i).eq.MSA) then !outlet subarea
          outlet_sub = i
          outlet_area = area_ha
        endif
        if(IEXT(I2)>0) then
	    II = IDOA(I2) !just the subarea interior contributions (the subarea has no contributing upstream subareas)
	    I3 = II
	  else
	    II = IDOA(I2) - 1 !the channel just upstream of the subarea
	    I3 = II+2 !after upstream loads have been added to the subarea
	  endif
        vol_runoff = QVOL(I3) !volume of surface runoff in the channel (m3)
        Q_runoff = vol_runoff / 86400.
        vol_total = WYLD(I3) !total volume of water in the channel (m3)
        Q_total = vol_total / 86400.
        !calculate concentration of each salt ion in the subarea channel
        load = 0.
        if(Q_total.gt.0.001) then
          do m=1,mion
            chan_csalt(i,m) = (SMQS_total(I3,m) * 1000.) / vol_total !kg --> g, / m3 = g/m3 = mg/L (i = actual subarea ID)
            load(m) = SMQS_total(I3,m)
            if(chan_csalt(i,m).lt.0) then
              chan_csalt(i,m) = 0.
            endif
            if(IBSA(i).eq.MSA) then !outlet subarea
              outlet_flow = Q_total
              outlet_csalt(m) = chan_csalt(i,m)
              outlet_load(m) = SMQS_total(I3,m)
            endif
          enddo
        endif
        
        !write out results for the current subarea
        write(2005,5000) i,IYR,IDA,area_ha,Q_runoff,Q_total,
     &                  (load(m),m=1,8), 
     &                  (chan_csalt(i,m),m=1,8) 
      
	enddo
      
      
      !write out results for the outlet subarea
      write(2006,5000) outlet_sub,IYR,IDA,outlet_area,outlet_flow,
     &                  (outlet_load(m),m=1,8), 
     &                  (outlet_csalt(m),m=1,8) 
      
      

      
 5000 format(I5,7x,I4,8x,I4,8x,500E16.8)
        
      return
      end


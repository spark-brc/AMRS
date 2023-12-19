      subroutine salt_budget !rtb salt
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    calculate and write out the salinity mass budget for the watershed
             
      use parm
      use GLOBAL,only:nrow,ncol,nlay,delr,delc,hnew,botm,ibound !MODFLOW
      use GWFRCHMODULE, only: RECH !MODFLOW
      use rt_global,only:cnew,prsity !RT3D
      use mf_rt_link,only:crch,rt_ncnh,CNH,DH !MODFLOW-RT3D link
      
      implicit none

      integer i,j,k,jj,num,m,sub,dum
      integer I1,I2,II,I3
      real    area_ha,salt_kg,volume,thickness,qss,flow_rate
      real    salt_sub(20,8),salt_subtot(20)
      real    salt_aqu(20),salt_wshd(20),salt_total(20)
         
          
      salt_total = 0.
      
      !subarea salt budget ----------------------------------------------------------------------------------------------------------------
      do i=1,MSA !loop through the subareas                                                                         
        
        salt_sub = 0.
        salt_subtot = 0.
        salt_kg = 0.
        
        sub = NBSA(i) !actual subarea
        area_ha = WSA(i) !subarea area in hectar
        
        !total salinity in soil profile (APEX)
        do jj=1,NBSL(i)
          ISL = LID(jj,i) !soil layer ID
          do m=1,mion !loop through the salt ions
            salt_sub(1,m) = salt_sub(1,m) + (wsalt(ISL,i,m)*WSA(i)) !kg/ha * ha = kg of salt
            salt_subtot(1)=salt_subtot(1) + (wsalt(ISL,i,m)*WSA(i))
          enddo
		enddo
        
        !sum up the salinity fluxes
        !loop through the salt ions
        do m=1,mion
          !track individual ions
          !(surfqsalt,surfsalt,latsalt,gwswsalt,swgwsalt not multipled by hectares because this was already done in BSUB, lines 1169-1172)
          salt_sub(2,m) = salt_sub(2,m) + surfqsalt_sub(sub,m)           !surface runoff 
          salt_sub(3,m) = salt_sub(3,m) + surfsalt_sub(sub,m)            !erosion runoff
          salt_sub(4,m) = salt_sub(4,m) + latsalt_sub(sub,m)             !lateral flow
          salt_sub(5,m) = salt_sub(5,m) + qrfsalt_sub(sub,m)             !quick return flow
          salt_sub(6,m) = salt_sub(6,m) + (irrsalt(i,m)*WSA(i))          !irrigation
          salt_sub(7,m) = salt_sub(7,m) +(saltdissolve_soil(i,m)*WSA(i)) !precip/dissolution
          salt_sub(8,m) = salt_sub(8,m) + (vsalt(i,m)*WSA(i))            !salinity flux: leaching from soil profile 
          salt_sub(9,m) = salt_sub(9,m) + gwswsalt_sub(sub,m)            !aquifer-->stream loading
          salt_sub(10,m) = salt_sub(10,m) + swgwsalt_sub(sub,m)            !stream-->aquifer loading
          salt_sub(11,m) = salt_sub(11,m) + pts_sub(sub,m)               !point source
          !track total salinity
          salt_subtot(2)=salt_subtot(2) + surfqsalt_sub(sub,m)           !surface runoff 
          salt_subtot(3)=salt_subtot(3) + surfsalt_sub(sub,m)            !erosion runoff
          salt_subtot(4)=salt_subtot(4) + latsalt_sub(sub,m)             !lateral flow
          salt_subtot(5)=salt_subtot(5) + qrfsalt_sub(sub,m)             !quick return flow
          salt_subtot(6)=salt_subtot(6) + (irrsalt(i,m)*WSA(i))          !irrigation
          salt_subtot(7)=salt_subtot(7) +(saltdissolve_soil(i,m)*WSA(i)) !precip/dissolution
          salt_subtot(8)=salt_subtot(8) + (vsalt(i,m)*WSA(i))            !salinity flux: leaching from soil profile 
          salt_subtot(9)=salt_subtot(9) + gwswsalt_sub(sub,m)            !aquifer-->stream loading
          salt_subtot(10)=salt_subtot(10) + swgwsalt_sub(sub,m)            !stream-->aquifer loading
          salt_subtot(11)=salt_subtot(11) + pts_sub(sub,m)               !point source
        enddo
        
        !write out salt ion mass for subarea
        write(2001,5000) sub,IYR,IDA,area_ha,(salt_sub(1,m),m=1,8),
     &                                       (salt_sub(2,m),m=1,8),
     &                                       (salt_sub(3,m),m=1,8),
     &                                       (salt_sub(4,m),m=1,8),
     &                                       (salt_sub(5,m),m=1,8),
     &                                       (salt_sub(6,m),m=1,8),
     &                                       (salt_sub(7,m),m=1,8),
     &                                       (salt_sub(8,m),m=1,8),
     &                                       (salt_sub(9,m),m=1,8),
     &                                       (salt_sub(10,m),m=1,8),
     &                                       (salt_sub(11,m),m=1,8)
        
        !write out salt budget array for subarea
        write(2002,5000) sub,IYR,IDA,area_ha,(salt_subtot(j),j=1,11)
        
        !add to the total salt array (summation of all subareas)
        salt_total(1) = salt_total(1) + salt_subtot(1) !soil salt
        salt_total(2) = salt_total(2) + salt_subtot(2) !surface runoff
        salt_total(3) = salt_total(3) + salt_subtot(3) !erosion runoff
        salt_total(4) = salt_total(4) + salt_subtot(4) !lateral flow
        salt_total(5) = salt_total(5) + salt_subtot(5) !quick return flow
        salt_total(6) = salt_total(6) + salt_subtot(6) !irrigation
        salt_total(7) = salt_total(7) + salt_subtot(7) !precipitation-dissolution
        salt_total(8) = salt_total(8) + salt_subtot(8) !leaching
        salt_total(9) = salt_total(9) + salt_subtot(9) !gw-sw discharge
        salt_total(10) = salt_total(10) + salt_subtot(10) !sw-gw seepage
        salt_total(11) = salt_total(11) + salt_subtot(11) !point sources

	enddo !go to next subarea
      
      
      
	
      !aquifer salt budget ----------------------------------------------------------------------------------------------------------------
      salt_aqu = 0.
      
      !total salt in aquifer (RT3D)
      do k=1,nlay
        do i=1,nrow
          do j=1,ncol
            thickness = 0
            if(ibound(j,i,k).ne.0) then
              volume = delr(j)*delc(i)*dh(j,i,k)*prsity(j,i,k) !m3 of groundwater
              do m=3,10
                salt_kg = volume * cnew(j,i,k,m) / 1000. !m3 * g/m3 / 1000 = kg
                salt_aqu(1) = salt_aqu(1) + salt_kg
              enddo
            endif
          enddo
        enddo
      enddo
      salt_total(12) = salt_aqu(1)
      
      !salinity flux: loading to water table
      do i=1,nrow
        do j=1,ncol
          !calculate total salt loading to the cell
          do m=3,10
            salt_kg = RECH(j,i) * CRCH(j,i,m) / 1000. !m3 * g/m3 / 1000 = kg
            salt_aqu(2) = salt_aqu(2) + salt_kg
          enddo
        enddo
      enddo
      salt_total(13) = salt_aqu(2)
      
      !salinity flux: aquifer<-->stream loading
      do i=1,MSA
        do m=1,mion
          if(sub_salt_exchange(i,m).gt.0) then
            salt_aqu(3) = salt_aqu(3) + sub_salt_exchange(i,m) !aquifer --> stream
          else
            salt_aqu(4) = salt_aqu(4) + (sub_salt_exchange(i,m)*-1) !stream --> aquifer
          endif
        enddo
      enddo
      salt_total(14) = salt_aqu(3) * (-1)
      salt_total(15) = salt_aqu(4)
      
      !salinity flux: boundary transport
      do num=1,rt_ncnh
        k = CNH(num,1)
        i = CNH(num,2)
        j = CNH(num,3)
        qss = CNH(num,4) !flow rate
        if(qss.lt.0) then 
          volume = delr(j)*delc(i)*dh(j,i,k) !bulk volume in m3
          flow_rate = qss * volume
          do m=1,mion
            salt_kg = flow_rate * cnew(j,i,k,m) / 1000. !m3 * g/m3 / 1000 = kg
            salt_aqu(5) = salt_aqu(5) + salt_kg
          enddo
        endif
      enddo
      salt_total(16) = salt_aqu(5)
      
      !write out salt_budget array for subarea
      write(2003,5001) IYR,IDA,(salt_aqu(j),j=1,5)
      
      
      
      
	!watershed salt budget --------------------------------------------------------------------------------------------------------------
      !make negative if the salt is leaving the system (soil profile or landscape)
      salt_total(2) = salt_total(2) * (-1) !surface runoff
      salt_total(3) = salt_total(3) * (-1) !erosion runoff
      salt_total(4) = salt_total(4) * (-1) !lateral flow
      salt_total(5) = salt_total(5) * (-1) !quick return flow
      salt_total(8) = salt_total(8) * (-1) !leaching
      salt_total(9) = salt_total(9) * (-1) !groundwater discharge
      write(2004,5001) IYR,IDA,(salt_total(j),j=1,16)
      
      

      
      !zero out arrays for the next day
      gwsalt = 0.
      surfqsalt = 0.
      surfsalt = 0.
      irrsalt = 0.
      latsalt = 0.
      qrfsalt = 0.
      vsalt = 0. !array that holds salt leaching from the bottom of the soil profile
      surfqsalt_sub = 0.
      surfsalt_sub = 0.
      latsalt_sub = 0.
      qrfsalt_sub = 0.
      gwswsalt_sub = 0.
      swgwsalt_sub = 0.
      pts_sub = 0.
      SMQS_total = 0.
      saltdissolve_soil = 0.
      
      
 5000 format(I5,7x,I4,8x,I4,8x,500E16.8)
 5001 format(I5,7x,I4,8x,50E16.8) 
        
      return
      end


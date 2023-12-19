      subroutine amrt_mfRiver
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine converts the APEX subarea channel depth to the MODFLOW River package
      
      !Import variables
      use parm, only:MSA,LPYR !APEX
      use GLOBAL, only:LENUNI,NCOL !MODFLOW
      use amrt_parm !APEX-MODFLOW
      use GWFRIVMODULE, only:RIVR,MXRIVR
      implicit none

      !Define local variables
      integer mf_lengthUnit,mf_timeUnit,i,j,k,subIndex,
     &        IC,IR,mf_gridID,subarea
      real rivlen(1),rivdepth(1),depth,length
      mf_lengthUnit = LENUNI + 10
      
      !output the channel depth for each APEX subarea
      if(out_APEX_channel.eq.1) then
        write(30003,1900) (dep_chan(j),j=1,MSA) !Ali        
      endif

      !output the river stage for each MODFLOW River Cell
      if(out_MF_channel.eq.1 .and.
     &   day_total.eq.outapexmf(apexmf_out_ctr)) then
        write(30004,*)
        write(30004,*) 'Day:',day_total
        write(30004,*) 'Channel depth for each MODFLOW River Cell'
      endif

      !Map the APEX subarea channel depth to MODFLOW's River Package
      do i=1,nrivcells_subs
        
        !reset averaged river properties
        rivlen = 0.
        rivdepth = 0.
        
        !calculate the channel depth for the River Cell
        !(loop through the APEX subareas connected with the current River Cell)
        do j=1,riv_nsubs(i)   
          !Get the river's properties to be based on a weighted average with:
          !weights = subbasin's river segment length / total river length in grid cell (rivlen)
          subarea = grid2riv_id(i,j) !APEX subarea connected with the current River Cell
          depth = dep_chan(subarea) !retrieve the channel depth for the APEX subarea
          length = grid2riv_len(i,j) !length of channel in the current River Cell
          rivlen(1) = rivlen(1) + length
          rivdepth(1) = rivdepth(1) + (depth * length)
        enddo
        
        !get weighted average channel depth
        if(rivlen(1).gt.0) then
          rivdepth(1) = rivdepth(1) / rivlen(1)
        else
          rivdepth(1) = rivdepth(1) / 1.
        endif

        !Convert to MODFLOW units
        call units(rivdepth, 12, mf_lengthUnit, 1, 1, LPYR);! to convert length units (m to LENUNI)
        
        !Populate MODFLOW's RIVR variable with channel depth ( = River Stage)
        do j=1,MXRIVR
          IR = RIVR(2,j) !row of MODFLOW river cell
          IC = RIVR(3,j) !column of MODFLOW river cell
          mf_gridID = ((IR-1)*NCOL) + IC !grid ID of MODFLOW river cell
          if(mf_gridID.eq.grid2riv_cellID(i)) then
            RIVR(4,j) = RIVR(6,j) + rivdepth(1) !stage = bottom elevation + depth
          endif
        enddo
        
        !Write out values to a file
        if(out_MF_channel.eq.1 .and.
     &     day_total.eq.outapexmf(apexmf_out_ctr)) then
          write(30004,*) rivdepth(1)
        endif
      enddo


 1900 format(1000(f12.4))
      
      return
      end

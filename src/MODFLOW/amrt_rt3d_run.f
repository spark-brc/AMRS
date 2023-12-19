      subroutine amrt_rt3d_run

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine calls conversions from APEX to RT3D, then runs RT3D
!!    for one day, then the groundwater results are converted back from RT3D
!!    to APEX in amrt_conversion2swat
    
      !Convert APEX variables units to RT3D variables and units
      call amrt_conversion2rt3d
          
      !Prepare RT3D
      call rtmf_prepare
          
      !Call RT3D
      call rt_run
          
      !There is no call to conversion back to APEX here because this is handled
      !in amrt_converstion2apex

      end subroutine amrt_rt3d_run

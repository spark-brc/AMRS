      subroutine amrt_init_rt3d

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine initializes RT3D linking subroutines of APEX-MODFLOW-RT3D


      !set up RT3D data and allocate arrays
      print *, 'RT3D is being used for groundwater solute transport'
      call rt_read
        
      end subroutine amrt_init_rt3d

      subroutine amrt_rt3d_stress  !Ali	  

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine calls RT3D's stress period reader

      use mf_rt_link, only: rt_kstp !MODFLOW-RT3D
      implicit none
          
      rt_kstp = 0
      call rt_stress

      end subroutine amrt_rt3d_stress
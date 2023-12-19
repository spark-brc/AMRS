      module amrt_parm
!!
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This module contains variables used for linking APEX variables to MODFLOW variables
      
      !flag for whether modflow is active -----------------------------------------------------------------------------
      !integer :: mf_active <-IMF in modparm by Jaehak
      
      !Declare global variables for APEX-MODFLOW linkage
      real, dimension (:,:), allocatable :: d2h_id, d2h_area
      integer :: d2h_size
      integer, dimension (:), allocatable :: cell_subareas

      !APEX-MODFLOW River segment variables
      integer :: nrivcells_subs
      integer, dimension (:), allocatable :: grid2riv_cellID
      real, dimension (:,:), allocatable :: grid2riv_id,grid2riv_len
      integer, dimension (:), allocatable :: riv_nsubs
      real, dimension (:), allocatable :: sub_gw_exchange,
     &                                    sub_gwno3_exchange,
     &                                    sub_gwp_exchange
      real, dimension (:), allocatable :: dep_chan,no3_chan,p_chan

      !groundwater storage
      real, dimension (:), allocatable :: sub_gw_volume,
     &                                    sub_gwno3,sub_gwp

      !Subarea --> Grid cell mapping
      type subarea_grid_mapping
        integer, allocatable :: subarea_id(:)
        real, allocatable :: subarea_area(:)
      endtype subarea_grid_mapping
      type (subarea_grid_mapping), dimension(:), allocatable :: s2g_map

      !Grid cell --> Subarea mapping
      type grid_subarea_mapping
        integer, allocatable :: cell_row(:)
        integer, allocatable :: cell_col(:)
        real, allocatable :: cell_perc(:)
      endtype grid_subarea_mapping
      type (grid_subarea_mapping), dimension(:), allocatable  :: g2s_map

      !specific yield variables
      real, dimension (:,:,:), allocatable :: mf_SC2
      integer :: mf_NTOP

      !APEX subarea variables
      real, dimension (:), allocatable :: etremain,sepbtm,
     &                                    rchrg,rch_volume,gw_delaye,
     &                                    rchrg_n,rchrg_p
      real, dimension (:), allocatable :: percn,percp
      integer :: nsubs
      
      !Variables for DRAIN package linkage to APEX subareas
      integer, dimension (:,:), allocatable :: drn_subs !rtb drain
      integer :: ndrn_subs,mf_drain_subs
      real, dimension (:), allocatable :: sub_drn

      !depth to water table
      real, dimension (:,:), allocatable :: wt_depth_cell
      real, dimension (:), allocatable :: sub_wt_depth

!!    Flags for writing output data (APEX, MODFLOW)
      integer :: out_APEX_recharge,out_MF_recharge,
     &           out_APEX_channel,out_MF_channel,
     &           out_MODFLOW_gwsw,out_APEX_gwsw

!!    Variables and Arrays for storing APEX-MODFLOW variable averages (monthly, annual) (rtb avg)
      integer :: day_count_mo,day_count_yr,amrt_month_counter,
     &           apexmf_out_avg
      
      integer, dimension (:), allocatable :: month_days,
     &                                       month_days_leap
      real, dimension (:,:), allocatable :: amrt_RECH_tot_mo,
     &                                      amrt_RECH_tot_yr
      real, dimension (:), allocatable :: amrt_RECH_APEX_tot_mo,
     &                                    amrt_RECH_APEX_tot_yr
      real, dimension (:), allocatable :: amrt_GWSW_MF_tot_mo,
     &                                    amrt_GWSW_MF_tot_yr
      real, dimension (:), allocatable :: amrt_GWSW_APEX_tot_mo,
     &                                    amrt_GWSW_APEX_tot_yr
      real, dimension (:,:,:,:), allocatable :: amrt_csolute_tot_mo,
     &                                          amrt_csolute_avg_mo,
     &                                          amrt_csolute_tot_yr,
     &                                          amrt_csolute_avg_yr
     
!!    Array for storing MODFLOW observation cells
      integer :: num_MF_obs,mf_obs_flag
      integer, dimension (:,:), allocatable :: MF_obs

!!    Array for storing APEX-MODFLOW output times
      integer :: n_outapexmf,apexmf_out_ctr
      integer, dimension (:), allocatable :: outapexmf

      !total simulation day counter
      integer :: day_total

      end module amrt_parm

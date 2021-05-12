! module type_definitions
!
! These are derived type definitions that are used in globals, setting, and
! elsewhere.
!
! Note that each of these is intended to be dependent only on derived
! types within this module.
!
! Types that are defined outside of here or setting_definition should
! be confined to the module in which they are defined.
!
!==========================================================================
module type_definitions

    implicit none

    ! EXTERNAL INFLOW OBJECT
    type totalInflow
        integer :: node_id ! index to element thar receives inflow
        ! t_series*sfactor + base_pat*baseline
        real(8), dimension(2) :: ext_t_series = [-1, -1] ! time_series
        real(8), dimension(2) :: ext_base_pat = [-1, -1] ! baseline pattern
        real(8) :: ext_baseline = 0 ! constant baseline value
        real(8) :: ext_sfactor = 0 ! time series scaling factor
        ! ---------------------------------------------------------
        real(8) :: dwf_avgValue = 0 ! average inflow value
        real(8), dimension(2) :: dwf_monthly_pattern = [-1, -1]
        real(8), dimension(2) :: dwf_daily_pattern = [-1, -1]
        real(8), dimension(2) :: dwf_hourly_pattern = [-1, -1]
        real(8), dimension(2) :: dwf_weekend_pattern = [-1, -1]
    end type totalInflow

end module type_definitions
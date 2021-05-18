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

    ! Temporal Array object
    type Array3D
        integer :: current_idx
        real(8), allocatable :: val(:,:,:)
    end type Array3D



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
        ! ---------------------------------------------------------
        ! fetch keeps track of the vars in nodeInflow that need to be
        ! fetched from SWMM C. If fetch[i] == .true. then the variable
        ! needs to be updated from SWMM C. Update as more data is needed.
        ! fetch = [
        !    ext_t_series, ext_base_pat, ext_baseline, ext_sfactor,
        !    dwf_avgValue, dwf_monthly_pattern, dwf_daily_pattern,
        !    dwf_hourly_pattern, dwf_weekend_pattern
        ! ]
        logical, dimension(9) :: fetch = .true.
    end type nodeInflow

end module type_definitions
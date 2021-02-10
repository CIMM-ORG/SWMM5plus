module objects

    use errors
    use dynamic_array

    implicit none

    public

    ! TABLE OBJECT
    type real_table
        integer :: table_type
        integer :: dim
        integer, allocatable :: tsize(:)
        type(real_array), allocatable :: data(:)
    end type real_table

    ! TIME SERIES OBJECT
    type tseries
        integer :: current = 1
        type(real_table) :: table
    end type tseries

    ! PATTERN OBJECT
    type pattern
        integer :: ptype
        integer :: count
        real, dimension(24) :: factor
    end type pattern

    ! EXTERNAL INFLOW OBJECT
    ! t_series*sfactor + base_pat*baseline
    type extInflow
        integer :: node_id ! index to element thar receives inflow
        integer :: t_series ! time_series
        integer :: base_pat ! pattern
        real :: baseline ! constant baseline value
        real :: sfactor ! time series scaling factor
        real :: max_inflow
    end type

    ! DRY INFLOW OBJECT
    type dwfInflow
        integer :: node_id ! index to element thar receives inflow
        real :: avgValue ! average inflow value
        integer :: monthly_pattern
        integer :: daily_pattern
        integer :: hourly_pattern
        integer :: weekly_pattern
        real :: max_inflow
    end type

    integer, private :: debuglevel = 0

    type(tseries), allocatable :: all_tseries(:)
    type(pattern), allocatable :: all_patterns(:)
    type(real_table), allocatable :: total_inflows(:)
    ! ----------------------------------------------
    type(extInflow), allocatable :: ext_inflows(:)
    type(dwfInflow), allocatable :: dwf_inflows(:)

end module objects
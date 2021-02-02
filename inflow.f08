module inflow

    use errors
    use objects
    use interface
    use dynamic_array
    use array_index

    implicit none

    type extInflow
        integer :: node_id ! index to element thar receives inflow
        type(tseries) :: t_series ! time_series
        type(pattern) :: base_pat ! pattern
        real :: baseline ! constant baseline value
        real :: sfactor ! time series scaling factor
        real :: max_inflow
    end type

    type dwfInflow
        integer :: node_id ! index to element thar receives inflow
        real :: avgValue ! average inflow value
        type(pattern) :: monthly_pattern
        type(pattern) :: daily_pattern
        type(pattern) :: hourly_pattern
        type(pattern) :: weekly_pattern
        real :: max_inflow
    end type

    type(extInflow), allocatable :: ext_inflows(:)
    type(dwfInflow), allocatable :: dwf_inflows(:)

    ! Temporal dynamics arrays
    type(integer_array) :: nodes_with_extinflow
    type(integer_array) :: nodes_with_dwfinflow

contains

    function get_max_inflow(k)
        integer, intent(in) :: k
        integer :: i, j, k
        type(extinflow) :: ext
        type(dwfinflow) :: dwf
        real :: get_max_inflow

        do i = 1, num_nodes
            j = nodeI(i, ni_extinflow)
            k = nodeI(i, ni_dwfinflow)
            if (j == -1 .and. k == -1) then
                cycle
            else if (j > -1 .and. k == -1) then
                ext = ext_inflows(j)
                get_max_inflow = maxval(ext%t_series%table%y%array * ext%sfactor + ext%baseline)
            else
                print *, MSG_FEATURE_NOT_COMPATIBLE
                stop
            end if
        end do
    end function get_max_inflow

    subroutine inflow_populate_inflows()
        integer :: num_extinflows
        integer :: num_dwfinflows
        integer :: extinflow_tseries
        integer :: extinflow_base_pat
        integer(4) :: patterns_idx
        integer :: i, j

        do i = 1, num_nodes
            if (get_node_attribute(i, node_has_extInflow) == 1) then
                call dyna_integer_append(nodes_with_extinflow, i)
            end if
            if (get_node_attribute(i, node_has_dwfInflow) == 1) then
                call dyna_integer_append(nodes_with_dwfinflow, i)
            end if
        end do

        allocate(ext_inflows(nodes_with_extinflow%len))
        allocate(dwf_inflows(nodes_with_dwfinflow%len))

        nodeI(:, ni_extinflow) = -1
        do j = 1, nodes_with_extinflow%len
            i = nodes_with_extinflow%array(j)
            ext_inflows(j)%node_id = i
            extinflow_tseries = get_node_attribute(i, node_extInflow_tSeries)
            if (extinflow_tseries .ne. -1) then
                ext_inflows(j)%t_series = get_inflow_tseries(extinflow_tseries)
            end if
            extinflow_base_pat = get_node_attribute(i, node_extInflow_basePat)
            if (extinflow_base_pat .ne. -1) then
                ext_inflows(j)%base_pat = get_pattern_factors(extinflow_base_pat)
            end if
            ext_inflows(j)%baseline = get_node_attribute(i, node_extInflow_baseline)
            ext_inflows(j)%sfactor =  get_node_attribute(i, node_extInflow_sFactor)
            nodeI(i, ni_extinflow) = j
        end do

        nodeI(:, ni_dwfinflow) = -1
        do j = 1, nodes_with_dwfinflow%len
            i = nodes_with_dwfinflow%array(j)
            dwf_inflows(j)%node_id = i
            dwf_inflows(j)%avgValue = get_node_attribute(i, node_dwfInflow_avgvalue)
            dwf_inflows(j)%monthly_pattern = get_pattern_factors(int(get_node_attribute(i, node_dwfInflow_monthly_pattern)))
            dwf_inflows(j)%daily_pattern = get_pattern_factors(int(get_node_attribute(i, node_dwfInflow_daily_pattern)))
            dwf_inflows(j)%hourly_pattern = get_pattern_factors(int(get_node_attribute(i, node_dwfInflow_hourly_pattern)))
            dwf_inflows(j)%weekly_pattern = get_pattern_factors(int(get_node_attribute(i, node_dwfInflow_weekly_pattern)))
            nodeI(i, ni_dwfinflow) = i
        end do
    end subroutine

    subroutine inflow_free_inflows()
        if (allocated(ext_inflows)) then
            deallocate(ext_inflows)
        end if
        if (allocated(dwf_inflows)) then
            deallocate(dwf_inflows)
        end if
    end subroutine inflow_free_inflows
end module inflow
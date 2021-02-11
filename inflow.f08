module inflow

    use errors
    use objects
    use globals
    use interface
    use dynamic_array
    use array_index
    use data_keys
    use selectors
    use setting_definition
    use datetime

    implicit none

    integer, private :: debuglevel = 0
    private :: inflow_get_maxima

contains

    function inflow_get_pattern_factor_at(p, date) result(pfactor)
        integer, intent(in) :: p
        real(8), intent(in) :: date
        real :: d2, pfactor
        integer :: dw, yy, mm, dd, d, h, s
        type(pattern) :: pat
        character(64) :: subroutine_name = 'inflow_get_pattern_factor_at'

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        dw = datetime_dayofweek(date)
        pfactor = 0

        if (p .ne. -1) then
            pat = all_patterns(p)
            if (pat%ptype == SWMM_MONTHLY_PATTERN) then
                call datetime_decodedate(date, yy, mm, dd)
                pfactor = pat%factor(mm)
            else if (pat%ptype == SWMM_DAILY_PATTERN) then
                pfactor = pat%factor(dw)
            else if ((pat%ptype == SWMM_HOURLY_PATTERN) .or. &
                ((pat%ptype == SWMM_WEEKEND_PATTERN) .and. ((dw == 1) .or. (dw == 7)))) then
                call datetime_decodetime(date, d, h, s)
                pfactor = pat%factor(h+1)
            endif
        endif
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end function inflow_get_pattern_factor_at

    function inflow_get_extinflow_at(e, date, current) result(flow)
        integer, intent(in) :: e ! extinflow id
        real(8), intent(in) :: date
        integer, optional, intent(in) :: current ! closest to date

        integer :: t, i
        real :: flow
        real :: pfactor = 0
        real :: tfactor = 0
        real :: t1, t2
        real :: y1, y2
        type(extinflow) :: ext
        type(tseries) :: ts
        real, allocatable, dimension(:) :: ttime, tvalues
        character(64) :: subroutine_name = 'inflow_get_extinflow_at'

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        if (date < swmm_start_time) then
            print *, MSG_INCORRECT_PARAMETER, 'date < start time'
            stop
        end if

        flow = 0
        ext = ext_inflows(e)
        t = ext%t_series
        ! Update pattern extinflow
        pfactor = inflow_get_pattern_factor_at(ext%base_pat, date)
        flow = flow + pfactor*ext%baseline

        if (t .ne. -1) then
            ts = all_tseries(t)
            ttime = ts%table%data(1)%array
            tvalues = ts%table%data(2)%array
            if (present(current)) then
                i = current
            else
                i = tables_find_time(ts%table, date)
            endif
            if (date == ttime(i)) then
                tfactor = tvalues(i)
            else ! interpolation is needed
                if (ttime(i) > date) then
                    if (i == 1) then
                        t1 = setting%time%starttime
                        y1 = tvalues(1)
                    else
                        t1 = ttime(i-1)
                        y1 = tvalues(i-1)
                    endif
                    t2 = ttime(i)
                    y2 = tvalues(i)
                else
                    if (i == ts%table%tsize(2)) then
                        t2 = setting%time%endtime
                        y2 = tvalues(i)
                    else
                        t2 = ttime(i+1)
                        y2 = tvalues(i+1)
                    endif
                    t1 = ttime(i)
                    y1 = tvalues(i)
                endif
                tfactor = interpolate(date, dble(t1), dble(t2), dble(y1), dble(y2), INTERPOLATION_LINEAR)
            endif
        endif

        flow = flow + tfactor*ext%sfactor
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end function inflow_get_extinflow_at

    function inflow_get_pattern_next_time_between(p, t1, t2) result(next_time)
        integer, intent(in) :: p
        real(8), intent(in) :: t1, t2

        real(8) :: next_time
        integer :: h1, m1, s1, h2, m2, s2
        type(pattern) :: pat
        character(64) :: subroutine_name = 'inflow_get_pattern_next_time_between'

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        if (p == -1) then
            next_time = t2
            return
        endif

        pat = all_patterns(p)

        if (pat%ptype == SWMM_MONTHLY_PATTERN) then
            next_time = datetime_get_next_month(t1)
        else if (pat%ptype == SWMM_DAILY_PATTERN) then
            next_time = datetime_get_next_day(t1)
        else if (pat%ptype == SWMM_HOURLY_PATTERN) then
            next_time = datetime_get_next_hour(t1)
        else if (pat%ptype == SWMM_WEEKEND_PATTERN) then
            next_time = datetime_get_next_weekendday_hour(t1)
        endif

        if (next_time > t2) then
            next_time = t2
        endif

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end function inflow_get_pattern_next_time_between

    subroutine inflow_get_maxima(nodeI, nodeR)
        ! needs to be executed after (inflow_load_inflows)
        integer, intent(in) :: nodeI(:,:)
        real, intent(inout) :: nodeR(:,:)
        integer :: i, j, k, t, node_id, tt, size_t
        type(tseries) :: ts
        type(extinflow) :: ext
        real(8), allocatable :: ttime(:)

        character(64) :: subroutine_name = 'inflow_get_maxima'

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        nodeR(:, nr_maxinflow) = 0
        do i = 1, nodes_with_inflow%len
            node_id = nodes_with_inflow%array(i)
            j = nodeI(node_id, ni_extinflow)
            k = nodeI(node_id, ni_dwfinflow)
            if (j == -1 .and. k > 0) then
                ! TODO - handle dry weather inflows
                cycle
            else if (j > 0 .and. k == -1) then
                ! I am assuming that there's no pattern in the external inflow
                ext = ext_inflows(j)
                ts = all_tseries(ext%t_series)
                size_t = ts%table%tsize(1)
                nodeR(node_id, nr_maxinflow) = maxval(ts%table%data(2)%array(1:size_t))
            else
                print *, MSG_FEATURE_NOT_COMPATIBLE
                stop
            end if
        end do

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine inflow_get_maxima

    subroutine inflow_load_inflows(nodeI, nodeR)
        integer, intent(inout) :: nodeI(:,:)
        real, intent(inout) :: nodeR(:,:)
        integer :: i, j
        logical :: l1, l2
        character(64) :: subroutine_name = 'inflow_load_inflows'

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        do i = 1, num_nodes
            l1 = get_node_attribute(i, node_has_extInflow) == 1
            l2 = get_node_attribute(i, node_has_dwfInflow) == 1
            if (l1) then
                call dyna_integer_append(nodes_with_extinflow, i)
            endif
            if (l2) then
                call dyna_integer_append(nodes_with_dwfinflow, i)
            endif
            if (l1 .or. l2) then
                call dyna_integer_append(nodes_with_inflow, i)
            endif
        end do

        allocate(ext_inflows(nodes_with_extinflow%len))
        allocate(dwf_inflows(nodes_with_dwfinflow%len))

        nodeI(:, ni_extinflow) = -1
        print *, "Retrieving External Inflows"
        do j = 1, nodes_with_extinflow%len
            print *, "Extinflow", j, '/', nodes_with_extinflow%len
            i = nodes_with_extinflow%array(j)
            ext_inflows(j)%node_id = i
            ext_inflows(j)%t_series = get_node_attribute(i, node_extInflow_tSeries)
            ext_inflows(j)%base_pat = get_node_attribute(i, node_extInflow_basePat)
            ext_inflows(j)%baseline = get_node_attribute(i, node_extInflow_baseline)
            ext_inflows(j)%sfactor =  get_node_attribute(i, node_extInflow_sFactor)
            nodeI(i, ni_extinflow) = j
        end do

        nodeI(:, ni_dwfinflow) = -1
        print *, "Retrieving Dry Inflows"
        do j = 1, nodes_with_dwfinflow%len
            print *, "Dryinflow", j, '/', nodes_with_dwfinflow%len
            i = nodes_with_dwfinflow%array(j)
            dwf_inflows(j)%node_id = i
            dwf_inflows(j)%avgValue = get_node_attribute(i, node_dwfInflow_avgvalue)
            dwf_inflows(j)%monthly_pattern = get_node_attribute(i, node_dwfInflow_monthly_pattern)
            dwf_inflows(j)%daily_pattern = get_node_attribute(i, node_dwfInflow_daily_pattern)
            dwf_inflows(j)%hourly_pattern = get_node_attribute(i, node_dwfInflow_hourly_pattern)
            dwf_inflows(j)%weekly_pattern = get_node_attribute(i, node_dwfInflow_weekly_pattern)
            nodeI(i, ni_dwfinflow) = i
        end do

        allocate(total_inflows(nodes_with_inflow%len))

        do i = 1, nodes_with_inflow%len
            total_inflows(i) = new_real_table(tinflow, 2)
        enddo

        call inflow_get_maxima(nodeI, nodeR)

        if ((debuglevel > 0) .or. (debuglevelall > 0)) then
            ! do i = 1, num_tseries
            !     print *, "TSERIES", i
            !     do j = 1, all_tseries(i)%table%tsize(1)
            !         print *, all_tseries(i)%table%data(1)%array(j), all_tseries(i)%table%data(2)%array(j)
            !     enddo
            ! enddo
            print *, '*** leave ',subroutine_name
        endif

    end subroutine inflow_load_inflows

    subroutine inflow_free_inflows()
        character(64) :: subroutine_name = 'inflow_free_inflows'

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name
        if (allocated(ext_inflows)) then
            deallocate(ext_inflows)
        endif
        if (allocated(dwf_inflows)) then
            deallocate(dwf_inflows)
        endif
        if (allocated(total_inflows)) then
            deallocate(total_inflows)
        endif
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine inflow_free_inflows
end module inflow
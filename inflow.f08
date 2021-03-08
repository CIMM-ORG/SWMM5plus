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
    use bc
    use tables

    implicit none

    integer, private :: debuglevel = 0

contains

    function inflow_get_pattern_factor_at(p, date) result(pfactor)
        integer, intent(in) :: p
        real(8), intent(in) :: date
        real(8) :: d2, pfactor
        integer :: dw, yy, mm, dd, h, m, s
        character(64) :: subroutine_name = 'inflow_get_pattern_factor_at'

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        dw = datetime_dayofweek(date)
        pfactor = 0

        if (p .ne. -1) then
            if (all_patterns(p)%ptype == SWMM_MONTHLY_PATTERN) then
                call datetime_decodedate(date, yy, mm, dd)
                pfactor = all_patterns(p)%factor(mm)
            else if (all_patterns(p)%ptype == SWMM_DAILY_PATTERN) then
                pfactor = all_patterns(p)%factor(dw)
            else if ((all_patterns(p)%ptype == SWMM_HOURLY_PATTERN) .or. &
                ((all_patterns(p)%ptype == SWMM_WEEKEND_PATTERN) .and. ((dw == 1) .or. (dw == 7)))) then
                call datetime_decodetime(date, h, m, s)
                pfactor = all_patterns(p)%factor(h+1)
            endif
        endif
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end function inflow_get_pattern_factor_at

    subroutine inflow_load_inflows(nodeI, nodeR, bcdataDn, bcdataUp)
        integer, intent(inout) :: nodeI(:,:)
        real(8), intent(inout) :: nodeR(:,:)
        type(bctype), allocatable, intent(inout) :: bcdataDn(:)
        type(bctype), allocatable, intent(inout) :: bcdataUp(:)
        integer :: i ! node index
        integer :: j ! loop over nodes with external inflow
        integer :: k ! loop over nodes with dry inflow inly
        integer :: ii ! total inflow index
        integer :: jj ! general purpose index
        integer :: min_res, tsize, ptype
        logical :: l1, l2
        character(64) :: subroutine_name = 'inflow_load_inflows'

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        ! Allocate boundary conditions
        ! Upstream Boundary conditions are stored in the following order:
        ! 1st - Nodes with external and with/without dry inflows
        ! 2nd - Nodes with dry inflows only
        ! All the nodes are associated to a total_inflow object which is also associated
        ! with a bctype object in BCUpstream whose total inflow is described by vectors
        ! TimeArray and ValueArray which add all the inflow types in SWMM as a single one

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
                nodeI(i, ni_node_type) = nBCup
            endif
        end do

        nodeI(1:N_BCdnstream, ni_temp1) = pack(nodeI(:,ni_idx),nodeI(:,ni_node_type) == nBCdn)
        N_BCdnstream = count(nodeI(:,ni_node_type) == nBCdn)
        N_BCupstream = count(nodeI(:,ni_node_type) == nBCup)

        call bc_allocate(bcdataDn, bcdataUp)

        print *, "Setting up BC downstream"
        do ii = 1, N_BCdnstream
            print *, "BC dnstream", ii, '/', N_BCdnstream
            bcdataDn(ii)%NodeID = nodeI(ii, ni_temp1)
            allocate(bcdataDn(ii)%TimeArray(2))
            allocate(bcdataDn(ii)%ValueArray(2))
            bcdataDn(ii)%TimeArray = (/dble(0.0), dble(setting%time%endtime)/)
            bcdataDn(ii)%ValueArray = nodeR(bcdataDn(ii)%NodeID, nr_Zbottom) + dble(3.0)
        enddo

        if (nodes_with_inflow%len == 0) then
            print *, "Error - There are no inflows"
            stop
        endif

        allocate(total_inflows(nodes_with_inflow%len))
        nodeI(:, ni_total_inflow) = -1
        min_res = 0 ! minimum resolution
        ii = 1

        do j = 1, nodes_with_extinflow%len
            i = nodes_with_extinflow%array(j)
            l2 = get_node_attribute(i, node_has_dwfInflow) == 1
            bcdataUp(ii)%NodeID = i
            nodeI(i, ni_total_inflow) = ii
            total_inflows(ii)%node_id = i
            total_inflows(ii)%ext_t_series = get_node_attribute(i, node_extInflow_tSeries)
            total_inflows(ii)%ext_base_pat = get_node_attribute(i, node_extInflow_basePat)
            total_inflows(ii)%ext_baseline = get_node_attribute(i, node_extInflow_baseline)
            total_inflows(ii)%ext_sfactor =  get_node_attribute(i, node_extInflow_sFactor)
            if (l2) then
                total_inflows(ii)%dwf_avgValue = get_node_attribute(i, node_dwfInflow_monthly_pattern)
                total_inflows(ii)%dwf_monthly_pattern = get_node_attribute(i, node_dwfInflow_monthly_pattern)
                total_inflows(ii)%dwf_daily_pattern = get_node_attribute(i, node_dwfInflow_daily_pattern)
                total_inflows(ii)%dwf_hourly_pattern = get_node_attribute(i, node_dwfInflow_hourly_pattern)
                total_inflows(ii)%dwf_weekend_pattern = get_node_attribute(i, node_dwfInflow_weekend_pattern)
                if (total_inflows(ii)%dwf_hourly_pattern /= -1) then
                    min_res = hourly
                else
                    if (total_inflows(ii)%dwf_weekend_pattern /= -1) then
                        min_res = weekend
                        if (total_inflows(ii)%dwf_daily_pattern /= -1) then
                            min_res = -daily
                        else
                            if (total_inflows(ii)%dwf_monthly_pattern /= -1) then
                                min_res = -monthly
                            endif
                        endif
                    else
                        if (total_inflows(ii)%dwf_daily_pattern /= -1) then
                            min_res = daily
                        else
                            if (total_inflows(ii)%dwf_monthly_pattern /= -1) then
                                min_res = monthly
                            endif
                        endif
                    endif
                endif
            endif

            if (total_inflows(ii)%ext_base_pat /= -1) then
                ptype = all_patterns(total_inflows(ii)%ext_base_pat)%ptype
                if (ptype == SWMM_MONTHLY_PATTERN) then
                    if (min_res == 0) then
                        min_res = monthly
                    endif
                else if (ptype == SWMM_DAILY_PATTERN) then
                    if ((min_res <= 0) .or. (min_res == monthly) .or. (min_res == weekend)) then
                        min_res = sign(1, min_res) * daily
                    endif
                else if (ptype == SWMM_HOURLY_PATTERN) then
                    min_res = hourly
                else if (ptype == SWMM_WEEKEND_PATTERN) then
                    if ((min_res > 0) .and. (min_res /= hourly)) then
                        min_res = -min_res
                    else if (min_res == 0) then
                        min_res = weekend
                    endif
                endif
            endif

            if (total_inflows(ii)%ext_t_series /= -1) then ! EXT INFLOW WITH TSERIES
                tsize = all_tseries(total_inflows(ii)%ext_t_series)%tsize(1)
                if (all_tseries(total_inflows(ii)%ext_t_series)%data(1)%array(tsize) < &
                    swmm_end_time) then
                        call tables_add_entry &
                        (all_tseries(total_inflows(ii)%ext_t_series), &
                        (/swmm_end_time, &
                        all_tseries(total_inflows(ii)%ext_t_series)%data(2)%array(tsize)/))
                endif
                ! The tseries is resampled in place
                if (min_res == -daily) then
                    call table_resample(all_tseries(total_inflows(ii)%ext_t_series), weekend)
                    call table_resample(all_tseries(total_inflows(ii)%ext_t_series), daily)
                else if (min_res == -monthly) then
                    call table_resample(all_tseries(total_inflows(ii)%ext_t_series), weekend)
                    call table_resample(all_tseries(total_inflows(ii)%ext_t_series), monthly)
                else if (min_res == monthly) then
                    call table_resample(all_tseries(total_inflows(ii)%ext_t_series), monthly)
                else if (min_res == daily) then
                    call table_resample(all_tseries(total_inflows(ii)%ext_t_series), daily)
                else if (min_res == hourly) then
                    call table_resample(all_tseries(total_inflows(ii)%ext_t_series), hourly)
                else if (min_res == weekend) then
                    call table_resample(all_tseries(total_inflows(ii)%ext_t_series), weekend)
                endif

                tsize = all_tseries(total_inflows(ii)%ext_t_series)%tsize(1)

                all_tseries(total_inflows(ii)%ext_t_series)%data(2)%array(1:tsize) = &
                    all_tseries(total_inflows(ii)%ext_t_series)%data(2)%array(1:tsize) * &
                    total_inflows(ii)%ext_sfactor !
                if (total_inflows(ii)%ext_base_pat /= -1) then ! EXT INFLOW WITH TSERIES AND WITH PATTERN
                    do jj = 1, tsize
                        all_tseries(total_inflows(ii)%ext_t_series)%data(2)%array(jj) &
                          = all_tseries(total_inflows(ii)%ext_t_series)%data(2)%array(jj) + &
                          inflow_get_pattern_factor_at(total_inflows(ii)%ext_base_pat, &
                            all_tseries(total_inflows(ii)%ext_t_series)%data(1)%array(jj)) &
                              * total_inflows(ii)%ext_baseline
                    enddo
                endif
            else ! EXT INFLOW WITHOUT TSERIES AND WITH PATTERN
                total_inflows(ii)%xy = new_real_table(tinflow, 2)
                call tables_add_entry(total_inflows(ii)%xy, (/swmm_start_time, total_inflows(ii)%ext_baseline/))
                call tables_add_entry(total_inflows(ii)%xy, (/swmm_end_time, total_inflows(ii)%ext_baseline/))
                if (min_res == -daily) then
                    call table_resample(total_inflows(ii)%xy, weekend)
                    call table_resample(total_inflows(ii)%xy, daily)
                else if (min_res == -monthly) then
                    call table_resample(total_inflows(ii)%xy, weekend)
                    call table_resample(total_inflows(ii)%xy, monthly)
                else if (min_res == monthly) then
                    call table_resample(total_inflows(ii)%xy, monthly)
                else if (min_res == daily) then
                    call table_resample(total_inflows(ii)%xy, daily)
                else if (min_res == hourly) then
                    call table_resample(total_inflows(ii)%xy, hourly)
                else if (min_res == weekend) then
                    call table_resample(total_inflows(ii)%xy, weekend)
                endif
                tsize = total_inflows(ii)%xy%tsize(1)
                total_inflows(ii)%xy%data(2)%array(1:tsize) = &
                    total_inflows(ii)%xy%data(2)%array(1:tsize) * total_inflows(ii)%ext_sfactor
                if (total_inflows(ii)%ext_base_pat /= -1) then
                    do jj = 1, tsize
                        total_inflows(ii)%xy%data(2)%array(jj) = total_inflows(ii)%xy%data(2)%array(jj) * &
                            inflow_get_pattern_factor_at &
                                (total_inflows(ii)%ext_base_pat, total_inflows(ii)%xy%data(1)%array(jj))
                    enddo
                endif
            endif

            if (total_inflows(ii)%ext_t_series /= -1) then ! EXT INFLOW WITH TSERIES
                allocate(bcdataUp(ii)%TimeArray(tsize))
                allocate(bcdataUp(ii)%ValueArray(tsize))
                bcdataUp(ii)%TimeArray(:) = &
                    all_tseries(total_inflows(ii)%ext_t_series)%data(1)%array(1:tsize)
                bcdataUp(ii)%ValueArray(:) = &
                    all_tseries(total_inflows(ii)%ext_t_series)%data(2)%array(1:tsize)
            else ! EXT INFLOW WITHOUT TSERIES
                allocate(bcdataUp(ii)%TimeArray(tsize))
                allocate(bcdataUp(ii)%ValueArray(tsize))
                bcdataUp(ii)%TimeArray(:) = total_inflows(ii)%xy%data(1)%array(1:tsize)
                bcdataUp(ii)%ValueArray(:) = total_inflows(ii)%xy%data(2)%array(1:tsize)
            endif
            ii = ii + 1
        end do

        do k = 1, nodes_with_dwfinflow%len
            l1 = get_node_attribute(i, node_has_extInflow) == 1
            i = nodes_with_dwfinflow%array(j)
            ! data associated to dry inflow was already extracted in
            ! previous loop for nodes with external and dry inflow
            if (.not. l1) then ! DRY INFLOW ONLY
                bcdataUp(ii)%NodeID = i
                nodeI(i, ni_total_inflow) = ii
                total_inflows(ii)%node_id = i
                total_inflows(ii)%dwf_avgvalue = get_node_attribute(i, node_dwfInflow_avgvalue)
                total_inflows(ii)%dwf_monthly_pattern = get_node_attribute(i, node_dwfInflow_monthly_pattern)
                total_inflows(ii)%dwf_daily_pattern = get_node_attribute(i, node_dwfInflow_daily_pattern)
                total_inflows(ii)%dwf_hourly_pattern = get_node_attribute(i, node_dwfInflow_hourly_pattern)
                total_inflows(ii)%dwf_weekend_pattern = get_node_attribute(i, node_dwfInflow_weekend_pattern)

                total_inflows(ii)%xy = new_real_table(tinflow, 2)
                call tables_add_entry(total_inflows(ii)%xy, (/swmm_start_time, dble(0)/))
                call tables_add_entry(total_inflows(ii)%xy, (/swmm_end_time, dble(0)/))
                min_res = 0

                if (total_inflows(ii)%dwf_monthly_pattern /= -1) then
                    min_res = monthly
                endif
                if (total_inflows(ii)%dwf_weekend_pattern /= -1) then
                    if (min_res == monthly) then
                        min_res = -monthly
                    else
                        min_res = monthly
                    endif
                endif
                if (total_inflows(ii)%dwf_daily_pattern /= -1) then
                    if ((min_res < 0) .or. (min_res == weekend)) then
                        min_res = -daily
                    else
                        min_res = daily
                    endif
                endif
                if (total_inflows(ii)%dwf_hourly_pattern /= -1) then
                    min_res = hourly
                endif

                if (min_res == -daily) then
                    call table_resample(total_inflows(ii)%xy, weekend)
                    call table_resample(total_inflows(ii)%xy, daily)
                else if (min_res == daily) then
                    call table_resample(total_inflows(ii)%xy, daily)
                else if (min_res == -monthly) then
                    call table_resample(total_inflows(ii)%xy, weekend)
                    call table_resample(total_inflows(ii)%xy, monthly)
                else if (min_res == monthly) then
                    call table_resample(total_inflows(ii)%xy, monthly)
                else if (min_res == hourly) then
                    call table_resample(total_inflows(ii)%xy, hourly)
                else if (min_res == weekend) then
                    call table_resample(total_inflows(ii)%xy, weekend)
                endif
                tsize = total_inflows(ii)%xy%tsize(1)
                allocate(bcdataUp(ii)%TimeArray(tsize))
                allocate(bcdataUp(ii)%ValueArray(tsize))
                bcdataUp(ii)%TimeArray(:) = total_inflows(ii)%xy%data(1)%array(1:tsize)
                bcdataUp(ii)%ValueArray(:) = total_inflows(ii)%xy%data(2)%array(1:tsize)
                ii = ii + 1
            endif
        enddo

        ii = 1
        nodeR(:, nr_maxinflow) = 0
        do ii = 1, nodes_with_inflow%len
            i = bcdataUp(ii)%NodeID
            ! DRY INFLOW WITH / WITHOUT EXT INFLOW
            do jj = 1, size(bcdataUp(ii)%ValueArray)
                bcdataUp(ii)%ValueArray(jj) = bcdataUp(ii)%ValueArray(jj) + &
                    inflow_get_pattern_factor_at(total_inflows(ii)%dwf_monthly_pattern, bcdataUp(ii)%TimeArray(jj)) * &
                    inflow_get_pattern_factor_at(total_inflows(ii)%dwf_daily_pattern, bcdataUp(ii)%TimeArray(jj)) * &
                    inflow_get_pattern_factor_at(total_inflows(ii)%dwf_hourly_pattern, bcdataUp(ii)%TimeArray(jj)) * &
                    inflow_get_pattern_factor_at(total_inflows(ii)%dwf_weekend_pattern, bcdataUp(ii)%TimeArray(jj)) * &
                    total_inflows(ii)%dwf_avgvalue
                ! Compute maximum inflow value
                if (bcdataUp(ii)%ValueArray(jj) > nodeR(i, nr_maxinflow)) then
                    nodeR(i, nr_maxinflow) = bcdataUp(ii)%ValueArray(jj)
                endif
            enddo
            ! Convert all times to seconds
            bcdataUp(ii)%TimeArray = (bcdataUp(ii)%TimeArray - swmm_start_time) * dble(secsperday)
            ! all_tseries are deallocated in interface.f08
            call free_table(total_inflows(ii)%xy)
        enddo

        if ((debuglevel > 0) .or. (debuglevelall > 0)) then
            ! do i = 1, nodes_with_inflow%len
                ! print *, "BCUPSTREAM"
                ! call print_object_name(bcdataUp(i)%NodeID, SWMM_NODE)
                ! print *, "Time", bcdataUp(i)%TimeArray(:)
                ! print *, "Array", bcdataUp(i)%ValueArray(:)
            ! enddo
            ! print *, "max inflow", nodeR(:, nr_maxinflow)
            print *, '*** leave ',subroutine_name
        endif

    end subroutine inflow_load_inflows
end module inflow
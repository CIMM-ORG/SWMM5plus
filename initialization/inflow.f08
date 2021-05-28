module inflow

    use setting_definition, only: setting
    use datetime

    implicit none

contains

    subroutine inflow_fetch()
    !-----------------------------------------------------------------------------
        integer :: i, idx, node_idx
        integer, allocatable :: nodes_with_inflow
        type(nodeInflow), pointer :: tmp_ninflow
    !-----------------------------------------------------------------------------

        idx = 1
        tmp => temporal_nodeInflow

        allocate(nodes_with_inflow(N_inflow))

        nodes_with_inflow = pack(nodeI(:,ni_idx), nodeYN(:,nYN_has_inflow), nodes_with_inflow)

        do i = 1, N_inflow
            node_idx = nodes_with_inflow(i)
            do while ()
                ! Extract values from SWMM C
                tmp(i, )
                total_inflows(i, idx, 0) = ! time
                total_inflows(i, idx, 1) = ! value
                if (idx == setting%Limiter%ArraySize%TemporalBC) then
                    idx = 1
                end if
                idx = idx + 1
            end do
        end do

        deallocate(nodes_with_inflow)
    end subroutine inflow_fetch

    function inflow_get_next_total(node_id)

    end function inflow_get_next_total

    function inflow_get_pattern_factor_at(pattern_id, date) result(pfactor)
        integer, intent(in) :: pattern_id
        real(8), intent(in) :: date
        real(8) :: d2, pfactor
        integer :: dw, yy, mm, dd, h, m, s
        character(64) :: subroutine_name = 'inflow_get_pattern_factor_at'

        if (setting%Debug%File%inflow) print *, '*** enter ',subroutine_name

        dw = datetime_dayofweek(date)
        pfactor = 0

        if (pattern_id /= -1) then
            if (all_patterns(pattern_id)%ptype == SWMM_MONTHLY_PATTERN) then
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
        if (setting%Debug%File%inflow)  print *, '*** leave ',subroutine_name
    end function inflow_get_pattern_factor_at

end module inflow
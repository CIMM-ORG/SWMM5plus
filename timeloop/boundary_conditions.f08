module boundary_conditions

    use interface
    use define_indexes
    use define_globals
    use utility, only: util_print_warning
    use define_settings, only: setting

    implicit none

    private

    public :: bc_step

contains

     subroutine bc_step()
        integer :: ii, nidx
        integer :: tstep_larger_than_resolution
        real(8) :: ttime, tnow, tend
        character(64) :: subroutine_name = "bc_step"
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%boundary_conditions)  print *, '*** enter ', subroutine_name

        tnow = setting%Time%Hydraulics%timeNow
        tend = setting%Time%EndTime

        if (N_flowBC > 0) then
            do ii = 1, N_flowBC
                nidx = BC%flowI(ii, bi_node_idx)
                if (BC%flowIdx(ii) == 0) then ! First fetch
                    call bc_fetch_flow(ii)
                else
                    ttime = BC%flowR_timeseries(ii, BC%flowIdx(ii), br_time) ! Current time slot (upper bound of time interval)
                    if (tnow > ttime) then ! Needs update
                        print *, "tnow:", tnow, "ttime", ttime
                        if (BC%flowIdx(ii) == setting%BC%BCSlots) then
                            call bc_fetch_flow(ii)
                        else
                            tstep_larger_than_resolution = -1
                            do while((tnow > ttime) .and. (BC%flowIdx(ii) < setting%BC%BCSlots))
                                BC%flowIdx(ii) = BC%flowIdx(ii) + 1
                                ttime = BC%flowR_timeseries(ii, BC%flowIdx(ii), br_time)
                                if ((BC%flowIdx(ii) == setting%BC%BCSlots) .and. (tnow > ttime)) then
                                    call bc_fetch_flow(ii)
                                    ttime = BC%flowR_timeseries(ii, BC%flowIdx(ii), br_time)
                                end if
                                tstep_larger_than_resolution = tstep_larger_than_resolution + 1
                            end do
                            if ((tstep_larger_than_resolution > 0) .and. setting%Warning) then
                                call util_print_warning("Warning (bc.f08): The flow boundary condition for node " &
                                // node%Names(nidx)%str // " has a higher resolution than the time step")
                            end if
                        end if
                    end if
                end if
            end do
        end if
        if (N_headBC > 0) then
            do ii = 1, N_headBC
                nidx = BC%headI(ii, bi_node_idx)
                if (BC%headIdx(ii) == 0) then ! First fetch
                    call bc_fetch_head(ii)
                else
                    ttime = BC%headR_timeseries(ii, BC%headIdx(ii), br_time) ! Current time slot (upper bound of time interval)
                    if (tnow > ttime) then ! Needs update
                        if (BC%headIdx(ii) == setting%BC%BCSlots) then
                            call bc_fetch_head(ii)
                        else
                            tstep_larger_than_resolution = -1
                            do while((tnow > ttime) .and. (BC%headIdx(ii) < setting%BC%BCSlots))
                                BC%headIdx(ii) = BC%headIdx(ii) + 1
                                ttime = BC%headR_timeseries(ii, BC%headIdx(ii), br_time)
                                if ((BC%headIdx(ii) == setting%BC%BCSlots) .and. (tnow > ttime)) then
                                    call bc_fetch_head(ii)
                                    ttime = BC%headR_timeseries(ii, BC%headIdx(ii), br_time)
                                end if
                                tstep_larger_than_resolution = tstep_larger_than_resolution + 1
                            end do
                            if ((tstep_larger_than_resolution > 0) .and. setting%Warning) then
                                call util_print_warning("Warning (bc.f08): The head boundary condition for node " &
                                // node%Names(nidx)%str // " has a higher resolution than the time step")
                            end if
                        end if
                    end if
                end if
            end do
        end if

        if (setting%Debug%File%boundary_conditions) then
            print *, "BC times"
            do ii = 1, setting%BC%BCSlots
                print *, BC%flowR_timeseries(:, ii, br_time)
            end do
            print *, "BC values"
            do ii = 1, setting%BC%BCSlots
                print *, BC%flowR_timeseries(:, ii, br_value)
            end do
            print *, '*** leave ', subroutine_name
        end if

    end subroutine bc_step

    subroutine bc_fetch_flow(bc_idx)
        integer, intent(in) :: bc_idx
        integer             :: ii, NN
        real(8)             :: new_inflow_time
        real(8)             :: new_inflow_value
        character(64)       :: subroutine_name = "bc_fetch_flow"
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%boundary_conditions)  print *, '*** enter ', subroutine_name

        NN = setting%BC%BCSlots

        if (BC%flowIdx(bc_idx) == 0) then ! First fetch
            BC%flowR_timeseries(bc_idx, 1, br_time) = setting%Time%StartTime
            BC%flowR_timeseries(bc_idx, 1, br_value) = interface_get_flowBC(bc_idx, setting%Time%StartTime)
        else ! last value becomes first
            BC%flowR_timeseries(bc_idx, 1, br_time) = BC%flowR_timeseries(bc_idx, NN, br_time)
            BC%flowR_timeseries(bc_idx, 1, br_value) = BC%flowR_timeseries(bc_idx, NN, br_value)
        end if

        do ii = 2, NN
            new_inflow_time = min(setting%Time%EndTime, interface_get_next_inflow_time(bc_idx, setting%Time%StartTime))
            BC%flowR_timeseries(bc_idx, ii, br_time) = new_inflow_time
            BC%flowR_timeseries(bc_idx, ii, br_value) = interface_get_flowBC(bc_idx, new_inflow_time)
            if (new_inflow_time == setting%Time%EndTime) exit
        end do
        BC%flowIdx(bc_idx) = 2

        if (setting%Debug%File%boundary_conditions) print *, '*** leave ', subroutine_name

    end subroutine bc_fetch_flow

    subroutine bc_fetch_head(bc_idx)
        integer, intent(in) :: bc_idx
        ! integer             :: ii
        ! real(8)             :: new_head

        ! if (BC%headIdx(bc_idx) == 0) then ! First fetch
        !     BC%headR_timeseries(bc_idx, 1, br_time) = setting%Time%StartTime
        ! else ! last value becomes first
        !     BC%headR_timeseries(bc_idx, 1, br_time) = BC%headR_timeseries(bc_idx, setting%BC%BCSlots, br_time)
        !     BC%headR_timeseries(bc_idx, 1, br_value) = BC%headR_timeseries(bc_idx, setting%BC%BCSlots, br_value)
        ! end if

        ! do ii = 2, setting%BC%BCSlots
        !     new_head = interface_get_next_head_time(setting%Time%StartTime)
        !     BC%headR_timeseries(bc_idx, ii, br_time) = new_head
        !     BC%headR_timeseries(bc_idx, ii, br_value) = interface_get_headBC(bc_idx, new_head)
        ! end do
        ! BC%headIdx(bc_idx) = 2
    end subroutine bc_fetch_head


end module boundary_conditions
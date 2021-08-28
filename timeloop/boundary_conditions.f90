module boundary_conditions

    use interface
    use define_indexes
    use define_globals
    use utility, only: util_print_warning
    use utility_interpolate
    use define_settings, only: setting
    use face, only: face_interpolate_bc
    implicit none

    private

    public :: bc_step
    public :: bc_update

contains

    subroutine bc_update()
        integer :: ii
        character(64) :: subroutine_name = "bc_update"
        !%-----------------------------------------------------------------------------

        if (setting%Debug%File%boundary_conditions)  print *, '*** enter ', this_image(), subroutine_name

        call bc_step()

        !% face interpolation for all BC faces
        if (N_flowBC > 0 .or. N_headBC > 0) then
            call bc_interpolate() ! computes interpolation
            call face_interpolate_bc() ! broadcast interpolation to face & elem arrays
        endif

        if (setting%Debug%File%boundary_conditions) then
            print *, "INFLOW BC"
            print *, "BC times"
            do ii = 1, setting%BC%slots
                print *, BC%flowR_timeseries(:, ii, br_time)
            end do
            print *, "BC values"
            do ii = 1, setting%BC%slots
                print *, BC%flowR_timeseries(:, ii, br_value)
            end do
            print *, '*** leave ', this_image(), subroutine_name

            print *, "HEAD BC"
            print *, "BC times"
            do ii = 1, setting%BC%slots
                print *, BC%headR_timeseries(:, ii, br_time)
            end do
            print *, "BC values"
            do ii = 1, setting%BC%slots
                print *, BC%headR_timeseries(:, ii, br_value)
            end do
            print *, '*** leave ', this_image(), subroutine_name
        end if
    end subroutine bc_update

    subroutine bc_step()
    !%-----------------------------------------------------------------------------
        integer :: ii, nidx, lidx
        integer :: tstep_larger_than_resolution
        real(8) :: ttime, tnow, tend
        character(64) :: subroutine_name = "bc_step"
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%boundary_conditions)  print *, '*** enter ', this_image(), subroutine_name

        tnow = setting%Time%Now
        tend = setting%Time%End

        if (N_flowBC > 0) then
            do ii = 1, N_flowBC
                nidx = BC%flowI(ii, bi_node_idx)
                if (BC%flowIdx(ii) == 0) then ! First fetch
                    call bc_fetch_flow(ii)
                else
                    ttime = BC%flowR_timeseries(ii, BC%flowIdx(ii), br_time) ! Current time slot (upper bound of time interval)
                    if (tnow > ttime) then ! Needs update
                        if (BC%flowIdx(ii) == setting%BC%slots) then
                            call bc_fetch_flow(ii)
                        else
                            tstep_larger_than_resolution = -1
                            do while((tnow > ttime) .and. (BC%flowIdx(ii) < setting%BC%slots))
                                BC%flowIdx(ii) = BC%flowIdx(ii) + 1
                                ttime = BC%flowR_timeseries(ii, BC%flowIdx(ii), br_time)
                                if ((BC%flowIdx(ii) == setting%BC%slots) .and. (tnow > ttime)) then
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
                    !% HACK - we are assuming that outfalls can only have one link upstream
                    lidx = node%I(nidx, ni_Mlink_u1)
                    link%R(lidx, lr_InitialDnstreamDepth) = BC%headR_timeseries(ii, 1, br_value) - node%R(nidx,nr_Zbottom)
                else
                    ttime = BC%headR_timeseries(ii, BC%headIdx(ii), br_time) ! Current time slot (upper bound of time interval)
                    if (tnow > ttime) then ! Needs update
                        if (BC%headIdx(ii) == setting%BC%slots) then
                            call bc_fetch_head(ii)
                        else
                            tstep_larger_than_resolution = -1
                            do while((tnow > ttime) .and. (BC%headIdx(ii) < setting%BC%slots))
                                BC%headIdx(ii) = BC%headIdx(ii) + 1
                                ttime = BC%headR_timeseries(ii, BC%headIdx(ii), br_time)
                                if ((BC%headIdx(ii) == setting%BC%slots) .and. (tnow > ttime)) then
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
    end subroutine bc_step

    subroutine bc_fetch_flow(bc_idx)
    !%-----------------------------------------------------------------------------
        integer, intent(in) :: bc_idx
        integer             :: ii, NN
        real(8)             :: new_inflow_time
        real(8)             :: new_inflow_value
        character(64)       :: subroutine_name = "bc_fetch_flow"
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%boundary_conditions)  print *, '*** enter ', this_image(), subroutine_name

        NN = setting%BC%slots

        if (BC%flowIdx(bc_idx) == 0) then ! First fetch
            BC%flowR_timeseries(bc_idx, 1, br_time) = setting%Time%Start
            BC%flowR_timeseries(bc_idx, 1, br_value) = interface_get_flowBC(bc_idx, setting%Time%Start)
        else ! last value becomes first
            BC%flowR_timeseries(bc_idx, 1, br_time) = BC%flowR_timeseries(bc_idx, NN, br_time)
            BC%flowR_timeseries(bc_idx, 1, br_value) = BC%flowR_timeseries(bc_idx, NN, br_value)
        end if

        new_inflow_time = setting%Time%Start
        do ii = 2, NN
            new_inflow_time = min(setting%Time%End, interface_get_next_inflow_time(bc_idx, new_inflow_time))
            BC%flowR_timeseries(bc_idx, ii, br_time) = new_inflow_time
            BC%flowR_timeseries(bc_idx, ii, br_value) = interface_get_flowBC(bc_idx, new_inflow_time)
            if (new_inflow_time == setting%Time%End) exit
        end do
        BC%flowIdx(bc_idx) = 2

        if (setting%Debug%File%boundary_conditions) then
            do ii = 1, NN
                write(*, "(*(G0.4 : ','))") BC%flowR_timeseries(bc_idx, ii, :)
            end do
            print *, '*** leave ', this_image(), subroutine_name
        end if

    end subroutine bc_fetch_flow

    subroutine bc_fetch_head(bc_idx)
    !%-----------------------------------------------------------------------------
        integer, intent(in) :: bc_idx
        integer             :: ii, NN
        real(8)             :: new_head_time
        real(8)             :: new_head_value
        character(64)       :: subroutine_name = "bc_fetch_head"
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%boundary_conditions)  print *, '*** enter ', this_image(), subroutine_name

        NN = setting%BC%slots

        if (BC%headIdx(bc_idx) == 0) then ! First fetch
            BC%headR_timeseries(bc_idx, 1, br_time) = setting%Time%Start
            BC%headR_timeseries(bc_idx, 1, br_value) = interface_get_headBC(bc_idx, setting%Time%Start)
        else ! last value becomes first
            BC%headR_timeseries(bc_idx, 1, br_time) = BC%headR_timeseries(bc_idx, NN, br_time)
            BC%headR_timeseries(bc_idx, 1, br_value) = BC%headR_timeseries(bc_idx, NN, br_value)
        end if

        do ii = 2, NN
            new_head_time = min(setting%Time%End, interface_get_next_head_time(bc_idx, setting%Time%Start))
            BC%headR_timeseries(bc_idx, ii, br_time) = new_head_time
            BC%headR_timeseries(bc_idx, ii, br_value) = interface_get_headBC(bc_idx, new_head_time)
            if (new_head_time == setting%Time%End) exit
        end do
        BC%headIdx(bc_idx) = 2

        if (setting%Debug%File%boundary_conditions) print *, '*** leave ', this_image(), subroutine_name
    end subroutine bc_fetch_head


    subroutine bc_interpolate()
    !%-----------------------------------------------------------------------------
    !% Description:
    !% This subroutine is for boundary condition interpolation.
    !% Base on the time passed from the time loop, we interpolate (linear interpolation for now)
    !% the boundary condition to get the corresponding value.
    !%-----------------------------------------------------------------------------
        real(8) :: tnow
        integer :: ii, slot_idx, upper_idx, lower_idx
        character(64) :: subroutine_name = 'bc_interpolate'
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%boundary_conditions)  print *, '*** enter ', this_image(), subroutine_name

        tnow = setting%Time%Now

        do ii=1, N_flowBC
            slot_idx = BC%flowIdx(ii)
            upper_idx = slot_idx
            lower_idx = slot_idx - 1
            !% Find the cloest index first, assign it to lower_idx for now
            if (BC%flowR_timeseries(ii, lower_idx, br_time) == tnow) then
                !% no need to do the interpolation, directly take the existing BC data
                BC%flowRI(ii) = BC%flowR_timeseries(ii, lower_idx, br_value)
            else if (lower_idx > 0) then
                if ( BC%flowR_timeseries(ii, lower_idx, br_value) .eq. BC%flowR_timeseries(ii, upper_idx, br_value)) then
                    BC%flowRI(ii) = BC%flowR_timeseries(ii, lower_idx, br_value)
                    !% constant value, no need to do the interpolation
                else
                    !% interpolation step
                    BC%flowRI(ii) = util_interpolate_linear( &
                        tnow, &
                        BC%flowR_timeseries(ii, lower_idx, br_time), &
                        BC%flowR_timeseries(ii, upper_idx, br_time), &
                        BC%flowR_timeseries(ii, lower_idx, br_value), &
                        BC%flowR_timeseries(ii, upper_idx, br_value))
                end if
            end if
        end do

        do ii=1, N_headBC
            slot_idx = BC%headIdx(ii)
            upper_idx = slot_idx
            lower_idx = slot_idx - 1
            !% Find the cloest index first, assign it to lower_idx just for now
            if (BC%headR_timeseries(ii, lower_idx, br_time) == tnow) then
                !% no need to do the interpolation, directly take the existing BC data
                BC%headRI(ii) = BC%headR_timeseries(ii, lower_idx, br_value)
            else if (lower_idx .ne. 0) then
                if ( BC%headR_timeseries(ii, lower_idx, br_value) .eq. BC%headR_timeseries(ii, upper_idx, br_value)) then
                    BC%headRI(ii) = BC%headR_timeseries(ii, lower_idx, br_value)
                    !% constant value, no need to do the interpolation
                else
                    BC%headRI(ii) = util_interpolate_linear( &
                        tnow, &
                        BC%headR_timeseries(ii, lower_idx, br_time), &
                        BC%headR_timeseries(ii, upper_idx, br_time), &
                        BC%headR_timeseries(ii, lower_idx, br_value), &
                        BC%headR_timeseries(ii, upper_idx, br_value))
                end if
            end if
        end do

        if (setting%Debug%File%boundary_conditions) print *, '*** leave ', this_image(), subroutine_name

    end subroutine bc_interpolate
end module boundary_conditions
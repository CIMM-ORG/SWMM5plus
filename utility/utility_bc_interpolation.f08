module utility_bc_interpolation

    use define_indexes
    use define_indexes
    use define_keys
    use define_globals
    use define_settings
    use utility_allocate

    implicit none

    private

    public :: interpolation_BC


    contains

    subroutine interpolation_BC(BCflow_interp_output, BChead_interp_output)
    !% ------------------------------
    !% Description: 
    !% This subroutine is for boundary condition interpolation. 
    !% Base on the time passed from the time loop, we interpolate (linear interpolation for now)
    !% the boundary condition to get the corresponding value.
        real(8) :: tnow
        integer :: ii, nidx, upper_idx = 0, lower_idx = 0
        real(8), intent(inout) :: BCflow_interp_output(:), BChead_interp_output(:)
        character(64) :: subroutine_name = 'interpolation_BC'
        !%-----------------------------------------------------------------------------

        if (setting%Debug%File%boundary_conditions)  print *, '*** enter ', subroutine_name

        tnow = setting%Time%Hydraulics%timeNow

        do ii=1, N_flowBC
            nidx = BC%flowIdx(ii)
            upper_idx = nidx
            lower_idx = nidx - 1 
            !% Find the cloest index first, assign it to lower_idx for now
            if (BC%flowR_timeseries(ii, lower_idx, br_time) == tnow) then 
                !% no need to do the interpolation, directly take the existing BC data
                BCflow_interp_output(ii) = BC%flowR_timeseries(ii, lower_idx, br_value)
            else if (lower_idx .ne. 0) then 
                if ( BC%flowR_timeseries(ii, lower_idx, br_value) .eq. BC%flowR_timeseries(ii, upper_idx, br_value)) then 
                    BCflow_interp_output(ii) = BC%flowR_timeseries(ii, lower_idx, br_value) 
                    !% constant value, no need to do the interpolation
                else
                    BCflow_interp_output(ii) = BC_flow_linear_interpolation(lower_idx, upper_idx, tnow, ii)
                    !% interpolation step
                end if
            end if
        end do

        do ii=1, N_headBC
            nidx = BC%headIdx(ii)
            upper_idx = nidx
            lower_idx = nidx - 1
            !% Find the cloest index first, assign it to lower_idx just for now 
            if (BC%headR_timeseries(ii, lower_idx, br_time) == tnow) then
                !% no need to do the interpolation, directly take the existing BC data
                BChead_interp_output(ii) = BC%headR_timeseries(ii, lower_idx, br_value)
            else if (lower_idx .ne. 0) then 
                if ( BC%headR_timeseries(ii, lower_idx, br_value) .eq. BC%headR_timeseries(ii, upper_idx, br_value)) then 
                    BChead_interp_output(ii) = BC%headR_timeseries(ii, lower_idx, br_value) 
                    !% constant value, no need to do the interpolation
                else
                    BChead_interp_output(ii) = BC_head_linear_interpolation(lower_idx, upper_idx, tnow, ii)
                end if
            end if
        end do

        if (setting%Debug%File%boundary_conditions) print *, '*** leave ', subroutine_name

    end subroutine interpolation_BC

    real(8) function BC_flow_linear_interpolation (lower_idx, upper_idx, current_T, BC_idx) result (outvalue)
    !% This is a linear interpolation function for BC interpoliation
        integer, intent(in) :: lower_idx, upper_idx, BC_idx
        real(8), intent(in) :: current_T
        real(8) :: ratio
        real(8), pointer :: flowBC(:,:)

        !%-------------------------------------------------------
        flowBC => BC%flowR_timeseries(BC_idx,:,:)

        ratio = (flowBC(upper_idx, br_time) - current_T)/(flowBC(upper_idx, br_time) - flowBC(lower_idx, br_time))
        outvalue = flowBC(lower_idx, br_value) + ratio * (flowBC(upper_idx, br_value) - flowBC(lower_idx, br_value))
    end function BC_flow_linear_interpolation
    
    real(8) function BC_head_linear_interpolation (lower_idx, upper_idx, current_T, BC_idx) result (outvalue)
    !% This is a linear interpolation function for BC interpoliation
        integer, intent(in) :: lower_idx, upper_idx, BC_idx
        real(8), intent(in) :: current_T
        real(8) :: ratio
        real(8), pointer :: headBC(:,:)

        !%-------------------------------------------------------
        headBC => BC%headR_timeseries(BC_idx,:,:)

        ratio = (headBC(upper_idx, br_time) - current_T)/(headBC(upper_idx, br_time) - headBC(lower_idx, br_time))
        outvalue = headBC(lower_idx, br_value) + ratio * (headBC(upper_idx, br_value) - headBC(lower_idx, br_value))
    end function BC_head_linear_interpolation

end module utility_bc_interpolation
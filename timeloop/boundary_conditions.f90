module boundary_conditions
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Converts from EPA-SWMM-C boundary condition input to BC structure
    !% for SWMM5+
    !%==========================================================================
    use interface_
    use define_indexes
    use define_keys
    use define_globals
    use utility, only: util_print_warning
    use utility_interpolate
    use define_settings, only: setting
    use face, only: face_interpolate_bc
    use define_xsect_tables
    use geometry
    use xsect_tables
    use utility_crash
    use utility_datetime, only: util_datetime_seconds_precision

    implicit none

    private

    public :: bc_step
    public :: bc_update

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine bc_update()
        !%------------------------------------------------------------------
        !% Description:
        !% Updates the BC data storage
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: ii
            logical :: isBConly = .true.
            character(64) :: subroutine_name = "bc_update"
        !%------------------------------------------------------------------
        !% Preliminaries
            !% --- ignore this subroutine if there are no BC
            if ((N_flowBC == 0) .and. (N_headBC == 0)) return    

            if (setting%Debug%File%boundary_conditions)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% --- ensure that we have BC data that brackets the current time
        call bc_step()

        !% --- interpolate the BC to the present time
        call bc_interpolate_flow()

        call bc_interpolate_head()

        !% --- store the BC value in face arrays
        call face_interpolate_bc(isBConly) 

        !% --- check for inflow BC that imply a large velocity if the
        !%     cross-section is full
        if (setting%Limiter%Velocity%UseLimitMaxYN) then
            do ii=1,N_flowBC
                ! print *, '=============================='
                ! print *, ii
                ! print *, 'bi_idx         ',BC%flowI(ii,bi_idx) 
                ! print *, 'bi_elem_idx    ',BC%flowI(ii,bi_elem_idx)
                ! print *, 'br_value       ',BC%flowR(ii,br_value)
                ! print *, 'full area      ',elemR(BC%flowI(ii,bi_elem_idx),er_FullArea)
                ! print *, ''
                if ((BC%flowR(ii,br_value) / elemR(BC%flowI(ii,bi_elem_idx),er_FullArea)) &
                    .ge. setting%Limiter%Velocity%Maximum) then
                    print *, 'USER CONFIGURATION ERROR'
                    print *, 'Implied inflow velocity exceeds the setting.Limiter.Velocity.Maximum'
                    print *, 'for one or more inflow elements (usually due to small area)'
                    print *, 'Error at BC bi_idx ',BC%flowI(ii,bi_idx)
                    print *, 'Node name  ', trim(node%Names(BC%flowI(ii,bi_node_idx))%str)
                    print *, 'Inflow value (CMS) = ',BC%flowR(ii,br_value)
                    print *, 'cross-sectional area of section ',elemR(BC%flowI(ii,bi_elem_idx),er_FullArea)
                    print *, 'Velocity limiter to use this inflow must be > ',BC%flowR(ii,br_value) / elemR(BC%flowI(ii,bi_elem_idx),er_FullArea)
                    print *, 'Current value of limiter ',setting%Limiter%Velocity%Maximum
                    call util_crashpoint(66987236)
                end if
            end do
        end if

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%boundary_conditions) then
                print *, "INFLOW BC"
                print *, "BC times"
                do ii = 1, setting%BC%TimeSlotsStored
                    print *, BC%flowTimeseries(:, ii, brts_time)
                end do
                print *, "BC values"
                do ii = 1, setting%BC%TimeSlotsStored
                    print *, BC%flowTimeseries(:, ii, brts_value)
                end do
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

                print *, "HEAD BC"
                print *, "BC times"
                do ii = 1, setting%BC%TimeSlotsStored
                    print *, BC%headTimeseries(:, ii, brts_time)
                end do
                print *, "BC values"
                do ii = 1, setting%BC%TimeSlotsStored
                    print *, BC%headTimeseries(:, ii, brts_value)
                end do

                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
           end if

    end subroutine bc_update
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine bc_step()
        !%------------------------------------------------------------------
        !% Description:
        !% Cycles through all the BC to make sure that the present time
        !% is between the upper location: BC%flowI(:,bi_TS_upper_idx, 
        !% BC%headI(:,bi_TS_upper_idx) and the
        !% lower location (upper location -1) in the BC%flowTimeseries
        !% and BC%headTimeseries arrays. If needed, this will fetch more 
        !% data from the timeseries files and overwrite the timeseries arrays
        !%------------------------------------------------------------------
        !% Declarations
            integer :: ii, nidx, lidx, mm
            integer :: interval_counter
            integer, pointer :: TimeSlotsStored, TS_upper_idx
            real(8), pointer :: tnow, tend, ttime
            character(64) :: subroutine_name = "bc_step"
        !%-----------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%boundary_conditions)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------
        !% Aliases:
            tnow            => setting%Time%Now
            tend            => setting%Time%End
            TimeSlotsStored => setting%BC%TimeSlotsStored
        !%-----------------------------------------------------------------

        !% --- flow BC 
        if (N_flowBC > 0) then
            !% --- cycle through all the flow BC
            do ii = 1, N_flowBC
                !% --- get the node index
                nidx = BC%flowI(ii, bi_node_idx)
                !% --- get the current location of the upper index in time series
                TS_upper_idx => BC%flowI(ii,bi_TS_upper_idx)
                !% --- check that this BC has an input file
                if (BC%flowYN(ii,bYN_read_input_series)) then
                    if (TS_upper_idx == 0) then 
                        !% --- first fetch of data from file
                        call bc_fetch_flow(ii)

                    else
                        !% --- get the current upper bound of time interval
                        ttime => BC%flowTimeseries(ii, TS_upper_idx, brts_time) 
                        !% --- check to see if we need to move to the next level of the BC data
                        if (tnow > ttime) then 
                            if (TS_upper_idx == TimeSlotsStored) then
                                !% --- if we're at the end, we need to fetch more data from file
                                call bc_fetch_flow(ii)
        
                            elseif (TS_upper_idx < TimeSlotsStored) then
                                !% --- there's still more stored BC data to cycle through
                                !% --- create a counter
                                interval_counter = -1
                                do while((tnow > ttime) .and. (TS_upper_idx < TimeSlotsStored))
                                    !% --- increment to the next stored time level
                                    TS_upper_idx = TS_upper_idx + 1
                                    !% --- get the time at the next level
                                    ttime => BC%flowTimeseries(ii, TS_upper_idx, brts_time)
                                    !% --- if tnow is still too large and we're at the end of storage, then fetch more data
                                    if ((TS_upper_idx == TimeSlotsStored) .and. (tnow > ttime)) then
                                        !% --- fetch more data and overwrite old storage
                                        call bc_fetch_flow(ii)
                                        !% --- use the new upper bound time value
                                        ttime => BC%flowTimeseries(ii, TS_upper_idx, brts_time)
                                    end if
                                    !% increment the counter
                                    interval_counter = interval_counter+ 1
                                end do
                                !% --- check if we had to go more than a single interval and print warning
                                if ((interval_counter > 0) .and. setting%Output%Warning) then
                                    call util_print_warning("Warning (bc_step): The flow boundary condition for node " &
                                    // trim(node%Names(nidx)%str) // " has smaller time intervals than the present model time step")
                                end if
                            else 
                                print *, 'CODE ERROR: unexpected else'
                                call util_crashpoint(582973)
                            end if
                        else
                            !% --- no action, update of BC storage not needed   
                        end if
                    end if
                    !% --- get the size of the time interval
                    BC%flowR(ii, br_timeInterval) =   BC%flowTimeseries(ii, TS_upper_idx,   brts_time) &
                                                    - BC%flowTimeseries(ii, TS_upper_idx-1, brts_time)                            
                else
                    !% --- HACK-future expansions should include getting BC from a data structure
                    !%     or external code through API
                    call util_print_warning("Error (bc.f08): The flow boundary condition for node " &
                    // node%Names(nidx)%str // " should always read from a time series")
                    call util_crashpoint( 87453)
                end if                 
            end do
        end if

        !% --- Head BC (outfall)
        if (N_headBC > 0) then
            !% --- cycle through all the head BC
            do ii = 1, N_headBC
                !% --- check that this BC has an input file
                if (BC%headYN(ii,bYN_read_input_series)) then
                    !% --- get the node index
                    nidx = BC%headI(ii, bi_node_idx)
                    !% --- get the upper index of the time series
                    TS_upper_idx => BC%headI(ii,bi_TS_upper_idx)
                    if (TS_upper_idx == 0) then 
                        !% --- first fetch of data from file
                        call bc_fetch_head(ii)
                        !%---  HACK - we are assuming that outfalls can only have one link upstream
                        !%     IMPORTANT -- WE NEED AN ERROR CHECK TO MAKE SURE THIS CONDITION ISN'T VIOLATED.
                        lidx = node%I(nidx, ni_Mlink_u1)
                        link%R(lidx, lr_InitialDnstreamDepth) = BC%headTimeseries(ii, 1, brts_value) - node%R(nidx,nr_Zbottom)
                    else
                        !% --- get the current upper bound of time interval
                        ttime => BC%headTimeseries(ii,TS_upper_idx, brts_time)

                        !% --- check to see if we need to move to the next level of the BC data
                        if (tnow > ttime) then 
                            if (TS_upper_idx == TimeSlotsStored) then
                                !% --- if we're at the end, we need to fetch more data from file
                                call bc_fetch_head(ii)
                            else
                                !% --- there's still more stored BC data to cycle through
                                !% --- create a counter
                                interval_counter = -1
                                do while((tnow > ttime) .and. (TS_upper_idx < TimeSlotsStored))
                                    !% --- increment to the next stored time level
                                    TS_upper_idx = TS_upper_idx + 1
                                    !% --- use the new upper bound time value
                                    ttime => BC%headTimeseries(ii, TS_upper_idx, brts_time)
                                    !% --- if tnow is still too large and we're at the end of storage, then fetch more data
                                    if ((TS_upper_idx == TimeSlotsStored) .and. (tnow > ttime)) then
                                        !% --- fetch more data and overwrite old storage
                                        call bc_fetch_head(ii)
                                        ttime => BC%headTimeseries(ii, TS_upper_idx, brts_time)
                                    end if
                                    !% increment the counter
                                    interval_counter = interval_counter + 1
                                end do
                                !% --- check if we had to go more than a single interval and print warning
                                if ((interval_counter > 0) .and. setting%Output%Warning) then
                                    call util_print_warning("Warning (bc_step): The head boundary condition for node " &
                                    // trim(node%Names(nidx)%str) // " has smaller time intervals than the present model time step")
                                end if
                            end if
                        end if
                    end if  

                    !% --- get the size of the time interval
                    BC%headR(ii, br_timeInterval) =   BC%headTimeseries(ii, TS_upper_idx,   brts_time) &
                                                    - BC%headTimeseries(ii, TS_upper_idx-1, brts_time)
                else
                    !% --- no BC file for this outfall -- BCH_fixed,..normal...free    
                    !% --- HACK-should have error checking that BC has appropriate setting
                end if
            end do
        end if

        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%boundary_conditions)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine bc_step
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine bc_fetch_flow(bc_idx)
        !%-------------------------------------------------------------------
        !% Description:
        !% Reads and stores data from a BC file.
        !% To prevent having to read the inflow file at every time step
        !% we store a number of "slots" of data controlled by the setting
        !% parameter "setting.BC.TimeSlotsStored". This code reads in and
        !% stores from the inflows in the SWMM input files
        !%-------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: bc_idx
            integer             :: ii
            integer, pointer    ::  NN, bc_level
            real(8)             :: new_inflow_time
            real(8)             :: new_inflow_value
            real(8)             :: tdummy
            real(8), pointer    :: timeEnd, timeEndEpoch
            character(64)       :: subroutine_name = "bc_fetch_flow"
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%boundary_conditions)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases:        
            NN           => setting%BC%TimeSlotsStored
            bc_level     => BC%flowI(bc_idx,bi_TS_upper_idx)
            timeEnd      => setting%Time%End
            timeEndEpoch => setting%Time%EndEpoch
        !%-------------------------------------------------------------------
        !% 
        !% --- check to see if this is the first fetch or a subsequent fetch.
        if (bc_level == 0) then 
            !% --- first fetch is always the simulation start time with interpolated value 
            BC%flowTimeseries(bc_idx, 1, brts_time)  = setting%Time%Start
            BC%flowTimeseries(bc_idx, 1, brts_value) = interface_get_flowBC(bc_idx, setting%Time%Start)
        else 
            !% --- on subsequent fetches we set the last value as the new first value
            BC%flowTimeseries(bc_idx, 1, brts_time)  = BC%flowTimeseries(bc_idx, NN, brts_time)
            BC%flowTimeseries(bc_idx, 1, brts_value) = BC%flowTimeseries(bc_idx, NN, brts_value)
        end if

        !% --- read in additional data to fill the timeseries arrays
        do ii = 2, NN
            !% --- get the next inflow time from the Time Series and advance the Tseries.x1 and Tseries.x2 locations
            !%     This uses the Epoch time as the last possible time (EPA-SWMM indexes of epoch time)
            new_inflow_time = interface_get_next_inflow_time(bc_idx, new_inflow_time,  timeEndEpoch)

            !% --- truncate the time in the table to the minimum of the end time and the next time
            new_inflow_time = min(timeEnd,new_inflow_time)

            !% --- set the timeseries to the new inflow time
            BC%flowTimeseries(bc_idx, ii, brts_time) = new_inflow_time

            !% --- get the new inflow value
            BC%flowTimeseries(bc_idx, ii, brts_value) = interface_get_flowBC(bc_idx, new_inflow_time)

            !% --- exit the loop if we've reached the maximum time for the simulation
            if (new_inflow_time == setting%Time%End) exit
        end do

        !% set the current location of the upper bound for interpolation to location 2
        BC%flowI(bc_idx,bi_TS_upper_idx) = 2

        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%boundary_conditions) then
                do ii = 1, NN
                    Print*, 'Time ', BC%flowTimeseries(bc_idx, ii, brts_time), 'Flow ', BC%flowTimeseries(bc_idx, ii, brts_value)
                end do
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            end if

    end subroutine bc_fetch_flow
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine bc_fetch_head(bc_idx)
        !%------------------------------------------------------------------
        !% Descriptions
        !% gets the bc head from EPA-SWMM-C
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: bc_idx
            integer             :: ii
            integer, pointer    :: NN, bc_level
            real(8)             :: new_head_time
            real(8)             :: new_head_value
            real(8)             :: tdummy
            real(8), pointer    :: timeEnd, timeEndEpoch
            character(64)       :: subroutine_name = "bc_fetch_head"
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%boundary_conditions)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases
            NN => setting%BC%TimeSlotsStored
            bc_level     => BC%headI(bc_idx,bi_TS_upper_idx)
            timeEnd      => setting%Time%End
            timeEndEpoch => setting%Time%EndEpoch
        !%-------------------------------------------------------------------

        !% --- check to see if this is the first fetch or a subsequent fetch.
        if (bc_level == 0) then 
            !% --- first fetch is always the simulation start time with interpolated value 
            BC%headTimeseries(bc_idx, 1, brts_time)  = setting%Time%Start
            BC%headTimeseries(bc_idx, 1, brts_value) = interface_get_headBC(bc_idx, setting%Time%Start) &
                                                         - setting%Solver%ReferenceHead
        else 
            !% --- on subsequent fetches we set the last value as the new first value
            BC%headTimeseries(bc_idx, 1, brts_time)  = BC%headTimeseries(bc_idx, NN, brts_time)
            BC%headTimeseries(bc_idx, 1, brts_value) = BC%headTimeseries(bc_idx, NN, brts_value)
        end if

        !% --- read in additional data to fill the timeseries arrays
        do ii = 2, NN
            !% --- get the next head time from the Time Series and advance the Tseries.x1 and Tseries.x2 locations
            !%     This uses the Epoch time as the last possible time (EPA-SWMM indexes of epoch time)
            new_head_time = interface_get_next_head_time(bc_idx, new_head_time,  timeEndEpoch)
            new_head_time  = util_datetime_seconds_precision (new_head_time)

            !% --- truncate the time in the table to the minimum of the end time and the next time
            new_head_time = min(timeEnd,new_head_time)

            !% --- set the timeseries to the new head time
            BC%headTimeseries(bc_idx, ii, brts_time)  = new_head_time
            !% --- get the new head value
            BC%headTimeseries(bc_idx, ii, brts_value) = interface_get_headBC(bc_idx, new_head_time) &
                                                             - setting%Solver%ReferenceHead

            !% --- exit the loop if we've reached the maximum time for the simulation
            if (new_head_time == setting%Time%End) exit
        end do
        !% set the current location of the upper bound for interpolation to location 2
        BC%headI(bc_idx,bi_TS_upper_idx) = 2

        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%boundary_conditions) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine bc_fetch_head
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine bc_interpolate_flow()
        !%-------------------------------------------------------------------
        !% Description:
        !% This subroutine is for flow boundary condition interpolation.
        !% Base on the time passed from the time loop, we interpolate (linear interpolation for now)
        !% the boundary condition to get the corresponding value.
        !%-------------------------------------------------------------------
        !% Declarations:
            real(8) :: normDepth, critDepth
            real(8), pointer :: tnow, flowValue(:)
            integer :: ii,  lower_idx , mm
            integer :: thisBCtype = BCFlow
            integer, pointer :: nodeIdx, faceIdx, elemUpIdx, upper_idx(:)
            character(64) :: subroutine_name = 'bc_interpolate_flow'
        !%-------------------------------------------------------------------
        !% Preliminaries:   
            if (setting%Debug%File%boundary_conditions)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases
            tnow      => setting%Time%Now
            flowValue   => BC%flowR(:,br_value)
            upper_idx => BC%flowI(:,bi_TS_upper_idx)
        !%-------------------------------------------------------------------    

        do ii=1, N_flowBC
            !% --- get the index below the current upper index
            lower_idx = upper_idx(ii) - 1

            call bc_interpolate_timeseries ( &
                    flowValue(ii), BC%flowTimeSeries, tnow, ii, lower_idx, upper_idx(ii), thisBCtype )

            !% --- HACK: do not let BC value to get smaller than zero
            !%     the absolute value is needed because of how max() handles very small differences.
            flowValue(ii) = abs(max(flowValue(ii),zeroR))
          
            !% --- error checking
            if (flowValue(ii) < zeroR) then
                print *, ' '
                print *, 'ILLEGAL NEGATIVE INFLOW BC ENCOUNTERED'
                print *, 'BC',ii,flowValue(ii)
                print *, upper_idx(ii), lower_idx
                print *, 'a', (BC%flowTimeseries(ii, lower_idx, brts_time) == tnow)
                print *, 'b', (lower_idx > 0)
                print *, 'c', ( BC%flowTimeseries(ii, lower_idx, brts_value) == BC%flowTimeseries(ii, upper_idx(ii), brts_value))
                print *, 'd', (.not. setting%BC%disableInterpolationYN)
                print *, 'tnow', tnow
                print *, BC%flowTimeseries(ii, lower_idx    , brts_time)
                print *, BC%flowTimeseries(ii, upper_idx(ii), brts_time)
                print *, BC%flowTimeseries(ii, lower_idx    , brts_value)
                print *, BC%flowTimeseries(ii, upper_idx(ii), brts_value)
                print *, 'result ',flowValue(ii)
                print *, ' ' 
                do mm=1,size(BC%flowTimeseries,DIM=2)
                    print *, BC%flowTimeseries(ii, mm,  brts_value)
                end do
                    
                call util_crashpoint(72098273)
            end if
            
        end do

        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%boundary_conditions) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine bc_interpolate_flow
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine bc_interpolate_head()
        !%-------------------------------------------------------------------
        !% Description:
        !% This subroutine is for head boundary condition interpolation.
        !% Base on the time passed from the time loop, we interpolate (linear interpolation for now)
        !% the boundary condition to get the corresponding value.
        !%-------------------------------------------------------------------
        !% Declarations:
            real(8) :: normDepth, critDepth, thisDepth, smallDepth
            real(8), pointer :: tnow, headValue(:), zbottom
            integer :: ii,  lower_idx , mm
            integer :: thisBCtype = BCHead
            integer, pointer :: nIdx, fIdx, upper_idx(:), eIdx
            
            character(64) :: subroutine_name = 'bc_interpolate_head'
        !%-------------------------------------------------------------------
        !% Preliminaries:   
            if (setting%Debug%File%boundary_conditions)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases
            tnow        => setting%Time%Now
            headValue   => BC%headR(:,br_value)
            upper_idx   => BC%headI(:,bi_TS_upper_idx)
        !%-------------------------------------------------------------------    
        !% --- cycle throuhg the head BC
        do ii=1, N_headBC

            nIdx        => BC%headI(ii,bi_node_idx)
            fIdx        => BC%headI(ii,bi_face_idx)
            eIdx        => BC%headI(ii,bi_elem_idx) !% upstream or downstream element

            zbottom     => faceR(fIdx,fr_Zbottom)
       
            !% --- select outfall type
            select case (BC%headI(ii,bi_subcategory))

                case (BCH_tidal)
                    print *, 'USER CONFIGURATION ERROR '
                    print *, 'Tidal BC not supported in SWMM5+'
                    call util_crashpoint(90222387)

                case (BCH_tseries)
                    !% --- get the index below the current upper index
                    lower_idx   =  upper_idx(ii) - 1

                    call bc_interpolate_timeseries ( &
                        headValue(ii), BC%headTimeSeries, tnow, ii, lower_idx, upper_idx(ii), thisBCtype )

                    thisDepth = max(headValue(ii) - zbottom, setting%ZeroValue%Depth)  

                case (BCH_fixed)
                    !% --- Note that SWMM.inp for FIXED BC is the elevation (not depth)
                    headValue(ii) = interface_get_headBC(ii, setting%Time%Start)

                    thisDepth = max(headValue(ii) - zbottom,setting%ZeroValue%Depth)  !% 20230511brh

                case (BCH_normal)
                    !% --- normal boundary handled in face_interp 

                case (BCH_free)
                    !% --- free boundary handled in face_interpolation 

                case default
                    call util_print_warning("CODE ERROR (bc_interpolate): Unknown downstream boundary condition type at " &
                        // trim(node%Names(nIdx)%str) // " node")
                    call util_crashpoint(86474)

            end select

        end do

        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%boundary_conditions) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            
    end subroutine bc_interpolate_head
    !%
!%==========================================================================
!%==========================================================================
!%
    subroutine bc_interpolate_timeseries &
        (interpout, TimeSeries, tnow, bc_idx, lower_idx, upper_idx, thisBCtype )
        !%------------------------------------------------------------------
        !% Description:
        !% Interpolates the timeseries for a single BC of either flow or head
        !% The upper index and lower index in the Time Series must bracket
        !% the time (this is checked elsewhere)
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(inout) :: interpout
            real(8), intent(in)    :: TimeSeries(:,:,:)
            real(8), intent(in)    :: tnow
            integer, intent(in)    :: bc_idx, lower_idx, upper_idx, thisBCtype
            integer :: ii
            character(64) :: subroutine_name = 'bc_interpolate_timeseries'
        !%------------------------------------------------------------------

        !% --- error checking
        if (lower_idx <= 0) then 
            !% lower_idx <= 0 is an error condition
            write(*,*), 'CODE ERROR: unexpected lower_idx value in ',trim(subroutine_name)
            call util_crashpoint( 987098)   
            return
        end if

        !%--- disabled interpolation: take the upper index value
        if (setting%BC%disableInterpolationYN) then
            interpout = Timeseries(bc_idx, upper_idx, brts_value)  
            return
        end if

        !% --- constant value in time series, no need to do the interpolation
        if (   Timeseries(bc_idx, lower_idx, brts_value)        &
            == Timeseries(bc_idx, upper_idx, brts_value) ) then 
            interpout = Timeseries(bc_idx, lower_idx, brts_value)
            return
        end if

        !% --- handle exact matches to time
        if (Timeseries(bc_idx, lower_idx, brts_time) == tnow) then                            
            !% ---  take the existing lower index BC data, no need to do the interpolation
            interpout = Timeseries(bc_idx, lower_idx, brts_value)
            return
        elseif (BC%flowTimeseries(bc_idx, upper_idx, brts_time) == tnow) then  
            !% --- take the existing upper index BC data, no need to do the interpolation
            interpout = Timeseries(bc_idx, upper_idx, brts_value)
            return
        end if
         
        !% --- standard interpolation
        interpout = util_interpolate_linear(          &
            tnow,                                      &
            Timeseries(bc_idx, lower_idx, brts_time),  &
            Timeseries(bc_idx, upper_idx, brts_time),  &
            Timeseries(bc_idx, lower_idx, brts_value), &
            Timeseries(bc_idx, upper_idx, brts_value))    
            
    
    end subroutine bc_interpolate_timeseries
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module boundary_conditions
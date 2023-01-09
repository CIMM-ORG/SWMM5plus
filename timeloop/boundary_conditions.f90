module boundary_conditions

    use interface_
    use define_indexes
    use define_keys
    use define_globals
    use utility, only: util_print_warning
    use utility_interpolate
    use define_settings, only: setting
    use face, only: face_interpolate_bc
    !use rectangular_triangular_conduit, only: rectangular_triangular_area_from_depth_singular, rectangular_triangular_topwidth_from_depth_singular
    use define_xsect_tables
    use geometry
    use xsect_tables
    use utility_crash
    use utility_datetime, only: util_datetime_seconds_precision
    ! use utility_unit_testing, only: util_utest_CLprint

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
        ! print *, 'calling bc_step'
        call bc_step()

            ! call util_utest_CLprint ('0000 after bc_step')

        !% --- interpolate the BC to the present time
        ! print *, 'calling bc_interpolate_flow'
        call bc_interpolate_flow()

            ! call util_utest_CLprint ('1111 after bc_interpolate_flow')

        ! print *, 'calling bc_interpolate_head'
        call bc_interpolate_head()

            ! call util_utest_CLprint ('2222 after bc_interpolate_head')

        !% --- store the BC value in face arrays
        ! print *, 'calling face_interpolate_bc'
        call face_interpolate_bc(isBConly) 

            ! call util_utest_CLprint ('3333 after face_interpolate_bc')

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
                    print *, 'CODE CONFIGURATION ERROR'
                    print *, 'Implied inflow velocity exceeds the setting.Limiter.Velocity.Maximum'
                    print *, 'for one or more inflow elements'
                    print *, 'Error at BC bi_idx ',BC%flowI(ii,bi_idx)
                    print *, 'Node name  ', trim(node%Names(BC%headI(ii,bi_node_idx))%str)
                    print *, 'Inflow value (CMS) = ',BC%flowR(ii,br_value)
                    print *, 'cross-sectional area of section ',elemR(BC%flowI(ii,bi_elem_idx),er_FullArea)
                    print *, 'Velocity limiter to use this inflow must be > ',BC%flowR(ii,br_value) / elemR(BC%flowI(ii,bi_elem_idx),er_FullArea)
                    print *, 'Current value of limiter ',setting%Limiter%Velocity%Maximum
                    call util_crashpoint(6698723)
                end if
            end do
        end if

        !  do ii=1,N_headBC
        !         print *, '=============================='
        !         print *, 'ii = ',ii
        !         print *, 'bi_idx         ',BC%headI(ii,bi_idx) 
        !         print *, 'bi_node_idx    ',BC%headI(ii,bi_node_idx), trim(node%Names(BC%headI(ii,bi_node_idx))%str)
        !         print *, 'bi_elem_idx    ',BC%headI(ii,bi_elem_idx)
        !         print *, 'bi_face_idx    ',BC%headI(ii,bi_face_idx)
        !         print *, 'bi_category    ',BC%headI(ii,bi_category), trim(reverseKey(BC%flowI(1,bi_category)))
        !         print *, 'bi_subcategory ',BC%headI(ii,bi_subcategory), trim(reverseKey(BC%flowI(1,bi_subcategory)))
        !         print *, 'bi_fetch       ',BC%headI(ii,bi_fetch)
        !         print *, 'bi_TS_upper_idx',BC%headI(ii,bi_TS_upper_idx)
        !         print *, 'br_value       ',BC%headR(ii,br_value)
        !         print *, 'face head_u    ',faceR(BC%headI(ii,bi_face_idx),fr_Head_u)
        !         print *, 'face head_d    ',faceR(BC%headI(ii,bi_face_idx),fr_Head_d)
        !         print *, 'downstream elem',faceI(BC%headI(ii,bi_face_idx),fi_Melem_dL)
        !         print *, 'link           ',elemI(faceI(BC%headI(ii,bi_face_idx),fi_Melem_dL),ei_link_Gidx_SWMM)
        !         !print *, 'link name      ',trim(link%Names(elemI(faceI(BC%headI(ii,bi_face_idx),fi_Melem_dL),ei_link_Gidx_SWMM))%str)
        !     end do

        !stop 44872

        ! do ii=1,N_flowBC
        !     print *, '=============================='
        !     print *, 'ii = ',ii
        !     print *, 'bi_idx         ',BC%flowI(ii,bi_idx) 
        !     print *, 'bi_node_idx    ',BC%flowI(ii,bi_node_idx),    trim(node%Names(BC%flowI(ii,bi_node_idx))%str)
        !     print *, 'bi_elem_idx    ',BC%flowI(ii,bi_elem_idx)
        !     print *, 'bi_face_idx    ',BC%flowI(ii,bi_face_idx)
        !     print *, 'bi_category    ',BC%flowI(ii,bi_category),    trim(reverseKey(BC%flowI(1,bi_category)))
        !     print *, 'bi_subcategory ',BC%flowI(ii,bi_subcategory), trim(reverseKey(BC%flowI(1,bi_subcategory)))
        !     print *, 'bi_fetch       ',BC%flowI(ii,bi_fetch)
        !     print *, 'bi_TS_upper_idx',BC%flowI(ii,bi_TS_upper_idx)
        !     print *, 'br_value       ',BC%flowR(ii,br_value)
        !     if (BC%flowI(ii,bi_face_idx) .ne. nullvalueI) then 
        !         print *, 'face flowrate  ',faceR(BC%flowI(ii,bi_face_idx),fr_Flowrate)
        !         print *, 'downstream elem',faceI(BC%flowI(ii,bi_face_idx),fi_Melem_dL)
        !         print *, 'link           ',elemI(faceI(BC%flowI(ii,bi_face_idx),fi_Melem_dL),ei_link_Gidx_SWMM)
        !         print *, 'link name      ',trim(link%Names(elemI(faceI(BC%flowI(ii,bi_face_idx),fi_Melem_dL),ei_link_Gidx_SWMM))%str)
        !     end if
        ! end do

        ! print *, ' =============='
        ! do ii=1,N_link
        !     print *, ' '
        !     print *, 'link ii =   ',ii, ' Name = ',trim(link%Names(ii)%str)
        !     print *, 'node up     ', link%I(ii,li_Mnode_u), ' Name = ',trim(node%Names(link%I(ii,li_Mnode_u))%str)
        !     print *, 'node dn     ', link%I(ii,li_Mnode_d), ' Name = ',trim(node%Names(link%I(ii,li_Mnode_d))%str)
        !     print *, 'elem index  ', link%I(ii,li_first_elem_idx), link%I(ii,li_last_elem_idx)
        !     print *, 'num elements', link%I(ii,li_N_element)
        !     print *, 'face up     ', elemI(link%I(ii,li_first_elem_idx),ei_Mface_uL)
        !     print *, 'face dn     ', elemI(link%I(ii,li_last_elem_idx),ei_Mface_dL)
        ! end do

        !print *, 'calling closing'
        !stop 330987
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
            !real(8) :: ttime
            integer, pointer :: TimeSlotsStored, TS_upper_idx
            real(8), pointer :: tnow, tend, ttime
            character(64) :: subroutine_name = "bc_step"
        !%-----------------------------------------------------------------
        !% Preliminaries
        !if (crashYN) return
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
                !print *, 'BC ii ',ii
                !% --- get the node index
                nidx = BC%flowI(ii, bi_node_idx)
                !print *, 'node ',nidx
                !% --- get the current location of the upper index in time series
                TS_upper_idx => BC%flowI(ii,bi_TS_upper_idx)
                !% --- check that this BC has an input file
                !print *, 'is inputfile ',BC%flowYN(ii,bYN_read_input_series)
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
                                print *, 'CODE ERROR: unexpected else at 582973'
                                call util_crashpoint(582973)
                            end if
                        else
                            !% --- no action, update of BC storage not needed   
                        end if
                    end if
                    !% --- get the size of the time interval
                    BC%flowR(ii, br_timeInterval) =   BC%flowTimeseries(ii, TS_upper_idx,   brts_time) &
                                                    - BC%flowTimeseries(ii, TS_upper_idx-1, brts_time)
                    
                    !print *, 'BC value ', BC%flowTimeseries(ii, TS_upper_idx, brts_value),  BC%flowTimeseries(ii, TS_upper_idx-1, brts_value)
                    ! print *, ' '
                    ! print *, 'BC flowR time interval ', ii, BC%flowR(ii, br_timeInterval)      
                    ! print *, ' '                             
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
                !print *, 'ii in HeadBC bc_step', ii
                !% --- check that this BC has an input file
                !print *, 'has file ',(BC%headYN(ii,bYN_read_input_series))
                if (BC%headYN(ii,bYN_read_input_series)) then
                    !% --- get the node index
                    nidx = BC%headI(ii, bi_node_idx)
                    !% --- get the upper index of the time series
                    TS_upper_idx => BC%headI(ii,bi_TS_upper_idx)
                    !print *, 'TS_upper_idx ',TS_upper_idx
                    if (TS_upper_idx == 0) then 
                        !% --- first fetch of data from file
                        call bc_fetch_head(ii)
                        !% HACK - we are assuming that outfalls can only have one link upstream
                        !% IMPORTANT -- WE NEED AN ERROR CHECK TO MAKE SURE THIS CONDITION ISN'T VIOLATED.
                        lidx = node%I(nidx, ni_Mlink_u1)
                        link%R(lidx, lr_InitialDnstreamDepth) = BC%headTimeseries(ii, 1, brts_value) - node%R(nidx,nr_Zbottom)
                    else
                        !% --- get the current upper bound of time interval
                        ttime => BC%headTimeseries(ii,TS_upper_idx, brts_time)

                        ! print *, 'times in bc_step'
                        ! print *, ttime, BC%headTimeseries(ii,TS_upper_idx-1, brts_time), ttime- BC%headTimeseries(ii,TS_upper_idx-1, brts_time)
                        ! print *, ' '

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

                    ! print *, ' '
                    ! print *, 'BC HeadR time interval ', ii, BC%headR(ii, br_timeInterval)      
                    ! print *, TS_upper_idx, TS_upper_idx-1
                    ! print *, BC%headTimeseries(ii, TS_upper_idx,   brts_time)
                    ! print *, BC%headTimeseries(ii, TS_upper_idx-1, brts_time)
                    ! print *, ' '                         
                    !stop 29873

                    ! print *, 'Time Series for Head'
                    ! do mm=1,TimeSlotsStored
                    !     print *, BC%headTimeseries(ii,mm,brts_time), BC%headTimeseries(ii,mm,brts_value)
                    ! end do
                    !stop 2098732
                else
                    !% --- no BC file for this outfall -- BCH_fixed,..normal...free    
                    !% --- HACK-should have error checking that BC has appropriate setting
                end if
            end do
        end if

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

        ! print *, ' '
        ! print *, 'in bc_fetch_flow'
        ! print *, 'bc_idx = ',bc_idx

        !% --- read in additional data to fill the timeseries arrays
        do ii = 2, NN
            !% --- get the next inflow time from the Time Series and advance the Tseries.x1 and Tseries.x2 locations
            !%     This uses the Epoch time as the last possible time (EPA-SWMM indexes of epoch time)
            new_inflow_time = interface_get_next_inflow_time(bc_idx, new_inflow_time,  timeEndEpoch)

            !print *, 'inflow time ',new_inflow_time, timeEndEpoch

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
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: bc_idx
            integer             :: ii
            integer, pointer    :: NN, bc_level
            real(8)             :: new_head_time
            real(8)             :: new_head_value
            real(8)             :: tdummy
            real(8), pointer    :: timeEnd, timeEndEpoch
            ! normDepth, critDepth
            !integer, pointer    :: nodeIdx, faceIdx, elemUpIdx
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
            ! print *, ' '
            ! print *, 'in AAA ',trim(subroutine_name)
            ! print *, BC%headTimeseries(bc_idx, 1, brts_time), BC%headTimeseries(bc_idx, 1, brts_value)
        else 
            !% --- on subsequent fetches we set the last value as the new first value
            BC%headTimeseries(bc_idx, 1, brts_time)  = BC%headTimeseries(bc_idx, NN, brts_time)
            BC%headTimeseries(bc_idx, 1, brts_value) = BC%headTimeseries(bc_idx, NN, brts_value)

                ! print *, ' '
                ! print *, 'in BBB ',trim(subroutine_name)
                ! print *, BC%headTimeseries(bc_idx, 1, brts_time), BC%headTimeseries(bc_idx, 1, brts_value)    
        end if

        !% --- read in additional data to fill the timeseries arrays
        do ii = 2, NN
            !% --- get the next head time from the Time Series and advance the Tseries.x1 and Tseries.x2 locations
            !%     This uses the Epoch time as the last possible time (EPA-SWMM indexes of epoch time)
            new_head_time = interface_get_next_head_time(bc_idx, new_head_time,  timeEndEpoch)
            new_head_time  = util_datetime_seconds_precision (new_head_time)
            !print *, 'ttime ',ttime

            !print *, ii, 'new_head_time', new_head_time, timeEndEpoch

            !new_head_time = min(setting%Time%End, interface_get_next_head_time(bc_idx, setting%Time%Start))

            !% --- truncate the time in the table to the minimum of the end time and the next time
            new_head_time = min(timeEnd,new_head_time)

            !print *, ii, 'limited new_head_time ',new_head_time

            !% --- set the timeseries to the new head time
            BC%headTimeseries(bc_idx, ii, brts_time)  = new_head_time
            !% --- get the new head value
            BC%headTimeseries(bc_idx, ii, brts_value) = interface_get_headBC(bc_idx, new_head_time) &
                                                             - setting%Solver%ReferenceHead

            ! print *, ' '
            ! print *, 'in CCC ',trim(subroutine_name)
            ! print *, BC%headTimeseries(bc_idx, ii, brts_time), BC%headTimeseries(bc_idx, ii, brts_value)

            !% --- exit the loop if we've reached the maximum time for the simulation
            if (new_head_time == setting%Time%End) exit
        end do
        !% set the current location of the upper bound for interpolation to location 2
        BC%headI(bc_idx,bi_TS_upper_idx) = 2


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

            !% HACK: do not let BC value to get smaller than zero
            !%       the absolute value is needed because of how max() handles very small differences.
            flowValue(ii) = abs(max(flowValue(ii),zeroR))

            !print *, 'in BC FLOW ',flowValue(ii)
            
            !% --- error checking
            ! if (lower_idx <= 0) then 
            !     !% lower_idx <= 0 is an error condition
            !     write(*,*), 'CODE ERROR: unexpected lower_idx value in ',trim(subroutine_name)
            !     call util_crashpoint( 987098)   
            !     return
            ! end if

            ! !%--- disabled interpolation: take the upper index value
            ! if (setting%BC%disableInterpolationYN) then
            !     flowValue(ii) = BC%flowTimeseries(ii, upper_idx(ii), brts_value)  
            !     return
            ! end if

            ! !% --- constant value in time series, no need to do the interpolation
            ! if (   BC%flowTimeseries(ii, lower_idx    , brts_value)        &
            !     == BC%flowTimeseries(ii, upper_idx(ii), brts_value) ) then 
            !     flowValue(ii) = BC%flowTimeseries(ii, lower_idx, brts_value)
            !     return
            ! end if

            ! !% --- handle exact matches to time
            ! if (BC%flowTimeseries(ii, lower_idx, brts_time) == tnow) then                            
            !     !% ---  take the existing lower index BC data, no need to do the interpolation
            !     flowValue(ii) = BC%flowTimeseries(ii, lower_idx, brts_value)
            !     return
            ! elseif (BC%flowTimeseries(ii, upper_idx(ii), brts_time) == tnow) then  
            !     !% --- take the existing upper index BC data, no need to do the interpolation
            !     flowValue(ii) = BC%flowTimeseries(ii, upper_idx(ii), brts_value)
            !     return
            ! end if
             
            ! !% --- standard interpolation
            ! flowValue(ii) = util_interpolate_linear(          &
            !     tnow,                                       &
            !     BC%flowTimeseries(ii, lower_idx    , brts_time),  &
            !     BC%flowTimeseries(ii, upper_idx(ii), brts_time),  &
            !     BC%flowTimeseries(ii, lower_idx    , brts_value), &
            !     BC%flowTimeseries(ii, upper_idx(ii), brts_value))                   
          
            !% --- error checking
            if (flowValue(ii) < zeroR) then
                print *, ' '
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
                    
                stop 449873
            end if
            
        end do
    
     !   stop 203987

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

            !% --- Error check: fail if offset present -- need to evaluate
            !%     Look at how Outfall offsets are used to set Z in link.c/link_setOutfallDepth
            !%     and Z is used in node.c/outfall_setOutletDepth
            !%  20221121 -- the use of zBottom from api_nodef_invertElev should solve this
            ! if ((link%R(node%I(nIdx,ni_Mlink_u1),lr_OutletOffset) > zeroR) .or. &
            !     (link%R(node%I(nIdx,ni_Mlink_u1),lr_InletOffset)  > zeroR) ) then
            !     print *, 'CODE ERROR: OUTFALL WITH OFFSET NEEDS ALGORITHM DESIGN'
            !     print *, 'outlet offset', link%R(node%I(nIdx,ni_Mlink_u1),lr_OutletOffset)
            !     print *, 'inlet offset' , link%R(node%I(nIdx,ni_Mlink_u1),lr_InletOffset)
            !     print *, 'this Node idx = ',nIdx, '; name =', trim(node%Names(nIdx)%str)
            !     print *, 'Link upstream idx =',node%I(nIdx,ni_Mlink_u1),'; name = ', trim(link%Names(node%I(nIdx,ni_Mlink_u1))%str)
            !     print *, 'Previously this error was seen when upstream of outfall was a weir or orifice'
            !     call util_crashpoint(6098734)
            ! end if

            critDepth = geo_criticaldepth_singular(BC%HeadI(ii,bi_UTidx))
            normDepth = geo_normaldepth_singular  (BC%HeadI(ii,bi_UTidx))
            
            !% --- select outfall type
            select case (BC%headI(ii,bi_subcategory))

            case (BCH_tidal)
                print *, 'ALGORITHM NEEDED FOR BCH_tidal ',trim(subroutine_name)
                call util_crashpoint(90222387)

            case (BCH_tseries)
                !% --- get the index below the current upper index
                lower_idx   =  upper_idx(ii) - 1

                call bc_interpolate_timeseries ( &
                    headValue(ii), BC%headTimeSeries, tnow, ii, lower_idx, upper_idx(ii), thisBCtype )

                thisDepth = max(headValue(ii) - zbottom, setting%SmallDepth%DepthCutoff)

                ! !% --- the time series are elevation (not depth)
                ! headValue(ii) = headValue(ii) - setting%Solver%ReferenceHead

                ! !% --- error check
                ! if (depth.le. zbottom + setting%SmallDepth%DepthCutoff) then
                !     print *, ' '
                !     print *, '*** WARNING ***'
                !     print *, '*** time series outfall stage is near the bottom at OUTFALL ',trim(node%Names(BC%headI(ii,bi_node_idx))%str)
                !     print *, '*** value is reset to smaller of normal or critical depth or'
                !     print *, '*** to the setting%SmallDepth%DepthCutoff'
                !     print *, ' '
                !     print *, 'head before fixing ',headValue(ii)
                !     print *, ' '
                !     ! critDepth = geo_criticaldepth_singular(BC%HeadI(ii,bi_UTidx))
                !     ! normDepth = geo_normaldepth_singular  (BC%HeadI(ii,bi_UTidx))

                !     !% --- BC head is the depth + Zbottom - referencehead
                !     headValue(ii) = zbottom        &
                !                   + max( min(critDepth,normDepth), setting%SmallDepth%DepthCutoff)      &
                !                   - setting%Solver%ReferenceHead

                !     !% -- note that time series head value is allowed to be larger than the
                !     !%    upstream value        

                ! end if

                ! print *, '*** head value at outfall ',headValue(ii)

            case (BCH_fixed)
                !% --- Note that SWMM.inp for FIXED BC is the elevation (not depth)
                headValue(ii) = interface_get_headBC(ii, setting%Time%Start)

                thisDepth = max(headValue(ii) - zbottom,setting%SmallDepth%DepthCutoff)

                ! !% --- BC head is the depth + Zbottom - referencehead
                ! headValue(ii) = zbottom  + thisDepth - setting%Solver%ReferenceHead                
                   
                ! !% --- error check
                ! if ((headValue(ii) + setting%Solver%ReferenceHead) .le. &
                !     (faceR(fIdx,fr_Zbottom) + setting%SmallDepth%DepthCutoff)) then
                !     print *, 'CONFIGURATION ERROR: a FIXED OUTFALL must have...'
                !     print *, '... a STAGE elevation greater than '
                !     print *, ' Zbottom + setting.SmallDepth.DepthCutff of the outfall'
                !     print *, 'Problem for Outfall ',trim(node%Names(BC%headI(ii,bi_node_idx))%str)
                !     call util_crashpoint(566823)
                ! end if

            case (BCH_normal)
                !% --- Error check, normal depth is infinite for a reverse flow, so a flap gate is needed
                if (.not. BC%headYN(ii,bYN_hasFlapGate)) then
                    print *, 'CONFIGURATION ERROR: a NORMAL OUTFALL must have flap gate set to YES'
                    print *, 'Problem for Outfall ',trim(node%Names(BC%headI(ii,bi_node_idx))%str)
                    call util_crashpoint(5668663)
                end if

                if (elemI(eIdx,ei_elementType) == CC) then
                    !% --- Error check, normal depth is infinite for adverse slope
                    !%     Use setting%Eps%Machine so that slope must be greater than precision 
                    if (elemR(eIdx,er_BottomSlope) .le. onehundredR*setting%Eps%Machine) then
                        print *, 'CONFIGURATION ERROR: a NORMAL OUTFALL must be connected to an...'
                        print *, '...conduit/channel element with non-zero, positive bottom slope.'
                        print *, 'Problem for Outfall ',trim(node%Names(BC%headI(ii,bi_node_idx))%str)
                        print *, 'Connected to element ', eIdx
                        print *, 'Part of link ',trim(  link%Names(elemI(eIdx,ei_link_Gidx_BIPquick))%str)
                        print *, 'Bottom slope is ',elemR(eIdx,er_BottomSlope)
                        call util_crashpoint(728474)
                    end if

                    !% --- depth is the larger of the normal depth or the small depth cutoff
                    thisDepth = max(normDepth, setting%SmallDepth%DepthCutoff)

                    ! !% --- BC head is the normal depth + Zbottom - referencehead
                    ! headValue(ii) = faceR(fIdx,fr_Zbottom)                   &
                    !               + thisDepth                                &
                    !               - setting%Solver%ReferenceHead

                    ! !% --- Outfall head should not be larger than the upstream
                    ! headValue(ii) = min(headValue(ii),elemR(eIdx,er_Head))    

                else
                    print *, 'CODE ERROR: NEED ALGORITHM DESIGN FOR OUTFALL WITH UPSTREAM DIAGNOSTIC ELEMENT'
                    call util_crashpoint(792873)
                    !% for free dnBC, if the upstream link is not CC (i.e. weir, orifice etc)
                    !% the depth in the node is zero
                    !headValue(ii) =  faceR(fIdx,fr_Zbottom)
                end if

                ! % --- Set final BC elevation following EPA SWMM in node.c/outfall_setOutletDepth
                ! %     Note that outfall z should already be included in zbottom
                ! %     and depth is already bounded on the low end by a min value    
                ! %     The BC must be larger than the smaller of the critical depth or normal depth

                ! print *, 'in ',trim(subroutine_name)
                ! print *, thisDepth, critDepth, normDepth

                smallDepth = min(critDepth,normDepth)
                if (thisDepth < smallDepth) then
                    !% --- Head should not be larger than the upstream
                    !%     to prohibit backwater and inflow when critical or normal depth controls.
                    headValue(ii) = min(smallDepth + zbottom - setting%Solver%ReferenceHead, &
                                        elemR(eIdx,er_Head)) 
                else
                    !% --- otherwise head based on thisDepth
                    headValue(ii) = thisDepth + zbottom - setting%Solver%ReferenceHead
                end if
    
                ! print *, 'head value ',headValue(ii)  !%- zbottom + setting%Solver%ReferenceHead

            case (BCH_free)
                !% --- Error check: Free outfall needs a flap gate, otherwise the negative flowrate into
                !%     the domain determines the outfall height, which can set up an instability
                if (.not. BC%headYN(ii,bYN_hasFlapGate)) then
                    print *, 'Problem for Outfall ',trim(node%Names(BC%headI(ii,bi_node_idx))%str)
                    write(*,*) '***************************************************************'
                    write(*,*) '** WARNING -- a FREE OUTFALL must have flap gate to prevent, **'
                    write(*,*) '**    back water flow. Thus, flap gate has been set to YES   **'
                    write(*,*) '***************************************************************'
                    BC%headYN(ii,bYN_hasFlapGate) = .true.
                    ! call util_crashpoint(566823)
                end if

                if (elemI(eIdx,ei_elementType) == CC) then
                    !% --- Use the smaller of critical depth or normal depth,
                    !%     but always use the depth cutoff as the smallest
                    !%     value allowed to prevent absurdly small depths
                    thisDepth = max(min(critDepth,normDepth), setting%SmallDepth%DepthCutoff)
                

                    ! !% --- BC head is the depth + Zbottom - referencehead
                    ! headValue(ii) = faceR(fIdx,fr_Zbottom)        &
                    !               + thisDepth     &
                    !               - setting%Solver%ReferenceHead

                    ! !% --- Head should not be larger than the upstream
                    ! headValue(ii) = min(headValue(ii),elemR(eIdx,er_Head))    

                else
                    print *, 'CODE ERROR: NEED ALGORITHM FOR OUTFALL WITH UPSTREAM DIAGNOSTIC ELEMENT'
                    call util_crashpoint(792873)
                    !% for free dnBC, if the upstream link is not CC (i.e. weir, orifice etc)
                    !% the depth in the node is zero
                    !headValue(ii) =  faceR(fIdx,fr_Zbottom)
                end if

            case default
                call util_print_warning("CODE ERROR (bc_interpolate): Unknown downstream boundary condition type at " &
                    // trim(node%Names(nIdx)%str) // " node")
                call util_crashpoint(86474)
            end select

            ! print *, 'here in bc_interpolate_head '
            ! print *, 'headvalue ',headValue(ii)

        end do

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

        ! print *, 'in ',trim(subroutine_name)
        ! print *, 'thisBC type ',trim(reverseKey(thisBCtype))
        ! print *, 'is disabled ',setting%BC%disableInterpolationYN

        !%--- disabled interpolation: take the upper index value
        if (setting%BC%disableInterpolationYN) then
            interpout = Timeseries(bc_idx, upper_idx, brts_value)  
            return
        end if

        !% --- constant value in time series, no need to do the interpolation
        if (   Timeseries(bc_idx, lower_idx, brts_value)        &
            == Timeseries(bc_idx, upper_idx, brts_value) ) then 
            interpout = Timeseries(bc_idx, lower_idx, brts_value)
            !print *, 'is constant value ',interpout
            return
        end if

        !% --- handle exact matches to time
        if (Timeseries(bc_idx, lower_idx, brts_time) == tnow) then                            
            !% ---  take the existing lower index BC data, no need to do the interpolation
            interpout = Timeseries(bc_idx, lower_idx, brts_value)
            ! print *, 'is exact  lower time match ',interpout
            return
        elseif (BC%flowTimeseries(bc_idx, upper_idx, brts_time) == tnow) then  
            !% --- take the existing upper index BC data, no need to do the interpolation
            interpout = Timeseries(bc_idx, upper_idx, brts_value)
            ! print *, 'is exact upper time match ',interpout
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
!%==========================================================================
!%
! function bc_get_CC_critical_depth(elemIdx) result (criticalDepth)
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% OBSOLETE 20220726 replaeced with geo_criticaldepth_from_velocity
    !     !% gets the critical depth of a CC element
    !     !% this code has been adapted from SWMM5-C code
    !     !%-----------------------------------------------------------------------------
    !         integer, intent(in)  :: elemIdx
    !         real(8)              :: criticalDepth
    !         integer, pointer     :: elemGeometry 
    !         real(8), pointer     :: Flowrate, FullDepth, BreadthMax, grav, FullArea
    !         real(8), pointer     :: alpha, bottomWidth
    !         real(8)              :: Q2g, critDepthEstimate, ratio, sideSlope, epsilon_c, tc0
    !         character(64) :: subroutine_name = 'bc_get_CC_critical_depth'
    !     !%-----------------------------------------------------------------------------
    !         !if (crashYN) return
    !         if (setting%Debug%File%boundary_conditions)  &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !         elemGeometry => elemI(elemIdx,ei_geometryType)
    !         Flowrate     => elemR(elemIdx,er_Flowrate)
    !         FullDepth    => elemR(elemIdx,er_FullDepth)
    !         FullArea     => elemR(elemIdx,er_FullArea)
    !         BreadthMax   => elemR(elemIdx,er_BreadthMax)
    !         grav         => setting%constant%gravity
    !         alpha        => setting%constant%energy_correction_factor

    !         Q2g = (Flowrate**twoR) / grav

    !         if (Q2g == zeroR) then
    !             criticalDepth =  zeroR
    !         else
    !             select case (elemGeometry)

    !             case (rectangular, rectangular_closed)
    !                 criticalDepth = (Q2g/BreadthMax**twoR)**onethirdR

    !             case (triangular)
    !                 sideSlope = elemSGR(elemIdx,esgr_Triangular_Slope)
    !                 criticalDepth = (twoR * Q2g / sideSlope ** twoR) ** onefifthR
    !             case (trapezoidal)
    !                 !% use the average side slope
    !                 sideSlope = onehalfR * (  elemSGR(elemIdx,esgr_Trapezoidal_LeftSlope)   &
    !                                         + elemSGR(elemIdx,esgr_Trapezoidal_RightSlope))
    !                 bottomWidth => elemSGR(elemIdx,esgr_Trapezoidal_Breadth)
    !                 !% Using approach of Vatakhah (2013)
    !                 !% non-dimensional discharge (epsilon_c), eq. 14
    !                 epsilon_c = fourR * sideSlope * ( alpha * (Flowrate**2) / (grav * (bottomWidth**5))  )**(onethirdR)
    !                 !% non-dimensional critical flow estimat (tc0), eq. 23
    !                 tc0 = (oneR + 1.161d0 * epsilon_c * (oneR + 0.666d0 * (epsilon_c**(1.041d0)) )**0.374d0 )**0.144d0
    !                 !% critical depth, eq. 22
    !                 criticalDepth = -onehalfR + onehalfR             &
    !                           * (                                    &                     
    !                                 (fiveR * (tc0**6) + oneR       ) &
    !                               / ( sixR * (tc0**5) - epsilon_c  ) &
    !                             )**3
        
    !             case (circular)  
    !                 !% first estimate Critical Depth for an equivalent circular conduit
    !                 critDepthEstimate = min(1.01*(Q2g/FullDepth)**onefourthR, FullDepth)

    !                 !% find ratio of conduit area to equiv. circular area
    !                 ratio = FullArea/(pi/fourR*(FullDepth**twoR))

    !                 if ((ratio >= onehalfR) .and. (ratio <= twoR)) then
    !                     criticalDepth = bc_critDepth_enum (elemIdx, Flowrate, critDepthEstimate)
    !                 else
    !                     criticalDepth = bc_critDepth_ridder (elemIdx, Flowrate, critDepthEstimate)
    !                 end if

    !             case default
    !                 print *, 'in ',trim(subroutine_name)
    !                 print *, 'CODE ERROR: unknown geometry of # ',elemGeometry 
    !                 print *, 'which has key ',trim(reverseKey(elemGeometry))
    !                 !stop 
    !                 call util_crashpoint( 389753)
    !                 !return
    !             end select
    !         end if

    !         if (setting%Debug%File%boundary_conditions) &
    !             write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

!e nd function bc_get_CC_critical_depth
!%
!%==========================================================================
!%==========================================================================
!%
! function bc_critDepth_enum (elemIdx, Q, yEstimate) result (criticalDepth)
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% gets the critical depth of a CC by enumeration method
    !     !% this code has been adapted from SWMM5-C code
    !     !%-----------------------------------------------------------------------------
    !         integer, intent(in)  :: elemIdx
    !         real(8), intent(in)  :: Q, yEstimate
    !         real(8)              :: criticalDepth
    !         real(8), pointer     :: FullDepth
    !         real(8)              :: dY, Y, Q0, QC
    !         integer              :: i1, ii
    !         character(64)        :: subroutine_name = 'bc_critDepth_enum'
    !     !%-----------------------------------------------------------------------------

    !     FullDepth => elemR(elemIdx,er_FullDepth)

    !     !% find the critical depth from SWMM5 ennumeration method
    !     dY = FullDepth/25.0
    !     i1 = int(yEstimate/dY)
    !     Q0 = bc_get_critical_flow(elemIdx,i1*dY,zeroR)

    !     if (Q0 < Q) then
    !         criticalDepth = FullDepth
    !         do  ii = i1+1,25
    !             QC = bc_get_critical_flow(elemIdx,ii*dY,zeroR)
    !             if (QC > Q) then
    !                 criticalDepth = ((Q-Q0)/(QC-Q0)+ real((ii-1),8))*dY
    !                 return
    !             end if
    !             Q0 = QC
    !         end do
    !     else
    !         criticalDepth = zeroR
    !         do  ii = i1-1,0,-1
    !             QC = bc_get_critical_flow(elemIdx,ii*dY,zeroR)
    !             if (QC < Q) then
    !                 criticalDepth = ((Q-QC)/(Q0-QC)+ real((ii),8))*dY
    !                 return
    !             end if
    !             Q0 = QC
    !         end do
    !     end if

! end function bc_critDepth_enum
!%
!%==========================================================================
!%==========================================================================
!%
! function bc_critDepth_ridder (elemIdx, Q, yEstimate) result (criticalDepth)
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% gets the critical depth of a CC by ridder method
    !     !% this code has been adapted from SWMM5-C code
    !     !%-----------------------------------------------------------------------------
    !         integer, intent(in)  :: elemIdx
    !         real(8), intent(in)  :: Q, yEstimate
    !         real(8)              :: criticalDepth, Y1, Y2
    !         real(8), pointer     :: FullDepth
    !         real(8)              :: Q0, Q1, Q2, Tol
    !         character(64)        :: subroutine_name = 'bc_critDepth_ridder'
    !     !%-----------------------------------------------------------------------------

    !     FullDepth => elemR(elemIdx,er_FullDepth)

    !     Y1 = zeroR
    !     Y2 = 0.99d0*FullDepth

    !     !% check if critical flow at full depth < target flow
    !     Q2 = bc_get_critical_flow(elemIdx,Y2,zeroR)
    !     if (Q2 < Q) then
    !         criticalDepth = FullDepth

    !     !% otherwise evaluate critical flow at initial estimate
    !     !% and 1/2 of fullDepth
    !     else
    !         Q0 = bc_get_critical_flow(elemIdx,yEstimate,zeroR)
    !         Q1 = bc_get_critical_flow(elemIdx, onehalfR*FullDepth,zeroR)

    !         if (Q0 > Q) then
    !             Y2 = yEstimate
    !             if (Q1 < Q) Y1 = onehalfR*FullDepth
    !         else
    !             Y1 = yEstimate
    !             if (Q1 > Q) Y2 = onehalfR*FullDepth
    !         end if
    !         Tol = 0.001
    !         criticalDepth = bc_findRoot_ridder_for_critDepth(elemIdx,Q,Y1,Y2,Tol)
    !     end if

! end function bc_critDepth_ridder
!%
!%==========================================================================
!%==========================================================================
!%
! function bc_findRoot_ridder_for_critDepth(elemIdx,Q,Y1,Y2,Tol) result (criticalDepth)
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% find root for the critical depth of a CC by ridder method
    !     !% this code has directly been adapted from SWMM5-C code
    !     !%-----------------------------------------------------------------------------
    !         integer, intent(in)  :: elemIdx
    !         real(8), intent(in)  :: Q, Y1, Y2, Tol
    !         real(8)              :: criticalDepth
    !         integer              :: ii, MaxIter
    !         real(8)              :: ans, fhi, flo, fm, fnew, s, Yhi, Ylo, Ym, Ynew

    !         character(64)        :: subroutine_name = 'bc_findRoot_ridder_for_critDepth'
    !     !%-----------------------------------------------------------------------------
    !         MaxIter = 60
    !         flo = bc_get_critical_flow(elemIdx,Y1,Q)
    !         fhi = bc_get_critical_flow(elemIdx,Y2,Q)

    !         if (flo == zeroR) then
    !             criticalDepth = Y1
    !             return
    !         end if

    !         if (fhi == zeroR) then
    !             criticalDepth = Y2
    !             return
    !         end if

    !         ans = onehalfR*(Y1+Y2)

    !         if ((flo > zeroR .and. fhi < zeroR) .or. (flo < zeroR .and. fhi > zeroR)) then
    !             Ylo = Y1
    !             Yhi = Y2
    !             do ii = 1,MaxIter
    !                 Ym = onehalfR*(Ylo+Yhi)
    !                 fm = bc_get_critical_flow(elemIdx,Ym,Q)
    !                 s = sqrt(fm*fm - flo*fhi)
    !                 if (s == zeroR) then
    !                     criticalDepth = ans
    !                     return
    !                 end if

    !                 if (flo >= fhi) then
    !                     Ynew = Ym + (Ym-Ylo)*fm/s
    !                 else
    !                     Ynew = Ym + (Ym-Ylo)*(-oneI*fm/s)
    !                 endif

    !                 if(abs(Ynew-ans) <= Tol) then 
    !                     criticalDepth = Ynew
    !                     return
    !                 end if
    !                 ans = Ynew
    !                 fnew = bc_get_critical_flow(elemIdx,Ynew,Q)

    !                 if (sign(fm,fnew) /= fm) then
    !                     Ylo = Ym
    !                     flo = fm
    !                     Yhi = ans
    !                     fhi = fnew
    !                 else if (sign(flo,fnew) /= flo) then
    !                     Yhi = ans
    !                     fhi = fnew
    !                 else if (sign(fhi,fnew) /= fhi) then
    !                     Ylo = ans
    !                     flo = fnew
    !                 else
    !                     criticalDepth = ans
    !                     return
    !                 end if

    !                 if ( abs(Yhi - Ylo) <= Tol) then
    !                     criticalDepth = ans
    !                     return
    !                 end if
    !             end do
    !             criticalDepth = ans
    !         end if

! end function bc_findRoot_ridder_for_critDepth
!%
!%==========================================================================
!%==========================================================================
!%
! function bc_get_critical_flow (elemIdx, Depth, Flow_0) result (Flow)
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% gets the critical flow of a CC element for a given depth
    !     !%-----------------------------------------------------------------------------
    !         integer, intent(in)  :: elemIdx
    !         real(8), intent(in)  :: Depth, Flow_0 
    !         real(8)              :: Flow
    !         integer, pointer     :: elemGeometry
    !         real(8), pointer     :: fullDepth, fullArea, grav
    !         real(8)              :: YoverYfull, Area, Topwidth
    !         character(64) :: subroutine_name = 'bc_get_critical_flow'
    !     !%-----------------------------------------------------------------------------

    !     elemGeometry => elemI(elemIdx,ei_geometryType)
    !     fullDepth    => elemR(elemIdx,er_FullDepth)
    !     fullArea     => elemR(elemIdx,er_FullArea)
    !     grav         => setting%constant%gravity

    !     select case (elemGeometry)

    !     case (circular)
    !         YoverYfull = Depth/fullDepth
    !         Area     = fullArea * xsect_table_lookup_singular (YoverYfull, ACirc ) !% 20220506brh removed NACirc
    !         Topwidth = max(fullDepth * xsect_table_lookup_singular (YoverYfull, TCirc), &   !% 20220506brh removed NTirc
    !                     setting%ZeroValue%Topwidth)

    !         Flow  = Area * sqrt(grav*Area/Topwidth) - Flow_0
        
    !     case (rect_triang)
    !         Area     = rectangular_triangular_area_from_depth_singular (elemIdx, Depth)
    !         Topwidth = rectangular_triangular_topwidth_from_depth_singular (elemIdx, Depth)
    !         Flow     = Area * sqrt(grav*Area/Topwidth) - Flow_0

    !     case default
    !         print *, 'in ',trim(subroutine_name)
    !         print *, 'CODE ERROR: unknown geometry of # ',elemGeometry 
    !         print *, 'which has key ',trim(reverseKey(elemGeometry))
    !         !stop 
    !         call util_crashpoint( 84198)
    !         !return
    !      end select

! end function bc_get_critical_flow
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module boundary_conditions
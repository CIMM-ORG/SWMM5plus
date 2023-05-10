module timeloop

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use output, only: outputML_store_data
    use pack_mask_arrays
    use runge_kutta2
    use hydrology
    use utility
    use utility_output
    use boundary_conditions
    use utility_profiler
    use utility_datetime
    use interface_, &
        only: interface_export_link_results, &
              interface_get_subcatch_runoff, &
              interface_call_runoff_execute, &
              interface_get_newRunoffTime,   &
              interface_export_runon_volume, &
              interface_call_climate_setState, &
              interface_get_evaporation_rate, &
              interface_get_RDII_inflow, &
              interface_get_groundwater_inflow
    use utility_crash
    use control_hydraulics, only: control_update
    ! use utility_unit_testing, only: util_utest_CLprint 

    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% top-level time-looping of simulation
    !%
    !% METHOD:
    !% Calls hydrology outer loop and then hydraulics inner loop (substeps)
    !%

    private
    public  :: timeloop_toplevel


contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine timeloop_toplevel()
        !%-------------------------------------------------------------------
        !% Description:
        !%     Loops over all the major time-stepping routines
        !%-------------------------------------------------------------------
        !% Declarations
            real(8)          :: startTime, endTime, reportStart
            logical          :: doHydraulicsStepYN, doHydrologyStepYN, SpinUpOutputYN
            real(8), pointer :: dtTol
            integer(kind=8)  :: cval, crate, cmax
            character(64)    :: subroutine_name = 'timeloop_toplevel'

        !%--------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !%--------------------------------------------------------------------
        !% Aliases
            dtTol   => setting%Time%DtTol
        !%--------------------------------------------------------------------
        !% --- store the start time so that we can reset after spin-up
        startTime   = setting%Time%Start
        endTime     = setting%Time%End
        reportStart = setting%Output%Report%StartTime 

        !print *, 'Starting timeloop toplevel '

        !% --- set spinup controls and call spinup
        call tl_spinup()

        !% --- if stopping after spin-up, end here
        if (setting%Simulation%stopAfterSpinUp) return

        !% --- reset the start time after spin-up
        setting%Time%Start = startTime
        setting%Time%Now   = startTime
        setting%Time%End   = endTime
        setting%Output%Report%StartTime = reportStart

        ! print *, 'BBB '

        sync all
        if (this_image()==1) then
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%TimeMarchStart = cval
            setting%Time%WallClock%LastTimeStored = cval
            setting%Time%WallClock%LastStepStored = setting%Time%Step
        end if 

        !print *, 'CCC '
        !% --- initialize the time settings for hydraulics and hydrology steps
        call tl_initialize_loop (doHydraulicsStepYN, doHydrologyStepYN, .false.)

        !print *, 'DDD '
        !-- perform the time-marching loop
        call tl_outerloop (doHydrologyStepYN, doHydraulicsStepYN, .false., .false.)

       ! print *, 'EEE '

        sync all
        !% --- close the timemarch time tick
        if (this_image() == 1) then
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%TimeMarchEnd= cval
        end if

        !print *, 'FFF '

        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine timeloop_toplevel
!% 
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine tl_spinup ()
        !% -----------------------------------------------------------------
        !% Description:
        !% Conducts spin-up simulations holding the time=0 BC fixed for the
        !% number of days set in setting%Simulation%SpinUpDays
        !% -----------------------------------------------------------------
        !% Declarations:
            logical :: inSpinUpYN, SpinUpOutputYN, doHydraulicsStepYN, doHydrologyStepYN
        !% -----------------------------------------------------------------

        if (.not. setting%Simulation%useSpinUp) return

        inSpinUpYN = .true.        
        setting%Time%Now = setting%Time%Start
        setting%Time%End = setting%Simulation%SpinUpDays * seconds_per_day

        !% --- only provide output during spinup if simulation stops after spinup  
        if (setting%Simulation%stopAfterSpinUp) then     
            SpinUpOutputYN = .true.
            !% default to report for all spinup
            setting%Output%Report%StartTime = setting%Time%Start
        else
            SpinUpOutputYN = .false.
            setting%Output%Report%StartTime = setting%Time%End + oneR
        end if

        !% --- initialize the time variables
        call tl_initialize_loop (doHydraulicsStepYN, doHydrologyStepYN, inSpinUpYN)
        
        !% --- perform the time loop for spin-up
        call tl_outerloop (doHydrologyStepYN, doHydraulicsStepYN, inSpinUpYN, SpinUpOutputYN)

    end subroutine tl_spinup
!% 
!%==========================================================================
!%==========================================================================
!%      
    subroutine tl_initialize_loop ( &
        doHydraulicsStepYN, doHydrologyStepYN, inSpinUpYN)
        !%------------------------------------------------------------------
        !% Description
        !% initialize the times before a time loop
        !%------------------------------------------------------------------
        !% Declarations
            logical, intent(inout) :: doHydrologyStepYN, doHydraulicsStepYN
            logical, intent(in)    :: inSpinUpYN
            real(8), pointer :: nextHydrologyTime, nextHydraulicsTime, nextControlRuleTime
            real(8), pointer :: lastHydrologyTime, lastHydraulicsTime, lastControlRuleTime, dtTol
            character(64)    :: subroutine_name = "tl_initialize_loop"    
        !%------------------------------------------------------------------
        !% Aliases  
            nextControlRuleTime => setting%Time%ControlRule%NextTime
            nextHydrologyTime   => setting%Time%Hydrology%NextTime
            nextHydraulicsTime  => setting%Time%Hydraulics%NextTime
            lastControlRuleTime => setting%Time%ControlRule%LastTime
            lastHydrologyTime   => setting%Time%Hydrology%LastTime
            lastHydraulicsTime  => setting%Time%Hydraulics%LastTime 
            dtTol               => setting%Time%DtTol 
        !%------------------------------------------------------------------

        !% --- set the last time storage
        lastHydrologyTime   = setting%Time%Start
        lastHydraulicsTime  = setting%Time%Start
        lastControlRuleTime = setting%Time%Start

        !% local T/F that changes with each time step depending on whether or 
        !% not the hydrology or hydraulics are conducted in that step
        doHydrologyStepYN  = setting%Simulation%useHydrology
        doHydraulicsStepYN = setting%Simulation%useHydraulics

        !% get the next hydrology time
        if (setting%Simulation%useHydrology) then      
            nextHydrologyTime  = interface_get_NewRunoffTime()
            ! print *, 'nextHydrologyTime from runoff = ',nextHydrologyTime
        else  
            !% set a large dummy time for hydrology if not used
            nextHydrologyTime  = setting%Time%End + onethousandR * dtTol
            ! print *, 'nextHydrologyTime from dummy = ',nextHydrologyTime
        end if
        
        if ((nextHydrologyTime == lastHydrologyTime) .and. &
            (setting%Simulation%useHydrology) ) then
            !% --- call the first hydrology step
            call tl_hydrology()
            nextHydrologyTime  = interface_get_NewRunoffTime()
        end if

        !% get the initial dt and the next hydraulics time
        if (setting%Simulation%useHydraulics) then
            call tl_smallestBC_timeInterval ()
            call tl_update_hydraulics_timestep(inSpinUpYN)
            call util_crashstop(229873)
        else
            !% NOTE -- WORKING WITHOUT SWMM5+ HYDRAULICS IS NOT SUPPORTED 
            !% the following is stub routines.
            !% set a large dummy time for hydraulics if not used
            nextHydraulicsTime = setting%Time%End + onethousandR * dtTol
            !% suppress the multi-level hydraulics output (otherwise seg faults)
            setting%Output%Report%suppress_MultiLevel_Output = .true.
            print *, 'CODE ERROR: SWMM5+ does not operate without hydraulics'
            call util_crashpoint(268743)
            call util_crashstop(268743)
        end if

        !% --- set the next control rule evaluation time
        nextControlRuleTime = lastControlRuleTime + real(setting%SWMMinput%ControlRuleStep,8)

        !% check to see if there is an initial step to hydrology
        !% if not, then skip the initial tl_hydrology
        if (( abs(nextHydrologyTime - setting%Time%Start) < dtTol) .and. (doHydrologyStepYN) ) then
            doHydrologyStepYN = .true.
        else
            doHydrologyStepYN = .false.
        endif 

    end subroutine tl_initialize_loop
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_outerloop ( &
        doHydrologyStepYN, doHydraulicsStepYN, inSpinUpYN, SpinUpOutputYN)
        !%------------------------------------------------------------------
        !% Description:
        !% do while loop for time marching
        !%-------------------------------------------------------------------
            logical, intent(inout)    :: doHydrologyStepYN, doHydraulicsStepYN
            logical, intent(in)       :: inSpinUpYN, SpinUpOutputYN
            real(8), pointer          :: dtTol
            integer(kind=8), pointer  :: thisStep
            integer(kind=8)           :: cval, crate, cmax
            integer :: ii
            logical :: BCupdateYN
        !%-------------------------------------------------------------------
        !% Aliases
            thisStep => setting%Time%Step
            dtTol    => setting%Time%DtTol 
        !%-------------------------------------------------------------------
        !% --- initialize the time step counter
        thisStep = 1

        !print *, 'starting tl_outerloop'

        do while (setting%Time%Now <= setting%Time%End - dtTol)
                ! print *, ' '
                ! print * , '==========================================================================='
                !print *, 'top of tl_outerloop loop at time ',setting%Time%Now
                ! print *, ' '
                    ! ! ! call util_utest_CLprint ('top of tl_outerloop')

                !% --- set the controls for using spin-up time
                if ((inSpinUpYN) .and. (thisStep > 1)) then
                    !% --- skip BC update during spin-up after first step
                    BCupdateYN = .false.
                else
                    BCupdateYN = .true.
                end if
    
                !% --- push the old values down the stack 
                !%     e.g. Flowrate to Flowrate_N0  and Flowrate_N0 to Flowrate_N1
                call tl_save_previous_values()

                    ! ! ! call util_utest_CLprint ('AAAAA before hydrology step in timeloop---------------------------')
    
                !% --- store the runoff from hydrology on a hydrology step
                if ((.not. inSpinUpYN) .or. (inSpinUpYN .and. BCupdateYN) ) then
                    !% --- note that hydrology in SWMM includes climate update
                    if (doHydrologyStepYN) then 
                        call tl_hydrology()
                        setting%Climate%EvapRate = interface_get_evaporation_rate()
                    else 
                        !% --- update climate if evaporation is needed in hydraulics-only simulations
                        if (setting%Climate%useHydraulicsEvaporationTF) then 
                            !print *, 'time now            ',setting%Time%Now
                            !print *, 'climate next update ',setting%Climate%NextTimeUpdate
                            if (setting%Time%Now .ge. setting%Climate%NextTimeUpdate) then 
                                !% --- update the time for the next climate call
                                setting%Climate%LastTimeUpdate = setting%Time%Now
                                setting%Climate%NextTimeUpdate = setting%Climate%LastTimeUpdate &
                                    + setting%Climate%HydraulicsOnlyIntervalHours * 3600.d0
                                call interface_call_climate_setState()    
                            end if
                        else 
                            !% --- no update to evaporation rate
                        end if
                    end if
                end if

                    ! ! ! call util_utest_CLprint ('BBBBB of tl_outerloop')

                !% --- main hydraulics time step
                if (doHydraulicsStepYN) then    

                    ! print *, 'in hydraulics step'
                    !% --- set a clock tick for hydraulic loop evaluation
                    if ((this_image()==1) .and. (.not. inSpinUpYN)) then
                        call system_clock(count=cval,count_rate=crate,count_max=cmax)
                        setting%Time%WallClock%HydraulicsStart = cval
                    end if 

                        ! ! ! call util_utest_CLprint ('CCCCC of tl_outerloop')

                    !% --- get updated boundary conditions
                    if (BCupdateYN) then
                        ! ! ! call util_utest_CLprint ('CCCCC_00 of tl_outerloop')

                        call bc_update() 

                            ! call util_utest_CLprint ('CCCCC_01 of tl_outerloop ')

                        call tl_lateral_inflow()
                        call tl_smallestBC_timeInterval ()
                    end if

                        ! call util_utest_CLprint ('DDDDD of tl_outerloop')

                    !% --- perform control rules
                    if ((.not. inSpinUpYN) .and. (setting%SWMMinput%N_control > 0)) then
                        if (setting%Time%Now .ge. setting%Time%ControlRule%NextTime) then
                            !% --- evaluate all the controls
                            !%     required for all images because monitorI data are spread
                            !%     across all images.
                            call control_update (oneI)
                            !% --- set the next time the controls will be evaluated
                            setting%Time%ControlRule%NextTime = setting%Time%Now  &
                                + real(setting%SWMMinput%ControlRuleStep,8)
                        end if
                    end if

                        ! call util_utest_CLprint ('EEEEE of tl_outerloop')

                    !% --- add subcatchment inflows
                    !%     note, this has "useHydrology" and not "doHydrologyStepYN" because the
                    !%     former is whether or not hydrology is used, and the latter is
                    !%     whether or not it is computed in this time step
                    if (setting%Simulation%useHydrology .and. BCupdateYN) then 
                        call tl_subcatchment_lateral_inflow () 
                    end if

                        ! call util_utest_CLprint ('FFFFF of tl_outerloop')

                    !% --- add RDII inflows
                    if (.not. setting%SWMMinput%IgnoreRDII) then 
                        call interface_get_RDII_inflow ()
                    end if
    
                    if ((.not. setting%SWMMinput%IgnoreGroundwater) .and. &
                              (setting%SWMMinput%N_groundwater > 0) ) then 
                        call interface_get_groundwater_inflow ()
                    end if 
    
                        ! call util_utest_CLprint ('GGGGG of tl_outerloop')

                    !% --- perform hydraulic routing
                    !  print *, '3A nextHydrologyTime ', setting%Time%Hydrology%NextTime
                    !  print *, 'calling tl_hydraulics'

                    
                    !% --- set hydraulics time step to handle inflow
                    call tl_update_hydraulics_timestep(.false.)


                    call tl_hydraulics()
                
                    !print *, 'out of tl_hydraulics'

                    ! print *, '4 nextHydrologyTime ', setting%Time%Hydrology%NextTime
                    ! print *, 'doHydraulicsStepYN',doHydraulicsStepYN
                    ! if (.not.doHydraulicsStepYN) then 
                    !     stop 4
                    ! end if

                        ! ! ! call util_utest_CLprint ('HHHHH in time_loop after tl_hydraulics')

                    !% --- accumulate RunOn from hydraulic elements to subcatchments
                    if ((setting%Simulation%useHydrology) .and. (any(subcatchYN(:,sYN_hasRunOn)))) then 
                        call tl_subcatchment_accumulate_runon ()
                    end if

                    ! print *, '5 nextHydrologyTime ', setting%Time%Hydrology%NextTime
                    ! print *, 'doHydraulicsStepYN',doHydraulicsStepYN
                    ! if (.not.doHydraulicsStepYN) then 
                    !     stop 5
                    ! end if

                        ! ! ! call util_utest_CLprint ('IIIII in time_loop after subcatchment_accumulate_runon')
    
                    !% --- close the clock tick for hydraulic loop evaluation
                    if ((this_image()==1) .and. (.not. inSpinUpYN)) then
                        call system_clock(count=cval,count_rate=crate,count_max=cmax)
                        setting%Time%WallClock%HydraulicsStop = cval
                        setting%Time%WallClock%HydraulicsCumulative &
                            = setting%Time%WallClock%HydraulicsCumulative &
                            + setting%Time%WallClock%HydraulicsStop &
                            - setting%Time%WallClock%HydraulicsStart
                    end if 

                    ! print *, '6 nextHydrologyTime ', setting%Time%Hydrology%NextTime
                    ! print *, 'doHydraulicsStepYN',doHydraulicsStepYN
                    ! if (.not.doHydraulicsStepYN) then 
                    !     stop 6
                    ! end if
    
                end if         
    
                ! print *, '7 nextHydrologyTime ', setting%Time%Hydrology%NextTime
                ! print *, 'doHydraulicsStepYN',doHydraulicsStepYN
                ! if (.not.doHydraulicsStepYN) then 
                !     stop 7
                ! end if

                !% --- handle output reporting
                if (setting%Output%Report%provideYN) then 
                    !% --- only provide output for spinup time if stopping after spinup
                     if ( (.not. inSpinUpYN)  &
                         .or. ((inSpinUpYN) .and. (SpinUpOutputYN)) ) then
    
                        !% set a time tick for output timing
                        if ((this_image()==1) .and. (.not. inSpinUpYN)) then
                            call system_clock(count=cval,count_rate=crate,count_max=cmax)
                            setting%Time%WallClock%LoopOutputStart = cval
                        end if
    
                        !% Results must be reported before the "do"counter increments
                        !call util_output_report()  --- this is a stub for future use
    
                        !% Multilevel time step output
                        if ((util_output_must_report()) .and. &
                            (.not. setting%Output%Report%suppress_MultiLevel_Output) ) then
                            call outputML_store_data (.false.)
                        end if
                        sync all
    
                        !% close the time tick for output timing
                        if ((this_image()==1) .and. (.not. inSpinUpYN)) then
                            call system_clock(count=cval,count_rate=crate,count_max=cmax)
                            setting%Time%WallClock%LoopOutputStop = cval
                            setting%Time%WallClock%LoopOutputCumulative &
                                = setting%Time%WallClock%LoopOutputCumulative &
                                + setting%Time%WallClock%LoopOutputStop  &
                                - setting%Time%WallClock%LoopOutputStart
                        end if
                    end if  
                end if

                ! print *, '8 nextHydrologyTime ', setting%Time%Hydrology%NextTime
                ! print *, 'doHydraulicsStepYN',doHydraulicsStepYN
                ! if (.not.doHydraulicsStepYN) then 
                !     stop 8
                ! end if
    
                call util_crashstop(13978)
    
                !% --- restart the hydraulics time tick
                sync all
                if ((this_image()==1) .and. (.not. inSpinUpYN) ) then
                    call system_clock(count=cval,count_rate=crate,count_max=cmax)
                    setting%Time%WallClock%HydraulicsStart = cval
                end if 

                ! print *, '9 nextHydrologyTime ', setting%Time%Hydrology%NextTime
                ! print *, 'doHydraulicsStepYN',doHydraulicsStepYN
                ! if (.not.doHydraulicsStepYN) then 
                !     stop 9
                ! end if
    
                !! ! call util_utest_CLprint('before time step change')
        
                sync all
                !% ---increment the time step and counters for the next time loop
                call tl_increment_timestep_and_counters(doHydraulicsStepYN, doHydrologyStepYN, inSpinUpYN)

                ! print *, '10 nextHydrologyTime ', setting%Time%Hydrology%NextTime
                ! print *, 'doHydraulicsStepYN',doHydraulicsStepYN
                ! if (.not.doHydraulicsStepYN) then 
                !     stop 10
                ! end if
    
                !% --- close the hydraulics time tick
                sync all
                if ((this_image()==1) .and. (.not. inSpinUpYN)) then
                    call system_clock(count=cval,count_rate=crate,count_max=cmax)
                    setting%Time%WallClock%HydraulicsStop = cval
                    setting%Time%WallClock%HydraulicsCumulative &
                        = setting%Time%WallClock%HydraulicsCumulative &
                        + setting%Time%WallClock%HydraulicsStop &
                        - setting%Time%WallClock%HydraulicsStart
                end if

                
                ! print *, 'doHydraulicsStepYN',doHydraulicsStepYN
                ! if (.not.doHydraulicsStepYN) then 
                !     stop 11
                ! end if

                !% --- check for blowup conditions
                ! call util_crashcheck (773623)
                if (crashI == 1) exit 

                ! print *, 'bottom of tl_outerloop'
    
            end do  !% end of time loop

    end subroutine tl_outerloop
!% 
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_hydrology()
        !%------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step
        !%------------------------------------------------------------------
        !% Declarations
            integer :: ii
            real(8), pointer :: sRunoff(:), dt
            integer(kind=8) :: crate, cmax, cval
            character(64) :: subroutine_name = 'tl_hydrology'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            !% wall clock tick
            if (this_image()==1) then
                call system_clock(count=cval,count_rate=crate,count_max=cmax)
                setting%Time%WallClock%HydrologyStart = cval
            end if 
        !%------------------------------------------------------------------              
        !% Aliases
            sRunoff => subcatchR(:,sr_RunoffRate_baseline)  
            dt      => setting%Time%Hydrology%Dt
        !%------------------------------------------------------------------   
            
        if (setting%Simulation%useSpinUp) then
            print *, 'ERROR setting.simulation.useSpinUp = true is not supported for hydrology at this time'
            call util_crashpoint(4482333)
        end if

        !% --- set the runon flowrate (if any)
        if (any(subcatchYN(:,sYN_hasRunOn))) then
            call tl_subcatchment_set_runon_volume ()
        end if
        
        ! print *, 'NewRunoffTime ',setting%Time%Hydrology%NextTime
        ! print *, 'Total Duration',setting%SWMMinput%TotalDuration
        
        !% --- only execute runoff if the next time is sufficently small
        if (setting%Time%Hydrology%NextTime < (setting%SWMMinput%TotalDuration - setting%Time%DtTol)) then
            
            ! print *, 'calling hydrology in tl_hydrology **********************************************'

            !% --- execute the EPA SWMM runoff for the next interval 
            call interface_call_runoff_execute()

            !% --- get the next runoff time
            setting%Time%Hydrology%NextTime = interface_get_NewRunoffTime()  
            ! print *, 'setting next hydrology time in tl_hydrology = ',setting%Time%Hydrology%NextTime

            !% --- cycle through the subcatchments to get the runoff 
            !%     ii-1 required in arg as C arrays start from 0
           ! print *, 'in tl_hydrology getting runoff'
            do ii = 1,setting%SWMMinput%N_subcatch
                sRunoff(ii) = interface_get_subcatch_runoff(ii-1)  
                !print *, 'ii,runoff ',ii,sRunoff(ii)
            end do
        else 
            sRunoff(:) = zeroR
        end if

        !%------------------------------------------------------------------   
        !% --- wall clock tick
        if (this_image()==1) then
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%HydrologyStop = cval
            setting%Time%WallClock%HydrologyCumulative  &
                = setting%Time%WallClock%HydrologyCumulative &
                + setting%Time%WallClock%HydrologyStop &
                - setting%Time%WallClock%HydrologyStart
        end if 
        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine tl_hydrology
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_hydraulics()
        !%------------------------------------------------------------------
        !% Description:
        !% Top level hydraulic solver for a single time step
        !%-------------------------------------------------------------------
        !% Declarations:
            integer, allocatable :: tempP(:)
            integer :: ii, kk
            real(8) :: tvolume
            character(64)    :: subroutine_name = 'tl_hydraulics'
        !%-------------------------------------------------------------------
        !% Preliminaries:
            !if (crashYN) return
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% --- repack all the dynamic arrays
        call pack_dynamic_arrays()

        !% --- ensure that the conservative flux terms are exactly zero in the entire array
        !%     so that we can be confident of conservation computation. 
        faceR(:,fr_Flowrate_Conservative) = zeroR  

        !% call the RK2 time-march
        call rk2_toplevel_ETM()

        call util_accumulate_volume_conservation () 

        ! tvolume = sum(elemR(1:N_elem(1),er_Volume))
        ! print *, 'tvolume in tl_hydraulics ',tvolume - faceR(1,fr_Flowrate_Conservative) * setting%Time%Hydraulics%Dt

        !%-------------------------------------------------------------------
        !% closing
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine tl_hydraulics
!% 
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_lateral_inflow()
        !%------------------------------------------------------------------
        !% Description:
        !% Handles lateral inflow update during time loop
        !% Note this sets lateral inflow to zero so it must be done before
        !% lateral inflows from subcatchments are added.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: npack, thisP(:), thisBC(:), nBarrels(:)
            real(8), pointer :: Qlateral(:), SeepRate(:), BreadthMax(:)
            real(8), pointer :: TopWidth(:), Area(:), AreaBelowBreadthMax(:)
            real(8), pointer :: Length(:), Fevap(:), EvapRate
            real(8) :: thisEpochTime, confac
            integer :: year, month, day
        !%------------------------------------------------------------------
        !% Aliases:
            Qlateral            => elemR(:,er_FlowrateLateral)
            SeepRate            => elemR(:,er_SeepRate)
            BreadthMax          => elemR(:,er_BreadthMax)
            TopWidth            => elemR(:,er_TopWidth)
            Area                => elemR(:,er_Area)
            AreaBelowBreadthMax => elemR(:,er_AreaBelowBreadthMax)
            Length              => elemR(:,er_Length)
            nBarrels            => elemI(:,ei_barrels)
            Fevap               => elemSR(:,esr_Storage_FractionEvap)
            EvapRate            => setting%Climate%EvapRate
            
            thisBC   => BC%P%BClat
        !%------------------------------------------------------------------
        !% set lateral to zero for all cells
        Qlateral(:) = zeroR 

        !% --- if ep_BClat exist add lateral inflow BC to lateral inflow accumulator
        !%     note that thisP and thisBC must be the same size or there is something wrong  
        !%     For multi-barrel elements, divide lateral inflow evenly between barrels
        npack    => npack_elemP(ep_BClat)
        thisP    => elemP(1:npack,ep_BClat)
        if (npack > 0) then    
            Qlateral(thisP) = Qlateral(thisP) &
                + BC%flowR(thisBC,br_value) / real(nBarrels(thisP),8)
                ! print *, 'tl_lateral_inflow'
                ! print *, BC%flowR(thisBC,br_value)
        end if

        ! print *, 'Qlateral(thisP)',Qlateral(thisP)

        !% --- find the Adjust.Conductivity multiplier for the current month
        thisEpochTime = util_datetime_secs_to_epoch(setting%Time%Now)
        call util_datetime_decodedate(thisEpochTime, year, month, day)
        confac = setting%Adjust%Conductivity(month)


        !% --- Subtract the seepage rate for all conduit and channel cells
        !%     Convert m/s seep to m^3/s flux by multiplying by topwidth and length
        !%     Note that Topwidth is either the actual topwidth or the maximum value if the present
        !%     cross-sectional area is greater than the area below the depth of maximum breadth
        !%     (i.e., we are using the maximum hydraulic cross section)
        !%     Note that in comparison between SWMM5+ and EPA-SWMM the average flowrate
        !%     in the SWMM5+ set of elements will be higher than the EPA-SWMM value for
        !%     a single link because the latter is giving the flow rate at the end of the
        !%     link after all the seepage is removed (i.e, what is expected in the last elem
        !%     of the link.)
        npack => npack_elemP  (ep_CC_NOTzerodepth)
        thisP => elemP(1:npack,ep_CC_NOTzerodepth)
        where (Area(thisP) > AreaBelowBreadthMax(thisP))
            Qlateral(thisP) = Qlateral(thisP) - confac * SeepRate(thisP) * BreadthMax(thisP) * Length(thisP)
        elsewhere
            Qlateral(thisP) = Qlateral(thisP) - confac * SeepRate(thisP) * TopWidth(thisP)   * Length(thisP)
        endwhere

        !% --- subtract the evaporation rate 
        if (setting%Climate%useHydraulicsEvaporationTF) then 

            !% --- apply to open channels
            npack    => npack_elemP  (ep_CC_Open_Elements)
            thisP    => elemP(1:npack,ep_CC_Open_Elements)
            !% --- apply evaporation to open-channel CC elements that are larger than zero depth
            where (elemR(thisP,er_Depth) > setting%ZeroValue%Depth)
                Qlateral(thisP) = Qlateral(thisP) - EvapRate * BreadthMax(thisP) * Length(thisP)
            endwhere

            !% --- apply evaporation to storage elements
            !%     note that Implied Storage has Fevap = 0.0, so this has no effect
            npack    => npack_elemP  (ep_JM)
            thisP    => elemP(1:npack,ep_JM)
            where (elemR(thisP,er_Depth) > setting%ZeroValue%Depth)
                Qlateral(thisP) = Qlateral(thisP) - EvapRate * Fevap(thisP) * BreadthMax(thisP) * Length(thisP)
            endwhere

            ! BRHWORKING NEED EXFILTRATION FOR STORAGE 
            ! routing.c  calls node_getLosses
            ! which calls node.c/node_getLosses
            ! which calls node.c/storage_getLosses
            ! which compute exfilRate as a local variable -- which is what we need

            ! IF EVERYTHING WE NEED IS STORED, THEN WE SHOULD BE ABLE TO 
            ! CALL exfil_getLoss() from an api.
            ! Must cycle through all the storage nodes.


        else
            !% --- evaporation has no effect
        end if




        ! print *, 'confac ',confac
        ! print *, 'seepRate ',SeepRate(thisP)
        ! print *, 'topwidth ',Topwidth(thisP)
        ! print *, 'BreadthMax ',BreadthMax(thisP)
        ! print *, 'Qlateral(thisP)',Qlateral(thisP)

        ! npack    => npack_elemP(ep_BClat)
        ! thisP    => elemP(1:npack,ep_BClat)
        ! print *, 'Qlateral(thisP)',Qlateral(thisP)
        ! stop 5590873

        !% --- HACK TODO: need seepage for storage elements 20220907 brh

    end subroutine tl_lateral_inflow
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_subcatchment_lateral_inflow ()
        !%------------------------------------------------------------------
        !% Description:
        !% gets the subcatchment inflows and adds to the hydraulics lateral
        !% inflow. Must be done AFTER the hydraulics lateral inflows are set
        !% Note: divide the inflow over the barrels in a multi-barrel element
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: sImage(:), eIdx(:), nBarrels(:)
            real(8), pointer :: Qrate(:), Qlateral(:)
            integer :: mm
        !%------------------------------------------------------------------
        !% Aliases
            sImage    => subcatchI(:,si_runoff_P_image)
            eIdx      => subcatchI(:,si_runoff_elemIdx)
            !% --- using the full runoff rate for this period
            Qrate    => subcatchR(:,sr_RunoffRate_baseline)
            Qlateral => elemR(:,er_FlowrateLateral)
            nBarrels => elemI(:,ei_barrels)
        !%------------------------------------------------------------------
        ! print *, 'in tl_subcatchment_lateral_inflow'
        ! print *, 'subcatch elem '
        ! print *, subcatchI(:,si_runoff_elemIdx)
            do mm = 1,setting%SWMMinput%N_subcatch
            !% --- only if this image holds this node
            
            ! print *, mm, eIdx(mm), Qlateral(eIdx(mm)), Qrate(mm)
            
            if (this_image() .eq. sImage(mm)) then
                Qlateral(eIdx(mm)) = Qlateral(eIdx(mm)) &
                     + Qrate(mm) / real(nBarrels(mm),8)
            end if
        end do



    end subroutine tl_subcatchment_lateral_inflow
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_subcatchment_accumulate_runon ()
        !%------------------------------------------------------------------
        !% Description
        !% accumulates the runon volume from hydraulic element to a 
        !% subcatchment
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: mm, kk
            integer, pointer :: nRunOn, fIdx
            real(8), pointer :: dt
        !%------------------------------------------------------------------
        !% Aliases
            dt       => setting%Time%Hydraulics%Dt
        !%------------------------------------------------------------------

        !% --- cycle through the subcatchments
        do mm = 1,setting%SWMMinput%N_subcatch
            !% --- get the count of runons
            nRunOn => subcatchI(mm,si_RunOn_count)
            if (nRunOn == 0) cycle

            !% --- cycle through the runons to this subcatchment
            do kk=1,nRunOn
                !% --- get runons for the image associated with the outfall node
                if (subcatchI(mm,si_RunOn_P_image(kk)) .eq. this_image()) then 
                    !% --- identify the face
                    fIdx => subcatchI(mm,si_RunOn_faceIdx(kk))
                    !% --- store the volume associated with the face flowrate
                    !%     this uses the conservative flowrate
                    subcatchR(mm,sr_RunOn_Volume(kk)) &
                        = subcatchR(mm,sr_RunOn_Volume(kk))               &
                          + faceR(fIdx,fr_Flowrate_Conservative) * dt
                else 
                    subcatchR(mm,sr_RunOn_Volume(kk)) = zeroR
                end if
            end do
        end do

        
    end subroutine tl_subcatchment_accumulate_runon  
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_subcatchment_set_runon_volume()
        !%------------------------------------------------------------------
        !% Description: 
        !% computes the runon volumes for all outfalls that end on
        !% subcatchments
        !%------------------------------------------------------------------
        !% Declarations
            integer :: mm, kk
            integer, pointer :: nRunOn 
            real(8), pointer :: dt
        !%------------------------------------------------------------------
        !% Aliases:
            dt   => setting%Time%Hydrology%Dt
        !%------------------------------------------------------------------
        
        !% --- cycle through the subcatchments
        do mm=1,setting%SWMMinput%N_subcatch
            nRunOn => subcatchI(mm,si_RunOn_count)
            if (nRunOn == 0) cycle
            
            !print *, 'dt hydrology in tl_subcatchment_set_runon_volume',dt

            do kk = 1,nRunOn               
                if (subcatchI(mm,si_RunOn_P_image(kk)) == this_image()) then
                    !% --- change all values to a flowrate
                    subcatchR(mm,sr_RunOn_Volume(kk)) = subcatchR(mm,sr_RunOn_Volume(kk))

                    ! print *, 'runon volume',subcatchR(mm,sr_RunOn_Volume(kk))
                else 
                    subcatchR(mm,sr_RunOn_Volume(kk)) = zeroR
                end if

                !% --- disgribute the runon volume to all images
                call co_max(subcatchR(mm,sr_RunOn_Volume(kk)))
            end do

            !% --- export the volume to EPASWMM outfall.vRouted storage
            call interface_export_runon_volume ()

            !% --- reset the runon values to zero
            subcatchR(mm,sr_RunOn_Volume(:)) = zeroR

        end do

    end subroutine tl_subcatchment_set_runon_volume
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_increment_timestep_and_counters ( &
        doHydraulicsStepYN, doHydrologyStepYN, inSpinUpYN)
        !%-------------------------------------------------------------------
        !% Description:
        !% increments the hydrology and hydraulics step counters and
        !%-------------------------------------------------------------------
        !% Declarations
            logical, intent(inout) :: doHydraulicsStepYN, doHydrologyStepYN
            logical, intent(in)    :: inSpinUpYN
            logical, pointer       :: useHydrology, useHydraulics
            real(8), pointer       :: nextHydraulicsTime, nextHydrologyTime
            real(8), pointer       :: lastHydraulicsTime, lastHydrologyTime
            real(8)                :: nextTime, nextTimeLocal
            real(8), pointer       :: dtTol, timeNow, dtHydraulics, dtHydrology
            integer                :: minImg
            integer, pointer       :: reportStep
            integer(kind=8), pointer    :: hydraulicStep, hydrologyStep, step
            character(64)          :: subroutine_name = 'tl_increment_counters'
        !%-------------------------------------------------------------------
        !% Preliminaries
            !if (crashYN) return
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%--------------------------------------------------------------------
        !% Aliases
            hydraulicStep => setting%Time%Hydraulics%Step
            hydrologyStep => setting%Time%Hydrology%Step
            dtHydraulics  => setting%Time%Hydraulics%Dt
            dtHydrology   => setting%Time%Hydrology%Dt
            timeNow       => setting%Time%Now
            step          => setting%Time%Step
            reportStep    => setting%Output%Report%ThisStep
            dtTol         => setting%Time%DtTol

            nextHydraulicsTime => setting%Time%Hydraulics%NextTime
            nextHydrologyTime  => setting%Time%Hydrology%NextTime

            lastHydraulicsTime => setting%Time%Hydraulics%LastTime
            lastHydrologyTime  => setting%Time%Hydrology%LastTime

            useHydrology  => setting%Simulation%useHydrology
            useHydraulics => setting%Simulation%useHydraulics
        !%--------------------------------------------------------------------

        !% --- get the timestep and the next time for hydraulics
        if (doHydraulicsStepYN) then
            call tl_update_hydraulics_timestep(inSpinUpYN)
            call util_crashstop(449873)
        else
            nextHydraulicsTime = setting%Time%End + tenR*DtTol
        end if    

        !print *, 'DT   AAAA       ',dtHydraulics 
        ! print *, 'next Hydraulics Time AAA ',nextHydraulicsTime

        !% --- The NextHydrologyTime is updated in SWMM, here we just need to
        !%     provide a large number if hydrology isn't used
        if (.not. useHydrology) then
            nextHydrologyTime = max(setting%Time%End + tenR*DtTol, &
                                    setting%Simulation%SpinUpDays * seconds_per_day)
           ! print *, 'setting NextHydrologyTime without hydrology ',nextHydrologyTime

        else
            !% --- SWMM-C will not return a next hydrology time that is
            !%     beyond the end of the simulation, so we need a special
            !%     treatment for this case
            if ((nextHydrologyTime .eq. lastHydrologyTime) .and. &
                (timeNow > setting%Time%Start)) then 
                nextHydrologyTime = max(setting%Time%End + tenR*DtTol , &
                                        setting%Simulation%SpinUpDays * seconds_per_day)

                !print *, 'setting last hydrology time = ',nextHydrologyTime
            end if
        end if
        ! print *, 'next Hydrology  Time BBB ',nextHydrologyTime
        ! print *, 'next Hydraulics Time BBB',nextHydraulicsTime

        !% --- Ensure that all processors use the same time step.
        !% --- find the minimum hydraulics time and store accross all processors
        call co_min(nextHydraulicsTime)
        !% --- take the nextTime as the minimum of either the Hydrology or Hydraulics time
        !%     done on a single processor because they all should have identical nextHydrologyTIme
        nextTime = min(nextHydraulicsTime, nextHydrologyTime)

        ! print *, 'next Time CCC ',nextTime

        !% --- note there is no need to broadcast across processors since co_min
        !%     ensures that each processor gets the min value.

        !% --- update the time step for the local processor (all have the same nextTime and lastHydraulicsTime)
        dtHydraulics = nextTime - lastHydraulicsTime

        ! print *, 'dtHydarulics, next time ',dtHydraulics, nextTime

        ! if (dtHydraulics == zeroR) stop 4498723

        doHydrologyStepYN  = (abs(nextTime - nextHydrologyTime)  <= dtTol) .and. useHydrology
        doHydraulicsStepYN = (abs(nextTime - nextHydraulicsTime) <= dtTol) .and. useHydraulics

        if (doHydrologyStepYN) then 
            dtHydrology = nextTime - lastHydrologyTime
        end if

        if ((.not. inSpinUpYN) .or. & 
            ((inSpinUpYN) .and. (setting%Simulation%stopAfterSpinUp)) ) then
            if (util_output_must_report()) reportStep = reportStep + 1
        end if
        if (doHydraulicsStepYN) hydraulicStep = hydraulicStep + 1
        if (doHydrologyStepYN)  hydrologyStep = hydrologyStep + 1

        step    = step + 1
        timeNow = nextTime !timeNow + dt

        ! print *, 'timeNow  in tl_increment_counters',timeNow

        call tl_command_line_step_output(inSpinUpYN)

        if (doHydraulicsStepYN) LastHydraulicsTime = NextHydraulicsTime
        if (doHydrologyStepYN)  LastHydrologyTime  = NextHydrologyTime

        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        end subroutine tl_increment_timestep_and_counters
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_update_hydraulics_timestep(inSpinUpYN)
        !%------------------------------------------------------------------
        !% Description:
        !% updates the timestep (dt) for hydraulics and computes the
        !% setting.Time.Hydraulics.NextTime for the current processor
        !%------------------------------------------------------------------
        !% Declarations
            logical, intent(in) :: inSpinUpYN
            logical, pointer :: matchHydrologyStep, useHydrology
            logical          :: do_roundoff
            real(8)          :: oldDT, reportDt, maxVelocity
            real(8)          :: timeleft, thisCFL, minCFL
            real(8), pointer :: targetCFL, maxCFL, maxCFLlow, timeNow, dtTol
            real(8), pointer :: increaseFactor
            real(8), pointer :: newDT
            !rm velocity(:), wavespeed(:), length(:), PCelerity(:)
            real(8), pointer :: nextHydrologyTime, nextHydraulicsTime
            real(8), pointer :: lastHydrologyTime, lastHydraulicsTime
            integer          :: ii, neededSteps, pindex(1)
            integer, pointer ::   checkStepInterval
            !rm integer, pointer :: thisCol, Npack, thisP(:)
            integer(kind=8), pointer :: stepNow, lastCheckStep
            character(64)    :: subroutine_name = 'tl_update_hydraulics_timestep'
            !integer :: kk !temporary
        !%-------------------------------------------------------------------
        !% Preliminaries
            !if (crashYN) return
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases:
            maxCFL             => setting%VariableDT%CFL_hi_max
            targetCFL          => setting%VariableDT%CFL_target
            maxCFLlow          => setting%VariableDT%CFL_lo_max
            increaseFactor     => setting%VariableDT%increaseFactor
            checkStepInterval  => setting%VariableDT%NstepsForCheck

            useHydrology       => setting%Simulation%useHydrology

            matchHydrologyStep => setting%Time%matchHydrologyStep

            nextHydrologyTime  => setting%Time%Hydrology%NextTime
            nextHydraulicsTime => setting%Time%Hydraulics%NextTime

            lastHydrologyTime  => setting%Time%Hydrology%LastTime
            lastHydraulicsTime => setting%Time%Hydraulics%LastTime        

            timeNow            => setting%Time%Now
            dtTol              => setting%Time%DtTol

            stepNow            => setting%Time%Hydraulics%Step
            lastCheckStep      => setting%VariableDT%LastCheckStep
        !%----------------------------------------------------------------------
        oldDT    =  setting%Time%Hydraulics%Dt  ! not a pointer (important!)
        newDT    => setting%Time%Hydraulics%Dt
        reportDt = setting%Output%Report%TimeInterval

            ! print *, '============================================================='
            ! print *, 'in ', trim(subroutine_name) 
            ! print *, 'starting dt', oldDT, setting%Time%Hydraulics%Dt
            ! print *, 'Hydrology DT',setting%Time%Hydrology%Dt

        !% --- set the minimum CFL, used to detect near zero flow conditions
        if (setting%Limiter%Dt%UseLimitMaxYN) then
            minCFL = setting%Eps%Velocity * oldDT / setting%Discretization%NominalElemLength
        end if

            ! print *, 'minCFL ',minCFL
            ! print *, ' '
            ! print *, 'in ',trim(subroutine_name)
            ! print *, matchHydrologyStep, useHydrology, inSpinUpYN
            ! print *, ' '

        if ((matchHydrologyStep) .and. (useHydrology) .and. (.not. inSpinUpYN)) then 

            print *, 'in matching hydrology step'

            !% --- for combined hydrology and hydraulics compute the CFL if we take a single
            !%     step to match the next hydrology time
            timeLeft = nextHydrologyTime - lastHydraulicsTime

            ! print *, 'time left in Hydrology timestep',timeleft

            if (timeLeft .le. dtTol) timeLeft = oldDT

            !% --- get the CFL if a single step is taken using only CC elements
            !if (setting%SmallDepth%useMomentumCutoffYN) then
            !    thisCFL = tl_get_max_CFL_CC(ep_CC_NOTsmalldepth,timeleft)  
            !else
                thisCFL = tl_get_max_CFL_CC(ep_CCJM_NOTzerodepth,timeleft)
            !end if

                ! print *, 'thisCFL to reach hydrology step ',thisCFL

            !% --- check to see if a single time step to match the hydrology time is possible
            if (thisCFL .le. maxCFL) then
                neededSteps = 1 !% only 1 hydraulic step to match hydrology step

                
                    ! print *, 'one time step timeleft, oldDT ',timeLeft, oldDT

                if (timeLeft > dtTol) then
                    !% --- check to be sure there is significant time left
                    ! print *, 'timeleft /oldDT, increasefactor ',timeLeft /oldDT , increaseFactor

                    !% --- check increase that is implied with a single time step to the next hydrology
                    if (timeleft / oldDT < increaseFactor) then
                        !% --- use a single hydraulics step for the remaining hydrology time
                        newDT = timeleft
                            ! print *, 'time left, newDT ',timeLeft, newDT
                    else
                        !% --- increase the oldDT by the maximum allowed
                        newDT = increaseFactor * oldDT
                            ! print *, 'newDT',newDT
                        neededSteps = ceiling(timeLeft/newDT)
                            ! print *, 'needed steps HERE ',neededSteps
                    end if
                else
                    !% --- not significant time left

                    ! print *, 'in ELSE section AAA timeleft ',timeleft 
                        
                    !% --- negligible time left, don't bother with it, go back to the old dt
                    newDT = oldDT
                    !% --- check that resetting to oldDT didn't cause a problem with CC elements
                    !if (setting%SmallDepth%useMomentumCutoffYN) then
                    !    thisCFL = tl_get_max_CFL_CC(ep_CCJM_NOTsmalldepth,newDT)
                    !else 
                        thisCFL = tl_get_max_CFL_CC(ep_CCJM_NOTzerodepth,newDT)
                    !end if
                    if (thisCFL > maxCFL) then 
                        !% --- if CFL too large, set the time step based on the target CFL
                        newDT = newDT * targetCFL / thisCFL
                        neededSteps = ceiling(timeLeft/newDT)
                    end if

                    ! print *, 'CFL, DT ',thisCFL, newDT

                end if
            else
                ! print *, 'in ELSE section BBB'
                ! print *, 'for required CFL > max CFL'

                !% --- if more than one time step is needed, compute the needed steps 
                !%     to break up the large CFL and time step size.
                if (thisCFL/targetCFL .ge. huge(neededSteps)) then
                    !% --- check to see if the implied time step is too small for the integer size
                    !%     this is likely a blow-up condition
                    newDT = setting%Limiter%Dt%Minimum + setting%Time%DtTol
                    neededSteps = 1000
                    write(*,*) 'warning -- really high velocity, setting dt to minimum to prevent overflow'
                    write(*,*), 'DT = ', newDT
                else
                    !% ---  for a moderate case where multiple hydraulic steps are needed for a hydrology step
                    if (thisCFL > minCFL) then
                        ! print *, 'thisCFL > minCFL'

                        !% --- note that neededSteps will be 2 or larger else thisCFL < maxCFL
                        neededSteps = ceiling( thisCFL / targetCFL )
                        !% --- the provisional time step that would get exactly to the hydrology time (if CFL didn't change)
                        newDT = timeleft / real(neededSteps,8)
                        
                        ! print *, 'timeleft, steps ',timeleft, neededSteps
                        ! print *, 'newDT', newDT
                        ! print *, 'increaseFactor ',increaseFactor

                        !% --- limit the change in the dt by the increase factor
                        if (newDT / oldDT > increaseFactor) then 
                            newDT = oldDT * increaseFactor
                        else
                            !% --- accept the provisional newDT
                        end if

                        ! print *, 'newDT after limiting by increaseFactor',newDT

                        !% --- check that the newDT didn't cause a CFL violation
                        ! if (setting%SmallDepth%useMomentumCutoffYN) then
                        !     thisCFL = tl_get_max_CFL_CC(ep_CC_NOTsmalldepth,newDT)
                        ! else
                            thisCFL = tl_get_max_CFL_CC(ep_CCJM_NOTzerodepth,newDT)
                        ! end if
                        if (thisCFL > maxCFL) then 
                            !% --- if CFL to large, set the time step based on the target CFL
                            newDT = newDT * targetCFL / thisCFL
                        end if

                        ! print *, 'newDT after CC check', newDT
                    else
                        !% --- Miniscule CFL is typically due to small volumes in the start-up condition
                        !%    we can use a larger time step based on volumes to fill
                        call tl_dt_vanishingCFL(newDT)

                        ! print *, 'Vanishing dt ',newDT
                    end if

                end if
            end if

            ! print *, 'newDT based on Hydrology step', newDT

        else
 
                ! print *, 'CFL with no hydrology '

            if (useHydrology) then 
                print *, 'USER CONFIGURATION ERROR: At this time, useHydrology requires '
                print *, 'setting%Time%matchHydrologyStep = .true.'
                print *, 'the setup for non-match has not been tested'
                call util_crashpoint(509872)
            end if

            !% --- for hydraulics without hydrology or if we are not matching the hydrology step
            neededSteps = 3 !% forces rounding check

            !% --- allowing hydrology and hydraulics to occur at different times
            ! if (setting%SmallDepth%useMomentumCutoffYN) then
            !     thisCFL = tl_get_max_CFL_CC(ep_CC_NOTsmalldepth,oldDT)
            ! else 
                thisCFL = tl_get_max_CFL_CC(ep_CCJM_NOTzerodepth,oldDT)
            ! end if
        
                ! print *, ' '
                ! print *, 'baseline CFL, minCFL, this step: '
                ! print *, thisCFL, minCFL, stepNow

            if (thisCFL .ge. maxCFL) then
                !% --- always decrease DT for exceeding the max CFL
                !%     new value is based on the target CFL
                newDT = oldDT * targetCFL / thisCFL
                lastCheckStep = stepNow

                    ! print *, 'Adjust DT for high CFL    ',newDT   

            else 
                !% --- if CFL is less than max, see if it can be raised (only checked at intervals)
                if (stepNow >= lastCheckStep + checkStepInterval) then
                !   print *, '------------------------------------------------------------------------------'
                !   print *, 'checking step: ',stepNow ,lastCheckStep, checkStepInterval

                    !% --- check for low CFL only on prescribed intervals and increase time step
                    if (thisCFL .le. minCFL) then
                        !% --- for really small CFL, the new DT could be unreasonably large (or infinite)
                        !%     so use a value based on inflows to fill to small volume
                            ! print *, 'vanishing CFL'
                        call tl_dt_vanishingCFL(newDT)

                    elseif ((minCFL < thisCFL) .and. (thisCFL .le. maxCFLlow)) then
                        !% --- increase the time step and reset the checkStep Counter
                            !    print *, 'lowCFL'
                        newDT = OldDT * targetCFL / thisCFL 
                        thisCFL = tl_get_max_CFL_CC(ep_CCJM_NOTzerodepth,newDT)
                    else
                        !% -- for maxCFLlow < thisCFL < maxCFL do nothing
                    end if
                    lastCheckStep = stepNow
                end if
            end if
        end if

        !% --- prevent large increases in the time step
        newDT = min(newDT,OldDT * increaseFactor)

        !% saz 20230328
        !% --- HACK: limit the newDT to be equal or smaller than the reportDT
        !% this ensures we get consistant output at every reporting steps
        !% however, this will slow down the timeloop if reportDT causes
        !% lower cfl values. Needs to revisit later. 
        newDT = min(newDT,reportDt)
            ! print *, 'newDT, oldDT :       ',newDT, OldDT

        !% 20220328brh time step limiter for inflows into small or zero volumes
        call tl_limit_BCinflow_dt (newDT)
            !  print *, 'dt BC flow limit     ',newDT

        call tl_limit_BChead_inflow_dt (newDT)
            !  print *, 'dt BC head limit     ',newDT

        call tl_limit_LatInflow_dt (newDT)
            !  print *, 'dt Qlat limit   ',newDt

        if ((matchHydrologyStep) .and. (useHydrology)) then
            neededSteps = ceiling(timeLeft/newDT)
        else 
            neededSteps = 1000
        end if


        !% --- limit by inflow/head external boundary conditions time intervals
        if (setting%VariableDT%limitByBC_YN) then
            ! print *, 'here limiting DT by BC ',newDT, setting%BC%smallestTimeInterval
            newDT = min(setting%BC%smallestTimeInterval,newDT)
        end if

        select case (neededSteps)
            case(1)
                !% -- accept dt as is
                do_roundoff = .false.
            case(2,3)
                do_roundoff = .true.
                !% -- set the smaller value to round off
                newDT = min(newDT, timeleft/real(neededSteps,8))
            case default
                do_roundoff = .true.
        end select

        if (do_roundoff) then 
            if (newDT > fiveR) then
                !% --- round down larger dt to counting numbers
                newDT = real(floor(newDT),8)
                    ! print *, 'rounding large    ',newDT
            elseif ((newDT > oneR) .and. (newDT .le. fiveR)) then
                !% --- round down smaller dt to two decimal places
                newDT = real(floor(newDT * tenR),8) / tenR
                    ! print *, 'rounding small    ',newDT   
            elseif ((newDT > onetenthR) .and. (newDT .le. oneR)) then
                !% --- round down smaller dt to three decimal places
                newDT = real(floor(newDT * onehundredR),8) / onehundredR
                    ! print *, 'rounding smaller  ',newDT
            else 
                !% do not round smaller values
            end if
        else 
            !% do not round
        end if


        ! ! ! print *, 'neededSteps ',neededSteps
        ! ! !% --- if dt is large and there is more than 2 steps, then round to an integer number
        ! if ((matchHydrologyStep) .and. (useHydrology) .and. (neededSteps .le. 2)) then
        !     !% don't round the dt
        ! else
        !     if (newDT > fiveR) then
        !         !% --- round larger dt to counting numbers
        !         newDT = real(floor(newDT),8)

        !             ! print *, 'rounding large    ',newDT
        !     elseif (newDT > oneR) then
        !         !% --- round smaller dt to two decimal places
        !         newDT = real(floor(newDT * onehundredR),8) / onehundredR

        !             ! print *, 'rounding small    ',newDT   
        !     else
        !         !%  HACK -- should round to 3 places for smaller numbers
        !     end if
        ! end if

            ! print *, 'final DT: ',newDT
            ! print *, ' '

        ! if (newDT < 0.01d0) then 
        !     print *, 'ERROR: time step has dropped below 0.01 s. need to investigate!'
        !     call util_crashpoint(119874)
        ! end if

        !% increment the hydraulics time clock
        nextHydraulicsTime = lastHydraulicsTime + newDT
        
        !% find the cfl for reporting
       ! cfl_max = tl_get_max_CFL_CC(ep_CCJBJM_NOTsmalldepth,newDT)
        ! if (setting%SmallDepth%useMomentumCutoffYN) then 
        !     cfl_max = tl_get_max_CFL_CC(ep_CCJM_NOTsmalldepth,newDT)
        ! else 
            cfl_max = tl_get_max_CFL_CC(ep_CCJM_NOTzerodepth,newDT)
        ! end if
        call co_max(cfl_max)


        !stop 29873

        !%----------------------------------------------------------------------
        !% closing
            if ((setting%Limiter%Dt%UseLimitMinYN) .and. (newDT .le. setting%Limiter%Dt%Minimum)) then
                print *, ' '
                print *, 'EXITING ON TIME STEP ERROR -- PROBABLY BLOWING UP DUE TO EXCESSIVE HEAD'
                print *,'timestep= ', setting%Time%Step
                print*, 'timeNow = ', timeNow, ' seconds'
                print*, 'dt = ', newDT, 'minDt = ',  setting%Limiter%Dt%Minimum
                ! if (setting%SmallDepth%useMomentumCutoffYN) then
                !     print*, 'max velocity  ', maxval( &
                !         elemR(elemP(1:npack_elemP(ep_CC_NOTsmalldepth),ep_CC_NOTsmalldepth),er_Velocity) )
                !     print*, 'max wavespeed ', maxval( &
                !         elemR(elemP(1:npack_elemP(ep_CC_NOTsmalldepth),ep_CC_NOTsmalldepth),er_WaveSpeed) )
                !     print*, 'warning: the dt value is smaller than the user supplied min dt value'
                ! else
                    print*, 'max velocity  ', maxval( &
                        elemR(elemP(1:npack_elemP(ep_CCJM_NOTzerodepth),ep_CCJM_NOTzerodepth),er_Velocity) )
                    print*, 'max wavespeed ', maxval( &
                        elemR(elemP(1:npack_elemP(ep_CCJM_NOTzerodepth),ep_CCJM_NOTzerodepth),er_WaveSpeed) )
                    print*, 'warning: the dt value is smaller than the user supplied min dt value'

                    print *, ' '
                    print *, 'element index location of max velocity'
                    pindex = maxloc(elemR(elemP(1:npack_elemP(ep_CCJM_NOTzerodepth),ep_CCJM_NOTzerodepth),er_Velocity))
                    print *, elemP(pindex,ep_CCJM_NOTzerodepth)
                    print *, 'element index location of max wavespeed'
                    pindex = maxloc(  elemR(elemP(1:npack_elemP(ep_CCJM_NOTzerodepth),ep_CCJM_NOTzerodepth),er_WaveSpeed))
                    print *, elemP(pindex,ep_CCJM_NOTzerodepth)
                ! end if
                !stop 1123938
                call util_crashpoint(1123938)
            end if

            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine tl_update_hydraulics_timestep
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_dt_vanishingCFL (dt)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes a time step when the maximum CFL is near zero. Here we look at the 
        !% maximum flow rates on any face or the element itself and the time it
        !% takes for the flow rate achieve a depth deeper than the small volume limit.
        !% This should only be used if EVERYWHERE has a very small CFL. We're
        !% looking for the inflow place that will control the time step
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), intent(inout) :: dt
            integer, pointer :: thisP(:), fDn(:), fUp(:), npack, ncol
            real(8), pointer :: Volume(:), SmallVolume(:), Area(:)
            real(8), pointer :: Qface(:), Qelem(:), Qlat(:)
            real(8), pointer :: flowrate(:), volumeDelta(:), timeToDepth(:)
            real(8) :: dt1, dt2
        !%------------------------------------------------------------------
        !% Aliases
            ncol        => col_elemP(ep_CC) !% cannot use JM elements
            npack       => npack_elemP(ncol)
            if (npack < 1) return
            thisP       => elemP(1:npack,ncol)
            flowrate    => elemR(1:npack,er_Temp01)
            volumeDelta => elemR(1:npack,er_Temp02)
            timeToDepth => elemR(1:npack,er_Temp03)
            fDn         => elemI(:,ei_Mface_dL)
            fUp         => elemI(:,ei_MFace_uL)
            Area        => elemR(:,er_Area)
            Volume      => elemR(:,er_Volume)
            SmallVolume => elemR(:,er_SmallVolume)
            Qface       => faceR(:,fr_Flowrate)
            Qelem       => elemR(:,er_Flowrate)
            Qlat        => elemR(:,er_FlowrateLateral)
        !%------------------------------------------------------------------
        !% --- initialize temporary arrays
        flowrate(:)    = nullvalueR
        volumeDelta(:) = nullvalueR
        timeToDepth(:) = nullvalueR

        !% --- maximum flowrate for any element
        flowrate(:)     = max( abs(Qface(fUp(thisP)) ),  &
                               abs(Qface(fDn(thisP)) ),  &
                               abs(Qelem    (thisP)  ),  &
                               abs(Qlat     (thisP)  ) ) 
                          
        !% --- set up minimum flowrate for dt purposes
        flowrate = max(flowrate, (Area(thisP) * setting%ZeroValue%Velocity))     
           
        !% --- volume delta to fill to small volume level
        volumeDelta(:) = SmallVolume(thisP) - Volume(thisP)
        
        !% --- replace negative and very small volume delta with zero volume value
        volumeDelta(:) = max(volumeDelta, setting%ZeroValue%Volume)
        

        !% --- time to reach 110% of small volume 
        timeToDepth = (volumeDelta + onetenthR * SmallVolume(thisP) )  / flowrate

        !% --- set timeToDepth for small volume delta to a large number
        !%     so that they will not play a role
        where (volumeDelta == setting%ZeroValue%Volume)
            timeToDepth = abs(nullValueR)
        endwhere

        !% --- set time step limiter as the minimum value 
        dt = minval(timeToDepth)

        !% --- use maximum limiter on the time step
        if (setting%Limiter%Dt%UseLimitMaxYN) then
            dt = min(dt,setting%Limiter%Dt%Maximum)
        end if

        !% --- check to see if result is null value
        !%     if so, set dt = 60.0  (HACK)
        if (dt == nullvalueR) then
            dt = sixtyR
        end if
    
        !% --- truncate extraneous subseconds
        if (dt > twoR) then
            dt = real(floor(dt),8)
        end if
        
        !print *, 'DT for new_zero subr ',dt

        !% --- reset the temporary arrays
        flowrate(:)    = nullvalueR
        volumeDelta(:) = nullvalueR
        timeToDepth(:) = nullvalueR

    end subroutine tl_dt_vanishingCFL
!%
!%==========================================================================
!%==========================================================================
!%       
!%
!%==========================================================================
!%==========================================================================
!%    
    ! subroutine tl_solver_select()
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% For ETM_AC dual method, this sets the elemI(:,ei_tmType) to the type of solver
    !     !% needed depending on the volume and the volume cutoffs.
    !     !% Should only be called if setting%Solver%SolverSelect == ETM_AC
    !     !%-----------------------------------------------------------------
    !     !% Delcarations
    !         integer :: thisCol
    !         integer, pointer :: Npack, tmType(:), thisP(:)
    !         real(8), pointer :: sfup, sfdn
    !         real(8), pointer :: volume(:), FullVolume(:)
    !         character(64) :: subroutine_name = 'tl_solver_select'
    !     !%-------------------------------------------------------------------
    !     !% Preliminaries
    !         if (setting%Solver%SolverSelect .ne. ETM_AC) return
    !         if (setting%Debug%File%timeloop) &
    !             write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     !%-------------------------------------------------------------------
    !     !% Aliases:       
    !         thiscol = ep_CCJM
    !         Npack => npack_elemP(thisCol)
    !         thisP => elemP(1:Npack,thisCol)

    !         tmType     => elemI(:,ei_tmType)
    !         volume     => elemR(:,er_Volume)
    !         FullVolume => elemR(:,er_FullVolume)

    !         sfup => setting%Solver%SwitchFractionUp
    !         sfdn => setting%Solver%SwitchFractionDn
    !     !%-------------------------------------------------------------------
    !     !% Look for ETM elements that are above the cutoff for going to AC and set
    !     ! !% these to AC
    !     ! where ( ( (volume(thisP) / FullVolume(thisP) ) > sfup ) .and. (tmType(thisP) == ETM) )
    !     !     tmType(thisP) = AC
    !     ! endwhere

    !     !% AC IS OBSOLETE

    !     !% Look for AC elements that are below the cutoff for going back to ETM and
    !     !% set these to ETM
    !     where ( ( (volume(thisP) / FullVolume(thisP) ) < sfdn) .and. (tmType(thisP) == AC) )
    !         tmType(thisP) = ETM
    !     endwhere

    !     !%-------------------------------------------------------------------
    !         if (setting%Debug%File%timeloop) &
    !             write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine tl_solver_select
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_save_previous_values()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Pushes the time N values into time N-1 storage, and the time N+1 values into
        !% the time N storage.
        !% HACK -- 20210809 brh This would be better done by changing the indexes without
        !% moving the data, but we will wait on doing this until we have the code
        !% debugged.
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'tl_save_previous_values'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%  push the old values down the stack
        !%  N values is the present, N0 is the last time step, and N1
        !%  is the timestep before (needed only for backwards 3rd in velocity and volume)
        elemR(:,er_Flowrate_N0)  = elemR(:,er_Flowrate)
        elemR(:,er_Head_N0)      = elemR(:,er_Head)
        elemR(:,er_Velocity_N1)  = elemR(:,er_Velocity_N0)
        elemR(:,er_Velocity_N0)  = elemR(:,er_Velocity)
        elemR(:,er_Volume_N1)    = elemR(:,er_Volume_N0)
        elemR(:,er_Volume_N0)    = elemR(:,er_Volume)
        elemR(:,er_Area_N1)      = elemR(:,er_Area_N0)
        elemR(:,er_Area_N0)      = elemR(:,er_Area)
        elemR(:,er_SlotVolume_N0) = elemR(:,er_SlotVolume)
        elemR(:,er_SlotDepth_N0)  = elemR(:,er_SlotDepth)

        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine tl_save_previous_values
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_command_line_step_output (inSpinUpYN)
        !%-----------------------------------------------------------------------------
        !% Description:
        !%
        !%-----------------------------------------------------------------------------
            logical, intent(in) :: inSpinUpYN
            character(64) :: subroutine_name = 'tl_command_line_step_output'
            real (8), pointer :: dt, timeNow, timeEnd
            real (8) :: thistime
            integer, pointer :: interval
            integer (kind=8), pointer :: step
            real(8) :: execution_realtime
            real(8) execution_realtime_per_step, steps_to_finish
            integer(kind=8) :: cval, crate, cmax
            real(8) :: simulation_fraction, seconds_to_completion, time_to_completion
            character(8) :: timeunit
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        dt            => setting%Time%Hydraulics%Dt
        timeNow       => setting%Time%Now
        timeEnd       => setting%Time%End
        step          => setting%Time%Step
        interval      => setting%Output%CommandLine%interval

        if (this_image() == 1) then
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%Now = cval

            ! --- estimate the remaining time based on total time from start
            execution_realtime = real((setting%Time%WallClock%Now - setting%Time%WallClock%TimeMarchStart),kind=8) &
                             /  (real(setting%Time%WallClock%CountRate,kind=8)) 
            execution_realtime_per_step = execution_realtime / real(1+step,kind=8)

            ! execution_realtime = real((setting%Time%WallClock%Now - setting%Time%WallClock%LastTimeStored),kind=8) &
            !                   /  (real(setting%Time%WallClock%CountRate,kind=8)) 

            ! seconds_to_completion =  (    (real(execution_realtime,kind=8))                     &
            !                            /  (real(setting%Time%WallClock%CountRate,kind=8))  )    &
            !                        * (    (setting%Time%End - setting%Time%Now)                 &
            !                            /  (setting%Time%Now - setting%Time%Start) )
            !execution_realtime_per_step = execution_realtime / real(1+step-setting%Time%WallClock%LastStepStored,kind=8)
            

        
            steps_to_finish = (setting%Time%End - setting%Time%Now) / dt
            seconds_to_completion = execution_realtime_per_step * steps_to_finish



            ! if (.not. inSpinUpYN) then
            ! print *, ' '
            ! print *, setting%Time%WallClock%Now, setting%Time%WallClock%TimeMarchStart
            ! print *, setting%Time%WallClock%Now - setting%Time%WallClock%TimeMarchStart
            ! print *, execution_realtime
            ! print *, execution_realtime_per_step
            ! print *, steps_to_finish
            ! print *, seconds_to_completion
            ! print *, setting%Time%End, setting%Time%Now, dt
            ! print *, step
            ! print *, ' '
            ! stop 4098734
            ! end if
        end if

        if (setting%Output%Verbose) then
            if (this_image() == 1) then
                if (mod(step,interval) == 0) then
                    thistime = timeNow
                    ! call util_datetime_display_time (thistime, timeunit)

                    ! ! write a time counter
                    ! if (.not. inSpinUpYN) then
                    !     if (dt > oneR) then
                    !         write(*,"(A12,i8,a17,F9.2,a1,a8,a6,f9.2,a3,a8,f9.2,a11,f9.2,a13,f9.2)") &
                    !             'time step = ',step,'; model time = ',thistime, &
                    !             ' ',trim(timeunit),'; dt = ',dt,' s', '; cfl = ',cfl_max!, &
                    !         !   '; cfl_CC = ',cfl_max_CC,'; cfl_JBJM = ',cfl_max_JBJM 
                    !     else
                    !         write(*,"(A12,i8,a17,F9.4,a1,a8,a6,f9.8,a3,a8,f9.2,a11,f9.2,a13,f9.2)") &
                    !             'time step = ',step,'; model time = ',thistime, &
                    !             ' ',trim(timeunit),'; dt = ',dt,' s', '; cfl = ',cfl_max
                    !     end if
                    ! else
                    !     write(*,"(A15,i8,a17,f9.2,a1,a8,a6,f9.2,a3,a8,f9.2)") &
                    !         'spin-up step = ',step,'; model time = ',thistime, &
                    !       ' ',trim(timeunit),'; dt = ',dt,' s', '; cfl = ',cfl_max!, &
                    !     !   '; cfl_CC = ',cfl_max_CC,'; cfl_JBJM = ',cfl_max_JBJM 
                    ! end if
                    ! if (.not. inSpinUpYN) then
                    !     ! write estimate of time remaining
                    !     thistime = seconds_to_completion
                    !     call util_datetime_display_time (thistime, timeunit)
                    !     write(*,"(A9,F10.2,A1,A3,A)") 'estimate ',thistime,' ',timeunit,' wall clock time until completion'
                    !     !write(*,"(A9,F6.2,A1,A3,A)") 'execution time ',thistime,' ',timeunit,' wall clock time thus far'
                    ! end if    
                    ! print *
                endif
            endif
        endif

        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine tl_command_line_step_output
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function tl_get_max_CFL_CC(thisCol,dt) result (outvalue)
        !%-------------------------------------------------------------------
        !% Description:
        !% computes the maximum CFL for the cells masked by "thisCol" in
        !% npack_elem(thisCol)
        !%-------------------------------------------------------------------
            integer, intent(in) :: thisCol
            real(8), intent(in), target :: dt
            real(8), pointer    :: thisDT
            integer, pointer :: Npack, thisP(:)
            real(8), pointer :: velocity(:), wavespeed(:), length(:), PCelerity(:)
            integer :: itemp(1), ip, fup, fdn, eup, edn, ii

            real(8), parameter :: smallDepth = 0.0001
        !%-------------------------------------------------------------------
            Npack              => npack_elemP(thisCol)
            thisP              => elemP(1:Npack,thisCol)
            velocity           => elemR(:,er_Velocity)
            wavespeed          => elemR(:,er_WaveSpeed)
            length             => elemR(:,er_Length)
            PCelerity          => elemR(:,er_Preissmann_Celerity)
        !%------------------------------------------------------------------- 
        if (dt .le. zeroR) then 
            thisDT => setting%Time%Hydraulics%Dt
        else    
            thisDT => dt
        end if

        ! do ii=1,N_elem(1)
        !     print *, ii, trim(reverseKey(elemI(ii,ei_elementType))), elemR(ii,er_WaveSpeed)
        ! end do

            ! print *, ' '
            ! print *, 'in tl_get_max_CFL_CC, Npack = ', Npack

            ! print *, ' element type '
            ! print *, elemI(thisP,ei_elementType)
            ! print *, trim(reverseKey(elemI(thisP(1),ei_elementType)))

            ! print *, 'depth '
            ! print *, elemR(thisP,er_Depth)
            ! do ii=1,size(thisP)
            !     print *, thisP(ii),elemR(thisP(ii),er_Depth)
            ! end do

            ! print *, 'head '
            ! print *, elemR(thisP,er_Head)

            ! print *, 'thisP '
            ! print *, thisP
            ! print *, ' '
            ! print *, 'Velocity  '
            ! print *, velocity(thisP)
            ! print *, ' '
            ! print *, 'WaveSpeed '
            ! print *, ' '
            ! print *, wavespeed(thisP)
            ! print *, ' '
            ! print *, 'Preiss C  '
            ! print *, PCelerity(thisP)
            ! print *, ' '
            ! print *, 'max values'
            ! print *, maxval(velocity(thisP)), maxval(wavespeed(thisP)) , maxval(PCelerity(thisP))

            ! print *,  ' '
            ! print *, 'wavespeed from iet'
            ! print *, wavespeed(iet)

            ! write(*,"(A,30f12.5)") 'Vcfl ' , velocity(iet) * thisDT / length(iet)
            ! write(*,"(A,30f12.5)") 'Hcfl ' , wavespeed(iet) * thisDT / length(iet)
            ! write(*,"(A,30f12.5)") 'Ccfl ' , PCelerity(iet) * thisDT / length(iet)
            ! print * ,' '
            ! print *, 'thisP '
            ! print *, thisP
            ! write(*,"(A,30f12.5)") 'Vcfl ' , velocity(thisP) * thisDT / length(thisP)
            ! write(*,"(A,30f12.5)") 'Hcfl ' , wavespeed(thisP) * thisDT / length(thisP)
            ! print *, 'thisDT = ',thisDT
            ! print *, length(thisP)
            ! print *, wavespeed(thisP)
            ! write(*,"(A,30f12.5)") 'Ccfl ' , PCelerity(thisP) * thisDT / length(thisP)

        !% HACK experiment to evaluate effect of JM wavespeed on time step
        elemR(:,er_Temp01) = wavespeed
        where (elemI(thisP,ei_elementType) .eq. JM)
            wavespeed(thisP) = zeroR
        endwhere

        !% --- set the outvalue
        if (Npack > 0) then 
            if (setting%Solver%PreissmannSlot%useSlotTF) then
                !% --- choose between maximum of the advective+wavespeed CFL or the 
                !%     Preissmann Slot celerity
                outvalue = max (maxval((abs(velocity(thisP)) + abs(wavespeed(thisP))) * thisDT / length(thisP)), &
                                maxval((abs(velocity(thisP)) + abs(PCelerity(thisP))) * thisDT / length(thisP)))
                ! print *, 'CFL options ', &
                !     maxval((abs(velocity(thisP)) + abs(wavespeed(thisP))) * thisDT / length(thisP)), &
                !     maxval((abs(PCelerity(thisP))) * thisDT / length(thisP))              
            else 
                outvalue = maxval((abs(velocity(thisP)) + abs(wavespeed(thisP))) * thisDT / length(thisP))
            end if
            ! print *, ' '
            ! print *, 'max vel, wave ', maxval(abs(velocity(thisP))), maxval(abs(wavespeed(thisP)))  
            ! print *, ' '      
        else
            outvalue = zeroR
        end if

        !% HACK EXPERIMENT CHANGING JM WaveSpeed -- KEEP THIS 20230505
        where (elemI(thisP,ei_elementType) .eq. JM)
            wavespeed(thisP) = elemR(thisP,er_Temp01)
        endwhere
        elemR(:,er_Temp01) = zeroR

        ! if (Npack == zeroI) then 
        !     !% --- no nonzero depths, so control time step by inflows
        !     outvalue = maxval(elemR(:,er_FlowrateLateral) * thisDT / (elemR(:,er_Length) * elemR(:,er_Breadthmax) * smallDepth))
        ! end if

        !     print *, 'max Q lateral ',maxval(elemR(:,er_FlowrateLateral) * thisDT / (elemR(:,er_Length) * elemR(:,er_Breadthmax) * smallDepth))

                ! print *, ' '
                ! print *, 'ThisDT ',thisDT
                ! print *, 'outvalue CFL in tl_get_max_CFL_CC',outvalue
                ! print *, ' '

            ! print *, ' '
            ! print *, 'in tl_get_max_CFL_CC'
            ! print *, length(thisP)
            ! print *, abs(wavespeed(iet(1))) * thisDT / length(iet(1))
            ! print *, outvalue
            ! print *, ' '

            !stop 298734

    end function tl_get_max_CFL_CC   
 !%
!%==========================================================================
!%==========================================================================
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_limit_BCinflow_dt (thisDT)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the maximum dt for inflows (BCup) into cells that
        !% have small depths
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), intent(inout) :: thisDT

            integer, pointer :: BCelemIdx(:)
            real(8), pointer :: depthLimit(:), DTlimit(:), depthScale(:)
            real(8), pointer :: BCQ(:), topwidth(:), length(:), depth(:)
            real(8), pointer :: alpha, gravity, smallDepth
            real(8) :: newDTlimit
        !%------------------------------------------------------------------
        !% Preliminaries:
            !if (crashYN) return
            temp_BCupI(:,:) = nullValueI
            temp_BCupR(:,:) = nullValueR
        !%------------------------------------------------------------------
        !% Aliases:
            BCelemIdx  => temp_BCupI(:,1)

            depthLimit => temp_BCupR(:,1)
            DTlimit    => temp_BCupR(:,2)
            BCQ        => temp_BCupR(:,3)
            depthScale => temp_BCupR(:,4)

            topwidth   => elemR(:,er_TopWidth)
            length     => elemR(:,er_Length)
            depth      => elemR(:,er_Depth)

            alpha      => setting%VariableDT%CFL_inflow_max
            gravity    => setting%Constant%gravity
            smallDepth => setting%SmallDepth%MomentumDepthCutoff
        !%------------------------------------------------------------------
        
        !% set the timestep limiter for upstream inflow boundary conditions (BCup)
        ! print *, ' '
        ! print *, 'this dt into limiter ',thisDT
        ! print *, 'size  in tl_limit_BCinflow_dt',size(BC%P%BCup)

        if (size(BC%P%BCup) > 0) then

            !% ensure flowrate used for limiter is  positive and non-zero
            BCQ(:) = max( abs(BC%flowR(BC%P%BCup,br_value)), setting%Eps%Velocity)

            !    print *, 'BCQ ',BCQ(:)

            !% store the element indexes downstream of a BCup face
            BCelemIdx =  faceI(BC%flowI(BC%P%BCup,bi_face_idx), fi_Melem_dL)

            !% store the scaling depth limit defined by when the induced velocity
            !% from an inflow is similar to the gravity wave speed of the BC inflow
            depthScale = ( (alpha**4) * (BCQ**2) / (gravity * (topwidth(BCelemIdx)**2) ) )**onethirdR

                ! print *, ' '
                ! print *, 'alpha ',alpha
                ! print *, 'BCQ   ',BCQ 
                ! print *, BCelemIdx, topwidth(BCelemIdx)

                ! print *, 'depthScale ',depthScale
                ! print *, 'depth      ',elemR( BCelemIdx,er_Depth)

            !% get the depth limit for maximum depth that the time step will be limited as 
            !% either the depth scale or the specified smallDepth limit
            depthLimit = max(smallDepth, depthScale)

                ! print *, 'depthLimit ',depthLimit
            
            !% time step limit based on inflow flux
            ! print *, 'length  ', length(BCelemIdx)
            ! print *, 'alpha   ', alpha
            ! print *, 'topwidth', topwidth(BCelemIdx) 
            ! print *, 'BCQ     ', BCQ

            DTlimit = length(BCelemIdx) * ( ( alpha * topwidth(BCelemIdx) / (gravity * BCQ) )**onethirdR) 

                ! print *, ' '
                ! print *, 'DTlimit ',DTlimit
                ! print *, 'length  ',length(BCelemIdx)
                ! print *, 'term    ',( ( alpha * topwidth(BCelemIdx) / (gravity * BCQ) )**onethirdR) 
                ! print *, ' '
        
            !% where the depth is greater than the depthlimit the DT inflow limiter
            !% is not needed, and we can use the existing DT value
            depthLimit = depth(BCelemIdx) - depthLimit
            where (depthLimit .ge. zeroR)
                DTlimit = thisDT
            endwhere

                ! print *, 'new 1 DTlimit ',DTlimit

            !% get the smallest DT in the limiter array
            newDTlimit = minval(DTlimit)

                ! print *, 'new 2 DTlimit ',DTlimit

            !% use the smaller value of the new limit or the input
            thisDT = min(newDTlimit,thisDT) 

            !% return to null value storage
            temp_BCupI(:,:) = nullValueI
            temp_BCupR(:,:) = nullValueR
        else
            !% continue -- no DT limitation if there are no BCup faces
        end if

        !%------------------------------------------------------------------
    end subroutine tl_limit_BCinflow_dt
!%
!%==========================================================================
!%==========================================================================
!%   
    subroutine tl_limit_BChead_inflow_dt (thisDT)
        !%------------------------------------------------------------------
        !% Description: 
        !% Computes the maximum DT allowed for a head BC that is providing
        !% an inflow
        !%------------------------------------------------------------------
        !% Declarations
            real(8), intent(inout) :: thisDt
            integer :: ii
            integer, pointer :: fidx, elemUp
            real(8), pointer :: alpha, gravity, smallDepth
            real(8), pointer :: Qface, length, depth, topWidth
            real(8) :: depthScale, depthLimit, DTlimit
        !%------------------------------------------------------------------
        !% Aliases    
            alpha      => setting%VariableDT%CFL_inflow_max
            gravity    => setting%Constant%gravity
            smallDepth => setting%SmallDepth%MomentumDepthCutoff
        !%------------------------------------------------------------------    

        do ii=1, N_headBC

            ! print *, ' '
            ! print *, 'BC Head inflow limit'

            !% face index
            fidx     => BC%headI(ii,bi_face_idx)
            !% face flowrate
            Qface    => faceR(fidx,fr_Flowrate)
            !% upstream element index
            elemUp   => faceI(fidx, fi_Melem_uL)
            !% topwidth of upstream element
            topWidth => elemR(elemUp,er_TopWidth)
            !% length of upstream element
            length   => elemR(elemUp,er_Length)
            !% depth of upstream element
            depth    => elemR(elemUp,er_Depth)

            if (Qface .ge. zeroR) exit  !% no dt limit if this is an outflow

            !% scaling depth limit defined by when the induced velocity
            !% from an inflow is similar to the gravity wave speed of the BC inflow
            depthScale = ( (alpha**4) * (Qface**2) / (gravity * (topWidth**2) ) )**onethirdR

            ! print *, 'depthScale ',depthScale
            ! print *, 'depth      ',elemR( BCelemIdx,er_Depth)

            !% get the depth limit for maximum depth that the time step will be limited as 
            !% either the depth scale or the specified smallDepth limit
            depthLimit = max(smallDepth, depthScale)

            ! print *, 'depthLimit ',depthLimit
            
            !% time step limit based on inflow flux
            ! print *, 'length  ', length(BCelemIdx)
            ! print *, 'alpha   ', alpha
            ! print *, 'topwidth', topwidth(BCelemIdx) 
            ! print *, 'BCQ     ', BCQ

            DTlimit = length * ( ( alpha * topwidth / (gravity * (-Qface)) )**onethirdR) 

            DTlimit = min(DTlimit, thisDT)

            ! print *, 'DTlimit ',DTlimit
        
            !% where the depth is greater than the depthlimit the DT inflow limiter
            !% is not needed, and we can use the existing DT value
            depthLimit = depth - depthLimit
            if (depthLimit .ge. zeroR) then
                DTlimit = thisDT
            end if

            thisDT = DTlimit

        end do

    end subroutine tl_limit_BChead_inflow_dt
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_limit_LatInflow_dt (thisDT)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes the maximum dt for inflows (Qlat) into cells that
        !% have small depths
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), intent(inout) :: thisDT
            real(8), pointer :: Qlat(:), PlanArea(:), HeightRise(:) !, depthScale(:), depthLimit(:)
            real(8), pointer :: DTlimit(:) !, topwidth(:), length(:), depth(:)
            !real(8), pointer :: alpha, gravity, smallDepth
            real(8) :: newDTlimit

            integer, pointer :: thisP(:), Npack
            integer :: thisCol, ii
        !%------------------------------------------------------------------
        !% Preliminaries:
            !if (crashYN) return
        !%------------------------------------------------------------------
        !% Aliases:
            Qlat       => elemR(:,er_Temp01)
            PlanArea   => elemR(:,er_Temp02)
            HeightRise => elemR(:,er_Temp03)
            DTLimit    => elemR(:,er_Temp04)

            ! !% --- use the largest scale of topwidth
            ! !%     otherwise, as D => 0 and T => 0 we have issues
            ! topwidth   => elemR(:,er_BreadthMax)
            ! length     => elemR(:,er_Length)
            ! depth      => elemR(:,er_Depth)

            ! alpha      => setting%VariableDT%CFL_inflow_max
            ! gravity    => setting%Constant%gravity
            !smallDepth => setting%SmallDepth%MomentumDepthCutoff

        !%------------------------------------------------------------------
        !% ---- use the pack for CC and JM with H time march 
        !%      (no lateral flows into JB or diagnostic)


        thisCol = ep_CCJM_H
        Npack => npack_elemP(thisCol)   

        !print *, 'Npack here ',Npack

        if (Npack < 1) return
        thisP => elemP(1:Npack,thisCol)

        !% --- initialize temporary arrays
        PlanArea   = zeroR
        HeightRise = zeroR
        Qlat       = zeroR
        DTLimit    = huge(zeroR)

        !% --- store positive flowrate for limiter
        Qlat(thisP) = abs( elemR(thisP,er_FlowrateLateral)) 

        ! print *, 'Qlat ',Qlat(thisP)

        !% --- get the geometry associated with a simple inflow
        where (Qlat > zeroR) 
            PlanArea   = elemR(:,er_Length) * elemR(:,er_BreadthMax)
            HeightRise = Qlat * thisDT / PlanArea
        endwhere

        ! print *, ' '
        ! print *, 'PlanArea ',PlanArea(thisP)
        ! print *, ' '
        ! print *, 'HeightRise ',HeightRise(thisP)

        !% --- if th depth is small and the heighrise is large, then limit the time step
        !%     to a height rise of twice the local small depth value
        where ((HeightRise > setting%SmallDepth%LateralInflowSmallDepth)  &
                .and.                                                     &
                (elemR(:,er_Depth) < fourR * setting%SmallDepth%LateralInflowSmallDepth))

            DTlimit = twoR * setting%SmallDepth%LateralInflowSmallDepth * PlanArea / Qlat
        endwhere

        ! print *, ' '
        ! print *, 'DTlimit ',DTlimit(thisP)

        newDTlimit = minval(DTlimit)

        ! print *, 'new dt limit ', newDTlimit

        thisDT = min(thisDT, newDTlimit)

        ! print *, ' '
        ! print *, 'limited DT ', thisDT

        ! if (newDTlimit < 10000.0) then 
        !     stop 6098723
        ! end if


        !% OBSOLETE BELOW 20230406 brh
        ! print *
        ! print *, 'topwidth ',topwidth(thisP)

        ! !% store the scaling depth limit defined by when the induced velocity
        ! !% from an inflow is similar to the gravity wave speed of the BC inflow
        ! depthScale(thisP) = ( (alpha**4) * (Qlat(thisP)**2) / (gravity * (topwidth(thisP)**2) ) )**onethirdR

        ! print *
        ! print *, 'depthScale ',depthScale(thisP)

        ! !% get the depth limit for maximum depth that the time step will be limited as 
        ! !% either the depth scale or the specified smallDepth limit
        ! depthLimit(thisP) = max(smallDepth, depthScale(thisP))

        ! print *
        ! print *, 'depthLimit ',depthLimit(thisP)
        
        ! !% time step limit based on inflow flux
        ! where (Qlat(thisP) > setting%Eps%Velocity)
        !     DTlimit(thisP) = length(thisP) * ( ( alpha * topwidth(thisP) / (gravity * Qlat(thisP)) )**onethirdR) 
        ! elsewhere
        !     DTlimit(thisP) = nullvalueR
        ! endwhere

        ! print *
        ! print *, 'DTlimit ',DTlimit(thisP)

        ! do ii=1,size(thisP)
        !     if (Qlat(thisP(ii)) > 0.d0) then
        !         write(*,"(i5,4(A,e12.5))") thisP(ii), ' ',Qlat(thisP(ii)), ' ',topwidth(thisP(ii)), ' ',depthLimit(thisP(ii)),' ',DTlimit(thisP(ii))
        !     end if
        ! end do
    
        ! !% where the depth is greater than the depthlimit the DT inflow limiter
        ! !% is not needed, and we can use the existing DT value
        ! depthLimit(thisP) = depth(thisP) - depthLimit(thisP)
        ! where (depthLimit .ge. zeroR)
        !     DTlimit = thisDT
        ! endwhere

        ! !print *
        ! !print *, 'new DTlimit ',DTlimit(thisP)

        ! !% get the smallest DT in the limiter array
        ! newDTlimit = minval(DTlimit(thisP))

        ! !% use the smaller value of the new limit or the input DT
        ! thisDT = min(newDTlimit,thisDT)

        ! !print *
        ! !print *, 'thisDT ',thisDT

        ! !%------------------------------------------------------------------
        ! !% Closing
        !     elemR(:,er_Temp01) = nullvalueR
        !     elemR(:,er_Temp02) = nullvalueR
        !     elemR(:,er_Temp03) = nullvalueR
        !     elemR(:,er_Temp04) = nullvalueR

    end subroutine tl_limit_LatInflow_dt
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_smallestBC_timeInterval ()
        !%------------------------------------------------------------------
        !% Description
        !% gets the smallest time interval in the latest BC data
        !%------------------------------------------------------------------
        !% Declarations
            real(8) :: smallHead, smallFlow
        !%------------------------------------------------------------------

        if (N_headBC > 0) then
            smallHead = minval(BC%headR(:,br_timeInterval))
        else
            smallHead = abs(nullvalueR)
        end if
        
        if (N_flowBC > 0) then
            smallFlow = minval(BC%flowR(:,br_timeInterval))
        else    
            smallFlow = abs(nullvalueR)
        end if

        setting%BC%smallestTimeInterval = min(smallFlow, smallHead)

    end subroutine tl_smallestBC_timeInterval
!% 
!%==========================================================================
!%==========================================================================
!%
    !% OBSOLETE APPROACH
    ! subroutine control_evaluate()
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Evaluates control and updates elemR(:,er_Setting) column
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer          :: ii, loc
    !         integer, pointer :: nControls, eIdx
    !         real(8), pointer :: TimeNow, TargetSetting(:), TimeArray(:), SettingArray(:)
    !     !%------------------------------------------------------------------
    !     !% Preliminaries:
    !     ! if (crashYN) return
    !     !%------------------------------------------------------------------
    !     !% Aliases:
    !     TimeNow   => setting%Time%Now 
    !     nControls => setting%Control%NumControl
    !     TargetSetting  => elemR(:,er_TargetSetting)

    !     !% only use controls if it is present in the settings file
    !     if (nControls > zeroI) then
    !         do ii = 1,nControls
    !             eIdx          => setting%Control%ElemIdx(ii)
    !             TimeArray     => setting%Control%TimeArray(:,ii)
    !             SettingArray  => setting%Control%SettingsArray(:,ii)

    !             !% now find where the current time (TImeNow) falls between
    !             !% the control time array (TimeArray)
    !             !% here, it is found using the maxloc Intrinsic function
    !             !% Returns the location of the minimum value of all elements 
    !             !% in an array, a set of elements in an array, or elements in 
    !             !% a specified dimension of an array. This function works only when
    !             !% the current time is above the minimum value in the TimeArray array
    !             if (TimeNow > minval(TimeArray)) then
    !                 loc = maxloc(TimeArray, 1, TimeArray <= TimeNow)
    !                 !% setting can not be greater than 1
    !                 TargetSetting(eIdx) = min(SettingArray(loc), oneR)
    !             end if
    !         end do
    !     end if

    ! end subroutine control_evaluate 
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module timeloop

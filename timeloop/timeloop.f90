module timeloop
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Provides outer loop of time march
    !%
    !% Methods:
    !% Outer loop calls both hydrology and hydraulics
    !% Determines the variable time step needed for the next hydraulics step
    !% Provides hydraulic spin-up simulations for steady-flow initial conditions
    !% Note that hydrologic spin-up is not presently supported
    !%==========================================================================
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


    implicit none
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
        !% Loops over all the major time-stepping routines
        !% Exit returns control to Main.f90
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

        !% --- set spinup controls and call spinup
        call tl_spinup()

        !% --- if stopping after spin-up, end here
        if (setting%Simulation%stopAfterSpinUp) return

        !% --- reset the start time after spin-up
        setting%Time%Start = startTime
        setting%Time%Now   = startTime
        setting%Time%End   = endTime
        setting%Output%Report%StartTime = reportStart

        sync all
        if (this_image()==1) then
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%TimeMarchStart = cval
            setting%Time%WallClock%LastTimeStored = cval
            setting%Time%WallClock%LastStepStored = setting%Time%Step
        end if 

        !% --- initialize the time settings for hydraulics and hydrology steps
        call tl_initialize_loop (doHydraulicsStepYN, doHydrologyStepYN, .false.)

        !-- perform the time-marching loop
        call tl_outerloop (doHydrologyStepYN, doHydraulicsStepYN, .false., .false.)

        sync all
        !% --- close the timemarch time tick
        if (this_image() == 1) then
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%TimeMarchEnd= cval
        end if

        !%--------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine timeloop_toplevel
!% 
!%==========================================================================
!% PRIVATE - 1st Level
!%==========================================================================
!%
    subroutine tl_spinup ()
        !% -----------------------------------------------------------------
        !% Description:
        !% Conducts spin-up simulations holding the time=0 BC fixed for the
        !% number of days set in setting%Simulation%SpinUpDays
        !% At finish, control returns to timeloop_toplevel() for continued
        !% simulation
        !% -----------------------------------------------------------------
        !% Declarations:
            logical :: inSpinUpYN, SpinUpOutputYN, doHydraulicsStepYN, doHydrologyStepYN
        !% -----------------------------------------------------------------
        !% Preliminaries
            if (.not. setting%Simulation%useSpinUp) return
        !% -----------------------------------------------------------------

        inSpinUpYN = .true.        
        setting%Time%Now = setting%Time%Start
        setting%Time%End = setting%Simulation%SpinUpDays * seconds_per_day

        !% --- only provide output during spinup if simulation stops after spinup  
        if (setting%Simulation%stopAfterSpinUp) then     
            SpinUpOutputYN = .true.
            !% --- default to report for all spinup
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
        !% initialize the time varaibles before a simulation
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

        !% --- local T/F that changes with each time step depending on whether or 
        !%     not the hydrology or hydraulics are conducted in that step
        doHydrologyStepYN  = setting%Simulation%useHydrology
        doHydraulicsStepYN = setting%Simulation%useHydraulics

        !% --- get the next hydrology time
        if (setting%Simulation%useHydrology) then      
            nextHydrologyTime  = interface_get_NewRunoffTime()
        else  
            !% set a large dummy time for hydrology if not used
            nextHydrologyTime  = setting%Time%End + onethousandR * dtTol
        end if
        
        if ((nextHydrologyTime == lastHydrologyTime) .and. &
            (setting%Simulation%useHydrology) ) then
            !% --- call the first hydrology step
            call tl_hydrology()
            nextHydrologyTime  = interface_get_NewRunoffTime()
        end if

        !% --- get the initial dt and the next hydraulics time
        if (setting%Simulation%useHydraulics) then
            call tl_smallestBC_timeInterval ()
            call tl_update_hydraulics_timestep(inSpinUpYN)
            
        else
            !% ---- HACK - WORKING WITHOUT SWMM5+ HYDRAULICS IS NOT SUPPORTED 
            !%      the following are stub routines for future design
            !%      set a large dummy time for hydraulics if not used
            nextHydraulicsTime = setting%Time%End + onethousandR * dtTol
            !% --- suppress the multi-level hydraulics output (otherwise seg faults)
            setting%Output%Report%suppress_MultiLevel_Output = .true.
            print *, 'CODE ERROR: SWMM5+ does not operate without hydraulics'
            print *, 'Ensure that setting.Simulation.useHydraulics = .true.'
            call util_crashpoint(268743)
        end if

        call util_crashstop(229873)

        !% --- set the next control rule evaluation time
        nextControlRuleTime = lastControlRuleTime + real(setting%SWMMinput%ControlRuleStep,8)

        !% --- check to see if there is an initial step to hydrology
        !%     if not, then skip the initial tl_hydrology ()
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
        !% The outer do while loop that links hydrology and hydraulics
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

        !% --- outer loop of the time-march
        do while (setting%Time%Now <= setting%Time%End - dtTol)

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

            !%------------------
            !%  HYDROLOGY
            !%------------------
            if ((.not. inSpinUpYN) .or. (inSpinUpYN .and. BCupdateYN) ) then
                !% --- note that hydrology in SWMM includes climate update
                if (doHydrologyStepYN) then 
                    !% --- store the runoff from hydrology on a hydrology step
                    call tl_hydrology()
                    setting%Climate%EvapRate = interface_get_evaporation_rate()
                else 
                    !% --- update climate if evaporation is needed in hydraulics-only simulations
                    if (setting%Climate%useHydraulicsEvaporationTF) then 
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

            !%-----------------
            !%  HYDRAULICS
            !%----------------
            if (doHydraulicsStepYN) then    

                !% --- set a clock tick for hydraulic loop evaluation
                if ((this_image()==1) .and. (.not. inSpinUpYN)) then
                    call system_clock(count=cval,count_rate=crate,count_max=cmax)
                    setting%Time%WallClock%HydraulicsStart = cval
                end if 

                !% --- get updated boundary conditions
                if (BCupdateYN) then
                    call bc_update() 
                    call tl_lateral_inflow()
                    call tl_smallestBC_timeInterval ()
                end if

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

                !% --- add subcatchment inflows
                !%     note, this has "useHydrology" and not "doHydrologyStepYN" because the
                !%     former is whether or not hydrology is used, and the latter is
                !%     whether or not it is computed in this time step
                if (setting%Simulation%useHydrology .and. BCupdateYN) then 
                    call tl_subcatchment_lateral_inflow () 
                end if

                !% --- add RDII inflows
                if (.not. setting%SWMMinput%IgnoreRDII) then 
                    call interface_get_RDII_inflow ()
                end if

                !% --- add groundwater inflows
                if ((.not. setting%SWMMinput%IgnoreGroundwater) .and. &
                            (setting%SWMMinput%N_groundwater > 0) ) then 
                    call interface_get_groundwater_inflow ()
                end if 

                !% --- set hydraulics time step to handle inflow
                call tl_update_hydraulics_timestep(.false.)

                !% --- perform one hydraulic routing step
                call tl_hydraulics()
            
                !% --- accumulate RunOn from hydraulic elements to subcatchments
                if ((setting%Simulation%useHydrology) .and. (any(subcatchYN(:,sYN_hasRunOn)))) then 
                    call tl_subcatchment_accumulate_runon ()
                end if

                !% --- close the clock tick for hydraulic loop evaluation
                if ((this_image()==1) .and. (.not. inSpinUpYN)) then
                    call system_clock(count=cval,count_rate=crate,count_max=cmax)
                    setting%Time%WallClock%HydraulicsStop = cval
                    setting%Time%WallClock%HydraulicsCumulative &
                        = setting%Time%WallClock%HydraulicsCumulative &
                        + setting%Time%WallClock%HydraulicsStop &
                        - setting%Time%WallClock%HydraulicsStart
                end if 

            end if         

            !%---------------
            !% OUTPUT
            !%-------------
            if (setting%Output%Report%provideYN) then 
                !% --- only provide output for spinup time if stopping after spinup
                    if ( (.not. inSpinUpYN)  &
                        .or. ((inSpinUpYN) .and. (SpinUpOutputYN)) ) then

                    !% set a time tick for output timing
                    if ((this_image()==1) .and. (.not. inSpinUpYN)) then
                        call system_clock(count=cval,count_rate=crate,count_max=cmax)
                        setting%Time%WallClock%LoopOutputStart = cval
                    end if

                    !% --- Report summary results in report HACK --NOT IMPLEMENTED
                    !%     Note must be reported before the "do"counter increments
                    !%     ARCHIVE FOR FUTURE USE
                    !call util_output_report()  --- this is a stub for future use

                    !%--- Multilevel time step output
                    if ((util_output_must_report()) .and. &
                        (.not. setting%Output%Report%suppress_MultiLevel_Output) ) then
                        call outputML_store_data (.false.)
                    end if
                    sync all

                    !% --- close the time tick for output timing
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

            !% --- check for crash conditions
            call util_crashstop(13978)

            !% --- restart the hydraulics time tick
            sync all
            if ((this_image()==1) .and. (.not. inSpinUpYN) ) then
                call system_clock(count=cval,count_rate=crate,count_max=cmax)
                setting%Time%WallClock%HydraulicsStart = cval
            end if 
    
            sync all
            !% ---increment the time step and counters for the next time loop
            call tl_increment_timestep_and_counters(doHydraulicsStepYN, doHydrologyStepYN, inSpinUpYN)

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

            !% --- check for crash conditions
            call util_crashstop(63978)
    
        end do  !% --- end of time loop

    end subroutine tl_outerloop
!% 
!%==========================================================================
!% PRIVATE - 2nd Level
!%==========================================================================
!%
    subroutine tl_hydrology()
        !%------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step
        !%------------------------------------------------------------------
        !% Declarations
            integer          :: ii
            real(8), pointer :: sRunoff(:), dt
            integer(kind=8)  :: crate, cmax, cval
            character(64)    :: subroutine_name = 'tl_hydrology'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            !% --- wall clock tick
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
            !% --- HACK hydrology spinup nto supported
            print *, 'CODE ERROR: setting.simulation.useSpinUp = true is not supported for hydrology at this time'
            call util_crashpoint(4482333)
        end if

        !% --- get the runon flowrate (if any) from hydraulics
        !%     HACK -- runon needs further testing
        if (any(subcatchYN(:,sYN_hasRunOn))) then
            call tl_subcatchment_set_runon_volume ()
        end if

        !% --- only execute runoff if the next time is sufficently small
        if (setting%Time%Hydrology%NextTime < (setting%SWMMinput%TotalDuration - setting%Time%DtTol)) then

            !% --- execute the EPA SWMM hydrology for the next interval 
            call interface_call_runoff_execute()

            !% --- get the next runoff time
            setting%Time%Hydrology%NextTime = interface_get_NewRunoffTime()  

            !% --- cycle through the subcatchments to get the runoff 
            !%     ii-1 required in arg as C arrays start from 0
            do ii = 1,setting%SWMMinput%N_subcatch
                sRunoff(ii) = interface_get_subcatch_runoff(ii-1)  
            end do
        else 
            sRunoff(:) = zeroR
        end if

        !%------------------------------------------------------------------ 
        !% Closing  
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
            integer              :: ii, kk
            real(8)              :: tvolume
            character(64)        :: subroutine_name = 'tl_hydraulics'
        !%-------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% --- repack all the dynamic arrays for new conditions
        call pack_dynamic_arrays()

        !% --- ensure that the conservative flux terms are exactly zero in the entire array
        !%     so that we can be confident of conservation computation. 
        faceR(:,fr_Flowrate_Conservative) = zeroR  

        !% --- call the RK2 time march
        call rk2_toplevel ()

        !% --- add non-conservation in this step to accumulator
        call util_accumulate_volume_conservation () 

        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine tl_hydraulics
!% 
!%==========================================================================
!% PRIVATE - 3rd Level
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
            real(8)          :: thisEpochTime, confac
            integer          :: year, month, day
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
        !% --- set lateral inflow to zero for all cells
        Qlateral(:) = zeroR 

        !%-----------------
        !% LATERAL INFLOW
        !%-----------------
        !% --- if ep_BClat exist add lateral inflow BC to lateral inflow accumulator
        !%     note that thisP and thisBC must be the same size or there is something wrong  
        !%     For multi-barrel elements, divide lateral inflow evenly between barrels
        npack    => npack_elemP(ep_BClat)
        thisP    => elemP(1:npack,ep_BClat)
        if (npack > 0) then    
            Qlateral(thisP) = Qlateral(thisP) &
                + BC%flowR(thisBC,br_value) / real(nBarrels(thisP),8)
        end if

        !%---------------
        !% CONDUCTIVITY
        !%---------------
        !% --- find the Adjust.Conductivity multiplier for the current month
        !%     The setting.Adjust.Conductivity is from EPA SWMM ADJUSTMENTS
        thisEpochTime = util_datetime_secs_to_epoch(setting%Time%Now)
        call util_datetime_decodedate(thisEpochTime, year, month, day)
        confac = setting%Adjust%Conductivity(month)

        !%--------------
        !% SEEPAGE
        !%--------------
        !% --- Subtract the seepage rate for all conduit and channel cells
        !%     HACK -- instead of using perimeter, we compute seepage based on 
        !%     element length and nominal topwidth. If we want to convert to a
        !%     Perimeter-based solution we need new geometry to store perimeter.
        !%
        !%     Convert m/s seep to m^3/s flux by multiplying by topwidth and length
        !%     Note that "Topwidth" here is either the actual topwidth or the maximum value if the present
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

        !%---------------------
        !% EXFILTRATION from STORAGE
        !%--------------------
        !% HACK - NEED EXFILTRATION FOR STORAGE 
        ! routing.c  calls node_getLosses
        ! which calls node.c/node_getLosses
        ! which calls node.c/storage_getLosses
        ! which compute exfilRate as a local variable -- which is what we need

        ! IF EVERYTHING WE NEED IS STORED, THEN WE SHOULD BE ABLE TO 
        ! CALL exfil_getLoss() from an api.
        ! Must cycle through all the storage nodes.

        !%------------------
        !% EVAPORATION
        !%------------------
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

        else
            !% --- evaporation has no effect
        end if

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
        !% Note: All images run the entire hydrology, but only the image 
        !% holding the runoff node is affected.
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

        do mm = 1,setting%SWMMinput%N_subcatch
            !% --- only if this image holds this node
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
            dt => setting%Time%Hydraulics%Dt
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

            do kk = 1,nRunOn               
                if (subcatchI(mm,si_RunOn_P_image(kk)) == this_image()) then
                    !% --- change all values to a flowrate
                    subcatchR(mm,sr_RunOn_Volume(kk)) = subcatchR(mm,sr_RunOn_Volume(kk))
                else 
                    subcatchR(mm,sr_RunOn_Volume(kk)) = zeroR
                end if

                !% --- distribute the runon volume to all images
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
        !% time counters
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
        else
            nextHydraulicsTime = setting%Time%End + tenR*DtTol
        end if   
        
        call util_crashstop(449873)

        !% --- The NextHydrologyTime is updated in EPA SWMM-C, here we just need to
        !%     provide a large number if hydrology isn't used
        if (.not. useHydrology) then
            nextHydrologyTime = max(setting%Time%End + tenR*DtTol, &
                                    setting%Simulation%SpinUpDays * seconds_per_day)
        else
            !% --- EPA SWMM-C will not return a next hydrology time that is
            !%     beyond the end of the simulation, so we need a special
            !%     treatment for this case
            if ((nextHydrologyTime .eq. lastHydrologyTime) .and. &
                (timeNow > setting%Time%Start)) then 
                nextHydrologyTime = max(setting%Time%End + tenR*DtTol , &
                                        setting%Simulation%SpinUpDays * seconds_per_day)
            end if
        end if

        !% --- Ensure that all processors use the same time step.
        !%     Find the minimum hydraulics time and store accross all processors.
        call co_min(nextHydraulicsTime)
        !% --- Take the nextTime as the minimum of either the Hydrology or Hydraulics time
        !%     Done on a single processor because they all should have identical nextHydrologyTIme
        nextTime = min(nextHydraulicsTime, nextHydrologyTime)

        !% --- Note there is no need to broadcast across processors since co_min
        !%     ensures that each processor gets the min value.

        !% --- Update the time step for the local processor 
        !%     (all have the same nextTime and lastHydraulicsTime)
        dtHydraulics = nextTime - lastHydraulicsTime

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

        step    = step + oneI
        timeNow = nextTime !timeNow + dt

        !% --- provide output to the command line
        call tl_command_line_step_output(inSpinUpYN)

        !% --- push old times down to last times
        !%     HACK -- need to check the logic for this
        if (doHydraulicsStepYN) LastHydraulicsTime = NextHydraulicsTime
        if (doHydrologyStepYN)  LastHydrologyTime  = NextHydrologyTime

        !%--------------------------------------------------------------------
        !% Closing
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
            logical, pointer    :: matchHydrologyStep, useHydrology

            real(8)             :: oldDT, oldCFL
            real(8), pointer    :: newDT, timeNow
            real(8), pointer    :: nextHydraulicsTime, lastHydraulicsTime

            integer             :: neededSteps, pindex(1)
    
            character(64)        :: subroutine_name = 'tl_update_hydraulics_timestep'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases:
            useHydrology       => setting%Simulation%useHydrology
            matchHydrologyStep => setting%Time%matchHydrologyStep
            nextHydraulicsTime => setting%Time%Hydraulics%NextTime
            lastHydraulicsTime => setting%Time%Hydraulics%LastTime        
            timeNow            => setting%Time%Now
        !%----------------------------------------------------------------------
        !% --- note oldDT is NOT a alias as we want it fixed while newDT
        !%     is the alias for the next time step, which will change
        oldDT    =  setting%Time%Hydraulics%Dt  
        newDT    => setting%Time%Hydraulics%Dt

        !% --- get the CFL based on last time step as a control
        oldCFL = tl_get_max_CFL(ep_CCJM_NOTzerodepth,oldDT)

        
        if ((matchHydrologyStep) .and. (useHydrology) .and. (.not. inSpinUpYN)) then 
            !%-----------------------------------------------------------
            !% IF HYDRAULICS STEP SHOULD END EXACTLY ON A HYDROLOGY STEP
            !%-----------------------------------------------------------
            !% --- get the time step that is less than or equal to the next
            !%     hydrology time
            call tl_DT_hydrology_match (newDT, oldDT, oldCFL, neededSteps)
        else 
            !%-----------------------------------------------------------------------
            !% IF HYDRAULICS STEP IS NOT RELATED TO HYDROLOGY TIME
            !%-----------------------------------------------------------------------               
            if (useHydrology) then 
                print *, 'USER CONFIGURATION ERROR:' 
                print *, 'At this time, useHydrology requires '
                print *, 'setting%Time%matchHydrologyStep = .true.'
                print *, 'the setup for non-match has not been tested'
                call util_crashpoint(509872)
            end if
            !% --- get the dt
            call tl_DT_standard(newDT, oldDT, oldCFL)  

            !% --- neededSteps is irrelevant without hydrology matching,
            !%     but this using 3 forces a rounding operation below
            neededSteps = 3  
        end if

        !% --- invoke other time step limiters
        call tl_limit_DT (newDT)

        !% --- round off of time steps with too many digits
        call tl_roundoff_DT (newDT, neededSteps)
       
        !% increment the hydraulics time clock
        nextHydraulicsTime = lastHydraulicsTime + newDT

        !% --- store the global value of maximum CFL
        cfl_max = tl_get_max_CFL(ep_CCJM_NOTzerodepth,newDT)
    
        sync all
        !% --- distribute across images
        call co_max(cfl_max)

        !%----------------------------------------------------------------------
        !% Closing
            if ((setting%Limiter%Dt%UseLimitMinYN)      &
                .and.                                   &
                (newDT .le. setting%Limiter%Dt%Minimum) &
                .and.                                   &
                (setting%Limiter%Dt%FailOnMinYN)        &
                ) then

                print *, ' '
                print *, 'EXITING ON TIME STEP ERROR -- PROBABLY BLOWING UP DUE TO EXCESSIVE HEAD'
                print *,'timestep= ', setting%Time%Step
                print*, 'timeNow = ', timeNow, ' seconds'
                print*, 'dt = ', newDT, 'minDt = ',  setting%Limiter%Dt%Minimum
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

        !% --- reset the temporary arrays
        flowrate(:)    = nullvalueR
        volumeDelta(:) = nullvalueR
        timeToDepth(:) = nullvalueR

    end subroutine tl_dt_vanishingCFL
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_save_previous_values()
        !%------------------------------------------------------------------
        !% Description:
        !% Pushes the time N values into time N-1 storage, and the time N+1 values into
        !% the time N storage.
        !% HACK -- 20210809 brh This would be better done by changing the indexes without
        !% moving the data, but we will wait on doing this until we have the code
        !% debugged.
        !%------------------------------------------------------------------
        !% Declarations
            character(64) :: subroutine_name = 'tl_save_previous_values'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------     
        !%  push the old values down the stack
        !%  N values is the present, N0 is the last time step, and N1
        !%  is the timestep before (needed only for backwards 3rd in velocity and volume)
        elemR(:,er_Flowrate_N0)   = elemR(:,er_Flowrate)
        elemR(:,er_Head_N0)       = elemR(:,er_Head)
        elemR(:,er_Velocity_N1)   = elemR(:,er_Velocity_N0)
        elemR(:,er_Velocity_N0)   = elemR(:,er_Velocity)
        elemR(:,er_Volume_N1)     = elemR(:,er_Volume_N0)
        elemR(:,er_Volume_N0)     = elemR(:,er_Volume)
        elemR(:,er_Area_N1)       = elemR(:,er_Area_N0)
        elemR(:,er_Area_N0)       = elemR(:,er_Area)
        elemR(:,er_SlotVolume_N0) = elemR(:,er_SlotVolume)
        elemR(:,er_SlotDepth_N0)  = elemR(:,er_SlotDepth)
        !%------------------------------------------------------------------
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine tl_save_previous_values
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_command_line_step_output (inSpinUpYN)
        !%------------------------------------------------------------------
        !% Description:
        !% provides output at the commandline
        !%------------------------------------------------------------------
            logical, intent(in) :: inSpinUpYN
        
            integer,         pointer :: interval
            integer(kind=8), pointer :: step
            integer(kind=8)          :: cval, crate, cmax

            real (8), pointer :: dt, timeNow, timeEnd

            real(8) :: execution_realtime
            real(8) :: execution_realtime_per_step, steps_to_finish
            real(8) :: simulation_fraction, seconds_to_completion, time_to_completion
            real(8) :: thistime

            character(8) :: timeunit
            character(64) :: subroutine_name = 'tl_command_line_step_output'
        !%------------------------------------------------------------------
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
            dt            => setting%Time%Hydraulics%Dt
            timeNow       => setting%Time%Now
            timeEnd       => setting%Time%End
            step          => setting%Time%Step
            interval      => setting%Output%CommandLine%interval
        !%------------------------------------------------------------------
        if (this_image() == 1) then
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%Now = cval

            ! --- estimate the remaining time based on total time from start
            execution_realtime = real((setting%Time%WallClock%Now - setting%Time%WallClock%TimeMarchStart),kind=8) &
                             /  (real(setting%Time%WallClock%CountRate,kind=8)) 
            execution_realtime_per_step = execution_realtime / real(1+step,kind=8)
        
            steps_to_finish       = (setting%Time%End - setting%Time%Now) / dt
            seconds_to_completion = execution_realtime_per_step * steps_to_finish
        end if

        if (setting%Output%Verbose) then
            if (this_image() == 1) then
                if (mod(step,interval) == 0) then
                    thistime = timeNow
                    call util_datetime_display_time (thistime, timeunit)

                    ! write a time counter
                    if (.not. inSpinUpYN) then
                        if (dt > oneR) then
                            write(*,"(A12,i8,a17,F9.2,a1,a8,a6,f9.2,a3,a8,f9.2,a11,f9.2,a13,f9.2)") &
                                'time step = ',step,'; model time = ',thistime, &
                                ' ',trim(timeunit),'; dt = ',dt,' s', '; cfl = ',cfl_max
                        else
                            write(*,"(A12,i8,a17,F9.4,a1,a8,a6,f9.8,a3,a8,f9.2,a11,f9.2,a13,f9.2)") &
                                'time step = ',step,'; model time = ',thistime, &
                                ' ',trim(timeunit),'; dt = ',dt,' s', '; cfl = ',cfl_max
                        end if
                    else
                        write(*,"(A15,i8,a17,f9.2,a1,a8,a6,f9.2,a3,a8,f9.2)") &
                            'spin-up step = ',step,'; model time = ',thistime, &
                          ' ',trim(timeunit),'; dt = ',dt,' s', '; cfl = ',cfl_max 
                    end if
                    if (.not. inSpinUpYN) then
                        ! write estimate of time remaining
                        thistime = seconds_to_completion
                        call util_datetime_display_time (thistime, timeunit)
                        write(*,"(A9,F10.2,A1,A3,A)") 'estimate ',thistime,' ',timeunit,' wall clock time until completion'
                        !write(*,"(A9,F6.2,A1,A3,A)") 'execution time ',thistime,' ',timeunit,' wall clock time thus far'
                    end if    
                    print *, ' '
                endif
            endif
        endif

        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine tl_command_line_step_output
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function tl_get_max_CFL(thisCol,dt) result (outvalue)
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
            integer(8), pointer :: step

            real(8), parameter :: smallDepth = 0.0001
        !%-------------------------------------------------------------------
            Npack              => npack_elemP(thisCol)
            thisP              => elemP(1:Npack,thisCol)
            velocity           => elemR(:,er_Velocity)
            wavespeed          => elemR(:,er_WaveSpeed)
            length             => elemR(:,er_Length)
            PCelerity          => elemR(:,er_Preissmann_Celerity)
        !%------------------------------------------------------------------- 
        step => setting%Time%Step
        if (dt .le. zeroR) then 
            thisDT => setting%Time%Hydraulics%Dt
        else    
            thisDT => dt
        end if

        !% --- Remove wavespeed from JM for CFL computation 
        !%     Note that storage is returned at end of this procedure
        elemR(:,er_Temp01) = wavespeed
        elemR(:,er_Temp02) = PCelerity
        where (elemI(thisP,ei_elementType) .eq. JM)
            wavespeed(thisP) = zeroR
            Pcelerity(thisP) = zeroR
        endwhere
    
        !% --- set the outvalue
        if (Npack > 0) then 
            if (setting%Solver%PreissmannSlot%useSlotTF) then
                !% --- choose between maximum of the advective+wavespeed CFL or the 
                !%     Preissmann Slot celerity
                outvalue = max (maxval((abs(velocity(thisP)) + abs(wavespeed(thisP))) * thisDT / length(thisP)), &
                                maxval((abs(velocity(thisP)) + abs(PCelerity(thisP))) * thisDT / length(thisP)))
                !% DEBUGGING OUTPUT                
                ! print *, 'CFL options ', &
                !     maxval((abs(velocity(thisP)) + abs(wavespeed(thisP))) * thisDT / length(thisP)), &
                !     maxval((abs(PCelerity(thisP))) * thisDT / length(thisP))              
            else 
                outvalue = maxval((abs(velocity(thisP)) + abs(wavespeed(thisP))) * thisDT / length(thisP))
            end if
            !% DEBUGGING OUTPUT
            ! print *, ' '
            ! write(*,"(A,i6,A,f8.4,A)") 'step = ',step,'  oldDT = ',thisDT,'      U+Wc     U+Pc'
            ! write(*,"(A, 3f11.3)") 'CFL using old DT                ', &     
            !                 maxval((abs(velocity(thisP))+abs(wavespeed(thisP))) * thisDT / length(thisP)), &
            !                 maxval((abs(velocity(thisP))+abs(PCelerity(thisP))) * thisDT / length(thisP))  
            ! write(*,"(A, 3f11.3)") 'max vel+wave/L, vel+Pcelerity/L ', &
            !                 maxval((abs(velocity(thisP))+abs(wavespeed(thisP))) * thisDT / length(thisP)), &
            !                 maxval((abs(velocity(thisP))+abs(PCelerity(thisP))) * thisDT / length(thisP))  

            ! write(*,"(A, 3i11)") 'max loc                     ', &
            !                 thisP(maxloc( (abs(velocity(thisP))+abs(wavespeed(thisP))) / length(thisP) )), &
            !                 thisP(maxloc( (abs(velocity(thisP))+abs(PCelerity(thisP))) / length(thisP) ))     
                            
            ! write(*,"(A, 3f11.3)") 'length at max loc               ', &
            !                 length(thisP(maxloc( (abs(velocity(thisP))+abs(wavespeed(thisP))) / length(thisP) ))), &
            !                 length(thisP(maxloc( (abs(velocity(thisP))+abs(PCelerity(thisP))) / length(thisP) )))     
            ! write(*,"(A, 3f11.3)") 'velocity at max loc             ', &
            !                 velocity(thisP(maxloc( (abs(velocity(thisP))+abs(wavespeed(thisP))) / length(thisP) ))), &
            !                 velocity(thisP(maxloc( (abs(velocity(thisP))+abs(PCelerity(thisP))) / length(thisP) )))  
            ! write(*,"(A, 3f11.3)") 'wavespeed at max loc            ', &
            !                 wavespeed(thisP(maxloc( (abs(velocity(thisP))+abs(wavespeed(thisP))) / length(thisP) ))), &
            !                 wavespeed(thisP(maxloc( (abs(velocity(thisP))+abs(PCelerity(thisP))) / length(thisP) )))                  
            ! write(*,"(A, 3f11.3)") 'Pcelerity at max loc            ', &
            !                 PCelerity(thisP(maxloc( (abs(velocity(thisP))+abs(wavespeed(thisP))) / length(thisP) ))), &
            !                 PCelerity(thisP(maxloc( (abs(velocity(thisP))+abs(PCelerity(thisP))) / length(thisP) )))   
                                       
            ! print *, ' '      
        else
            outvalue = zeroR
        end if

        !% --- returning JM wavespeed to regular storage
        where (elemI(thisP,ei_elementType) .eq. JM)
            wavespeed(thisP) = elemR(thisP,er_Temp01)
            Pcelerity(thisP) = elemR(thisP,er_Temp02)
        endwhere
        elemR(:,er_Temp01) = zeroR
        elemR(:,er_Temp02) = zeroR

    end function tl_get_max_CFL   
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
        
        !% --- set the timestep limiter for upstream inflow boundary conditions (BCup)
        if (size(BC%P%BCup) > 0) then

            !% --- ensure flowrate used for limiter is  positive and non-zero
            BCQ(:) = max( abs(BC%flowR(BC%P%BCup,br_value)), setting%Eps%Velocity)

            !% --- store the element indexes downstream of a BCup face
            BCelemIdx =  faceI(BC%flowI(BC%P%BCup,bi_face_idx), fi_Melem_dL)

            !% --- store the scaling depth limit defined by when the induced velocity
            !%     from an inflow is similar to the gravity wave speed of the BC inflow
            depthScale = ( (alpha**4) * (BCQ**2) / (gravity * (topwidth(BCelemIdx)**2) ) )**onethirdR

            !% --- get the depth limit for maximum depth that the time step will be limited as 
            !%     either the depth scale or the specified smallDepth limit
            depthLimit = max(smallDepth, depthScale)

            DTlimit = length(BCelemIdx) * ( ( alpha * topwidth(BCelemIdx) / (gravity * BCQ) )**onethirdR) 
        
            !% --- where the depth is greater than the depthlimit the DT inflow limiter
            !%     is not needed, and we can use the existing DT value
            depthLimit = depth(BCelemIdx) - depthLimit
            where (depthLimit .ge. zeroR)
                DTlimit = thisDT
            endwhere

            !% --- get the smallest DT in the limiter array
            newDTlimit = minval(DTlimit)

            !% --- use the smaller value of the new limit or the input
            thisDT = min(newDTlimit,thisDT) 

            !% --- return to null value storage
            temp_BCupI(:,:) = nullValueI
            temp_BCupR(:,:) = nullValueR
        else
            !% --- continue -- no DT limitation if there are no BCup faces
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
            !% --- face index
            fidx     => BC%headI(ii,bi_face_idx)
            !% --- face flowrate
            Qface    => faceR(fidx,fr_Flowrate)
            !% --- upstream element index
            elemUp   => faceI(fidx, fi_Melem_uL)
            !% --- topwidth of upstream element
            topWidth => elemR(elemUp,er_TopWidth)
            !% --- length of upstream element
            length   => elemR(elemUp,er_Length)
            !% --- depth of upstream element
            depth    => elemR(elemUp,er_Depth)

            if (Qface .ge. zeroR) exit  !% --- no dt limit if this is an outflow

            !% --- scaling depth limit defined by when the induced velocity
            !%     from an inflow is similar to the gravity wave speed of the BC inflow
            depthScale = ( (alpha**4) * (Qface**2) / (gravity * (topWidth**2) ) )**onethirdR

            !% --- get the depth limit for maximum depth that the time step will be limited as 
            !%     either the depth scale or the specified smallDepth limit
            depthLimit = max(smallDepth, depthScale)
            
            !% --- time step limit based on inflow flux
            DTlimit = length * ( ( alpha * topwidth / (gravity * (-Qface)) )**onethirdR) 

            DTlimit = min(DTlimit, thisDT)

            !% --- where the depth is greater than the depthlimit the DT inflow limiter
            !%     is not needed, and we can use the existing DT value
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
            real(8), pointer :: Qlat(:), PlanArea(:), HeightRise(:) 
            real(8), pointer :: DTlimit(:) 
            real(8)          :: newDTlimit

            integer, pointer :: thisP(:), Npack
            integer :: thisCol, ii
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
            Qlat       => elemR(:,er_Temp01)
            PlanArea   => elemR(:,er_Temp02)
            HeightRise => elemR(:,er_Temp03)
            DTLimit    => elemR(:,er_Temp04)
        !%------------------------------------------------------------------
        !% ---- use the pack for CC and JM with H time march 
        !%      (no lateral flows into JB or diagnostic)
        thisCol =  ep_CCJM_H
        Npack   => npack_elemP(thisCol)   

        if (Npack < 1) return
        thisP => elemP(1:Npack,thisCol)

        !% --- initialize temporary arrays
        PlanArea   = zeroR
        HeightRise = zeroR
        Qlat       = zeroR
        DTLimit    = huge(zeroR)

        !% --- store positive flowrate for limiter
        Qlat(thisP) = abs( elemR(thisP,er_FlowrateLateral)) 

        !% --- get the geometry associated with a simple inflow
        where (Qlat > zeroR) 
            PlanArea   = elemR(:,er_Length) * elemR(:,er_BreadthMax)
            HeightRise = Qlat * thisDT / PlanArea
        endwhere

        !% --- if th depth is small and the heighrise is large, then limit the time step
        !%     to a height rise of twice the local small depth value
        where ((HeightRise > setting%SmallDepth%LateralInflowSmallDepth)  &
                .and.                                                     &
                (elemR(:,er_Depth) < fourR * setting%SmallDepth%LateralInflowSmallDepth))

            DTlimit = twoR * setting%SmallDepth%LateralInflowSmallDepth * PlanArea / Qlat
        endwhere

        newDTlimit = minval(DTlimit)

        thisDT = min(thisDT, newDTlimit)

        !%------------------------------------------------------------------
        !% Closing
            elemR(:,er_Temp01) = nullvalueR
            elemR(:,er_Temp02) = nullvalueR
            elemR(:,er_Temp03) = nullvalueR
            elemR(:,er_Temp04) = nullvalueR

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
    subroutine tl_DT_hydrology_match (newDT, oldDT, oldCFL, neededSteps)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes new hydraualic dt with a goal of matching the end of
        !% the next hydrology step where possible.
        !% Also computes the number of steps needed to reach Hydrology Step
        !% neededSteps >= 3 is treated as a large number 
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), intent(inout) :: newDT
            real(8), intent(in)    :: oldDT, oldCFL
            integer, intent(inout) :: neededSteps

            real(8), pointer :: nextHydrologyTime, LastHydraulicsTime
            real(8), pointer :: CFL_hi, dtTol
            real(8) :: timeLeft, timeLeftCFL
        !%------------------------------------------------------------------
        !% Aliases
            CFL_hi             => setting%VariableDT%CFL_hi
            dtTol              => setting%Time%DtTol
            nextHydrologyTime  => setting%Time%Hydrology%NextTime
            lastHydraulicsTime => setting%Time%Hydraulics%LastTime 
        !%------------------------------------------------------------------

        !% --- for combined hydrology and hydraulics 
        !%     we would like to take a single step, or uniform steps towards
        !%     the next hydrology time
        timeLeft = nextHydrologyTime - lastHydraulicsTime

        !% --- small time left will simply step beyond the nextHydrologyTime
        !%     this ensures we don't try to match micro-seconds
        !%     in the hydraulics/hydrology time stepping
        if (timeLeft > dtTol) then
            !% --- get the CFL of a single step using the present time step
            !%     argument determines the packed set of elements
            !%     used for the CFL computation
            timeLeftCFL = tl_get_max_CFL(ep_CCJM_NOTzerodepth,timeleft)

            !% --- time step to hydrology is less than max CFL
            !%     implies one or a few steps required
            if (timeLeftCFL .le. CFL_hi) then 
                call tl_DT_close_to_hydrology_time &
                    (newDT, neededSteps, timeLeft, timeLeftCFL, oldDT, oldCFL)
            else 
                !% --- cannot reach hydrology time step without
                !%     violating CFL_hi
                call tl_DT_standard (newDT, oldDT, oldCFL)
                neededSteps = ceiling(timeLeft/newDT)
                !% --- reset time if a few steps will get to the Hydrology time
                if (neededSteps < 3) then 
                    newDT = timeLeft / real(neededSteps,8)
                end if
            end if
        else 
            !% --- small time remaining to hydrology, so treat as
            !%     next step
            call tl_DT_standard (newDT, oldDT, oldCFL)
            neededSteps = 3
        end if

  end subroutine tl_DT_hydrology_match
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_DT_close_to_hydrology_time &
        (newDT, neededSteps, timeLeft, timeLeftCFL, oldDT, oldCFL)
        !%------------------------------------------------------------------
        !% Description
        !% When the next time step is close to the next hydrology time, this
        !% computes the one or few time steps needed that are consistent
        !% with CFL limitation
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(inout) :: neededSteps
            real(8), intent(inout) :: newDT
            real(8), intent(in)    :: timeLeft, timeLeftCFL, oldDT, oldCFL
            real(8), pointer       :: increaseFactor
        !%------------------------------------------------------------------
        !% Aliases
            increaseFactor     => setting%VariableDT%increaseFactor
        !%------------------------------------------------------------------
        if (timeLeftCFL < oldCFL*increaseFactor) then
            !% --- reach hydrology time in 1 step
            newDT = timeLeft
            neededSteps = oneI
        else 
            !% --- need multiple steps to reach hydrology
            !%     maximum allowed time step
            newDT = oldDT*increaseFactor
            !% --- number of steps to hydrology
            neededSteps = ceiling(timeLeft/newDT)
            newDT = timeLeft / real(neededSteps,8)
        end if

    end subroutine tl_DT_close_to_hydrology_time
!% 
!%==========================================================================
!%==========================================================================
!%   
    subroutine tl_DT_standard (newDT, oldDT, oldCFL)
        !%------------------------------------------------------------------
        !% Description
        !% Standard approach to setting time step
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), intent(inout) :: newDT
            real(8), intent(in)    :: oldDT, oldCFL

            real(8), pointer :: CFL_hi, targetCFL, CFL_lo
            real(8), pointer :: increaseFactor, decreaseFactor

            integer, pointer :: checkStepInterval

            integer(kind=8), pointer :: stepNow, lastCheckStep

            real(8) :: minCFL
        !%------------------------------------------------------------------
            CFL_hi             => setting%VariableDT%CFL_hi
            targetCFL          => setting%VariableDT%CFL_target
            CFL_lo             => setting%VariableDT%CFL_lo

            increaseFactor     => setting%VariableDT%increaseFactor
            decreaseFactor     => setting%VariableDT%decreaseFactor
            checkStepInterval  => setting%VariableDT%NstepsForCheck

            stepNow            => setting%Time%Hydraulics%Step
            lastCheckStep      => setting%VariableDT%LastCheckStep
        !%------------------------------------------------------------------

        !% --- set the minimum CFL, used to detect near zero flow conditions
        if (setting%Limiter%Dt%UseLimitMaxYN) then
            minCFL = setting%Eps%Velocity * oldDT / setting%Discretization%NominalElemLength
        end if

        if (oldCFL .ge. CFL_hi) then 
            !% --- always reduce large CFL to target
            newDT = oldDT * targetCFL / oldCFL
            lastCheckStep = stepNow
        else
            !% --- CFL is below CFL_hi, only modify on check intervals
            !%     and if value is less than CFL_lo or between targetCFL
            !%     and CFL_hi (i.e., we leave time step if between 
            !%     target and CFL_lo)
            if (stepNow >= lastCheckStep + checkStepInterval) then
                if (oldCFL .le. minCFL) then 
                    !% --- vanishing CFL
                    !%     For a really small CFL, the new DT could be unreasonably large 
                    !%     (or infinite) so use a value based on inflows to fill to
                    !%     small volume
                    call tl_dt_vanishingCFL (newDT)
                else
                    !% --- increase DT if below CFL_lo limit
                    if ((oldCFL <  CFL_lo)) then 
                        newDT = min(OldDT * targetCFL / oldCFL, oldDT*increaseFactor)
                    elseif ((oldCFL > targetCFL) .and. oldCFL < CFL_hi) then
                    !% --- decrease DT if above target
                        newDT = max(OldDT * targetCFL / oldCFL, OldDT*decreaseFactor )
                    else
                        !% --- CFL is between CFL_lo and targetCFL
                    end if
                end if
                lastCheckStep = stepNow
            else
                !% --- no change in DT
            end if
        endif

        !% --- check for minimum limit
        if (setting%Limiter%Dt%UseLimitMinYN) then
            if (newDT < setting%Limiter%Dt%Minimum) then 
                print *, 'WARNING time step is set to minimum ', &
                     setting%Limiter%Dt%Minimum, ' seconds'
                newDT = setting%Limiter%Dt%Minimum
            endif 
        end if

        !% --- check for maximum
        if (setting%Limiter%Dt%UseLimitMaxYN) then 
            if (newDT > setting%Limiter%Dt%Maximum) then 
                print *, 'WARNING time step set to maximum ', &
                    setting%Limiter%Dt%Maximum, ' seconds'
                newDT = setting%Limiter%Dt%Maximum
            end if
        end if

    end subroutine tl_DT_standard
!% 
!%==========================================================================
!%==========================================================================
!%       
    subroutine tl_limit_DT (newDT)
        !%------------------------------------------------------------------
        !% Description:
        !% Apply various limiters to time step
        !% This is called after dt is set by hydraulic flows in channels
        !% and junctions to account for other destabilizing conditions
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), intent(inout) :: newDT

            real(8), pointer :: reportDT
        !%------------------------------------------------------------------
        !% Aliases
            reportDt => setting%Output%Report%TimeInterval
        !%------------------------------------------------------------------

        !% --- limit dt for reporting
        !%     HACK: limit the newDT to be equal or smaller than the reportDT
        !%     this ensures we get consistant output at every reporting steps
        !%     however, this will slow down the timeloop if reportDT causes
        !%     lower cfl values. Needs to revisit later. 
        newDT = min(newDT,reportDt)
    
        !% --- time step limiter for inflows into small or zero volumes
        call tl_limit_BCinflow_dt (newDT)

        call tl_limit_BChead_inflow_dt (newDT)

        call tl_limit_LatInflow_dt (newDT)

        !% --- limit by inflow/head external boundary conditions time intervals
        if (setting%VariableDT%limitByBC_YN) then
            newDT = min(setting%BC%smallestTimeInterval,newDT)
        end if

    end subroutine tl_limit_DT
!% 
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_roundoff_DT (newDT, neededSteps)
        !%------------------------------------------------------------------
        !% Description:
        !% Rounds off dt to reduce the number of significant digits
        !%------------------------------------------------------------------
        !% Declarations:
            real(8), intent(inout) :: newDT 
            integer, intent(in)    :: neededSteps

            logical :: do_roundoff = .false.
        !%------------------------------------------------------------------

        select case (neededSteps)
            case(1,2)
                !% -- accept dt as is
                do_roundoff = .false.
            case default !%-- for 3 steps or higher, always round off
                do_roundoff = .true.
        end select

        if (do_roundoff) then 
            if (newDT > fiveR) then
                !% --- round down larger dt to counting numbers
                newDT = real(floor(newDT),8)
            elseif ((newDT > oneR) .and. (newDT .le. fiveR)) then
                !% --- round down smaller dt to two decimal places
                newDT = real(floor(newDT * tenR),8) / tenR
            elseif ((newDT > onetenthR) .and. (newDT .le. oneR)) then
                !% --- round down smaller dt to three decimal places
                newDT = real(floor(newDT * onehundredR),8) / onehundredR
            else 
                !% do not round smaller values
            end if
        else 
            !% do not round
        end if

    end subroutine tl_roundoff_DT
!% 
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module timeloop

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
              interface_get_newRunoffTime
    use utility_crash
    use control_hydraulics, only: control_update

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
        end if 

        !% --- initialize the time settings for hydraulics and hydrology steps
        call tl_initialize_loop (doHydraulicsStepYN, doHydrologyStepYN)

        !-- perform the time-marching loop
        call tl_outerloop (doHydrologyStepYN, doHydraulicsStepYN, .false., .false.)

        sync all
        !% --- close the timemarch time tick
        if (this_image() == 1) then
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%TimeMarchEnd= cval
        end if

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
        call tl_initialize_loop (doHydraulicsStepYN, doHydrologyStepYN)
        
        !% --- perform the time loop for spin-up
        call tl_outerloop (doHydrologyStepYN, doHydraulicsStepYN, inSpinUpYN, SpinUpOutputYN)

    end subroutine tl_spinup
!% 
!%==========================================================================
!%==========================================================================
!%      
    subroutine tl_initialize_loop (doHydraulicsStepYN, doHydrologyStepYN)
        !%------------------------------------------------------------------
        !% Description
        !% initialize the times before a time loop
        !%------------------------------------------------------------------
        !% Declarations
            logical, intent(inout) :: doHydrologyStepYN, doHydraulicsStepYN
            real(8), pointer :: nextHydrologyTime, nextHydraulicsTime, nextControlRuleTime
            real(8), pointer :: lastHydrologyTime, lastHydraulicsTime, lastControlRuleTime, dtTol
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
        else   
            !% set a large dummy time for hydrology if not used
            nextHydrologyTime  = setting%Time%End + onethousandR * dtTol
        end if

        !% get the initial dt and the next hydraulics time
        if (setting%Simulation%useHydraulics) then
            call tl_smallestBC_timeInterval ()
            call tl_update_hydraulics_timestep()
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
        nextControlRuleTime = lastControlRuleTime + real(setting%SWMMinput%RuleStep,8)

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
            logical :: BCupdateYN
        !%-------------------------------------------------------------------
        !% Aliases
            thisStep => setting%Time%Step
            dtTol    => setting%Time%DtTol 
        !%-------------------------------------------------------------------
        !% --- initialize the time step counter
        thisStep = 1

        do while (setting%Time%Now <= setting%Time%End - dtTol)
                !% --- set the controls for using spin-up time
                if ((inSpinUpYN) .and. (thisStep > 1)) then
                    !% --- skip BC update during spin-up after first step
                    BCupdateYN = .false.
                else
                    BCupdateYN = .true.
                end if
                    ! print *, ' '
                    ! print *, ' '
                    ! print *, '*******************************************************************************'
                    ! print *, '*******************************************************************************'
                    ! write(6,"(A,f12.5,A,f12.5)") ' ... beginning time loop ======== time (h):',&
                    !    setting%Time%Now/3600.d0, ';  DT (s) =',setting%Time%Hydraulics%Dt
                   
                    ! call util_CLprint ('at start of time loop')
    
                !% --- push the old values down the stack 
                call tl_save_previous_values()
    
                !% store the runoff from hydrology on a hydrology step
                if ((.not. inSpinUpYN) .or. (inSpinUpYN .and. BCupdateYN) ) then
                    if (doHydrologyStepYN) call tl_hydrology()
                end if
    
                !% --- main hydraulics time steop
                if (doHydraulicsStepYN) then    

                    !% --- set a clock tick for hydraulic loop evaluation
                    if ((this_image()==1) .and. (.not. inSpinUpYN)) then
                        call system_clock(count=cval,count_rate=crate,count_max=cmax)
                        setting%Time%WallClock%HydraulicsStart = cval
                    end if 
    
                    !print *, 'about to update BC in timeloop'

                    !% --- get updated boundary conditions
                    if (BCupdateYN) then
                        call bc_update() 
                        call tl_lateral_inflow()
                        call tl_smallestBC_timeInterval ()
                    end if

                    !print *, 'about to perform control rules '

                    !% --- perform control rules
                    if ((.not. inSpinUpYN) .and. (setting%SWMMinput%N_control > 0)) then
                        if (setting%Time%Now .ge. setting%Time%ControlRule%NextTime) then
                            !% --- evaluate all the controls
                            !%     required for all images because monitorI data are spread
                            !%     across all images.
                            call control_update ()
                            !% --- set the next time the controls will be evaluated
                            setting%Time%ControlRule%NextTime = setting%Time%Now + real(setting%SWMMinput%RuleStep,8)
                        end if
                        ! print *, 'CONTROL----------------------------'
                        ! print *, 'orifice setting ',elemR(iet(3),er_Setting)
                    end if

                    !print *, 'about to call tl_subcatchment_lateral_inflow'
    
                    !% --- add subcatchment inflows
                    !%     note, this has "useHydrology" and not "doHydrologyStepYN" because the
                    !%     former is whether or not hydrology is used, and the latter is
                    !%     whether or not it is computed in this time step
                    if (setting%Simulation%useHydrology .and. BCupdateYN) call tl_subcatchment_lateral_inflow ()
    
                    !% --- perform hydraulic routing
                    call tl_hydraulics()

    
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
    
                call util_crashstop(13978)
    
                !% --- restart the hydraulics time tick
                sync all
                if ((this_image()==1) .and. (.not. inSpinUpYN) ) then
                    call system_clock(count=cval,count_rate=crate,count_max=cmax)
                    setting%Time%WallClock%HydraulicsStart = cval
                end if 
    
                !call util_CLprint('before time step change')
        
                sync all
                !% ---increment the time step and counters for the next time loop
                call tl_increment_timestep_and_counters(doHydraulicsStepYN, doHydrologyStepYN)
    
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

                
    
                !% --- check for blowup conditions
                call util_crashcheck (773623)
                if (crashI == 1) exit 
    
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
            real(8), pointer :: sRunoff(:)
            integer(kind=8) :: crate, cmax, cval
            character(64) :: subroutine_name = 'tl_hydrology'
        !%------------------------------------------------------------------
        !% Preliminaries
            !if (crashYN) return
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
        !%------------------------------------------------------------------   
            
        if (setting%Simulation%useSpinUp) then
            print *, 'ERROR setting.simulation.useSpinUp = true is not supported for hydrology at this time'
            call util_crashpoint(4482333)
        end if

        !% execute the SWMM-C runoff for the next interval 
        call interface_call_runoff_execute()

        !% get the next runoff time
        setting%Time%Hydrology%NextTime = interface_get_NewRunoffTime()   

        !% cycle through the subcatchments to get the runoff 
        !% ii-1 required in arg as C arrays start from 0
        do ii = 1,setting%SWMMinput%N_subcatch
            sRunoff(ii) = interface_get_subcatch_runoff(ii-1)  
        end do

        !%------------------------------------------------------------------   
        !% wall clock tick
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
            character(64)    :: subroutine_name = 'tl_hydraulics'
        !%-------------------------------------------------------------------
        !% Preliminaries:
            !if (crashYN) return
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------

        !% ---check for where solver needs to switch in dual-solver model
        if (setting%Solver%SolverSelect == ETM_AC) call tl_solver_select()

        !% --- repack all the dynamic arrays
        call pack_dynamic_arrays()

        !% --- ensure that the conservative flux terms are exactly zero in the entire array
        !%     so that we can be confident of conservation computation. 
        faceR(:,fr_Flowrate_Conservative) = zeroR  

        !% call the RK2 time-march, depending on the type of solver
        select case (setting%Solver%SolverSelect)
        case (ETM_AC)
            !% dual-method, ETM for open channel, AC for surcharge
            call rk2_toplevel_ETMAC() 
        case (ETM)
        !    print *, ' '
        !    print *, 'here calling rk2 toplevel------------------------------------------',this_image()
        !    call util_CLprint ()

                !outstring = '    tl_hydraulics 111 '
                !call util_syncwrite

            !% ETM with Preissmann slot for surcharge
            call rk2_toplevel_ETM()
            
                !outstring = '    tl_hydraulics 222 '
                !call util_syncwrite

            !call util_CLprint ()
        case (AC)
            !% AC for both open-channel and surcharge
            call rk2_toplevel_AC()
        case DEFAULT
            print *, 'CODE ERROR setting.Solver.SolverSelect type unknown for # ', setting%Solver%SolverSelect
            print *, 'which has key ',trim(reverseKey(setting%Solver%SolverSelect))
            stop 497895
        end select    

        call util_accumulate_volume_conservation () 

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
            integer, pointer :: npack, thisP(:), thisBC(:)
            real(8), pointer :: Qlateral(:)
        !%------------------------------------------------------------------
        !% Aliases:
            Qlateral => elemR(:,er_FlowrateLateral)
            npack    => npack_elemP(ep_BClat)
            thisP    => elemP(1:npack,ep_BClat)
            thisBC   => BC%P%BClat
        !%------------------------------------------------------------------
        !% set lateral to zero
        Qlateral(:) = zeroR 

        !% --- add lateral inflow BC to lateral inflow accumulator
        !%     note that thisP and thisBC must be the same size or there is something wrong       
        Qlateral(thisP) = Qlateral(thisP) + BC%flowR(thisBC,br_value) 

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
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: sImage(:), eIdx(:)
            real(8), pointer :: Qrate(:), Qlateral(:)
            integer :: mm
        !%------------------------------------------------------------------
        !% Aliases
            sImage => subcatchI(:,si_runoff_P_image)
            eIdx   => subcatchI(:,si_runoff_elemIdx)
            !% --- using the full runoff rate for this period
            Qrate => subcatchR(:,sr_RunoffRate_baseline)
            Qlateral => elemR(:,er_FlowrateLateral)
        !%------------------------------------------------------------------
        do mm = 1,setting%SWMMinput%N_subcatch
            !% --- only if this image holds this node
            !print *, mm, eIdx(mm), Qlateral(eIdx(mm)), Qrate(mm)
            if (this_image() .eq. sImage(mm)) then
                Qlateral(eIdx(mm)) = Qlateral(eIdx(mm)) + Qrate(mm)
            end if
        end do

    end subroutine tl_subcatchment_lateral_inflow
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_increment_timestep_and_counters(doHydraulicsStepYN, doHydrologyStepYN)
        !%-------------------------------------------------------------------
        !% Description:
        !% increments the hydrology and hydraulics step counters and
        !%-------------------------------------------------------------------
        !% Declarations
            logical, intent(inout) :: doHydraulicsStepYN, doHydrologyStepYN
            logical, pointer       :: useHydrology, useHydraulics
            real(8), pointer       :: nextHydraulicsTime, nextHydrologyTime
            real(8), pointer       :: lastHydraulicsTime, lastHydrologyTime
            real(8)                :: nextTime, nextTimeLocal
            real(8), pointer       :: dtTol, timeNow, dtHydraulics
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
            call tl_update_hydraulics_timestep()
            call util_crashstop(449873)
        else
            nextHydraulicsTime = setting%Time%End + tenR*DtTol
        end if    

        ! print *, 'DT            ',dtHydraulics 
        ! print *, 'next Hydraulics Time AAA ',nextHydraulicsTime

        !% --- The NextHydrologyTime is updated in SWMM, here we just need to
        !%     provide a large number if hydrology isn't used
        if (.not. useHydrology) then
            nextHydrologyTime = max(setting%Time%End + tenR*DtTol, &
                                    setting%Simulation%SpinUpDays * seconds_per_day)
        else
            !% --- SWMM-C will not return a next hydrology time that is
            !%     beyond the end of the simulation, so we need a special
            !%     treatment for this case
            if ((nextHydrologyTime .eq. lastHydrologyTime) .and. &
                (timeNow > setting%Time%Start)) then 
                nextHydrologyTime = max(setting%Time%End + tenR*DtTol , &
                                        setting%Simulation%SpinUpDays * seconds_per_day)
            end if
        end if

        ! print *, 'next Hydraulics Time BBB',nextHydraulicsTime

        !% --- Ensure that all processors use the same time step.
        !% --- find the minimum hydraulics time and store accross all processors
        call co_min(nextHydraulicsTime)
        !% --- take the nextTime as the minimum of either the Hydrology or Hydraulics time
        !%     done on a single processor because they all should have identical nextHydrologyTIme
        nextTime = min(nextHydraulicsTime, nextHydrologyTime)

        !print *, 'next Time CCC ',nextTime

        !% --- note there is no need to broadcast across processors since co_min
        !%     ensures that each processor gets the min value.

        !% --- update the time step for the local processor (all have the same nextTime and lastHydraulicsTime)
        dtHydraulics = nextTime - lastHydraulicsTime

        ! print *, 'dtHydarulics, next time ',dtHydraulics, nextTime

        ! if (dtHydraulics == zeroR) stop 4498723

        doHydrologyStepYN  = (abs(nextTime - nextHydrologyTime)  <= dtTol) .and. useHydrology
        doHydraulicsStepYN = (abs(nextTime - nextHydraulicsTime) <= dtTol) .and. useHydraulics

        !nextTimeLocal = nextTime
        !% Get the smallest nextTime across all the processors, and make that the
        !% nextTime on all processors
        !call co_min(dt)
        !call co_min(nextTime)

        ! if (nextTime == nextTimeLocal) then
        !     minImg = this_image()
        ! else
        !     minImg = -1
        ! end if
        ! call co_max(minImg)
        ! call co_broadcast(doHydraulicsStepYN, minImg)
        ! call co_broadcast(doHydrologyStepYN, minImg)

        if ((.not. inSpinUpYN) .or. & 
            ((inSpinUpYN) .and. (setting%Simulation%stopAfterSpinUp)) ) then
            if (util_output_must_report()) reportStep = reportStep + 1
        end if
        if (doHydraulicsStepYN) hydraulicStep = hydraulicStep + 1
        if (doHydrologyStepYN)  hydrologyStep = hydrologyStep + 1

        step    = step + 1
        timeNow = nextTime !timeNow + dt

        ! print *, 'timeNow  in tl_increment_counters',timeNow

        call tl_command_line_step_output()

        if (doHydraulicsStepYN) LastHydraulicsTime = NextHydraulicsTime
        if (doHydrologyStepYN)  LastHydrologyTime  = NextHydrologyTime

        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        end subroutine tl_increment_timestep_and_counters
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_update_hydraulics_timestep()
        !%------------------------------------------------------------------
        !% Description:
        !% updates the timestep (dt) for hydraulics and computes the
        !% setting.Time.Hydraulics.NextTime for the current processor
        !%------------------------------------------------------------------
        !% Declarations
            logical, pointer :: matchHydrologyStep, useHydrology
            real(8)          :: oldDT, maxVelocity
            real(8)          :: timeleft, thisCFL, minCFL
            real(8), pointer :: targetCFL, maxCFL, maxCFLlow, timeNow, dtTol
            real(8), pointer :: increaseFactor
            real(8), pointer :: newDT
            !rm velocity(:), wavespeed(:), length(:), PCelerity(:)
            real(8), pointer :: nextHydrologyTime, nextHydraulicsTime
            real(8), pointer :: lastHydrologyTime, lastHydraulicsTime
            integer          :: ii, neededSteps
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
        oldDT =  setting%Time%Hydraulics%Dt  ! not a pointer (important!)
        newDT => setting%Time%Hydraulics%Dt

        !% --- set the minimum CFL, used to detect near zero flow conditions
        if (setting%Limiter%Dt%UseLimitMaxYN) then
            minCFL = setting%Eps%Velocity * oldDT / setting%Discretization%NominalElemLength
        end if

        if ((matchHydrologyStep) .and. (useHydrology) .and. (.not. inSpinUpYN)) then 
            !% --- for combined hydrology and hydraulics compute the CFL if we take a single
            !%     step to match the next hydrology time
            timeLeft = nextHydrologyTime - lastHydraulicsTime
            if (timeLeft .le. dtTol) timeLeft = oldDT
            !% --- get the CFL if a single step is taken
            !thisCFL = tl_get_max_cfl(ep_CCJBJM_NOTsmalldepth,timeleft)  
            thisCFL = tl_get_max_cfl(ep_CCJM_NOTsmalldepth,timeleft)  

            !% --- check to see if a single time step to match the hydrology time is possible
            if (thisCFL .le. maxCFL) then
                neededSteps = 1 !% only 1 hydraulic step to match hydrology step
                !% --- check to be sure there is significant time left
                if (timeLeft > dtTol) then
                    !% --- check increase that is implied with a single time step to the next hydrology
                    if (timeleft / oldDT < increaseFactor) then
                        !% --- use a single hydraulics step for the remaining hydrology time
                        newDT = timeleft
                    else
                        !% --- increase the oldDT by the maximum allowed
                        newDT = increaseFactor * oldDT
                    end if
                else
                    !% --- negligible time left, don't bother with it, go back to the old dt
                    newDT = oldDT
                    !% --- check that resetting to oldDT didn't cause a problem
                    !thisCFL = tl_get_max_cfl(ep_CCJBJM_NOTsmalldepth,newDT)
                    thisCFL = tl_get_max_cfl(ep_CCJM_NOTsmalldepth,newDT)
                    if (thisCFL > maxCFL) then 
                        !% --- if CFL too large, set the time step based on the target CFL
                        newDT = newDT * targetCFL / thisCFL
                    end if

                end if
            else
                !% --- if more than one time step is needed, compute the needed steps 
                !%     to break up the large CFL and time step size.
                if (thisCFL/targetCFL .ge. huge(neededSteps)) then
                    !% --- check to see if the implied time step is too small for the integer size
                    !%     this is likely a blow-up condition
                    write(*,*) 'warning -- really high velocity, setting dt to minimum to prevent overflow'
                    newDT = setting%Limiter%Dt%Minimum + setting%Time%DtTol
                    neededSteps = 1000
                else
                    !% ---  for a moderate case where multiple hydraulic steps are needed for a hydrology step
                    if (thisCFL > minCFL) then
                        !% --- note that neededSteps will be 2 or larger else thisCFL < maxCFL
                        neededSteps = ceiling( thisCFL / targetCFL )
                        !% --- the provisional time step that would get exactly to the hydrology time (if CFL didn't change)
                        newDT = timeleft / real(neededSteps,8)
                        !% --- limit the change in the dt by the increase factor
                        if (newDT / oldDT > increaseFactor) then 
                            newDT = oldDT * increaseFactor
                        else
                            !% --- accept the provisional newDT
                        end if
                        !% --- check that the newDT didn't cause a CFL violation
                        !thisCFL = tl_get_max_cfl(ep_CCJBJM_NOTsmalldepth,newDT)
                        thisCFL = tl_get_max_cfl(ep_CCJM_NOTsmalldepth,newDT)
                        if (thisCFL > maxCFL) then 
                            !% --- if CFL to large, set the time step based on the target CFL
                            newDT = newDT * targetCFL / thisCFL
                        end if
                    else
                        !% --- Miniscule CFL is typically due to small volumes in the start-up condition
                        !%    we can use a larger time step based on volumes to fill
                        call tl_dt_vanishingCFL(newDT)
                    end if

                end if
            end if
        else
            !% --- for hydraulics without hydrology or if we are not matching the hydrology step
            neededSteps = 3 !% forces rounding check

            !% --- allowing hydrology and hydraulics to occur at different times
            !thisCFL = tl_get_max_cfl(ep_CC_NOTsmalldepth,oldDT)
            !thisCFL = tl_get_max_cfl(ep_CCJBJM_NOTsmalldepth,oldDT)
            thisCFL = tl_get_max_cfl(ep_CCJM_NOTsmalldepth,oldDT)

                !print *, 'baseline CFL, minCFL, this step: '
                !print *, thisCFL, minCFL, stepNow

            if (thisCFL .ge. maxCFL) then
                !% --- always decrease DT for exceeding the max CFL
                !%     new value is based on the target CFL
                newDT = oldDT * targetCFL / thisCFL
                lastCheckStep = stepNow

                !print *, 'Adjust DT for high CFL    ',newDT   

            else 
                !% --- if CFL is less than max, see if it can be raised (only checked at intervals)
                if (stepNow >= lastCheckStep + checkStepInterval) then
                   ! print *, '------------------------------------------------------------------------------'
                   ! print *, 'checking step: ',stepNow ,lastCheckStep, checkStepInterval

                    !% --- check for low CFL only on prescribed intervals and increase time step
                    if (thisCFL .le. minCFL) then
                        !% --- for really small CFL, the new DT could be unreasonably large (or infinite)
                        !%     so use a value based on inflows to fill to small volume
                        !print *, 'vanishing CFL'
                        call tl_dt_vanishingCFL(newDT)

                    elseif ((minCFL < thisCFL) .and. (thisCFL .le. maxCFLlow)) then
                        !% --- increase the time step and reset the checkStep Counter
                       ! print *, 'lowCFL'
                        newDT = OldDT * targetCFL / thisCFL 
                    else
                        !% -- for maxCFLlow < thisCFL < maxCFL do nothing
                    end if
                    lastCheckStep = stepNow

                end if
            end if
        end if

        !% --- prevent large increases in the time step
        newDT = min(newDT,OldDT * increaseFactor)
           ! print *, 'before adjust newDT, oldDT : ',newDT, OldDT
            ! print *, ' '

        !% 20220328brh time step limiter for inflows into small or zero volumes
            ! print *, 'dt before limit ', newDT
        call tl_limit_BCinflow_dt (newDT)
            !print *, 'dt BC flow limit     ',newDT

        call tl_limit_BChead_inflow_dt (newDT)
            !print *, 'dt BC head limit     ',newDT

        call tl_limit_LatInflow_dt (newDT)
            !print *, 'dt Qlat limit   ',newDt

        !% --- limit by inflow/head external boundary conditions time intervals
        if (setting%VariableDT%limitByBC_YN) then
            !print *, 'here limiting DT by BC ',newDT, setting%BC%smallestTimeInterval
            newDT = min(setting%BC%smallestTimeInterval,newDT)
        end if

        !% --- if dt is large and there is more than 2 steps, then round to an integer number
        if ((matchHydrologyStep) .and. (useHydrology) .and. (neededSteps .le. 2)) then
            !% don't round the dt
        else
            if (newDT > fiveR) then
                !% --- round larger dt to counting numbers
                newDT = real(floor(newDT),8)

                    !print *, 'rounding large    ',newDT
            elseif (newDT > oneR) then
                !% --- round smaller dt to two decimal places
                newDT = real(floor(newDT * onehundredR),8) / onehundredR

                    !print *, 'rounding small    ',newDT   
            else
                !%  HACK -- should round to 3 places for smaller numbers
            end if
        end if

         !print *, 'final DT: ',newDT
        !  print *, ' '

        !% increment the hydraulics time clock
        nextHydraulicsTime = lastHydraulicsTime + newDT
        
        !% find the cfl for reporting
       ! cfl_max = tl_get_max_cfl(ep_CCJBJM_NOTsmalldepth,newDT)
        cfl_max = tl_get_max_cfl(ep_CCJM_NOTsmalldepth,newDT)
        call co_max(cfl_max)


        !stop 29873

        !%----------------------------------------------------------------------
        !% closing
            if ((setting%Limiter%Dt%UseLimitMinYN) .and. (newDT .le. setting%Limiter%Dt%Minimum)) then
                print*, 'timeNow = ', timeNow
                print*, 'dt = ', newDT, 'minDt = ',  setting%Limiter%Dt%Minimum
                print*, 'max velocity  ', maxval( &
                    elemR(elemP(1:npack_elemP(ep_CC_NOTsmalldepth),ep_CC_NOTsmalldepth),er_Velocity) )
                print*, 'max wavespeed ', maxval( &
                    elemR(elemP(1:npack_elemP(ep_CC_NOTsmalldepth),ep_CC_NOTsmalldepth),er_WaveSpeed) )
                print*, 'warning: the dt value is smaller than the user supplied min dt value'
                !stop 1123938
                call util_crashpoint(1123938)
            end if

            ! print*, 'timeNow = ', timeNow
            ! print*, 'dt = ', newDT, 'minDt = ',  setting%Limiter%Dt%Minimum
            ! maxVelocity = maxval( &
            !     elemR(elemP(1:npack_elemP(ep_CC_NOTsmalldepth),ep_CC_NOTsmalldepth),er_Velocity) )
            ! print*, 'max velocity  ', maxVelocity
            ! print*, 'max V locate  ', findloc(elemR(:,er_Velocity),maxVelocity)
            ! print*, 'max wavespeed ', maxval( &
            !     elemR(elemP(1:npack_elemP(ep_CC_NOTsmalldepth),ep_CC_NOTsmalldepth),er_WaveSpeed) )

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
            ncol        => col_elemP(ep_CC_Alltm) !% cannot use JM elements
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
            Qelem       => elemR(:,er_FlowrateLateral)
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
    subroutine tl_solver_select()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% For ETM_AC dual method, this sets the elemI(:,ei_tmType) to the type of solver
        !% needed depending on the volume and the volume cutoffs.
        !% Should only be called if setting%Solver%SolverSelect == ETM_AC
        !%-----------------------------------------------------------------------------
        integer :: thisCol
        integer, pointer :: Npack, tmType(:), thisP(:)
        real(8), pointer :: sfup, sfdn
        real(8), pointer :: volume(:), FullVolume(:)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'tl_solver_select'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        thiscol = ep_ALLtm
        Npack => npack_elemP(thisCol)
        thisP => elemP(1:Npack,thisCol)

        tmType     => elemI(:,ei_tmType)
        volume     => elemR(:,er_Volume)
        FullVolume => elemR(:,er_FullVolume)

        sfup => setting%Solver%SwitchFractionUp
        sfdn => setting%Solver%SwitchFractionDn
        !%-----------------------------------------------------------------------------
        !% Look for ETM elements that are above the cutoff for going to AC and set
        !% these to AC
        where ( ( (volume(thisP) / FullVolume(thisP) ) > sfup ) .and. (tmType(thisP) == ETM) )
            tmType(thisP) = AC
        endwhere

        !% Look for AC elements that are below the cutoff for going back to ETM and
        !% set these to ETM
        where ( ( (volume(thisP) / FullVolume(thisP) ) < sfdn) .and. (tmType(thisP) == AC) )
            tmType(thisP) = ETM
        endwhere

        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine tl_solver_select
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

        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine tl_save_previous_values
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_command_line_step_output ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !%
        !%-----------------------------------------------------------------------------
            character(64) :: subroutine_name = 'tl_command_line_step_output'
            real (8), pointer :: dt, timeNow, timeEnd
            real (8) :: thistime
            integer, pointer :: interval
            integer (kind=8), pointer :: step
            integer (kind=8) :: execution_realtime
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

            ! estimate the remaining time
            execution_realtime = (setting%Time%WallClock%Now - setting%Time%WallClock%TimeMarchStart)
            seconds_to_completion =  (    (real(execution_realtime,kind=8))                     &
                                       /  (real(setting%Time%WallClock%CountRate,kind=8))  )    &
                                   * (    (setting%Time%End - setting%Time%Now)                 &
                                       /  (setting%Time%Now - setting%Time%Start) )
        end if

        if (setting%Output%Verbose) then
            if (this_image() == 1) then
                if (mod(step,interval) == 0) then
                    thistime = timeNow
                    call util_datetime_display_time (thistime, timeunit)

                    ! write a time counter
                    if (.not. inSpinUpYN) then
                        write(*,"(A12,i8,a17,F9.2,a1,a8,a6,f9.2,a3,a8,f9.2,a11,f9.2,a13,f9.2)") &
                            'time step = ',step,'; model time = ',thistime, &
                            ' ',trim(timeunit),'; dt = ',dt,' s', '; cfl = ',cfl_max!, &
                        !   '; cfl_CC = ',cfl_max_CC,'; cfl_JBJM = ',cfl_max_JBJM 
                    else
                        write(*,"(A15,i8,a17,f9.2,a1,a8,a6,f9.2,a3,a8,f9.2)") &
                            'spin-up step = ',step,'; model time = ',thistime, &
                          ' ',trim(timeunit),'; dt = ',dt,' s', '; cfl = ',cfl_max!, &
                        !   '; cfl_CC = ',cfl_max_CC,'; cfl_JBJM = ',cfl_max_JBJM 
                    end if
                    if (.not. inSpinUpYN) then
                        ! write estimate of time remaining
                        thistime = seconds_to_completion
                        call util_datetime_display_time (thistime, timeunit)
                        write(*,"(A9,F6.2,A1,A3,A)") 'estimate ',thistime,' ',timeunit,' wall clock time until completion'
                    end if    
                    print *
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
    real(8) function tl_get_max_cfl(thisCol,dt) result (outvalue)
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
            integer :: itemp(1), ip, fup, fdn, eup, edn
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

        ! print *, ' '
        ! print *, 'in tl_get_max_cfl, Npack = ', Npack
        ! write(*,"(A,10f12.5)") 'Vcfl ' , velocity(iet) * thisDT / length(iet)
        ! write(*,"(A,10f12.5)") 'Hcfl ' , wavespeed(iet) * thisDT / length(iet)
        ! write(*,"(A,10f12.5)") 'Ccfl ' , PCelerity(iet) * thisDT / length(iet)
        ! print * ,' '
        ! print *, 'thisP '
        ! print *, thisP
        ! write(*,"(A,10f12.5)") 'Vcfl ' , velocity(thisP) * thisDT / length(thisP)
        ! write(*,"(A,10f12.5)") 'Hcfl ' , wavespeed(thisP) * thisDT / length(thisP)
        ! print *, thisDT
        ! print *, length(thisP)
        ! print *, wavespeed(thisP)
        !write(*,"(A,10f12.5)") 'Ccfl ' , PCelerity(thisP) * thisDT / length(thisP)

        if (Npack > 0) then 
            outvalue = max (maxval((abs(velocity(thisP)) + abs(wavespeed(thisP))) * thisDT / length(thisP)), &
                            maxval((abs(PCelerity(thisP))) * thisDT / length(thisP)))
           ! print *, 'max vel, wave ', maxval(abs(velocity(thisP))), maxval(abs(wavespeed(thisP)))        
        else
            outvalue = zeroR
        end if
        ! print *, 'outvalue CFL in tl_get_max_cfl',outvalue
        ! print *, ' '

        ! print *, ' '
        ! print *, 'in tl_get_max_cfl'
        ! print *, length(iet(1))
        ! print *, abs(wavespeed(iet(1))) * thisDT / length(iet(1))
        ! print *, outvalue
        ! print *, ' '

    end function tl_get_max_cfl    
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
            smallDepth => setting%SmallDepth%DepthCutoff
        !%------------------------------------------------------------------
        
        !% set the timestep limiter for upstream inflow boundary conditions (BCup)
        ! print *, ' '
        ! print *, 'this dt into limiter ',thisDT
        ! print *, 'size  in tl_limit_BCinflow_dt',size(BC%P%BCup)
        if (size(BC%P%BCup) > 0) then

            !% ensure flowrate used for limiter is  positive and non-zero
            !tinyQ(:) = setting%Eps%Velocity
            BCQ(:) = max( abs(BC%flowR(BC%P%BCup,br_value)), setting%Eps%Velocity)

           ! print *, 'BCQ ',BCQ(:)

            !% store the element indexes downstream of a BCup face
            BCelemIdx =  faceI(BC%flowI(BC%P%BCup,bi_face_idx), fi_Melem_dL)

            !% store the scaling depth limit defined by when the induced velocity
            !% from an inflow is similar to the gravity wave speed of the BC inflow
            depthScale = ( (alpha**4) * (BCQ**2) / (gravity * (topwidth(BCelemIdx)**2) ) )**onethirdR

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

            ! print *, 'DTlimit ',DTlimit
        
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
            smallDepth => setting%SmallDepth%DepthCutoff
        !%------------------------------------------------------------------    

        do ii=1, N_headBC

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
            real(8), pointer :: Qlat(:), depthScale(:), depthLimit(:)
            real(8), pointer :: DTlimit(:), topwidth(:), length(:), depth(:)
            real(8), pointer :: alpha, gravity, smallDepth
            real(8) :: newDTlimit

            integer, pointer :: thisP(:), Npack
            integer :: thisCol
        !%------------------------------------------------------------------
        !% Preliminaries:
            !if (crashYN) return
        !%------------------------------------------------------------------
        !% Aliases:
            Qlat       => elemR(:,er_Temp01)
            depthScale => elemR(:,er_Temp02)
            depthLimit => elemR(:,er_Temp03)
            DTLimit    => elemR(:,er_Temp04)

            topwidth   => elemR(:,er_TopWidth)
            length     => elemR(:,er_Length)
            depth      => elemR(:,er_Depth)

            alpha      => setting%VariableDT%CFL_inflow_max
            gravity    => setting%Constant%gravity
            smallDepth => setting%SmallDepth%DepthCutoff
        !%------------------------------------------------------------------
        !% ---- use the pack for CC and JM with H time march 
        !%      (no lateral flows into JB or diagnostic)
        thisCol = ep_CCJM_H_ETM
        Npack => npack_elemP(thisCol)   
        if (Npack < 1) return
        thisP => elemP(1:Npack,thisCol)
        
        !% ensure flowrate used for limiter is  positive and non-zero
        !tinyQ(:) = setting%Eps%Velocity
        Qlat(thisP) = max( abs( elemR(thisP,er_FlowrateLateral)) , &
                           setting%Eps%Velocity)

        !print *, 'Qlat ',Qlat(thisP)

        !print *
        !print *, 'topwidth ',topwidth(thisP)

        !% store the scaling depth limit defined by when the induced velocity
        !% from an inflow is similar to the gravity wave speed of the BC inflow
        depthScale(thisP) = ( (alpha**4) * (Qlat(thisP)**2) / (gravity * (topwidth(thisP)**2) ) )**onethirdR

        !print *
        !print *, 'depthScale ',depthScale(thisP)

        !% get the depth limit for maximum depth that the time step will be limited as 
        !% either the depth scale or the specified smallDepth limit
        depthLimit(thisP) = max(smallDepth, depthScale(thisP))

        !print *
        !print *, 'depthLimit ',depthLimit(thisP)
        
        !% time step limit based on inflow flux
        DTlimit(thisP) = length(thisP) * ( ( alpha * topwidth(thisP) / (gravity * Qlat(thisP)) )**onethirdR) 

        !print *
        !print *, 'DTlimit ',DTlimit(thisP)
    
        !% where the depth is greater than the depthlimit the DT inflow limiter
        !% is not needed, and we can use the existing DT value
        depthLimit(thisP) = depth(thisP) - depthLimit(thisP)
        where (depthLimit .ge. zeroR)
            DTlimit = thisDT
        endwhere

        !print *
        !print *, 'new DTlimit ',DTlimit(thisP)

        !% get the smallest DT in the limiter array
        newDTlimit = minval(DTlimit(thisP))

        !% use the smaller value of the new limit or the input DT
        thisDT = min(newDTlimit,thisDT)

        !print *
        !print *, 'thisDT ',thisDT

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

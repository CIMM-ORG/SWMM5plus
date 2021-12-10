module timeloop

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use output, only: outputML_store_data
    use pack_mask_arrays
    use runge_kutta2
    use utility_output
    use boundary_conditions
    use utility_profiler
    use interface, &
        only: interface_export_link_results, &
              interface_get_subcatch_runoff, &
              interface_call_runoff_execute, &
              interface_get_newRunoffTime

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
            integer          :: ii, additional_rows
            logical          :: isTLfinished
            logical          :: doHydraulics, doHydrology
            real(8), pointer :: nextHydrologyTime, nextHydraulicsTime
            real(8), pointer :: lastHydrologyTime, lastHydraulicsTime
            real(8), pointer :: timeEnd, timeNow, dtTol, dtHydraulics
            character(64)    :: subroutine_name = 'timeloop_toplevel'
            integer :: kk !temporary
        !%--------------------------------------------------------------------
        !% Preliminaries
            if (icrash) return
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

            if (setting%Output%Verbose) &
                write(*,"(2A,i5,A)") new_line(" "), 'begin timeloop [Processor ', this_image(), "]"
            call system_clock(count=setting%Time%Real%EpochTimeLoopStartSeconds)
        !%--------------------------------------------------------------------
        !% Aliases
            nextHydrologyTime  => setting%Time%Hydrology%NextTime
            nextHydraulicsTime => setting%Time%Hydraulics%NextTime
            lastHydrologyTime  => setting%Time%Hydrology%LastTime
            lastHydraulicsTime => setting%Time%Hydraulics%LastTime        
            timeEnd => setting%Time%End
            timeNow => setting%Time%Now
            dtTol   => setting%Time%DtTol
            dtHydraulics => setting%Time%Hydraulics%Dt
        !%--------------------------------------------------------------------
        !% initialize the times
        setting%Time%Hydrology%LastTime  = timeNow  
        setting%Time%Hydraulics%LastTime = timeNow

        !% local T/F that changes with each time step depending on whether or 
        !% not the hydrology or hydraulics are conducted in that step
        doHydrology  = setting%simulation%useHydrology 
        doHydraulics = setting%simulation%useHydraulics

        !% get the next hydrology time
        if (doHydrology) then      
            nextHydrologyTime = interface_get_NewRunoffTime()
        else   
            !% set a large dummy time for hydrology if not used
            nextHydrologyTime = timeEnd + onethousandR * dtTol
        end if

        !% get the initial dt and the next hydraulics time
        if (doHydraulics) then
            call tl_update_hydraulics_timestep()
        else
            !% set a large dummy time for hydraulics if not used
            nextHydraulicsTime = timeEnd + onethousandR * dtTol
            !% suppress the multi-level hydraulics output (otherwise seg faults)
            setting%Output%suppress_MultiLevel_Output = .true.
        end if

        !% check to see if there is an initial step to hydrology
        !% if not, then skip the initial tl_hydrology
        if (( abs(nextHydrologyTime - timeNow) < dtTol) .and. (doHydrology) ) then
            doHydrology = .true.
        else
            doHydrology = .false.
        endif 

        !% Combined hydrology (SWMM-C) and hydraulics simulation
        !% The loop starts at t = setting%Time%Start
        !% Adust the end time for the dtTol (precision error in epoch to seconds conversion)
        ii=0
        do ii=1,5
        !do while (setting%Time%Now <= (setting%Time%End - dtTol))
        !    ii = ii+1
            print *, ii, timeNow, setting%Time%End, doHydrology, doHydraulics
            !print *, 'here 3198755'
            !do kk=1,size(elemR,2)
            !    print *, kk, elemR(kk,er_Depth)
            !end do

            if (doHydrology) call tl_hydrology()
            if (doHydraulics) then
                call bc_update()
                call tl_hydraulics()
            end if
            call util_output_report() !% Results must be reported before counter increment
            !% Multilevel time step output
            if ( (setting%Output%report) .and. &
                 (util_output_must_report()) .and. &
                 (.not. setting%Output%suppress_MultiLevel_Output) ) then
                call outputML_store_data (.false.)
            end if
            !% increment the time step and counters for the next time loop
            call tl_increment_counters(doHydraulics, doHydrology)
            if (icrash) then
                if (.not. setting%Output%suppress_MultiLevel_Output) then
                    call outputML_store_data (.true.)
                end if
                exit
            end if

            !print *, 'here 47685'
            !do kk=1,size(elemR,2)
            !    print *, kk, elemR(kk,er_Velocity)
            !end do

        end do

        !% >>> BEGIN HACK
        !%     Temporary for debugging (can be deleted for deployment)
        ! if (setting%Debug%Output) then
        !     !% Write .out in readable .csv
        !     if (this_image() == 1) then
        !         additional_rows = num_images() - 1
        !         do ii = 1, size(Link%I,1) - additional_rows
        !             call interface_export_link_results(ii)
        !         end do
        !     end if
        ! end if
        !% >>> END HACK

        print *, 'finished timeloop'

        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine timeloop_toplevel
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine tl_hydrology()
        ! %-----------------------------------------------------------------------------
        ! % Description:
        ! % Performs a single hydrology step
        ! %-----------------------------------------------------------------------------
        integer :: ii
        character(64) :: subroutine_name = 'tl_hydrology'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !rm !% return if not at the next hydrology time
        !rm if (setting%Time%Now > (setting.Time.NextHydrologyTime - setting.Eps.seconds)) return

        !% execute the EPA SWMM runoff for the next interval
        print *, '...inside tl_hydrology for total subcatchments =', SWMM_N_subcatch    
        call interface_call_runoff_execute()

        !% runoff_execute updates the next runoff time, get the value and store
        setting%Time%Hydrology%NextTime = interface_get_NewRunoffTime()   

        !% cycle through the subcatchments to get the runoff 
        !% starts at 0 because this is index into EPA SWMM storage
        do ii = 0,SWMM_N_subcatch-1  
            call interface_get_subcatch_runoff(ii)  
            print *, '...need to get runoff out of subcatchments'
        end do  

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
            character(64)    :: subroutine_name = 'tl_hydraulics'
        !%-------------------------------------------------------------------
            if (icrash) return
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------

        !% check for where solver needs to switch in dual-solver model
        if (setting%Solver%SolverSelect == ETM_AC) then
            call tl_solver_select()
        end if

        !% repack all the dynamic arrays
        !% FUTURE 20210609 brh need to decide where this goes
        call pack_dynamic_arrays()
        ! print *, "Need to decide on pack_dynamic_arrays 94837"

        !%  push the old values down the stack for AC solver
        call tl_save_previous_values()

        !%  Reset the flowrate adhoc detection before flowrates are updated.
        !%  Note that we do not reset the small volume detection here -- that should
        !%  be in geometry routines.
        elemYN(:,eYN_IsAdhocFlowrate) = .false.

        select case (setting%Solver%SolverSelect)
            case (ETM_AC)
                call rk2_toplevel_ETMAC()
            case (ETM)
                call rk2_toplevel_ETM()
            case (AC)
                call rk2_toplevel_AC()
            case DEFAULT
                print *, 'error, code should not be reached.'
                STOP 1001 !% need error statement
        end select

        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine tl_hydraulics
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_update_hydraulics_timestep()
        !%------------------------------------------------------------------
        !% Description:
        !% updates the timestep (dt) for hydraulics and computes the
        !% setting.Time.Hydraulics.NextTime
        !%------------------------------------------------------------------
            logical, pointer :: matchHydrologyStep, useHydrology
            real(8)          :: oldDT
            real(8)          :: timeleft, thisCFL
            real(8), pointer :: targetCFL, maxCFL, maxCFLlow, timeNow, dtTol
            real(8), pointer :: decreaseFactor, increaseFactor
            real(8), pointer :: newDT
            !rm velocity(:), wavespeed(:), length(:), PCelerity(:)
            real(8), pointer :: nextHydrologyTime, nextHydraulicsTime
            real(8), pointer :: lastHydrologyTime, lastHydraulicsTime
            integer          :: ii, neededSteps
            integer, pointer :: stepNow, stepNext, stepfinal, lastCheckStep, checkStepInterval
            !rm integer, pointer :: thisCol, Npack, thisP(:)
            character(64)    :: subroutine_name = 'tl_update_hydraulics_timestep'

            integer :: kk !temporary
        !%-------------------------------------------------------------------
            if (icrash) return
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases:
            maxCFL             => setting%VariableDT%CFL_hi_max
            targetCFL          => setting%VariableDT%CFL_target
            maxCFLlow          => setting%VariableDT%CFL_lo_max
            decreaseFactor     => setting%VariableDT%decreaseFactor
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
        !% brh2021210 removing JB JM and small volumes from CFL computation
        !thisCol            => col_elemP(ep_CC_Q_NOTsmallvolume)  !col_elemP(ep_CC_ALLtm)
        !Npack              => npack_elemP(thisCol)
        !thisP              => elemP(1:Npack,thisCol)

        !velocity           => elemR(:,er_Velocity)
        !wavespeed          => elemR(:,er_WaveSpeed)
        !length             => elemR(:,er_Length)
        !PCelerity          => elemR(:,er_Preissmann_Celerity)

        oldDT =  setting%Time%Hydraulics%Dt  ! not a pointer
        newDT => setting%Time%Hydraulics%Dt

        !print *, 'here ============================== in tl_update_hydraulics_timestep'

        if ((matchHydrologyStep) .and. (useHydrology)) then
    
            !% brh20211209 
            !rm nextTimeHydrology = (setting%Time%Hydrology%Step + 1) * setting%Time%Hydrology%Dt
            
            !% For combined hydrology and hydraulics compute the CFL if we take a single
            !% step to match the next hydrology time
            timeLeft = nextHydrologyTime - lastHydraulicsTime
            !rm thisCFL = max (maxval((abs(velocity(thisP)) + abs(wavespeed(thisP))) * timeleft / length(thisP)), &
            !rm               maxval((abs(velocity(thisP)) + abs(PCelerity(thisP))) * timeleft / length(thisP)))
            thisCFL = tl_get_max_cfl(ep_CC_Q_NOTsmallvolume,timeleft)

            !print *, 'timeLeft ',timeLeft    
            !print *, 'thisCFL  ',thisCFL    
            !print *, 'maxCFL   ',maxCFL
            !print *, 'oldDT    ',oldDT       

            !% check to see if a single time step to match the hydrology time is possible
            if (thisCFL < maxCFL) then
                !% check to be sure there is significant time left
                if (timeLeft > dtTol) then
                    !% check increase that is implied with a single time step to the next hydrology
                    if (timeleft / oldDT < increaseFactor) then
                        !% use a single hydraulics step for the remaining hydrology time
                        newDT = timeleft
                    else
                        !% increase the oldDT by the maximum allowed
                        newDT = increaseFactor * oldDT
                    end if
                else
                    !% small time left, don't bother with it, go back to the old dt
                    newDT = oldDT
                    !% check that resetting to oldDT didn't cause a problem
                    thisCFL = tl_get_max_cfl(ep_CC_Q_NOTsmallvolume,newDT)
                    if (thisCFL > maxCFL) then 
                        !% if CFL to large, set the time step based on the target CFL
                        newDT = newDT * targetCFL / thisCFL
                    end if
                end if
            else
                !% if more than one time step is needed
                !% compute the needed steps and time step size
                !% first check to see if the implied time step is too small for the integer size
                if (thisCFL/targetCFL .ge. huge(neededSteps)) then
                    write(*,*) 'warning -- really high CFL, setting dt to minimum to prevent overflow'
                    newDT = setting%Limiter%Dt%Minimum + setting%Time%DtTol
                else
                    !% note that neededSteps will be 2 or larger else thisCFL < maxCFL
                    neededSteps = ceiling( thisCFL / targetCFL )
                    !% the provisional time step that would get exactly to the hydrology time (if CFL didn't change)
                    newDT = timeleft / real(neededSteps,8)
                    !% limit the change in the dt by the increase and decrease factor
                    if (newDT / oldDT > increaseFactor) then 
                        newDT = oldDT * increaseFactor
                    elseif (newDT / oldDT < decreaseFactor) then 
                        newDT = oldDT * decreaseFactor
                        !% check that the DT still isn't too large, if so, then take half the timeleft
                        if (newDT > timeleft / twoR) then 
                            newDT = timeleft / twoR 
                            neededSteps = twoI
                        end if
                    else
                        !% accept the provisional newDT
                    end if
                    !% check thatthe newDT didn't cause a CFL violation
                    thisCFL = tl_get_max_cfl(ep_CC_Q_NOTsmallvolume,newDT)
                    if (thisCFL > maxCFL) then 
                        !% if CFL to large, set the time step based on the target CFL
                        newDT = newDT * targetCFL / thisCFL
                    end if
                end if
            end if
        else
            !% Allowing hydrology and hydraulics to occur at different times
            thisCFL = tl_get_max_cfl(ep_CC_Q_NOTsmallvolume,oldDt)
            !rm thisCFL = max (maxval((abs(velocity(thisP)) + abs(wavespeed(thisP))) * oldDT / length(thisP)), &
            !rm               maxval((abs(velocity(thisP)) + abs(PCelerity(thisP))) * oldDT / length(thisP)))

 
            !print *, 'thisCFL  ',thisCFL    
            !print *, 'oldDT    ',oldDT  
            !do kk = 1,size(thisP)
            !    print *, kk, velocity(thisP(kk)), elemR(thisP(kk),er_Depth)
            !end do
            !print *, maxval(abs(velocity(thisP))), maxval(abs(wavespeed(thisP))), maxval(elemR(thisP,er_ell))

            if (thisCFL > maxCFL) then
                !% decrease the time step to the target CFL and reset the checkStep counter
                newDT = oldDT * targetCFL / thisCFL
                lastCheckStep = stepNow
            else
                !% if CFL is less than max, see if it can be raised (only checked at intervals)
                if (stepNow > lastCheckStep + checkStepInterval) then
                    !% check for low CFL only on prescribed intervals and increase time step
                    if (thisCFL < maxCFLlow) then
                        !% increase the time step and reset the checkStep Counter
                        newDT = oldDT * increaseFactor
                        lastCheckStep = stepNow
                    end if
                end if
            end if
        end if

        !% if dt is large and there is more than 2 steps, then round to an integer number
        if (neededSteps > 2) then
            if (newDT > fiveR) then
                !% round larger dt to counting numbers
                newDT = real(floor(newDT),8)
            elseif (newDT > oneR) then
                newDT = real(floor(newDT * onehundredR),8) / onehundredR
            else
                !%  HACK NEED SOMETHING TO ROUND TO 3 DIGITS FOR SMALLER NUMBERS
            end if
        end if

        !print *, 'timeNow,newDT ',timeNow, newDT
        !print *, 'new CFL ', tl_get_max_cfl(ep_CC_Q_NOTsmallvolume,newDt)

        nextHydraulicsTime = lastHydraulicsTime + newDT

        if ((setting%Limiter%Dt%UseLimitMin) .and. (newDT .le. setting%Limiter%Dt%Minimum)) then
            print*, 'timeNow = ', timeNow
            print*, 'dt = ', newDT, 'minDt = ',  setting%Limiter%Dt%Minimum
            print*, 'max velocity  ', maxval( &
                elemR(elemP(1:npack_elemP(ep_CC_Q_NOTsmallvolume),ep_CC_Q_NOTsmallvolume),er_Velocity) )
            print*, 'max wavespeed ', maxval( &
                elemR(elemP(1:npack_elemP(ep_CC_Q_NOTsmallvolume),ep_CC_Q_NOTsmallvolume),er_WaveSpeed) )
            print*, 'warning: the dt value is smaller than the user supplied min dt value'
            stop 1123938
        end if

        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine tl_update_hydraulics_timestep
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine tl_increment_counters(doHydraulics, doHydrology)
        !%-------------------------------------------------------------------
        !% Description:
        !% increments the hydrology and hydraulics step counters and
        !%-------------------------------------------------------------------
        !% Declarations
            logical, intent(inout) :: doHydraulics, doHydrology
            logical, pointer       :: useHydrology, useHydraulics
            real(8), pointer       :: nextHydraulicsTime, nextHydrologyTime
            real(8), pointer       :: lastHydraulicsTime, lastHydrologyTime
            real(8)                :: nextTime, nextTimeLocal
            real(8), pointer       :: dtTol, timeNow, dtHydraulics
            integer                :: minImg
            integer, pointer       :: hydraulicStep, hydrologyStep, step, reportStep
            character(64)          :: subroutine_name = 'tl_increment_counters'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (icrash) return
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%--------------------------------------------------------------------
        !% Aliases
            hydraulicStep => setting%Time%Hydraulics%Step
            hydrologyStep => setting%Time%Hydrology%Step
            dtHydraulics  => setting%Time%Hydraulics%Dt
            timeNow       => setting%Time%Now
            step          => setting%Time%Step
            reportStep    => setting%Output%reportStep
            dtTol         => setting%Time%DtTol

            nextHydraulicsTime => setting%Time%Hydraulics%NextTime
            nextHydrologyTime  => setting%Time%Hydrology%NextTime

            lastHydraulicsTime => setting%Time%Hydraulics%LastTime
            lastHydrologyTime  => setting%Time%Hydrology%LastTime

            useHydrology  => setting%Simulation%useHydrology
            useHydraulics => setting%Simulation%useHydraulics
        !%--------------------------------------------------------------------

        !% get the timestep and the next time for hydraulics
        if (doHydraulics) then
            call tl_update_hydraulics_timestep()
        else
            nextHydraulicsTime = setting%Time%End + tenR*DtTol
        end if    

        !% brh20211209s the hydrology is updated as a setting in tl_hydrology
        !rm nextTimeHydrology = (hydrologyStep + 1) * setting%Time%Hydrology%Dt

        !% The NextHydrologyTime is updated in SWMM, here we just need to
        !% provide a large number if hydrology isn't used
        if (.not. useHydrology) then
            nextHydrologyTime = setting%Time%End + tenR*DtTol
        else
            !% EPA SWMM will not return a next hydrology time that is
            !% beyond the end of the simulation, so we need a special
            !% treatment for this case
            if ((nextHydrologyTime .eq. lastHydrologyTime) .and. &
                (timeNow > setting%Time%Start)) then 
                nextHydrologyTime = setting%Time%End + tenR*DtTol
            end if
        end if

        !print *, 'last times ',lastHydraulicsTime, lastHydrologyTime
        !print *, 'next times ',nextHydraulicsTime, nextHydrologyTime

        !% find the minimum hydraulics time and store accross all processors
        call co_min(nextHydraulicsTime)
        !% take the nextTime as the minimum of either the Hydrology or Hydraulics time
        nextTime = min(nextHydraulicsTime, nextHydrologyTime)

        doHydrology  = (abs(nextTime - nextHydrologyTime)  <= dtTol) .and. useHydrology
        doHydraulics = (abs(nextTime - nextHydraulicsTime) <= dtTol) .and. useHydraulics

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
        ! call co_broadcast(doHydraulics, minImg)
        ! call co_broadcast(doHydrology, minImg)

        if (util_output_must_report()) reportStep = reportStep + 1
        if (doHydraulics) hydraulicStep = hydraulicStep + 1
        if (doHydrology)  hydrologyStep = hydrologyStep + 1

        step    = step + 1
        timeNow = nextTime !timeNow + dt

        !print *, 'timeNow  in tl_increment_counters',timeNow

        call tl_command_line_step_output()

        if (doHydraulics) LastHydraulicsTime = NextHydraulicsTime
        if (doHydrology)  LastHydrologyTime  = NextHydrologyTime

        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        end subroutine tl_increment_counters
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
        if (icrash) return
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
        if (icrash) return
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
            integer, pointer :: step, interval
            integer :: execution_realtime
            real(8) :: simulation_fraction, seconds_to_completion, time_to_completion
            character(3) :: timeunit
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%timeloop) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        dt            => setting%Time%Hydraulics%Dt
        timeNow       => setting%Time%Now
        timeEnd       => setting%Time%End
        step          => setting%Time%Step
        interval      => setting%Output%CommandLine%interval

        call system_clock(count=setting%Time%Real%EpochNowSeconds) ! Fortran function returns real epoch time

        ! estimate the remaining time
        execution_realtime = setting%Time%Real%EpochNowSeconds - setting%Time%Real%EpochTimeLoopStartSeconds
        seconds_to_completion = execution_realtime * (setting%Time%End - setting%Time%Now) &
                                                   / (setting%Time%Now - setting%Time%Start)

        if (setting%Output%Verbose) then
            if (this_image() == 1) then
                if (mod(step,interval) == 0) then
                    ! translate time in seconds into something more useful
                    if (timeNow  < sixtyR) then
                        thistime = timeNow
                        timeunit = 's  '
                    elseif (timeNow >= sixtyR .and. timeNow < seconds_per_hour) then
                        thistime = timeNow / sixtyR
                        timeunit = 'min'
                    elseif (timeNow >= seconds_per_hour .and. timeNow < 3.0*seconds_per_day) then
                        thistime = timeNow / seconds_per_hour
                        timeunit = 'hr '
                    elseif (timeNow >= 3.0 * seconds_per_day) then
                        thistime = timeNow / seconds_per_day
                        timeunit = 'day'
                    endif

                    ! write a time counter
                    write(*,"(A12,i8,a17,F9.2,a1,a3,a6,f5.2)") &
                        'time step = ',step,'; model time = ',thistime, &
                        ' ',timeunit,'; dt = ',dt

                    ! write estimate of time remaining
                    if (seconds_to_completion < sixtyR) then
                        timeunit = 's  '
                        time_to_completion = seconds_to_completion
                    elseif (seconds_to_completion >=sixtyR .and. seconds_to_completion < seconds_per_hour ) then
                        timeunit = 'min'
                        time_to_completion = seconds_to_completion / sixtyR
                    elseif (seconds_to_completion >=seconds_per_hour .and. seconds_to_completion < seconds_per_day) then
                        timeunit = 'hr '
                        time_to_completion = seconds_to_completion / seconds_per_hour
                    else
                        timeunit = 'day'
                        time_to_completion = seconds_to_completion / seconds_per_day
                    endif
                    write(*,"(A9,F6.2,A1,A3,A28)") 'estimate ',time_to_completion,' ',timeunit,' clock time until completion'
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
            real(8), intent(in) :: dt
            integer, pointer :: Npack, thisP(:)
            real(8), pointer :: velocity(:), wavespeed(:), length(:), PCelerity(:)
        !%-------------------------------------------------------------------
            Npack              => npack_elemP(thisCol)
            thisP              => elemP(1:Npack,thisCol)
            velocity           => elemR(:,er_Velocity)
            wavespeed          => elemR(:,er_WaveSpeed)
            length             => elemR(:,er_Length)
            PCelerity          => elemR(:,er_Preissmann_Celerity)
        !%-------------------------------------------------------------------

        outvalue = max (maxval((abs(velocity(thisP)) + abs(wavespeed(thisP))) * dt/ length(thisP)), &
                        maxval((abs(velocity(thisP)) + abs(PCelerity(thisP))) * dt / length(thisP)))

    end function tl_get_max_cfl    
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module timeloop

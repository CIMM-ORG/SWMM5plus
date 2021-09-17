module timeloop

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use pack_mask_arrays
    use runge_kutta2
    use utility_output
    use boundary_conditions
    use utility_profiler
    use utility_prof_jobcount

    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% top-level time-looping of simulation
    !%
    !% METHOD:
    !% Calls hydrology outer loop and then hydraulics inner loop (substeps)
    !%

    private

    public :: timeloop_toplevel

    integer :: itemp

    contains
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    subroutine timeloop_toplevel()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Loops over all the major time-stepping routines
        !%-----------------------------------------------------------------------------
        logical :: isTLfinished
        logical, pointer :: useHydrology, useHydraulics
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'timeloop_toplevel'
        type(wall_clk) :: timer

        real(8) :: start, intermediate, finish
        call cpu_time(start)

        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name
        if (setting%Profile%File%timeloop) call util_tic(timer, 3)

        !%-----------------------------------------------------------------------------
        useHydrology  => setting%Simulation%useHydrology
        useHydraulics => setting%simulation%useHydraulics
        !%-----------------------------------------------------------------------------
        !% logical to detect end of time loop computations
        isTLfinished = .false.
        !%
        
        !% Combined hydrology (SWMM-C) and hydraulics simulation
        !%
        if (useHydrology .and. useHydraulics) then
            !% set the counters used for outer loop iteration
            call tl_setup_counters(hydrology)
            call bc_update()
            !% outer loop (Hydrology) time stepping
            do while (.not. isTLfinished)

                ! print *, '--- in ',trim(subroutine_name),'    AAA'
                ! print *, setting%Time%Hydrology%timeNow, setting%Time%Hydrology%timeNext, &
                !     setting%Time%Hydrology%timeFinal, 'Hydrology timeNow timeNext, timeFinal'
                ! print *, setting%Time%Hydrology%stepNow, setting%Time%Hydrology%stepNext, &
                !     setting%Time%Hydrology%stepFinal, 'Hydrology stepNow, stepNext StepFinal'
                ! print *, setting%Time%Hydrology%Dt, 'Hydrology Dt'
                ! print *, setting%Time%EndTime , 'end time'

                !% Perform one time step of hydrology
                call tl_hydrology()
                !% Call inner loop (multiple subtime steps) of hydraulics
                call tl_hydraulics()
                call tl_increment_counters(hydrology)
                
                call bc_update()
                call tl_check_finish_status(isTLfinished)
                
            !% HACK to prevent infinite loop in testing
            ! print *, "HACK hard-code stopping time loop  39872"
            ! isTLfinished = .true.

            end do !% (while not isTLfinished)
        !%
        !% Hydrology only simulation
        !%
        elseif (useHydrology .and. .not. useHydraulics) then
            call tl_setup_counters(hydrology)
            do while (.not. isTLfinished)
                !% Perform 1 time step of hydrology
                call tl_hydrology()
                call tl_increment_counters(hydrology)
                call tl_check_finish_status(isTLfinished)

            print *, 'error, the useHydrology only does not work.'
            stop 37864

            !% HACK to prevent infinite loop in testing
            print *, "HACK hard-code stopping time loop  93785"
            isTLfinished = .true.

            end do !% (while not isTLfinished)
        !%
        !% Hydraulics only simulation
        !%
        elseif (useHydraulics .and. .not. useHydrology) then
            !% time-loop for hydraulics only is self-contained and doesn't
            !% require an external loop
            call tl_hydraulics()

            print *, 'error, the useHydraulics only does not work.'
            stop 83789

        else
            print *, 'error, condition that should not occur.'
            stop 76408
        endif

        if (setting%Profile%File%timeloop) then
            call util_toc(timer, 3)
            print *, '** time', this_image(),subroutine_name, ' = ', duration(timer%jobs(3))
            ! call util_free_jobs(timer)
        end if

        if (setting%Debug%File%timeloop)  print *, '*** leave ', this_image(), subroutine_name
    end subroutine timeloop_toplevel
    !%
    !%==========================================================================
    !% PRIVATE
    !%==========================================================================
    !%
    subroutine tl_setup_counters(timeloop_type)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Sets up the counters that are used to handle time loops
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: timeloop_type
        integer, pointer :: stepNow, stepNext
        real(8), pointer :: timeNow, timeNext, timeStart, timeEnd, loopTimeFinal, dt
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'tl_setup_counters'
        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------
        select case (timeloop_type)
        case (hydrology)
            stepNow       => setting%Time%Hydrology%stepNow
            stepNext      => setting%Time%Hydrology%stepNext
            timeNow       => setting%Time%Hydrology%timeNow
            timeNext      => setting%Time%Hydrology%timeNext
            dt            => setting%Time%Hydrology%Dt
            loopTimeFinal => setting%Time%Hydrology%timeFinal
            !% Note the timeStart for hydrology (outer loop) is simulation start time
            timeStart => setting%Time%StartTime
            !% Note the timeEnd for hydrology (outer loop) is simulation end time
            timeEnd => setting%Time%EndTime
        case (hydraulics)
            stepNow       => setting%Time%Hydraulics%stepNow
            stepNext      => setting%Time%Hydraulics%stepNext
            timeNow       => setting%Time%Hydraulics%timeNow
            timeNext      => setting%Time%Hydraulics%timeNext
            dt            => setting%Time%Hydraulics%Dt
            loopTimeFinal => setting%Time%Hydraulics%timeFinal
            if (setting%Simulation%useHydrology) then
                !% Combined hydraulics and hydrology
                !% The time start and end for hydraulics are the present hydrology step
                timeStart   => setting%Time%Hydrology%timeNow
                timeEnd     => setting%Time%Hydrology%timeNext
            else
                !% For hydraulics only simulation
                !% The time start and end are the entire simulation
                timeStart   => setting%Time%StartTime
                timeEnd     => setting%Time%EndTime
            endif

        case default
            print *, 'error -- should be unreachable'
            stop 1001
        end select
       !%-----------------------------------------------------------------------------
        stepNow  = 0
        stepNext = 1

        timeNow = timeStart
        timeNext = timeNow + dt
        loopTimeFinal = timeEnd

        if (setting%Debug%File%timeloop) print *, '*** leave ', this_image(), subroutine_name
    end subroutine tl_setup_counters
    !
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine tl_hydrology()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'tl_hydrology'
        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------

        !% need to execute a hydrology step and extract boundary conditions for
        !% upstream, downstream, and lateral inflows.
        ! print *, "Hydrology calls to SWMM-C are needed 8473"
        !% stop 8473

        if (setting%Debug%File%timeloop) print *, '*** leave ', this_image(), subroutine_name
    end subroutine tl_hydrology

    !%==========================================================================
    !%==========================================================================
    !
    subroutine tl_hydraulics()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs inner loop of hydraulics for a single hydrology step
        !%-----------------------------------------------------------------------------
        integer :: ii
        integer, pointer :: stepnext, stepfinal
        real(8), pointer :: timeNext, timeFinal
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'tl_hydraulics'
        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------
        timeNext  => setting%Time%Hydraulics%timeNext
        timeFinal => setting%Time%Hydraulics%timeFinal
        !%-----------------------------------------------------------------------------
        !% set the counters used for inner loop iteration
        call tl_setup_counters(hydraulics)

    
        !% set the expected number of substeps for hydraulics given the present CFL
        call tl_set_hydraulic_substep()

        !% setup for checking volume conservation during the hydraulic steps.
        ! FUTURE 20210609 brh need to decide where to place this and pull code from old version
        ! call diagnostic_volume_conservation

        !% these are hydraulics substeps within a single hydrology step
        do while (timeNext <= timeFinal)
            call tl_update_hydraulic_BC()
            call tl_hydraulic_solver()
            call tl_increment_counters(hydraulics)
            call tl_set_hydraulic_substep()

        !% HACK to prevent infinite loop in testing
        ! print *, "Hard-code hydrualic subtime-step loop exit for testing 7647"
        ! timeNext = timeFinal+1.0

            !print *, '--- in ',trim(subroutine_name),' in hydraulics loop'
            !print *, setting%Time%Hydraulics%timeNow, setting%Time%Hydraulics%timeFinal, 'Hydraulics timeNow, timeFinal'
            !print *, setting%Time%Hydraulics%stepNow, setting%Time%Hydraulics%stepFinal, 'Hydraulics stepNow, StepFinal'    
            !print *, '------------------------------------------'
            
        end do

        if (setting%Debug%File%timeloop) print *, '*** leave ', this_image(), subroutine_name
    end subroutine tl_hydraulics

    !%==========================================================================
    !%==========================================================================
    !%
    subroutine tl_set_hydraulic_substep()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% For combined hydrology and hydraulics simulations, this sets the size of the
        !% hydraulic substep (dt) and the total number of hydraulic substeps in the
        !% inner loop for the present CFL conditions.
        !% For a hydraulics-only simulation, this adjust the time steps up or down from
        !% its present value depending on the brackets of setting.VariableDT.CFL_hi_max and
        !% setting.VariableDT.CFL_lo_max. Note that the increase of dt for low CFL will
        !% only occur every N time steps, as set by setting.VariableDT.NstepsForCheck
        !%-----------------------------------------------------------------------------
        real(8) :: timeleft, thisCFL
        integer :: neededSteps, ii

        real(8), pointer :: dt, maxCFL, maxCFLlow, targetCFL
        real(8), pointer :: timeNow, timeNext, timeFinal, decreaseFactor, increaseFactor
        real(8), pointer :: velocity(:), wavespeed(:), length(:)
        integer, pointer :: stepNow, stepNext, stepfinal, checkStepInterval, lastCheckStep
        integer, pointer :: thisCol, Npack, thisP(:)
        logical, pointer :: useHydrology, useHydraulics
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'tl_set_hydraulic_substep'
        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------
        useHydrology  => setting%Simulation%useHydrology
        useHydraulics => setting%Simulation%useHydraulics
        dt        => setting%Time%Hydraulics%Dt
        timeNow   => setting%Time%Hydraulics%timeNow
        timeNext  => setting%Time%Hydraulics%timeNext
        timeFinal => setting%Time%Hydraulics%timeFinal
        stepNow   => setting%Time%Hydraulics%stepNow
        stepNext  => setting%Time%Hydraulics%stepNext
        stepFinal => setting%Time%Hydraulics%stepFinal
        !% Note that the timeFinal above assumes that this correctly stores either the
        !% total simulation time (for a hydraulics-only simulation) or the next hydrology
        !% step time (for a combined hydrology-hydraulics simulation)
        maxCFL            => setting%VariableDT%CFL_hi_max
        targetCFL         => setting%VariableDT%CFL_target
        maxCFLlow         => setting%VariableDT%CFL_lo_max
        decreaseFactor    => setting%VariableDT%decreaseFactor
        increaseFactor    => setting%VariableDT%increaseFactor
        checkStepInterval => setting%VariableDT%NstepsForCheck
        lastCheckStep     => setting%VariableDT%LastCheckStep

        thisCol => col_elemP(ep_CC_ALLtm)
        Npack   => npack_elemP(thisCol)
        thisP   => elemP(1:Npack,thisCol)

        velocity  => elemR(:,er_Velocity)
        wavespeed => elemR(:,er_WaveSpeed)
        length    => elemR(:,er_Length)
        !%-----------------------------------------------------------------------------
        !% how much time is remaining in the inner loop (or entire simulation)
        timeleft = timeFinal - timeNow

        !print *, '---- in ',trim(subroutine_name),'    z01'
        !print *, timeleft, ' timeleft'

        if (timeleft > zeroR) then
            !% compute the maximum CFL if a single step is taken

            if (useHydrology) then
                !% For combined hydrology and hydraulics use the hydrology time step
                !% as the target CFL.
                thisCFL = maxval( (velocity(thisP) + wavespeed(thisP)) * timeleft / length(thisP) )

                !% check to see if max CFL is exceeded
                if (thisCFL < maxCFL) then
                    !% use a single hydraulics step for the remaining hydrology time
                    dt = timeleft
                    stepFinal = stepNext
                else
                    !% compute the needed steps and time step size
                    neededSteps = ceiling(thisCFL / targetCFL)
                    stepFinal = stepNow + neededSteps
                    dt = timeleft / real(neededSteps,8)
                    if ((dt > fiveR) .and. (neededSteps > 2)) then
                        !% round larger dt to integer values
                        dt = real(floor(dt),8)
                    endif


                endif
            else
                !% For hydraulics only, keep the timestep stable unless it
                !% exceeds CFL limits (both high and low limits).
                thisCFL = maxval( (velocity(thisP) + wavespeed(thisP)) * dt / length(thisP) )

                if (thisCFL > maxCFL) then
                    !% decrease the time step and reset the checkStep counter
                    dt = dt * decreaseFactor * maxCFL /thisCFL
                    lastCheckStep = stepNow
                else
                    if (stepNow > lastCheckStep + checkStepInterval) then
                        !% check for low CFL only on prescribed intervals and increase time step
                        if (thisCFL < maxCFLlow) then
                            !% increase the time step and reset the checkStep Counter
                            dt = dt * increaseFactor
                            lastCheckStep = stepNow
                        endif
                    endif
                endif  
            endif
        else
            !% for timeleft <= 0 there is no change as the hydraulics loop should exit
        endif

        !% set the smallest time step on any processor as the time step
        sync all
        call co_min(dt)
        timeNext = timeNow + dt

        if ((setting%Limiter%Dt%UseLimitMin) .and. (dt <= setting%Limiter%Dt%Minimum)) then
            print*, 'timeNow = ', timeNow
            print*, 'dt = ', dt, 'minDt = ',  setting%Limiter%Dt%Minimum
            print*, 'warning: the dt value is smaller than the user supplied min dt value'
            stop 3369
        endif

        if (setting%Debug%File%timeloop) print *, '*** leave ', this_image(), subroutine_name
    end subroutine tl_set_hydraulic_substep
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine tl_update_hydraulic_BC()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Updates (if needed) BC to hydraulic solver, including upstream, downstream
        !% and lateral inflows.
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'tl_update_hydraulic_BC'
        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------
        !% This needs to take the BC from the hydrology step (obtained in tl_hydrology)
        !% and subdivide for the subtime stepping of the hydraulics. For flow rates this
        !% is a simple task -- take the flow over the hydrology timestep and subdivide
        !% by the hydraulics time step (however ,this gets a little tricky if the hydraulic
        !% timestep is allowed to change during the overarching hydrology timestep, which
        !% may be necessary for stability). For the elevation BC we will need to think
        !% more carefully. Let's imagine that we have a hydrology time step of 15 minutes.
        !% Does the hydrology elevation BC represent the average water surface elevation
        !% during those 15 minutes? or does it represent the instantaneous water surface
        !% elevation at the start of the 15 minute step? Let us assume that it represents
        !% average elevation over the hydrology time step "m".  In which case, we can make
        !% an estimate of the  surface elevation at the start of the time step as
        !% H^{m-1/2} = (H^{m} + H^{m-1}) / 2
        !% Furthermore, the rate of change of the water surface elevation can be estimated
        !% for the 15 minute time step as a simple difference
        !% dH/dt = ( H^{m} - H^{m-1} ) / (15 * 60)
        !% Let us assume that we use a 3 minute subtime step for hydraulics, in which case
        !% for the n subtime step we have
        !% H^{n} = H^{m-1/2} + n(3 * 60) dH/dt
        !% The above should be written out more clearly in the SWMM5+ Code Narration.

        ! print *, "Need tl_updated_hydraulic_BC to be written 38972"

        if (setting%Debug%File%timeloop) print *, '*** leave ', this_image(), subroutine_name
    end subroutine tl_update_hydraulic_BC
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine tl_hydraulic_solver()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Top level hydraulic solver for a single time step
        !%-----------------------------------------------------------------------------
        real(8) :: thisCFL

        real(8), pointer :: dt, timeNow
        real(8), pointer :: velocity(:), wavespeed(:), length(:)
        integer, pointer :: thisCol, Npack, thisP(:)
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'tl_hydraulic_solver'
        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------
        dt        => setting%Time%Hydraulics%Dt
        timeNow   => setting%Time%Hydraulics%timeNow
        thisCol   => col_elemP(ep_CC_ALLtm)
        Npack     => npack_elemP(thisCol)
        thisP     => elemP(1:Npack,thisCol)

        velocity  => elemR(:,er_Velocity)
        wavespeed => elemR(:,er_WaveSpeed)
        length    => elemR(:,er_Length)

        ! BRH useful for debugging -- comment out if you want quiet
        !print *, '** in ',trim(subroutine_name), ' timeNow=',timeNow,' dt=', dt

        !% check for where solver needs to switch in dual-solver model
        if (setting%Solver%SolverSelect == ETM_AC) then
            call tl_solver_select()
        endif

        !% repack all the dynamic arrays
        !% FUTURE 20210609 brh need to decide where this goes
        call pack_dynamic_arrays()
        ! print *, "Need to decide on pack_dynamic_arrays 94837"

        !%  push the old values down the stack for AC solver
        call tl_save_previous_values()

        !% print the cfl to check for model blowout
        call util_output_report_summary()

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

        !% report timestep
        if (setting%Output%report) call util_output_report()

        if (setting%Debug%File%timeloop) print *, '*** leave ', this_image(), subroutine_name
    end subroutine tl_hydraulic_solver
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine tl_increment_counters(timeloop_type)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Increments the counters that are used to monitor hydrology (outerloop)
        !%-----------------------------------------------------------------------------

        integer, intent(in) :: timeloop_type

        integer, pointer :: thisstep, nextstep
        real(8), pointer :: thistime, nexttime, dt
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'tl_increment_counters'
        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------

        select case (timeloop_type)
        case (hydrology)
            thistime => setting%Time%Hydrology%timeNow
            nexttime => setting%Time%Hydrology%timeNext
            dt       => setting%Time%Hydrology%Dt
            thisstep => setting%Time%Hydrology%stepNow
            nextstep => setting%Time%Hydrology%stepNext
        case (hydraulics)
            thistime => setting%Time%Hydraulics%timeNow
            nexttime => setting%Time%Hydraulics%timeNext
            dt       => setting%Time%Hydraulics%Dt
            thisstep => setting%Time%Hydraulics%stepNow
            nextstep => setting%Time%Hydraulics%stepNext
        case default
            print *, 'error should be unreachable'
            stop 1001
        end select
        !%-----------------------------------------------------------------------------
        !% increment time
        thistime = nexttime
        nexttime = nexttime + dt

        !% increment step
        thisstep = nextstep
        nextstep = nextstep + 1

        !% HACK to prevent inifinite looping
        !print *, '======================= in ',trim(subroutine_name),'; nexstep = ',nextstep
        if (nextstep > 100000) then
            print *, 'stopping because too many time steps (HACK) in ',subroutine_name
            stop 2987
        endif

        if (setting%Debug%File%timeloop) print *, '*** leave ', subroutine_name
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
        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------

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

        if (setting%Debug%File%timeloop) print *, '*** leave ', this_image(), subroutine_name
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
        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------
        !%  push the old values down the stack
        !%  N values is the present, N0 is the last time step, and N1
        !%  is the timestep before (needed only for backwards 3rd in velocity and volume)
        elemR(:,er_Flowrate_N0)  = elemR(:,er_Flowrate)
        elemR(:,er_Head_N0)      = elemR(:,er_Head)
        elemR(:,er_Velocity_N1)  = elemR(:,er_Velocity_N0)
        elemR(:,er_Velocity_N0)  = elemR(:,er_Velocity)
        elemR(:,er_Volume_N1)    = elemR(:,er_Volume_N0)
        elemR(:,er_Volume_N0)    = elemR(:,er_Volume)

        if (setting%Debug%File%timeloop) print *, '*** leave ', this_image(), subroutine_name
    end subroutine tl_save_previous_values
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine tl_check_finish_status(isTLfinished)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Signasl the timeloop is finished when next time step would be beyond the
        !% end time.
        !%-----------------------------------------------------------------------------
        logical, intent(inout) :: isTLfinished

        real(8), pointer :: endtime, thistime
        integer, pointer :: thisstep, finalstep
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'tl_check_finish_status'
        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------

        endtime   => setting%Time%EndTime
        thistime  => setting%Time%Hydrology%timeNow
        thisstep  => setting%Time%Hydrology%stepNow
        finalstep => setting%Time%Hydrology%stepFinal
        !%-----------------------------------------------------------------------------

        if ((thistime >= endtime) .or. (thisstep >= finalstep)) then
            isTLfinished = .true.
        endif

        !% BRHbugfix 20210813 start
        if (isTLfinished) then
            print *, '** in ',subroutine_name
            print *,  thistime, ' = thistime'
            print *, endtime, ' = planned endtime'
        endif
        !% BRHbugfix 20210813 end
        
        !% FUTURE brh 20210607 Need a control to exit on error
        if (setting%Debug%File%timeloop) print *, '*** leave ', this_image(), subroutine_name
    end subroutine tl_check_finish_status
    !%
    !%==========================================================================
    !% END OF MODULE
    !%+=========================================================================
end module timeloop

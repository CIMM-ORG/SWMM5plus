module timeloop

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use pack_mask_arrays
    use runge_kutta2
    use utility_output
    use boundary_conditions
    use interface, only: interface_export_link_results

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

    !%==========================================================================
    !% PUBLIC
    !%==========================================================================

    subroutine timeloop_toplevel()
    !%-----------------------------------------------------------------------------
    !% Description:
    !%     Loops over all the major time-stepping routines
    !%-----------------------------------------------------------------------------
        integer          :: ii, additional_rows
        logical          :: isTLfinished
        logical          :: doHydraulics, doHydrology
        character(64)    :: subroutine_name = 'timeloop_toplevel'
    !%-----------------------------------------------------------------------------
        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name

        doHydraulics = setting%simulation%useHydraulics
        doHydrology = setting%simulation%useHydrology

        !% Combined hydrology (SWMM-C) and hydraulics simulation
        !% The loop starts at t = setting%Time%Start
        do while (setting%Time%Now <= setting%Time%End)
            if (doHydrology) call tl_hydrology()
            if (doHydraulics) then
                call bc_update()
                call tl_hydraulics()
            end if
            call util_output_report() !% Results must be reported before counter increment
            call tl_increment_counters(doHydraulics, doHydrology)
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

        if (setting%Debug%File%timeloop) print *, '*** leave ', this_image(), subroutine_name
    end subroutine timeloop_toplevel

    !%==========================================================================
    !% PRIVATE
    !%==========================================================================

    subroutine tl_hydrology()
    !%-----------------------------------------------------------------------------
    !% Description:
    !% Performs a single hydrology step
    !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'tl_hydrology'
    !%-----------------------------------------------------------------------------
        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name


        if (setting%Debug%File%timeloop) print *, '*** leave ', this_image(), subroutine_name
    end subroutine tl_hydrology
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine tl_hydraulics()
    !%-----------------------------------------------------------------------------
    !% Description:
    !% Top level hydraulic solver for a single time step
    !%-----------------------------------------------------------------------------
        character(64)    :: subroutine_name = 'tl_hydraulics'
    !%-----------------------------------------------------------------------------
        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name

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

        if (setting%Debug%File%timeloop) print *, '*** leave ', this_image(), subroutine_name
    end subroutine tl_hydraulics
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine tl_update_hydraulic_step()
    !%-----------------------------------------------------------------------------
    !% Description:
    !%-----------------------------------------------------------------------------

        logical          :: matchHydrologyStep
        real(8)          :: timeleft, timeNow, thisCFL, targetCFL, maxCFL, maxCFLlow
        real(8)          :: decreaseFactor, increaseFactor, nextTimeHydrology
        real(8), pointer :: dt, velocity(:), wavespeed(:), length(:)
        integer          :: ii, neededSteps, checkStepInterval
        integer, pointer :: stepNow, stepNext, stepfinal, lastCheckStep
        integer, pointer :: thisCol, Npack, thisP(:)
        character(64)    :: subroutine_name = 'tl_update_hydraulic_step'

    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name

        maxCFL             = setting%VariableDT%CFL_hi_max
        targetCFL          = setting%VariableDT%CFL_target
        maxCFLlow          = setting%VariableDT%CFL_lo_max
        decreaseFactor     = setting%VariableDT%decreaseFactor
        increaseFactor     = setting%VariableDT%increaseFactor
        checkStepInterval  = setting%VariableDT%NstepsForCheck
        matchHydrologyStep = setting%Time%matchHydrologyStep
        timeNow            = setting%Time%Now

        dt                => setting%Time%Hydraulics%Dt
        stepNow           => setting%Time%Step
        lastCheckStep     => setting%VariableDT%LastCheckStep

        thisCol           => col_elemP(ep_CC_ALLtm)
        Npack             => npack_elemP(thisCol)
        thisP             => elemP(1:Npack,thisCol)

        velocity          => elemR(:,er_Velocity)
        wavespeed         => elemR(:,er_WaveSpeed)
        length            => elemR(:,er_Length)

        if (matchHydrologyStep) then
            !% For combined hydrology and hydraulics use the hydrology time step
            !% as the target CFL.
            nextTimeHydrology = (setting%Time%Hydrology%Step + 1) * setting%Time%Hydrology%Dt
            timeLeft = nextTimeHydrology - timeNow
            thisCFL = maxval( (abs(velocity(thisP)) + abs(wavespeed(thisP))) * timeleft / length(thisP) )

            !% check to see if max CFL is exceeded
            if (thisCFL < maxCFL) then
                !% use a single hydraulics step for the remaining hydrology time
                dt = timeleft
            else
                !% compute the needed steps and time step size
                neededSteps = ceiling( thisCFL / targetCFL )
                dt = timeleft / real(neededSteps,8)
                if ((dt > fiveR) .and. (neededSteps > 2)) then
                    !% round larger dt to integer values
                    dt = real(floor(dt),8)
                end if
            end if
        else
            !% For hydraulics only, keep the timestep stable unless it
            !% exceeds CFL limits (both high and low limits).
            thisCFL = maxval( (abs(velocity(thisP)) + abs(wavespeed(thisP))) * dt / length(thisP) )

            if (thisCFL > maxCFL) then
                !% decrease the time step and reset the checkStep counter
                dt = dt * decreaseFactor * maxCFL / thisCFL
                lastCheckStep = stepNow
            else
                if (stepNow > lastCheckStep + checkStepInterval) then
                    !% check for low CFL only on prescribed intervals and increase time step
                    if (thisCFL < maxCFLlow) then
                        !% increase the time step and reset the checkStep Counter
                        dt = dt * increaseFactor
                        lastCheckStep = stepNow
                    end if
                end if
            end if
        end if

        if ((setting%Limiter%Dt%UseLimitMin) .and. (dt <= setting%Limiter%Dt%Minimum)) then
            print*, 'timeNow = ', timeNow
            print*, 'dt = ', dt, 'minDt = ',  setting%Limiter%Dt%Minimum
            print*, 'max velocity', maxval(velocity(thisP)), 'max wavespeed', maxval(wavespeed(thisP))
            print*, 'warning: the dt value is smaller than the user supplied min dt value'
            stop 1123
        end if

        if (setting%Debug%File%timeloop) print *, '*** leave ', this_image(), subroutine_name
    end subroutine tl_update_hydraulic_step
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine tl_increment_counters(doHydraulics, doHydrology)
    !%-----------------------------------------------------------------------------
    !% Description:
    !%-----------------------------------------------------------------------------
        logical, intent(inout) :: doHydraulics, doHydrology
        logical                :: useHydrology, useHydraulics
        real(8)                :: nextTimeHydraulics, nextTimeHydrology, nextTime, dtTol
        real(8), pointer       :: timeNow, dt
        integer                :: minImg
        integer, pointer       :: hydraulicStep, hydrologyStep, step, reportStep
        character(64)          :: subroutine_name = 'tl_increment_counters'
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name

        hydraulicStep => setting%Time%Hydraulics%Step
        hydrologyStep => setting%Time%Hydrology%Step
        dt            => setting%Time%Dt
        timeNow       => setting%Time%Now
        step          => setting%Time%Step
        reportStep    => setting%Output%reportStep

        dtTol         = setting%Time%DtTol
        useHydrology  = setting%Simulation%useHydrology
        useHydraulics = setting%Simulation%useHydraulics

        !% Check if CFL > CFLmax condition
        if (useHydraulics) call tl_update_hydraulic_step()

        nextTimeHydraulics = timeNow + setting%Time%Hydraulics%Dt
        nextTimeHydrology = (hydrologyStep + 1) * setting%Time%Hydrology%Dt

        !% Check what needs to be computed next (Hydrology | Hydraulics)

        !% 1. Update simulation time step
        nextTime = min(nextTimeHydraulics, nextTimeHydrology)
        doHydrology  = (abs(nextTime - nextTimeHydrology) <= dtTol) .and. useHydrology
        doHydraulics = (abs(nextTime - nextTimeHydraulics) <= dtTol) .and. useHydraulics
        dt = nextTime - timeNow

        !% 2. Communicate new Dt (ALL the processors compute the same step)
        call co_min(dt)
        if (dt == (nextTime - timeNow)) then
            minImg = this_image()
        else
            minImg = -1
        end if
        call co_max(minImg)
        call co_broadcast(doHydraulics, minImg)
        call co_broadcast(doHydrology, minImg)

        if (util_output_must_report()) reportStep = reportStep + 1
        if (doHydraulics) hydraulicStep = hydraulicStep + 1
        if (doHydrology) hydrologyStep = hydrologyStep + 1

        step    = step + 1
        timeNow = timeNow + dt

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
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name
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
    !% END OF MODULE
    !%+=========================================================================
end module timeloop

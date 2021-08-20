module timeloop

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use pack_mask_arrays
    use runge_kutta2
    use utility_output
    use boundary_conditions
    use, intrinsic :: ISO_FORTRAN_ENV, only: team_type

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
        logical          :: isTLfinished
        logical          :: doHydraulics = .true.
        logical, pointer :: useHydrology, useHydraulics
        character(64)    :: subroutine_name = 'timeloop_toplevel'
    !%-----------------------------------------------------------------------------
        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name

        useHydrology  => setting%Simulation%useHydrology
        useHydraulics => setting%simulation%useHydraulics

        !% Combined hydrology (SWMM-C) and hydraulics simulation
        do while (setting%Time%Now <= setting%Time%End)

            if (.not. useHydraulics) stop "ERROR - Hydraulics solver disabled"

            call bc_update()
            if (doHydraulics) call tl_hydraulics()
            doHydraulics = tl_increment_counters()
        end do

        if (setting%Debug%File%timeloop)  print *, '*** leave ', this_image(), subroutine_name
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

        call bc_update()

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

        call tl_update_hydraulic_step()

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
    end subroutine tl_hydraulics
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine tl_update_hydraulic_step()
    !%-----------------------------------------------------------------------------
    !% Description:
    !%-----------------------------------------------------------------------------

        real(8) :: timeleft, thisCFL
        integer :: neededSteps, ii
        real(8), pointer :: dt, maxCFL, maxCFLlow, targetCFL
        real(8), pointer :: timeNow, timeNext, timeFinal, decreaseFactor, increaseFactor
        real(8), pointer :: velocity(:), wavespeed(:), length(:)
        integer, pointer :: stepNow, stepNext, stepfinal, checkStepInterval, lastCheckStep
        integer, pointer :: thisCol, Npack, thisP(:)
        character(64)    :: subroutine_name = 'tl_update_hydraulic_step'

    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name

        dt                => setting%Time%HydraulicStep
        maxCFL            => setting%VariableDT%CFL_hi_max
        targetCFL         => setting%VariableDT%CFL_target
        maxCFLlow         => setting%VariableDT%CFL_lo_max
        decreaseFactor    => setting%VariableDT%decreaseFactor
        increaseFactor    => setting%VariableDT%increaseFactor
        checkStepInterval => setting%VariableDT%NstepsForCheck
        stepNow           => setting%Time%Step
        lastCheckStep     => setting%VariableDT%LastCheckStep

        thisCol           => col_elemP(ep_CC_ALLtm)
        Npack             => npack_elemP(thisCol)
        thisP             => elemP(1:Npack,thisCol)

        velocity          => elemR(:,er_Velocity)
        wavespeed         => elemR(:,er_WaveSpeed)
        length            => elemR(:,er_Length)

        !% For hydraulics only, keep the timestep stable unless it
        !% exceeds CFL limits (both high and low limits).
        thisCFL = maxval( (velocity(thisP) + wavespeed(thisP)) * dt / length(thisP) )

        if (thisCFL > maxCFL) then
            !% decrease the time step and reset the checkStep counter
            dt = dt * decreaseFactor * maxCFL / thisCFL
            lastCheckStep = stepNow
        else
            if (dt > lastCheckStep + checkStepInterval) then
                !% check for low CFL only on prescribed intervals and increase time step
                if (thisCFL < maxCFLlow) then
                    !% increase the time step and reset the checkStep Counter
                    dt = dt * increaseFactor
                    lastCheckStep = stepNow
                endif
            endif
        endif

        print *, dt, maxval(velocity(thisP))
        if ((setting%Limiter%Dt%UseLimitMin) .and. (dt <= setting%Limiter%Dt%Minimum)) then
            print*, 'timeNow = ', timeNow
            print*, 'dt = ', dt, 'minDt = ',  setting%Limiter%Dt%Minimum
            print*, 'warning: the dt value is smaller than the user supplied min dt value'
            stop 3369
        endif

        if (setting%Debug%File%timeloop) print *, '*** leave ', this_image(), subroutine_name
    end subroutine tl_update_hydraulic_step
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    function tl_increment_counters() result(doHydraulics)
    !%-----------------------------------------------------------------------------
    !% Description:
    !%-----------------------------------------------------------------------------
        logical         :: doHydraulics
        real(8)         :: newDt
        character(64)   :: subroutine_name     = 'tl_increment_counters'
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%timeloop) print *, '*** enter ', this_image(), subroutine_name

        !% Check what needs to be computed next (Hydrology | Hydraulics)
        !% 1. Update simulation time step
        if (setting%Simulation%useHydrology) then
            setting%Time%Dt = min(setting%Time%HydraulicStep, setting%Time%HydrologyStep)
        else
            setting%Time%Dt = setting%Time%HydraulicStep
        end if

        !% 2. Enable routing if Hydraulics is next
        doHydraulics = (setting%Time%Dt == setting%Time%HydraulicStep)
        !% Only the processors that run hydraulics next, exchange the minimum time step
        newDt = setting%Time%Dt
        call co_min(newDt)
        if (doHydraulics) setting%Time%Dt = newDt

        setting%Time%Step = setting%Time%Step + 1
        setting%Time%Now = setting%Time%Now + setting%Time%Dt

        if (setting%Debug%File%timeloop) print *, '*** leave ', subroutine_name
    end function tl_increment_counters
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

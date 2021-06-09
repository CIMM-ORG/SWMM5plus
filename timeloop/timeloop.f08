module timeloop

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use runge_kutta2

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
        useHydrology  => setting%Simulation%useHydrology
        useHydraulics => setting%simulation%useHydraulics
        !%-----------------------------------------------------------------------------
        isTLfinished = .false.
        !%
        !% Combined hydrology and hydraulics simulation
        !%
        !% HACK 20210608 BRH commented out until further testing
        if (useHydrology .and. useHydraulics) then
            !% set the counters used for outer loop iteration
            call tl_setup_counters(hydrology)     
            !% outer loop (Hydrology) time stepping

            do while (.not. isTLfinished)
                !% Perform 1 time step of hydrology
                call tl_hydrology()
                !% Call inner loop (multiple subtime steps) of hydraulics
                call tl_hydraulics()
                call tl_increment_counters(hydrology)
                call tl_check_finish_status(isTLfinished)
            !% HACK to prevent infinite loop in testing
            isTLfinished = .true.                
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
            !% HACK to prevent infinite loop in testing
            isTLfinished = .true.                    
            end do !% (while not isTLfinished)
        !%
        !% Hydraulics only simulation
        !%    
        elseif (useHydraulics .and. .not. useHydrology) then
            !% time-loop for hydraulics is self-contained
            call tl_hydraulics()
        else
            print *, 'error, condition that should not occur.'  
        endif  
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
                timeStart =>  setting%Time%StartTime    
                timeEnd => setting%Time%EndTime
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

    end subroutine tl_hydrology    
    
    !%==========================================================================  
    !%==========================================================================  
    !
    subroutine tl_hydraulics()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs inner loop of hydraulics for a single hydrology step 
        !%-----------------------------------------------------------------------------
        integer, pointer :: stepnext, stepfinal
        real(8), pointer :: timeNext, timeFinal
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
        timeNext = timeFinal+1.0            
        end do

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
        integer :: neededSteps
        
        real(8), pointer :: dt, maxCFL, maxCFLlow, targetCFL 
        real(8), pointer :: timeNow, timeFinal, decreaseFactor, increaseFactor
        real(8), pointer :: velocity(:), wavespeed(:), length(:)
        integer, pointer :: stepNow, stepNext, stepfinal, checkStepInterval, lastCheckStep
        integer, pointer :: thisCol, Npack, thisP(:)
        logical, pointer :: useHydrology, useHydraulics
        !%-----------------------------------------------------------------------------
        useHydrology  => setting%Simulation%useHydrology
        useHydraulics => setting%Simulation%useHydraulics
        dt        => setting%Time%Hydraulics%Dt 
        timeNow   => setting%Time%Hydraulics%timeNow
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

        thisCol => col_elemP(ep_ALLtm)
        Npack   => col_elemP(thisCol)
        thisP   => elemP(1:Npack,thisCol)

        velocity  => elemR(:,er_Velocity)
        wavespeed => elemR(:,er_WaveSpeed)
        length    => elemR(:,er_Length)
        !%-----------------------------------------------------------------------------
        !% how much time is remaining in the inner loop (or entire simulation)
        timeleft = timeFinal - timeNow

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

        !% check for where solver needs to switch
        if (setting%Solver%SolverSelect == ETM_AC) then
            call tl_solver_select ()
        endif    
        
        !% repack all the dynamic arrays
        !% FUTURE 20210609 brh need to decide where this goes
        !call pack_dynamic_arrays
        
        !%  push the old values down the stack for AC solver
        call tl_save_previous_values ()
        
        !%  Reset the flowrate adhoc detection before flowrates are updated.
        !%  Note that we do not reset the small volume detection here -- that should
        !%  be in geometry routines.
        elemYN(:,eYN_IsAdhocFlowrate) = .false.
        
        select case (setting%Solver%SolverSelect)
            case (ETM_AC)
                call rk2_ETMAC_toplevel ()
                call rk2_AC_toplevel ()
            case (ETM)
                call rk2_ETM_toplevel ()
            case (AC)
                call rk2_AC_toplevel ()
            case DEFAULT
                print *, 'error, code should not be reached.'
                STOP 1001 !% need error statement
        end select
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

    end subroutine tl_solver_select
    !%
    !%==========================================================================  
    !%==========================================================================  
    !%   
    subroutine tl_save_previous_values ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Pushes the time N values into time N-1 storage, and the time N+1 values into
        !% the time N storage. 
        !% HACK -- 20210809 brh This would be better done by changing the indexes without
        !% moving the data, but we will wait on doing this until we have the code
        !% debugged.
        !%-----------------------------------------------------------------------------
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
        endtime   => setting%Time%EndTime
        thistime  => setting%Time%Hydrology%timeNow
        thisstep  => setting%Time%Hydrology%stepNow
        finalstep => setting%Time%Hydrology%stepFinal
        !%-----------------------------------------------------------------------------

        if ((thistime > endtime) .or. (thisstep > finalstep)) then
            isTLfinished = .true.
        endif

        !% FUTURE brh 20210607 Need a control to exit on error
    
    end subroutine tl_check_finish_status 
    !%  
    !%==========================================================================
    !% END OF MODULE
    !%+=========================================================================
end module timeloop
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
    use interface, &
        only: interface_export_link_results, &
              interface_get_subcatch_runoff, &
              interface_call_runoff_execute, &
              interface_get_newRunoffTime
    use utility_crash

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
            integer          :: ii, additional_rows, iblank
            logical          :: isTLfinished
            logical          :: doHydraulics, doHydrology
            logical, pointer :: useHydraulics, useHydrology
            real(8), pointer :: nextHydrologyTime, nextHydraulicsTime
            real(8), pointer :: lastHydrologyTime, lastHydraulicsTime
            real(8), pointer :: timeEnd, timeNow, dtTol, dtHydraulics
            character(64)    :: subroutine_name = 'timeloop_toplevel'
            integer(kind=8) :: cval, crate, cmax
            integer :: kk !temporary

            !% temporary for lateral flow testing
            integer :: mm
            real(8), pointer :: Qlateral(:), Qrate(:)
            integer, pointer :: thisCol, npack, thisP(:), thisBC(:)
            integer, pointer :: sImage(:), eIdx(:)

            integer, allocatable :: tempP(:)

            logical :: writeYN

        !%--------------------------------------------------------------------
        !% Preliminaries
            !if (crashYN) return
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

            !if (setting%Output%Verbose) &
            !    write(*,"(2A,i5,A)") new_line(" "), 'begin timeloop [Processor ', this_image(), "]"
            !if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin timeloop'
        !%--------------------------------------------------------------------
        !% Aliases
            useHydrology       => setting%Simulation%useHydrology
            useHydraulics      => setting%Simulation%useHydraulics
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
        doHydrology  = useHydrology 
        doHydraulics = useHydraulics

        !% get the next hydrology time
        if (useHydrology) then      
            nextHydrologyTime = interface_get_NewRunoffTime()
        else   
            !% set a large dummy time for hydrology if not used
            nextHydrologyTime = timeEnd + onethousandR * dtTol
        end if

        !% get the initial dt and the next hydraulics time
        if (useHydraulics) then
            call tl_update_hydraulics_timestep()
            call util_crashstop(229873)
        else
            !% NOTE -- WORKING WITHOUT SWMM5+ HYDRAULICS IS NOT SUPPORTED 
            !% the following is stub routines.
            !% set a large dummy time for hydraulics if not used
            nextHydraulicsTime = timeEnd + onethousandR * dtTol
            !% suppress the multi-level hydraulics output (otherwise seg faults)
            setting%Output%Report%suppress_MultiLevel_Output = .true.
        end if

        !% check to see if there is an initial step to hydrology
        !% if not, then skip the initial tl_hydrology
        if (( abs(nextHydrologyTime - timeNow) < dtTol) .and. (doHydrology) ) then
            doHydrology = .true.
        else
            doHydrology = .false.
        endif 

        if (this_image()==1) then
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%TimeMarchStart = cval
        end if 
        
        !% ========================================================================================
        !% BEGIN DO LOOP
        !%
        !% Combined hydrology (SWMM-C) and hydraulics simulation
        !% The loop starts at t = setting%Time%Start
        !% Adust the end time for the dtTol (precision error in epoch to seconds conversion)
        do while (setting%Time%Now <= (setting%Time%End - dtTol))

            !% --- check for blowup conditions
            call util_crashcheck (44804)
            call util_crashstop (44804)

                ! sync all
                ! do ii=1,num_images()
                !     sync all
                !     if ((this_image()==1) .and. (ii==1)) then
                !         write(6,*) ' ... beginning time loop ===================',this_image(),setting%Time%Now
                !         flush(6)
                !         !wait(6)
                !         !call sleep(1)
                !     end if
                !     sync all
                ! end do 
                ! sync all

            !% --- push the old values down the stack 
            call tl_save_previous_values()

                ! write(6,*) 'timeloop 01 ',this_image(), crashI
                ! call flush(6)

            !% store the runoff from hydrology on a hydrology step
            if (doHydrology) call tl_hydrology()

                !call util_syncwrite('timeloop 02 ',crashI)


            !% --- main hydraulics time steop
            if (doHydraulics) then       

                !% --- set a clock tick for hydraulic loop evaluation
                if (this_image()==1) then
                    call system_clock(count=cval,count_rate=crate,count_max=cmax)
                    setting%Time%WallClock%HydraulicsStart = cval
                end if 

                !% HACK -- much of the following needs to be combined into an upper level BC subroutine

                !% --- get updated boundary conditions
                !% ***** BUGCHECK -- the lateral flowrate is cumulative of BC and hydrology, 
                !% ***** so check that it is zeroed before the first BC added.
                call bc_update() 

                
                    !call util_syncwrite('    timeloop aaa ',crashI)


                if (crashI==1) then
                    call util_crashpoint (884782)
                end if
            

                !% HACK put lateral here for now -- moved from face_interpolation.
                !% set lateral to zero
                Qlateral => elemR(:,er_FlowrateLateral)
                Qlateral(:) = zeroR 
        
                !% --- add lateral inflow BC to lateral inflow accumulator
                npack   => npack_elemP(ep_BClat)
                !% --- note that thisP and thisBC must be the same size or there is something wrong
                thisP   => elemP(1:npack,ep_BClat)
                thisBC  => BC%P%BClat
                Qlateral(thisP) = Qlateral(thisP) + BC%flowRI(thisBC)  

                !print *, 'QLATERAL before hydrology '
                !print *, Qlateral(780)

                !write(6,*) '   timeloop bbb ',this_image()
                ! flush(6)

                !% --- add subcatchment inflows
                if (useHydrology) then 
                    sImage => subcatchI(:,si_runoff_P_image)
                    eIdx   => subcatchI(:,si_runoff_elemIdx)
                    !% --- using the full runoff rate for this period
                    Qrate => subcatchR(:,sr_RunoffRate_baseline)

                    do mm = 1,SWMM_N_subcatch
                        !% --- only if this image holds this node
                        !print *, mm, eIdx(mm), Qlateral(eIdx(mm)), Qrate(mm)
                        if (this_image() .eq. sImage(mm)) then
                            Qlateral(eIdx(mm)) = Qlateral(eIdx(mm)) + Qrate(mm)
                        end if
                    end do
                else
                    !% continue
                end if


                    !outstring = '    timeloop ccc '
                    !call util_syncwrite

              
                write (*,*) ' '
                write(*,"(A,f12.2,A)") 'time loop ',timeNow/3600.0, ' ===================================================================='

                !% --- perform hydraulic routing
                call tl_hydraulics()

                    !outstring = '    timeloop ddd '
                    !call util_syncwrite

                !% --- close the clock tick for hydraulic loop evaluation
                sync all
                if (this_image()==1) then
                    call system_clock(count=cval,count_rate=crate,count_max=cmax)
                    setting%Time%WallClock%HydraulicsStop = cval
                    setting%Time%WallClock%HydraulicsCumulative &
                        = setting%Time%WallClock%HydraulicsCumulative &
                        + setting%Time%WallClock%HydraulicsStop &
                        - setting%Time%WallClock%HydraulicsStart
                end if 

            end if

                !call util_syncwrite('timeloop 03 ',crashI)
                

            !% --- handle output reporting
            if (setting%Output%Report%provideYN) then 

                sync all

                !% set a time tick for output timing
                if (this_image()==1) then
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
                if (this_image()==1) then
                    call system_clock(count=cval,count_rate=crate,count_max=cmax)
                    setting%Time%WallClock%LoopOutputStop = cval
                    setting%Time%WallClock%LoopOutputCumulative &
                         = setting%Time%WallClock%LoopOutputCumulative &
                         + setting%Time%WallClock%LoopOutputStop  &
                         - setting%Time%WallClock%LoopOutputStart
                end if
            end if    

                ! write(6,*) 'timeloop 04 ',this_image(), crashI
                ! call flush(6)

            call util_crashstop(13978)

            !% --- restart the hydraulics time tick
            sync all
            if (this_image()==1) then
                call system_clock(count=cval,count_rate=crate,count_max=cmax)
                setting%Time%WallClock%HydraulicsStart = cval
            end if 

                ! write(6,*) 'timeloop 05 ',this_image(), crashI
                ! call flush(6)

            sync all
            !% ---increment the time step and counters for the next time loop
            call tl_increment_counters(doHydraulics, doHydrology)

                ! write(6,*) 'timeloop 06 ',this_image(), crashI
                ! call flush(6)
 
            ! !% --- check for a crash
            ! !if (crashYN) then
            !     if ((setting%Output%Report%provideYN) .and. &
            !     (.not. setting%Output%Report%suppress_MultiLevel_Output)) then
            !         call outputML_store_data (.true.)
            !     end if
            !     exit ! leave the time-step do loop
            ! end if

            !% --- close the hydraulics time tick
            sync all
            if (this_image()==1) then
                call system_clock(count=cval,count_rate=crate,count_max=cmax)
                setting%Time%WallClock%HydraulicsStop = cval
                setting%Time%WallClock%HydraulicsCumulative &
                    = setting%Time%WallClock%HydraulicsCumulative &
                    + setting%Time%WallClock%HydraulicsStop &
                    - setting%Time%WallClock%HydraulicsStart
            end if

 
                ! sync all
                ! do ii=1,num_images()
                !     sync all
                !     flush(6)
                !     if ((this_image()==1) .and. (ii==1)) then 
                !         flush(6)
                !         write(6,*) ' ... end of time loop ======================',this_image(), crashI
                !         flush(6)
                !         !wait(6)
                !         !call sleep(1)
                !     else
                !         flush(6)
                !     end if
                !     sync all
                ! end do
                ! sync all


            call util_crashstop (773623)


        end do  !% end of time loop

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
        !% execute the SWMM-C runoff for the next interval 
        call interface_call_runoff_execute()

        !% get the next runoff time
        setting%Time%Hydrology%NextTime = interface_get_NewRunoffTime()   

        !% cycle through the subcatchments to get the runoff 
        !% ii-1 required in arg as C arrays start from 0
        do ii = 1,SWMM_N_subcatch
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
        if (doHydraulics) then
            call tl_update_hydraulics_timestep()
            call util_crashstop(449873)
        else
            nextHydraulicsTime = setting%Time%End + tenR*DtTol
        end if    

        !% --- The NextHydrologyTime is updated in SWMM, here we just need to
        !%     provide a large number if hydrology isn't used
        if (.not. useHydrology) then
            nextHydrologyTime = setting%Time%End + tenR*DtTol
        else
            !% --- SWMM-C will not return a next hydrology time that is
            !%     beyond the end of the simulation, so we need a special
            !%     treatment for this case
            if ((nextHydrologyTime .eq. lastHydrologyTime) .and. &
                (timeNow > setting%Time%Start)) then 
                nextHydrologyTime = setting%Time%End + tenR*DtTol
            end if
        end if

        !% --- Ensure that all processors use the same time step.
        !% --- find the minimum hydraulics time and store accross all processors
        call co_min(nextHydraulicsTime)
        !% --- take the nextTime as the minimum of either the Hydrology or Hydraulics time
        !%     done on a single processor because they all should have identical nextHydrologyTIme
        nextTime = min(nextHydraulicsTime, nextHydrologyTime)

        !% --- note there is no need to broadcast across processors since co_min
        !%     ensures that each processor gets the min value.

        !% --- update the time step for the local processor (all have the same nextTime and lastHydraulicsTime)
        dtHydraulics = nextTime - lastHydraulicsTime

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

        !% find the cfl for reporting
        cfl_max = tl_get_max_cfl(ep_CC_Q_NOTsmalldepth,dtHydraulics)
        call co_max(cfl_max)

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
    subroutine tl_update_hydraulics_timestep()
        !%------------------------------------------------------------------
        !% Description:
        !% updates the timestep (dt) for hydraulics and computes the
        !% setting.Time.Hydraulics.NextTime for the current processor
        !%------------------------------------------------------------------
        !% Declarations
            logical, pointer :: matchHydrologyStep, useHydrology
            real(8)          :: oldDT, maxVelocity
            real(8)          :: timeleft, thisCFL
            real(8), pointer :: targetCFL, maxCFL, maxCFLlow, timeNow, dtTol
            real(8), pointer :: increaseFactor
            real(8), pointer :: newDT
            !rm velocity(:), wavespeed(:), length(:), PCelerity(:)
            real(8), pointer :: nextHydrologyTime, nextHydraulicsTime
            real(8), pointer :: lastHydrologyTime, lastHydraulicsTime
            integer          :: ii, neededSteps
            integer, pointer :: stepNow, stepNext, stepfinal, lastCheckStep, checkStepInterval
            !rm integer, pointer :: thisCol, Npack, thisP(:)
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
            !rm 20220209brh decreaseFactor     => setting%VariableDT%decreaseFactor
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
        
        ! print *, '=============================================='
        ! print *, 'oldDT ',oldDT
        ! print *, 'newDT ',newDT

        if ((matchHydrologyStep) .and. (useHydrology)) then 
            !% --- for combined hydrology and hydraulics compute the CFL if we take a single
            !%     step to match the next hydrology time
            timeLeft = nextHydrologyTime - lastHydraulicsTime
            if (timeLeft .le. dtTol) timeLeft = oldDT
            thisCFL = tl_get_max_cfl(ep_CC_Q_NOTsmalldepth,timeleft)  

            ! print *, 'this CFL 01 ',thisCFL

            !% --- check to see if a single time step to match the hydrology time is possible
            if (thisCFL < maxCFL) then
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

                    ! print *, 'newDT 02   ', newDT
                else
                    !% --- small time left, don't bother with it, go back to the old dt
                    newDT = oldDT
                    !% --- check that resetting to oldDT didn't cause a problem
                    thisCFL = tl_get_max_cfl(ep_CC_Q_NOTsmalldepth,newDT)
                    if (thisCFL > maxCFL) then 
                        !% --- if CFL to large, set the time step based on the target CFL
                        newDT = newDT * targetCFL / thisCFL
                    end if

                    ! print *, 'this CFL 03 ',thisCFL
                    ! print *, 'newDT 03    ',newDT
                end if
            else
                !% --- if more than one time step is needed, compute the needed steps 
                !%     to break up the large CFL and time step size.
                !%     First check to see if the implied time step is too small for the integer size
                if (thisCFL/targetCFL .ge. huge(neededSteps)) then
                    write(*,*) 'warning -- really high CFL, setting dt to minimum to prevent overflow'
                    newDT = setting%Limiter%Dt%Minimum + setting%Time%DtTol
                    neededSteps = 1000

                    ! print *, 'newDT 04    ',newDT
                else
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
                    thisCFL = tl_get_max_cfl(ep_CC_Q_NOTsmalldepth,newDT)
                    if (thisCFL > maxCFL) then 
                        !% --- if CFL to large, set the time step based on the target CFL
                        newDT = newDT * targetCFL / thisCFL
                    end if

                    ! print *, 'this CFL 05 ',thisCFL
                    ! print *, 'newDT 05    ',newDT
                end if
            end if
        else
            neededSteps = 3 !% used to control rounding
            !% --- allowing hydrology and hydraulics to occur at different times
            thisCFL = tl_get_max_cfl(ep_CC_Q_NOTsmalldepth,oldDT)

            ! print *, 'thisCFL 06 ',thisCFL

            if (thisCFL > maxCFL) then
                !% --- decrease the time step to the target CFL and reset the checkStep counter
                newDT = oldDT * targetCFL / thisCFL
                lastCheckStep = stepNow

                ! print *, 'newDT 07    ',newDT
            else
                !% --- if CFL is less than max, see if it can be raised (only checked at intervals)
                if (stepNow > lastCheckStep + checkStepInterval) then
                    !% --- check for low CFL only on prescribed intervals and increase time step
                    if (thisCFL < maxCFLlow) then
                        !% --- increase the time step and reset the checkStep Counter
                        !%newDT = oldDT * increaseFactor
                        newDT = OldDT * targetCFL / thisCFL  ! 20220214brh
                        lastCheckStep = stepNow
                    end if

                    ! print *, 'thisCFL 08 ',thisCFL
                    ! print *, 'newDT   08 ',newDT
                end if
            end if
        end if

        !% 20220328brh time step limiter for inflows into small or zero volumes
        !print *, 'dt before limit ', newDT
        call tl_limit_BCinflow_dt (newDT)
        ! print *, 'dt BC limit     ',newDT

        call tl_limit_LatInflow_dt (newDT)
        ! print *, 'dt Qlat limit   ',newDt

        !% --- if dt is large and there is more than 2 steps, then round to an integer number
        if ((matchHydrologyStep) .and. (useHydrology) .and. (neededSteps .le. 2)) then
            !% don't round the dt
        else
            if (newDT > fiveR) then
                !% --- round larger dt to counting numbers
                newDT = real(floor(newDT),8)

                ! print *, 'newDT 10    ',newDT
            elseif (newDT > oneR) then
                !% --- round smaller dt to two decimal places
                newDT = real(floor(newDT * onehundredR),8) / onehundredR

                ! print *, 'newDT 11    ',newDT   
            else
                !%  HACK -- should round to 3 places for smaller numbers
            end if
        end if

        !% increment the hydraulics time clock
        nextHydraulicsTime = lastHydraulicsTime + newDT

        !%----------------------------------------------------------------------
        !% closing
            if ((setting%Limiter%Dt%UseLimitMinYN) .and. (newDT .le. setting%Limiter%Dt%Minimum)) then
                print*, 'timeNow = ', timeNow
                print*, 'dt = ', newDT, 'minDt = ',  setting%Limiter%Dt%Minimum
                print*, 'max velocity  ', maxval( &
                    elemR(elemP(1:npack_elemP(ep_CC_Q_NOTsmalldepth),ep_CC_Q_NOTsmalldepth),er_Velocity) )
                print*, 'max wavespeed ', maxval( &
                    elemR(elemP(1:npack_elemP(ep_CC_Q_NOTsmalldepth),ep_CC_Q_NOTsmalldepth),er_WaveSpeed) )
                print*, 'warning: the dt value is smaller than the user supplied min dt value'
                !stop 1123938
                call util_crashpoint(1123938)
            end if

            ! print*, 'timeNow = ', timeNow
            ! print*, 'dt = ', newDT, 'minDt = ',  setting%Limiter%Dt%Minimum
            ! maxVelocity = maxval( &
            !     elemR(elemP(1:npack_elemP(ep_CC_Q_NOTsmalldepth),ep_CC_Q_NOTsmalldepth),er_Velocity) )
            ! print*, 'max velocity  ', maxVelocity
            ! print*, 'max V locate  ', findloc(elemR(:,er_Velocity),maxVelocity)
            ! print*, 'max wavespeed ', maxval( &
            !     elemR(elemP(1:npack_elemP(ep_CC_Q_NOTsmalldepth),ep_CC_Q_NOTsmalldepth),er_WaveSpeed) )

            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine tl_update_hydraulics_timestep
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
            integer, pointer :: step, interval
            integer (kind=8) :: execution_realtime
            integer(kind=8) :: cval, crate, cmax
            real(8) :: simulation_fraction, seconds_to_completion, time_to_completion
            character(8) :: timeunit
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
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
            
        !rint *,'++++++++++++++++++++++'
        !print *, (setting%Time%WallClock%ClockNow - setting%Time%WallClock%ClockLoopStart ) / setting%Time%WallClock%ClockCountRate  
        !print *,'++++++++++++++++++++++'                                         
        !print *, 'setting%Time%WallClock%EpochTimeLoopStartSeconds ',setting%Time%WallClock%ClockLoopStart     
        !print *, 'setting%Time%WallClock%EpochNowSeconds           ',setting%Time%WallClock%ClockNow                                  
        !print *, 'execution_realtime ',execution_realtime    
        !print *, 'seconds_to_completion ',  seconds_to_completion      
        !stop 487566

        if (setting%Output%Verbose) then
            if (this_image() == 1) then
                if (mod(step,interval) == 0) then
                    thistime = timeNow
                    call util_datetime_display_time (thistime, timeunit)
                    ! ! translate time in seconds into something more useful
                    ! if (timeNow  < sixtyR) then
                    !     thistime = timeNow
                    !     timeunit = 's  '
                    ! elseif (timeNow >= sixtyR .and. timeNow < seconds_per_hour) then
                    !     thistime = timeNow / sixtyR
                    !     timeunit = 'min'
                    ! elseif (timeNow >= seconds_per_hour .and. timeNow < 3.0*seconds_per_day) then
                    !     thistime = timeNow / seconds_per_hour
                    !     timeunit = 'hr '
                    ! elseif (timeNow >= 3.0 * seconds_per_day) then
                    !     thistime = timeNow / seconds_per_day
                    !     timeunit = 'day'
                    ! endif

                    ! write a time counter
                    write(*,"(A12,i8,a17,F9.2,a1,a8,a6,f9.2,a3,a8,f9.2)") &
                        'time step = ',step,'; model time = ',thistime, &
                        ' ',trim(timeunit),'; dt = ',dt,' s', '; cfl = ',cfl_max

                    ! write estimate of time remaining
                    thistime = seconds_to_completion
                    call util_datetime_display_time (thistime, timeunit)
                    ! if (seconds_to_completion < sixtyR) then
                    !     timeunit = 's  '
                    !     time_to_completion = seconds_to_completion
                    ! elseif (seconds_to_completion >=sixtyR .and. seconds_to_completion < seconds_per_hour ) then
                    !     timeunit = 'min'
                    !     time_to_completion = seconds_to_completion / sixtyR
                    ! elseif (seconds_to_completion >=seconds_per_hour .and. seconds_to_completion < seconds_per_day) then
                    !     timeunit = 'hr '
                    !     time_to_completion = seconds_to_completion / seconds_per_hour
                    ! else
                    !     timeunit = 'day'
                    !     time_to_completion = seconds_to_completion / seconds_per_day
                    ! endif
                    write(*,"(A9,F6.2,A1,A3,A)") 'estimate ',thistime,' ',timeunit,' wall clock time until completion'
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

        !print *, 'in tl_get_max_CFL'
        !print *, size(thisP), Npack

        if (Npack > 0) then 
            outvalue = max (maxval((abs(velocity(thisP)) + abs(wavespeed(thisP))) * thisDT / length(thisP)), &
                            maxval((abs(PCelerity(thisP))) * thisDT / length(thisP)))
        !    print *, 'max vel, wave ', maxval(abs(velocity(thisP))), maxval(abs(wavespeed(thisP))) 
        !    itemp = maxloc(abs(velocity(thisP))) 
        !    ip = thisP(itemp(1))
        !    fup = elemI(ip,ei_Mface_uL)
        !    fdn = elemI(ip,ei_Mface_dL)
        !    eup = faceI(fup,fi_Melem_uL)
        !    edn = faceI(fdn,fi_Melem_dL)
           !print *, ' at max loc '
           !print *, trim(reverseKey(elemI(ip,ei_elementType))), ' ', trim(reverseKey(elemI(ip,ei_geometryType)))
           !print *, trim(reverseKey(elemI(eup,ei_elementType))), ' ', trim(reverseKey(elemI(eup,ei_geometryType)))
           !print *, trim(reverseKey(elemI(edn,ei_elementType))), ' ', trim(reverseKey(elemI(edn,ei_geometryType)))
           !print *, velocity(ip), wavespeed(ip), Pcelerity(ip)
           !print *, elemR(ip,er_Depth), elemR(ip,er_Area), elemR(ip,er_Length)
           !print *, faceR(fup,fr_Flowrate), faceR(fdn,fr_Flowrate)
           
        !    print *, elemR(eup,er_Head), elemR(ip,er_Head), elemR(edn,er_Head)
        !    print *, elemR(eup,er_Zbottom), elemR(ip,er_Zbottom), elemR(ip,er_Zbottom)
        !    print *, elemR(eup,er_Depth), elemR(ip,er_Depth), elemR(edn,er_Depth)
        !    print *, elemI(eup,ei_link_Gidx_SWMM), elemI(ip,ei_link_Gidx_SWMM), elemI(edn,ei_link_Gidx_SWMM)
        !    print *, eup, ip, edn
        !    stop 29873

                     
        else
            outvalue = zeroR
        end if

        ! print * , ' '
        ! print *, velocity(thisP)
        ! print *, ' '
        ! print *, wavespeed(thisP)
        ! print *, ' '
        ! print *, elemR(thisP,er_Area)
        ! print *, ' '
        ! print *, elemR(thisP,er_Head)
        ! print *, ' '
        ! print *, elemR(thisP,er_WaveSpeed)
        ! print *, ' '
        ! print *, elemR(thisP,er_Preissmann_Celerity)
        ! print *, ' '
        ! print *, outvalue
        ! print *, thisDT
        ! print *, 'exiting tl_get_max_cfl'
        ! !stop 39875

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
        
        !% set the inflow limiter for upstream boundary conditions (BCup)
        if (size(BC%P%BCup) > 0) then

            !% ensure flowrate used for limiter is  positive and non-zero
            !tinyQ(:) = setting%Eps%Velocity
            BCQ(:) = max( abs(BC%FlowRI(BC%P%BCup)), setting%Eps%Velocity)

            !print *, 'BCQ', BCQ

            !% store the element indexes downstream of a BCup face
            BCelemIdx =  faceI(BC%flowI(BC%P%BCup,bi_face_idx), fi_Melem_dL)

            !% store the scaling depth limit defined by when the induced velocity
            !% from an inflow is similar to the gravity wave speed of the BC inflow
            depthScale = ( (onehalfR**4) * (BCQ**2) / (gravity * (topwidth(BCelemIdx)**2) ) )**onethirdR

            !print *, 'depthScale ',depthScale

            !% get the depth limit for maximum depth that the time step will be limited as 
            !% either the depth scale or the specified smallDepth limit
            depthLimit = max(smallDepth, depthScale)

            !print *, 'depthLimit ',depthLimit
            
            !% time step limit based on inflow flux
            DTlimit = length(BCelemIdx) * ( ( alpha * topwidth(BCelemIdx) / (gravity * BCQ) )**onethirdR) 

            !print *, 'DTlimit ',DTlimit
        
            !% where the depth is greater than the depthlimit the DT inflow limiter
            !% is not needed, and we can use the existing DT value
            depthLimit = depth(BCelemIdx) - depthLimit
            where (depthLimit .ge. zeroR)
                DTlimit = thisDT
            endwhere

            !print *, 'new DTlimit ',DTlimit

            !% get the smallest DT in the limiter array
            newDTlimit = minval(DTlimit)

            !% use the smaller value of the new limit or the input
            thisDT = min(newDTlimit,thisDT)

            !print *, 'thisDT ',thisDT

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
!% END OF MODULE
!%+=========================================================================
end module timeloop

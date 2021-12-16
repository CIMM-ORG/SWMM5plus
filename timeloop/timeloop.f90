module timeloop

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use output, only: outputML_store_data
    use pack_mask_arrays
    use runge_kutta2
    use hydrology
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
        !%--------------------------------------------------------------------
        !% Preliminaries
            if (icrash) return
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

            if (setting%Output%Verbose) &
                write(*,"(2A,i5,A)") new_line(" "), 'begin timeloop [Processor ', this_image(), "]"
            call system_clock(setting%Time%Real%ClockLoopStart)
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

        print *, 'here 37685 ',dtHydraulics

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

        print *, 'here 37686 ',dtHydraulics

        !% get the initial dt and the next hydraulics time
        if (useHydraulics) then
            call tl_update_hydraulics_timestep()
            !print *, dtHydraulics
            !stop 398705
        else
            !% set a large dummy time for hydraulics if not used
            nextHydraulicsTime = timeEnd + onethousandR * dtTol
            !% suppress the multi-level hydraulics output (otherwise seg faults)
            setting%Output%Report%suppress_MultiLevel_Output = .true.
        end if

        print *, 'here 3875 ',dtHydraulics

        !% check to see if there is an initial step to hydrology
        !% if not, then skip the initial tl_hydrology
        if (( abs(nextHydrologyTime - timeNow) < dtTol) .and. (doHydrology) ) then
            doHydrology = .true.
        else
            doHydrology = .false.
        endif 

        print *, 'here 9875 ',dtHydraulics

        !% Combined hydrology (SWMM-C) and hydraulics simulation
        !% The loop starts at t = setting%Time%Start
        !% Adust the end time for the dtTol (precision error in epoch to seconds conversion)
        ii=0
        !do ii=1,10
        do while (setting%Time%Now <= (setting%Time%End - dtTol))
            ii = ii+1
            !print *, '----------------------'
            !print *, ii, timeNow, setting%Time%End, doHydrology, doHydraulics
            !print *

            !% store the runoff from hydrology on a hydrology step
            if (doHydrology) call tl_hydrology()

            if (doHydraulics) then        

                print *, 'in 870533 ',trim(subroutine_name)
                write(*,'(8f9.2)') faceR(1:N_elem(1)+1,fr_Flowrate)
                !% get updated boundary conditions
                !% ***** BUGCHECK -- the lateral flowrate is cumulative of BC and hydrology, 
                !% ***** so check that it is zeroed before the first BC added.
                call bc_update()

                print *, 'in 879533 ',trim(subroutine_name)
                write(*,'(8f9.2)') faceR(1:N_elem(1)+1,fr_Flowrate)

                !% HACK put lateral here for now -- cut out of face_interpolation.
                !% set lateral to zero
                Qlateral => elemR(:,er_FlowrateLateral)
                Qlateral(:) = zeroR 
        
                !% add inflow BC to lateral inflow
                npack   => npack_elemP(ep_BClat)
                !% note that thisP and thisBC must be the same size or there is something wrong
                thisP   => elemP(1:npack,ep_BClat)
                thisBC  => BC%P%BClat
                Qlateral(thisP) = Qlateral(thisP) + BC%flowRI(thisBC)  

                print *, 'lateral inflows'
                print *, Qlateral(thisP)

                !% add subcatchment inflows
                if (useHydrology) then 
                    sImage => subcatchI(:,si_runoff_P_image)
                    eIdx   => subcatchI(:,si_runoff_elemIdx)
                    !% using the full runoff rate for this period
                    Qrate => subcatchR(:,sr_RunoffRate_baseline)

                    !print *, 'all of Qrate ', Qrate(:)

                    do mm = 1,SWMM_N_subcatch
                        !% only if this image holds this node
                        if (this_image() .eq. sImage(mm)) then
                            Qlateral(eIdx(mm)) = Qlateral(eIdx(mm)) + Qrate(mm)
                        end if
                    end do
                else
                    !% continue
                end if

                print *, 'here 22355 ',dtHydraulics

                print *, 'in 771552 ',trim(subroutine_name)
                write(*,'(8f9.2)') faceR(1:N_elem(1)+1,fr_Flowrate)

                !% perform hydraulic routing
                call tl_hydraulics()

                print *, 'in 93705 ',trim(subroutine_name)
                write(*,'(8f9.2)') faceR(1:N_elem(1)+1,fr_Flowrate)
                !do kk=1,size(elemI,DIM=1)
                !    write(*,'(i4, 3f8.3)') kk,elemR(kk,er_FlowrateLateral), elemR(kk,er_Flowrate), elemR(kk,er_Depth)
                !end do

                ! kk=24
                
                ! write(*,'(A,i4, 5f8.5)') 'Head    ', kk, &
                !      elemR(kk-1,er_Head)-1514.0, &
                !      elemR(kk  ,er_Head)-1514.0, &
                !      elemR(kk+1,er_Head)-1514.0
                ! write(*,'(A,i4, 5f8.5)') 'depth   ', kk, &
                !      elemR(kk-1,er_Depth), &
                !      elemR(kk  ,er_Depth), &
                !      elemR(kk+1,er_Depth)
                ! write(*,'(A,i4, 5f8.5)') 'Flowrate', kk, &
                !      elemR(kk-1,er_Flowrate), &
                !      elemR(kk  ,er_Flowrate), &
                !      elemR(kk+1,er_Flowrate)

                ! write(*,'(A,i4, 5f12.5)') 'FaceQ  ', kk, &
                !      faceR(faceI(kk,fi_Melem_uL),fr_Flowrate), &
                !      faceR(faceI(kk,fi_Melem_dL),fr_Flowrate)
                ! !write(*,'(A,i4, 3f8.3)') 'lateral', kk, elemR(kk,er_FlowrateLateral), elemR(kk+1,er_FlowrateLateral)
                ! !write(*,'(A,i4, 3f8.3)') 'flow   ', kk, elemR(kk,er_Flowrate), elemR(kk+1,er_Flowrate)
                ! write(*,*) faceI(kk,fi_Melem_uL), kk, faceI(kk,fi_Melem_dL)
            end if

            if (setting%Output%Report%provideYN) then 
                !% Results must be reported before the "do"counter increments
                !call util_output_report()  --- this is a stub for future use

                !% Multilevel time step output
                if ((util_output_must_report()) .and. &
                    (.not. setting%Output%Report%suppress_MultiLevel_Output) ) then
                    call outputML_store_data (.false.)
                end if
            end if    

            !% increment the time step and counters for the next time loop
            call tl_increment_counters(doHydraulics, doHydrology)
            if (icrash) then
                if ((setting%Output%Report%provideYN) .and. &
                (.not. setting%Output%Report%suppress_MultiLevel_Output)) then
                    call outputML_store_data (.true.)
                end if
                exit
            end if
        end do

        !% >>> BEGIN HACK
        !%     Temporary for debugging (can be deleted for deployment)
        ! if (setting%Debug%OutputYN) then
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
        !%------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step
        !%------------------------------------------------------------------
        !% Declarations
            integer :: ii
            real(8), pointer :: sRunoff(:)
            character(64) :: subroutine_name = 'tl_hydrology'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (icrash) return
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
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
            integer, allocatable :: tempP(:)
            character(64)    :: subroutine_name = 'tl_hydraulics'
        !%-------------------------------------------------------------------
            if (icrash) return
            if (setting%Debug%File%timeloop) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------

        print *, 'in 77155 ',trim(subroutine_name)
        write(*,'(8f9.2)') faceR(1:N_elem(1)+1,fr_Flowrate)

        print *, 'here aaa75',setting%Time%Hydraulics%Dt
        tempP = pack(elemI(:,ei_Lidx),elemI(:,ei_elementType)== CC)
        print *, elemR(tempP,er_Velocity)
        deallocate(tempP)
        print *, ' '

        !% check for where solver needs to switch in dual-solver model
        if (setting%Solver%SolverSelect == ETM_AC) then
            call tl_solver_select()
        end if

        print *, 'here bbb75',setting%Time%Hydraulics%Dt

        print *, 'in 8701875 ',trim(subroutine_name)
        write(*,'(8f9.2)') faceR(1:N_elem(1)+1,fr_Flowrate)

        !% repack all the dynamic arrays
        !% FUTURE 20210609 brh need to decide where this goes
        call pack_dynamic_arrays()
        ! print *, "Need to decide on pack_dynamic_arrays 94837"

        print *, 'here ccc75',setting%Time%Hydraulics%Dt

        print *, 'in 3705 ',trim(subroutine_name)
        write(*,'(8f9.2)') faceR(1:N_elem(1)+1,fr_Flowrate)

        !%  push the old values down the stack for AC solver
        call tl_save_previous_values()

        print *, 'here ddd75',setting%Time%Hydraulics%Dt

        print *, 'in 9775 ',trim(subroutine_name)
        write(*,'(8f9.2)') faceR(1:N_elem(1)+1,fr_Flowrate)

        !%  Reset the flowrate adhoc detection before flowrates are updated.
        !%  Note that we do not reset the small volume detection here -- that should
        !%  be in geometry routines.
        elemYN(:,eYN_IsAdhocFlowrate) = .false.

        print *,"****************** velocity before"
        print *, elemR(:,er_Velocity)

        print *, 'in 870955 ',trim(subroutine_name)
        write(*,'(8f9.2)') faceR(1:N_elem(1)+1,fr_Flowrate)

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

        print *,"****************** velocity after "
        print *, elemR(:,er_Velocity)
        !tempP = pack(elemI(:,ei_Lidx),elemI(:,ei_elementType)== CC)
        !print *, elemR(tempP,er_Velocity)
        !deallocate(tempP)
        print *, ' '

        print *, 'in 38705 ',trim(subroutine_name)
        write(*,'(8f9.2)') faceR(1:N_elem(1)+1,fr_Flowrate)

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

        if ((matchHydrologyStep) .and. (useHydrology)) then 
            !% For combined hydrology and hydraulics compute the CFL if we take a single
            !% step to match the next hydrology time
            timeLeft = nextHydrologyTime - lastHydraulicsTime
            if (timeLeft .le. dtTol) timeLeft = oldDT
            thisCFL = tl_get_max_cfl(ep_CC_Q_NOTsmallvolume,timeleft)  

            !% check to see if a single time step to match the hydrology time is possible
            if (thisCFL < maxCFL) then
                neededSteps = 1
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
                !% compute the needed steps to break up the large CFL and time step size
                !% first check to see if the implied time step is too small for the integer size
                if (thisCFL/targetCFL .ge. huge(neededSteps)) then
                    write(*,*) 'warning -- really high CFL, setting dt to minimum to prevent overflow'
                    newDT = setting%Limiter%Dt%Minimum + setting%Time%DtTol
                    neededSteps = 1000
                else
                    !% note that neededSteps will be 2 or larger else thisCFL < maxCFL
                    neededSteps = ceiling( thisCFL / targetCFL )
                    !print *, 'neededSteps ', neededSteps
                    !% the provisional time step that would get exactly to the hydrology time (if CFL didn't change)
                    newDT = timeleft / real(neededSteps,8)
                    !% limit the change in the dt by the increase factor
                    if (newDT / oldDT > increaseFactor) then 
                        newDT = oldDT * increaseFactor
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
            neededSteps = 3 !% used to control rounding
            !% Allowing hydrology and hydraulics to occur at different times
            thisCFL = tl_get_max_cfl(ep_CC_Q_NOTsmallvolume,oldDT)

            !print *,'here 111 ',thisCFL
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
        
        !print *, matchHydrologyStep, useHydrology, neededSteps

        !% if dt is large and there is more than 2 steps, then round to an integer number
        if ((matchHydrologyStep) .and. (useHydrology) .and. (neededSteps .le. 2)) then
            !% don't round the dt
        else
            !print *, 'new dt ', newDT
            if (newDT > fiveR) then
                !% round larger dt to counting numbers
                newDT = real(floor(newDT),8)
            elseif (newDT > oneR) then
                newDT = real(floor(newDT * onehundredR),8) / onehundredR
            else
                !%  HACK NEED SOMETHING TO ROUND TO 3 DIGITS FOR SMALLER NUMBERS
            end if
        end if

        !print *, ' '
        !print *, 'in here 397057 dt ',newDT
        !print *, ' '
        !% increment the hydraulics time clock
        nextHydraulicsTime = lastHydraulicsTime + newDT

        if ((setting%Limiter%Dt%UseLimitMinYN) .and. (newDT .le. setting%Limiter%Dt%Minimum)) then
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
            reportStep    => setting%Output%Report%ThisStep
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
            !% SWMM-C will not return a next hydrology time that is
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
            integer (kind=8) :: execution_realtime
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

        call system_clock(count=setting%Time%Real%ClockNow) ! Fortran function returns real time

        ! estimate the remaining time
        execution_realtime = (setting%Time%Real%ClockNow - setting%Time%Real%ClockLoopStart)
        seconds_to_completion = real(execution_realtime,kind=8) * (setting%Time%End - setting%Time%Now) &
                                                   / (setting%Time%Now - setting%Time%Start)

        seconds_to_completion = execution_realtime / real(setting%Time%Real%ClockCountRate) 
                                                   
        !print *, 'setting%Time%Real%EpochTimeLoopStartSeconds ',setting%Time%Real%ClockLoopStart     
        !print *, 'setting%Time%Real%EpochNowSeconds           ',setting%Time%Real%ClockNow                                  
        !print *, 'execution_realtime ',execution_realtime    
        !print *, 'seconds_to_completion ',  seconds_to_completion      
        !stop 487566

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
                    write(*,"(A12,i8,a17,F9.2,a1,a3,a6,f9.2,a3)") &
                        'time step = ',step,'; model time = ',thistime, &
                        ' ',timeunit,'; dt = ',dt,' s'

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
            real(8), intent(in), target :: dt
            real(8), pointer    :: thisDT
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
        if (dt .le. zeroR) then 
            thisDT => setting%Time%Hydraulics%Dt
        else    
            thisDT => dt
        end if
        
        outvalue = max (maxval((abs(velocity(thisP)) + abs(wavespeed(thisP))) * thisDT / length(thisP)), &
                        maxval((abs(velocity(thisP)) + abs(PCelerity(thisP))) * thisDT / length(thisP)))

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
!% END OF MODULE
!%+=========================================================================
end module timeloop

! module time_loop
!
! top-level iteration for continuity and momentum solutions
!
!==========================================================================
!
module time_loop

    use ac_convergence_loop
    use array_index
    use bc
    use control
    use data_keys
    use debug
    use diagnostic
    use explicit_euler
    use face_values
    use globals
    use output
    use runge_kutta
    use setting_definition
    use utility

    implicit none

    private

    public  :: time_marching

    integer :: debuglevel = 0

    integer,   parameter :: idummy = 0

contains
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine time_marching &
        (elem2R, elemMR, faceR, elem2I, elemMI, faceI, elem2YN, elemMYN, &
        faceYN, bcdataDn, bcdataUp,  gateSetting, linkI, nodeI, linkR,   &
        nodeR, debugfile, diagnostic, threadedfile, ID, numberPairs,     &
        ManningsN, Length, zBottom, xDistance, Breadth, widthDepthData,  &
        cellType)

        !
        ! top-level iteration for continuity and momentum solution
        !
        character(64) :: subroutine_name = 'time_marching'

        real(8),      target, intent(in out) :: elem2R(:,:),  elemMR(:,:),  faceR(:,:)
        integer,   target, intent(in out) :: elem2I(:,:),  elemMI(:,:),  faceI(:,:)
        logical,   target, intent(in out) :: elem2YN(:,:), elemMYN(:,:), faceYN(:,:)

        type(bcType),                  intent(in out) :: bcdataDn(:), bcdataUp(:)
        type(controlType),             intent(in out) :: gateSetting(:)
        type(debugfileType),  target,  intent(in)     :: debugfile(:)
        type(diagnosticType), target,  intent(in out) :: diagnostic(:)
        type(threadedfileType),        intent(in)     :: threadedfile(:)

        integer,                       intent(in)     :: linkI(:,:), nodeI(:,:)
        real(8),                          intent(in out) :: linkR(:,:), nodeR(:,:)

        real(8), pointer :: rkVol(:), rkU(:)
        real(8), pointer :: fQ(:), fUdn(:), fUup(:), fAdn(:), fAup(:)
        real(8), pointer :: fEdn(:), fEup(:), eE(:)

        real(8)    :: AnormH(3), AnormQ(3), AnormHlast(3), AnormQlast(3)
        real(8)    :: TnormH(3), TnormQ(3), RTnormH(3), RTnormQ(3), RLnormH(3), RLnormQ(3)

        real(8), pointer :: thistime, nexttime

        integer, pointer :: thisStep, restartStep

        logical :: realLoop
        integer :: ii

        character(len=32) :: outdataName

        integer, intent(in out)    :: ID(:)
        integer, intent(in out)    :: numberPairs(:)
        real(8),    intent(in out)    :: ManningsN(:)
        real(8),    intent(in out)    :: Length(:)
        real(8),    intent(in out)    :: zBottom(:)
        real(8),    intent(in out)    :: xDistance(:)
        real(8),    intent(in out)    :: Breadth(:)
        real(8),    intent(in out)    :: widthDepthData(:,:,:)
        type(string), intent(in out)   :: cellType(:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        thisStep    => setting%Step%Current
        !%  Note that for a restart, the simulation may begin from a step other than 1
        !%  In which case, diagnostic storage indexed off of the step are using the
        !%  step "FromRestart" index
        restartStep => setting%Step%FromRestart
        thistime    => setting%Time%ThisTime
        nexttime    => setting%Time%NextTime

        !% get the upper time
        nexttime = thistime + dt

        !%  initialization of local storage the Lnorms for head and flowrate convergence
        !%  these ares stored as L1, L2, Linf
        !%  Note that these are different than the trun. norms that are used
        !%  for steady-state checking
        !%  Absolute norms are the dimensional rates of change dpsi/dtau
        AnormH = zeroR
        AnormQ = zeroR  
        AnormHlast = zeroR
        AnormQlast = zeroR   
        
        !%  Total norm is the dimensional rate of change from the N0 values
        TnormH = zeroR
        TnormQ = zeroR      
        
        !%  Relative norms are the relative rates of change. 
        !%  The RT norm is the recent norm compared to the norm for total rate 
        !%  of change since N0
        RTnormH = zeroR
        RTnormQ = zeroR     
                
        !%  Relative last norm is the recen Anorm compared to AnormLast
        RLnormH = zeroR
        RLnormQ = zeroR
        

        call bc_applied_onelement &
            (elem2R, bcdataDn, bcdataUp, thistime, bc_category_inflowrate, e2r_Velocity)

        call bc_applied_onelement &
            (elem2R, bcdataDn, bcdataUp, thistime, bc_category_elevation, idummy)

        call bc_applied_onface (faceR, faceI, elem2R, elem2I, bcdataDn, bcdataUp, e2r_Velocity, thistime)
                  
        call diagnostic_volume_conservation &
            (diagnostic, elem2R, elem2I, elemMR, elemMI, faceR, bcdataUp, bcdataDn, restartStep,  1)
            
        restartStep = restartStep + 1   ! HACK required for diagnostics so end of first step is index 2

        call  output_all_threaded_data_by_link &
            (threadedfile, elem2R, elem2I, elemMR, elemMI, faceR, faceI, linkI, &
            bcdataUp, bcdataDn, 0)
            
        
        !% Iterate for a fixed number of steps
        !% HACK - need to develop better time and iteration controls
        do while (thisstep <= setting%step%final)

            !% display to the screen
            if (setting%Debugout%DisplayInterval > 0) then
                if (mod(thisstep,setting%Debugout%DisplayInterval) == 0) then
                    print *, '--------------------------------------'
                    print *, thisstep,'=current step; ', thistime, '=this time', &
                        diagnostic_CFL(elem2R, e2r_Timescale_Q_u, e2r_Timescale_Q_d),'=CFL max' & !Im putting the timescale for Q here. Might needed to be changed later
                        ,maxval(abs(elem2R(2:size(elem2R,1)-1,e2r_Velocity))),'=max velocity'
                    !print *, thisstep,'=current step; ', &
                    !         maxval(elem2R(2:size(elem2R,1)-1,e2r_Flowrate))

                endif
            endif

            !%  Reset the flowrate adhoc detection before flowrates are updated.
            !%  Note that we do not reset the small volume detection here - that is done
            !%  in the element_geometry module before the volumes are updated.
            elem2YN(:,e2YN_IsAdhocFlowrate) = .false.
            elemMYN(:,eMYN_IsAdhocFlowrate) = .false.

            !%  push the old values down the stack for AC solver
            call save_previous_values (elem2R, elemMR, faceR)

            !% Runge-Kutta 2nd-order advance
            !% realLoop mask ensures steptime is advanced for timeloop
            realLoop = .true.

            !% HACK: only works for trajkovic cases
            ! call control_evaluate &
            !     (elem2I, elem2R, gateSetting, N_Gates, thistime, .false.)

            call rk2 &
                (elem2R, elemMR, elem2I, elemMI, faceR, faceI, elem2YN, elemMYN, faceYN, &
                bcdataDn, bcdataUp, thistime, dt, ID, numberPairs, ManningsN, Length,   &
                zBottom, xDistance, Breadth, widthDepthData, cellType, realLoop)
                
            if (  count(elem2I(:,e2i_solver) == AC) &
                + count(elemMI(:,eMi_solver) == AC)> zeroI ) then               
                !% Artifical compressibility convergence pseudo time marching loop[    
                call pseudo_time_marching &
                    (elem2R, elemMR, faceR, elem2I, elemMI, faceI, elem2YN, elemMYN,    &
                    faceYN, bcdataDn, bcdataUp, linkI, thisStep, thisTime, AnormH,      &
                    AnormQ, AnormHlast, AnormQlast, TnormH, TnormQ, RTnormH, RTnormQ,   &
                    RLnormH,  RLnormQ, debugfile, diagnostic, threadedfile, ID,         &
                    numberPairs, ManningsN, Length, zBottom, xDistance, Breadth,        &
                    widthDepthData, cellType)
            endif

            !% compute the element froude number (diagnostic only)
            call diagnostic_froude_number (elem2R, elem2I, elemMR, elemMI)

            !% compute the volume conservation
            call diagnostic_volume_conservation &
                (diagnostic, elem2R, elem2I, elemMR, elemMI, faceR, bcdataUp, bcdataDn, restartStep, 2)
                

            if (setting%Debugout%DisplayInterval > 0) then
                if (mod(thisstep,setting%Debugout%DisplayInterval) == 0) then
                    print *, 'Volume Conservation (this step and total) = ', &
                        diagnostic(restartStep)%Volume%ConservationThisStep,  &
                        diagnostic(restartStep)%Volume%ConservationTotal
                endif
            endif

            call output_translation_from_elements_to_link_node &
                (elem2I, elem2R, elem2YN, elemMI, elemMR, elemMYN, faceI, faceR, &
                linkI, linkR, nodeI, nodeR, bcdataUp, bcdataDn, thisstep)
                
            call debug_output &
                (debugfile, nodeR, linkR, elem2R, elem2I, elem2YN, elemMR,  &
                elemMI, elemMYN, faceR, faceI, faceYN,bcdataUp, bcdataDn, thisstep)

            call  output_all_threaded_data_by_link &
                (threadedfile, elem2R, elem2I, elemMR, elemMI, faceR, faceI, linkI, &
                bcdataUp, bcdataDn, thisstep)

            !%  change solver based on A/Afull if both SVE and AC solver is used
            if (setting%Solver%SolverSelect == 'SVE-AC') then
                !%  assign the solver for the next step depending on area
                call assign_solver &
                    (elem2I, elem2R, e2r_Area, e2r_FullArea, e2i_elem_type, ePipe,  &
                    e2i_solver, e2r_Temp, e2r_n_temp, next_e2r_temparray)

                ! !%  HACK: AC solver for Junction Pipe has not derived yet
                call assign_solver &
                    (elemMI, elemMR, eMr_Area, eMr_FullArea, eMi_elem_type, eJunctionPipe, &
                    eMi_solver, eMr_Temp, eMr_n_temp, next_eMr_temparray)
            endif

            !% TEST ROUTINES
            !    call explicit_euler_advance &
            !        (elem2R, elem2I, elem2YN, elemMR, elemMI, elemMYN, faceR, faceI, faceYN, &
            !         bcdataDn, bcdataUp, thistime, dt)

            !    call explicit_test_advance &
            !        (elem2R, elem2I, elem2YN, elemMR, elemMI, elemMYN, faceR, faceI, faceYN, &
            !         bcdataDn, bcdataUp, thistime, dt)
            !
            !% increment the counters
            thisstep    = thisstep + 1
            restartStep = restartStep + 1
            thistime    = nexttime
            nexttime    = thistime + dt
        enddo

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine time_marching
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine save_previous_values &
        (elem2R, elemMR, faceR)
        !
        character(64) :: subroutine_name = 'save_previous_values'

        real(8),   intent(inout)  :: elem2R(:,:), elemMR(:,:), faceR(:,:)
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !%  push the old values down the stack
        !%  N values is the present, N0 is the last time step, and N1
        !%  is the timestep before (needed only for backwards 3rd)
        elem2R(:,e2r_Volume_N1)    = elem2R(:,e2r_Volume_N0)
        elem2R(:,e2r_Volume_N0)    = elem2R(:,e2r_Volume)
        elem2R(:,e2r_Area_N1)      = elem2R(:,e2r_Area_N0)
        elem2R(:,e2r_Area_N0)      = elem2R(:,e2r_Area)
        elem2R(:,e2r_Flowrate_N1)  = elem2R(:,e2r_Flowrate_N0)
        elem2R(:,e2r_Flowrate_N0)  = elem2R(:,e2r_Flowrate)
        elem2R(:,e2r_Eta_N0)       = elem2R(:,e2r_Eta)

        elemMR(:,eMr_Volume_N1)    = elemMR(:,eMr_Volume_N0)
        elemMR(:,eMr_Volume_N0)    = elemMR(:,eMr_Volume)
        elemMR(:,eMr_Area_N1)      = elemMR(:,eMr_Area_N0)
        elemMR(:,eMr_Area_N0)      = elemMR(:,eMr_Area)
        elemMR(:,eMr_Flowrate_N1)  = elemMR(:,eMr_Flowrate_N0)
        elemMR(:,eMr_Flowrate_N0)  = elemMR(:,eMr_Flowrate)
        elemMR(:,eMr_Eta_N0)       = elemMR(:,eMr_Eta)

        faceR(:,fr_Flowrate_N0)    = faceR(:,fr_Flowrate)
        
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine save_previous_values
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine assign_solver &
        (elemI, elemR, er_Area, er_FullArea, ei_elem_type, ThisElemType,  &
        ei_solver, er_Temp, er_n_temp, next_er_temparray)
        !
        character(64) :: subroutine_name = 'assign_solver'

        real(8),      target, intent(inout)  :: elemR(:,:)
        integer,   target, intent(inout)  :: elemI(:,:)

        integer,   intent(in)      ::  er_Area, er_FullArea, er_n_temp
        integer,   intent(in)      ::  ei_elem_type, ei_solver, ThisElemType
        integer,   intent(in)      ::  er_Temp(:)
        integer,   intent(inout)   ::  next_er_temparray

        real(8),      pointer  :: AoverAfull(:), Area(:), FullArea(:)
        real(8)                :: switchBufferP, switchBufferM
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


        !%  temporary space for pipe elements
        AoverAfull => elemR(:,er_Temp(next_er_temparray))
        next_er_temparray = utility_advance_temp_array (next_er_temparray,er_n_temp)

        Area     => elemR(:,er_Area)
        FullArea => elemR(:,er_FullArea)

        !%  values +/- buffer for the solver switch
        switchBufferP = setting%DefaultAC%Switch%Area + setting%DefaultAC%Switch%Buffer
        switchBufferM = setting%DefaultAC%Switch%Area - setting%DefaultAC%Switch%Buffer

        where ( elemI(:,ei_elem_type) == ThisElemType )
            AoverAfull = Area / FullArea
        endwhere

        ! selecting appropriate solver for pipe
        where ( (elemI(:,ei_elem_type) == ThisElemType) .and. &
                (elemI(:,ei_solver) == SVE)             .and. &
                (AoverAfull .GE. switchBufferP) )

            elemI(:,ei_solver) = AC

        elsewhere( (elemI(:,ei_elem_type) == ThisElemType) .and. &
                   (elemI(:,ei_solver) == AC)              .and. &
                   (AoverAfull .LE. switchBufferM) )

            elemI(:,ei_solver) = SVE
        endwhere

        ! print*, trim(subroutine_name)
        ! print*, 'AoverAfull', AoverAfull
        ! print*, 'Selected solver', elemI(:,ei_solver)
        
        AoverAfull = nullvalueR
        nullify(AoverAfull)
        next_er_temparray = next_er_temparray - 1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine assign_solver
    !
    !==========================================================================
    !==========================================================================
    !
    ! subroutine dynamic_time_step_size &
    !     (elem2R, elemMR, faceR)
    !     ! 
    !     ! Increases the model time step if the CFL everywhere is low, and 
    !     ! decreases the model time step if the CFL everywhere is high.
    !     ! To prevent this from oscillating between increase and decrease on
    !     ! successive time steps, there is a stepcounter that keeps track of the
    !     ! number of time steps since the last decrease. There must be a number
    !     ! of steps (cfl_increase_stepinterval_from_decrease) without a decrease
    !     ! before the step size can be increased.
    !     ! 
    !     character(64) :: subroutine_name = 'dynamic_time_step_size'

    !     real(8),   intent(inout)  :: elem2R(:,:), elemMR(:,:), faceR(:,:)
    !     real(8)                   :: maxcfl
    !     !--------------------------------------------------------------------------
    !     if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


    !     maxcfl = diagnostic_CFL(elem2R, e2r_Timescale_Q_u, e2r_Timescale_Q_d)
        
    !     if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    ! end subroutine dynamic_time_step_size
    !
    !==========================================================================
    ! END OF MODULE time_loop
    !==========================================================================
end module time_loop

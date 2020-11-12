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
        (elem2R, elemMR, faceR, elem2I, elemMI, faceI, elem2YN, elemMYN,   &
        faceYN, bcdataDn, bcdataUp, linkI, debugfile, diagnostic,         &
        threadedfile, ID, numberPairs, ManningsN, Length, zBottom,        &
        xDistance, Breadth, widthDepthData, cellType)
        !
        ! top-level iteration for continuity and momentum solution
        !
        character(64) :: subroutine_name = 'time_marching'

        real,      target, intent(in out) :: elem2R(:,:),  elemMR(:,:),  faceR(:,:)
        integer,   target, intent(in out) :: elem2I(:,:),  elemMI(:,:),  faceI(:,:)
        logical,   target, intent(in out) :: elem2YN(:,:), elemMYN(:,:), faceYN(:,:)

        type(bcType),                  intent(in out) :: bcdataDn(:), bcdataUp(:)
        type(debugfileType),  target,  intent(in)     :: debugfile(:)
        type(diagnosticType), target,  intent(in out) :: diagnostic(:)
        type(threadedfileType),        intent(in)     :: threadedfile(:)

        integer,                       intent(in)     :: linkI(:,:)

        real, pointer :: rkVol(:), rkU(:)
        real, pointer :: fQ(:), fUdn(:), fUup(:), fAdn(:), fAup(:)
        real, pointer :: fEdn(:), fEup(:), eE(:)

        real    :: AnormH(3), AnormQ(3), AnormHlast(3), AnormQlast(3)
        real    :: TnormH(3), TnormQ(3), RTnormH(3), RTnormQ(3), RLnormH(3), RLnormQ(3)

        real, pointer :: thistime, nexttime

        integer, pointer :: thisStep, restartStep

        logical :: rkCycle(2)
        integer :: ii

        character(len=32) :: outdataName

        integer, intent(in out)    :: ID(:)
        integer, intent(in out)    :: numberPairs(:)
        real,    intent(in out)    :: ManningsN(:)
        real,    intent(in out)    :: Length(:)
        real,    intent(in out)    :: zBottom(:)
        real,    intent(in out)    :: xDistance(:)
        real,    intent(in out)    :: Breadth(:)
        real,    intent(in out)    :: widthDepthData(:,:,:)
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
                    print *, thisstep,'=current step; ', &
                        diagnostic_CFL(elem2R, e2r_Timescale_Q_u, e2r_Timescale_Q_d),'=CFL max' & !Im putting the timescale for Q here. Might needed to be changed later
                        ,maxval(abs(elem2R(2:size(elem2R,1)-1,e2r_Velocity))),'=max velocity'
                    !print *, thisstep,'=current step; ', &
                    !         maxval(elem2R(2:size(elem2R,1)-1,e2r_Flowrate))

                endif
            endif

            call debug_output &
                (debugfile, elem2R, elem2I, elem2YN, elemMR,  &
                elemMI, elemMYN, faceR, faceI, faceYN,bcdataUp, bcdataDn, thisstep)

            !%  Reset the flowrate adhoc detection before flowrates are updated.
            !%  Note that we do not reset the small volume detection here - that is done
            !%  in the element_geometry module before the volumes are updated.
            elem2YN(:,e2YN_IsAdhocFlowrate) = .false.
            elemMYN(:,eMYN_IsAdhocFlowrate) = .false.

            !%  push the old values down the stack
            call save_previous_values (elem2R, elemMR, faceR)

            !% Runge-Kutta 2nd-order advance
            !% rkCycle mask ensures both the RK steps are taken for time loop
            rkCycle(1) = .true.
            rkCycle(2) = .true.

            ! call rk2 &
            !     (elem2R, elemMR, elem2I, elemMI, faceR, faceI, elem2YN, elemMYN, faceYN, &
            !     bcdataDn, bcdataUp, thistime, dt, ID, numberPairs, ManningsN, Length,   &
            !     zBottom, xDistance, Breadth, widthDepthData, cellType, rkCycle)

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
                    ! print*, 'Flowrate ==>'
                    ! call utility_print_values_by_link &
                    !     (elem2R, elem2I, elemMR, elemMI, faceR, faceI, 1, &
                    !     fr_Flowrate, 0, e2r_Flowrate, eMr_Flowrate, eMr_AreaDn, eMr_AreaDn)
                    ! print*, 'Eta ==>'
                    ! call utility_print_values_by_link &
                    !     (elem2R, elem2I, elemMR, elemMI, faceR, faceI, 1, &
                    !     fr_Eta_d, fr_Eta_u, e2r_Eta, eMr_Eta, eMr_EtaDn, eMr_EtaUp)
                    ! print*, 'Area ==>'
                    ! call utility_print_values_by_link &
                    !     (elem2R, elem2I, elemMR, elemMI, faceR, faceI, 1, &
                    !     fr_Area_d, fr_Area_u, e2r_Area, eMr_Area, eMr_AreaDn, eMr_AreaUp)
                    ! print*, 'Depth ==>'
                    ! call utility_print_values_by_link &
                    !     (elem2R, elem2I, elemMR, elemMI, faceR, faceI, 1, &
                    !     fr_HydDepth_d, fr_HydDepth_u, e2r_Depth, eMr_Depth, eMr_AreaDn, eMr_AreaUp)
                    ! print*, 'TopWidth ==>'
                    ! call utility_print_values_by_link &
                    !     (elem2R, elem2I, elemMR, elemMI, faceR, faceI, 1, &
                    !     fr_Topwidth, 0, e2r_Topwidth, eMr_Topwidth, eMr_AreaDn, eMr_AreaUp)
                    ! print*, 'Solver ==>', elem2I(:,e2i_solver)
                endif
            endif

            call  output_all_threaded_data_by_link &
                (threadedfile, elem2R, elem2I, elemMR, elemMI, faceR, faceI, linkI, &
                bcdataUp, bcdataDn, thisstep)

            !% TEST ROUTINES
            !    call explicit_euler_advance &
            !        (elem2R, elem2I, elem2YN, elemMR, elemMI, elemMYN, faceR, faceI, faceYN, &
            !         bcdataDn, bcdataUp, thistime, dt)

            !    call explicit_test_advance &
            !        (elem2R, elem2I, elem2YN, elemMR, elemMI, elemMYN, faceR, faceI, faceYN, &
            !         bcdataDn, bcdataUp, thistime, dt)
            !
            ! increment the counters
            ! print*, 'after pseudo loop'
            ! print*, 'Flowrate ', elem2R(997:1001,e2r_Flowrate)
            ! print*, 'Velocity ', elem2R(997:1001,e2r_Velocity)
            ! print*, 'Eta      ', elem2R(997:1001,e2r_Eta)
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

        real,   intent(inout)  :: elem2R(:,:), elemMR(:,:), faceR(:,:)
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !%  push the old values down the stack
        !%  here elN is the ell (length scale) used in AC derivation
        !%  N values is the present, N0 is the last time step, and N1
        !%  is the timestep before (needed only for backwards 3rd)

        elem2R(:,e2r_elN)          = elem2R(:,e2r_HydDepth)
        elem2R(:,e2r_Area_N1)      = elem2R(:,e2r_Area_N0)
        elem2R(:,e2r_Area_N0)      = elem2R(:,e2r_Area)
        elem2R(:,e2r_Flowrate_N1)  = elem2R(:,e2r_Flowrate_N0)
        elem2R(:,e2r_Flowrate_N0)  = elem2R(:,e2r_Flowrate)
        elem2R(:,e2r_Eta_N0)       = elem2R(:,e2r_Eta)

        elemMR(:,eMr_elN)          = elemMR(:,eMr_HydDepth)
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
    ! END OF MODULE time_loop
    !==========================================================================
end module time_loop

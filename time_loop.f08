! module time_loop
!
! top-level iteration for continuity and momentum solutions
!
!==========================================================================
!
 module time_loop

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

    implicit none

    private

    public :: time_marching

    integer :: debuglevel = 0

    integer,   parameter :: idummy = 0

 contains
!
!==========================================================================
!==========================================================================
!
 subroutine time_marching &
    (elem2R, elemMR, faceR, elem2I, elemMI, faceI, elem2YN, elemMYN, faceYN, &
     bcdataDn, bcdataUp, linkI, debugfile, diagnostic, threadedfile)
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

 real, pointer :: thistime, nexttime

 integer, pointer :: thisStep, restartStep

 integer :: ii

 character(len=32) :: outdataName


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
                     diagnostic_CFL(elem2R, e2r_Timescale_u, e2r_Timescale_d),'=CFL max' &
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

    !% Runge-Kutta 2nd-order advance
    call rk2 &
        (elem2R, elemMR, elem2I, elemMI, faceR, faceI, elem2YN, elemMYN, faceYN, &
         bcdataDn, bcdataUp, thistime, dt)

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
    thisstep    = thisstep + 1
    restartStep = restartStep + 1
    thistime    = nexttime
    nexttime    = thistime + dt
 enddo

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine time_marching
!
!==========================================================================
! END OF MODULE time_loop
!==========================================================================
 end module time_loop
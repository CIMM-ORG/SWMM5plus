module update

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry
    use adjust
    use storage_geometry, only: storage_plan_area_from_volume
    use utility_profiler
    use utility_crash
    !use utility, only: util_syncwrite 
    ! use utility_unit_testing, only: util_utest_CLprint

    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% Updates values during timeloop of hydraulics.
    !%

    private

    public :: update_auxiliary_variables
    public :: update_auxiliary_variables_CC
    public :: update_auxiliary_variables_JMJB
    public :: update_interpweights_JB
    public :: update_element_psi_CC
    !public :: update_element_energyHead_CC
    !public :: update_Froude_number_junction_branch

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine update_auxiliary_variables_CC ( &
        pCol, pCol_Open, pCol_Closed, isAllYN, isSingularYN, idx)
        !%------------------------------------------------------------------
        !% Description:
        !% Updates the variables dependent on the solution of Volume and
        !% Velocity in the RK step
        !% if isSingular is true then only a single point is being updated
        !% if isAllYN is true then all the CC points are being updated
        !% if both are false then some subset is being used
        !% The isALLYN is used for efficient processing of geometry
        !%------------------------------------------------------------------
        !% Declarations
            !% packed columns
            integer,           intent(in) :: pCol, pCol_Open, pCol_Closed 
            !% singular element
            integer,           intent(in) :: idx
            logical,           intent(in) :: isSingularYN, isALLYN

            integer, target, dimension(1) :: zeroIdx, aIdx

            integer, pointer :: thisP(:), thisP_Closed(:), thisP_Open(:)
            !integer, pointer :: thisCol, thisCol_Open, thisCol_Closed
            integer :: npackP, npackP_Closed, npackP_Open
            character(64) :: subroutine_name = 'update_auxiliary_variables_CC'
        !%------------------------------------------------------------------
        !% Aliases
            if (isSingularYN) then 
                !% --- only a single element with index idx is handled
                npackP  = oneI
                aIdx(1) = idx  !% convert the scalar idx to an array
                thisP => Aidx
                if (elemYN(idx,eYN_canSurcharge)) then 
                    thisP_Closed  => aIdx
                    thisP_Open    => zeroIdx
                    npackP_Closed =  oneI
                    npackP_Open   =  zeroI
                else
                    thisP_Closed  => zeroIdx
                    thisP_Open    => aIdx
                    npackP_Open   =  oneI 
                    npackP_Closed =  zeroI
                end if
            else
                !% --- handling a packed set
                npackP        = npack_elemP(pCol)
                npackP_Open   = npack_elemP(pCol_Open)
                npackP_Closed = npack_elemP(pCol_Closed) 

                if (npackP < 1) return 
                thisP => elemP(1:npackP,pCol)
                
                if (npackP_Open > 0) then 
                    thisP_Open => elemP(1:npackP_Open,pCol_Open)
                else 
                    thisP_Open => zeroIdx
                end if

                if (npackP_Closed > 0) then 
                    thisP_Closed => elemP(1:npackP_Closed,pCol_Closed)
                else
                    thisP_Closed => zeroIdx
                end if
            end if
        !%------------------------------------------------------------------

        !% --- update the head (non-surcharged) and geometry
        call geometry_toplevel_CC ( &
            thisP, npackP, thisP_Open, npackP_Open, thisP_Closed, npackP_Closed, &
             isSingularYN, isAllYN)
        
        if (npackP > 0) then
            !% --- Compute the flowrate on CC.
            call update_flowrate_CC (thisP)

            !% --- compute element Froude numbers for CC
            call update_Froude_number_element (thisP)

            !% --- compute the wave speeds
            call update_wavespeed_element(thisP)

            !% --- compute element-face interpolation weights on CC
            call update_interpweights_CC(thisP)

        end if    

    end subroutine update_auxiliary_variables_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine update_element_psi_CC (thisCol) 
        !%------------------------------------------------------------------
        !% Description:
        !% Updates the 2 * beta * psi * L for conduit/channel elements adjacent
        !% to junction. 
        !% thisCol is packed ep_CC_DownstreamOfJunction or ep_CC_UpstreamOfJunction
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisCol
            integer, pointer :: Npack, thisP(:), fup(:), fdn(:)
            real(8), pointer :: hL(:), Hup(:), Hdn(:), Vup(:), Vdn(:), grav
            real(8), pointer :: psiL(:), beta(:), Qelem(:)
            real(8) :: Qepsilon
            integer :: ii
        !%------------------------------------------------------------------
        !% Preliminaries
            Npack => npack_elemP(thisCol)
            if (Npack < 1) return 
            Qepsilon = setting%ZeroValue%Velocity * setting%ZeroValue%Area
        !%------------------------------------------------------------------
        !% Aliases
            thisP => elemP(1:Npack,thisCol)
            fup   => elemI(:,ei_Mface_uL)
            fdn   => elemI(:,ei_Mface_dl)

            hL    => elemR(:,er_Temp01)  !% head loss
            beta  => elemR(:,er_Temp03)  !% +-1 depending on branch type
            Qelem => elemR(:,er_Flowrate)
            !psiL  => elemR(:,er_2B_psiL)  !% 2 * beta * psi * L term

            Hdn   => faceR(:,fr_Head_d) !% Head downstream face
            Hup   => faceR(:,fr_Head_u) !% Head upstream face
            Vdn   => faceR(:,fr_Velocity_d) !% Velocity downstream face
            Vup   => faceR(:,fr_Velocity_u) !% Velocity upstream face

            grav => setting%Constant%gravity
        !%------------------------------------------------------------------

        !% --- ensure Temp arrays are initialized
        hL(thisP)   = zeroR
        beta(thisP) = zeroR

        !% --- set beta for an upstream (+1) or downstream (-1) branch
        select case (thisCol)
        case (ep_CC_DownstreamOfJunction)
            beta(thisP) = -oneR
        case (ep_CC_UpstreamOfJunction)
            beta(thisP) = +oneR
        case default 
            print *, 'CODE ERROR'
            print *, 'Unexpected case default'
            call util_crashpoint(729873)
        end select

        !% ---compute the head loss over the element
        hL(thisP) = Hdn(fup(thisP)) - Hup(fdn(thisP)) &
            + onehalfR * ( (Vdn(fup(thisP))**2) - (Vup(fdn(thisP))**2) ) / grav

        ! print *, ' '
        ! print *, 'headloss ',hL(2:3)

        !% --- error checking
        !%     head loss, which is + for nominal downstream flow,
        !%     should also be the same sign for Q, so their product should be 
        !%     greater than zero. In some dynamic cases this may not be true.
        !%     When this happens we set hL = 0 so that psiL = 0
        if (any((Qelem(thisP) * hL(thisP)) < zeroR)) then 
            where ((Qelem(thisP) * hL(thisP)) < zeroR)
                hL(thisP) = zeroR
            endwhere
            ! print *, 'CODE ERROR'
            ! print *, 'Unexpected mismatch between flowrate direction and head loss'
            ! print *, 'Likely bug in code'
            ! call util_crashpoint(6298723)
        end if

        ! !% --- compute the 2 \beta * psi * L term
        ! where (Qelem(thisP) > Qepsilon)
        !     psiL(thisP) = twoR * beta(thisP) * hL(thisP) / (Qelem(thisP)**2)
        ! elsewhere
        !     psiL(thisP) = zeroR 
        ! endwhere

        ! print *, ' '
        ! print *, 'thisP ',thisP
        ! print *, 'beta   ',beta(2:3)
        ! print *, 'psiL   ',psiL(2:3)
        ! print *, 'Qelem2 ',Qelem(2:3)**2

    end subroutine update_element_psi_CC    
!%
!%==========================================================================
!%==========================================================================
!% 
    ! subroutine update_element_energyHead_CC (thisCol)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Updates the energy head on elements
    !     !% Must be done after all zerodepth and smalldepth adjustments
    !     !% thisCol should be (e.g.) ep_CC_ETM
    !     !%------------------------------------------------------------------
    !     !% Declarations
    !         integer, intent(in) :: thisCol
    !         integer, pointer    :: thisP(:), Npack
    !         real(8), pointer    :: Ehead(:), Head(:), Velocity(:), grav
    !     !%------------------------------------------------------------------
    !     !% Aliases
    !         Npack => npack_elemP(thisCol)
    !         if (Npack < 1) return
    !         thisP => elemP(1:Npack,thisCol)
    !         grav  => setting%Constant%gravity

    !         Ehead    => elemR(:,er_EnergyHead)
    !         Head     => elemR(:,er_Head)
    !         Velocity => elemR(:,er_Velocity)
    !     !%------------------------------------------------------------------

    !     Ehead(thisP) = Head(thisP) + onehalfR * (Velocity(thisP)**2) / grav
        

    ! end subroutine update_element_energyHead_CC
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine update_auxiliary_variables_JMJB  (forceJByn)
        !%------------------------------------------------------------------
        !% Description:
        !% Updates the variables for JM junctions after their solution
        !%------------------------------------------------------------------
        !% Declarations
            logical, intent(in) :: forceJByn !% .true. forces interpweightQ to favor JB
            integer, pointer    :: Npack, thisCol, thisCol_JB, thisP(:), thisCol_JM
            !integer, pointer    :: elemPGx(:,:), npack_elemPGx(:), col_elemPGx(:)
            integer             :: mm
            character(64) :: subroutine_name = 'update_auxiliary_variables_JMJB'
        !%------------------------------------------------------------------
        !% Aliases
            !% --- set packed column for updated elements
           ! elemPGx                => elemPGetm(:,:)
           ! npack_elemPGx          => npack_elemPGetm(:)
           ! col_elemPGx            => col_elemPGetm(:)
            thisCol_JM             => col_elemP(ep_JM)
            thisCol_JB             => col_elemP(ep_JB)
          
        !%------------------------------------------------------------------

        !% --- geometry for both JM and JB
        call geometry_toplevel_JMJB ()  
       
            ! ! ! call util_utest_CLprint ('------- in update after geometry_toplevel_JMJB')
        
        
        !% --- Froude number and wavespeed on JB
        Npack => npack_elemP(thisCol_JB)
        if (Npack > 0) then 
            thisP => elemP(1:Npack, thisCol_JB)
            call update_Froude_number_element (thisP) 
            call update_wavespeed_element(thisP)
        end if

        !% --- wave speed on JM (no Fr since no velocity on JM)
        Npack => npack_elemP(thisCol_JM)
        if (Npack > 0) then
            thisP => elemP(1:Npack, thisCol_JM)
            call update_wavespeed_element(thisP)
        end if

        !% --- interpolation weights on JB
        Npack => npack_elemp(thisCol_JB) 
        if (Npack > 0) then  
            thisP => elemP(1:Npack,thisCol_JB)
            call update_interpweights_JB (thisP, Npack, forceJByn)
        end if

        
    end subroutine update_auxiliary_variables_JMJB
! !%
! !%==========================================================================
!%==========================================================================
!%    
    subroutine update_auxiliary_variables (whichTM)
        ! !%------------------------------------------------------------------
        ! !% Description:
        ! !% Updates the variables dependent on the TM solution of volume
        ! !% and velocity
        ! !%------------------------------------------------------------------
        ! !% Declarations
             integer, intent(in) :: whichTM  !% indicates which Time marching sets (ALLtm, AC, ETM)
        !     integer, pointer :: thisCol_CC, thisCol_JM
        !     character(64) :: subroutine_name = 'update_auxiliary_variables_OLD'

             print *, 'OBSOLETE'
             stop 6908743
        ! !%------------------------------------------------------------------
        ! !% Preliminaries:
        !     !if (crashYN) return
        !     if (setting%Debug%File%update) &
        !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !     if (setting%Profile%useYN) call util_profiler_start (pfc_update_auxiliary_variables)    
        ! !%------------------------------------------------------------------
        ! !%

        !     print *, 'OBSOLETE' 
        !     stop 298734
        !      ! ! ! ! call util_utest_CLprint ('in update before geometry toplevel')
        
        ! ! !% --- update the head (non-surcharged) and geometry
        ! ! call geometry_toplevel (whichTM)

        ! !      ! ! ! ! call util_utest_CLprint ('in update before adjust_limit_velocity_max')

        ! !      !stop 1098734

        ! ! !% --- adjust velocity with limiters
        ! ! call adjust_limit_velocity_max_CC (whichTM)
        ! ! call util_crashstop(21987)

        ! !     ! ! ! ! call util_utest_CLprint ('in update before update_CC_element_flowrate')

        ! ! !% --- set packed column for updated elements
        ! ! select case (whichTM)
        ! !     case (ALLtm)
        ! !         thisCol_CC  => col_elemP(ep_CC_ALLtm)
        ! !         thisCol_JM  => col_elemP(ep_JM_ALLtm)
        ! !     case (ETM)
        ! !         thisCol_CC  => col_elemP(ep_CC_ETM)
        ! !         thisCol_JM  => col_elemP(ep_JM_ETM)
        ! !     case (AC)
        ! !         thisCol_CC  => col_elemP(ep_CC_AC)
        ! !         thisCol_JM  => col_elemP(ep_JM_AC)
        ! !     case default
        ! !         print *, 'CODE ERROR: time march type unknown for # ', whichTM
        ! !         print *, 'which has key ',trim(reverseKey(whichTM))
        ! !         call util_crashpoint(45834)
        ! ! end select

        ! ! !% --- Compute the flowrate on CC.
        ! ! !%     Note that JM should have 0 flowrate and JB has lagged flowrate at this point.
        ! ! !%     The JB flowrate is not updated until after face interpolation
        ! ! call update_element_flowrate (thisCol_CC)

        ! !     ! ! ! ! call util_utest_CLprint ('in update before update_Froude_number_element')

        ! ! !% --- compute element Froude numbers for CC
        ! ! call update_Froude_number_element (thisCol_CC)

        ! !      ! ! ! ! call util_utest_CLprint ('in update before CC wavespeed')

        ! ! !% --- compute the wave speeds
        ! ! call update_wavespeed_element(thisCol_CC)

        ! !     ! ! ! ! call util_utest_CLprint ('in update before JM wavespeed')

        ! ! call update_wavespeed_element(thisCol_JM)

        ! !     ! ! ! ! call util_utest_CLprint ('in update before CC interpweights')

        ! ! !% --- compute element-face interpolation weights on CC
        ! ! call update_interpweights_CC(thisCol_CC, whichTM)

        ! !     ! ! ! ! call util_utest_CLprint ('in update before JB interpweights')

        ! ! !% --- compute element-face interpolation weights on JB
        ! ! !%     .false. as the general call to update aux does not force JB
        ! ! call update_interpweights_JB (thisCol_JM, .false.)

        
        ! !     ! ! ! ! call util_utest_CLprint ('in update before update Froude Number Junction Branch')

        ! ! !% --- compute element Froude number for JB
        ! ! call update_Froude_number_JB (thisCol_JM) 

        ! !     ! ! ! ! call util_utest_CLprint ('in update before update BCoutlet_flowrate')

        ! ! !% --- not needed 20220716brh
        ! ! !% --- flow values on an BC outlet face 20220714brh
        ! ! !%     required so that an inflow to a zero or small depth will not be lost
        ! ! ! call update_BCoutlet_flowrate ()

        ! ! !%------------------------------------------------------------------
        ! ! !% Closing:
        ! !     if (setting%Profile%useYN) call util_profiler_stop (pfc_update_auxiliary_variables)

        ! !      if (setting%Debug%File%update)  &
        ! !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]" 
    end subroutine update_auxiliary_variables
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================

    subroutine update_flowrate_CC (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !%
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisP(:)
            real(8), pointer :: flowrate(:), velocity(:), area(:), Qmax(:)
            character(64) :: subroutine_name = 'update_element_flowrate'
        !%------------------------------------------------------------------
        !% Aliases
            flowrate => elemR(:,er_Flowrate)
            velocity => elemR(:,er_Velocity)
            area     => elemR(:,er_Area)
            Qmax     => elemR(:,er_FlowrateLimit)
        !%-----------------------------------------------------------------------------
 
        flowrate(thisP) = area(thisP) * velocity(thisP)

        ! ! ! call util_utest_CLprint ('in update element Flowrate B')

        !% --- limit flowrate by the full value (if it exists)
        where ((Qmax(thisP) > zeroR) .and. (abs(flowrate(thisP)) > Qmax(thisP)))
            flowrate(thisP) = sign(Qmax(thisP), flowrate(thisP))
        end where
        
        ! ! ! call util_utest_CLprint ('in update element Flowrate C')

        ! print *, flowrate(139), area(139), velocity(139)
        ! print*, flowrate(thisP), 'flowrate(thisP)'

    end subroutine update_flowrate_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine update_Froude_number_element (thisP)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes Froude number on each element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisP(:)
        real(8), pointer :: Froude(:), velocity(:), ellDepth(:), grav
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        Froude   => elemR(:,er_FroudeNumber)
        velocity => elemR(:,er_Velocity)
        ellDepth    => elemR(:,er_EllDepth)  !% Use the ell value (modified hydraulic depth)
        grav     => setting%constant%gravity
        !%-----------------------------------------------------------------------------

        Froude(thisP) = velocity(thisP) / sqrt(grav * ellDepth(thisP))
      

    end subroutine update_Froude_number_element
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine update_Froude_number_JB (thisCol_JM)
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% computes Froude number on each junction branch element
    !     !% BRHbugfix 20210812
    !     !%-----------------------------------------------------------------------------
    !     character(64) :: subroutine_name = 'update_Froude_number_JB'
    !     integer, intent(in) :: thisCol_JM
    !     integer, pointer :: Npack, thisP(:), tM, BranchExists(:)
    !     real(8), pointer :: Froude(:), velocity(:), Depth(:), grav
    !     integer :: ii, kk, tB
    !     !%-----------------------------------------------------------------------------
    !     !if (crashYN) return
    !     Froude   => elemR(:,er_FroudeNumber)
    !     velocity => elemR(:,er_Velocity)
    !     !ellDepth    => elemR(:,er_EllDepth)  !% Use the ell value (modified hydraulic depth) NOT AVAILALBE
    !     Depth       => elemR(:,er_Depth)
    !     BranchExists => elemSI(:,esi_JunctionBranch_Exists)
    !     grav     => setting%constant%gravity
    !     !%-----------------------------------------------------------------------------
    !     if (setting%Debug%File%update) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !     Npack => npack_elemP(thisCol_JM)
    !     if (Npack > 0) then
    !         thisP => elemP(1:Npack,thisCol_JM)
    !         do ii=1,Npack
    !             tM => thisP(ii)
    !             do kk=1,max_branch_per_node
    !                 tB = tM + kk
    !                 if (BranchExists(tB)==oneI) then
    !                     !Froude(tB) = velocity(tB) / sqrt(grav * ellDepth(tB))
    !                     Froude(tB) = velocity(tB) / sqrt(grav * Depth(tB))
    !                     !print *, kk, tB, Froude(tB), velocity(tB),'  Froude JB'
    !                 end if
    !             end do
    !         end do
    !     end if

    !     if (setting%Debug%File%update)  &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine update_Froude_number_JB
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine update_wavespeed_element (thisP)
        !%------------------------------------------------------------------
        !% Description
        !% computes the wavespeed on a CC or JB element
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisP(:)
            real(8), pointer :: wavespeed(:), ellDepth(:), grav
        !%------------------------------------------------------------------
        !% Aliases:
            wavespeed => elemR(:,er_WaveSpeed)
            ellDepth  => elemR(:,er_EllDepth)
            grav      => setting%constant%gravity
        !%------------------------------------------------------------------

       ! print *, 'in update wavespeed element'    
        
        !% wavespeed at modified hydraulic depth (ell) 
        wavespeed(thisP) = sqrt(grav * ellDepth(thisP))

    end subroutine update_wavespeed_element    
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine update_interpweights_CC (thisP)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes the interpolation weights on each element for CC
        !% tim-marching elements
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'update_interpweights_CC'
        integer, intent(in) :: thisP(:)
        integer, pointer :: fUp(:), fDn(:)
        real(8), pointer :: velocity(:), wavespeed(:), ellDepth(:), length(:), QLateral(:)
        real(8), pointer :: PCelerity(:), SlotVolume(:),SlotWidth(:), fullArea(:)
        real(8), pointer :: w_uQ(:), w_dQ(:),  w_uG(:), w_dG(:),  w_uH(:), w_dH(:), w_uP(:), w_dP(:), Area(:)
        real(8), pointer :: Fr(:), grav !BRHbugfix20210811 test
        logical, pointer :: isSlot(:), fSlot(:)
        integer :: ii
        !%-----------------------------------------------------------------------------
        ! if (crashYN) return
        if (setting%Debug%File%update) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        Qlateral  => elemR(:,er_FlowrateLateral)
        velocity  => elemR(:,er_Velocity)
        wavespeed => elemR(:,er_WaveSpeed)
        ellDepth  => elemR(:,er_EllDepth)  !% modified hydraulic depth!
        length    => elemR(:,er_Length)
        w_uQ      => elemR(:,er_InterpWeight_uQ)
        w_dQ      => elemR(:,er_InterpWeight_dQ)
        w_uG      => elemR(:,er_InterpWeight_uG)
        w_dG      => elemR(:,er_InterpWeight_dG)
        w_uH      => elemR(:,er_InterpWeight_uH)
        w_dH      => elemR(:,er_InterpWeight_dH)
        w_uP      => elemR(:,er_InterpWeight_uP)
        w_dP      => elemR(:,er_InterpWeight_dP)
        Fr        => elemR(:,er_FroudeNumber)  !BRHbugfix20210811 test
        isSlot    => elemYN(:,eYN_isPSsurcharged)  !% Preissmann

        fSlot    => faceYN(:,fYN_isPSsurcharged)  !% Preissmann
        fUp      => elemI(:,ei_Mface_uL)
        fDn      => elemI(:,ei_Mface_dL)

        PCelerity  => elemR(:,er_Preissmann_Celerity)
        SlotVolume => elemR(:,er_SlotVolume) !% Preissmann
        SlotWidth  => elemR(:,er_SlotWidth)  !% Preissmann
        fullArea   => elemR(:,er_FullArea)
        grav       => setting%constant%gravity


        Area => faceR(:,er_Area)
        !%-----------------------------------------------------------------------------
        !% 2nd cases needed for handling surcharged AC elements and using the celerity
        !% multiplier of the AC method for the wavespeed
        ! select case (whichTM)
        !     case (ALLtm)
        !         !thisCol_AC          =>  col_elemP(ep_ACsurcharged)
        !         print *, 'CODE ERROR: ALLtm not complete'
        !         call util_crashpoint(5598723) 
        !     case (ETM)
        !         thisCol_ClosedElems =>  col_elemP(ep_CC_Closed_Elements)
        !     case (AC)
        !         !thisCol_AC          =>  col_elemP(ep_ACsurcharged)
        !         print *, 'CODE ERROR: AC not complete'
        !         call util_crashpoint(55987233) 
        !     case default
        !         print *, 'CODE ERROR: time march type unknown for # ', whichTM
        !         print *, 'which has key ',trim(reverseKey(whichTM))
        !         stop 3987
        ! end select

        ! Npack => npack_elemP(thisCol)
        ! if (Npack < 1) return

        ! thisP => elemP(1:Npack,thisCol)

        !% wavespeed at modified hydraulic depth (ell)
        wavespeed(thisP) = sqrt(grav * EllDepth(thisP))

        ! print *, ' '
        ! print *, 'in update interpweights CC'
        ! print *, wavespeed(49),velocity(49), Fr(49)
    
        ! !% modify wavespeed for surcharged AC cells
        ! if (whichTM .ne. ETM) then
        !     Npack2 => npack_elemP(thisCol_AC)
        !     if (Npack2 > 0) then
        !         thisP2 => elemP(1:Npack2,thisCol_AC)
        !         wavespeed(thisP2) = wavespeed(thisP2) * setting%ACmethod%Celerity%RC
        !     end if
        ! end if

        where (.not. isSlot(thisP))
            w_uQ(thisP) = - onehalfR * length(thisP)  / (abs(Fr(thisp)**0) * velocity(thisP) - wavespeed(thisP)) !bugfix SAZ 09212021 
            w_dQ(thisP) = + onehalfR * length(thisP)  / (abs(Fr(thisp)**0) * velocity(thisP) + wavespeed(thisP)) !bugfix SAZ 09212021 
        elsewhere (isSlot(thisP))
            !% --- Preissmann slot
            w_uQ(thisP) = - onehalfR * length(thisP)  / (abs(Fr(thisp)**0) * velocity(thisP) - PCelerity(thisP)) !bugfix SAZ 23022022 
            w_dQ(thisP) = + onehalfR * length(thisP)  / (abs(Fr(thisp)**0) * velocity(thisP) + PCelerity(thisP)) !bugfix SAZ 23022022 
        end where

        !% apply limiters to timescales
        where (w_uQ(thisP) < zeroR)
            w_uQ(thisP) = setting%Limiter%InterpWeight%Maximum
        endwhere
        where (w_uQ(thisP) < setting%Limiter%InterpWeight%Minimum)
            w_uQ(thisP) = setting%Limiter%InterpWeight%Minimum
        endwhere
        where (w_uQ(thisP) > setting%Limiter%InterpWeight%Maximum)
            w_uQ(thisP) = setting%Limiter%InterpWeight%Maximum
        endwhere

        where (w_dQ(thisP) < zeroR)
            w_dQ(thisP) = setting%Limiter%InterpWeight%Maximum
        endwhere
        where (w_dQ(thisP) < setting%Limiter%InterpWeight%Minimum)
            w_dQ(thisP) = setting%Limiter%InterpWeight%Minimum
        endwhere
        where (w_dQ(thisP) > setting%Limiter%InterpWeight%Maximum)
            w_dQ(thisP) = setting%Limiter%InterpWeight%Maximum
        endwhere

        !% timescale interpolation for geometry are identical to flowrate
        !% but may be modified elsewhere
        w_uG(thisP) = w_uQ(thisP)
        w_dG(thisP) = w_dQ(thisP)
        w_uP(thisP) = w_uQ(thisP)
        w_dP(thisP) = w_dQ(thisP)

        !% head uses length scale interpolation
        !% This shouldn't need limiters.
        w_uH(thisP) = onehalfR * length(thisP)
        w_dH(thisP) = onehalfR * length(thisP)

        ! print *, ' '
        ! print *, 'in update interpweights CC'
        ! print *, 'G ',elemR(49,er_InterpWeight_uG),elemR(49,er_InterpWeight_dG)
        ! print *, ' '

        !% adjust upstream interpolation weights for downstream flow in presence of lateral inflows
        !% so that upstream interpolation is used
        !% HACK -- this probably could use an approach with some kind of ad hoc blend -- needs work

        !% 20220817brh REMOVING LATERAL RESET AS IT IS CAUSING OSCILLATIONS IN HIGH INSTREAM FLOW CONDITIONS
        !% MAY NEED TO PUT IT BACK IN FOR CASES WHERE LATERAL FLOWRATE IS LARGER THAN DOWNSTREAM FLOW

        ! where ( (velocity(thisP) > zeroR) .and. (Qlateral(thisP) > zeroR) )
        !     w_uQ(thisP) =  setting%Limiter%InterpWeight%Maximum
        !     w_uG(thisP) =  setting%Limiter%InterpWeight%Maximum
        ! endwhere

        ! ! !% adjust downstream interpolation weights for upstream flow in presence of lateral inflow
        ! where ( (velocity(thisP) < zeroR) .and. (Qlateral(thisP) > zeroR) )
        !     w_dQ(thisP) = setting%Limiter%InterpWeight%Maximum
        !     w_dG(thisP) = setting%Limiter%InterpWeight%Maximum
        ! endwhere

        if (setting%Debug%File%update)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine update_interpweights_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine update_interpweights_JB (thisP, Npack, forceJBQyn)
        !%------------------------------------------------------------------
        !% Description:
        !% compute the interpolation weights for junction branches
        !% inpute is the set of all JB
        !%------------------------------------------------------------------
            integer, intent(in) :: thisP(:), Npack
            logical, intent(in) :: forceJBQyn !% --- if true then forces JB weight to max
            integer             :: ii, mm, jB
            real(8), pointer    :: grav, wavespeed(:), PCelerity(:), velocity(:), length(:), depth(:)
            real(8), pointer    :: w_uQ(:), w_dQ(:), w_uG(:), w_dG(:), w_uH(:), w_dH(:), w_uP(:), w_dP(:)
            logical, pointer    :: isSlot(:)
        !%------------------------------------------------------------------
        !% Aliases
            grav => setting%Constant%gravity
            velocity  => elemR(:,er_Velocity)
            wavespeed => elemR(:,er_WaveSpeed)
            PCelerity => elemR(:,er_Preissmann_Celerity)
            depth     => elemR(:,er_EllDepth)  !% modified hydraulic depth!
            length    => elemR(:,er_Length)
            w_uQ      => elemR(:,er_InterpWeight_uQ)
            w_dQ      => elemR(:,er_InterpWeight_dQ)
            w_uG      => elemR(:,er_InterpWeight_uG)
            w_dG      => elemR(:,er_InterpWeight_dG)
            w_uH      => elemR(:,er_InterpWeight_uH)
            w_dH      => elemR(:,er_InterpWeight_dH)
            w_uP      => elemR(:,er_InterpWeight_uP)
            w_dP      => elemR(:,er_InterpWeight_dP)
            isSlot    => elemYN(:,eYN_isPSsurcharged)  !% Preissmann
        !%------------------------------------------------------------------

        !% cycle through the branches to compute weights
        do mm=1,Npack
            !% --- JB index
            jB = thisP(mm)
        
            !% --- cycle if not valid
            if (elemSI(jB,esi_JunctionBranch_Exists) .ne. oneI) cycle

            if (.not. forceJBQyn) then
                !% --- flowrate interpweight as time scaled
                if (.not. isSlot(jB)) then
                    w_uQ(jB) = - onehalfR * length(jB)  / (velocity(jB) - wavespeed(jB))
                    w_dQ(jB) = + onehalfR * length(jB)  / (velocity(jB) + wavespeed(jB))
                else
                    !% --- Preissmann slot
                    w_uQ(jB) = - onehalfR * length(jB)  / (velocity(jB) - PCelerity(jB))
                    w_dQ(jB) = + onehalfR * length(jB)  / (velocity(jB) + PCelerity(jB))
                end if
            else
                w_uQ(jB) = setting%Limiter%InterpWeight%Minimum
                w_dQ(jB) = setting%Limiter%InterpWeight%Minimum
            end if 

            !% --- geometry interpweight as time scaled
            if (.not. isSlot(jB)) then
                w_uG(jB) = - onehalfR * length(jB)  / (velocity(jB) - wavespeed(jB))
                w_dG(jB) = + onehalfR * length(jB)  / (velocity(jB) + wavespeed(jB))
            else
                !% --- Preissmann slot
                w_uG(jB) = - onehalfR * length(jB)  / (velocity(jB) - PCelerity(jB))
                w_dG(jB) = + onehalfR * length(jB)  / (velocity(jB) + PCelerity(jB))
            end if

            !% apply upstream limiters to timescales for geometry
            if (w_uG(jB) < zeroR) then
                w_uG(jB) = setting%Limiter%InterpWeight%Maximum
            end if
            if (w_uG(jB) < setting%Limiter%InterpWeight%Minimum) then
                w_uG(jB) = setting%Limiter%InterpWeight%Minimum
            end if
            if (w_uG(jB) > setting%Limiter%InterpWeight%Maximum) then
                w_uG(jB) = setting%Limiter%InterpWeight%Maximum
            end if

            !% apply downstream limiters to timescales for geometry
            if (w_dG(jB) < zeroR) then
                w_dG(jB) = setting%Limiter%InterpWeight%Maximum
            end if
            if (w_dG(jB) < setting%Limiter%InterpWeight%Minimum) then
                w_dG(jB) = setting%Limiter%InterpWeight%Minimum
            end if
            if (w_dG(jB) > setting%Limiter%InterpWeight%Maximum) then
                w_dG(jB) = setting%Limiter%InterpWeight%Maximum
            end if

            !% appl upstream limiters to timescales for flowrate
            if (w_uQ(jB) < zeroR) then
                w_uQ(jB) = setting%Limiter%InterpWeight%Maximum
            end if
            if (w_uQ(jB) < setting%Limiter%InterpWeight%Minimum) then
                w_uQ(jB) = setting%Limiter%InterpWeight%Minimum
            end if
            if (w_uQ(jB) > setting%Limiter%InterpWeight%Maximum) then
                w_uQ(jB) = setting%Limiter%InterpWeight%Maximum
            end if

            !% apply downstream limiters to timescales for flowrate
            if (w_dQ(jB) < zeroR) then
                w_dQ(jB) = setting%Limiter%InterpWeight%Maximum
            end if
            if (w_dQ(jB) < setting%Limiter%InterpWeight%Minimum) then
                w_dQ(jB) = setting%Limiter%InterpWeight%Minimum
            end if
            if (w_dQ(jB) > setting%Limiter%InterpWeight%Maximum) then
                w_dQ(jB) = setting%Limiter%InterpWeight%Maximum
            end if

            !% OBSOLETE -- the decision on where to set the Preissmann interp
            !% is made in face.f90 in the choice of whether it is in the G, H or Q
            !% interpolation sets
            !% set the Preissman interp the same as geometry interp
            !w_uG(jB) = w_uQ(jB)
            !w_dG(jB) = w_dQ(jB)
            !w_uP(jB) = w_uG(jB)
            !w_dP(jB) = w_dG(jB)

            !% --- set head interp as length-scaled (always)
            w_uH(jB) = onehalfR * length(jB)
            w_dH(jB) = onehalfR * length(jB)

            !print *, 'JB, iw ',JB, w_dH(jB)
        end do

        ! print *, ' '
        ! print *, 'in update interpweights_JB ',w_dH(13)
        ! print *, 'length ',length(13)
        ! print *, ' '
        !stop 5098723

    end subroutine update_interpweights_JB
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine update_BCoutlet_flowrate ()
        !%------------------------------------------------------------------
        !% Description:
        !% sets the outlet (face) flowrate equal to the interior element
        !% flowrate to ensure
        !%------------------------------------------------------------------
        !% Declarations
            integer, pointer :: eup(:), idx_fBC(:)
        !%------------------------------------------------------------------
        !% Aliases
            if (npack_faceP(fp_BCdn) < 1) return
            eup       => faceI(:,fi_Melem_uL)
            idx_fBC   => faceP(1:npack_faceP(fp_BCdn),fp_BCdn)
        !%------------------------------------------------------------------    
        !%    
        faceR(idx_fBC,fr_Flowrate) = elemR(eup(idx_fBC),er_Flowrate)

    end subroutine update_BCoutlet_flowrate    
!%
!%==========================================================================
!%==========================================================================
!%
    ! subroutine update_interpolation_weights_ds_JB ()
    !     !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% This is a test subroutine that violates no neighbour algorithm
    !     !% This subroutine sets the interpolation wights in ds JB to its
    !     !% conneceted link element
    !     !%-----------------------------------------------------------------------------
    !     character(64) :: subroutine_name = 'update_interpolation_weights_ds_JB'
    !     integer, pointer :: thisColP_dsJB, thisColP_ds_of_JB
    !     integer, pointer :: Npack1, Npack2,  thisP1(:), thisP2(:)
    !     real(8), pointer :: w_uQ(:), w_dQ(:),  w_uG(:), w_dG(:),  w_uH(:), w_dH(:)
    !     !%-----------------------------------------------------------------------------
    !     !if (crashYN) return
    !     if (setting%Debug%File%update)  &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    !     w_uQ      => elemR(:,er_InterpWeight_uQ)
    !     w_dQ      => elemR(:,er_InterpWeight_dQ)
    !     w_uG      => elemR(:,er_InterpWeight_uG)
    !     w_dG      => elemR(:,er_InterpWeight_dG)
    !     w_uH      => elemR(:,er_InterpWeight_uH)
    !     w_dH      => elemR(:,er_InterpWeight_dH)
    !     !%-----------------------------------------------------------------------------

    !     !% replace the interpolation weights for downstream JB
    !     thisColP_dsJB  => col_elemP(ep_JB_Downstream)
    !     Npack1         => npack_elemP(thisColP_dsJB)

    !     if (Npack1 > 0) then
    !         thisP1 => elemP(1:Npack1,thisColP_dsJB)
    !         w_dQ(thisP1) = oneR
    !         w_dG(thisP1) = oneR
    !         w_dH(thisP1) = oneR
    !     end if

    !     thisColP_ds_of_JB => col_elemP(ep_CC_DownstreamOfJunction)
    !     Npack2            => npack_elemP(ep_CC_DownstreamOfJunction)

    !     !% replace the interpolation weights for elements downstream of dn JB
    !     if (Npack2 > 0) then
    !         thisP2 => elemP(1:Npack2,thisColP_ds_of_JB)
    !         w_uQ(thisP2) = oneR
    !         w_uG(thisP2) = oneR
    !         w_uH(thisP2) = oneR
    !     end if

    !     if (setting%Debug%File%update) &
    !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end subroutine update_interpolation_weights_ds_JB
    !%
    !%==========================================================================
    !% END OF MODULE
    !%==========================================================================
end module update
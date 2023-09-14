module update
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Updates depth, head, flowrate based on geometry during time-march
    !%
    !% Methods:
    !% Element updates based on new volume and geometry
    !%==========================================================================

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry
    use adjust
    use storage_geometry, only: storage_plan_area_from_volume
    use utility_profiler
    use utility_crash

    use utility_unit_testing, only: util_utest_CLprint

    implicit none

    private

    public :: update_auxiliary_variables_CC
    public :: update_Froude_number_element
    public :: update_wavespeed_element
    public :: update_interpweights_JB
    public :: update_interpweights_Diag
    

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

            if (.not. isSingularYN) call util_utest_CLprint('    aaa update - - - - - - - - - - ')

        !% --- update the head (non-surcharged) and geometry
        call geometry_toplevel_CC ( &
            thisP, npackP, thisP_Open, npackP_Open, thisP_Closed, npackP_Closed, &
             isSingularYN, isAllYN)

             if (.not. isSingularYN)  call util_utest_CLprint('    bbb update - - - - - - - - - - ')
        
        if (npackP > 0) then
            !% --- Compute the flowrate on CC.
            call update_flowrate_CC (thisP)

            ! if (.not. isSingularYN)  call util_utest_CLprint('    ccc update - - - - - - - - - - ')

            !% --- compute element Froude numbers for CC
            call update_Froude_number_element (thisP)

            ! if (.not. isSingularYN) call util_utest_CLprint('    ddd update - - - - - - - - - - ')

            !% --- compute the wave speeds
            call update_wavespeed_element(thisP)

            ! if (.not. isSingularYN)  call util_utest_CLprint('    eee update - - - - - - - - - - ')

            !% --- compute element-face interpolation weights on CC
            call update_interpweights_CC(thisP)

            ! if (.not. isSingularYN) call util_utest_CLprint('    fff update - - - - - - - - - - ')

            !% --- compute element total energyhead 
            call update_energyhead_CC(thisP)

        end if    

    end subroutine update_auxiliary_variables_CC
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
            real(8), pointer    :: Froude(:), velocity(:), ellDepth(:), grav
        !%-----------------------------------------------------------------------------
        !% Aliases
            Froude   => elemR(:,er_FroudeNumber)
            velocity => elemR(:,er_Velocity)
            ellDepth => elemR(:,er_EllDepth)  !% Use the ell value (modified hydraulic depth)
            grav     => setting%constant%gravity
        !%-----------------------------------------------------------------------------

        Froude(thisP) = velocity(thisP) / sqrt(grav * ellDepth(thisP))
      
    end subroutine update_Froude_number_element
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
        
        !% --- wavespeed at modified hydraulic depth (ell) 
        wavespeed(thisP) = sqrt(grav * ellDepth(thisP))

    end subroutine update_wavespeed_element    
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
            integer             :: mm, jB
            real(8), pointer    :: grav, wavespeed(:), PCelerity(:), velocity(:), length(:), depth(:)
            real(8), pointer    :: w_uQ(:), w_dQ(:), w_uG(:), w_dG(:), w_uH(:), w_dH(:)
            logical, pointer    :: isSlot(:)
        !%------------------------------------------------------------------
        !% Aliases
            grav      => setting%Constant%gravity
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
            isSlot    => elemYN(:,eYN_isPSsurcharged)  !% Preissmann
        !%------------------------------------------------------------------

        !% --- cycle through the branches to compute weights
        do mm=1,Npack
            !% --- JB index
            jB = thisP(mm)
        
            !% --- cycle if not valid
            if (elemSI(jB,esi_JunctionBranch_Exists) .ne. oneI) cycle

            !% --- if zero depth, then JB is maximum for all interp
            if (depth(jB) .le. setting%ZeroValue%Depth) then 
                w_uQ(jB) = setting%Limiter%InterpWeight%Maximum
                w_dQ(jB) = setting%Limiter%InterpWeight%Maximum
                w_uG(jB) = setting%Limiter%InterpWeight%Maximum
                w_dG(jB) = setting%Limiter%InterpWeight%Maximum
                w_uH(jB) = setting%Limiter%InterpWeight%Maximum
                w_dH(jB) = setting%Limiter%InterpWeight%Maximum
                cycle
            end if
            
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

            !% --- set head interp as length-scaled (always)
            w_uH(jB) = onehalfR * length(jB)
            w_dH(jB) = onehalfR * length(jB)

        end do

    end subroutine update_interpweights_JB
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine update_interpweights_Diag (thisP, Npack)
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets the interpolation weights for diagnostic elements
        !% only called during initialization
        !%-----------------------------------------------------------------
        !% Declarations
          integer, intent(in) :: thisP(:), Npack
          real(8), pointer    :: w_uQ(:), w_dQ(:), w_uG(:), w_dG(:), w_uH(:), w_dH(:)
        !%-----------------------------------------------------------------
        !% Aliases
          w_uQ      => elemR(:,er_InterpWeight_uQ)
          w_dQ      => elemR(:,er_InterpWeight_dQ)
          w_uG      => elemR(:,er_InterpWeight_uG)
          w_dG      => elemR(:,er_InterpWeight_dG)
          w_uH      => elemR(:,er_InterpWeight_uH)
          w_dH      => elemR(:,er_InterpWeight_dH)
        !%-----------------------------------------------------------------
        if (Npack >0) then
          w_uQ(thisP) = setting%Limiter%Interpweight%Minimum
          w_dQ(thisP) = setting%Limiter%Interpweight%Minimum

          w_uG(thisP) = setting%Limiter%Interpweight%Maximum
          w_dG(thisP) = setting%Limiter%Interpweight%Maximum

          w_uH(thisP) = setting%Limiter%Interpweight%Maximum
          w_dH(thisP) = setting%Limiter%Interpweight%Maximum
        end if

    end subroutine update_interpweights_Diag
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine update_flowrate_CC (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% Updates flowrate from velocity and cross-sectional area
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisP(:)
            real(8), pointer    :: flowrate(:), velocity(:), area(:), Qmax(:)
            character(64) :: subroutine_name = 'update_element_flowrate'
        !%------------------------------------------------------------------
        !% Aliases
            flowrate => elemR(:,er_Flowrate)
            velocity => elemR(:,er_Velocity)
            area     => elemR(:,er_Area)
            Qmax     => elemR(:,er_FlowrateLimit)
        !%------------------------------------------------------------------
 
        flowrate(thisP) = area(thisP) * velocity(thisP)

        !% --- limit flowrate by the full value (if it exists)
        where ((Qmax(thisP) > zeroR) .and. (abs(flowrate(thisP)) > Qmax(thisP)))
            flowrate(thisP) = sign(Qmax(thisP), flowrate(thisP))
        end where
        
        !% --- small values are set to zero
        where ((flowrate(thisP) < 1.0d-10) .and. (flowrate(thisP) > -1.0d-10))
            flowrate(thisP) = zeroR
        endwhere

    end subroutine update_flowrate_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine update_interpweights_CC (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% computes the interpolation weights on each element for CC
        !% tim-marching elements
        !%------------------------------------------------------------------
        !% Declarations
            character(64)       :: subroutine_name = 'update_interpweights_CC'
            integer, intent(in) :: thisP(:)
            integer, pointer    :: fUp(:), fDn(:)
            real(8), pointer    :: velocity(:), wavespeed(:), ellDepth(:), length(:), QLateral(:)
            real(8), pointer    :: PCelerity(:), SlotVolume(:),SlotWidth(:), fullArea(:)
            real(8), pointer    :: w_uQ(:), w_dQ(:),  w_uG(:), w_dG(:),  w_uH(:), w_dH(:), w_uP(:), w_dP(:), Area(:)
            real(8), pointer    :: Fr(:), grav
            logical, pointer    :: isSlot(:), fSlot(:)
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%update) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
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

            Area       => faceR(:,er_Area)
        !%------------------------------------------------------------------
        !% --- wavespeed at modified hydraulic depth (ell)
        wavespeed(thisP) = sqrt(grav * EllDepth(thisP))

        !% --- limiters below zero depth
        where (elemR(thisP,er_Depth) .le. setting%ZeroValue%Depth)
            w_uQ(thisP) = setting%Limiter%InterpWeight%Maximum
            w_dQ(thisP) = setting%Limiter%InterpWeight%Maximum
            w_uG(thisP) = setting%Limiter%InterpWeight%Maximum
            w_dG(thisP) = setting%Limiter%InterpWeight%Maximum
            w_uH(thisP) = setting%Limiter%InterpWeight%Maximum
            w_dQ(thisP) = setting%Limiter%InterpWeight%Maximum
        elsewhere

            ! --- free surface uses wave speed, Preissmann Slot use Preissmann Celerity
            where (.not. isSlot(thisP)) 
                w_uQ(thisP) = - onehalfR * length(thisP)  / (abs(Fr(thisp)**0) * velocity(thisP) - wavespeed(thisP)) !bugfix SAZ 09212021 
                w_dQ(thisP) = + onehalfR * length(thisP)  / (abs(Fr(thisp)**0) * velocity(thisP) + wavespeed(thisP)) !bugfix SAZ 09212021 
            elsewhere (isSlot(thisP))
                !% --- Preissmann slot
                w_uQ(thisP) = - onehalfR * length(thisP)  / (abs(Fr(thisp)**0) * velocity(thisP) - PCelerity(thisP)) !bugfix SAZ 23022022 
                w_dQ(thisP) = + onehalfR * length(thisP)  / (abs(Fr(thisp)**0) * velocity(thisP) + PCelerity(thisP)) !bugfix SAZ 23022022 
            end where

            !% --- apply limiters to timescales
            !% --- negative weight indicates supercritical downstream flow
            where (w_uQ(thisP) < zeroR)
                w_uQ(thisP) = setting%Limiter%InterpWeight%Maximum
            endwhere
            where (w_uQ(thisP) < setting%Limiter%InterpWeight%Minimum)
                w_uQ(thisP) = setting%Limiter%InterpWeight%Minimum
            endwhere
            where (w_uQ(thisP) > setting%Limiter%InterpWeight%Maximum)
                w_uQ(thisP) = setting%Limiter%InterpWeight%Maximum
            endwhere

            !% --- negative weight indicates supercritical
            where (w_dQ(thisP) < zeroR)
                w_dQ(thisP) = setting%Limiter%InterpWeight%Maximum
            endwhere
            where (w_dQ(thisP) < setting%Limiter%InterpWeight%Minimum)
                w_dQ(thisP) = setting%Limiter%InterpWeight%Minimum
            endwhere
            where (w_dQ(thisP) > setting%Limiter%InterpWeight%Maximum)
                w_dQ(thisP) = setting%Limiter%InterpWeight%Maximum
            endwhere

            !% --- timescale interpolation for geometry are identical to flowrate
            !%     but may be modified elsewhere
            w_uG(thisP) = w_uQ(thisP)
            w_dG(thisP) = w_dQ(thisP)
            w_uP(thisP) = w_uQ(thisP)
            w_dP(thisP) = w_dQ(thisP)

            !% --- head uses length scale interpolation
            !%     This shouldn't need limiters.
            w_uH(thisP) = onehalfR * length(thisP)
            w_dH(thisP) = onehalfR * length(thisP)

        endwhere

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%update)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine update_interpweights_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine update_energyhead_CC (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% computes the energy head (H + v^2/2g) on each element for CC
        !% tim-marching elements
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisP(:)
        !%------------------------------------------------------------------

        elemR(thisP,er_EnergyHead) = elemR(thisP,er_Head) &
            + (elemR(thisP,er_Velocity)**2) / (twoR * setting%Constant%gravity)    

    end subroutine update_energyhead_CC
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module update
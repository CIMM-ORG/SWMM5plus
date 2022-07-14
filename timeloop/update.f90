module update

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry
    use adjust
    use utility_profiler
    use utility_crash
    use utility, only: util_CLprint, util_syncwrite

    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% Updates values during timeloop of hydraulics.
    !%

    private

    public :: update_auxiliary_variables
    !public :: update_Froude_number_junction_branch

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
    !%
    subroutine update_auxiliary_variables (whichTM)
        !%------------------------------------------------------------------
        !% Description:
        !% Updates the variables dependent on the TM solution of volume
        !% and velocity
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: whichTM  !% indicates which Time marching sets (ALLtm, AC, ETM)
            integer, pointer :: thisCol_CC, thisCol_JM
            character(64) :: subroutine_name = 'update_auxiliary_variables'
        !%------------------------------------------------------------------
        !% Preliminaries:
            !if (crashYN) return
            if (setting%Debug%File%update) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            if (setting%Profile%useYN) call util_profiler_start (pfc_update_auxiliary_variables)    
        !%------------------------------------------------------------------
        !%
            ! print *, this_image(),'    update 000 before geom TL',setting%Time%Step
            ! call util_CLprint ('in update before geometry toplevel')
        
        !% --- update the head (non-surcharged) and geometry
        call geometry_toplevel (whichTM)

            !print *, this_image(),'    update 001 after geom TL',setting%Time%Step
            ! call util_CLprint ('in update before adjust_limit_velocity_max')

        !% --- adjust velocity with limiters
        call adjust_limit_velocity_max (whichTM)
        call util_crashstop(21987)

            ! print *, this_image(),'    update 002 after adjust limit velocity',this_image()
            ! call util_CLprint ('in update before update_CC_element_flowrate')

        !% --- set packed column for updated elements
        select case (whichTM)
            case (ALLtm)
                thisCol_CC  => col_elemP(ep_CC_ALLtm)
                thisCol_JM  => col_elemP(ep_JM_ALLtm)
            case (ETM)
                thisCol_CC  => col_elemP(ep_CC_ETM)
                thisCol_JM  => col_elemP(ep_JM_ETM)
            case (AC)
                thisCol_CC  => col_elemP(ep_CC_AC)
                thisCol_JM  => col_elemP(ep_JM_AC)
            case default
                print *, 'CODE ERROR: time march type unknown for # ', whichTM
                print *, 'which has key ',trim(reverseKey(whichTM))
                call util_crashpoint(45834)
        end select

        !% --- Compute the flowrate on CC.
        !%     Note that JM should have 0 flowrate and JB has lagged flowrate at this point.
        !%     The JB flowrate is not updated until after face interpolation
        call update_CC_element_flowrate (thisCol_CC)

            ! print *, 'update 003 after element flowrate',this_image()
            ! call util_CLprint ('in update before update_Froude_number_element')

        !% --- compute element Froude numbers for CC
        call update_Froude_number_element (thisCol_CC)

            !print *, 'update 004 after Froude Number',this_image()
            !  call util_CLprint ('in update before CC interpweights in update')

        !% --- compute the wave speeds
        call update_wavespeed_element(thisCol_CC)
        call update_wavespeed_element(thisCol_JM)

        !% --- compute element-face interpolation weights on CC
        call update_interpweights_CC(thisCol_CC, whichTM)

            ! print *, 'update 005 after CC interpweights',this_image()
            ! call util_CLprint ('in update before JB interpweights')

        !% --- compute element-face interpolation weights on JB
        call update_interpweights_JB (thisCol_JM)

            ! print *, 'update 006 after JB interpweights',this_image()
            ! call util_CLprint ('in update before update Froude Number Junction Branch')

        !% --- compute element Froude number for JB
        call update_Froude_number_JB (thisCol_JM) 

            ! print *, 'update 007 at end after update Froud number JB',this_image()
            ! call util_CLprint ('in update before update BCoutlet_flowrate')

        !% --- flow values on an BC outlet face 20220714brh
        !%     required so that an inflow to a zero or small depth will not be lost
        call update_BCoutlet_flowrate ()

            ! call util_CLprint ('in update at end')

        !%------------------------------------------------------------------
        !% Closing:
            if (setting%Profile%useYN) call util_profiler_stop (pfc_update_auxiliary_variables)

             if (setting%Debug%File%update)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]" 
    end subroutine update_auxiliary_variables
    !%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine update_CC_element_flowrate (thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !%
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisCol
        !%-----------------------------------------------------------------------------
        integer, pointer ::  Npack, thisP(:)
        real(8), pointer :: flowrate(:), velocity(:), area(:)
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        flowrate => elemR(:,er_Flowrate)
        velocity => elemR(:,er_Velocity)
        area     => elemR(:,er_Area)
        !%-----------------------------------------------------------------------------
        Npack => npack_elemP(thisCol)
        if (Npack > 0) then
            thisP    => elemP(1:Npack,thisCol)
            flowrate(thisP) = area(thisP) * velocity(thisP)
        end if
        ! print*, flowrate(thisP), 'flowrate(thisP)'
    end subroutine update_CC_element_flowrate
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine update_Froude_number_element (thisCol)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes Froude number on each element
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisCol
        integer, pointer :: Npack, thisP(:)
        real(8), pointer :: Froude(:), velocity(:), depth(:), grav
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        Froude   => elemR(:,er_FroudeNumber)
        velocity => elemR(:,er_Velocity)
        depth    => elemR(:,er_ell)  !% Use the ell value (modified hydraulic depth)
        grav     => setting%constant%gravity
        !%-----------------------------------------------------------------------------

        Npack => npack_elemP(thisCol)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisCol)
            Froude(thisP) = velocity(thisP) / sqrt(grav * depth(thisP))
        end if

    end subroutine update_Froude_number_element
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine update_Froude_number_JB (thisCol_JM)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes Froude number on each junction branch element
        !% BRHbugfix 20210812
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'update_Froude_number_JB'
        integer, intent(in) :: thisCol_JM
        integer, pointer :: Npack, thisP(:), tM, BranchExists(:)
        real(8), pointer :: Froude(:), velocity(:), depth(:), grav
        integer :: ii, kk, tB
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        Froude   => elemR(:,er_FroudeNumber)
        velocity => elemR(:,er_Velocity)
        depth    => elemR(:,er_ell)  !% Use the ell value (modified hydraulic depth)
        BranchExists => elemSI(:,esi_JunctionBranch_Exists)
        grav     => setting%constant%gravity
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%update) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        Npack => npack_elemP(thisCol_JM)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisCol_JM)
            do ii=1,Npack
                tM => thisP(ii)
                do kk=1,max_branch_per_node
                    tB = tM + kk
                    if (BranchExists(tB)==1) then
                        Froude(tB) = velocity(tB) / sqrt(grav * depth(tB))
                        !print *, kk, tB, Froude(tB), velocity(tB),'  Froude JB'
                    end if
                end do
            end do
        end if

        if (setting%Debug%File%update)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine update_Froude_number_JB
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine update_wavespeed_element(thisCol)
        !%------------------------------------------------------------------
        !% Description
        !% computes the wavespeed on a CC or JM element
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisCol
            integer, pointer :: Npack, thisP(:)
            real(8), pointer :: wavespeed(:), ell(:), grav
        !%------------------------------------------------------------------
        !% Aliases:
            Npack     => npack_elemP(thisCol)
            if (Npack < 1) return
            thisP     => elemP(1:Npack,thisCol)
            wavespeed => elemR(:,er_WaveSpeed)
            ell       => elemR(:,er_ell)
            grav      => setting%constant%gravity
        !%------------------------------------------------------------------

        !% wavespeed at modified hydraulic depth (ell) 
        wavespeed(thisP) = sqrt(grav * ell(thisP))

    end subroutine update_wavespeed_element    
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine update_interpweights_CC (thisCol, whichTM)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes the interpolation weights on each element for CC
        !% tim-marching elements
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'update_interpweights_CC'
        integer, intent(in) :: thisCol, whichTM
        integer, pointer :: Npack, Npack2, thisCol_AC,  thisCol_ClosedElems, thisP(:), thisP2(:)
        real(8), pointer :: velocity(:), wavespeed(:), depth(:), length(:), QLateral(:)
        real(8), pointer :: PCelerity(:), SlotVolume(:),SlotWidth(:), fullArea(:)
        real(8), pointer :: w_uQ(:), w_dQ(:),  w_uG(:), w_dG(:),  w_uH(:), w_dH(:), w_uP(:), w_dP(:), Area(:)
        real(8), pointer :: Fr(:), grav !BRHbugfix20210811 test
        logical, pointer :: isSlot(:)
        integer :: ii
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%update) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        Qlateral  => elemR(:,er_FlowrateLateral)
        velocity  => elemR(:,er_Velocity)
        wavespeed => elemR(:,er_WaveSpeed)
        depth     => elemR(:,er_ell)  !% modified hydraulic depth!
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
        isSlot    => elemYN(:,eYN_isSlot)

        PCelerity  => elemR(:,er_Preissmann_Celerity)
        SlotVolume => elemR(:,er_SlotVolume)
        SlotWidth  => elemR(:,er_SlotWidth)
        fullArea   => elemR(:,er_FullArea)
        grav       => setting%constant%gravity


        Area => faceR(:,er_Area)
        !%-----------------------------------------------------------------------------
        !% 2nd cases needed for handling surcharged AC elements and using the celerity
        !% multiplier of the AC method for the wavespeed
        select case (whichTM)
            case (ALLtm)
                thisCol_AC          =>  col_elemP(ep_Surcharged_AC)
            case (ETM)
                thisCol_ClosedElems =>  col_elemP(ep_Closed_Elements)
            case (AC)
                thisCol_AC          =>  col_elemP(ep_Surcharged_AC)
            case default
                print *, 'CODE ERROR: time march type unknown for # ', whichTM
                print *, 'which has key ',trim(reverseKey(whichTM))
                stop 3987
        end select

        Npack => npack_elemP(thisCol)
        if (Npack < 1) return

        thisP => elemP(1:Npack,thisCol)

        ! !% wavespeed at modified hydraulic depth (ell) ! moved to update_wavespeed_element()
        ! wavespeed(thisP) = sqrt(grav * depth(thisP))
        ! ! PCelerity(thisP) = zeroR !% initialize to zero

        ! if (setting%Time%Now/3600.0 > 388.0) then
        !     print *, 'in ',trim(subroutine_name), wavespeed(ietU1(2))
        ! end if
    
        !% modify wavespeed for surcharged AC cells
        if (whichTM .ne. ETM) then
            Npack2 => npack_elemP(thisCol_AC)
            if (Npack2 > 0) then
                thisP2 => elemP(1:Npack2,thisCol_AC)
                wavespeed(thisP2) = wavespeed(thisP2) * setting%ACmethod%Celerity%RC
            end if
        ! else if (whichTM .eq. ETM) then
        !     Npack2 => npack_elemP(thisCol_ClosedElems)
        !     if (Npack2 > 0) then
        !         thisP2 => elemP(1:Npack2,thisCol_ClosedElems)
        !         !% initialize preissmann slot celerity
        !         PCelerity(thisP2) = zeroR
        !         where (SlotVolume(thisP2) .gt. zeroR) 
        !             PCelerity(thisP2) = sqrt(grav * FullArea(thisP2)/SlotWidth(thisP2)) 
        !             ! PCelerity(thisP2) = sqrt(grav * Area(thisP2)/SlotWidth(thisP2))        
        !         end where
        !     end if
        end if

        ! print *, 'in update_interpweights_CC'
        ! print *, '*** AAA vel - wave ',Velocity(iet) - wavespeed(iet)
        ! print *, '*** AAA vel + wave ',Velocity(iet) + wavespeed(iet)

        !% timescale interpolation weights for flowrate
        !% Modified from original approach by Froude number weighting
        !% Note that Fr is +/- depending on flow direction, so if the Fr is an odd power
        !% it needs to have an abs() e.g, abs(Fr(thisp)**3) *
        ! w_uQ(thisP) = - onehalfR * length(thisP)  / ( abs(Fr(thisp)**10) * velocity(thisP) - wavespeed(thisP)) !BRHbugfix 20210813 testing Fr
        ! w_dQ(thisP) = + onehalfR * length(thisP)  / ( abs(Fr(thisp)**10) * velocity(thisP) + wavespeed(thisP)) !BRHbugfix 20210813 testing Fr

        !% brh 20220104 -- for now, leave the structure of the froude number weighting, but make irrelevant with 0 power
        where (.not. isSlot(thisP))
            w_uQ(thisP) = - onehalfR * length(thisP)  / (abs(Fr(thisp)**0) * velocity(thisP) - wavespeed(thisP)) !bugfix SAZ 09212021 
            w_dQ(thisP) = + onehalfR * length(thisP)  / (abs(Fr(thisp)**0) * velocity(thisP) + wavespeed(thisP)) !bugfix SAZ 09212021 
        elsewhere (isSlot(thisP))
            w_uQ(thisP) = - onehalfR * length(thisP)  / (abs(Fr(thisp)**0) * velocity(thisP) - PCelerity(thisP)) !bugfix SAZ 23022022 
            w_dQ(thisP) = + onehalfR * length(thisP)  / (abs(Fr(thisp)**0) * velocity(thisP) + PCelerity(thisP)) !bugfix SAZ 23022022 
        end where

        ! print *, '*** BBB  w_uQ    ',w_uQ(iet)
        ! print *, '*** CCC  w_dQ    ',w_dQ(iet)

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

        ! print *, '*** DDD  w_uQ    ',w_uQ(iet)
        ! print *, '*** EEE  w_dQ    ',w_dQ(iet)

        where (w_dQ(thisP) < zeroR)
            w_dQ(thisP) = setting%Limiter%InterpWeight%Maximum
        endwhere
        where (w_dQ(thisP) < setting%Limiter%InterpWeight%Minimum)
            w_dQ(thisP) = setting%Limiter%InterpWeight%Minimum
        endwhere
        where (w_dQ(thisP) > setting%Limiter%InterpWeight%Maximum)
            w_dQ(thisP) = setting%Limiter%InterpWeight%Maximum
        endwhere

        ! print *, '*** FFF  w_uQ    ',w_uQ(iet)
        ! print *, '*** GGG  w_dQ    ',w_dQ(iet)

        !% timescale interpolation for geometry are identical to flowrate
        !% but may be modified elsewhere
        w_uG(thisP) = w_uQ(thisP)
        w_dG(thisP) = w_dQ(thisP)

        !% timescale interpolation for the preissmann number only depends on the preissmann celerity
        w_uP(thisP) = w_uQ(thisP) 
        w_dP(thisP) = w_dQ(thisP) 

        !% head uses length scale interpolation
        !% This shouldn't need limiters.
        w_uH(thisP) = onehalfR * length(thisP)
        w_dH(thisP) = onehalfR * length(thisP)

        !% adjust upstream interpolation weights for downstream flow in presence of lateral inflows
        !% so that upstream interpolation is used
        !% HACK -- this probably could use an approach with some kind of ad hoc blend -- needs work
        where ( (velocity(thisP) > zeroR) .and. (Qlateral(thisP) > zeroR) )
            w_uQ(thisP) =  setting%Limiter%InterpWeight%Maximum
            w_uG(thisP) =  setting%Limiter%InterpWeight%Maximum
            !w_uH(thisP) = setting%Limiter%InterpWeight%Minimum !do not use!
        endwhere

        ! print *, '*** HHH  w_uQ    ',w_uQ(iet)
        ! print *, '*** III  w_dQ    ',w_dQ(iet)

        ! !% adjust downstream interpolation weights for upstream flow in presence of lateral inflow
        where ( (velocity(thisP) < zeroR) .and. (Qlateral(thisP) > zeroR) )
            w_dQ(thisP) = setting%Limiter%InterpWeight%Maximum
            w_dG(thisP) = setting%Limiter%InterpWeight%Maximum
            !w_dH(thisP) = setting%Limiter%InterpWeight%Minimum ! do not use!
        endwhere

        ! print *, '*** JJJ  w_uQ    ',w_uQ(iet)
        ! print *, '*** KKK  w_dQ    ',w_dQ(iet)

        if (setting%Debug%File%update)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine update_interpweights_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine update_interpweights_JB (thisCol)
        !%------------------------------------------------------------------
        !% Description:
        !% compute the interpolation weights for junction branches
        !%------------------------------------------------------------------
            integer, intent(in) :: thisCol
            integer, pointer    :: npack, thisP(:)
            integer             :: ii
            real(8), pointer    :: grav, wavespeed(:), velocity(:), length(:), depth(:)
            real(8), pointer    :: w_uQ(:), w_dQ(:), w_uG(:), w_dG(:), w_uH(:), w_dH(:)
        !%------------------------------------------------------------------
        !% Aliases
            npack => npack_elemP(thisCol)
            if (npack < 1) return
            thisP => elemP(1:npack,thisCol)
            grav => setting%Constant%gravity
            velocity  => elemR(:,er_Velocity)
            wavespeed => elemR(:,er_WaveSpeed)
            depth     => elemR(:,er_ell)  !% modified hydraulic depth!
            length    => elemR(:,er_Length)
            w_uQ      => elemR(:,er_InterpWeight_uQ)
            w_dQ      => elemR(:,er_InterpWeight_dQ)
            w_uG      => elemR(:,er_InterpWeight_uG)
            w_dG      => elemR(:,er_InterpWeight_dG)
            w_uH      => elemR(:,er_InterpWeight_uH)
            w_dH      => elemR(:,er_InterpWeight_dH)
        !%------------------------------------------------------------------
        !% cycle through the branches to compute weights
        !print *, 'here in JB update '
        do ii=1,max_branch_per_node
            wavespeed(thisP+ii) = sqrt(grav * depth(thisP+ii))
            w_uQ(thisP+ii) = - onehalfR * length(thisP+ii)  / (velocity(thisP+ii) - wavespeed(thisP+ii))
            w_dQ(thisP+ii) = + onehalfR * length(thisP+ii)  / (velocity(thisP+ii) + wavespeed(thisP+ii))
            

            ! if (setting%Time%Now/3600.0 > 388.0) then
            !     write(*,"(A,10f16.9)") 'interp before', w_dQ(ietU1(1)), w_uQ(ietU1(2))
            !     print *, 'depth, wavespeed JB ',depth(ietU1(2)),wavespeed(ietU1(2))
            ! end if

            !% apply limiters to timescales
            where (w_uQ(thisP+ii) < zeroR)
                w_uQ(thisP+ii) = setting%Limiter%InterpWeight%Maximum
            endwhere
            where (w_uQ(thisP+ii) < setting%Limiter%InterpWeight%Minimum)
                w_uQ(thisP+ii) = setting%Limiter%InterpWeight%Minimum
            endwhere
            where (w_uQ(thisP+ii) > setting%Limiter%InterpWeight%Maximum)
                w_uQ(thisP+ii) = setting%Limiter%InterpWeight%Maximum
            endwhere

            where (w_dQ(thisP+ii) < zeroR)
                w_dQ(thisP+ii) = setting%Limiter%InterpWeight%Maximum
            endwhere
            where (w_dQ(thisP+ii) < setting%Limiter%InterpWeight%Minimum)
                w_dQ(thisP+ii) = setting%Limiter%InterpWeight%Minimum
            endwhere
            where (w_dQ(thisP+ii) > setting%Limiter%InterpWeight%Maximum)
                w_dQ(thisP+ii) = setting%Limiter%InterpWeight%Maximum
            endwhere

            ! if (setting%Time%Now/3600.0 > 388.0) then
            !     write(*,"(A,10f16.9)") 'interp after ', w_dQ(ietU1(1)), w_uQ(ietU1(2))
            ! end if

            !% set the geometry interp the same as flow interp
            w_uG(thisP+ii) = w_uQ(thisP+ii)
            w_dG(thisP+ii) = w_dQ(thisP+ii)

            !% use head interp as length-scaled
            w_uH(thisP+ii) = onehalfR * length(thisP+ii)
            w_dH(thisP+ii) = onehalfR * length(thisP+ii)  !% 20220224brh
        end do

        ! print *, ' '
        ! print *, 'in update_JB_Interpweight'
        ! write(*,"(A,10f12.5)") 'wavespeed ',wavespeed(iet)
        ! write(*,"(A,10f12.5)") 'ell       ',depth(iet)

    end subroutine update_interpweights_JB
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine update_BCoutlet_flowrate ()
        !%------------------------------------------------------------------
        !% Description:
        !% sets the outlet (face) flowrate equal to the interior element
        !% flowrate
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
    !     thisColP_dsJB  => col_elemP(ep_JB_DownStreamJB)
    !     Npack1         => npack_elemP(thisColP_dsJB)

    !     if (Npack1 > 0) then
    !         thisP1 => elemP(1:Npack1,thisColP_dsJB)
    !         w_dQ(thisP1) = oneR
    !         w_dG(thisP1) = oneR
    !         w_dH(thisP1) = oneR
    !     end if

    !     thisColP_ds_of_JB => col_elemP(ep_CC_DownstreamJbAdjacent)
    !     Npack2            => npack_elemP(ep_CC_DownstreamJbAdjacent)

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
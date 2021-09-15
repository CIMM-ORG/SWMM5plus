module update

    use define_settings, only: setting
    use define_globals
    use define_indexes
    use define_keys
    use geometry
    use adjust

    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% Updates values during timeloop of hydraulics.
    !%

    private

    public :: update_auxiliary_variables
    public :: update_Froude_number_junction_branch

    real(8), pointer :: grav => setting%constant%gravity

    contains
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    subroutine update_auxiliary_variables (whichTM)
        !%-----------------------------------------------------------------------------
        !% Description:
        !%
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: whichTM  !% indicates which Time marching method
        integer, pointer :: thisCol_all
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'update_auxiliary_variables'
        if (setting%Debug%File%update) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        !%
        !% update the head (non-surcharged) and geometry

        !print *, '---- in ',subroutine_name,'   y01'
        !write(*,'(7F9.4,A15)') elemR(ietmp,er_Head),' Head elem '

        call geometry_toplevel (whichTM)

        !print *, '---- in ',subroutine_name,'   y02'
        !write(*,'(7F9.4,A15)') elemR(ietmp,er_Head),' Head elem '

        !% adjust velocity with limiters and small volume treatment
        call adjust_velocity (whichTM, er_Velocity, er_Volume)

        !print *, '---- in ',subroutine_name,'   y03'
        !write(*,'(7F9.4,A15)') elemR(ietmp,er_Head),' Head elem '

        !% set packed column for updated elements
        select case (whichTM)
            case (ALLtm)
                thisCol_all => col_elemP(ep_CC_ALLtm)
            case (ETM)
                thisCol_all => col_elemP(ep_CC_ETM)
            case (AC)
                thisCol_all => col_elemP(ep_CC_AC)
            case default
                print *, 'error, default case should not be reached'
                stop 7489
        end select

        !% Compute the flowrate on CC.
        !% Note that JM should have 0 flowrate and JB has lagged flowrate at this point.
        !% The JB flowrate is not updated until after face interpolation
        call update_CC_element_flowrate (thisCol_all)

        !print *, '---- in ',subroutine_name,'   y04'
        !write(*,'(7F9.4,A15)') elemR(ietmp,er_Head),' Head elem '

        !% compute element Froude numbers for CC, JM
        call update_Froude_number_element (thisCol_all)

        !print *, '---- in ',subroutine_name,'   y05'
        !write(*,'(7F9.4,A15)') elemR(ietmp,er_Head),' Head elem '

        !print *, '---- in ',subroutine_name,'   y07'
        !write(*,'(7F9.4,A15)') elemR(ietmp,er_Head),' Head elem '

        !% compute element face interpolation weights on CC, JM
        call update_interpolation_weights_element (thisCol_all, whichTM)

        !print *, '---- in ',subroutine_name,'   y06'
        !write(*,'(7F9.4,A15)') elemR(ietmp,er_Head),' Head elem '

        if (setting%Debug%File%update)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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
        flowrate => elemR(:,er_Flowrate)
        velocity => elemR(:,er_Velocity)
        area     => elemR(:,er_Area)
        !%-----------------------------------------------------------------------------
        Npack => npack_elemP(thisCol)
        if (Npack > 0) then
            thisP    => elemP(1:Npack,thisCol)
            flowrate(thisP) = area(thisP) * velocity(thisP)
        end if

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
        real(8), pointer :: Froude(:), velocity(:), depth(:)
        !%-----------------------------------------------------------------------------
        Froude   => elemR(:,er_FroudeNumber)
        velocity => elemR(:,er_Velocity)
        depth    => elemR(:,er_ell)  !% Use the ell value (modified hydraulic depth)
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
    subroutine update_Froude_number_junction_branch (thisCol_JM)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes Froude number on each junction branch element
        !% BRHbugfix 20210812
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'update_Froude_number_junction_branch'
        integer, intent(in) :: thisCol_JM
        integer, pointer :: Npack, thisP(:), tM, BranchExists(:)
        real(8), pointer :: Froude(:), velocity(:), depth(:)
        integer :: ii, kk, tB
        !%-----------------------------------------------------------------------------
        Froude   => elemR(:,er_FroudeNumber)
        velocity => elemR(:,er_Velocity)
        depth    => elemR(:,er_ell)  !% Use the ell value (modified hydraulic depth)
        BranchExists => elemSI(:,eSI_JunctionBranch_Exists)
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%update) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

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
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine update_Froude_number_junction_branch
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine update_interpolation_weights_element (thisCol, whichTM)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% computes the interpolation weights on each element form CC, JM
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'update_interpolation_weights_element'
        integer, intent(in) :: thisCol, whichTM
        integer, pointer :: Npack, Npack2, thisCol_AC,  thisCol_ClosedElems, thisP(:), thisP2(:)
        real(8), pointer :: velocity(:), wavespeed(:), depth(:), length(:)
        real(8), pointer :: PCelerity(:), SlotVolume(:),SlotWidth(:), fullArea(:)
        real(8), pointer :: w_uQ(:), w_dQ(:),  w_uG(:), w_dG(:),  w_uH(:), w_dH(:)
        real(8), pointer :: Fr(:) !BRHbugfix20210811 test
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%update) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

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
        Fr        => elemR(:,er_FroudeNumber)  !BRHbugfix20210811 test

        PCelerity  => elemR(:,er_Preissmann_Celerity)
        SlotVolume => elemR(:,er_SlotVolume)
        SlotWidth  => elemR(:,er_SlotWidth)
        fullArea   => elemR(:,er_FullArea)
        !%-----------------------------------------------------------------------------
        !% 2nd cases needed for handling surcharged AC elements and using the celerity
        !% multiplier of the AC method for the wavespeed
        select case (whichTM)
            case (ALLtm)
                thisCol_AC =>  col_elemP(ep_Surcharged_AC)
            case (ETM)
                thisCol_ClosedElems => col_elemP(ep_Closed_Elements)
            case (AC)
                thisCol_AC =>  col_elemP(ep_Surcharged_AC)
            case default
                print *, 'error, case default should not be reached.'
                stop 3987
        end select

        Npack => npack_elemP(thisCol)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisCol)

            !% wavespeed at modified hydraulic depth (ell)
            wavespeed(thisP) = sqrt(grav * depth(thisP))

            !% modify wavespeed for surcharged AC cells
            if (whichTM .ne. ETM) then
                Npack2 => npack_elemP(thisCol_AC)
                if (Npack2 > 0) then
                    thisP2 => elemP(1:Npack2,thisCol_AC)
                    wavespeed(thisP2) = wavespeed(thisP2) * setting%ACmethod%Celerity%RC
                end if
            else if (whichTM .eq. ETM) then
                Npack2 => npack_elemP(thisCol_ClosedElems)
                if (Npack2 > 0) then
                    thisP2 => elemP(1:Npack2,thisCol_ClosedElems)
                    where(SlotVolume(thisP2) .gt. zeroR) 
                        PCelerity(thisP2) = sqrt(grav * fullArea(thisP2)/SlotWidth(thisP2))
                    end where
                end if
            end if


            !% timescale interpolation weights for flowrate
            !% Modified from original approach by Froude number weighting
            !% Note that Fr is +/- depending on flow direction, so if the Fr is an odd power
            !% it needs to have an abs() e.g, abs(Fr(thisp)**3) *
            w_uQ(thisP) = - onehalfR * length(thisP)  / ( abs(Fr(thisp)**10) * velocity(thisP) - wavespeed(thisP)) !BRHbugfix 20210813 testing Fr
            w_dQ(thisP) = + onehalfR * length(thisP)  / ( abs(Fr(thisp)**10) * velocity(thisP) + wavespeed(thisP)) !BRHbugfix 20210813 testing Fr

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

            !BRHbugfix 20210829
            where (w_dQ(thisP) < zeroR)
                w_dQ(thisP) = setting%Limiter%InterpWeight%Maximum
            endwhere
            where (w_dQ(thisP) < setting%Limiter%InterpWeight%Minimum)
                w_dQ(thisP) = setting%Limiter%InterpWeight%Minimum
            endwhere
            where (w_dQ(thisP) > setting%Limiter%InterpWeight%Maximum)
                w_dQ(thisP) = setting%Limiter%InterpWeight%Maximum
            endwhere
            !BRHbugfix 20210829

            !% timescale interpolation for geometry are identical to flowrate
            !% but may be modified elsewhere
            w_uG(thisP) = w_uQ(thisP)
            w_dG(thisP) = w_dQ(thisP)

            !% head uses length scale interpolation
            !% This shouldn't need limiters.

            w_uH(thisP) = onehalfR * length(thisP)
            w_dH(thisP) = onehalfR * length(thisP)

        end if

        if (setting%FaceInterp%DownJBFaceInterp == dynamic) then
            !% testin a new branch interp technique
            call update_interpolation_weights_ds_JB ()
        endif

        !print *
        !print *,'--- in ',trim(subroutine_name),' ----------------------------------------- end'
        !write(*,'(7e11.4,A15)') elemR(ietmp,er_InterpWeight_dQ),' InterpWeight_dQ'
        !write(*,'(7e11.4,A15)') elemR(ietmp,er_InterpWeight_uQ),' InterpWeight_uQ'
        ! print *, elemR(ietmp(1), er_InterpWeight_dQ)
        ! print *, elemR(ietmp(2), er_InterpWeight_dQ)
        ! print *, elemR(ietmp(3), er_InterpWeight_dQ)
        ! print *, elemR(ietmp(4), er_InterpWeight_dQ)
        ! print *, elemR(ietmp(5), er_InterpWeight_dQ)
        ! print *, elemR(ietmp(6), er_InterpWeight_dQ)
        ! print *, elemR(ietmp(7), er_InterpWeight_dQ)

        if (setting%Debug%File%update)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine update_interpolation_weights_element
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine update_interpolation_weights_ds_JB ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% This is a test subroutine that violates no neighbour algorithm
        !% This subroutine sets the interpolation wights in ds JB to its
        !% conneceted link element
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'update_interpolation_weights_ds_JB'
        integer, pointer :: thisColP_dsJB, thisColP_ds_of_JB
        integer, pointer :: Npack1, Npack2,  thisP1(:), thisP2(:)
        real(8), pointer :: w_uQ(:), w_dQ(:),  w_uG(:), w_dG(:),  w_uH(:), w_dH(:)
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%update)  &
        write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"
        w_uQ      => elemR(:,er_InterpWeight_uQ)
        w_dQ      => elemR(:,er_InterpWeight_dQ)
        w_uG      => elemR(:,er_InterpWeight_uG)
        w_dG      => elemR(:,er_InterpWeight_dG)
        w_uH      => elemR(:,er_InterpWeight_uH)
        w_dH      => elemR(:,er_InterpWeight_dH)
        !%-----------------------------------------------------------------------------

        !% replace the interpolation weights for downstream JB
        thisColP_dsJB  => col_elemP(ep_JB_DownStreamJB)
        Npack1         => npack_elemP(thisColP_dsJB)

        if (Npack1 > 0) then
            thisP1 => elemP(1:Npack1,thisColP_dsJB)
            w_dQ(thisP1) = oneR
            w_dG(thisP1) = oneR
            w_dH(thisP1) = oneR
        end if

        thisColP_ds_of_JB => col_elemP(ep_CC_DownstreamJbAdjacent)
        Npack2            => npack_elemP(ep_CC_DownstreamJbAdjacent)

        !% replace the interpolation weights for elements downstream of dn JB
        if (Npack2 > 0) then
            thisP2 => elemP(1:Npack2,thisColP_ds_of_JB)
            w_uQ(thisP2) = oneR
            w_uG(thisP2) = oneR
            w_uH(thisP2) = oneR
        end if

        if (setting%Debug%File%update) &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine update_interpolation_weights_ds_JB
    !%
    !%==========================================================================
    !% END OF MODULE
    !%==========================================================================
end module update
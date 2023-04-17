module diagnostic_elements

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use face
    use weir_elements
    use pump_elements
    use orifice_elements
    use outlet_elements
    use adjust
    !use utility, only: util_CLprint
    use utility_profiler
    use utility_crash, only: util_crashpoint
    ! use utility_unit_testing, only: util_utest_CLprint

    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% Computes diagnostic elements
    !%
    !% METHOD:
    !%
    !%

    private

    public :: diagnostic_toplevel
    public :: diagnostic_fix_JB_adjacent
    public :: diagnostic_by_type 

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine diagnostic_toplevel (elemPCol, facePCol, istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: elemPCol, facePCol, istep
        integer, pointer    :: facePackCol
        
        character(64) :: subroutine_name = 'diagnostic_toplevel'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%diagnostic_elements) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if (setting%Profile%useYN) call util_profiler_start (pfc_diagnostic_toplevel)
        !%-----------------------------------------------------------------------------
        !%
        !thisCol => col_elemP(ep_Diag)
        !Npack   => npack_elemP(elemPCol)

        !% --- if (Npack > 0) was commented out by SS
        !%     within this conditional face update has been called.
        !%     the face update syncs all the images when running in parallel
        !%     if daignostic elements are presents in a system, and if one/some
        !%     does not include a diagnostic element, it will cause a race condition.
        !%     because the images with the diagnostic elements will wait for syncs 
        !%     inside the face update. however, images that do not include a diagnostic 
        !%     element will never reach that condition.

        ! if (Npack > 0) then

            ! ! ! call util_utest_CLprint ('in diagnostic_toplevel  AAAA')

        call diagnostic_by_type (elemPCol, istep)

            ! call util_utest_CLprint ('in diagnostic_toplevel  BBB')

        !% reset any face values affected
        call face_interpolation (facePCol,.true.,.true.,.true.,.true.,.true.)

            ! call util_utest_CLprint ('in diagnostic_toplevel  CCC')

        !% --- reset the zero and small depth fluxes
        call adjust_zero_and_small_depth_face (.false.)

            ! call util_utest_CLprint ('in diagnostic_toplevel  DDD')

        ! end if
       
        if (setting%Profile%useYN) call util_profiler_stop (pfc_diagnostic_toplevel)

        if (setting%Debug%File%diagnostic_elements)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine diagnostic_toplevel
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine diagnostic_fix_JB_adjacent ()
        !%------------------------------------------------------------------
        !% Description
        !% ensures that a diagnostic element adjacent to a JB junction branch
        !% has exactly the JB flowrate on its faces and element
        !% This assumes that JB flowrate has been forced to the faces of
        !% the JB
        !%------------------------------------------------------------------
        !% Declarations
            integer, pointer :: thisColP, Npack, thisE(:), fup(:), fdn(:)
            logical, pointer :: fFrozenYN(:), fJBupstreamYN(:), fJBdownstreamYN(:)
            integer :: mm, eIdx
        !%------------------------------------------------------------------
        !% Aliases
            thisColP =>   col_elemP(ep_Diag_JBadjacent)
            Npack    => npack_elemP(thisColP)
            if (Npack < 1) return

            thisE => elemP(:,thisColP)
            fup   => elemI(:,ei_Mface_uL)
            fdn   => elemI(:,ei_Mface_dL)

            fFrozenYN       => faceYN(:,fYN_isJB_QfrozenByDiag)
            fJBupstreamYN   => faceYN(:,fYN_isUpstreamJBFace)
            fJBdownstreamYN => faceYN(:,fYN_isDownstreamJBFace)
        !%------------------------------------------------------------------

        !print *, 'in Fix JB adjacent '
        !% --- cycle through diagnostic elements adjacent to JB
        do mm=1,Npack
            eIdx = thisE(mm)

            !% --- check if faces are not frozen
            if ( (.not. fFrozenYN(fup(eIdx))) .and. (.not. fFrozenYN(fdn(eIdx))) ) then
                if (fJBupstreamYN(fup(eIdx))) then 
                    !% --- if JB upstream then store that Q as the diagnostic element
                    !%     and the downstream face
                    elemR(eIdx,er_Flowrate)      = faceR(fup(eIdx),fr_Flowrate)
                    faceR(fdn(eIdx),fr_Flowrate) = faceR(fup(eIdx),fr_Flowrate)
                elseif (fJBdownstreamYN(fdn(eIdx))) then
                    !% --- if JB downstreamstream then store that Q as the diagnostic element
                    !%     and the upstream face
                    elemR(eIdx,er_Flowrate)      = faceR(fdn(eIdx),fr_Flowrate)
                    faceR(fup(eIdx),fr_Flowrate) = faceR(fdn(eIdx),fr_Flowrate)
                else
                    print *, 'CODE ERROR: unexpected else'
                    call util_crashpoint(219874)
                end if
            end if
        end do

    end subroutine diagnostic_fix_JB_adjacent
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine diagnostic_by_type (thisCol, istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Solves for flow/head on all the diagnostic elements.
        !%
        !% Because the diagnostic elements are not vectorized by type, we simply
        !% must step through the packed array and solve each element on an individual
        !% basis. Although this is not efficient as vectorizing, for our purposes
        !% the number of diagnostic elements is small and it simply isn't worth the
        !% difficulty in storing them in vector groupings.
        !%-----------------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisCol, istep
            integer, pointer    :: thisType, thisP(:), Npack
            real(8), pointer    :: FlowRate(:)
            real(8)             :: FlowRateOld
            integer :: ii
        !%-----------------------------------------------------------------------------
            Npack => npack_elemP(thisCol)
            if (Npack < 1) return
        !%-----------------------------------------------------------------------------
        !% Aliases
            FlowRate => elemR(:,er_Flowrate)    
            thisP    => elemP(1:Npack,thisCol)
        !%-----------------------------------------------------------------------------

        !% this cycles through the individual elements, but each
        !% cycle is entirely independent
        do ii=1,Npack
            !% replace with do concurrent if every procedure called in this loop can be PURE
            thisType => elemI(thisP(ii),ei_elementType)

            !print *, 'in diagnostic by type ',ii, thisP(ii)
            !print *, 'case ',reverseKey(thisType)

            !% -- store the old flowrate for use in first step of an RK2
            FlowRateOld = FlowRate(thisP(ii))

            select case (thisType)
                
            case (weir)
                call weir_toplevel (thisP(ii))

            case (orifice)
                call orifice_toplevel (thisP(ii))

            case (pump)
                call pump_toplevel (thisP(ii))

            case (outlet)
                call outlet_toplevel (thisP(ii))
                
            case default
                print *, 'CODE ERROR element type unknown for # ', thisType
                print *, 'which has key ',trim(reverseKey(thisType))
                call util_crashpoint( 9472)
                !return
            end select

            !% --- prevent an RK2 first step from setting the flowrate to zero
            !%     Otherwise the conservative flux is identically zero for the
            !%     entire time step
            if (((istep == oneI) .or. (istep == zeroI)) .and. (FlowRate(thisP(ii)) .eq. zeroR)) then
                FlowRate(thisP(ii)) = onehalfR * (FlowRate(thisP(ii)) + FlowRateOld)
            end if
        end do

    end subroutine diagnostic_by_type
!%
!%==========================================================================
!%==========================================================================
!%
!%==========================================================================
!%==========================================================================
!%
        !%-----------------------------------------------------------------------------
        !% Description:
        !%
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%
!%
!%==========================================================================
!%==========================================================================
!%
        !%-----------------------------------------------------------------------------
        !% Description:
        !%
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%
!%
!%==========================================================================
!%==========================================================================
!%
        !%-----------------------------------------------------------------------------
        !% Description:
        !%
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module diagnostic_elements
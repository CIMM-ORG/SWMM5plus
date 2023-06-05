module diagnostic_elements
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Calls routines for diagnostic (not time-marching) elements
    !%
    !% Methods:
    !% Varies depending on element type
    !%==========================================================================
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
    use utility_profiler
    use utility_crash, only: util_crashpoint

    implicit none

    private

    public :: diagnostic_fix_JB_adjacent
    public :: diagnostic_by_type 

    contains
!%==========================================================================
!% PUBLIC
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
        !% Preliminaries:
            Npack => npack_elemP(thisCol)
            if (Npack < 1) return
        !%-----------------------------------------------------------------------------
        !% Aliases:
            FlowRate => elemR(:,er_Flowrate)    
            thisP    => elemP(1:Npack,thisCol)
        !%-----------------------------------------------------------------------------

        !% this cycles through the individual elements, but each
        !% cycle is entirely independent
        do ii=1,Npack
            !% replace with do concurrent if every procedure called in this loop can be PURE
            thisType => elemI(thisP(ii),ei_elementType)

            !% -- store the old flowrate for use in first step of an RK2
            FlowRateOld = FlowRate(thisP(ii))

            select case (thisType)
                
            case (weir)
                call weir_toplevel (thisP(ii))

            case (orifice)
                call orifice_toplevel (thisP(ii))

            case (pump)
                call pump_toplevel (thisP(ii),istep)

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
!% END OF MODULE
!%+=========================================================================
end module diagnostic_elements
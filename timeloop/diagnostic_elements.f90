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
    use utility_profiler
    use utility_crash, only: util_crashpoint

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

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine diagnostic_toplevel
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step
        !%-----------------------------------------------------------------------------
        integer, pointer :: thisCol, Npack, facePackCol
        
        character(64) :: subroutine_name = 'diagnostic_toplevel'
        !%-----------------------------------------------------------------------------
        if (crashYN) return
        if (setting%Debug%File%diagnostic_elements) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if (setting%Profile%useYN) call util_profiler_start (pfc_diagnostic_toplevel)
        !%-----------------------------------------------------------------------------
        !%
        thisCol => col_elemP(ep_Diag)
        Npack   => npack_elemP(thisCol)

        if (Npack > 0) then
            call diagnostic_by_type (thisCol, Npack)
            if (crashI==1) return
            call face_interpolation (fp_Diag, dummy)
        end if

        if (setting%Profile%useYN) call util_profiler_stop (pfc_diagnostic_toplevel)

        if (setting%Debug%File%diagnostic_elements)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine diagnostic_toplevel
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine diagnostic_by_type (thisCol, Npack)
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
        integer, intent(in) :: Npack, thisCol
        integer, pointer :: thisType, thisP(:)
        integer :: ii
        !%-----------------------------------------------------------------------------
        if (crashYN) return
        thisP => elemP(1:Npack,thisCol)

        !% this cycles through the individual elements, but each
        !% cycle is entirely independent
        do ii=1,Npack
            !% replace with do concurrent if every procedure called in this loop can be PURE
            thisType => elemI(thisP(ii),ei_elementType)

            select case (thisType)
            case (weir)
                call weir_toplevel (thisP(ii))

            case (orifice)
                call orifice_toplevel (thisP(ii))

            case (pump)
                print *, 'CODE ERROR: Pump has not yet been developed'
                print *, 'which has key ',trim(reverseKey(thisType))
                call util_crashpoint(564321)
                return

            case (outlet)
                call outlet_toplevel (thisP(ii))
                
            case default
                print *, 'CODE ERROR element type unknown for # ', thisType
                print *, 'which has key ',trim(reverseKey(thisType))
                call util_crashpoint( 9472)
                return
            end select
        end do

    end subroutine diagnostic_by_type
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
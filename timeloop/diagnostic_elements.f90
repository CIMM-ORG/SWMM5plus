module diagnostic_elements

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use face
    use weir_elements
    use pump_elements
    use orifice_elements

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
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'diagnostic_toplevel'
        if (setting%Debug%File%diagnostic_elements) print *, '*** enter ', this_image(), subroutine_name
        !%-----------------------------------------------------------------------------
        !%
        thisCol => col_elemP(ep_Diag)
        Npack   => npack_elemP(thisCol)

        if (Npack > 0) then
            call diagnostic_by_type (thisCol, Npack)
            call face_interpolation (fp_Diag)
        endif

        if (setting%Debug%File%diagnostic_elements)  print *, '*** leave ', this_image(), subroutine_name
    end subroutine diagnostic_toplevel
    ! %
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
        thisP => elemP(1:Npack,thisCol)

        !% this cycles through the individual elements, but each
        !% cycle is entirely independent
        do ii=1,Npack
            !% replace with do concurrent if every procedure called in this loop can be PURE
            thisType => elemI(thisP(ii),ei_specificType)

            select case (thisType)
                case (weir)
                    call weir_toplevel (thisP(ii))

                case (orifice)
                    call orifice_toplevel (thisP(ii))

                case (pump)
                    ! call diagnostic_pump (thisP(ii))

                case default
                    print *, 'error, default case should not be reached'
                    stop 9472
            end select
        enddo

        !% HACK not sure what we need for diagnostic aux variables
        !% The weir geometry is set in weir routines, as is flowrate, head, and velocity
        call diagnostic_auxiliary_variables (thisCol, Npack)

    end subroutine diagnostic_by_type
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine diagnostic_auxiliary_variables (thisCol, Npack)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes auxiliary variables for diagnostic elements
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisCol, Npack
!% Not sure what we need here. The diagnostic should produce the flowrate and head of the
!% diagnostic element. For face interpolation we need arguabely area, topwidth, hydraulic depth
!% However, if the interpolation weighting makes these negligible, then maybe they can just use
!% arbitrary small values.

!% QUESTION -- how do we handle geometry for the diagnostic elements themselves?

        !%-----------------------------------------------------------------------------
        !%
    end subroutine diagnostic_auxiliary_variables
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
    !    !%
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
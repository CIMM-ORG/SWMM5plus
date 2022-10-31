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
    use utility, only: util_CLprint
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
    subroutine diagnostic_toplevel (isRKfirstStep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step
        !%-----------------------------------------------------------------------------
        logical, intent(in) :: isRKfirstStep
        integer, pointer :: thisCol, Npack, facePackCol
        
        character(64) :: subroutine_name = 'diagnostic_toplevel'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%diagnostic_elements) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if (setting%Profile%useYN) call util_profiler_start (pfc_diagnostic_toplevel)
        !%-----------------------------------------------------------------------------
        !%
        thisCol => col_elemP(ep_Diag)
        Npack   => npack_elemP(thisCol)

        if (Npack > 0) then
            !print *, 'calling diagnostic by type'
            call diagnostic_by_type (thisCol, Npack, isRKfirstStep)

            !% reset any face values affected
            call face_interpolation (fp_Diag, dummy)

            !% --- reset the zero and small depth fluxes
            call adjust_zero_and_small_depth_face (ETM, .false.)

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
    subroutine diagnostic_by_type (thisCol, Npack, isRKfirstStep)
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
            logical, intent(in) :: isRKfirstStep
            integer, intent(in) :: Npack, thisCol
            integer, pointer    :: thisType, thisP(:)
            real(8), pointer    :: FlowRate(:)
            real(8)             :: FlowRateOld
            integer :: ii
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
            if ((isRKfirstStep) .and. (FlowRate(thisP(ii)) .eq. zeroR)) then
                FlowRate(thisP(ii)) = onehalfR * (FlowRate(thisP(ii)) + FlowRateOld)
            end if
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
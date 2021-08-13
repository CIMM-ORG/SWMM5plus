module discretization

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting

    implicit none

    real(8), pointer :: elem_nominal_length => setting%Discretization%NominalElemLength
    real(8), pointer :: elem_shorten_cof    => setting%Discretization%LinkShortingFactor

contains
    
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_discretization_adjustlinklength()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine loans some of the length of a link to an adjacent nJm node.
    !   The purpose is to allow the nJm branch elements to have some volume.
    !
    ! HACK:
    !   this is a subroutine for adjusting the length of links.
    !   Put it here for now but can be moved to somewhere else
    !
    !-----------------------------------------------------------------------------
        integer :: ii, Adjustment_flag
        real(8) :: temp_length

        character(64) :: subroutine_name = 'init_discretization_adjustlinklength'
    !-----------------------------------------------------------------------------

        if (setting%Debug%File%discretization) print *, '*** enter ', subroutine_name

        do ii =1, N_link
            temp_length = link%R(ii,lr_Length) ! lenght of link ii
            Adjustment_flag = oneI

            if ( node%I(link%I(ii,li_Mnode_u), ni_node_type) == nJm ) then
                temp_length = temp_length - elem_shorten_cof * elem_nominal_length ! make a cut for upstream M junction
                Adjustment_flag = Adjustment_flag + oneI
            endif

            if ( node%I(link%I(ii,li_Mnode_d), ni_node_type) == nJm ) then
                temp_length = temp_length - elem_shorten_cof * elem_nominal_length ! make a cut for downstream M junction
                Adjustment_flag = Adjustment_flag + oneI
            endif

            link%R(ii,lr_AdjustedLength) = temp_length
            link%I(ii,li_length_adjusted) = Adjustment_flag
            !% set the new element length based on the adjusted link length
            link%R(ii,lr_ElementLength) = link%R(ii,lr_AdjustedLength)/link%I(ii,li_N_element)
        enddo

        if (setting%Debug%File%discretization)  print *, '*** leave ', subroutine_name
    end subroutine init_discretization_adjustlinklength
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_discretization_nominal()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine sets the number of elements per link.  The element length
    !   is adjusted so that an integer number of elements is assigned to each link.
    !
    !-----------------------------------------------------------------------------
        integer :: ii
        real(8) :: remainder
        character(64) :: subroutine_name = 'init_discretization_nominal'
               
    !-----------------------------------------------------------------------------

        if (setting%Debug%File%discretization) print *, '*** enter ', subroutine_name

        do ii = 1, N_link
            remainder = mod(link%R(ii,lr_Length), elem_nominal_length)
            if ( remainder == zeroR ) then
                link%I(ii, li_N_element) = int(link%R(ii, lr_Length)/elem_nominal_length)
                link%R(ii, lr_ElementLength) = link%R(ii, lr_Length)/link%I(ii, li_N_element)
            elseif ( remainder .ge. onehalfR * elem_nominal_length ) then
                link%I(ii, li_N_element) = ceiling(link%R(ii,lr_Length)/elem_nominal_length)
                link%R(ii, lr_ElementLength) = link%R(ii, lr_Length)/link%I(ii, li_N_element)
            else
                link%I(ii, li_N_element) = floor(link%R(ii,lr_Length)/elem_nominal_length)
                link%R(ii, lr_ELementLength) = link%R(ii, lr_Length)/link%I(ii, li_N_element)
            endif
        enddo

        if (setting%Debug%File%discretization)  print *, '*** leave ', subroutine_name

    end subroutine init_discretization_nominal
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_discretization_cfl()
    end subroutine init_discretization_cfl
    !
    !==========================================================================
    ! END OF MODULE
    !==========================================================================
    !
end module discretization

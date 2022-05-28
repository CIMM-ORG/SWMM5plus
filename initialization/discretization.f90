module discretization

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting

    implicit none

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
        integer          :: ii, Adjustment_flag
        real(8)          :: temp_length
        real(8), pointer :: elem_nominal_length
        real(8)          :: elem_shorten_cof
        character(64)    :: subroutine_name = 'init_discretization_adjustlinklength'
    !-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%discretization) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        elem_nominal_length => setting%Discretization%NominalElemLength
        
        do ii =1, N_link
            !% --- default shorting coefficient (reset for each link)
            elem_shorten_cof = setting%Discretization%LinkShortingFactor

            temp_length = link%R(ii,lr_Length) ! length of link ii
            Adjustment_flag = oneI
           
            !% --- adjust the shortening for small link lengths 20220520brh 
            if (temp_length < (oneR + twoR*elem_shorten_cof) * elem_nominal_length) then
                !% --- limit the shortening to 1/4 of the total element length (1/8 on either side)
                elem_shorten_cof = temp_length / (elem_nominal_length * eightR)
            end if

            if ( node%I(link%I(ii,li_Mnode_u), ni_node_type) == nJm ) then
                temp_length = temp_length - elem_shorten_cof * elem_nominal_length ! make a cut for upstream M junction
                Adjustment_flag = Adjustment_flag + oneI
            end if

            if ( node%I(link%I(ii,li_Mnode_d), ni_node_type) == nJm ) then
                temp_length = temp_length - elem_shorten_cof * elem_nominal_length ! make a cut for downstream M junction
                Adjustment_flag = Adjustment_flag + oneI
            end if

            ! if (ii==30) then
            !     print *, 'here in ',trim(subroutine_name)
            !     print *, temp_length
            !     if ((temp_length < onehalfR * link%R(ii,lr_Length)) &
            !        .and. (temp_length < elem_nominal_length) ) then
            !         print *, temp_length,  link%R(ii,lr_Length)
            !     end if
            !     print *, temp_length,  link%R(ii,lr_Length)
            !     stop 29873
            ! end if

            if ((link%I(ii,li_link_type) == lChannel) .or. (link%I(ii,li_link_type) == lPipe)) then
                link%R(ii,lr_AdjustedLength) = temp_length
                link%I(ii,li_length_adjusted) = Adjustment_flag
                !% set the new element length based on the adjusted link length
                link%R(ii,lr_ElementLength) = link%R(ii,lr_AdjustedLength)/link%I(ii,li_N_element)
            else
                link%R(ii,lr_AdjustedLength) = link%R(ii,lr_ElementLength)
                link%I(ii,li_length_adjusted) = DiagAdjust
            end if

        end do

        if (setting%Debug%File%discretization)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_discretization_adjustlinklength
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_discretization_nominal(link_idx)
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine sets the number of elements per link.  The element length
    !   is adjusted so that an integer number of elements is assigned to each link.
    !
    !-----------------------------------------------------------------------------
        integer, intent(in) :: link_idx
        real(8) :: remainder
        real(8), pointer :: elem_nominal_length
        character(64) :: subroutine_name = 'init_discretization_nominal'

    !-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%discretization) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        elem_nominal_length => setting%Discretization%NominalElemLength

        !% Adjusts the number of elements in a link based on the length
        remainder = mod(link%R(link_idx,lr_Length), elem_nominal_length)

        !% If the elements fit evenly into the length of link
        if ( remainder == zeroR ) then
            link%I(link_idx, li_N_element) = int(link%R(link_idx, lr_Length)/elem_nominal_length)
            link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length)/link%I(link_idx, li_N_element)

        !% If the remainder is greater than half an element length
        elseif ( remainder .ge. onehalfR * elem_nominal_length ) then
            link%I(link_idx, li_N_element) = ceiling(link%R(link_idx,lr_Length)/elem_nominal_length)
            link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length)/link%I(link_idx, li_N_element)

        !% If the remainder is less than half an element length
        else
            link%I(link_idx, li_N_element) = floor(link%R(link_idx,lr_Length)/elem_nominal_length)
            link%R(link_idx, lr_ELementLength) = link%R(link_idx, lr_Length)/link%I(link_idx, li_N_element)
        end if

        !% Additional check to ensure that every link has at least one element
        if ( link%R(link_idx, lr_Length) .le. elem_nominal_length ) then
            link%I(link_idx, li_N_element) = oneI
            link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length)
        end if

        !% treatment of for special links
        if ((link%I(link_idx,li_link_type) == lWeir)    .or. &
            (link%I(link_idx,li_link_type) == lOrifice) .or. &
            (link%I(link_idx,li_link_type) == lOutlet)  .or. &
            (link%I(link_idx,li_link_type) == lPump) ) then

            link%I(link_idx, li_N_element) = oneI
            link%R(link_idx, lr_ElementLength) = link%R(link_idx, lr_Length)

        end if 

        ! !% HARD CODED for trajkovic case A_nx40
        ! link%I(1, li_N_element) = 6
        ! link%I(2, li_N_element) = 40
        ! link%I(3, li_N_element) = 6
        ! link%I(4, li_N_element) = 20

        !% HARD CODED for trajkovic case A_nx80
        ! link%I(1, li_N_element) = 12
        ! link%I(2, li_N_element) = 80
        ! link%I(3, li_N_element) = 12
        ! link%I(4, li_N_element) = 40

        if (setting%Debug%File%discretization)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

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

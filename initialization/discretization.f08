module discretization

    use globals
    use array_index
    use setting_definition, only: setting

    implicit none

contains

    ! this is a subroutine for adjusting the length of links.
    ! Put it here for now but can be moved to somewhere else
    subroutine adjust_link_length()
        integer :: ii
        real(8) :: temp_length
        character(64) :: subroutine_name = 'adjust_link_length'

        if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name

        do ii =1, N_link
            temp_length = linkR(ii,lr_Length) ! lenght of link ii

            if ( nodeI(linkI(ii,li_Mnode_u), ni_node_type) .eq. nJm ) then
                temp_length = temp_length - elem_branch_factor * elem_nominal_length ! make a cut for upstream M junction
            endif

            if ( nodeI(linkI(ii,li_Mnode_d), ni_node_type) .eq. nJm ) then
                temp_length = temp_length - elem_branch_factor * elem_nominal_length ! make a cut for downstream M junction
            endif

            linkR(ii,lr_Length) = temp_length
        enddo

        if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name
    end subroutine adjust_link_length

    subroutine nominal_discretization()
        integer :: ii
        real(8) :: remainder
        character(64) :: subroutine_name = 'nominal_discretization'

        if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name

        do ii = 1, N_link
            remainder = mod(linkR(ii,lr_Length), elem_nominal_length)
            if ( remainder .eq. zeroR ) then
                linkI(ii, li_N_element) = int(linkR(ii, lr_Length)/elem_nominal_length)
                linkR(ii, lr_ElementLength) = linkR(ii, lr_Length)/linkI(ii, li_N_element)
            elseif ( remainder .ge. onehalfR * elem_nominal_length ) then
                linkI(ii, li_N_element) = ceiling(linkR(ii,lr_Length)/elem_nominal_length)
                linkR(ii, lr_ElementLength) = linkR(ii, lr_Length)/linkI(ii, li_N_element)
            else
                linkI(ii, li_N_element) = floor(linkR(ii,lr_Length)/elem_nominal_length)
                linkR(ii, lr_ELementLength) = linkR(ii, lr_Length)/linkI(ii, li_N_element)
            endif
        enddo

        if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name

    end subroutine nominal_discretization

    subroutine cfl_discretization()
    end subroutine cfl_discretization

end module discretization

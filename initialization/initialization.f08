module initialization
    use allocate_storage
    use array_index
    use data_keys
    use globals
    use interface
    use BIPquick
    use coarray_partition
    use utility, only: utility_export_linknode_csv
    use setting_definition, only: setting

    implicit none


!-----------------------------------------------------------------------------
!
! Description:
!   General initialization of data structures (not including network)
!
! Method:
!    Creates the arrays index structures that are used for accessing data.
!    Arguably, this could be done more simply but we want the fundamental
!    column indexes in array_index to be parameters rather than variables. By
!    using parameters we reduce the possibility of accidentally changing a
!    column definition.
!
!-----------------------------------------------------------------------------

    private

    public :: initialize_linknode_arrays, count_node_types

contains

    subroutine initialize_linknode_arrays()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   Retrieves data from SWMM C interface and populates link and node tables
    !
    !-----------------------------------------------------------------------------

        integer       :: i, total_n_links
        logical       :: l1, l2
        character(64) :: subroutine_name = 'initialize_arrays'

    !-----------------------------------------------------------------------------

        if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name

        if (.not. api_is_initialized) then
            print *, "ERROR: API is not initialized"
            stop
        end if

        ! Allocate storage for link & node tables
        call allocate_linknode_storage ()

        nodeI(:,ni_N_link_u) = 0
        nodeI(:,ni_N_link_d) = 0

        do i = 1, N_link
            linkI(i,li_idx) = i
            linkI(i,li_link_type) = get_link_attribute(i, link_type)
            linkI(i,li_geometry) = get_link_attribute(i, link_geometry)
            linkI(i,li_Mnode_u) = get_link_attribute(i, link_node1) + 1 ! node1 in C starts from 0
            linkI(i,li_Mnode_d) = get_link_attribute(i, link_node2) + 1 ! node2 in C starts from 0

            ! HACK This is a temporary hardcode until Gerardo can populate this column from the CFL condition
            linkI(i, li_N_element) = 10

            nodeI(linkI(i,li_Mnode_u), ni_N_link_u) = nodeI(linkI(i,li_Mnode_u), ni_N_link_u) + 1
            nodeI(linkI(i,li_Mnode_u), ni_idx_base1 + nodeI(linkI(i,li_Mnode_u), ni_N_link_u)) = i
            nodeI(linkI(i,li_Mnode_d), ni_N_link_d) = nodeI(linkI(i,li_Mnode_d), ni_N_link_d) + 1
            nodeI(linkI(i,li_Mnode_d), ni_idx_base2 + nodeI(linkI(i,li_Mnode_d), ni_N_link_d)) = i

            linkI(i,li_InitialDepthType) = 1 ! TODO - get from params file
            linkR(i,lr_Length) = get_link_attribute(i, conduit_length)

            ! linkR(i,lr_TopWidth): defined in network_define.f08
            linkR(i,lr_BreadthScale) = get_link_attribute(i, link_xsect_wMax)
            ! linkR(i,lr_Slope): defined in network_define.f08
            linkR(i,lr_LeftSlope) = get_link_attribute(i, link_left_slope)
            linkR(i,lr_RightSlope) = get_link_attribute(i, link_right_slope)
            linkR(i,lr_Roughness) = get_link_attribute(i, conduit_roughness)
            linkR(i,lr_InitialFlowrate) = get_link_attribute(i, link_q0)
            linkR(i,lr_InitialUpstreamDepth) = get_node_attribute(linkI(i,li_Mnode_u), node_initDepth)
            linkR(i,lr_InitialDnstreamDepth) = get_node_attribute(linkI(i,li_Mnode_d), node_initDepth)
            linkR(i,lr_InitialDepth) = (linkR(i,lr_InitialDnstreamDepth) + linkR(i,lr_InitialUpstreamDepth)) / 2.0
        end do

        do i = 1, N_node
            total_n_links = nodeI(i,ni_N_link_u) + nodeI(i,ni_N_link_d)
            nodeI(i, ni_idx) = i
            if (get_node_attribute(i, node_type) == 1) then ! OUTFALL
                nodeI(i, ni_node_type) = nBCdn
            else if (total_n_links == 2) then
                nodeI(i, ni_node_type) = nJ2
            else if (total_n_links > 2) then
                nodeI(i, ni_node_type) = nJm
            end if
            ! Nodes with nBCup are defined in inflow.f08 -> (inflow_load_inflows)
            l1 = get_node_attribute(i, node_has_extInflow) == 1
            l2 = get_node_attribute(i, node_has_dwfInflow) == 1
            if (l1 .or. l2) then
                !nodeYN(i, nYN_has_inflow) = .true.
                if (total_n_links == 1) then
                    nodeI(i, ni_node_type) = nBCup
                end if
            end if


            nodeR(i,nr_InitialDepth) = get_node_attribute(i, node_initDepth)
            nodeR(i,nr_Zbottom) = get_node_attribute(i, node_invertElev)
        end do

        ! adjust the length and calculate the number/length of elements in each link
        call link_length_adjust()
        call N_elem_assign()

        
        call initialize_partition_coarray()


        if (setting%Debug%File%initialization) then
            call utility_export_linknode_csv()
        end if

        if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name
    end subroutine initialize_linknode_arrays


    ! this is a subroutine for adjusting the length of links.
    ! Put it here for now but can be moved to somewhere else
    subroutine link_length_adjust()
        integer :: ii
        real(8) :: temp_length
        character(64) :: subroutine_name = 'link_length_adjust'
        
        if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name

        do ii =1, N_link
            temp_length = linkR(ii,lr_Length) ! lenght of link ii
            
            if ( nodeI(linkI(ii,li_Mnode_u), ni_node_type) .eq. nJm ) then
                temp_length = temp_length - elem_shorten_cof * element_length ! make a cut for upstream M junction
            endif

            if ( nodeI(linkI(ii,li_Mnode_d), ni_node_type) .eq. nJm ) then
                temp_length = temp_length - elem_shorten_cof * element_length ! make a cut for downstream M junction
            endif

            linkR(ii,lr_Length) = temp_length
        enddo

        if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name
    end subroutine link_length_adjust

    subroutine N_elem_assign()
        integer :: ii
        real(8) :: remainder
        character(64) :: subroutine_name = 'N_elem_assign'
        
        if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name

        do ii = 1, N_link
            remainder = mod(linkR(ii,lr_Length), element_length)
            if ( remainder .eq. zeroR ) then
                linkI(ii, li_N_element) = int(linkR(ii, lr_Length)/element_length)
                linkR(ii, lr_ElementLength) = linkR(ii, lr_Length)/linkI(ii, li_N_element)
            elseif ( remainder .ge. onehalfR * element_length ) then
                linkI(ii, li_N_element) = ceiling(linkR(ii,lr_Length)/element_length)
                linkR(ii, lr_ElementLength) = linkR(ii, lr_Length)/linkI(ii, li_N_element)
            else
                linkI(ii, li_N_element) = floor(linkR(ii,lr_Length)/element_length)
                linkR(ii, lr_ELementLength) = linkR(ii, lr_Length)/linkI(ii, li_N_element)
            endif
        enddo

        if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name

    end subroutine N_elem_assign

    subroutine initialize_partition_coarray()
        character(64) :: subroutine_name = 'initialize_partition'
        
        if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name

        !% In order to keep the main() clean, move the following two subroutines here, BIPquick can be removed 
        call BIPquick_YJunction_Hardcode()
        
        call coarray_length_calculation()
        
        call coarray_storage_allocation()  ! once we finish the image flag this is ready to use

        if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name

    end subroutine initialize_partition_coarray


    subroutine count_node_types(N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2)
        integer, intent(in out) :: N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2
        integer :: ii
    
        ! This subroutine uses the vectorized count() function to search the array for number of instances of each node type
        N_nBCup = count(nodeI(:, ni_node_type) == nBCup)
        N_nBCdn = count(nodeI(:, ni_node_type) == nBCdn)
        N_nJm = count(nodeI(:, ni_node_type) == nJM)
        N_nStorage = count(nodeI(:, ni_node_type) == nStorage)
        N_nJ2 = count(nodeI(:, ni_node_type) == nJ2)
    
        ! The nodes that correspond to having 7, 1, and 0 attributed elements are summed together
        ! num_nJm_nodes = N_nJm
        ! num_one_elem_nodes = N_nBCup + N_nBCdn + N_nStorage
        ! num_zero_elem_nodes = N_nJ2
    
    end subroutine count_node_types

end module initialization


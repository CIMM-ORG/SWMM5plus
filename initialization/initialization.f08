module initialization

    use allocate_storage
    use array_index
    use data_keys
    use globals
    use interface
    use utility, only: utility_export_linknode_csv
    use setting_definition, only: setting

    implicit none

!-----------------------------------------------------------------------------
!
! Description:
!    General initialization of data structures (not including network)
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

        integer       :: ii, total_n_links
        logical       :: l1, l2
        character(64) :: subroutine_name = 'initialize_arrays'

    !-----------------------------------------------------------------------------

        if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name

        if (.not. api_is_initialized) then
            print *, "ERROR: API is not initialized"
            stop
        end if

        nodeI(:,ni_N_link_u) = 0
        nodeI(:,ni_N_link_d) = 0

        do ii = 1, N_link
            linkI(ii,li_idx) = ii
            linkI(ii,li_link_type) = get_link_attribute(ii, link_type)
            linkI(ii,li_geometry) = get_link_attribute(ii, link_geometry)
            linkI(ii,li_Mnode_u) = get_link_attribute(ii, link_node1) + 1 ! node1 in C starts from 0
            linkI(ii,li_Mnode_d) = get_link_attribute(ii, link_node2) + 1 ! node2 in C starts from 0

            nodeI(linkI(ii,li_Mnode_u), ni_N_link_u) = nodeI(linkI(ii,li_Mnode_u), ni_N_link_u) + 1
            nodeI(linkI(ii,li_Mnode_u), ni_idx_base1 + nodeI(linkI(ii,li_Mnode_u), ni_N_link_u)) = ii
            nodeI(linkI(ii,li_Mnode_d), ni_N_link_d) = nodeI(linkI(ii,li_Mnode_d), ni_N_link_d) + 1
            nodeI(linkI(ii,li_Mnode_d), ni_idx_base2 + nodeI(linkI(ii,li_Mnode_d), ni_N_link_d)) = ii

            linkI(ii,li_InitialDepthType) = 1 ! TODO - get from params file
            linkR(ii,lr_Length) = get_link_attribute(ii, conduit_length)
            ! linkR(ii,lr_TopWidth): defined in network_define.f08
            linkR(ii,lr_BreadthScale) = get_link_attribute(ii, link_xsect_wMax)
            ! linkR(ii,lr_Slope): defined in network_define.f08
            linkR(ii,lr_LeftSlope) = get_link_attribute(ii, link_left_slope)
            linkR(ii,lr_RightSlope) = get_link_attribute(ii, link_right_slope)
            linkR(ii,lr_Roughness) = get_link_attribute(ii, conduit_roughness)
            linkR(ii,lr_InitialFlowrate) = get_link_attribute(ii, link_q0)
            linkR(ii,lr_InitialUpstreamDepth) = get_node_attribute(linkI(ii,li_Mnode_u), node_initDepth)
            linkR(ii,lr_InitialDnstreamDepth) = get_node_attribute(linkI(ii,li_Mnode_d), node_initDepth)
            linkR(ii,lr_InitialDepth) = (linkR(ii,lr_InitialDnstreamDepth) + linkR(ii,lr_InitialUpstreamDepth)) / 2.0
        end do

        do ii = 1, N_node
            total_n_links = nodeI(ii,ni_N_link_u) + nodeI(ii,ni_N_link_d)
            nodeI(ii, ni_idx) = ii
            if (get_node_attribute(ii, node_type) == 1) then ! OUTFALL
                nodeI(ii, ni_node_type) = nBCdn
            else if (total_n_links == 2) then
                nodeI(ii, ni_node_type) = nJ2
            else if (total_n_links > 2) then
                nodeI(ii, ni_node_type) = nJm
            end if

            ! Determine if node has inflow and is nBCup
            l1 = get_node_attribute(ii, node_has_extInflow) == 1
            l2 = get_node_attribute(ii, node_has_dwfInflow) == 1
            if (l1 .or. l2) then
                nodeYN(ii, nYN_has_inflow) = .true.
                if (total_n_links == 1) then
                    nodeI(ii, ni_node_type) = nBCup
                end if
            endif

            nodeR(ii,nr_InitialDepth) = get_node_attribute(ii, node_initDepth)
            nodeR(ii,nr_Zbottom) = get_node_attribute(ii, node_invertElev)
        end do

        ! Count number of nodes with inflow
        N_inflow = count(nodeYN(:, nYN_has_inflow))
        if (setting%Debug%File%initialization) then
            call utility_export_linknode_csv()
        end if

        if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name
    end subroutine initialize_linknode_arrays

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

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

    public :: initialize_linknode_arrays

contains

    subroutine initialize_linknode_arrays()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   Retrieves data from SWMM C interface and populates link and node tables
    !
    !-----------------------------------------------------------------------------

        integer       :: i, total_n_links
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
            nodeR(i,nr_InitialDepth) = get_node_attribute(i, node_initDepth)
            nodeR(i,nr_Zbottom) = get_node_attribute(i, node_invertElev)
        end do

        if (setting%Debug%File%initialization) then
            call utility_export_linknode_csv()
        end if

        if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name
    end subroutine initialize_linknode_arrays

end module initialization

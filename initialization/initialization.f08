! module initialization
!
! This creates the arrays index structures that are used for accessing
! data. Arguably, this could be done more simply but we want the fundamental
! column indexes in array_index to be parameters rather than variables. By
! using parameters we reduce the possibility of accidentally changing a
! column definition.
!
!==========================================================================
!
module initialization
    !
    ! general initialization of data structures (not including network)
    !
    use allocate_storage
    use array_index
    use data_keys
    use globals
    use interface

    implicit none
    private

    public :: initialize_arrayindex ! handles indexes for multiple face per element
    public :: initialize_arrayindex_status ! status check at end of simulation
    public :: initialize_array_zerovalues ! sets some elemMR values to zero
    public :: initialize_dummy_values
    public :: initialize_linknode_arrays ! Retrieves data from SWMM C interface and populates link and node tables
    integer, private :: debuglevel = 0

contains

    subroutine initialize_linknode_arrays()
        character(64) :: subroutine_name
        integer :: i, total_n_links

        subroutine_name = 'initialize_arrays'

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ', subroutine_name

        if (.not. api_is_initialized) then
            print *, MSG_API_NOT_INITIALIZED
            stop
        end if

        N_link = num_links
        N_node = num_nodes

        ! Allocate storage for link & node tables
        call allocate_linknode_storage &
            (linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName)

        nodeI(:,ni_N_link_u) = 0
        nodeI(:,ni_N_link_d) = 0

        do i = 1, N_link
            linkI(i,li_idx) = i
            linkI(i,li_link_type) = get_link_attribute(i, link_type)
            linkI(i,li_geometry) = get_link_attribute(i, link_geometry)
            linkI(i,li_roughness_type) = 1 ! TODO - get from params file
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

        if ((debuglevel > 0) .or. (debuglevel > 0)) then
            print*, "li_idx, ", "li_link_type, ", "li_weir_type, ", "li_orif_type, ", "li_pump_type, ", "li_geometry, ",&
                    "li_roughness_type, ", "li_N_element, ", "li_Mnode_u, ", "li_Mnode_d, ", "li_Melem_u, ", &
                    "li_Melem_d, ", "li_Mface_u, ", "li_Mface_d, ", "li_assigned, ", "li_InitialDepthType, ", &
                    "li_temp1, ", "li_idx_max"
            do i = 1, N_link
                print *, linkI(i,:)
            end do
            print *, "lr_Length, ", "lr_BreadthScale, ", "lr_TopWidth, ", "lr_ElementLength, ", "lr_Slope, ",&
                     "lr_LeftSlope, ", "lr_RightSlope, ", "lr_Roughness, ", "lr_InitialFlowrate, ", &
                     "lr_InitialDepth, ", "lr_InitialUpstreamDepth, ", "lr_InitialDnstreamDepth, ", &
                     "lr_ParabolaValue, ", "lr_SideSlope, ", "lr_InletOffset, ", "lr_DischargeCoeff1, ",&
                     "lr_DischargeCoeff2, ", "lr_FullDepth, ", "lr_EndContractions, ", "lr_temp1"
            do i = 1, N_link
                print *, linkR(i,:)
            end do
            print *, "ni_idx, ", "ni_node_type, ", "ni_N_link_u, ", "ni_N_link_d, ", "ni_curve_type, ", &
                     "ni_assigned, ", "ni_temp1, ", "ni_idx_base1, ", "ni_Mlink_u1, ", "ni_Mlink_u2, ", &
                     "ni_Mlink_u3, ", "ni_idx_base2, ", "ni_Mlink_d1, ", "ni_Mlink_d2, ", "ni_Mlink_d3"
            do i = 1, N_node
                print *, nodeI(i,:)
            end do
            print *, "nr_Zbottom, ", "nr_InitialDepth, ", "nr_FullDepth, ", "nr_StorageConstant, ", &
                     "nr_StorageCoeff, ", "nr_StorageExponent, ", "nr_PondedArea, ", "nr_SurchargeDepth, ", &
                     "nr_temp1, ", "nr_idx_base1, ", "nr_ElementLength_u1, ", "nr_ElementLength_u2, ", &
                     "nr_ElementLength_u3, ", "nr_idx_base2, ", "nr_ElementLength_d1, ", "nr_ElementLength_d2, ", &
                     "nr_ElementLength_d3"
            do i = 1, N_node
                print *, nodeR(i,:)
            end do
        end if

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name
    end subroutine initialize_linknode_arrays
end module initialization

! module case_simple_weir
!
! Custom test case for a simple weir created with with 4 nodes and 3 links
! Used in conjunction with test_cases module for custom setup.
!
!==========================================================================
!
module case_simple_weir
    !
    use allocate_storage
    use array_index
    use bc
    use data_keys
    use globals
    use setting_definition

    implicit none

    private

    public :: case_simple_weir_initialize

    integer :: debuglevel = 0

contains
    !
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine case_simple_weir_initialize &
        (channel_length, channel_breadth, subdivide_length, lowerZ, upperZ, &
        initial_flowrate, depth_upstream, depth_dnstream, initial_depth,    &
        side_slope, inlet_offset, discharge_coefficient1,                   &
        discharge_coefficient2, full_depth, end_contractions, ManningsN,    &
        roughness_type, idepth_type, linkR, nodeR, linkI, nodeI,linkYN,     &
        nodeYN, linkName, nodeName, bcdataDn, bcdataUp)
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   initialize the link-node system and boundary conditions for a simple channel
        !
        ! Method:
        !    
        !-----------------------------------------------------------------------------
        
        character(64) :: subroutine_name = 'case_simple_weir_initialize'

        real(8),  intent(in)  :: channel_length(:), channel_breadth(:), subdivide_length(:)
        real(8),  intent(in)  :: lowerZ(:), upperZ(:),  initial_flowrate(:)
        real(8),  intent(in)  :: depth_upstream(:), depth_dnstream(:), initial_depth(:)
        real(8),  intent(in)  :: side_slope(:), inlet_offset(:), end_contractions(:)
        real(8),  intent(in)  :: discharge_coefficient1(:), discharge_coefficient2(:)
        real(8),  intent(in)  :: full_depth(:), ManningsN(:)

        integer, intent(in):: roughness_type, idepth_type(:)

        integer,   dimension(:,:), allocatable, target, intent(out)    :: linkI
        integer,   dimension(:,:), allocatable, target, intent(out)    :: nodeI

        real(8),      dimension(:,:), allocatable, target, intent(out)    :: linkR
        real(8),      dimension(:,:), allocatable, target, intent(out)    :: nodeR

        logical,   dimension(:,:), allocatable, target, intent(out)    :: linkYN
        logical,   dimension(:,:), allocatable, target, intent(out)    :: nodeYN

        type(string), dimension(:), allocatable, target, intent(out)   :: linkName
        type(string), dimension(:), allocatable, target, intent(out)   :: nodeName

        type(bcType), dimension(:), allocatable, intent(out) :: bcdataUp, bcdataDn

        integer    :: ntimepoint, N_BCdnstream, N_BCupstream

        integer            :: allocation_status, ii
        character(len=99)  :: emsg

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        ! Boundary conditions
        ntimepoint = 2
        N_BCupstream = 1
        N_BCdnstream = 1

        ! check if

        call bc_allocate &
            (bcdataDn, bcdataUp)

        do ii = 1, N_BCdnstream
            allocate( bcdataDn(ii)%TimeArray(ntimepoint))
            allocate( bcdataDn(ii)%ValueArray(ntimepoint))
            bcdataDn(ii)%TimeArray      = nullvalueR
            bcdataDn(ii)%ValueArray     = nullvalueR
        end do
        do ii = 1, N_BCupstream
            allocate( bcdataUp(ii)%TimeArray(ntimepoint))
            allocate( bcdataUp(ii)%ValueArray(ntimepoint))
            bcdataUp(ii)%TimeArray      = nullvalueR
            bcdataUp(ii)%ValueArray     = nullvalueR
        end do

        ! assign values
        ! downstream is default to elevation
        bcdataDn(1)%NodeID = 1
        bcdataDn(1)%TimeArray(1)     = setting%Time%StartTime
        bcdataDn(1)%TimeArray(2)     = setting%Time%EndTime + 100.0 !s
        bcdataDn(1)%ValueArray(1)    = lowerZ(1) +  depth_dnstream(1)   ! m
        bcdataDn(1)%ValueArray(2)    = lowerZ(1) +  depth_dnstream(1) ! m

        ! upstream is default to flowrate
        bcdataUp(1)%NodeID = 4
        bcdataUp(1)%TimeArray(1)  = setting%Time%StartTime
        bcdataUp(1)%TimeArray(2)  = setting%Time%EndTime + 100.0 !s
        bcdataUp(1)%ValueArray(1) = initial_flowrate(3)  ! m^3/s
        bcdataUp(1)%ValueArray(2) = initial_flowrate(3)  ! m^3/2

        call case_simple_weir_and_nodes &
            (channel_length, channel_breadth, subdivide_length, lowerZ, upperZ, &
            initial_flowrate, depth_upstream, depth_dnstream, initial_depth,    &
            side_slope, inlet_offset, discharge_coefficient1,                   &
            discharge_coefficient2, full_depth, end_contractions, ManningsN,    &
            roughness_type, idepth_type, linkR, nodeR, linkI, nodeI,linkYN,     &
            nodeYN, linkName, nodeName)

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name

    end subroutine case_simple_weir_initialize
    !
    !==========================================================================
    !
    ! PRIVATE BELOW HERE
    !
    !==========================================================================
    subroutine case_simple_weir_and_nodes &
        (channel_length, channel_breadth, subdivide_length, lowerZ, upperZ,  &
        initial_flowrate, depth_upstream, depth_dnstream, initial_depth,    &
        side_slope, inlet_offset, discharge_coefficient1,                   &
        discharge_coefficient2, full_depth, end_contractions, ManningsN,    &
        roughness_type, idepth_type, linkR, nodeR, linkI, nodeI,linkYN,     &
        nodeYN, linkName, nodeName)
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   creates a weir in between two rectangular channels
        !
        ! Method:
        !    
        !-----------------------------------------------------------------------------
        
        character(64) :: subroutine_name = 'case_simple_weir_and_nodes'

        real(8),  intent(in)  :: channel_length(:), channel_breadth(:), subdivide_length(:)
        real(8),  intent(in)  :: lowerZ(:), upperZ(:), initial_flowrate(:)
        real(8),  intent(in)  :: depth_upstream(:), depth_dnstream(:), initial_depth(:)
        real(8),  intent(in)  :: side_slope(:), inlet_offset(:), end_contractions(:)
        real(8),  intent(in)  :: discharge_coefficient1(:), discharge_coefficient2(:)
        real(8),  intent(in)  :: full_depth(:), ManningsN(:)

        integer, intent(in):: roughness_type, idepth_type(:)

        integer,   dimension(:,:), allocatable, target, intent(out)    :: linkI
        integer,   dimension(:,:), allocatable, target, intent(out)    :: nodeI

        real(8),      dimension(:,:), allocatable, target, intent(out)    :: linkR
        real(8),      dimension(:,:), allocatable, target, intent(out)    :: nodeR

        logical,   dimension(:,:), allocatable, target, intent(out)    :: linkYN
        logical,   dimension(:,:), allocatable, target, intent(out)    :: nodeYN

        type(string), dimension(:), allocatable, target, intent(out)   :: linkName
        type(string), dimension(:), allocatable, target, intent(out)   :: nodeName

        integer :: ii, mm

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        N_link = 3
        N_node = 4

        call allocate_linknode_storage &
            (linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName)

        ! assign the indexes
        linkI(:,li_idx) = (/ (ii, ii=1,N_link) /)
        nodeI(:,ni_idx) = (/ (ii, ii=1,N_node) /)

        ! assign no names for links
        linkName(1)%str = 'Channel Dn'
        linkName(2)%str = 'Weir'
        linkName(3)%str = 'Channel Up'

        ! assign zeros for accumulators
        nodeI(:,ni_N_link_d) = 0
        nodeI(:,ni_N_link_u) = 0

        ! assign uniform physical data
        linkI(:,li_roughness_type)  = roughness_type
        linkR(:,lr_Roughness)       = ManningsN

        ! designate the upstream nodes
        nodeI(1,ni_node_type) = nBCdn
        nodeR(1,nr_Zbottom) = lowerZ(1)
        nodeName(1)%str = 'DownstreamBC'

        ! designate the downstream node
        nodeI(2,ni_node_type) = nJ2
        nodeR(2,nr_Zbottom) = lowerZ(2)
        nodeName(2)%str = 'JunctionDn'

        ! designate the downstream node
        nodeI(3,ni_node_type) = nJ2
        nodeR(3,nr_Zbottom) = lowerZ(3)
        nodeName(3)%str = 'JunctionUp'

        ! designate the downstream node
        nodeI(4,ni_node_type) = nBCup
        nodeR(4,nr_Zbottom) = upperZ(3)
        nodeName(4)%str = 'UpstreamBC'

        ! assign the link types
        linkI(1,li_link_type) = lChannel
        linkI(2,li_link_type) = lWeir
        linkI(3,li_link_type) = lChannel

        ! assign weir type
        linkI(2,li_weir_type) = lVnotchWeir

        ! assign link geometry
        linkI(1,li_geometry) = lRectangular
        linkI(2,li_geometry) = lTriangular
        linkI(3,li_geometry) = lRectangular

        ! assign the link position and mappings
        linkI(1,li_Mnode_d) = 1
        linkI(1,li_Mnode_u) = 2

        linkI(2,li_Mnode_d) = 2
        linkI(2,li_Mnode_u) = 3

        linkI(3,li_Mnode_d) = 3
        linkI(3,li_Mnode_u) = 4

        do mm = 1,N_link
            linkR(mm,lr_Length)                 = channel_length(mm)
            linkR(mm,lr_BreadthScale)           = channel_breadth(mm)
            linkR(mm,lr_ElementLength)          = subdivide_length(mm)
            linkR(mm,lr_InitialFlowrate)        = initial_flowrate(mm)
            linkR(mm,lr_InitialUpstreamDepth)   = depth_upstream(mm)
            linkR(mm,lr_InitialDnstreamDepth)   = depth_dnstream(mm)
            linkR(mm,lr_InitialDepth)           = initial_depth(mm)
            linkR(mm,lr_SideSlope)              = side_slope(mm)
            linkR(mm,lr_LeftSlope)              = side_slope(mm)
            linkR(mm,lr_RightSlope)             = side_slope(mm)
            linkR(mm,lr_InletOffset)            = inlet_offset(mm)
            linkR(mm,lr_OutletOffset)           = inlet_offset(mm)
            linkR(mm,lr_DischargeCoeff1)        = discharge_coefficient1(mm)
            linkR(mm,lr_DischargeCoeff2)        = discharge_coefficient2(mm)
            linkR(mm,lr_FullDepth)              = full_depth(mm)
            linkR(mm,lr_EndContractions)        = end_contractions(mm)
            linkI(mm,li_InitialDepthType)       = idepth_type(mm)
        end do

        if ((debuglevel > 0) .or. (debuglevelall > 0)) then
            print *
            print *, subroutine_name,'-----------------------------------'
            print *, 'link info'
            print *, linkI(:,li_idx), ' idx'
            print *, linkI(:,li_link_type), ' type'
            print *, linkI(:,li_Mnode_u) , ' upstream node'
            print *, linkI(:,li_Mnode_d) , ' downstream node'
            print *, linkR(:,lr_InitialDnstreamDepth), 'downstream depth'
            print *, linkR(:,lr_InitialUpstreamDepth), 'upstream depth'
            print *, linkR(:,lr_InletOffset), ' inletOffset'
            print *, linkR(:,lr_OutletOffset), ' outletOffset'
            print *, 'node info'
            print *, nodeI(:,ni_idx), ' idx'
            print *, nodeI(:,ni_node_type), ' type'
            print *, nodeR(:,nr_Zbottom),   ' zbottom'
            ! print *, nodeI(:,ni_N_link_d), 'number of downstream links'
            ! print *, nodeI(:,ni_Mlink_d1), 'downstream1 link'
            ! print *, nodeI(:,ni_N_link_u), 'number of upstream links'
            ! print *, nodeI(:,ni_Mlink_u1), 'upstream1 link'
        end if

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine case_simple_weir_and_nodes
    !
end module case_simple_weir
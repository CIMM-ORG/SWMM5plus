!
! Custom test case for swashes
!==========================================================================
!
module case_swashes
    !
    use allocate_storage
    use array_index
    use bc
    use data_keys
    use globals
    use setting_definition

    implicit none

    private

    public :: case_shashes_initialize

    integer :: debuglevel = 0

contains
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine case_shashes_initialize &
        (channel_length, channel_breadth, channel_topwidth, subdivide_length, &
        lowerZ, upperZ, initial_flowrate, init_depth, depth_upstream,         &
        depth_dnstream, ManningsN, roughness_type, idepth_type, full_depth,   &
        N_link, N_node, linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, &
        nodeName, bcdataDn, bcdataUp, zbottom)
        !
        ! initialize the link-node system and boundary conditions for a simple channel
        !
        character(64) :: subroutine_name = 'case_shashes_initialize'

        real,  intent(in)  :: channel_length(:), channel_breadth(:)
        real,  intent(in)  :: channel_topwidth(:), subdivide_length(:)
        real,  intent(in)  :: lowerZ(:), upperZ(:),  initial_flowrate(:)
        real,  intent(in)  :: depth_upstream(:), depth_dnstream(:), init_depth(:)
        real,  intent(in)  :: ManningsN(:), full_depth(:), zbottom(:)

        integer, intent(in):: roughness_type, idepth_type(:), N_link, N_node

        integer,   dimension(:,:), allocatable, target, intent(out)    :: linkI
        integer,   dimension(:,:), allocatable, target, intent(out)    :: nodeI

        real,      dimension(:,:), allocatable, target, intent(out)    :: linkR
        real,      dimension(:,:), allocatable, target, intent(out)    :: nodeR

        logical,   dimension(:,:), allocatable, target, intent(out)    :: linkYN
        logical,   dimension(:,:), allocatable, target, intent(out)    :: nodeYN

        type(string), dimension(:), allocatable, target, intent(out)   :: linkName
        type(string), dimension(:), allocatable, target, intent(out)   :: nodeName

        type(bcType), dimension(:), allocatable, intent(out) :: bcdataUp, bcdataDn

        integer    :: ntimepoint, ndnstreamBC, nupstreamBC


        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        ! Boundary conditions
        ntimepoint = 2
        nupstreamBC = 1
        ndnstreamBC = 1

        ! check if

        call bc_allocate &
            (bcdataDn, bcdataUp, ndnstreamBC, nupstreamBC, ntimepoint)

        ! assign values
        ! downstream is default to elevation
        bcdataDn(1)%NodeID = N_node
        bcdataDn(1)%TimeArray(1)     = setting%Time%StartTime
        bcdataDn(1)%TimeArray(2)     = setting%Time%EndTime + 100.0 !s
        bcdataDn(1)%ValueArray(1)    = lowerZ(N_link) +  depth_dnstream(N_link)  ! m
        bcdataDn(1)%ValueArray(2)    = lowerZ(N_link) +  depth_dnstream(N_link)  ! m

        ! upstream is default to flowrate
        bcdataUp(1)%NodeID = 1
        bcdataUp(1)%TimeArray(1)  = setting%Time%StartTime
        bcdataUp(1)%TimeArray(2)  = setting%Time%EndTime + 100.0 !s
        bcdataUp(1)%ValueArray(1) = initial_flowrate(1)  ! m^3/s
        bcdataUp(1)%ValueArray(2) = initial_flowrate(1) ! m^3/2

        call case_swashes_links_and_nodes &
            (channel_length, channel_breadth, channel_topwidth, subdivide_length, &
            lowerZ, upperZ, initial_flowrate, init_depth, depth_upstream,         &
            depth_dnstream, ManningsN, roughness_type, idepth_type, full_depth,   &
            linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName,       &
            N_link, N_node, zbottom)

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine case_shashes_initialize
    !
    !==========================================================================
    !
    ! PRIVATE BELOW HERE
    !
    !==========================================================================
    !
    subroutine case_swashes_links_and_nodes &
        (channel_length, channel_breadth, channel_topwidth, subdivide_length, &
        lowerZ, upperZ, initial_flowrate, init_depth, depth_upstream,         &
        depth_dnstream, ManningsN, roughness_type, idepth_type, full_depth,   &
        linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName,       &
        N_link, N_node, zbottom)
        !
        ! creates swashes test case
        !
        character(64) :: subroutine_name = 'case_swashes_links_and_nodes'

        real,  intent(in)  :: channel_length(:), channel_breadth(:), full_depth(:)
        real,  intent(in)  :: channel_topwidth(:), subdivide_length(:), zbottom(:)
        real,  intent(in)  :: lowerZ(:), upperZ(:), ManningsN(:), initial_flowrate(:)
        real,  intent(in)  :: depth_upstream(:), depth_dnstream(:), init_depth(:)

        integer, intent(in):: roughness_type, idepth_type(:), N_link, N_node

        integer,   dimension(:,:), allocatable, target, intent(out)    :: linkI
        integer,   dimension(:,:), allocatable, target, intent(out)    :: nodeI

        real,      dimension(:,:), allocatable, target, intent(out)    :: linkR
        real,      dimension(:,:), allocatable, target, intent(out)    :: nodeR

        logical,   dimension(:,:), allocatable, target, intent(out)    :: linkYN
        logical,   dimension(:,:), allocatable, target, intent(out)    :: nodeYN

        type(string), dimension(:), allocatable, target, intent(out)   :: linkName
        type(string), dimension(:), allocatable, target, intent(out)   :: nodeName

        integer :: ii, mm

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        call allocate_linknode_storage &
            (linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName)

        ! assign the indexes
        linkI(:,li_idx) = (/ (ii, ii=1,N_link) /)
        nodeI(:,ni_idx) = (/ (ii, ii=1,N_node) /)

        ! assign no names for links
        do ii=1,N_link
            linkName(ii)%str = "Pipe"
        end do

        ! assign zeros for accumulators
        nodeI(:,ni_N_link_d) = 0
        nodeI(:,ni_N_link_u) = 0

        ! assign uniform physical data
        linkI(:,li_roughness_type)  = roughness_type
        linkR(:,lr_Roughness)       = ManningsN

        
        ! print *
        ! print *, 'b)       ii  ,     xvalue  ,    lowerz    ,     upperz  '
        ! do ii = 1, N_link
        !     print *, ii, (ii)*channel_length(ii) - channel_length(ii)/2.0, &
        !     lowerz(ii), upperz(ii)
        ! enddo
        ! stop

        do ii = 2,N_node-1
            nodeI(ii,ni_node_type) = nJ2
            nodeR(ii,nr_Zbottom)   = lowerZ(ii) 
        end do

        ! designate the upstream nodes
        nodeI(1,ni_node_type) = nBCup
        nodeR(1,nr_Zbottom) = upperZ(1)
        nodeName(1)%str = 'UpstreamBC'

        ! designate the downstream node
        nodeI(N_node,ni_node_type) = nBCdn
        nodeR(N_node,nr_Zbottom) = lowerZ(N_link)
        nodeName(N_node)%str = 'DownstreamBC'

        ! print *
        ! print *, 'c)       ii  ,     xvalue  ,    zbottom    '
        ! do ii = 2, N_node-1
        !     print *, ii, (ii-1)*channel_length(ii), nodeR(ii,nr_Zbottom)
        ! enddo
        ! stop

        ! assign the link types
        linkI(:,li_link_type) = lPipe

        ! assign all as rectangular channels
        linkI(:,li_geometry) = lRectangular

        ! assign the link position and mappings
        do ii = 1, N_link
            linkI(ii,li_Mnode_u) = ii ! map to upstream node
            linkI(ii,li_Mnode_d) = ii + 1 ! map to downstream node

            linkR(ii,lr_Length)          = channel_length(ii)
            linkR(ii,lr_BreadthScale)    = channel_breadth(ii)
            linkR(ii,lr_TopWidth)        = channel_topwidth(ii)
            linkR(ii,lr_ElementLength)   = subdivide_length(ii)
            linkR(ii,lr_InitialFlowrate) = initial_flowrate(ii)
            linkR(ii,lr_InitialDepth)    = init_depth(ii)
            linkR(ii,lr_FullDepth)       = full_depth(ii)
            linkR(ii,lr_InitialUpstreamDepth)    = depth_upstream(ii)
            linkR(ii,lr_InitialDnstreamDepth)    = depth_dnstream(ii)
            linkI(ii,li_InitialDepthType)        = idepth_type(ii)

        enddo
        if ((debuglevel > 0) .or. (debuglevelall > 0)) then
            print *
            print *, subroutine_name,'-----------------------------------'
            print *, 'link info'
            print *, linkI(:,li_idx), ' idx'
            print *, linkI(:,li_link_type), ' type'
            print *, linkI(:,li_Mnode_u) , ' upstream node'
            print *, linkI(:,li_Mnode_d) , ' downstream node'
            print *, linkR(:,lr_InitialFlowrate), 'initial flowrate'
            print *, linkR(:,lr_InitialDepth), 'initial depth' 
            print *, ''
            print *, 'node info'
            print *, nodeI(:,ni_idx), ' idx'
            print *, nodeI(:,ni_node_type), ' type'
            print *, nodeI(:,ni_N_link_d), 'number of downstream links'
            !print *, nodeI(:,ni_Mlink_d1), 'downstream1 link'
            print *, nodeI(:,ni_N_link_u), 'number of upstream links'
            !print *, nodeI(:,ni_Mlink_u1), 'upstream1 link'
            print *, nodeR(:,nr_Zbottom), 'node zbottom'
        endif
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine case_swashes_links_and_nodes
    !
    !==========================================================================
    ! END OF MODULE case_simple_channel
    !==========================================================================
end module case_swashes
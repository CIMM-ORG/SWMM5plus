! module case_y_channel
!
! Custom setup for 3 links and 4 nodes arranged in a Y
!
!==========================================================================
!
module case_y_channel
    !
    use allocate_storage
    use array_index
    use bc
    use data_keys
    use globals
    use setting_definition
    use utility

    implicit none

    private

    public :: case_y_channel_initialize

    integer :: debuglevel = 1

contains
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine case_y_channel_initialize &
        (channel_length, channel_breadth, channel_topwidth, subdivide_length, &
        lowerZ, upperZ, initial_flowrate, init_depth, depth_upstream, &
        depth_dnstream, ManningsN, roughness_type, idepth_type,             &
        linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName,     &
        bcdataDn, bcdataUp)

        character(64) :: subroutine_name = 'case_y_channel_initialize'

        real(8),  intent(in)  :: channel_length(:), channel_breadth(:)
        real(8),  intent(in)  :: channel_topwidth(:), subdivide_length(:)
        real(8),  intent(in)  :: lowerZ(:), upperZ(:),  initial_flowrate(:)
        real(8),  intent(in)  :: depth_upstream(:), depth_dnstream(:), init_depth(:)
        real(8),  intent(in)  :: ManningsN(:)

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

        ntimepoint = 2
        N_BCupstream = 2
        N_BCdnstream = 1

        ! check if

        call bc_allocate &
            (bcdataDn, bcdataUp)

        do ii = 1, N_BCdnstream
            allocate( bcdataDn(ii)%TimeArray(ntimepoint), stat=allocation_status, errmsg=emsg)
            call utility_check_allocation (allocation_status, emsg)
            allocate( bcdataDn(ii)%ValueArray(ntimepoint), stat=allocation_status, errmsg=emsg)
            call utility_check_allocation (allocation_status, emsg)
            bcdataDn(ii)%TimeArray      = nullvalueR
            bcdataDn(ii)%ValueArray     = nullvalueR
        enddo
        do ii = 1, N_BCupstream
            allocate( bcdataUp(ii)%TimeArray(ntimepoint), stat=allocation_status, errmsg=emsg)
            call utility_check_allocation (allocation_status, emsg)
            allocate( bcdataUp(ii)%ValueArray(ntimepoint), stat=allocation_status, errmsg=emsg)
            call utility_check_allocation (allocation_status, emsg)
            bcdataUp(ii)%TimeArray      = nullvalueR
            bcdataUp(ii)%ValueArray     = nullvalueR
        enddo

        ! assign values
        ! downstream is default to elevation
        bcdataDn(1)%NodeID = 1
        bcdataDn(1)%TimeArray(1)     = setting%Time%StartTime
        bcdataDn(1)%TimeArray(2)     = setting%Time%EndTime + 100.0 !s
        bcdataDn(1)%ValueArray(1)    = lowerZ(1) +  depth_dnstream(1)  ! m
        bcdataDn(1)%ValueArray(2)    = lowerZ(1) +  depth_dnstream(1) ! m

        ! upstream is default to flowrate
        bcdataUp(1)%NodeID = 3
        bcdataUp(1)%TimeArray(1)  = setting%Time%StartTime
        bcdataUp(1)%TimeArray(2)  = setting%Time%EndTime + 100.0 !s
        bcdataUp(1)%ValueArray(1) = initial_flowrate(1)  ! m^3/s
        bcdataUp(1)%ValueArray(2) = initial_flowrate(1)  ! m^3/2

        bcdataUp(2)%NodeID = 4
        bcdataUp(2)%TimeArray(1)  = setting%Time%StartTime
        bcdataUp(2)%TimeArray(2)  = setting%Time%EndTime + 100.0 !s
        bcdataUp(2)%ValueArray(1) = initial_flowrate(2)  ! m^3/s
        bcdataUp(2)%ValueArray(2) = initial_flowrate(2)  ! m^3/2

        call case_y_channel_links_and_nodes &
            (channel_length, channel_breadth, channel_topwidth, subdivide_length, &
            lowerZ, upperZ, initial_flowrate, init_depth, depth_upstream, &
            depth_dnstream, ManningsN, roughness_type, idepth_type, &
            linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName)

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine case_y_channel_initialize
    !
    !==========================================================================
    !
    ! PRIVATE BELOW HERE
    !
    !==========================================================================
    !
    subroutine case_y_channel_links_and_nodes &
        (channel_length, channel_breadth, channel_topwidth, subdivide_length, &
        lowerZ, upperZ, initial_flowrate, init_depth, depth_upstream, &
        depth_dnstream, ManningsN, roughness_type, idepth_type, &
        linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName)
        !
        ! creates a simple rectangular channel with 1 link and 2 nodes
        !
        character(64) :: subroutine_name = 'case_simple_channel_links_and_nodes'

        real(8),  intent(in)  :: channel_length(:), channel_breadth(:)
        real(8),  intent(in)  :: channel_topwidth(:), subdivide_length(:)
        real(8),  intent(in)  :: lowerZ(:), upperZ(:), ManningsN(:), initial_flowrate(:)
        real(8),  intent(in)  :: depth_upstream(:), depth_dnstream(:), init_depth(:)

        integer, intent(in):: roughness_type, idepth_type(:)

        integer,   dimension(:,:), allocatable, target, intent(out)    :: linkI
        integer,   dimension(:,:), allocatable, target, intent(out)    :: nodeI

        real(8),      dimension(:,:), allocatable, target, intent(out)    :: linkR
        real(8),      dimension(:,:), allocatable, target, intent(out)    :: nodeR

        logical,   dimension(:,:), allocatable, target, intent(out)    :: linkYN
        logical,   dimension(:,:), allocatable, target, intent(out)    :: nodeYN

        type(string), dimension(:), allocatable, target, intent(out)   :: linkName
        type(string), dimension(:), allocatable, target, intent(out)   :: nodeName

        integer :: mm, ii


        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        N_link = 3
        N_node = 4

        call allocate_linknode_storage &
            (linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName)

        ! assign the indexes
        linkI(:,li_idx) = (/ (ii, ii=1,N_link) /)
        nodeI(:,ni_idx) = (/ (ii, ii=1,N_node) /)

        ! assign names for links
        linkName(1)%str = 'Channel_Downstream'
        linkName(2)%str = 'Channel_Upstream_Right'
        linkName(3)%str = 'Channel_Upstream_Left'

        ! assign zeros for accumulators
        nodeI(:,ni_N_link_d) = 0
        nodeI(:,ni_N_link_u) = 0

        ! assign uniform physical data
        linkI(:,li_roughness_type)  = roughness_type
        linkR(:,lr_Roughness)       = ManningsN

        ! designate the downstream node
        nodeI(1,ni_node_type) = nBCdn

        nodeR(1,nr_Zbottom) = lowerZ(1)

        nodeName(1)%str = 'DownstreamBC'

        ! designate the junction node
        nodeI(2,ni_node_type) = nJm

        nodeR(2,nr_Zbottom) = upperZ(1)

        nodeName(2)%str = 'Junction'

        ! designated the upstream nodes
        nodeI(3,ni_node_type) = nBCup
        nodeI(4,ni_node_type) = nBCup

        nodeR(3,nr_Zbottom) = upperZ(2)
        nodeR(4,nr_Zbottom) = upperZ(3)

        nodeName(3)%str = 'UpstreamBC_Right'
        nodeName(4)%str = 'UpstreamBC_Left'

        ! assign the link types
        linkI(:,li_link_type) = lChannel

        ! assign all as rectangular channels
        linkI(:,li_geometry) = lRectangular

        ! assign the link position and mappings

        linkI(1,li_Mnode_d) = 1
        linkI(1,li_Mnode_u) = 2

        linkI(2,li_Mnode_d) = 2
        linkI(2,li_Mnode_u) = 3

        linkI(3,li_Mnode_d) = 2
        linkI(3,li_Mnode_u) = 4

        do mm=1,N_link
            linkR(mm,lr_Length)          = channel_length(mm)
            linkR(mm,lr_BreadthScale)    = channel_breadth(mm)
            linkR(mm,lr_TopWidth)        = channel_topwidth(mm)
            linkR(mm,lr_ElementLength)   = subdivide_length(mm)
            linkR(mm,lr_InitialFlowrate) = initial_flowrate(mm)
            linkI(mm,li_InitialDepthType)= idepth_type(mm)
        enddo

        linkR(:  ,lr_InitialDepth)         = init_depth(:)
        linkR(:  ,lr_InitialDnstreamDepth) = depth_dnstream(:)
        linkR(:  ,lr_InitialUpstreamDepth) = depth_upstream(:)


        !print *, initial_flowrate
        !print *,trim(subroutine_name)
        !stop

        if ((debuglevel > 0) .or. (debuglevelall > 0)) then
            print *
            print *, subroutine_name,'-----------------------------------'
            print *, 'link info'
            print *, linkI(:,li_idx), ' idx'
            print *, linkI(:,li_link_type), ' type'
            print *, linkI(:,li_Mnode_u) , ' upstream node'
            print *, linkI(:,li_Mnode_d) , ' downstream node'
            print *, ''
            print *, 'node info'
            print *, nodeI(:,ni_idx), ' idx'
            print *, nodeI(:,ni_node_type), ' type'
            !print *, nodeI(:,ni_N_link_d), 'number of downstream links'
            !print *, nodeI(:,ni_Mlink_d1), 'downstream1 link'
            !print *, nodeI(:,ni_N_link_u), 'number of upstream links'
            !print *, nodeI(:,ni_Mlink_u1), 'upstream1 link'
        endif

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine case_y_channel_links_and_nodes
    !

    !==========================================================================
    ! END OF MODULE stub
    !==========================================================================
end module case_y_channel

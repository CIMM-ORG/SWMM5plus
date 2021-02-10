!==========================================================================
!
module custom_network

    use array_index
    use data_keys
    use globals

    implicit none
    private

    public :: custom_1link_network
    public :: custom_6link_1_line_network
    public :: custom_3link_Y_network
    public :: custom_6link_Y_network

    integer, private :: debuglevel = 0

contains
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine custom_1link_network &
        (linkR, nodeR, linkI, nodeI, linkName, nodeName)
        !
        ! one link that is a simple rectangular open channel
        ! Basic setup
        ! Overall length = 10,000 m
        ! Upstream Zbottom = 11 m
        ! Downstream Zbottom = 1 m
        ! effective slope is 0.001
        ! Channel breadth is 3 m
        ! Mannings n = 0.03
        ! Target depth = 0.5 m
        ! flow rate = (1/0.03) * 1.5 * (1.5 / 4)^(2/3) * (0.001)^(1/2) = 0.822225 m^3/s
        !
        character(64) :: subroutine_name = 'custom_1link_network'

        integer, intent(in out) :: linkI(:,:)
        integer, intent(in out) :: nodeI(:,:)

        real, intent(in out)    :: linkR(:,:)
        real, intent(in out)    :: nodeR(:,:)

        type(string), dimension(:), intent(in out)   :: nodeName, linkName

        integer :: ii
        integer :: thisnode

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        ! ERROR CHECK --------------------------
        if (N_link /= 1) then
            print *, 'error this network requires 1 links in ',subroutine_name
            STOP
        endif
        if (N_node /= 2) then
            print *, 'error this network requires 2 nodes in ',subroutine_name
            STOP
        endif

        ! SET UP CONNECTIVITY AND PHYSICS ------------------

        ! assign the indexes
        linkI(:,li_idx) = (/ (ii, ii=1,N_link) /)
        nodeI(:,ni_idx) = (/ (ii, ii=1,N_node) /)

        ! assign no names for links
        do ii=1,N_link
            linkName(ii)%str = ""
        end do

        ! assign zeros for accumulators
        nodeI(:,ni_N_link_d) = 0
        nodeI(:,ni_N_link_u) = 0

        ! assign uniform physical data
        linkI(:,li_roughness_type) = lManningsN

        ! designate the upstream nodes
        nodeI(1,ni_node_type) = nBCup

        nodeR(1,nr_Zbottom) = 10.0

        nodeName(1)%str = 'UpstreamBC'

        ! designate the downstream node
        nodeI(2,ni_node_type) = nBCdn

        nodeR(2,nr_Zbottom) = 1.0

        nodeName(2)%str = 'DownstreamBC'

        ! assign the link types
        linkI(:,li_link_type) = lChannel

        ! assign all as rectangular channels
        linkI(:,li_geometry) = lRectangular

        ! assign the link position and mappings

        ! link 1
        linkI(1,li_Mnode_u) = 1 ! map to upstream node
        linkI(1,li_Mnode_d) = 2 ! map to downstream node

        linkR(1,lr_Length) = 10000.0
        linkR(1,lr_Breadth) = 3.0

        ! cycle through links to assign nodes
        do ii=1,N_link
            ! get the upstream node of the link (if it exists)
            thisnode = linkI(ii,li_Mnode_u)
            ! look for next available downstream link position
            if (thisnode > 0) then
                if (nodeI(thisnode,ni_Mlink_d1) == nullvalueI) then
                    nodeI(thisnode,ni_Mlink_d1) = ii
                else
                    print *, 'node ',thisnode
                    print *, 'error: attempt to assign 2 downstream links to 1 node in ',subroutine_name
                    stop
                endif
                ! increment the downstream link counter
                nodeI(thisnode,ni_N_link_d) = nodeI(thisnode,ni_N_link_d)+1
            endif

            ! get the downstream node of the link (if it exists)
            thisnode = linkI(ii,li_Mnode_d)
            ! look for the next available upstream link position
            if (thisnode > 0) then
                if (nodeI(thisnode,ni_Mlink_u1) == nullvalueI) then
                    nodeI(thisnode,ni_Mlink_u1) = ii
                elseif (nodeI(thisnode,ni_Mlink_u2) == nullvalueI) then
                    nodeI(thisnode,ni_Mlink_u2) = ii
                else
                    print *, 'node ',thisnode
                    print *, 'error: attempt to assign 3 upstream links to 1 node in ',subroutine_name
                    stop
                endif
                ! increment the upstream link counter
                nodeI(thisnode,ni_N_link_u) = nodeI(thisnode,ni_N_link_u)+1
            endif
        enddo

        print *
        print *, 'link info'
        print *, linkI(:,li_idx), ' idx'
        print *, linkI(:,li_link_type), ' type'
        print *, linkI(:,li_Mnode_u) , ' upstream node'
        print *, linkI(:,li_Mnode_d) , ' downstream node'
        print *,
        print *, 'node info'
        print *, nodeI(:,ni_idx), ' idx'
        print *, nodeI(:,ni_node_type), ' type'
        print *, nodeI(:,ni_N_link_d), 'number of downstream links'
        print *, nodeI(:,ni_Mlink_d1), 'downstream1 link'
        print *, nodeI(:,ni_N_link_u), 'number of upstream links'
        print *, nodeI(:,ni_Mlink_u1), 'upstream1 link'

        !stop

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine custom_1link_network
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine custom_6link_1_line_network &
        (linkR, nodeR, linkI, nodeI, linkName, nodeName)

        character(64) :: subroutine_name = 'custom_6link_1_line_network'

        integer, intent(in out) :: linkI(:,:)
        integer, intent(in out) :: nodeI(:,:)

        real, intent(in out)    :: linkR(:,:)
        real, intent(in out)    :: nodeR(:,:)

        type(string), dimension(:), intent(in out)   :: linkName, nodeName

        integer :: ii
        integer :: thisnode

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        ! ERROR CHECK --------------------------
        if (N_link /= 6) then
            print *, 'error this network requires 6 links in ',subroutine_name
            STOP
        endif
        if (N_node /= 7) then
            print *, 'error this network requires 7 nodes in ',subroutine_name
            STOP
        endif

        ! SET UP CONNECTIVITY AND PHYSICS ------------------

        ! assign the indexes
        linkI(:,li_idx) = (/ (ii, ii=1,N_link) /)
        nodeI(:,ni_idx) = (/ (ii, ii=1,N_node) /)

        ! assign zeros for accumulators
        nodeI(:,ni_N_link_d) = 0
        nodeI(:,ni_N_link_u) = 0

        ! assign uniform physical data
        linkI(:,li_roughness_type) = lManningsN

        ! assign no names for links
        do ii=1,N_link
            linkName(ii)%str = ""
        end do

        !print *, linkI(:,li_idx)
        !print *, nodeI(:,ni_idx)

        ! designate the upstream nodes
        nodeI(1,ni_node_type) = nBCup

        nodeR(1,nr_Zbottom) = 10.0

        nodeName(1)%str = 'UpstreamBC'

        ! designate intermediate nodes
        nodeI(2:6,ni_node_type) = nJ2

        nodeR(2,nr_Zbottom) = 10.0
        nodeR(3,nr_Zbottom) = 9.0
        nodeR(4,nr_Zbottom) = 7.0
        nodeR(5,nr_Zbottom) = 6.0
        nodeR(6,nr_Zbottom) = 4.0

        nodeName(2)%str = 'inter01'
        nodeName(3)%str = 'inter02'
        nodeName(4)%str = 'inter03'
        nodeName(5)%str = 'inter04'
        nodeName(6)%str = 'inter05'

        ! designate the downstream node
        nodeI(7,ni_node_type) = nBCdn

        nodeR(7,nr_Zbottom) = 1.5

        nodeName(7)%str = 'DownstreamBC'

        ! assign the link types
        linkI(:,li_link_type) = lChannel

        ! assign all as rectangular channels
        linkI(:,li_geometry) = lRectangular

        ! assign the link position and mappings

        ! link 1
        linkI(1,li_Mnode_u) = 1 ! map to upstream node
        linkI(1,li_Mnode_d) = 2 ! map to downstream node

        linkR(1,lr_Length) = 155.0
        linkR(1,lr_Breadth) = 2.5

        ! link 2
        linkI(2,li_Mnode_u) = 2 ! map to upstream node
        linkI(2,li_Mnode_d) = 3 ! map to downstream node

        linkR(2,lr_Length) = 247.0
        linkR(2,lr_Breadth) = 3.2

        ! link 3
        linkI(3,li_Mnode_u) = 3 ! map to upstream node
        linkI(3,li_Mnode_d) = 4 ! map to downstream node

        linkR(3,lr_Length) = 323.0
        linkR(3,lr_Breadth) = 4.1

        ! link 4
        linkI(4,li_Mnode_u) = 4 ! map to upstream node
        linkI(4,li_Mnode_d) = 5 ! map to downstream node

        linkR(4,lr_Length) = 280.0
        linkR(4,lr_Breadth) = 3.9

        ! link 5
        linkI(5,li_Mnode_u) = 5 ! map to upstream node
        linkI(5,li_Mnode_d) = 6 ! map to downstream node

        linkR(5,lr_Length) = 256.0
        linkR(5,lr_Breadth) = 4.5

        ! link 6
        linkI(6,li_Mnode_u) = 6 ! map to upstream node
        linkI(6,li_Mnode_d) = 7 ! map to downstream node

        linkR(6,lr_Length) = 123.0
        linkR(6,lr_Breadth) = 4.7


        ! cycle through links to assign nodes
        do ii=1,N_link
            ! get the upstream node of the link (if it exists)
            thisnode = linkI(ii,li_Mnode_u)
            ! look for next available downstream link position
            if (thisnode > 0) then
                if (nodeI(thisnode,ni_Mlink_d1) == nullvalueI) then
                    nodeI(thisnode,ni_Mlink_d1) = ii
                else
                    print *, 'node ',thisnode
                    print *, 'error: attempt to assign 2 downstream links to 1 node in ',subroutine_name
                    stop
                endif
                ! increment the downstream link counter
                nodeI(thisnode,ni_N_link_d) = nodeI(thisnode,ni_N_link_d)+1
            endif

            ! get the downstream node of the link (if it exists)
            thisnode = linkI(ii,li_Mnode_d)
            ! look for the next available upstream link position
            if (thisnode > 0) then
                if (nodeI(thisnode,ni_Mlink_u1) == nullvalueI) then
                    nodeI(thisnode,ni_Mlink_u1) = ii
                elseif (nodeI(thisnode,ni_Mlink_u2) == nullvalueI) then
                    nodeI(thisnode,ni_Mlink_u2) = ii
                else
                    print *, 'node ',thisnode
                    print *, 'error: attempt to assign 3 upstream links to 1 node in ',subroutine_name
                    stop
                endif
                ! increment the upstream link counter
                nodeI(thisnode,ni_N_link_u) = nodeI(thisnode,ni_N_link_u)+1
            endif
        enddo

        print *
        print *, 'link info'
        print *, linkI(:,li_idx), ' idx'
        print *, linkI(:,li_link_type), ' type'
        print *, linkI(:,li_Mnode_u) , ' upstream node'
        print *, linkI(:,li_Mnode_d) , ' downstream node'
        print *,
        print *, 'node info'
        print *, nodeI(:,ni_idx), ' idx'
        print *, nodeI(:,ni_node_type), ' type'
        print *, nodeI(:,ni_N_link_d), 'number of downstream links'
        print *, nodeI(:,ni_Mlink_d1), 'downstream1 link'
        !print *, nodeI(:,ni_Mlink_d2), 'downstream2 link'
        !print *, nodeI(:,ni_Mlink_d3), 'downstream2 link'
        print *, nodeI(:,ni_N_link_u), 'number of upstream links'
        print *, nodeI(:,ni_Mlink_u1), 'upstream1 link'
        !print *, nodeI(:,ni_Mlink_u2), 'upstream2 link'
        !print *, nodeI(:,ni_Mlink_u3), 'upstream2 link'

        !print *, nodeI(2,8)
        !stop

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine custom_6link_1_line_network
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine custom_3link_Y_network &
        (linkR, nodeR, linkI, nodeI, linkName, nodeName)

        character(64) :: subroutine_name = 'custom_3link_Y_network'

        integer, intent(in out) :: linkI(:,:)
        integer, intent(in out) :: nodeI(:,:)

        real, intent(in out)    :: linkR(:,:)
        real, intent(in out)    :: nodeR(:,:)

        type(string), dimension(:), intent(in out)   :: linkName, nodeName

        integer :: ii
        integer :: thisnode

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        ! ERROR CHECK --------------------------
        if (N_link /= 3) then
            print *, 'error this network requires 3 links in ',subroutine_name
            STOP
        endif
        if (N_node /= 4) then
            print *, 'error this network requires 4 nodes in ',subroutine_name
            STOP
        endif

        ! SET UP CONNECTIVITY AND PHYSICS ------------------

        ! assign the indexes
        linkI(:,li_idx) = (/ (ii, ii=1,N_link) /)
        nodeI(:,ni_idx) = (/ (ii, ii=1,N_node) /)

        ! assign zeros for accumulators
        nodeI(:,ni_N_link_d) = 0
        nodeI(:,ni_N_link_u) = 0

        ! assign uniform physical data
        linkI(:,li_roughness_type) = lManningsN

        ! assign no names for links
        do ii=1,N_link
            linkName(ii)%str = ""
        end do

        !print *, linkI(:,li_idx)
        !print *, nodeI(:,ni_idx)

        ! designate the upstream nodes
        nodeI(1,ni_node_type) = nBCup
        nodeI(3,ni_node_type) = nBCup

        nodeR(1,nr_Zbottom) = 10.0
        nodeR(3,nr_Zbottom) = 10.0

        nodeName(1)%str = 'UpstreamBC01'
        nodeName(3)%str = 'UpstreamBC02'

        ! designate the junction node for multiple links
        nodeI(2,ni_node_type) = nJm

        nodeR(2,nr_Zbottom) = 5.0

        nodeName(2)%str = 'Junction'

        ! designate the downstream node
        nodeI(4,ni_node_type) = nBCdn

        nodeR(4,nr_Zbottom) = 0

        nodeName(4)%str = 'DownstreamBC'

        ! assign the link types
        linkI(:,li_link_type) = lChannel

        ! assign all as rectangular channels
        linkI(:,li_geometry) = lRectangular

        ! assign the link position and mappings

        ! link 1 - upstream Y branch 1
        linkI(1,li_Mnode_u) = 1 ! map to upstream node
        linkI(1,li_Mnode_d) = 2 ! map to downstream node

        linkR(1,lr_Length) = 10000.0
        linkR(1,lr_Breadth) = 2.5

        ! link 2 - upstream Y branch 2
        linkI(2,li_Mnode_u) = 3 ! map to upstream node
        linkI(2,li_Mnode_d) = 2 ! map to downstream node

        linkR(2,lr_Length) = 10000.0
        linkR(2,lr_Breadth) = 3.2

        ! link 3 - downstram branch
        linkI(3,li_Mnode_u) = 2 ! map to upstream node
        linkI(3,li_Mnode_d) = 4 ! map to downstream node

        linkR(3,lr_Length) = 10000.0
        linkR(3,lr_Breadth) = 4.1

        ! cycle through links to assign nodes
        do ii=1,N_link
            ! get the upstream node of the link (if it exists)
            thisnode = linkI(ii,li_Mnode_u)
            ! look for next available downstream link position
            if (thisnode > 0) then
                if (nodeI(thisnode,ni_Mlink_d1) == nullvalueI) then
                    nodeI(thisnode,ni_Mlink_d1) = ii
                else
                    print *, 'node ',thisnode
                    print *, 'ni_Mlink_d1 ', ni_Mlink_d1
                    print *, 'nodeI(thisnode,ni_Mlink_d1) ', nodeI(thisnode,ni_Mlink_d1)
                    print *, 'error: attempt to assign 2 downstream links to 1 node in ',subroutine_name
                    stop
                endif
                ! increment the downstream link counter
                nodeI(thisnode,ni_N_link_d) = nodeI(thisnode,ni_N_link_d)+1
            endif

            ! get the downstream node of the link (if it exists)
            thisnode = linkI(ii,li_Mnode_d)
            ! look for the next available upstream link position
            if (thisnode > 0) then
                if (nodeI(thisnode,ni_Mlink_u1) == nullvalueI) then
                    nodeI(thisnode,ni_Mlink_u1) = ii
                elseif (nodeI(thisnode,ni_Mlink_u2) == nullvalueI) then
                    nodeI(thisnode,ni_Mlink_u2) = ii
                else
                    print *, 'node ',thisnode
                    print *, 'error: attempt to assign 3 upstream links to 1 node in ',subroutine_name
                    stop
                endif
                ! increment the upstream link counter
                nodeI(thisnode,ni_N_link_u) = nodeI(thisnode,ni_N_link_u)+1
            endif
        enddo

        print *
        print *, 'link info'
        print *, linkI(:,li_idx), ' idx'
        print *, linkI(:,li_link_type), ' type'
        print *, linkI(:,li_Mnode_u) , ' upstream node'
        print *, linkI(:,li_Mnode_d) , ' downstream node'
        print *,
        print *, 'node info'
        print *, nodeI(:,ni_idx), ' idx'
        print *, nodeI(:,ni_node_type), ' type'
        print *, nodeI(:,ni_N_link_d), 'number of downstream links'
        print *, nodeI(:,ni_Mlink_d1), 'downstream1 link'
        !print *, nodeI(:,ni_Mlink_d2), 'downstream2 link'
        !print *, nodeI(:,ni_Mlink_d3), 'downstream2 link'
        print *, nodeI(:,ni_N_link_u), 'number of upstream links'
        print *, nodeI(:,ni_Mlink_u1), 'upstream1 link'
        print *, nodeI(:,ni_Mlink_u2), 'upstream2 link'
        !print *, nodeI(:,ni_Mlink_u3), 'upstream2 link'

        !print *, nodeI(2,8)
        !stop

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine custom_3link_Y_network
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine custom_6link_Y_network &
        (linkR, nodeR, linkI, nodeI, linkName, nodeName)

        character(64) :: subroutine_name = 'custom_6link_Y_network'

        integer, intent(in out) :: linkI(:,:)
        integer, intent(in out) :: nodeI(:,:)

        real, intent(in out)    :: linkR(:,:)
        real, intent(in out)    :: nodeR(:,:)

        type(string), dimension(:), intent(in out)   :: linkName, nodeName

        integer :: ii
        integer :: thisnode

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        ! ERROR CHECK --------------------------
        if (N_link /= 6) then
            print *, 'error this network requires 6 links in ',subroutine_name
            STOP
        endif
        if (N_node /= 7) then
            print *, 'error this network requires 7 nodes in ',subroutine_name
            STOP
        endif

        ! SET UP CONNECTIVITY AND PHYSICS ------------------

        ! assign the indexes
        linkI(:,li_idx) = (/ (ii, ii=1,N_link) /)
        nodeI(:,ni_idx) = (/ (ii, ii=1,N_node) /)


        ! assign zeros for accumulators
        nodeI(:,ni_N_link_d) = 0
        nodeI(:,ni_N_link_u) = 0

        ! assign uniform physical data
        linkI(:,li_roughness_type) = lManningsN


        ! assign no names for links
        do ii=1,N_link
            linkName(ii)%str = ""
        end do

        !print *, linkI(:,li_idx)
        !print *, nodeI(:,ni_idx)

        ! designate the upstream nodes
        nodeI(1,ni_node_type) = nBCup
        nodeI(4,ni_node_type) = nBCup

        nodeR(1,nr_Zbottom) = 10.0
        nodeR(4,nr_Zbottom) = 8.0

        nodeName(1)%str = 'UpstreamBC01'
        nodeName(4)%str = 'UpstreamBC02'

        ! designated the 2-element junctions
        nodeI(2,ni_node_type) = nJ2
        nodeI(5,ni_node_type) = nJ2
        nodeI(6,ni_node_type) = nJ2

        nodeR(2,nr_Zbottom) = 9.0
        nodeR(5,nr_Zbottom) = 7.0
        nodeR(6,nr_Zbottom) = 3.0

        nodeName(2)%str = 'inter01'
        nodeName(5)%str = 'inter02'
        nodeName(6)%str = 'inter03'

        ! designate the junction node for multiple links
        nodeI(3,ni_node_type) = nJm

        nodeR(3,nr_Zbottom) = 4.0

        ! designate the downstream node
        nodeI(7,ni_node_type) = nBCdn

        nodeR(7,nr_Zbottom) = 1.5

        nodeName(7)%str = 'DownstreamBC'

        ! assign the link types
        linkI(:,li_link_type) = lChannel

        ! assign all as rectangular channels
        linkI(:,li_geometry) = lRectangular

        ! assign the link position and mappings

        ! link 1 - upstream Y branch 1a
        linkI(1,li_Mnode_u) = 1 ! map to upstream node
        linkI(1,li_Mnode_d) = 2 ! map to downstream node

        linkR(1,lr_Length) = 155.0
        linkR(1,lr_Breadth) = 2.5

        ! link 2 - upstream Y branch 1b
        linkI(2,li_Mnode_u) = 2 ! map to upstream node
        linkI(2,li_Mnode_d) = 3 ! map to downstream node (Y)

        linkR(2,lr_Length) = 247.0
        linkR(2,lr_Breadth) = 3.2

        ! link 3 - upstream Y branch 2a
        linkI(3,li_Mnode_u) = 4 ! map to upstream node
        linkI(3,li_Mnode_d) = 5 ! map to downstream node

        linkR(3,lr_Length) = 192.0
        linkR(3,lr_Breadth) = 2.9

        ! link 4 - upstream Y branch 2b
        linkI(4,li_Mnode_u) = 5 ! map to upstream node
        linkI(4,li_Mnode_d) = 3 ! map to downstream node

        linkR(4,lr_Length) = 239.0
        linkR(4,lr_Breadth) = 3.7

        ! link 5 - downstram branch
        linkI(5,li_Mnode_u) = 3 ! map to upstream node
        linkI(5,li_Mnode_d) = 6 ! map to downstream node

        linkR(5,lr_Length) = 323.0
        linkR(5,lr_Breadth) = 4.1

        linkI(6,li_Mnode_u) = 6 ! map to upstream node
        linkI(6,li_Mnode_d) = 7 ! map to downstream node

        linkR(6,lr_Length) = 389.0
        linkR(6,lr_Breadth) = 4.6

        ! cycle through links to assign nodes
        do ii=1,N_link
            ! get the upstream node of the link (if it exists)
            thisnode = linkI(ii,li_Mnode_u)
            ! look for next available downstream link position
            if (thisnode > 0) then
                if (nodeI(thisnode,ni_Mlink_d1) == nullvalueI) then
                    nodeI(thisnode,ni_Mlink_d1) = ii
                else
                    print *, 'node ',thisnode
                    print *, 'error: attempt to assign 2 downstream links to 1 node in ',subroutine_name
                    stop
                endif
                ! increment the downstream link counter
                nodeI(thisnode,ni_N_link_d) = nodeI(thisnode,ni_N_link_d)+1
            endif

            ! get the downstream node of the link (if it exists)
            thisnode = linkI(ii,li_Mnode_d)
            ! look for the next available upstream link position
            if (thisnode > 0) then
                if (nodeI(thisnode,ni_Mlink_u1) == nullvalueI) then
                    nodeI(thisnode,ni_Mlink_u1) = ii
                elseif (nodeI(thisnode,ni_Mlink_u2) == nullvalueI) then
                    nodeI(thisnode,ni_Mlink_u2) = ii
                else
                    print *, 'node ',thisnode
                    print *, 'error: attempt to assign 3 upstream links to 1 node in ',subroutine_name
                    stop
                endif
                ! increment the upstream link counter
                nodeI(thisnode,ni_N_link_u) = nodeI(thisnode,ni_N_link_u)+1
            endif
        enddo

        print *
        print *, 'link info'
        print *, linkI(:,li_idx), ' idx'
        print *, linkI(:,li_link_type), ' type'
        print *, linkI(:,li_Mnode_u) , ' upstream node'
        print *, linkI(:,li_Mnode_d) , ' downstream node'
        print *,
        print *, 'node info'
        print *, nodeI(:,ni_idx), ' idx'
        print *, nodeI(:,ni_node_type), ' type'
        print *, nodeI(:,ni_N_link_d), 'number of downstream links'
        print *, nodeI(:,ni_Mlink_d1), 'downstream1 link'
        !print *, nodeI(:,ni_Mlink_d2), 'downstream2 link'
        !print *, nodeI(:,ni_Mlink_d3), 'downstream2 link'
        print *, nodeI(:,ni_N_link_u), 'number of upstream links'
        print *, nodeI(:,ni_Mlink_u1), 'upstream1 link'
        print *, nodeI(:,ni_Mlink_u2), 'upstream2 link'
        !print *, nodeI(:,ni_Mlink_u3), 'upstream2 link'

        !print *, nodeI(2,8)
        !stop

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine custom_6link_Y_network
    !
    !==========================================================================
    !==========================================================================
end module custom_network

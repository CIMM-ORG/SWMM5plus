!
!% module network_define
!
!% Handles relationship between coarse link-node network and high-resolution
!% element-face network. This module defines all the indexes and mappings
!
!% Sazzad Sharior 06/01/2021
!
!==========================================================================
!
module network_define
    !
    use interface
    use utility_allocate
    use discretization
    use define_indexes
    use define_keys
    use define_globals
    use define_settings

    implicit none

    private

    public :: init_network_define_toplevel

contains
    !
    !==========================================================================
    ! PUBLIC
    !==========================================================================
    !
    subroutine init_network_define_toplevel ()
    !--------------------------------------------------------------------------
    !
    !% Initializes a element-face network from a link-node network.
    !%   Requires network links and nodes before execution
    !
    !--------------------------------------------------------------------------

        integer :: ii, jj

        character(64) :: subroutine_name = 'init_network_define_toplevel'

    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name

        !% get the slope of each link given the node Z values
        call init_network_linkslope ()

        !% divide the link node networks in elements and faces
        call init_network_datacreate ()

        !% replaces ni_elemface_idx of nJ2 nodes for the upstream elem
        !% of the face associated with the node
        call init_network_update_nj2_elem ()

        sync all

        !% print result
        if (setting%Debug%File%network_define) then
            !% only using the first processor to print results
            if (this_image() == 1) then
                do ii = 1,num_images()
                    print*
                    print*, '===================================================================' //&
                    '==================================================================='
                    print*, 'image = ', ii
                    print*, '.......................Elements...............................'
                    print*
                    print*, '     ei_Lidx     ei_Gidx     ei_link    ei_node    Mface_uL    Mface_dL'
                    do jj = 1,N_elem(ii)
                        print*, elemI(jj,ei_Lidx)[ii], elemI(jj,ei_Gidx)[ii], elemI(jj,ei_link_Gidx_SWMM)[ii], &
                        elemI(jj,ei_node_Gidx_SWMM)[ii], elemI(jj,ei_Mface_uL)[ii], elemI(jj,ei_Mface_dL)[ii]
                    enddo
                    print*
                    print*, '.......................Faces.............................'
                    print*, 'a)  fi_Lidx     fi_Gidx     elem_uL     elem_dL     C_image' //&
                    '    GElem_up    GElem_dn     node_id     link_id'
                    do jj = 1,N_face(ii)
                        print*, faceI(jj,fi_Lidx)[ii],faceI(jj,fi_Gidx)[ii],faceI(jj,fi_Melem_uL)[ii], &
                        faceI(jj,fi_Melem_dL)[ii],faceI(jj,fi_Connected_image)[ii], faceI(jj,fi_GhostElem_uL)[ii],&
                        faceI(jj,fi_GhostElem_dL)[ii],faceI(jj,fi_node_idx)[ii],faceI(jj, fi_link_idx)[ii]
                    enddo

                    ! print*
                    ! print*, '.......................Faces..................................'
                    ! print*, 'b)     fi_Lidx     fi_BCtype   fYN_isInteriorFace    fYN_isSharedFace'//&
                    ! '    fYN_isnull    fYN_isUpGhost    fYN_isDnGhost'
                    ! do jj = 1,N_face(ii)
                    !     print*, faceI(jj,fi_Lidx)[ii],' ',faceI(jj,fi_BCtype)[ii],&
                    !     '           ',faceYN(jj,fYN_isInteriorFace)[ii], &
                    !     '                    ',faceYN(jj,fYN_isSharedFace)[ii], &
                    !     '              ',faceYN(jj,fYN_isnull)[ii], &
                    !     '            ',faceYN(jj,fYN_isUpGhost)[ii], &
                    !     '              ',faceYN(jj,fYN_isDnGhost)[ii]
                    ! enddo
                    print*, '===================================================================' //&
                    '==================================================================='
                    print*
                    call execute_command_line('')
                enddo
            endif
        endif
        
        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name
    end subroutine init_network_define_toplevel
    !
    !==========================================================================
    ! PRIVATE
    !==========================================================================
    !
    subroutine init_network_update_nj2_elem()
        integer, allocatable :: nJ2_nodes(:)
        integer              :: N_nJ2_nodes
        character(64) :: subroutine_name = 'init_network_update_nj2_elem'

        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name

        N_nJ2_nodes = count((node%I(:, ni_node_type) == nJ2) .and. (node%I(:, ni_P_image) == this_image()))

        if (N_nJ2_nodes > 0) then

            nJ2_nodes = pack(node%I(:, ni_idx), (node%I(:, ni_node_type) == nJ2) &
                        .and. (node%I(:, ni_P_image) == this_image()))

            node%I(nJ2_nodes, ni_elemface_idx) = faceI(node%I(nJ2_nodes, ni_elemface_idx), fi_Melem_uL)

        endif

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name
    end subroutine init_network_update_nj2_elem
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_linkslope()
    !--------------------------------------------------------------------------
    !
    !% compute the slope across each link
    !
    !--------------------------------------------------------------------------

        character(64) :: subroutine_name = 'init_network_linkslope'

        integer, pointer :: NodeUp, NodeDn
        real(8), pointer :: zUp, zDn, Slope, Length
        integer          :: mm

    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name

        do mm = 1,N_link
            !% Inputs
            NodeUp      => link%I(mm,li_Mnode_u)
            NodeDn      => link%I(mm,li_Mnode_d)
            zUp         => node%R(NodeUp,nr_Zbottom)
            zDn         => node%R(NodeDn,nr_Zbottom)
            Length      => link%R(mm,lr_Length)
            !%-------------------------------------------------------------------------
            !% HACK: Original length is used for slope calculation instead of adjusted
            !% length. Using adjusted lenghts will induce errors in slope calculations
            !% and will distort the original network.
            !%-------------------------------------------------------------------------
            !% Output
            Slope       => link%R(mm,lr_Slope)

            Slope = (zUp - zDn) / Length
        end do

        if (setting%Debug%File%network_define) then
            !% provide output for debugging
            print *, subroutine_name,'--------------------------------'
            print *, 'link ID,               Slope,             length'
            do mm=1,N_link
                print *, mm, link%R(mm,lr_Slope), link%R(mm,lr_Length)
            end do
        endif

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name

    end subroutine init_network_linkslope
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_datacreate()
    !--------------------------------------------------------------------------
    !
    !% creates the network of elements and faces from nodes and link
    !
    !--------------------------------------------------------------------------

        integer :: ii, image
        integer :: ElemGlobalIdx, FaceGlobalIdx
        integer :: ElemLocalIdx, FacelocalIdx
        integer :: unique_faces
        integer, pointer :: N_Images, Lidx
        integer, pointer :: NodeUp, NodeDn, NodeUpTyp, NodeDnTyp

        character(64) :: subroutine_name = 'init_network_datacreate'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name

        !% initializing global element and face id
        ElemGlobalIdx = first_elem_index
        FaceGlobalIdx = first_face_index - oneI

        !% initializing local element and face id
        ElemLocalIdx = first_elem_index
        FacelocalIdx = first_face_index - oneI

        !% Setting the local image value
        image = this_image()

        !% initialize the global indexes of elements and faces
        call init_network_set_global_indexes &
            (image, ElemGlobalIdx, FaceGlobalIdx)

        !% set the dummy element
        call init_network_set_dummy_elem ()

        !% handle all the links and nodes in a partition
        call init_network_handle_partition &
            (image, ElemLocalIdx, FacelocalIdx, ElemGlobalIdx, FaceGlobalIdx)

        !% finish mapping all the junction branch and faces that were not
        !% handeled in handle_link_nodes subroutine
        call init_network_map_nodes (image)

        !% set interior face logical
        call init_network_set_interior_faceYN ()

        !% shared faces are mapped by copying data from different images
        !% thus a sync all is needed
        sync all

        !% set the same global face idx for shared faces across images
        call init_network_map_shared_faces (image)

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name
    end subroutine init_network_datacreate
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_set_global_indexes &
        (image, ElemGlobalCounter, FaceGlobalCounter)
    !
    !--------------------------------------------------------------------------
    !
    !% This subroutine initaializes the global element and face indexes.
    !% N_elem and N_face are global vectors of whose row number corresponds
    !% to the number of elements and faces in a certain image. By looping
    !% through all the images the element and face global indexes are set
    !
    !--------------------------------------------------------------------------
    !
        integer, intent(in)     :: image
        integer, intent(inout)  :: ElemGlobalCounter, FaceGlobalCounter

        integer                 :: ii

        character(64) :: subroutine_name = 'init_network_set_global_indexes'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name

        if(image /= 1) then
           do ii=1, image-1
              ElemGlobalCounter = ElemGlobalCounter + N_elem(ii)
              FaceGlobalCounter = FaceGlobalCounter + N_unique_face(ii) + oneI
           end do
        end if

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name
    end subroutine init_network_set_global_indexes
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_set_dummy_elem ()
    !
    !--------------------------------------------------------------------------
    !
    !% initialize dummy row in the elem arrays
    !
    !--------------------------------------------------------------------------
    !
        integer       :: dummyIdx

        character(64) :: subroutine_name = 'init_network_set_dummy_elem'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name

        !% index for the dummy element
        dummyIdx = max_caf_elem_N + N_dummy_elem

        !% set the elem arrays with dummy values
        elemI(dummyIdx, ei_Lidx)        = dummyIdx
        elemI(dummyIdx,ei_elementType)  = dummy

        elemYN(dummyIdx,eYN_isDummy)    = .true.

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name
    end subroutine init_network_set_dummy_elem
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_handle_partition &
        (image, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, FaceGlobalCounter)
    !
    !--------------------------------------------------------------------------
    !
    !% Traverse through all the links and nodes in a partition and creates
    !% elements and faces. This subroutine assumes there will be at least
    !% one link in a partition.
    !
    !--------------------------------------------------------------------------
    !
        integer, intent(in)     :: image
        integer, intent(inout)  :: ElemLocalCounter, FaceLocalCounter
        integer, intent(inout)  :: ElemGlobalCounter, FaceGlobalCounter

        integer                 :: ii, pLink
        integer, pointer        :: thisLink, upNode, dnNode
        integer, dimension(:), allocatable, target :: packed_link_idx

        character(64) :: subroutine_name = 'init_network_handle_partition'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name

        !% pack all the link indexes in a partition
        packed_link_idx = pack(link%I(:,li_idx), (link%I(:,li_P_image) == image))

        !% find the number of links in a partition
        pLink = size(packed_link_idx)

        !% cycling through all the links in a partition
        do ii = 1,pLink
            !% necessary pointers to the link and its connected nodes
            thisLink => packed_link_idx(ii)
            upNode   => link%I(thisLink,li_Mnode_u)
            dnNode   => link%I(thisLink,li_Mnode_d)

            !% handle the upstream node of the link to create elements and faces
            call init_network_handle_upstreamnode &
                (image, thisLink, upNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                FaceGlobalCounter)

            !% handle the link to create elements and faces
            call init_network_handle_link &
                (image, thisLink, upNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                FaceGlobalCounter)

            !% handle the downstream node of the link to create elements and faces
            call init_network_handle_downstreamnode &
                (image, thisLink, dnNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                FaceGlobalCounter)
        end do

        !% deallocate temporary array
        deallocate(packed_link_idx)

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name
    end subroutine init_network_handle_partition
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_map_nodes (image)
    !
    !--------------------------------------------------------------------------
    !
    !% map all the interior junction faces in an image
    !
    !--------------------------------------------------------------------------
    !
        integer, intent(in)    :: image

        integer :: ii, pNodes
        integer, pointer :: thisJunctionNode, nodeType
        integer, dimension(:), allocatable, target :: packed_node_idx, JunctionElementIdx

        character(64) :: subroutine_name = 'init_network_map_nodes'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name


        !% pack all the interior node indexes in a partition to find face maps
        packed_node_idx = pack(node%I(:,ni_idx),                       &
                              ((node%I(:,ni_P_image)   == image) .and. &
                              ((node%I(:,ni_node_type) == nJm ) .or.   &
                               (node%I(:,ni_node_type) == nJ2 ) ) )    )

        !% number of interior nodes in a partition
        pNodes = size(packed_node_idx)


        !% cycle through all the nJm nodes and set the face maps
        do ii = 1,pNodes
            thisJunctionNode => packed_node_idx(ii)
            nodeType         => node%I(thisJunctionNode,ni_node_type)

            select case (nodeType)

                case (nJm)
                    JunctionElementIdx = pack( elemI(:,ei_Lidx), &
                                             ( elemI(:,ei_node_Gidx_SWMM) == thisJunctionNode) )

                    call init_network_map_nJm_branches (image, thisJunctionNode, JunctionElementIdx)

                    !% deallocate temporary array
                    deallocate(JunctionElementIdx)

                case (nJ2)
                    call init_network_map_nJ2 (image, thisJunctionNode)

            end select
        end do

        !% deallocate temporary array
        deallocate(packed_node_idx)

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name
    end subroutine init_network_map_nodes
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_map_shared_faces (image)
    !
    !--------------------------------------------------------------------------
    !
    !% set the global indexes for shared faces across images
    !
    !--------------------------------------------------------------------------
    !
        integer, intent(in)    :: image

        integer :: ii, NsharedFaces
        integer, pointer :: fLidx, nIdx, eUp, eDn, targetImage, nodeType
        logical, pointer :: isUpGhost, isDnGhost
        integer, dimension(:), allocatable, target ::  sharedFaces

        character(64) :: subroutine_name = 'init_network_map_shared_faces'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name

        !% pack the shared faces in an image
        sharedFaces = pack(faceI(:,fi_Lidx), faceYN(:,fYN_isSharedFace))

        !% fid the size of the pack
        NsharedFaces = size(sharedFaces)

        do ii = 1,NsharedFaces
            fLidx       => sharedFaces(ii)
            nIdx        => faceI(fLidx,fi_node_idx)
            nodeType    => node%I(nIdx,ni_node_type)

            select case (nodeType)

                case (nJ2)
                    call init_network_map_shared_nJ2_nodes (image, fLidx, nIdx)

                case (nJm)
                    call init_network_map_shared_nJm_nodes (image, fLidx, nIdx)

            end select
        enddo

        !% deallocate temporary array
        deallocate(sharedFaces)

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name
    end subroutine init_network_map_shared_faces
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_handle_upstreamnode &
        (image, thisLink, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
        FaceGlobalCounter)
    !
    !--------------------------------------------------------------------------
    !
    !% Handle the node upstream of a link to create elements and faces
    !
    !--------------------------------------------------------------------------
    !
        integer, intent(in)     :: image, thisLink, thisNode
        integer, intent(inout)  :: ElemLocalCounter, FaceLocalCounter
        integer, intent(inout)  :: ElemGlobalCounter, FaceGlobalCounter

        integer                 :: ii
        integer, pointer        :: nAssignStatus, nodeType, linkUp

        character(64) :: subroutine_name = 'init_network_handle_upstreamnode'
    !--------------------------------------------------------------------------
        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name

        !% Check 1: If the node is in the partition
        if (node%I(thisNode,ni_P_image) == image) then

            !% necessary pointers
            nAssignStatus => node%I(thisNode,ni_assigned)
            nodeType      => node%I(thisNode,ni_node_type)

            select case (nodeType)

                !% Handle upstream boundary nodes
                case(nBCup)
                    !% Check 2: If the node has already been assigned
                    if (nAssignStatus == nUnassigned) then

                        !% Advance the local and global counter for the upstream
                        !% boundary node
                        FaceLocalCounter  = FaceLocalCounter + oneI
                        FaceGlobalCounter = FaceGlobalCounter + oneI

                        !% integer data
                        faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                        faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                        faceI(FacelocalCounter,fi_node_idx) = thisNode
                        !% an upstream boundary face does not have any local upstream element
                        !% thus, it is mapped to the dummy element
                        faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem
                        faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter
                        faceI(FaceLocalCounter,fi_BCtype)   = BCup
                        node%I(thisNode,ni_elemface_idx)    = FaceLocalCounter
                        !% set zbottom
                        faceR(FaceLocalCounter,fr_Zbottom)  = node%R(thisNode,nr_Zbottom)
                        !% change the node assignmebt value
                        nAssignStatus =  nAssigned
                    endif

                !% Handle 2 branch junction nodes
                case (nJ2)
                    !% Check 2: If the node has already been assigned
                    if (nAssignStatus == nUnassigned) then

                        !% Advance the local and global counter for the upstream
                        !% nJ2 node
                        FaceLocalCounter  = FaceLocalCounter + oneI
                        FaceGlobalCounter = FaceGlobalCounter + oneI

                        !% integer data
                        faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                        faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                        faceI(FacelocalCounter,fi_node_idx) = thisNode
                        faceI(FacelocalCounter,fi_Melem_uL) = ElemLocalCounter - oneI
                        faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter
                        faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
                        !% set zbottom
                        faceR(FaceLocalCounter,fr_Zbottom)  = node%R(thisNode,nr_Zbottom)
                        !% First assign the face index to the nJ2 node, then will
                        !% update with the elem_uL index of the upstream element
                        node%I(thisNode,ni_elemface_idx)     = FaceLocalCounter

                        !% Check 3: If the node is an edge node (meaning this node is the
                        !% connecting node across partitions)
                        if (node%I(thisNode,ni_P_is_boundary) == EdgeNode) then

                            !% An upstream edge node indicates there are no local
                            !% elements upstream of that node. Thus it is mapped to
                            !% the dummy element
                            faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem

                            !% logical data
                            faceYN(FaceLocalCounter,fYN_isSharedFace) = .true.
                            faceYN(FaceLocalCounter,fYN_isUpGhost)    = .true.

                            !% find the connecting image to this face
                            linkUp  => node%I(thisNode,ni_Mlink_u1)

                            faceI(FaceLocalCounter,fi_Connected_image) = link%I(linkUp,li_P_image)

                            if (image < faceI(FaceLocalCounter,fi_Connected_image)) then
                                !% we only set the global indexes where the connection
                                !% is in higher order than the current image.
                                !% (for example if current image = 1 and connection is 2,
                                !% we set the global counter. But when the current image = 2 but
                                !% the connection is 1, we set it from init_network_map_shared_faces
                                !% subroutine)
                                faceI(FaceLocalCounter,fi_Gidx) = FaceGlobalCounter
                            else
                                !% set global index as nullvalue for shared faces.
                                !% these global indexes will be set later
                                faceI(FaceLocalCounter,fi_Gidx)     = nullvalueI
                            endif

                        endif

                        !% change the node assignmebt value
                        nAssignStatus =  nAssigned
                    endif

                !% Handle junction nodes with more than 2 branches (multi branch junction node).
                case (nJm)

                    !% Check 2: If the node has already been assigned
                    if (nAssignStatus == nUnassigned) then

                        !% multibranch junction nodes will have both elements and faces.
                        !% thus, a seperate subroutine is required to handle these nodes
                        call init_network_handle_nJm &
                            (image, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                            FaceGlobalCounter, nAssignStatus)
                    endif

                !% Handle storage nodes
                case (nStorage)
                    print*
                    print*, 'In ', subroutine_name
                    print*, 'error: storage node is not handeled yet'

                case default
                    print*
                    print*, 'In ', subroutine_name
                    print*, 'error: node ' // node%Names(thisNode)%str // &
                            ' has an unexpected nodeType', nodeType
                    stop
            end select

        !% handle the node if it is not in the partition
        else

            !% Advance the local and global counter for the upstream
            !% node outside the partition
            FaceLocalCounter  = FaceLocalCounter + oneI
            FaceGlobalCounter = FaceGlobalCounter + oneI

            !% integer data
            faceI(FacelocalCounter,fi_Lidx)     = FacelocalCounter
            faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
            faceI(FacelocalCounter,fi_node_idx) = thisNode
            faceI(FacelocalCounter,fi_link_idx) = thisLink
            faceI(FacelocalCounter,fi_Connected_image) = node%I(thisNode,ni_P_image)

            !% if the upstream node is not in the partiton,
            !% the face map to upstream is mapped to
            !% the dummy element
            faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem
            !% since no upstream node indicates start of a partiton,
            !% the downstream element will be initialized elem idx
            faceI(FacelocalCounter,fi_Melem_dL) = ElemLocalCounter

            !% since this is a shared face, it will have a copy in other image and they will
            !% both share same global index. so, the face immediately after this shared face
            !% will have the global index set from the init_network_set_global_indexes subroutine.
            !% However, since the init_network_handle_link subroutine will advance the global face
            !% count anyway, the count  here is needed to be adjusted by substracting one from the
            !% count.
            FaceGlobalCounter = FaceGlobalCounter - oneI

            !% logical data
            faceYN(FacelocalCounter,fYN_isSharedFace)   = .true.
            faceYN(FacelocalCounter,fYN_isUpGhost)      = .true.

            if (image < faceI(FaceLocalCounter,fi_Connected_image)) then
                !% we only set the global indexes where the connection
                !% is in higher order than the current image.
                !% (for example if current image = 1 and connection is 2,
                !% we set the global counter. But when the current image = 2 but
                !% the connection is 1, we set it from init_network_map_shared_faces
                !% subroutine)
                faceI(FaceLocalCounter,fi_Gidx) = FaceGlobalCounter
            else
                !% set global index as nullvalue for shared faces.
                !% these global indexes will be set later
                faceI(FaceLocalCounter,fi_Gidx)     = nullvalueI
            endif
        endif

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name
    end subroutine init_network_handle_upstreamnode
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_handle_link &
        (image, thisLink, upNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
        FaceGlobalCounter)
    !--------------------------------------------------------------------------
    !
    !% handle links in a partition
    !
    !--------------------------------------------------------------------------

        integer, intent(in)     :: image, thisLink, upNode
        integer, intent(inout)  :: ElemLocalCounter, FaceLocalCounter
        integer, intent(inout)  :: ElemGlobalCounter, FaceGlobalCounter

        integer                 :: ii
        real(8)                 :: zCenter, zDownstream
        integer, pointer        :: lAssignStatus, NlinkElem
        real(8), pointer        :: zUpstream

        character(64) :: subroutine_name = 'init_network_handle_link'

    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name

        !% necessary pointers
        lAssignStatus => link%I(thisLink,li_assigned)

        if (lAssignStatus == lUnAssigned) then
            NlinkElem     => link%I(thisLink,li_N_element)
            zUpstream     => node%R(upNode,nr_Zbottom)

            !% store the ID of the first (upstream) element in this link
            link%I(thisLink,li_first_elem_idx)   = ElemLocalCounter
            !% reference elevations at cell center
            zCenter     = zUpstream - 0.5 * link%R(thisLink,lr_ElementLength) * link%R(thisLink,lr_Slope)
            zDownstream = zUpstream - link%R(thisLink,lr_ElementLength) * link%R(thisLink,lr_Slope)
            do ii = 1, NlinkElem
                !%................................................................
                !% Element arrays update
                !%................................................................

                !% integer data
                elemI(ElemLocalCounter,ei_Lidx)             = ElemLocalCounter
                elemI(ElemLocalCounter,ei_Gidx)             = ElemGlobalCounter
                elemI(ElemLocalCounter,ei_elementType)      = link%I(thisLink,li_link_type)
                elemI(ElemLocalCounter,ei_link_Gidx_SWMM)   = thisLink
                elemI(ElemLocalCounter,ei_Mface_uL)         = FaceLocalCounter
                elemI(ElemLocalCounter,ei_Mface_dL)         = FaceLocalCounter + oneI
                elemI(ElemLocalCounter,ei_link_pos)         = ii

                !% real data
                elemR(ElemLocalCounter,er_Length)           = link%R(thisLink,lr_ElementLength)
                elemR(ElemLocalCounter,er_Zbottom)          = zCenter

                !%................................................................
                !% Face arrays update
                !%................................................................

                if (ii < NlinkElem) then
                !% advance only the downstream interior face counter of a link element
                    FaceLocalCounter  = FaceLocalCounter  + oneI
                    FaceGlobalCounter = FaceGlobalCounter + oneI

                    !% face integer data
                    faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                    faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                    faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter + oneI
                    faceI(FaceLocalCounter,fi_Melem_uL) = ElemLocalCounter
                    faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
                    faceI(FaceLocalCounter,fi_link_idx) = thisLink
                    faceR(FaceLocalCounter,fr_Zbottom)  = zDownstream
                endif

                !% counter for element z bottom calculation
                zCenter     = zCenter - link%R(thisLink,lr_ElementLength) * link%R(thisLink,lr_Slope)
                zDownstream = zDownstream - link%R(thisLink,lr_ElementLength) * link%R(thisLink,lr_Slope)

                !% Advance the element counter
                ElemLocalCounter  = ElemLocalCounter  + oneI
                ElemGlobalCounter = ElemGlobalCounter + oneI
            end do

            lAssignStatus = lAssigned
            link%I(thisLink,li_last_elem_idx)    = ElemLocalCounter - oneI

            !% re initialize the face local and global counters for the next expected node
            ! FaceLocalCounter  = FaceLocalCounter  - oneI
            ! FaceGlobalCounter = FaceGlobalCounter - oneI
        endif

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name
    end subroutine init_network_handle_link
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_handle_downstreamnode &
        (image, thisLink, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
        FaceGlobalCounter)
    !--------------------------------------------------------------------------
    !
    !% handle the node downstream of a link
    !
    !--------------------------------------------------------------------------
    !
        integer, intent(in)    :: image, thisLink, thisNode
        integer, intent(inout) :: ElemLocalCounter, FaceLocalCounter
        integer, intent(inout) :: ElemGlobalCounter, FaceGlobalCounter

        integer, pointer :: nAssignStatus, nodeType, linkDn

        character(64) :: subroutine_name = 'init_network_handle_downstreamnode'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name

        !% Check 1: Is the node is in the partition
        if (node%I(thisNode,ni_P_image) == image) then

            !% necessary pointers
            nAssignStatus => node%I(thisNode,ni_assigned)
            nodeType      => node%I(thisNode,ni_node_type)

            select case (nodeType)

                case(nBCdn)
                    !% Check 2: If the node has already been assigned
                    if (nAssignStatus == nUnassigned) then

                        !% Advance face local and global counters for BCdn node
                        FaceLocalCounter  = FaceLocalCounter  + oneI
                        FaceGlobalCounter = FaceGlobalCounter + oneI

                        !% integer data
                        !% an downstream boundary face does not have any local downstream element
                        !% thus, it is mapped to the dummy element
                        faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                        faceI(FacelocalCounter,fi_Gidx)     = FaceGlobalCounter
                        faceI(FacelocalCounter,fi_Melem_uL) = ElemLocalCounter - oneI
                        faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem
                        faceI(FaceLocalCounter,fi_BCtype)   = BCdn
                        faceI(FacelocalCounter,fi_node_idx) = thisNode
                        !% set zbottom
                        faceR(FaceLocalCounter,fr_Zbottom)  = node%R(thisNode,nr_Zbottom)
                        node%I(thisNode,ni_elemface_idx)    = FaceLocalCounter

                        !% change the node assignmebt value
                        nAssignStatus =  nAssigned
                    endif

                case (nJ2)
                    !% Check 2: If the node has already been assigned
                    if (nAssignStatus == nUnassigned) then

                        !% Advance face local and global counters for nJ2 node
                        FaceLocalCounter  = FaceLocalCounter  + oneI
                        FaceGlobalCounter = FaceGlobalCounter + oneI

                        !% integer data
                        faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                        faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                        faceI(FacelocalCounter,fi_node_idx) = thisNode
                        faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
                        faceI(FacelocalCounter,fi_Melem_uL) = ElemLocalCounter - oneI
                        faceI(FacelocalCounter,fi_Melem_dL) = ElemLocalCounter
                        !% set zbottom
                        faceR(FaceLocalCounter,fr_Zbottom)  = node%R(thisNode,nr_Zbottom)
                        !% First assign the face index to the nJ2 node, then will
                        !% update with the elem_uL index of the upstream element
                        node%I(thisNode,ni_elemface_idx)   = FaceLocalCounter

                        !% integer data
                        if (node%I(thisNode,ni_P_is_boundary) == EdgeNode) then

                            !% A downstream edge node indicates there are no local
                            !% elements downstream of that node
                            faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem

                            !% logical data
                            faceYN(FaceLocalCounter,fYN_isSharedFace) = .true.
                            faceYN(FaceLocalCounter,fYN_isDnGhost)    = .true.

                            !% find the connecting image to this face
                            linkDn  => node%I(thisNode,ni_Mlink_d1)

                            faceI(FaceLocalCounter,fi_Connected_image)    = link%I(linkDn,li_P_image)

                            if (image < faceI(FaceLocalCounter,fi_Connected_image)) then
                                !% we only set the global indexes where the connection
                                !% is in higher order than the current image.
                                !% (for example if current image = 1 and connection is 2,
                                !% we set the global counter. But when the current image = 2 but
                                !% the connection is 1, we set it from init_network_map_shared_faces
                                !% subroutine)
                                faceI(FaceLocalCounter,fi_Gidx) = FaceGlobalCounter
                            else
                                !% set global index as nullvalue for shared faces.
                                !% these global indexes will be set later
                                faceI(FaceLocalCounter,fi_Gidx)     = nullvalueI
                            endif
                        endif

                        !% change the node assignmebt value
                        nAssignStatus =  nAssigned
                    endif

                case (nJm)
                    !% Check 2: If the node has already been assigned
                    if (nAssignStatus == nUnassigned) then

                        call init_network_handle_nJm &
                            (image, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                            FaceGlobalCounter, nAssignStatus)

                    endif

                case (nStorage)
                    print*
                    print*, 'In ', subroutine_name
                    print*, 'error: storage node is not handeled yet'

                case default

                    print*
                    print*, 'In ', subroutine_name
                    print*, 'error: node ' // node%Names(thisNode)%str // &
                            ' has an unexpected nodeType', nodeType
                    stop
            end select
        else
            !% Advance face local and global counters for nodes outside of the partition
            FaceLocalCounter  = FaceLocalCounter  + oneI
            FaceGlobalCounter = FaceGlobalCounter + oneI
            !% if the downstream node is not in the partiton.
            !% through subdivide_link_going_downstream subroutine
            !% upstream map to the element has alrady been set.
            !% However, downstream map has set to wrong value.
            !% Thus, setting the map elem ds to dummy elem
            !% integer data
            faceI(FacelocalCounter,fi_Lidx)     = FaceLocalCounter
            faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem
            faceI(FacelocalCounter,fi_Melem_uL) = ElemLocalCounter - oneI
            faceI(FacelocalCounter,fi_node_idx) = thisNode
            faceI(FacelocalCounter,fi_link_idx) = thisLink
            faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
            faceI(FacelocalCounter,fi_Connected_image) = node%I(thisNode,ni_P_image)

            !% logical data
            faceYN(FacelocalCounter,fYN_isSharedFace) = .true.
            faceYN(FaceLocalCounter,fYN_isDnGhost)    = .true.

            if (image < faceI(FaceLocalCounter,fi_Connected_image)) then
                !% we only set the global indexes where the connection
                !% is in higher order than the current image.
                !% (for example if current image = 1 and connection is 2,
                !% we set the global counter. But when the current image = 2 but
                !% the connection is 1, we set it from init_network_map_shared_faces
                !% subroutine)
                faceI(FaceLocalCounter,fi_Gidx) = FaceGlobalCounter
            else
                !% set global index as nullvalue for shared faces.
                !% these global indexes will be set later
                faceI(FaceLocalCounter,fi_Gidx)     = nullvalueI
            endif
        endif

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name
    end subroutine init_network_handle_downstreamnode
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_handle_nJm &
        (image, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
        FaceGlobalCounter, nAssignStatus)
    !--------------------------------------------------------------------------
    !
    !% subdivides the multi branch junctions into elements and faces
    !
    !--------------------------------------------------------------------------

        integer, intent(in)    :: image, thisNode
        integer, intent(inout) :: ElemLocalCounter, FaceLocalCounter
        integer, intent(inout) :: ElemGlobalCounter, FaceGlobalCounter
        integer, intent(inout) :: nAssignStatus

        integer, pointer :: upBranchIdx, dnBranchIdx

        integer :: ii, upBranchSelector, dnBranchSelector

        character(64) :: subroutine_name = 'init_network_handle_nJm'

    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name

        !%................................................................
        !% Junction main
        !%................................................................

        !% Element Arrays
        !% integer data
        elemI(ElemLocalCounter,ei_Lidx)             = ElemLocalCounter
        elemI(ElemLocalCounter,ei_Gidx)             = ElemGlobalCounter
        elemI(ElemLocalCounter,ei_elementType)      = JM
        elemI(ElemLocalCounter,ei_node_Gidx_SWMM)   = thisNode
        !% Assign junction main element to node
        node%I(thisNode,ni_elemface_idx)             = ElemLocalCounter

        !% real data
        elemR(ElemLocalCounter,er_Zbottom) = node%R(thisNode,nr_zbottom)

        !% Advance the element counter to 1st upstream branch
        ElemLocalCounter  = ElemLocalCounter  + oneI
        ElemGlobalCounter = ElemGlobalCounter + oneI

        !%................................................................
        !% Handle Junction Branches
        !%................................................................

        !% initialize selecteros for upstream and downstream branches
        upBranchSelector = zeroI
        dnBranchSelector = zeroI

        !% loopthrough all the branches
        do ii = 1,max_branch_per_node

            !% common junction branch data
            !% element arrays
            !% integer data
            elemI(ElemLocalCounter,ei_Lidx)           = ElemLocalCounter
            elemI(ElemLocalCounter,ei_Gidx)           = ElemGlobalCounter
            elemI(ElemLocalCounter,ei_elementType)    = JB
            elemI(ElemLocalCounter,ei_node_Gidx_SWMM) = thisNode

            !% real data
            elemR(ElemLocalCounter,er_Zbottom) = node%R(thisNode,nr_zbottom)

            !% advance the face counters for the branch
            FaceLocalCounter  = FaceLocalCounter  + oneI
            FaceGlobalCounter = FaceGlobalCounter + oneI

            !%......................................................
            !% Upstream Branches
            !%......................................................
            ! if ((ii == 1) .or. (ii == 3) .or. (ii == 5)) then
            select case (mod(ii,2))
            case (1)
            !% finds odd number branches
            !% all the odd numbers are upstream branches
                upBranchSelector = upBranchSelector + oneI
                !% pointer to upstream branch
                upBranchIdx => node%I(thisNode,ni_idx_base1 + upBranchSelector)

                !% elem array
                !% map the upstream face of the branch element
                elemI(ElemLocalCounter,ei_Mface_uL) = FaceLocalCounter

                !% face array
                !% integer data
                faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                faceI(FacelocalCounter,fi_Gidx)     = FaceGlobalCounter
                faceI(FacelocalCounter,fi_node_idx) = thisNode
                faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter
                faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist

                !% real branches
                if (upBranchIdx /= nullvalueI) then
                    !% integer data
                    elemSI(ElemLocalCounter,eSI_JunctionBranch_Exists)           = oneI
                    elemSI(ElemLocalCounter,eSI_JunctionBranch_Link_Connection)  = upBranchIdx
                    elemR(ElemLocalCounter,er_Length) = init_network_nJm_branch_length(upBranchIdx)
                    faceI(FaceLocalCounter,fi_link_idx) = upBranchIdx
                    !% set zbottom
                    faceR(FaceLocalCounter,fr_Zbottom)  = node%R(thisNode,nr_Zbottom)
                    !% Check 4: this node is the connecting node across partitions
                    if ( (node%I(thisNode,ni_P_is_boundary) == EdgeNode)  .and. &
                         (link%I(upBranchIdx,li_P_image)    /= image   ) )  then

                        !% since this is a shared face in the upstream direction,
                        !% the up_map is set to dummy element
                        faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem
                        faceI(FaceLocalCounter,fi_Connected_image) = link%I(upBranchIdx,li_P_image)

                        !% logical data
                        faceYN(FacelocalCounter,fYN_isSharedFace) = .true.
                        faceYN(FaceLocalCounter,fYN_isUpGhost)    = .true.

                        if (image < faceI(FaceLocalCounter,fi_Connected_image)) then
                            !% we only set the global indexes where the connection
                            !% is in higher order than the current image.
                            !% (for example if current image = 1 and connection is 2,
                            !% we set the global counter. But when the current image = 2 but
                            !% the connection is 1, we set it from init_network_map_shared_faces
                            !% subroutine)
                            faceI(FaceLocalCounter,fi_Gidx) = FaceGlobalCounter
                        else
                            !% set global index as nullvalue for shared faces.
                            !% these global indexes will be set later
                            faceI(FaceLocalCounter,fi_Gidx)     = nullvalueI
                        endif
                    endif
                else

                    !% since this is a null face in the upstream direction,
                    !% the up_map is set to dummy element
                    faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem

                    call init_network_nullify_nJm_branch &
                        (ElemLocalCounter, FaceLocalCounter)
                endif
            !%......................................................
            !% Downstream Branches
            !%......................................................
            !else
            case (0)
            !% even number branches

                !% all the even numbers are downstream branches
                dnBranchSelector = dnBranchSelector + oneI
                !% pointer to upstream branch
                dnBranchIdx => node%I(thisNode,ni_idx_base2 + dnBranchSelector)

                !% elem array
                !% integer data
                elemI(ElemLocalCounter,ei_Mface_dL) = FaceLocalCounter

                !% face array
                !% integer data
                faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                faceI(FacelocalCounter,fi_Gidx)     = FaceGlobalCounter
                faceI(FacelocalCounter,fi_node_idx) = thisNode
                faceI(FaceLocalCounter,fi_Melem_uL) = ElemLocalCounter
                faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist

                !% Check 3: if the branch is a valid branch
                if (dnBranchIdx /= nullvalueI) then
                    !% integer data
                    elemSI(ElemLocalCounter,eSI_JunctionBranch_Exists)          = oneI
                    elemSI(ElemLocalCounter,eSI_JunctionBranch_Link_Connection) = dnBranchIdx
                    elemR(ElemLocalCounter,er_Length) = init_network_nJm_branch_length(dnBranchIdx)
                    faceI(FacelocalCounter,fi_link_idx) = dnBranchIdx
                    !% set zbottom
                    faceR(FaceLocalCounter,fr_Zbottom)  = node%R(thisNode,nr_Zbottom)
                    !% Check 4: if the link connecting this branch is a part of this partition and
                    !% the node is not an edge node (meaning this node is the connecting node
                    !% across partitions)
                    if ( (node%I(thisNode,ni_P_is_boundary) == EdgeNode)  .and. &
                         (link%I(dnBranchIdx,li_P_image)    /= image   ) )  then

                        !% since this is a shared face in the downstream direction,
                        !% the dn_map is set to dummy element
                        faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem
                        faceI(FaceLocalCounter,fi_Connected_image) = link%I(dnBranchIdx,li_P_image)

                        !% logical data
                        faceYN(FacelocalCounter,fYN_isSharedFace) = .true.
                        faceYN(FaceLocalCounter,fYN_isDnGhost)    = .true.

                        if (image < faceI(FaceLocalCounter,fi_Connected_image)) then
                            !% we only set the global indexes where the connection
                            !% is in higher order than the current image.
                            !% (for example if current image = 1 and connection is 2,
                            !% we set the global counter. But when the current image = 2 but
                            !% the connection is 1, we set it from init_network_map_shared_faces
                            !% subroutine)
                            faceI(FaceLocalCounter,fi_Gidx) = FaceGlobalCounter
                        else
                            !% set global index as nullvalue for shared faces.
                            !% these global indexes will be set later
                            faceI(FaceLocalCounter,fi_Gidx)     = nullvalueI
                        endif
                    endif

                else

                    !% since this is a null face in the downstream direction,
                    !% the dn_map is set to dummy element
                    faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem

                    call init_network_nullify_nJm_branch &
                        (ElemLocalCounter, FaceLocalCounter)
                endif
            end select

            !% Advance the element counter for next branch
            ElemLocalCounter  = ElemLocalCounter  + oneI
            ElemGlobalCounter = ElemGlobalCounter + oneI
        end do

        !% set status to assigned
        nAssignStatus = nAssigned

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name
    end subroutine init_network_handle_nJm
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_map_nJm_branches (image, thisJNode, JelemIdx)
    !
    !--------------------------------------------------------------------------
    !
    !% map all the multi branch junction elements
    !
    !--------------------------------------------------------------------------
    !
        integer, intent(in)                       :: image, thisJNode
        integer, dimension(:), target, intent(in) :: JelemIdx

        integer          :: ii, upBranchSelector, dnBranchSelector
        integer          :: LinkFirstElem, LinkLastElem
        integer, pointer :: upBranchIdx, dnBranchIdx
        integer, pointer :: eIdx, fLidx

        character(64) :: subroutine_name = 'init_network_map_nJm_branches'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name


        !% initialize selecteros for upstream and downstream branches
        upBranchSelector = zeroI
        dnBranchSelector = zeroI

        !% cycle through the junction elements of map faces
        do ii = 1,Nelem_in_Junction

            !% now we are considering all the junction elements including
            !% junction main.

            !% all the even numbers are upstream branch elements
            if ((ii == 2) .or. (ii == 4) .or. (ii == 6)) then

                upBranchSelector = upBranchSelector + oneI
                !% pointer to upstream branch
                upBranchIdx => node%I(thisJNode,ni_idx_base1 + upBranchSelector)

                !% condition for a link connecting this branch is valid and
                !% included in this partition.
                if ( (upBranchIdx /= nullvalueI)   .and.       &
                     (link%I(upBranchIdx,li_P_image) == image) ) then

                    !% find the last element index of the link
                    LinkLastElem = link%I(upBranchIdx,li_last_elem_idx)

                    !% pointer to the specific branch element
                    eIdx => JelemIdx(ii)

                    !% find the downstream face index of that last element
                    fLidx => elemI(eIdx,ei_Mface_uL)

                    !% if the face is a shared face across images,
                    !% it will not have any upstream local element
                    if ( .not. faceYN(fLidx,fYN_isSharedFace)) then

                        !% the upstream face of the upstream branch will be the
                        !% last downstream face of the connected link
                        !% here, one important thing to remember is that
                        !% the upstrem branch elements does not have any
                        !% downstream faces.

                        !% local d/s face map to element u/s of the branch
                        elemI(LinkLastElem,ei_Mface_dL) = fLidx

                        !% local u/s element of the face
                        faceI(fLidx,fi_Melem_uL) = LinkLastElem
                    endif
                endif

            !% all odd numbers starting from 3 are downstream branch elements
            elseif ((ii == 3) .or. (ii == 5) .or. (ii == 7)) then

                dnBranchSelector = dnBranchSelector + oneI
                !% pointer to upstream branch
                dnBranchIdx => node%I(thisJNode,ni_idx_base2 + dnBranchSelector)

                !% condition for a link connecting this branch is valid and
                !% included in this partition.
                if ( (dnBranchIdx /= nullvalueI)   .and.       &
                     (link%I(dnBranchIdx,li_P_image) == image) ) then

                    !% find the first element index of the link
                    LinkFirstElem = link%I(dnBranchIdx,li_first_elem_idx)

                    !% pointer to the specific branch element
                    eIdx => JelemIdx(ii)

                    !% find the downstream face index of that last element
                    fLidx => elemI(eIdx,ei_Mface_dL)

                    !% if the face is a shared face across images,
                    !% it will not have any upstream local element
                    !% (not sure if we need this condition)
                    if ( .not. faceYN(fLidx,fYN_isSharedFace)) then

                        !% the downstream face of the downstream branch will be the
                        !% first upstream face of the connected link
                        !% here, one important thing to remember is that
                        !% the downstream branch elements does not have any
                        !% upstream faces.

                        !% local map to upstream face for elemI
                        elemI(LinkFirstElem,ei_Mface_uL) = fLidx

                        !% local downstream element of the face
                        faceI(fLidx,fi_Melem_dL) = LinkFirstElem
                    endif
                endif
            endif
        enddo

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name

    end subroutine init_network_map_nJm_branches
        !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_map_shared_nJm_nodes (image, fLidx, nIdx)
    !
    !--------------------------------------------------------------------------
    !
    !% set the global index, map, and ghost element for nJm nodes
    !
    !--------------------------------------------------------------------------
    !
        integer, intent(in) :: image, fLidx, nIdx

        integer             :: ii
        integer, pointer    :: fGidx, eUp, eDn, targetImage, branchIdx

        character(64) :: subroutine_name = 'init_network_map_shared_nJm_nodes'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name

        !% necessary pointers
        fGidx       => faceI(fLidx,fi_Gidx)
        eUp         => faceI(fLidx,fi_Melem_uL)
        eDn         => faceI(fLidx,fi_Melem_dL)
        targetImage => faceI(fLidx,fi_Connected_image)
        branchIdx   => faceI(fLidx,fi_link_idx)

        do ii = 1,N_face(targetImage)

            if ((faceI(ii,fi_Connected_image)[targetImage] == image) .and. &
                (faceI(ii,fi_node_idx)[targetImage] == nIdx        ) .and. &
                (faceI(ii,fi_link_idx)[targetImage] == branchIdx   ))   then

                !% find the local ghost element index of the connected image
                if (faceYN(ii,fYN_isUpGhost)[targetImage]) then
                    faceI(ii,fi_GhostElem_uL)[targetImage] = eUp

                elseif (faceYN(ii,fYN_isDnGhost)[targetImage]) then
                    faceI(ii,fi_GhostElem_dL)[targetImage] = eDn
                endif

                !% find the global index and set to target image
                if(faceI(ii,fi_Gidx)[targetImage] == nullvalueI) then
                    faceI(ii,fi_Gidx)[targetImage] = fGidx
                endif
            endif
        enddo

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name
    end subroutine init_network_map_shared_nJm_nodes
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_map_nJ2 (image, thisJNode)
    !
    !--------------------------------------------------------------------------
    !
    !% map all the nJ2 nodes. All the nJ2 node maps are handeled in the partition
    !% this is for the special cases where a disconnected nJ2 has not been mapped
    !% properly
    !
    !--------------------------------------------------------------------------
    !
        integer, intent(in) :: image, thisJNode

        integer             :: ii, upBranchSelector, dnBranchSelector
        integer             :: LinkFirstElem, LinkLastElem
        integer, pointer    :: upBranchIdx, dnBranchIdx
        integer, pointer    :: eIdx, fLidx

        character(64) :: subroutine_name = 'init_network_map_nJ2'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name


        !% initialize selecteros for upstream and downstream branches
        upBranchSelector = zeroI
        dnBranchSelector = zeroI

        upBranchSelector = upBranchSelector + oneI
        !% pointer to upstream branch
        upBranchIdx => node%I(thisJNode,ni_idx_base1 + upBranchSelector)

        !% condition for a link included in this partition.
        if (link%I(upBranchIdx,li_P_image) == image) then

            !% find the last element index of the link
            LinkLastElem = link%I(upBranchIdx,li_last_elem_idx)

            !% find the downstream face index of that last element
            fLidx => node%I(thisJNode,ni_elemface_idx)

            !% if the face is a shared face across images,
            !% it will not have any upstream local element
            if ( .not. faceYN(fLidx,fYN_isSharedFace)) then
                !% local d/s face map to element u/s of the branch
                elemI(LinkLastElem,ei_Mface_dL) = fLidx

                !% local u/s element of the face
                faceI(fLidx,fi_Melem_uL) = LinkLastElem
            endif
        endif

        dnBranchSelector = dnBranchSelector + oneI
        !% pointer to upstream branch
        dnBranchIdx => node%I(thisJNode,ni_idx_base2 + dnBranchSelector)

        !% condition for a link included in this partition.
        if (link%I(dnBranchIdx,li_P_image) == image) then

            !% find the first element index of the link
            LinkFirstElem = link%I(dnBranchIdx,li_first_elem_idx)

            !% find the downstream face index of that last element
            fLidx => node%I(thisJNode,ni_elemface_idx)

            !% if the face is a shared face across images,
            !% it will not have any downstream local element
            !% (not sure if we need this condition)
            if ( .not. faceYN(fLidx,fYN_isSharedFace)) then

                !% local map to upstream face for elemI
                elemI(LinkFirstElem,ei_Mface_uL) = fLidx

                !% local downstream element of the face
                faceI(fLidx,fi_Melem_dL) = LinkFirstElem
            endif
        endif

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name

    end subroutine init_network_map_nJ2
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_map_shared_nJ2_nodes (image, fLidx, nIdx)
    !
    !--------------------------------------------------------------------------
    !
    !% set the global index, map, and ghost element for nJ2 nodes
    !
    !--------------------------------------------------------------------------
    !
        integer, intent(in) :: image, fLidx, nIdx

        integer             :: ii
        integer, pointer    :: fGidx, eUp, eDn, targetImage
        logical, pointer    :: isUpGhost, isDnGhost

        character(64) :: subroutine_name = 'init_network_map_shared_nJ2_nodes'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name

        !% necessary pointers
        fGidx       => faceI(fLidx,fi_Gidx)
        eUp         => faceI(fLidx,fi_Melem_uL)
        eDn         => faceI(fLidx,fi_Melem_dL)
        targetImage => faceI(fLidx,fi_Connected_image)

        do ii = 1,N_face(targetImage)

            if ((faceI(ii,fi_Connected_image)[targetImage] == image) .and. &
                (faceI(ii,fi_node_idx)[targetImage] == nIdx)) then

                !% find the local ghost element index of the connected image
                if (faceYN(ii,fYN_isUpGhost)[targetImage]) then
                    faceI(ii,fi_GhostElem_uL)[targetImage] = eUp

                elseif (faceYN(ii,fYN_isDnGhost)[targetImage]) then
                    faceI(ii,fi_GhostElem_dL)[targetImage] = eDn
                endif

                !% find the global index and set to target image
                if(faceI(ii,fi_Gidx)[targetImage] == nullvalueI) then
                    faceI(ii,fi_Gidx)[targetImage] = fGidx
                endif
            endif
        enddo

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name
    end subroutine init_network_map_shared_nJ2_nodes
    !
    !==========================================================================
    !==========================================================================
    !
    function init_network_nJm_branch_length (LinkIdx) result (BranchLength)
    !--------------------------------------------------------------------------
    !
    !% compute the length of a junction branch
    !
    !--------------------------------------------------------------------------

        integer, intent(in)  :: LinkIdx
        real(8)              :: BranchLength

        character(64) :: subroutine_name = 'init_network_nJm_branch_length'
    !--------------------------------------------------------------------------
        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name

        !% find the length of the junction branch
        if (link%I(LinkIdx,li_length_adjusted) == OneSideAdjust) then
            BranchLength = link%R(LinkIdx,lr_Length) - link%R(LinkIdx,lr_AdjustedLength)
        elseif (link%I(LinkIdx,li_length_adjusted) == BothSideAdjust) then
            BranchLength = (link%R(LinkIdx,lr_Length) - link%R(LinkIdx,lr_AdjustedLength))/twoR
        endif

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name
    end function init_network_nJm_branch_length
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_nullify_nJm_branch (ElemIdx, FaceIdx)
    !--------------------------------------------------------------------------
    !
    !% set all the values to zero for a null junction
    !
    !--------------------------------------------------------------------------

        integer, intent(in)  :: ElemIdx, FaceIdx

        character(64) :: subroutine_name = 'init_network_nullify_nJm_branch'
    !--------------------------------------------------------------------------
        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name

        !% set everything to zero for a non existant branch
        elemR(ElemIdx,:)                            = zeroR
        elemSR(ElemIdx,:)                           = zeroR
        elemSGR(ElemIdx,:)                          = zeroR
        elemSI(ElemIdx,eSI_JunctionBranch_Exists)   = zeroI
        faceR(FaceIdx,:)                            = zeroR
        faceYN(FaceIdx,fYN_isnull)                  = .true.

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name
    end subroutine init_network_nullify_nJm_branch
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_set_interior_faceYN ()
    !
    !--------------------------------------------------------------------------
    !
    !% set the logicals of fYN_isInteriorFace
    !
    !--------------------------------------------------------------------------
        character(64) :: subroutine_name = 'init_network_set_interior_faceYN'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ', this_image(),subroutine_name

        where ( (faceI(:,fi_BCtype)         ==  doesnotexist) &
                .and. &
                (faceYN(:,fYN_isnull)       .eqv. .false.     ) &
                .and. &
                (faceYN(:,fYN_isSharedFace) .eqv. .false.     ) )

            faceYN(:,fYN_isInteriorFace) = .true.
        endwhere

        if (setting%Debug%File%network_define) print *, '*** leave ', this_image(),subroutine_name
    end subroutine init_network_set_interior_faceYN
    !
    !==========================================================================
    ! END OF MODULE
    !==========================================================================
    !
end module network_define

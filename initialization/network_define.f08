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

    public :: init_network

contains
    !
    !==========================================================================
    ! PUBLIC
    !==========================================================================
    !
    subroutine init_network ()
    !--------------------------------------------------------------------------
    !
    !% Initializes a element-face network from a link-node network.
    !%   Requires network links and nodes before execution
    !
    !--------------------------------------------------------------------------

        integer :: ii

        character(64) :: subroutine_name = 'init_network'

    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name


        !% get the slope of each link given the node Z values
        call init_network_linkslope()

        !% divide the link node networks in elements and faces
        call init_network_datacreate()

        sync all
        !% print result
        if (setting%Debug%File%network_define) then

            !% only using the first processor to print results
            if (this_image() == 1) then

                do ii = 1,num_images()
                   print*, '----------------------------------------------------'
                   print*, 'image = ', ii
                   print*, '..................elements..........................'
                   print*, elemI(:,ei_Lidx)[ii], 'Lidx'
                   print*, elemI(:,ei_Gidx)[ii], 'Gidx'
                   print*, elemI(:,ei_elementType)[ii], 'elementType'
                   print*, elemI(:,ei_geometryType)[ii], 'geometryType'
                   print*, elemI(:,ei_link_Gidx_SWMM)[ii], 'link_Gidx_SWMM'
                   print*, elemI(:,ei_node_Gidx_SWMM)[ii], 'node_Gidx_SWMM'
                   print*, elemI(:,ei_Mface_uL)[ii],'Mface_uL'
                   print*, elemI(:,ei_Mface_dL)[ii],'Mface_dL'
                   print*, '..................faces.............................'
                   print*, faceI(:,fi_Lidx)[ii], 'face Lidx'
                   print*, faceI(:,fi_Gidx)[ii], 'face Gidx'
                   print*, faceI(:,fi_Melem_dL)[ii], 'face Melem_dL'
                   print*, faceI(:,fi_Melem_uL)[ii], 'face Melem_uL'
                   print*, faceI(:,fi_Connected_image)[ii], 'fi_Connected_image'
                   print*, faceYN(:,fYN_isSharedFace)[ii], 'face is shared face'
                   print*, '----------------------------------------------------'
                   call execute_command_line('')
                enddo

            endif
        endif

        sync all

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine init_network
    !
    !==========================================================================
    ! PRIVATE
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

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        do mm = 1,N_link
            !% Inputs
            NodeUp      => linkI(mm,li_Mnode_u)
            NodeDn      => linkI(mm,li_Mnode_d)
            zUp         => nodeR(NodeUp,nr_Zbottom)
            zDn         => nodeR(NodeDn,nr_Zbottom)
            Length      => linkR(mm,lr_Length)
            !%-------------------------------------------------------------------------
            !% HACK: Original length is used for slope calculation instead of adjusted
            !% length. Using adjusted lenghts will induce errors in slope calculations
            !% and will distort the original network.
            !%-------------------------------------------------------------------------
            !% Output
            Slope       => linkR(mm,lr_Slope)

            Slope = (zUp - zDn) / Length
        end do

        if (setting%Debug%File%network_define) then
            !% provide output for debugging
            print *, subroutine_name,'--------------------------------'
            print *, 'link ID,               Slope,             length'
            do mm=1,N_link
                print *, mm, linkR(mm,lr_Slope), linkR(mm,lr_Length)
            end do
        endif

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

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

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% initializing global element and face id
        ElemGlobalIdx = first_elem_index
        FaceGlobalIdx = first_face_index

        !% initializing local element and face id
        ElemLocalIdx = first_elem_index
        FacelocalIdx = first_face_index

        !% Setting the local image value
        image = this_image()

        !% initalizing the global element and face id by looping through N_elem and N_face for images not equal to one

        if(this_image() /= 1) then
           do ii=1, this_image()-1
              ElemGlobalIdx = ElemGlobalIdx + N_elem(ii)
              !% we have to subtract one from the global face id such that faces along the image boundaries are shared.
              FaceGlobalIdx = FaceGlobalIdx + N_unique_face(ii)
           end do
        end if

        !% handle all the links and nodes in a partition
        call init_network_map_linknode &
            (image, ElemLocalIdx, FacelocalIdx, ElemGlobalIdx, FaceGlobalIdx)

        !% finish mapping all the junction branch and faces that were not
        !% handeled in handle_link_nodes subroutine
        call init_network_map_nJm (image)

        !% shared faces are mapped by copying data from different images
        !% thus a sync all is needed
        sync all

        !% finally set the same global face idx for shared faces across images
        call init_network_map_shared_faces (image)

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine init_network_datacreate
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_map_linknode &
        (image, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, FaceGlobalCounter)
    !
    !--------------------------------------------------------------------------
    !
    !% Traverse through all the links and nodes in a partition and creates
    !%   elements and faces.
    !
    !--------------------------------------------------------------------------
    !
        integer, intent(in)    :: image
        integer, intent(inout) :: ElemLocalCounter, FaceLocalCounter
        integer, intent(inout) :: ElemGlobalCounter, FaceGlobalCounter

        integer :: ii, pLink, pNode
        integer, pointer :: thisLink, upNode, dnNode
        integer, dimension(:), allocatable, target :: packed_link_idx, packed_node_idx

        character(64) :: subroutine_name = 'init_network_map_linknode'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% pack all the link indexes in an image
        packed_link_idx = pack(linkI(:,li_idx), (linkI(:,li_P_image) .eq. image))

        !% find the number of links in an image
        pLink = size(packed_link_idx)

        !% pack all the node indexes in a partition to determine which nodes are in the partition
        packed_node_idx = pack(nodeI(:,ni_idx), (nodeI(:,ni_P_image) .eq. image))

        !% number of nodes in a partition
        pNode = size(packed_node_idx)

        !% cycle through the links in an image
        do ii = 1,pLink
            !% necessary pointers to the links and connected nodes
            thisLink => packed_link_idx(ii)
            upNode   => linkI(thisLink,li_Mnode_u)
            dnNode   => linkI(thisLink,li_Mnode_d)

            call init_network_map_upstreamnode &
                (image, upNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                FaceGlobalCounter)

            call init_network_map_subdividelink &
                (thisLink, upNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                FaceGlobalCounter)

            call init_network_map_downstreamnode &
                (image, dnNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                FaceGlobalCounter)
        end do

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine init_network_map_linknode
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_map_nJm (image)
    !
    !--------------------------------------------------------------------------
    !
    !% map all the multi branch junction faces in an image
    !
    !--------------------------------------------------------------------------
    !
        integer, intent(in)    :: image

        integer :: ii, pnJm
        integer, pointer :: thisJunctionNode
        integer, dimension(:), allocatable, target :: packed_nJm_idx, JunctionElementIdx



        character(64) :: subroutine_name = 'init_network_map_nJm'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name


        !% pack all the nJm node indexes in a partition to find face maps
        packed_nJm_idx = pack(nodeI(:,ni_idx),                       &
                             ((nodeI(:,ni_P_image) .eq. image) .and. &
                              (nodeI(:,ni_node_type) .eq. nJm ) ) )

        !% number of nodes in a partition
        pnJm = size(packed_nJm_idx)

        !% cycle through the nJm nodes to find the face maps
        do ii = 1,pnJm
            thisJunctionNode => packed_nJm_idx(ii)
            JunctionElementIdx = pack( elemI(:,ei_Lidx), &
                                     ( elemI(:,ei_node_Gidx_SWMM) .eq. thisJunctionNode) )

            call init_network_map_nJm_branches (image, thisJunctionNode, JunctionElementIdx)

        end do

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine init_network_map_nJm
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

        integer :: ii, jj, targetImage, NsharedFaces, fGidx
        integer, pointer :: fLidx
        integer, dimension(:), allocatable, target ::  sharedFaces

        character(64) :: subroutine_name = 'init_network_map_shared_faces'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% pack the shared faces in an image
        sharedFaces = pack(faceI(:,fi_Lidx), faceYN(:,fYN_isSharedFace))

        !% fid the size of the pack
        NsharedFaces = size(sharedFaces)


        !% HACK: bellow is absolutely rubbish code
        !% it can be written in a better way
        !% but it works for now

        do ii = 1,NsharedFaces
            fLidx => sharedFaces(ii)
            !% find the target image
            targetImage = faceI(fLidx,fi_Connected_image)

            !% find the global index and set to target image
            if(faceI(fLidx,fi_Gidx) .ne. nullvalueI) then

                fGidx = faceI(fLidx,fi_Gidx)

                do jj = 1, size(faceI(:,fi_Lidx))

                    if (faceI(jj,fi_Connected_image)[targetImage] .eq. image) then
                        faceI(jj,fi_Gidx)[targetImage] = fGidx
                    endif
                enddo
            endif
        enddo

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine init_network_map_shared_faces
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_map_upstreamnode &
        (image, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
        FaceGlobalCounter)
    !
    !--------------------------------------------------------------------------
    !
    !% Handle the node upstream of a link
    !
    !--------------------------------------------------------------------------
    !
        integer, intent(in)    :: image, thisNode
        integer, intent(inout) :: ElemLocalCounter, FaceLocalCounter
        integer, intent(inout) :: ElemGlobalCounter, FaceGlobalCounter

        integer, pointer :: nAssignStatus, nodeType, linkUp

        integer :: ii

        character(64) :: subroutine_name = 'init_network_map_upstreamnode'
    !--------------------------------------------------------------------------
        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% check 1: Is the node is in the partition
        if (nodeI(thisNode,ni_P_image) .eq. image) then

            !% necessary pointers
            nAssignStatus => nodeI(thisNode,ni_assigned)
            nodeType      => nodeI(thisNode,ni_node_type)

            select case (nodeType)

                case(nBCup)
                    !% check if the node has already been assigned
                    if (nAssignStatus .eq. nUnassigned) then
                        !% integer data
                        faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                        faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                        faceI(FaceLocalCounter,fi_Melem_uL) = nullvalueI
                        faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter
                        faceI(FaceLocalCounter,fi_BCtype)   = fBCup

                        !% change the node assignmebt value
                        nAssignStatus =  nAssigned
                    endif

                case (nJ2)
                    !% check if the node has already been assigned
                    if (nAssignStatus .eq. nUnassigned) then
                        !% integer data
                        faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                        faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter

                        if (nodeI(thisNode,ni_P_is_boundary) .eq. EdgeNode) then

                            !% An upstream edge node indicates there are no local
                            !% elements upstream of that node
                            faceI(FaceLocalCounter,fi_Melem_uL) = nullvalueI

                            !% logical data
                            faceYN(FaceLocalCounter,fYN_isSharedFace) = .true.

                            !% find the connecting image to this face
                            linkUp  => nodeI(thisNode,ni_Mlink_u1)

                            faceI(FaceLocalCounter, fi_Connected_image) = linkI(linkUp,li_P_image)
                        else
                            !% if the node is not an edge node, there is an upstream
                            !% local link element which has already been handeled
                            faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                            faceI(FaceLocalCounter,fi_Melem_uL) = ElemLocalCounter - oneI
                        endif

                        !% change the node assignmebt value
                        nAssignStatus =  nAssigned
                    endif

                case (nJm)

                    !% check if the node has already been assigned
                    if (nAssignStatus .eq. nUnassigned) then

                        call init_network_map_subdivide_nJm &
                            (image, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                            FaceGlobalCounter, nAssignStatus)

                    endif

                case default

                    print*, 'In ', subroutine_name
                    print*, 'error: unexpected node, ', nodeType,'  at upstream boundary'
                    stop
            end select


        else

            !% if the upstream node is not in the partiton,
            !% the face map to upstream nullvalue

            !% integer data
            faceI(FacelocalCounter,fi_Lidx) = FacelocalCounter
            faceI(FacelocalCounter,fi_Gidx) = nullvalueI
            faceI(FacelocalCounter,fi_Connected_image) = nodeI(thisNode,ni_P_image)

            !% since no upstream node indicates start of a partiton,
            !% the downstream element will be initialized elem idx
            faceI(FacelocalCounter,fi_Melem_dL) = ElemLocalCounter

            !% since the first face is a shared face, the global
            !% index will be updated later.
            !% However, subdivide_link_going_downstream will advance
            !% the global face counter without considering this update
            !% later. Thus, this index should be adjusted by subtracting
            !% one from the FaceGlobalCounter

            !% HACK: i dont know if it will work for different networks

            FaceGlobalCounter = FaceGlobalCounter - oneI

            !% logical data
            faceYN(FacelocalCounter,fYN_isSharedFace) = .true.
        endif

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine init_network_map_upstreamnode
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_map_subdividelink &
        (thisLink, upNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, FaceGlobalCounter)
    !--------------------------------------------------------------------------
    !
    !% subdivides the links into elements going downstream
    !
    !--------------------------------------------------------------------------

        integer, intent(in)     :: thisLink, upNode
        integer, intent(inout)  :: ElemLocalCounter, FaceLocalCounter
        integer, intent(inout)  :: ElemGlobalCounter, FaceGlobalCounter

        integer                 :: ii, NlinkElem
        integer, pointer        :: lAssignStatus, NodeUp
        real                    :: zUpstream, zCenter

        character(64) :: subroutine_name = 'init_network_map_subdividelink'

    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% necessary pointers
        !% link assign status
        lAssignStatus => linkI(thisLink,li_assigned)

        NlinkElem  = linkI(thisLink,li_N_element)
        zUpstream  =  nodeR(upNode,nr_Zbottom)

        if (lAssignStatus .eq. lUnassigned) then


            !%  store the ID of the first (upstream) element in this link
            linkI(thisLink,li_first_elem_idx)   = ElemLocalCounter
            linkI(thisLink,li_Melem_u)          = ElemLocalCounter - oneI ! the new element is the 1st (downstream)
            linkI(thisLink,li_Mface_u)          = FaceLocalCounter        ! the old face is the 1st

            if (linkI(thisLink,li_Melem_u) .eq. zeroI) then
                !% HACK: needs further testing
                !% a zero value of Melem_u indicates there isnt
                !% any local element upstream of that link
                linkI(thisLink,li_Melem_u) = nullvalueI
            endif

            !%  reference elevations at cell center and cell face
            zCenter = zUpstream - 0.5 * linkR(thisLink,lr_ElementLength) * linkR(thisLink,lr_Slope)

            do ii = 1,NlinkElem
                !%................................................................
                !% Element arrays update
                !%................................................................

                !% integer data
                elemI(ElemLocalCounter,ei_Lidx)             = ElemLocalCounter
                elemI(ElemLocalCounter,ei_Gidx)             = ElemGlobalCounter
                elemI(ElemLocalCounter,ei_geometryType)     = linkI(thisLink,li_geometry)
                elemI(ElemLocalCounter,ei_elementType)      = linkI(thisLink,li_link_type)
                elemI(ElemLocalCounter,ei_link_Gidx_SWMM)   = thisLink
                elemI(ElemLocalCounter,ei_Mface_uL)         = FaceLocalCounter
                elemI(ElemLocalCounter,ei_Mface_dL)         = FaceLocalCounter + oneI

                !% real data
                elemR(ElemLocalCounter,er_Length) = linkR(thisLink,lr_ElementLength)
                elemR(ElemLocalCounter,er_Zcrown) = zCenter

                !%................................................................
                !% Face arrays update
                !%................................................................

                !% advance the face counters
                FaceLocalCounter  = FaceLocalCounter  + oneI
                FaceGlobalCounter = FaceGlobalCounter + oneI

                !% face integer data
                faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter + oneI
                faceI(FaceLocalCounter,fi_Melem_uL) = ElemLocalCounter

                !% counter for element z bottom calculation
                zCenter = zCenter - linkR(thisLink,lr_ElementLength) * linkR(thisLink,lr_Slope)

                !% Advance the element counter
                ElemLocalCounter  = ElemLocalCounter  + oneI
                ElemGlobalCounter = ElemGlobalCounter + oneI
            end do

            lAssignStatus = lAssigned


            linkI(thisLink,li_last_elem_idx)    = ElemLocalCounter - oneI

            ! linkI(thisLink,li_Melem_d)          = ElemLocalCounter        ! HACK: This may not be right in some cases
            ! linkI(thisLink,li_Mface_d)          = FaceLocalCounter + oneI ! HACK: This may not be right in some cases
            !% HACK:
            !% the last face idx may need to be corrected based on junction element
        endif

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine init_network_map_subdividelink
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_map_downstreamnode &
        (image, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
        FaceGlobalCounter)
    !--------------------------------------------------------------------------
    !
    !% handle the node downstream of a link
    !
    !--------------------------------------------------------------------------
    !
        integer, intent(in)    :: image, thisNode
        integer, intent(inout) :: ElemLocalCounter, FaceLocalCounter
        integer, intent(inout) :: ElemGlobalCounter, FaceGlobalCounter

        integer, pointer :: nAssignStatus, nodeType, linkDn

        character(64) :: subroutine_name = 'init_network_map_downstreamnode'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% check 1: Is the node is in the partition
        if (nodeI(thisNode,ni_P_image) .eq. image) then

            !% necessary pointers
            nAssignStatus => nodeI(thisNode,ni_assigned)
            nodeType      => nodeI(thisNode,ni_node_type)

            select case (nodeType)

                case(nBCdn)
                    !% check if the node has already been assigned
                    if (nAssignStatus .eq. nUnassigned) then
                        !% by the time we reach a downstream boundary
                        !% node, all the indexes have already been set
                        !% from the subdivide_link_going_downstream.
                        !% only the map to downstream element is
                        !% needed to be fixed.

                        !% integer data
                        faceI(FaceLocalCounter,fi_Melem_dL) = nullvalueI
                        faceI(FaceLocalCounter,fi_BCtype)   = fBCup

                        !% change the node assignmebt value
                        nAssignStatus =  nAssigned
                    endif

                case (nJ2)
                    !% check if the node has already been assigned
                    if (nAssignStatus .eq. nUnassigned) then
                        !% by the time we reach a downstream nJ2 node,
                        !% all the indexes have already been set from
                        !% subdivide_link_going_downstream. only the map
                        !% to downstream element is needed to be fixed if
                        !% the node is an edge node.

                        !% integer data
                        if (nodeI(thisNode,ni_P_is_boundary) .eq. EdgeNode) then

                            !% An downstream edge node indicates there are no local
                            !% elements downstream of that node
                            faceI(FaceLocalCounter,fi_Melem_dL) = nullvalueI

                            !% logical data
                            faceYN(FaceLocalCounter,fYN_isSharedFace) = .true.

                            !% find the connecting image to this face
                            linkDn  => nodeI(thisNode,ni_Mlink_d1)

                            faceI(FaceLocalCounter, fi_Connected_image) = linkI(linkDn,li_P_image)
                        endif

                        !% change the node assignmebt value
                        nAssignStatus =  nAssigned
                    endif

                case (nJm)
                    !% check if the node has already been assigned
                    if (nAssignStatus .eq. nUnassigned) then

                        call init_network_map_subdivide_nJm &
                            (image, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                            FaceGlobalCounter, nAssignStatus)

                    endif

                case default

                    print*, 'In ', subroutine_name
                    print*, 'error: unexpected node, ', nodeType,'  at upstream boundary'
                    stop
            end select
        else
            !% if the downstream node is not in the partiton.
            !% through subdivide_link_going_downstream subroutine
            !% upstream map to the element has alrady been set.
            !% However, downstream map has set to wrong value.
            !% Thus, setting the map elem ds as nullvaleI
            !% integer data
            faceI(FaceLocalCounter,fi_Melem_dL) = nullvalueI
            faceI(FacelocalCounter,fi_Gidx)     = nullvalueI
            faceI(FacelocalCounter,fi_Connected_image) = nodeI(thisNode,ni_P_image)

            !% logical data
            faceYN(FacelocalCounter,fYN_isSharedFace) = .true.
        endif

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine init_network_map_downstreamnode
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_map_subdivide_nJm &
        (image, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
        FaceGlobalCounter, nAssignStatus)
    !--------------------------------------------------------------------------
    !
    !% subdivides the links into elements going downstream
    !
    !--------------------------------------------------------------------------

        integer, intent(in)    :: image, thisNode
        integer, intent(inout) :: ElemLocalCounter, FaceLocalCounter
        integer, intent(inout) :: ElemGlobalCounter, FaceGlobalCounter
        integer, intent(inout) :: nAssignStatus

        integer, pointer :: upBranchIdx, dnBranchIdx

        integer :: ii, upBranchSelector, dnBranchSelector

        character(64) :: subroutine_name = 'init_network_map_subdivide_nJm'

    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !%................................................................
        !% Junction main
        !%................................................................
        !% Element Arrays
        !% integer data
        elemI(ElemLocalCounter,ei_Lidx) = ElemLocalCounter
        elemI(ElemLocalCounter,ei_Gidx) = ElemGlobalCounter
        elemI(ElemLocalCounter,ei_elementType) = eJunctionMain
        elemI(ElemLocalCounter,ei_node_Gidx_SWMM) = thisNode

        !% real data
        elemR(ElemLocalCounter,er_Zbottom) = nodeR(thisNode,nr_zbottom)

        !% Advance the element counter to 1st upstream branch
        ElemLocalCounter  = ElemLocalCounter  + oneI
        ElemGlobalCounter = ElemGlobalCounter + oneI

        !%................................................................
        !% Handle Junctin Branches
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
            elemI(ElemLocalCounter,ei_elementType)    = eJunctionBranch
            elemI(ElemLocalCounter,ei_node_Gidx_SWMM) = thisNode

            !% real data
            elemR(ElemLocalCounter,er_Zbottom) = nodeR(thisNode,nr_zbottom)

            !% face arrays
            !% integer data
            !% HACK: for non-edge node, the first face indexes
            !% should already be updated

            !% Now moving into branch specific calculations

            !%......................................................
            !% Upstream Branches
            !%......................................................
            if ((ii .eq. 1) .or. (ii .eq. 3) .or. (ii .eq. 5)) then
            !% all the odd numbers are upstream branches
                upBranchSelector = upBranchSelector + oneI
                !% pointer to upstream branch
                upBranchIdx => nodeI(thisNode,ni_idx_base1 + upBranchSelector)

                if (upBranchIdx .ne. nullvalueI) then
                    !% integer data
                    elemSI(ElemLocalCounter,eSI_JunctionBranch_Exists)  = oneI
                    elemR(ElemLocalCounter,er_Length) = get_junction_length(upBranchIdx)

                    !% check if the link connecting this branch
                    !% is a part of this partition
                    if ( (nodeI(thisNode,ni_P_is_boundary) .eq. EdgeNode)  .and. &
                         (linkI(upBranchIdx,li_P_image) .ne. image) )  then

                        !% HACK:
                        !% face counters are always advanced by link
                        !% elements unless a branch does not have any
                        !% elements associated with it. if the links
                        !% connected to the junction branch is in a different
                        !% partition, the local face count will not be advanced.
                        !% thus to keep the count consistant, face is advanced

                        !% advance the face counters for next branch
                        FaceLocalCounter  = FaceLocalCounter  + oneI
                        FaceGlobalCounter = FaceGlobalCounter + oneI

                        !% elem array
                        elemI(ElemLocalCounter,ei_Mface_uL) = FaceLocalCounter

                        !% face array
                        !% integer data
                        faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                        faceI(FacelocalCounter,fi_Gidx)     = FaceGlobalCounter
                        faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter
                        faceI(FaceLocalCounter,fi_Connected_image) = linkI(upBranchIdx,li_P_image)


                        !% logical data
                        faceYN(FacelocalCounter,fYN_isSharedFace) = .true.

                    endif
                else

                    !% HACK:
                    !% face counters are always advanced by link
                    !% elements unless a branch does not have any
                    !% elements associated with it. for null branches,
                    !% no links are connected to the junction branch.
                    !% to keep the count consistant, face is advanced

                    !% advance the face counters for next branch
                    FaceLocalCounter  = FaceLocalCounter  + oneI
                    FaceGlobalCounter = FaceGlobalCounter + oneI

                    !% elem array
                    elemI(ElemLocalCounter,ei_Mface_uL) = FaceLocalCounter

                    !% face array
                    !% integer data
                    faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                    faceI(FacelocalCounter,fi_Gidx)     = FaceGlobalCounter
                    faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter

                    call nullify_junction_branch &
                        (ElemLocalCounter, FaceLocalCounter)

                endif
            !%......................................................
            !% Downstream Branches
            !%......................................................
            else
                !% all the even numbers are downstream branches
                dnBranchSelector = dnBranchSelector + oneI
                !% pointer to upstream branch
                dnBranchIdx => nodeI(thisNode,ni_idx_base2 + dnBranchSelector)

                if (dnBranchIdx .ne. nullvalueI) then
                    !% integer data
                    elemSI(ElemLocalCounter,eSI_JunctionBranch_Exists)  = oneI
                    elemR(ElemLocalCounter,er_Length) = get_junction_length(dnBranchIdx)

                    !% check if the link connecting this branch
                    !% is a part of this partition
                    if ( (nodeI(thisNode,ni_P_is_boundary) .eq. EdgeNode)  .and. &
                         (linkI(dnBranchIdx,li_P_image) .ne. image) )  then

                        !% HACK:
                        !% faces are always advanced by link elements
                        !% however, if there aren't any link element
                        !% connected to the junction branch in a partition,
                        !% the face is advances

                        !% advance the face counters for next branch
                        FaceLocalCounter  = FaceLocalCounter  + oneI
                        FaceGlobalCounter = FaceGlobalCounter + oneI

                        !% elem array
                        !% integer data
                        elemI(ElemLocalCounter,ei_Mface_dL) = FaceLocalCounter

                        !% face array
                        !% integer data
                        faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                        faceI(FacelocalCounter,fi_Gidx)     = FaceGlobalCounter
                        faceI(FaceLocalCounter,fi_Melem_uL) = ElemLocalCounter
                        faceI(FaceLocalCounter,fi_Connected_image) = linkI(dnBranchIdx,li_P_image)

                        !% logical data
                        faceYN(FacelocalCounter,fYN_isSharedFace) = .true.

                    endif

                else
                    !% HACK:
                    !% faces are always advanced by link elements
                    !% however, for null branches, no links are
                    !% connected to the junction branch. To keep
                    !% the count consistant, face is advanced

                    !% advance the face counters for next branch
                    FaceLocalCounter  = FaceLocalCounter  + oneI
                    FaceGlobalCounter = FaceGlobalCounter + oneI

                    !% elem array
                    !% integer data
                    elemI(ElemLocalCounter,ei_Mface_dL) = FaceLocalCounter

                    !% face array
                    !% integer data
                    faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                    faceI(FacelocalCounter,fi_Gidx)     = FaceGlobalCounter
                    faceI(FaceLocalCounter,fi_Melem_uL) = ElemLocalCounter

                    call nullify_junction_branch &
                        (ElemLocalCounter, FaceLocalCounter)
                endif
            endif

            !% Advance the element counter for next branch
            ElemLocalCounter  = ElemLocalCounter  + oneI
            ElemGlobalCounter = ElemGlobalCounter + oneI
        end do

        !% set status to assigned
        nAssignStatus = nAssigned


        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine init_network_map_subdivide_nJm
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

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name


        !% initialize selecteros for upstream and downstream branches
        upBranchSelector = zeroI
        dnBranchSelector = zeroI

        !% cycle through the junction elements of map faces
        do ii = 1,Nelem_in_Junction

            !% now we are considering all the junction elements including
            !% junction main.

            !% all the even numbers are upstream branch elements
            if ((ii .eq. 2) .or. (ii .eq. 4) .or. (ii .eq. 6)) then

                upBranchSelector = upBranchSelector + oneI
                !% pointer to upstream branch
                upBranchIdx => nodeI(thisJNode,ni_idx_base1 + upBranchSelector)

                !% condition for a link connecting this branch is valid and
                !% included in this partition.
                if ( (upBranchIdx .ne. nullvalueI)   .and.       &
                     (linkI(upBranchIdx,li_P_image) .eq. image) ) then

                    !% find the last element index of the link
                    LinkLastElem = linkI(upBranchIdx,li_last_elem_idx)

                    !% find the downstream face index of that last element
                    fLidx => elemI(LinkLastElem,ei_Mface_dL)

                    !% pointer to the specific branch element
                    eIdx => JelemIdx(ii)

                    !% if the face is a shared face across images,
                    !% it will not have any upstream local element
                    if ( .not. faceYN(fLidx,fYN_isSharedFace)) then

                        !% the upstream face of the upstream branch will be the
                        !% last downstream face of the connected link
                        !% here, one important thing to remember is that
                        !% the upstrem branch elements does not have any
                        !% downstream faces.

                        !% local map to upstream face for elemI
                        elemI(eIdx,ei_Mface_uL) = fLidx

                        !% local downstream element of the face
                        faceI(fLidx,fi_Melem_dL) = eIdx

                    endif

                endif

            !% all odd numbers starting from 3 are downstream branch elements
            elseif ((ii .eq. 3) .or. (ii .eq. 5) .or. (ii .eq. 7)) then

                dnBranchSelector = dnBranchSelector + oneI
                !% pointer to upstream branch
                dnBranchIdx => nodeI(thisJNode,ni_idx_base2 + dnBranchSelector)

                !% condition for a link connecting this branch is valid and
                !% included in this partition.
                if ( (dnBranchIdx .ne. nullvalueI)   .and.       &
                     (linkI(dnBranchIdx,li_P_image) .eq. image) ) then

                    !% find the last element index of the link
                    LinkFirstElem = linkI(dnBranchIdx,li_first_elem_idx)

                    !% find the downstream face index of that last element
                    fLidx => elemI(LinkFirstElem,ei_Mface_dL)

                    !% pointer to the specific branch element
                    eIdx => JelemIdx(ii)

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
                        elemI(eIdx,ei_Mface_dL) = fLidx

                        !% local downstream element of the face
                        faceI(fLidx,fi_Melem_uL) = eIdx

                    endif

                endif

            endif

        enddo

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine init_network_map_nJm_branches
    !
    !==========================================================================
    !==========================================================================
    !
    function get_junction_length (LinkIdx) result (BranchLength)
    !--------------------------------------------------------------------------
    !
    !% compute the length of a junction branch
    !
    !--------------------------------------------------------------------------

        integer, intent(in)  :: LinkIdx
        real(8)              :: BranchLength

        character(64) :: subroutine_name = 'get_junction_length'
    !--------------------------------------------------------------------------
        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

        !% find the length of the junction branch
        if (linkI(LinkIdx,li_length_adjusted) .eq. OneSideAdjust) then

            BranchLength = linkR(LinkIdx,lr_Length) - linkR(LinkIdx,lr_AdjustedLength)

        elseif (linkI(LinkIdx,li_length_adjusted) .eq. BothSideAdjust) then

            BranchLength = (linkR(LinkIdx,lr_Length) - linkR(LinkIdx,lr_AdjustedLength))/twoR
        endif

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name
    end function get_junction_length
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine nullify_junction_branch (ElemIdx, FaceIdx)
    !--------------------------------------------------------------------------
    !
    !% set all the values to zero for a null junction
    !
    !--------------------------------------------------------------------------

        integer, intent(in)  :: ElemIdx, FaceIdx

        character(64) :: subroutine_name = 'nullify_junction_branch'
    !--------------------------------------------------------------------------
        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

        !% set everything to zero for a non existant branch
        elemR(ElemIdx,:)                            = zeroR
        elemSR(ElemIdx,:)                           = zeroR
        elemSI(ElemIdx,eSI_JunctionBranch_Exists)   = zeroI
        faceR(FaceIdx,:)                            = zeroR
        faceYN(FaceIdx,fYN_isnull)                  = .true.

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name
    end subroutine nullify_junction_branch
    !
    !==========================================================================
    ! END OF MODULE
    !==========================================================================
    !
end module network_define

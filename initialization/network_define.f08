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

        !% print results
        if (setting%Debug%File%network_define) then
            !% only using the first processor to print results
            if (this_image() == 1) then

                do ii = 1,num_images()
                   print*, '----------------------------------------------------'
                   print*, 'image = ', ii
                   print*, '..................elements..........................'
                   print*, elemI(:,ei_Lidx)[ii], 'Lidx'
                   print*, elemI(:,ei_Gidx)[ii], 'Gidx'
                   print*, elemI(:,ei_link_Gidx_SWMM)[ii], 'link_Gidx_SWMM'
                   print*, elemI(:,ei_node_Gidx_SWMM)[ii], 'node_Gidx_SWMM'
                   print*, elemI(:,ei_Mface_uL)[ii],'Mface_uL'
                   print*, elemI(:,ei_Mface_dL)[ii],'Mface_dL'
                   print*, '..................faces.............................'
                   print*, faceI(:,fi_Lidx)[ii], 'face Lidx'
                   print*, faceI(:,fi_Gidx)[ii], 'face Gidx'
                   print*, faceI(:,fi_Melem_uL)[ii], 'face Melem_uL'
                   print*, faceI(:,fi_Melem_dL)[ii], 'face Melem_dL'
                   print*, faceI(:,fi_Connected_image)[ii], 'fi_Connected_image'
                   print*, faceI(:,fi_BCtype)[ii], 'fi_BCtype'
                   print*, faceI(:,fi_GhostElem_uL)[ii], 'fi_GhostElem_uL'
                   print*, faceI(:,fi_GhostElem_dL)[ii], 'fi_GhostElem_dL'
                   print*, faceYN(:,fYN_isInteriorFace)[ii], 'face is interior'
                   print*, faceYN(:,fYN_isSharedFace)[ii], 'face is shared'
                   print*, '----------------------------------------------------'
                   call execute_command_line('')
                enddo
            endif
        endif

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

        !% initialize the global indexes of elements and faces
        call init_network_set_global_indexes &
            (image, ElemGlobalIdx, FaceGlobalIdx)

        !% set the dummy element
        call init_network_set_dummy_elem()

        !% handle all the links and nodes in a partition
        call init_network_handle_partition &
            (image, ElemLocalIdx, FacelocalIdx, ElemGlobalIdx, FaceGlobalIdx)

        !% finish mapping all the junction branch and faces that were not
        !% handeled in handle_link_nodes subroutine
        call init_network_map_nJm (image)

        !% set interior face logical
        call init_network_set_interior_faceYN ()

        !% shared faces are mapped by copying data from different images
        !% thus a sync all is needed
        sync all

        !% set the same global face idx for shared faces across images
        call init_network_map_shared_faces (image)

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name
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

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name
 
        if(image /= 1) then
           do ii=1, image-1
              ElemGlobalCounter = ElemGlobalCounter + N_elem(ii)
              FaceGlobalCounter = FaceGlobalCounter + N_unique_face(ii)
           end do
        end if

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name
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

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% index for the dummy element
        dummyIdx = max_caf_elem_N + N_dummy_elem

        !% set the elem arrays with dummy values
        elemI(dummyIdx, ei_Lidx)        = dummyIdx
        elemI(dummyIdx,ei_elementType)  = dummy

        elemYN(dummyIdx,eYN_isDummy)    = .true.

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name
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
        logical                 :: firstUpBcHandeled
        integer, dimension(:), allocatable, target :: packed_link_idx

        character(64) :: subroutine_name = 'init_network_handle_partition'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% pack all the link indexes in a partition
        packed_link_idx = pack(linkI(:,li_idx), (linkI(:,li_P_image) .eq. image))

        !% find the number of links in a partition
        pLink = size(packed_link_idx)

        !% initializing first upstream boundary node handeled as false.
        !% if a partition has multiple upstream boundary node, this flag
        !% will keep the count consistant.
        firstUpBcHandeled = .false.

        !% cycling through all the links in a partition
        do ii = 1,pLink
            !% necessary pointers to the link and its connected nodes
            thisLink => packed_link_idx(ii)
            upNode   => linkI(thisLink,li_Mnode_u)
            dnNode   => linkI(thisLink,li_Mnode_d)

            !% handle the upstream node of the link to create elements and faces
            call init_network_handle_upstreamnode &
                (image, upNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                FaceGlobalCounter,firstUpBcHandeled)

            !% handle the link to create elements and faces
            call init_network_handle_link &
                (image, thisLink, upNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                FaceGlobalCounter)

            !% handle the downstream node of the link to create elements and faces
            call init_network_handle_downstreamnode &
                (image, dnNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                FaceGlobalCounter)
        end do

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name
    end subroutine init_network_handle_partition
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
        packed_nJm_idx = pack(nodeI(:,ni_idx),                         &
                             ((nodeI(:,ni_P_image)   .eq. image) .and. &
                              (nodeI(:,ni_node_type) .eq. nJm  ) ) )

        !% number of nJm nodes in a partition
        pnJm = size(packed_nJm_idx)

        !% cycle through all the nJm nodes and set the face maps
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

        integer :: ii, jj, NsharedFaces
        integer, pointer :: fLidx, fGidx, eUp, eDn, targetImage
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
            fGidx => faceI(fLidx,fi_Gidx)
            eUp   => faceI(fLidx,fi_Melem_uL)
            eDn   => faceI(fLidx,fi_Melem_dL)
            !% find the target image
            targetImage => faceI(fLidx,fi_Connected_image)
            
            do jj = 1, size(faceI(:,fi_Lidx))
                if (faceI(jj,fi_Connected_image)[targetImage] .eq. image) then

                    !% find the local ghost element index of the connected image
                    if (faceYN(jj,fYN_isUpGhost)[targetImage]) then
                        faceI(jj,fi_GhostElem_uL)[targetImage] = eUp

                    elseif (faceYN(jj,fYN_isDnGhost)[targetImage]) then
                        faceI(jj,fi_GhostElem_dL)[targetImage] = eDn
                    endif

                    !% find the global index and set to target image
                    if(faceI(fLidx,fi_Gidx) .ne. nullvalueI) then
                        faceI(jj,fi_Gidx)[targetImage] = fGidx
                    endif
                endif
            enddo
        enddo


        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name
    end subroutine init_network_map_shared_faces
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_handle_upstreamnode &
        (image, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
        FaceGlobalCounter,firstUpBcHandeled)
    !
    !--------------------------------------------------------------------------
    !
    !% Handle the node upstream of a link to create elements and faces
    !
    !--------------------------------------------------------------------------
    !
        integer, intent(in)     :: image, thisNode
        integer, intent(inout)  :: ElemLocalCounter, FaceLocalCounter
        integer, intent(inout)  :: ElemGlobalCounter, FaceGlobalCounter
        logical, intent(inout)  :: firstUpBcHandeled

        integer                 :: ii
        integer, pointer        :: nAssignStatus, nodeType, linkUp

        character(64) :: subroutine_name = 'init_network_handle_upstreamnode'
    !--------------------------------------------------------------------------
        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% Check 1: If the node is in the partition
        if (nodeI(thisNode,ni_P_image) .eq. image) then

            !% necessary pointers
            nAssignStatus => nodeI(thisNode,ni_assigned)
            nodeType      => nodeI(thisNode,ni_node_type)

            select case (nodeType)

                !% Handle upstream boundary nodes
                case(nBCup)
                    !% Check 2: If the node has already been assigned
                    if (nAssignStatus .eq. nUnassigned) then

                        !% Check 3: If there are multiple UpBc 
                        if (firstUpBcHandeled) then
                            !% if first boundary condition already assigned this condition will be
                            !% true. Thus, if the handler finds other upstream boundary node in the
                            !% images, and the face counder will be advanced.
                            FaceLocalCounter  = FaceLocalCounter + oneI
                            FaceGlobalCounter = FaceGlobalCounter + oneI
                        endif

                        !% integer data
                        faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                        faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                        !% an upstream boundary face does not have any local upstream element
                        !% thus, it is mapped to the dummy element
                        faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem
                        faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter
                        faceI(FaceLocalCounter,fi_BCtype)   = BCup

                        !% change the node assignmebt value
                        nAssignStatus =  nAssigned
                        firstUpBcHandeled = .true.
                    endif

                !% Handle 2 branch junction nodes
                case (nJ2)
                    !% Check 2: If the node has already been assigned
                    if (nAssignStatus .eq. nUnassigned) then
                        !% integer data
                        faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                        faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter
                        faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist

                        !% Check 3: If the node is an edge node (meaning this node is the
                        !% connecting node across partitions)
                        if (nodeI(thisNode,ni_P_is_boundary) .eq. EdgeNode) then

                            !% An upstream edge node indicates there are no local
                            !% elements upstream of that node. Thus it is mapped to 
                            !% the dummy element
                            faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem

                            !% logical data
                            faceYN(FaceLocalCounter,fYN_isSharedFace) = .true.
                            faceYN(FaceLocalCounter,fYN_isUpGhost)    = .true.

                            !% find the connecting image to this face
                            linkUp  => nodeI(thisNode,ni_Mlink_u1)

                            faceI(FaceLocalCounter, fi_Connected_image) = linkI(linkUp,li_P_image)
                        else
                            !% if the node is not an edge node, there is an upstream
                            !% local link element which has already been handeled
                            faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                            faceI(FaceLocalCounter,fi_Melem_uL) = ElemLocalCounter - oneI

                            faceYN(FaceLocalCounter,fYN_isInteriorFace) = .true.
                        endif

                        !% change the node assignmebt value
                        nAssignStatus =  nAssigned
                    endif

                !% Handle junction nodes with more than 2 branches (multi branch junction node).
                case (nJm)

                    !% Check 2: If the node has already been assigned
                    if (nAssignStatus .eq. nUnassigned) then

                        !% multibranch junction nodes will have both elements and faces.
                        !% thus, a seperate subroutine is required to handle these nodes
                        call init_network_handle_nJm &
                            (image, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                            FaceGlobalCounter, nAssignStatus)
                    endif

                !% Handle storage nodes
                case (nStorage)

                    print*, 'In ', subroutine_name 
                    print*, 'error: storage node is not handeled yet'

                case default

                    print*, 'In ', subroutine_name
                    print*, 'error: unexpected node, ', nodeType,'  at upstream boundary'
                    stop
            end select

        !% handle the node if it is not in the partition
        else

            !% integer data
            faceI(FacelocalCounter,fi_Lidx)     = FacelocalCounter
            faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
            faceI(FacelocalCounter,fi_Connected_image) = nodeI(thisNode,ni_P_image)


            !% if the upstream node is not in the partiton,
            !% the face map to upstream is mapped to 
            !% the dummy element
            faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem
            !% since no upstream node indicates start of a partiton,
            !% the downstream element will be initialized elem idx
            faceI(FacelocalCounter,fi_Melem_dL) = ElemLocalCounter

            !% since this will be a shared face, the global counter will be set from
            !% init_network_map_shared_faces subroutine
            faceI(FacelocalCounter,fi_Gidx) = nullvalueI

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
        endif

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name
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
        real(8)                 :: zCenter
        integer, pointer        :: lAssignStatus, NlinkElem
        real(8), pointer        :: zUpstream

        character(64) :: subroutine_name = 'init_network_handle_link'

    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% necessary pointers
        lAssignStatus => linkI(thisLink,li_assigned)

        if (lAssignStatus .eq. lUnAssigned) then
            NlinkElem     => linkI(thisLink,li_N_element)
            zUpstream     => nodeR(upNode,nr_Zbottom)

            !% store the ID of the first (upstream) element in this link
            linkI(thisLink,li_first_elem_idx)   = ElemLocalCounter
            !% reference elevations at cell center
            zCenter = zUpstream - 0.5 * linkR(thisLink,lr_ElementLength) * linkR(thisLink,lr_Slope)

            !% HACK CODE
            !% the node handler usually dont advance face counter.
            !% thus, if an assigned junction node in the same partition
            !% upstream of the link is encountered, the face counters
            !% needed to be advanced (this will be the upstream face of 
            !% the link element). This is necessary because the do loop
            !% here will always advance the downstream face of a link element

            !% the mapping of this face will be carried out later
            if ( (nodeI(upNode,ni_P_image)   .eq. image    )  .and. &
                 (nodeI(upNode,ni_assigned)  .eq. nAssigned)  .and. &
                 (nodeI(upNode,ni_node_type) .eq. nJm      ) ) then

                !% advancing face counter
                FaceLocalCounter  = FaceLocalCounter  + oneI
                FaceGlobalCounter = FaceGlobalCounter + oneI

                !% face integer data
                faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
                faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter

            endif

            do ii = 1, NlinkElem
                !%................................................................
                !% Element arrays update
                !%................................................................

                !% integer data
                elemI(ElemLocalCounter,ei_Lidx)             = ElemLocalCounter
                elemI(ElemLocalCounter,ei_Gidx)             = ElemGlobalCounter
                elemI(ElemLocalCounter,ei_elementType)      = linkI(thisLink,li_link_type)
                elemI(ElemLocalCounter,ei_link_Gidx_SWMM)   = thisLink
                elemI(ElemLocalCounter,ei_Mface_uL)         = FaceLocalCounter
                elemI(ElemLocalCounter,ei_Mface_dL)         = FaceLocalCounter + oneI
                elemI(ElemLocalCounter,ei_link_pos)         = ii

                !% real data
                elemR(ElemLocalCounter,er_Length)           = linkR(thisLink,lr_ElementLength)
                elemR(ElemLocalCounter,er_Zbottom)          = zCenter

                !%................................................................
                !% Face arrays update
                !%................................................................

                !% advance the downstream face counter of a link element
                FaceLocalCounter  = FaceLocalCounter  + oneI
                FaceGlobalCounter = FaceGlobalCounter + oneI

                !% face integer data
                faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter + oneI
                faceI(FaceLocalCounter,fi_Melem_uL) = ElemLocalCounter
                faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist

                !% counter for element z bottom calculation
                zCenter = zCenter - linkR(thisLink,lr_ElementLength) * linkR(thisLink,lr_Slope)

                !% Advance the element counter
                ElemLocalCounter  = ElemLocalCounter  + oneI
                ElemGlobalCounter = ElemGlobalCounter + oneI
            end do

            lAssignStatus = lAssigned
            linkI(thisLink,li_last_elem_idx)    = ElemLocalCounter - oneI

        endif

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name
    end subroutine init_network_handle_link
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_network_handle_downstreamnode &
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

        character(64) :: subroutine_name = 'init_network_handle_downstreamnode'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% Check 1: Is the node is in the partition
        if (nodeI(thisNode,ni_P_image) .eq. image) then

            !% necessary pointers
            nAssignStatus => nodeI(thisNode,ni_assigned)
            nodeType      => nodeI(thisNode,ni_node_type)

            select case (nodeType)

                case(nBCdn)
                    !% Check 2: If the node has already been assigned
                    if (nAssignStatus .eq. nUnassigned) then
                        !% by the time we reach a downstream boundary
                        !% node, all the indexes have already been set
                        !% from the subdivide_link_going_downstream.
                        !% only the map to downstream element is
                        !% needed to be fixed.

                        !% integer data
                        !% an downstream boundary face does not have any local downstream element
                        !% thus, it is mapped to the dummy element
                        faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem
                        faceI(FaceLocalCounter,fi_BCtype)   = BCdn

                        !% change the node assignmebt value
                        nAssignStatus =  nAssigned
                    endif

                case (nJ2)
                    !% Check 2: If the node has already been assigned
                    if (nAssignStatus .eq. nUnassigned) then
                        !% by the time we reach a downstream nJ2 node,
                        !% all the face indexes have already been set from
                        !% subdivide_link_going_downstream. only the map
                        !% to downstream element is needed to be fixed if
                        !% the node is an edge node.
                        faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist

                        !% Check 3: If the node is an edge node (meaning this node is the
                        !% connecting node across partitions)
                        if (nodeI(thisNode,ni_P_is_boundary) .eq. EdgeNode) then

                            !% An downstream edge node indicates there are no local
                            !% elements downstream of that node. thus, it is mapped to the dummy element
                            faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem

                            !% logical data
                            faceYN(FaceLocalCounter,fYN_isSharedFace) = .true.
                            faceYN(FaceLocalCounter,fYN_isDnGhost)    = .true.

                            !% find the connecting image to this face
                            linkDn  => nodeI(thisNode,ni_Mlink_d1)

                            faceI(FaceLocalCounter, fi_Connected_image) = linkI(linkDn,li_P_image)
                        endif

                        !% change the node assignmebt value
                        nAssignStatus =  nAssigned
                    endif

                case (nJm)
                    !% Check 2: If the node has already been assigned
                    if (nAssignStatus .eq. nUnassigned) then

                        call init_network_handle_nJm &
                            (image, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                            FaceGlobalCounter, nAssignStatus)

                    endif

                case (nStorage)

                    print*, 'In ', subroutine_name 
                    print*, 'error: storage node is not handeled yet'

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
            !% Thus, setting the map elem ds to dummy elem
            !% integer data
            faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem
            faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
            !% since this will be a shared face, the global counter will be set from
            !% init_network_map_shared_faces subroutine
            faceI(FacelocalCounter,fi_Gidx)     = nullvalueI
            faceI(FacelocalCounter,fi_Connected_image) = nodeI(thisNode,ni_P_image)

            !% logical data
            faceYN(FacelocalCounter,fYN_isSharedFace) = .true.
            faceYN(FaceLocalCounter,fYN_isDnGhost)    = .true.
        endif

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name
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

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !%................................................................
        !% Junction main
        !%................................................................
        !% Element Arrays
        !% integer data
        elemI(ElemLocalCounter,ei_Lidx)             = ElemLocalCounter
        elemI(ElemLocalCounter,ei_Gidx)             = ElemGlobalCounter
        elemI(ElemLocalCounter,ei_elementType)      = JM
        elemI(ElemLocalCounter,ei_node_Gidx_SWMM)   = thisNode

        !% real data
        elemR(ElemLocalCounter,er_Zbottom) = nodeR(thisNode,nr_zbottom)

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
            elemR(ElemLocalCounter,er_Zbottom) = nodeR(thisNode,nr_zbottom)
            elemR(ElemLocalCounter,er_Depth)   = nodeR(thisNode,nr_InitialDepth)

            !%......................................................
            !% Upstream Branches
            !%......................................................
            ! if ((ii .eq. 1) .or. (ii .eq. 3) .or. (ii .eq. 5)) then
            select case (mod(ii,2))   
            case (1)    
            !% finds odd number branches    
            !% all the odd numbers are upstream branches
                upBranchSelector = upBranchSelector + oneI
                !% pointer to upstream branch
                upBranchIdx => nodeI(thisNode,ni_idx_base1 + upBranchSelector)

                !% real branches
                if (upBranchIdx .ne. nullvalueI) then
                    !% integer data
                    elemSI(ElemLocalCounter,eSI_JunctionBranch_Exists)           = oneI
                    elemSI(ElemLocalCounter,eSI_JunctionBranch_Link_Connection)  = upBranchIdx
                    elemR(ElemLocalCounter,er_Length) = init_network_nJm_branch_length(upBranchIdx)

                    !% Check 4: if the link connecting this branch is a part of this partition and 
                    !% the node is not an edge node (meaning this node is the connecting node 
                    !% across partitions)
                    if ( (nodeI(thisNode,ni_P_is_boundary) .eq. EdgeNode)  .and. &
                         (linkI(upBranchIdx,li_P_image)    .ne. image   ) )  then

                        !% faces are always advanced by link elements unless it is
                        !% a null-branch or the branch is in a different image.

                        !% advance the face counters for next branch
                        FaceLocalCounter  = FaceLocalCounter  + oneI
                        FaceGlobalCounter = FaceGlobalCounter + oneI

                        !% elem array
                        !% map the face if the branch is in a different image
                        elemI(ElemLocalCounter,ei_Mface_uL) = FaceLocalCounter

                        !% face array
                        !% integer data
                        faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                        faceI(FacelocalCounter,fi_Gidx)     = FaceGlobalCounter
                        faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter
                        !% since this is a shared face in the upstream direction,
                        !% the up_map is set to dummy element 
                        faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem
                        faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
                        faceI(FaceLocalCounter,fi_Connected_image) = linkI(upBranchIdx,li_P_image)

                        !% logical data
                        faceYN(FacelocalCounter,fYN_isSharedFace) = .true.
                        faceYN(FaceLocalCounter,fYN_isUpGhost)    = .true.
                    endif
                else
                    
                !% null branches require a valid face row
                !% face counters are always advanced by link
                !% elements unless a branch does not have any
                !% elements associated with it. For null branches,
                !% no links are connected to the junction branch,
                !% but a face row is defined and connected.

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
                    !% since this is a null face in the upstream direction,
                    !% the up_map is set to dummy element 
                    faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem
                    faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist

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
                dnBranchIdx => nodeI(thisNode,ni_idx_base2 + dnBranchSelector)

                !% Check 3: if the branch is a valid branch
                if (dnBranchIdx .ne. nullvalueI) then
                    !% integer data
                    elemSI(ElemLocalCounter,eSI_JunctionBranch_Exists)          = oneI
                    elemSI(ElemLocalCounter,eSI_JunctionBranch_Link_Connection) = dnBranchIdx
                    elemR(ElemLocalCounter,er_Length) = init_network_nJm_branch_length(dnBranchIdx)

                    !% Check 4: if the link connecting this branch is a part of this partition and 
                    !% the node is not an edge node (meaning this node is the connecting node 
                    !% across partitions)
                    if ( (nodeI(thisNode,ni_P_is_boundary) .eq. EdgeNode)  .and. &
                         (linkI(dnBranchIdx,li_P_image)    .ne. image   ) )  then

                        !% faces are always advanced by link elements unless it is
                        !% a null-branch or the branch is in a different image.

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
                        !% since this is a shared face in the downstream direction,
                        !% the dn_map is set to dummy element 
                        faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem
                        faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
                        faceI(FaceLocalCounter,fi_Connected_image) = linkI(dnBranchIdx,li_P_image)

                        !% logical data
                        faceYN(FacelocalCounter,fYN_isSharedFace) = .true.
                        faceYN(FaceLocalCounter,fYN_isDnGhost)    = .true.
                    endif

                else

                    !% face counters are always advanced for null branches
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
                    !% since this is a null face in the downstream direction,
                    !% the dn_map is set to dummy element 
                    faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem
                    faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist

                    call init_network_nullify_nJm_branch &
                        (ElemLocalCounter, FaceLocalCounter)
                endif
            end select !% case (mod(ii,2))  

            !% Advance the element counter for next branch
            ElemLocalCounter  = ElemLocalCounter  + oneI
            ElemGlobalCounter = ElemGlobalCounter + oneI
        end do

        !% set status to assigned
        nAssignStatus = nAssigned

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name
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

                    !% find the first element index of the link
                    LinkFirstElem = linkI(dnBranchIdx,li_first_elem_idx)

                    !% find the downstream face index of that last element
                    fLidx => elemI(LinkFirstElem,ei_Mface_uL)

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
        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

        !% find the length of the junction branch
        if (linkI(LinkIdx,li_length_adjusted) .eq. OneSideAdjust) then

            BranchLength = linkR(LinkIdx,lr_Length) - linkR(LinkIdx,lr_AdjustedLength)

        elseif (linkI(LinkIdx,li_length_adjusted) .eq. BothSideAdjust) then

            BranchLength = (linkR(LinkIdx,lr_Length) - linkR(LinkIdx,lr_AdjustedLength))/twoR
        endif

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name
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
        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

        !% set everything to zero for a non existant branch
        elemR(ElemIdx,:)                            = zeroR
        elemSR(ElemIdx,:)                           = zeroR
        elemSI(ElemIdx,eSI_JunctionBranch_Exists)   = zeroI
        faceR(FaceIdx,:)                            = zeroR
        faceYN(FaceIdx,fYN_isnull)                  = .true.

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name
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
    !
        character(64) :: subroutine_name = 'init_network_set_interior_faceYN'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        where ( (faceI(:,fi_BCtype)         .eq.  doesnotexist) &
                .and. &
                (faceYN(:,fYN_isnull)       .eqv. .false.     ) &
                .and. &
                (faceYN(:,fYN_isSharedFace) .eqv. .false.     ) )

            faceYN(:,fYN_isInteriorFace) = .true.
        endwhere

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name
    end subroutine init_network_set_interior_faceYN
    !
    !==========================================================================
    ! END OF MODULE
    !==========================================================================
    !
end module network_define

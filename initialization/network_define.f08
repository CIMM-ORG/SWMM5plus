!
! module network_define
!
! Handles relationship between coarse link-node network and high-resolution
! element-face network.
!
! This module defines all the indexes and mappings
!
!==========================================================================
!
module network_define
    !
    use allocate_storage
    use assign_index
    use initialization
    use data_keys
    use globals
    use interface
    use setting_definition

    implicit none

    private

    public :: network_initiation

contains
    !
    !==========================================================================
    ! PUBLIC
    !==========================================================================
    !
    subroutine network_initiation ()
    !--------------------------------------------------------------------------
    !
    !% Initializes a element-face network from a link-node network.
    !%   Requires network links and nodes before execution 
    !
    !--------------------------------------------------------------------------

        character(64) :: subroutine_name = 'network_initiation'

    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

    
        !% get the slope of each link given the node Z values
        call network_get_link_slope()

        !% divide the link node networks in elements and faces 
        call network_data_create()

        sync all 
        !% print result
        if (setting%Debug%File%network_define) then
           ! image = this_image()
           if (this_image() == 1) then

            do image = 1,3
               print*, '----------------------------------------------------'
               print*, 'image = ', image
               print*, '..................elements..........................'
               print*, elemI(:,ei_Lidx)[image], 'Lidx'
               print*, elemI(:,ei_Gidx)[image], 'Gidx'
               print*, elemI(:,ei_elementType)[image], 'elementType'
               print*, elemI(:,ei_geometryType)[image], 'geometryType'
               print*, elemI(:,ei_link_Gidx_SWMM)[image], 'link_Gidx_SWMM'
               print*, elemI(:,ei_node_Gidx_SWMM)[image], 'node_Gidx_SWMM'
               print*, elemI(:,ei_Mface_uL)[image],'Mface_uL'
               print*, elemI(:,ei_Mface_dL)[image],'Mface_dL'
               print*, '..................faces.............................'
               print*, faceI(:,fi_Lidx)[image], 'face Lidx'
               print*, faceI(:,fi_Gidx)[image], 'face Gidx'
               print*, faceI(:,fi_Melem_dL)[image], 'face Melem_dL'
               print*, faceI(:,fi_Melem_uL)[image], 'face Melem_uL'
               print*, faceYN(:,fYN_isSharedFace)[image], 'face is shared face'
               print*, '----------------------------------------------------'
               call execute_command_line('')
           enddo
            endif
        endif

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine network_initiation
    !
    !==========================================================================
    ! PRIVATE
    !==========================================================================
    !
    subroutine network_get_link_slope
    !-------------------------------------------------------------------------- 
    !
    !% compute the slope across each link
    !
    !--------------------------------------------------------------------------

        character(64) :: subroutine_name = 'network_get_link_slope'

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

    end subroutine network_get_link_slope
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine network_data_create()
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
        integer, dimension(:), allocatable, target :: pack_link_idx, pack_node_idx

        character(64) :: subroutine_name = 'network_data_create'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% Number of Images
        !% HACK: This is only works on every single image, so we will need to change when we move to a control and computation node structure.
        N_Images => setting%Partitioning%N_Image

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
              print*, '--------------------------------------'
              print*, image, 'image'
              print*, ii, 'ii'
              print*, FaceGlobalIdx, 'FaceGlobalIdx'
              print*, '--------------------------------------'
           end do
        end if        

        sync all

        !% handle all the links and nodes in a partition
        call handle_partition &
            (image, ElemLocalIdx, FacelocalIdx, ElemGlobalIdx, FaceGlobalIdx)

        sync all

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine network_data_create
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine handle_partition &
        (image, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, FaceGlobalCounter)
    !    
    !-------------------------------------------------------------------------- 
    !
    !% Handle all the links in a partition
    !
    !--------------------------------------------------------------------------
    ! 
        integer, intent(in)    :: image
        integer, intent(inout) :: ElemLocalCounter, FaceLocalCounter
        integer, intent(inout) :: ElemGlobalCounter, FaceGlobalCounter

        integer :: ii, pLink, pNode
        integer, pointer :: thisLink, upNode, dnNode
        integer, dimension(:), allocatable, target :: packed_link_idx, packed_node_idx



        character(64) :: subroutine_name = 'handle_partition'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% pack all the link indexes in an image
        packed_link_idx = pack(linkI(:,li_idx), (linkI(:,li_BQ_image) == image))

        !% find the number of links in an image
        pLink = size(packed_link_idx)

        !% pack all the node indexes in a partition to determine which nodes are in the partition
        packed_node_idx = pack(nodeI(:,ni_idx), (nodeI(:,ni_BQ_image) == image))
        
        !% number of nodes in a partition
        pNode = size(packed_node_idx)

        !% cycle through the links in an image
        do ii = 1,pLink
            !% necessary pointers to the links and connected nodes
            thisLink => packed_link_idx(ii)
            upNode   => linkI(thisLink,li_Mnode_u)
            dnNode   => linkI(thisLink,li_Mnode_d)

            call handle_upstream_node &
                (image, upNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                FaceGlobalCounter)

            call subdivide_link_going_downstream &
                (thisLink, upNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                FaceGlobalCounter)

            call handle_downstream_node &
                (image, dnNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                FaceGlobalCounter)
        end do

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine handle_partition
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine handle_upstream_node &
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

        integer, pointer :: nAssignStatus, nodeType

        integer :: ii
        
        character(64) :: subroutine_name = 'handle_upstream_node'
    !--------------------------------------------------------------------------
        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% check 1: Is the node is in the partition
        if (nodeI(thisNode,ni_BQ_image) .eq. image) then

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

                        if (nodeI(thisNode,ni_BQ_edge) .eq. EdgeNode) then

                            !% An upstream edge node indicates there are no local
                            !% elements upstream of that node
                            faceI(FaceLocalCounter,fi_Melem_uL) = nullvalueI

                            !% logical data
                            faceYN(FaceLocalCounter,fYN_isSharedFace) = .true.
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

                    print*, 'In ', subroutine_name
                    print*, 'Upstream multi branch junction node will be coded later'
                    stop

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

            !% since no upstream node indicates start of a partiton,
            !% the downstream element will be initialized elem idx
            faceI(FacelocalCounter,fi_Melem_dL) = ElemLocalCounter

            !% since the first face is a shared face, the global
            !% index will be updated later.
            !% However, subdivide_link_going_downstream will advance
            !% the global face counter without considering this update 
            !% later. Thus, this index should be adjusted by subtracting 
            !% one from the FaceGlobalCounter

            FaceGlobalCounter = FaceGlobalCounter - oneI
            
            !% logical data
            faceYN(FacelocalCounter,fYN_isSharedFace) = .true.
        endif
    
        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine handle_upstream_node
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine subdivide_link_going_downstream &
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

        character(64) :: subroutine_name = 'subdivide_link_going_downstream'

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

    end subroutine subdivide_link_going_downstream
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine handle_downstream_node &
        (image, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
        FaceGlobalCounter)
    !-------------------------------------------------------------------------- 
    !
    !% handle the node dwonstream of a link
    !
    !--------------------------------------------------------------------------
    ! 
        integer, intent(in)    :: image, thisNode
        integer, intent(inout) :: ElemLocalCounter, FaceLocalCounter 
        integer, intent(inout) :: ElemGlobalCounter, FaceGlobalCounter

        integer, pointer :: nAssignStatus, nodeType
        integer, pointer :: upBranchIdx, dnBranchIdx

        integer :: ii, upBranchSelector, dnBranchSelector
        
        character(64) :: subroutine_name = 'handle_downstream_node'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% check 1: Is the node is in the partition
        if (nodeI(thisNode,ni_BQ_image) .eq. image) then

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
                        if (nodeI(thisNode,ni_BQ_edge) .eq. EdgeNode) then

                            !% An downstream edge node indicates there are no local
                            !% elements downstream of that node
                            faceI(FaceLocalCounter,fi_Melem_dL) = nullvalueI

                            !% logical data
                            faceYN(FaceLocalCounter,fYN_isSharedFace) = .true.
                        endif

                        !% change the node assignmebt value
                        nAssignStatus =  nAssigned
                    endif

                case (nJm)
                    !% check if the node has already been assigned
                    if (nAssignStatus .eq. nUnassigned) then

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
                                    if ( (nodeI(thisNode,ni_BQ_edge) .eq. EdgeNode)  .and. &
                                         (linkI(upBranchIdx,li_BQ_image) .ne. image) )  then

                                        !% HACK:
                                        !% faces are always advanced by link elements
                                        !% however, if there aren't any link element
                                        !% connected to the junction branch in a partition,
                                        !% the face is advanced

                                        !% advance the face counters for next branch
                                        FaceLocalCounter  = FaceLocalCounter  + oneI
                                        FaceGlobalCounter = FaceGlobalCounter + oneI

                                        !% integer data
                                        faceI(FaceLocalCounter,fi_Lidx) = FaceLocalCounter
                                        faceI(FacelocalCounter,fi_Gidx) = FaceGlobalCounter

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

                                    !% integer data
                                    faceI(FaceLocalCounter,fi_Lidx) = FaceLocalCounter
                                    faceI(FacelocalCounter,fi_Gidx) = FaceGlobalCounter

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
                                    if ( (nodeI(thisNode,ni_BQ_edge) .eq. EdgeNode)  .and. &
                                         (linkI(dnBranchIdx,li_BQ_image) .ne. image) )  then

                                        !% HACK:
                                        !% faces are always advanced by link elements
                                        !% however, if there aren't any link element
                                        !% connected to the junction branch in a partition,
                                        !% the face is advances

                                        !% advance the face counters for next branch
                                        FaceLocalCounter  = FaceLocalCounter  + oneI
                                        FaceGlobalCounter = FaceGlobalCounter + oneI

                                        !% integer data
                                        faceI(FaceLocalCounter,fi_Lidx) = FaceLocalCounter
                                        faceI(FacelocalCounter,fi_Gidx) = FaceGlobalCounter

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

                                    !% integer data
                                    faceI(FaceLocalCounter,fi_Lidx) = FaceLocalCounter
                                    faceI(FacelocalCounter,fi_Gidx) = FaceGlobalCounter

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

            !% logical data
            faceYN(FacelocalCounter,fYN_isSharedFace) = .true.
        endif

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine handle_downstream_node
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
    !==========================================================================
    !
    !
    !==========================================================================
    ! END OF MODULE
    !==========================================================================
    !
end module network_define

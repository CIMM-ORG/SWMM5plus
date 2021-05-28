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
               print*, faceYN(:,fYN_isnull)[image], 'face is_null'
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
        
        integer :: ii, image, pLink, pNode
        integer :: ElemGlobalIdx, FaceGlobalIdx 
        integer :: ElemLocallIdx, FacelocallIdx
        integer :: lastElem, lastFace, P_elem, P_face
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
        ElemLocallIdx = first_elem_index
        FacelocallIdx = first_face_index

        !% Setting the local image value
        image = this_image()

        !% initalizing the global element and face id by looping through N_elem and N_face for images not equal to one
        
        if(this_image() /= 1) then

           do ii=1, this_image()-1
              ElemGlobalIdx = ElemGlobalIdx + N_elem(ii)
              !% we have to subtract one from the global face id such that faces along the image boundaries are shared.
              FaceGlobalIdx = FaceGlobalIdx + N_face(ii)-1
           end do
        end if        

        !% pack all the link indexes in a partition to cycle through the links
        pack_link_idx = pack(linkI(:,li_idx), (linkI(:,li_BQ_image) == image))

        !% number of links in a partition
        pLink = size(pack_link_idx)

        !% pack all the node indexes in a partition to determine which nodes are in the partition
        pack_node_idx = pack(nodeI(:,ni_idx), (nodeI(:,ni_BQ_image) == image))
        
        !% number of nodes in a partition
        pNode = size(pack_node_idx)

        do ii = 1, pLink
           !% cycle through link indexs in a partition 
           Lidx      => pack_link_idx(ii)
           NodeUp    => linkI(Lidx,li_Mnode_u)
           NodeDn    => linkI(Lidx,li_Mnode_d)
           NodeUpTyp => nodeI(NodeUp,ni_node_type)
           NodeDnTyp => nodeI(NodeDn,ni_node_type)
           
           sync all

            !% Search if the link in a partition has a upstream node
            if (any(pack_node_idx .eq. NodeUp)) then

                call handle_upstream_node &
                    (NodeUp, NodeUpTyp, ElemLocallIdx, FacelocallIdx, &
                    ElemGlobalIdx, FaceGlobalIdx)
            endif

            call subdivide_link_going_downstream &
                    (Lidx, ElemLocallIdx, FacelocallIdx, ElemGlobalIdx, FaceGlobalIdx)

           if (any(pack_node_idx .eq. NodeDn)) then
              call handle_downstream_node &
                    (NodeDn, NodeDnTyp, ElemLocallIdx, FacelocallIdx, &
                    ElemGlobalIdx, FaceGlobalIdx)
           endif

           !% at this point we have the information of which elemI array idxs
           !% are connected to which links/nodes. Now mapping out the faces
           !% by cycling through the nodes should be simple

        enddo
        sync all

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine network_data_create
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine handle_upstream_node &
        (NodeIdx, NodeTyp, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, FaceGlobalCounter)
    !-------------------------------------------------------------------------- 
    !
    !% handle the node upstream of a link
    !
    !--------------------------------------------------------------------------

        integer, pointer       :: nAssignStatus, usBranch, dsBranch
        integer, intent(in)    :: NodeIdx, NodeTyp
        integer, intent(inout) :: ElemLocalCounter, FaceLocalCounter
        integer, intent(inout) :: ElemGlobalCounter, FaceGlobalCounter
        integer                :: ii, jj, kk

        character(64) :: subroutine_name = 'handle_upstream_node'

    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% node assign status
        nAssignStatus => nodeI(NodeIdx,ni_assigned)

        select case (NodeTyp)

            case (nBCup)
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

                    !% advance the face counters
                    FaceLocalCounter = FaceLocalCounter + oneI
                    FaceGlobalCounter = FaceGlobalCounter + oneI

                endif

            case (nJ2)
                !% check if the node has already been assigned
                if (nAssignStatus .eq. nUnassigned) then
                    !% integer data
                    faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                    faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                    !%-----------------------------------------------------------
                    !% HACK: we still dont know the upstream element to this face
                    !%-----------------------------------------------------------
                    faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter

                    !% change the node assignmebt value
                    nAssignStatus =  nAssigned

                    !% advance the face counters
                    FaceLocalCounter = FaceLocalCounter + oneI
                    FaceGlobalCounter = FaceGlobalCounter + oneI 
                else
                    !% figure out a way to map the upstream element
                    stop
                endif

            case (nJm)
                !% check if the node has already been assigned
                if (nAssignStatus .eq. nUnassigned) then
                    !%-----------------------------------------------------------------
                    !% Junction main 
                    !% integer data
                    elemI(ElemLocalCounter,ei_Lidx) = ElemLocalCounter
                    elemI(ElemLocalCounter,ei_Gidx) = ElemGlobalCounter
                    elemI(ElemLocalCounter,ei_elementType) = eJunctionMain
                    elemI(ElemLocalCounter,ei_node_Gidx_SWMM) = NodeIdx

                    !% real data
                    elemR(ElemLocalCounter,er_Zbottom) = nodeR(NodeIdx,nr_zbottom)

                    !% Advance the element counter to 1st upstream branch
                    ElemLocalCounter = ElemLocalCounter + oneI
                    ElemGlobalCounter = ElemGlobalCounter + oneI

                    !%-----------------------------------------------------------------
                    !% Handling all the junction branches
                    !%-----------------------------------------------------------------
                    jj = 0
                    kk = 0
                    do ii = 1,max_branch_per_node
                        !% condition for upstrem branches
                        if ((ii .eq. 1) .or. (ii .eq. 2) .or. (ii .eq. 3)) then
                            jj = jj + 1
                            usBranch => nodeI(NodeIdx,ni_idx_base1 + jj)

                            !%-----------------------------------------------------------------
                            !% Junction Upstream Branch 
                            !%-----------------------------------------------------------------
                            elemI(ElemLocalCounter,ei_Lidx)         = ElemLocalCounter
                            elemI(ElemLocalCounter,ei_Gidx)         = ElemGlobalCounter
                            elemI(ElemLocalCounter,ei_elementType)  = eJunctionBranch
                            elemI(ElemLocalCounter,ei_node_Gidx_SWMM) = NodeIdx

                            !% face data
                            faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                            faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                            faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter
                            
                            !% right now, we dont know any information about the upstream
                            !% element of this branch
                            if (usBranch .ne. nullvalueI) then
                                !% integer data
                                !% HACK: Not sure about this. Needs furtehr testing
                                elemI(ElemLocalCounter,ei_Mface_uL)                 = FaceLocalCounter
                                elemSI(ElemLocalCounter,eSI_JunctionBranch_Exists)  = oneI

                                !% real data
                                elemR(ElemLocalCounter,er_Zbottom) = nodeR(NodeIdx,nr_zbottom)

                                !% find the length of the junction branch
                                if (linkI(usBranch,li_length_adjusted) .eq. OneSideAdjust) then
                                    elemR(ElemLocalCounter,er_Length) = linkR(usBranch,lr_Length) - &
                                                linkR(usBranch,lr_AdjustedLength)
                                elseif (linkI(usBranch,li_length_adjusted) .eq. BothSideAdjust) then
                                    elemR(ElemLocalCounter,er_Length) = (linkR(usBranch,lr_Length) - &
                                                linkR(usBranch,lr_AdjustedLength))/twoR
                                else
                                    print*, 'error in, ', subroutine_name
                                    print*, 'link connected to junction has not been shortened'
                                    stop
                                endif

                            else
                                !% set everything to zero for a non existant branch
                                elemR(ElemLocalCounter,:) = zeroR
                                elemSR(ElemLocalCounter,:) = zeroR
                                elemSI(ElemLocalCounter,eSI_JunctionBranch_Exists) = zeroI
                                faceR(FaceLocalCounter,:) = zeroR
                                faceYN(FaceLocalCounter,fYN_isnull) = .true.
                            endif

                            !% Advance the element counter to 2nd downstream branch
                            ElemLocalCounter = ElemLocalCounter + oneI
                            ElemGlobalCounter = ElemGlobalCounter + oneI 

                            !% advance the face counters
                            FaceLocalCounter = FaceLocalCounter + oneI
                            FaceGlobalCounter = FaceGlobalCounter + oneI

                        !% condition for downstream branch
                        else      
                            kk = kk + 1
                            dsBranch => nodeI(NodeIdx,ni_idx_base2 + kk)

                            !%-----------------------------------------------------------------
                            !% Junction Downstream Branch
                            !%-----------------------------------------------------------------
                            elemI(ElemLocalCounter,ei_Lidx) = ElemLocalCounter
                            elemI(ElemLocalCounter,ei_Gidx) = ElemGlobalCounter
                            elemI(ElemLocalCounter,ei_elementType) = eJunctionBranch
                            elemI(ElemLocalCounter,ei_node_Gidx_SWMM) = NodeIdx

                            !% face data
                            faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                            faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                            faceI(FaceLocalCounter,fi_Melem_uL) = ElemLocalCounter

                            !% right now, we dont know any information about the downstream
                            !% element of this branch
                            if (dsBranch .ne. nullvalueI) then
                                !% integer data
                                !% HACK: Not sure about this. Needs furtehr testing
                                elemI(ElemLocalCounter,ei_Mface_dL)                 = FaceLocalCounter
                                elemSI(ElemLocalCounter,eSI_JunctionBranch_Exists)  = oneI
                                !% real data
                                elemR(ElemLocalCounter,er_Zbottom) = nodeR(NodeIdx,nr_zbottom)
                                !% find the length of the junction branch
                                if (linkI(dsBranch,li_length_adjusted) .eq. OneSideAdjust) then
                                    elemR(ElemLocalCounter,er_Length) = linkR(dsBranch,lr_Length) - &
                                            linkR(dsBranch,lr_AdjustedLength)
                                elseif (linkI(dsBranch,li_length_adjusted) .eq. BothSideAdjust) then
                                    elemR(ElemLocalCounter,er_Length) = (linkR(dsBranch,lr_Length) - &
                                            linkR(dsBranch,lr_AdjustedLength))/twoR
                                else
                                    print*, 'error in, ', subroutine_name
                                    print*, 'link connected to junction has not been shortened'
                                    stop
                                endif

                            else
                                !% set everything to zero for a non existant branch
                                elemR(ElemLocalCounter,:) = zeroR
                                elemSR(ElemLocalCounter,:) = zeroR
                                elemSI(ElemLocalCounter,eSI_JunctionBranch_Exists) = zeroI
                                faceR(FaceLocalCounter,:) = zeroR
                                faceYN(FaceLocalCounter,fYN_isnull) = .true.
                            endif

                            !% change the node assignmebt value
                            nAssignStatus =  nAssigned

                            !% Advance the element counter to 3rd upstream branch
                            ElemLocalCounter = ElemLocalCounter + oneI
                            ElemGlobalCounter = ElemGlobalCounter + oneI 

                            !% advance the face counters
                            FaceLocalCounter = FaceLocalCounter + oneI
                            FaceGlobalCounter = FaceGlobalCounter + oneI

                        endif
                    enddo
                                                  
                endif

            case default
                print*, 'error: unexpected node, ', NodeTyp,'  at upstream boundary'
                stop
        end select

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine handle_upstream_node
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine subdivide_link_going_downstream &
        (Lidx, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, FaceGlobalCounter)
    !-------------------------------------------------------------------------- 
    !
    !% subdivides the links into elements going downstream
    !
    !--------------------------------------------------------------------------

        integer, pointer        :: NlinkElem, lAssignStatus, NodeUp
        integer, intent(in)     :: Lidx
        integer, intent(inout)  :: ElemLocalCounter, FaceLocalCounter
        integer, intent(inout)  :: ElemGlobalCounter, FaceGlobalCounter
        integer                 :: ii

        real                    :: zUpstream, zCenter

        character(64) :: subroutine_name = 'subdivide_link_going_downstream'

    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% necessary pointers
        NlinkElem  => linkI(Lidx,li_N_element)
        ! LinkTyp       => linkI(Lidx,li_link_type)
        nodeUp     => linkI(Lidx,li_Mnode_u)
        zUpstream  =  nodeR(NodeUp,nr_Zbottom)

        !% link assign status
        lAssignStatus => linkI(Lidx,li_assigned)

        !%  store the ID of the first (upstream) element in this link
        linkI(Lidx,li_Melem_u) = ElemLocalCounter - oneI ! the new element is the 1st (downstream)
        linkI(Lidx,li_Mface_u) = FaceLocalCounter - oneI ! the old face is the 1st

        !%  reference elevations at cell center and cell face
        zCenter = zUpstream - 0.5 * linkR(Lidx,lr_ElementLength) * linkR(Lidx,lr_Slope)

        do ii = 1,NlinkElem
            !% integer data
            elemI(ElemLocalCounter,ei_Lidx)             = ElemLocalCounter
            elemI(ElemLocalCounter,ei_Gidx)             = ElemGlobalCounter
            elemI(ElemLocalCounter,ei_geometryType)     = linkI(Lidx,li_geometry)
            elemI(ElemLocalCounter,ei_elementType)      = linkI(Lidx,li_link_type)
            elemI(ElemLocalCounter,ei_link_Gidx_SWMM)   = Lidx
            elemI(ElemLocalCounter,ei_Mface_uL)         = FaceLocalCounter - oneI
            elemI(ElemLocalCounter,ei_Mface_dL)         = FaceLocalCounter

            !% real data
            elemR(ElemLocalCounter,er_Length) = linkR(Lidx,lr_ElementLength)
            elemR(ElemLocalCounter,er_Zcrown) = zCenter

            !% face integer data
            faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
            faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
            faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter + oneI
            faceI(FaceLocalCounter,fi_Melem_uL) = ElemLocalCounter
            

            !% face logical data
            faceYN(FaceLocalCounter,fYN_isnull) = .false.

            !% counter for element z bottom calculation
            zCenter = zCenter - linkR(Lidx,lr_ElementLength) * linkR(Lidx,lr_Slope)

            !% Advance the element counter to 3rd upstream branch
            ElemLocalCounter = ElemLocalCounter + oneI
            ElemGlobalCounter = ElemGlobalCounter + oneI 

            !% advance the face counters
            FaceLocalCounter = FaceLocalCounter + oneI
            FaceGlobalCounter = FaceGlobalCounter + oneI
        end do

        lAssignStatus = lAssigned

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine subdivide_link_going_downstream
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine handle_downstream_node &
        (NodeIdx, NodeTyp, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, FaceGlobalCounter)
    !-------------------------------------------------------------------------- 
    !
    !% handle the node dwonstream of a link
    !
    !--------------------------------------------------------------------------

        integer, pointer       :: nAssignStatus, usBranch, dsBranch
        integer, intent(in)    :: NodeIdx, NodeTyp
        integer, intent(inout) :: ElemLocalCounter, FaceLocalCounter
        integer, intent(inout) :: ElemGlobalCounter, FaceGlobalCounter
        integer                :: ii, jj, kk

        character(64) :: subroutine_name = 'handle_downstream_node'

    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        nAssignStatus => nodeI(NodeIdx,ni_assigned)

        select case (NodeTyp)

            case (nBCdn)
                !% check if the node has already been assigned
                if (nAssignStatus .eq. nUnassigned) then
                    !% integer data
                    faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                    faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                    faceI(FaceLocalCounter,fi_Melem_uL) = ElemLocalCounter - oneI
                    faceI(FaceLocalCounter,fi_Melem_dL) = nullvalueI
                    faceI(FaceLocalCounter,fi_BCtype)   = fBCup

                    !% change the node assignmebt value
                    nAssignStatus =  nAssigned

                    !% advance the face counters
                    FaceLocalCounter  = FaceLocalCounter + oneI
                    FaceGlobalCounter = FaceGlobalCounter + oneI 
                endif

            case (nJ2)
                !% check if the node has already been assigned
                if (nAssignStatus .eq. nUnassigned) then
                    !% integer data
                    faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                    faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                    !%-----------------------------------------------------------
                    !% HACK: we still dont know the upstream element to this face
                    !%-----------------------------------------------------------
                    faceI(FaceLocalCounter,fi_Melem_uL) = ElemLocalCounter - oneI

                    !% change the node assignmebt value
                    nAssignStatus =  nAssigned

                    !% advance the face counters
                    FaceLocalCounter = FaceLocalCounter + oneI
                    FaceGlobalCounter = FaceGlobalCounter + oneI 
                else
                    !% figure out a way to map the upstream element
                    stop
                endif

            case (nJm)
                !% check if the node has already been assigned
                if (nAssignStatus .eq. nUnassigned) then
                    !%-----------------------------------------------------------------
                    !% Junction main 
                    !% integer data
                    elemI(ElemLocalCounter,ei_Lidx) = ElemLocalCounter
                    elemI(ElemLocalCounter,ei_Gidx) = ElemGlobalCounter
                    elemI(ElemLocalCounter,ei_elementType) = eJunctionMain
                    elemI(ElemLocalCounter,ei_node_Gidx_SWMM) = NodeIdx

                    !% real data
                    elemR(ElemLocalCounter,er_Zbottom) = nodeR(NodeIdx,nr_zbottom)

                    !% Advance the element counter to 1st upstream branch
                    ElemLocalCounter = ElemLocalCounter + oneI
                    ElemGlobalCounter = ElemGlobalCounter + oneI

                    !%-----------------------------------------------------------------
                    !% Handling all the junction branches
                    !%-----------------------------------------------------------------
                    jj = 0
                    kk = 0
                    do ii = 1,max_branch_per_node
                        !% condition for upstrem branches
                        if ((ii .eq. 1) .or. (ii .eq. 3) .or. (ii .eq. 5)) then
                            jj = jj + 1
                            print*, 'jj', jj
                            print*, 'BRANCH', usBranch
                            usBranch => nodeI(NodeIdx,ni_idx_base1 + jj)

                            !%-----------------------------------------------------------------
                            !% Junction Upstream Branch 
                            !%-----------------------------------------------------------------
                            elemI(ElemLocalCounter,ei_Lidx)           = ElemLocalCounter
                            elemI(ElemLocalCounter,ei_Gidx)           = ElemGlobalCounter
                            elemI(ElemLocalCounter,ei_elementType)    = eJunctionBranch
                            elemI(ElemLocalCounter,ei_node_Gidx_SWMM) = NodeIdx
                            elemI(ElemLocalCounter,ei_Mface_uL)       = ElemLocalCounter - oneI

                            !% face data
                            faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                            faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                            faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter
                            
                            !% right now, we dont know any information about the upstream
                            !% element of this branch
                            if (usBranch .ne. nullvalueI) then
                                
                                !% integer data
                                !% HACK: Not sure about this. Needs furtehr testing
                                
                                elemSI(ElemLocalCounter,eSI_JunctionBranch_Exists)  = oneI

                                !% real data
                                elemR(ElemLocalCounter,er_Zbottom) = nodeR(NodeIdx,nr_zbottom)

                                !% find the length of the junction branch
                                if (linkI(usBranch,li_length_adjusted) .eq. OneSideAdjust) then
                                    elemR(ElemLocalCounter,er_Length) = linkR(usBranch,lr_Length) - &
                                                linkR(usBranch,lr_AdjustedLength)
                                elseif (linkI(usBranch,li_length_adjusted) .eq. BothSideAdjust) then
                                    elemR(ElemLocalCounter,er_Length) = (linkR(usBranch,lr_Length) - &
                                                linkR(usBranch,lr_AdjustedLength))/twoR
                                else
                                    print*, 'error in, ', subroutine_name
                                    print*, 'link connected to junction has not been shortened'
                                endif

                            else
                                !% set everything to zero for a non existant branch
                                elemR(ElemLocalCounter,:) = zeroR
                                elemSR(ElemLocalCounter,:) = zeroR
                                elemSI(ElemLocalCounter,eSI_JunctionBranch_Exists) = zeroI
                                faceR(FaceLocalCounter,:) = zeroR
                                faceYN(FaceLocalCounter,fYN_isnull) = .true.
                            endif

                            !% Advance the element counter to 2nd downstream branch
                            ElemLocalCounter = ElemLocalCounter + oneI
                            ElemGlobalCounter = ElemGlobalCounter + oneI 

                            !% advance the face counters
                            FaceLocalCounter = FaceLocalCounter + oneI
                            FaceGlobalCounter = FaceGlobalCounter + oneI

                        else 
                            !% condition for downstream branch
                            kk = kk + 1
                            dsBranch => nodeI(NodeIdx,ni_idx_base2 + kk)

                            !%-----------------------------------------------------------------
                            !% Junction Downstream Branch
                            !%-----------------------------------------------------------------
                            elemI(ElemLocalCounter,ei_Lidx)             = ElemLocalCounter
                            elemI(ElemLocalCounter,ei_Gidx)             = ElemGlobalCounter
                            elemI(ElemLocalCounter,ei_elementType)      = eJunctionBranch
                            elemI(ElemLocalCounter,ei_node_Gidx_SWMM)   = NodeIdx
                            elemI(ElemLocalCounter,ei_Mface_dL)         = ElemLocalCounter - oneI

                            !% face data
                            faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                            faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                            faceI(FaceLocalCounter,fi_Melem_uL) = ElemLocalCounter

                            !% right now, we dont know any information about the downstream
                            !% element of this branch
                            if (dsBranch .ne. nullvalueI) then
                                !% integer data
                                !% HACK: Not sure about this. Needs furtehr testing
                                elemSI(ElemLocalCounter,eSI_JunctionBranch_Exists)  = oneI
                                !% real data
                                elemR(ElemLocalCounter,er_Zbottom) = nodeR(NodeIdx,nr_zbottom)
                                !% find the length of the junction branch
                                if (linkI(dsBranch,li_length_adjusted) .eq. OneSideAdjust) then
                                    elemR(ElemLocalCounter,er_Length) = linkR(dsBranch,lr_Length) - &
                                                linkR(dsBranch,lr_AdjustedLength)
                                elseif (linkI(dsBranch,li_length_adjusted) .eq. BothSideAdjust) then
                                    elemR(ElemLocalCounter,er_Length) = (linkR(dsBranch,lr_Length) - &
                                                linkR(dsBranch,lr_AdjustedLength))/twoR
                                else
                                    print*, 'error in, ', subroutine_name
                                    print*, 'link connected to junction has not been shortened'
                                endif

                            else
                                !% set everything to zero for a non existant branch
                                elemR(ElemLocalCounter,:) = zeroR
                                elemSR(ElemLocalCounter,:) = zeroR
                                elemSI(ElemLocalCounter,eSI_JunctionBranch_Exists) = zeroI
                                faceR(FaceLocalCounter,:) = zeroR
                                faceYN(FaceLocalCounter,fYN_isnull) = .true.
                            endif

                            !% Advance the element counter to 3rd upstream branch
                            ElemLocalCounter = ElemLocalCounter + oneI
                            ElemGlobalCounter = ElemGlobalCounter + oneI 

                            !% advance the face counters
                            FaceLocalCounter = FaceLocalCounter + oneI
                            FaceGlobalCounter = FaceGlobalCounter + oneI
                        endif
                    enddo
                                                  
                endif

            case default
                print*, 'error: unexpected node, ', NodeTyp,'  at upstream boundary'
                stop
        end select

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine handle_downstream_node
    !
    !==========================================================================
    ! END OF MODULE
    !==========================================================================
    !
end module network_define

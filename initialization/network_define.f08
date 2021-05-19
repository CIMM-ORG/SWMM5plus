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

        !% print result
        if (setting%Debug%File%network_define) then
            do image = 1, num_images()
                print*, '----------------------------------------------------'
                print*, 'image = ', image
                print*, '..................elements..........................'
                print*, elemI(:,ei_Lidx)[image], 'ei_Lidx'
                print*, elemI(:,ei_Gidx)[image], 'ei_Gidx'
                print*, elemI(:,ei_elementType)[image], 'ei_elementType'
                print*, elemI(:,ei_geometryType)[image], 'ei_geometryType'
                print*, elemI(:,ei_link_Gidx_SWMM)[image], 'ei_link_Gidx_SWMM'
                print*, elemI(:,ei_node_Gidx_SWMM)[image], 'ei_node_Gidx_SWMM'
                print*, '..................faces.............................'
                print*, faceI(:,fi_Lidx)[image], 'fi_Lidx'
                print*, faceI(:,fi_Gidx)[image], 'fi_Gidx'
                print*, '----------------------------------------------------'  
            enddo
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
        
        integer :: ii, jj, image, firstIdx
        integer :: ElemGlobalIdx, FaceGlobalIdx 
        integer :: ElemLocallIdx, FacelocallIdx
        integer :: lastElem, lastFace 
        integer, pointer :: N_Images, P_elem, P_face, Lidx
        integer, pointer :: NodeUp, NodeDn, NodeUpTyp, NodeDnTyp
        integer, dimension(:), allocatable, target :: pack_link_idx, pack_node_idx

        character(64) :: subroutine_name = 'network_data_create'
    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        !% Number of Images
        !% HACK: This is hardcoded for now. We have to come up with a way 
        !% to use it as a user defined parameters
        N_Images => setting%Partitioning%N_Image

        !% initializing global element and face id
        ElemGlobalIdx = first_elem_index
        FaceGlobalIdx = first_face_index

        !% distribute the network across images
        do image = 1, num_images()

            !% initializing local element and face id
            ElemLocallIdx = first_elem_index
            FacelocallIdx = first_face_index
            !% initializing the first element number of link elements in a partition
            firstIdx = oneI

            !% pointer to number of elements and faces to this image
            P_elem => N_elem(image)
            P_face => N_face(image)


            elemI(ElemLocallIdx:P_elem, ei_Lidx)[image] = [(jj,jj=ElemLocallIdx,P_elem)]
            elemI(ElemLocallIdx:P_elem, ei_Gidx)[image] = [(jj,jj=ElemGlobalIdx,P_elem + ElemGlobalIdx - 1)]
            faceI(ElemLocallIdx:P_face, fi_Lidx)[image] = [(jj,jj=FacelocallIdx,P_face)]
            faceI(ElemLocallIdx:P_face, fi_Gidx)[image] = [(jj,jj=FaceGlobalIdx,P_face + FaceGlobalIdx -1)]

            ElemGlobalIdx = ElemGlobalIdx + P_elem
            FaceGlobalIdx = FaceGlobalIdx + P_face

            !% pack all the link indexes in a partition to cycle through the links
            pack_link_idx = pack(linkI(:,li_idx), (linkI(:,li_BQ_image) == image))
            !% pack all the node indexes in a partition to determine which nodes are in the partition
            pack_node_idx = pack(nodeI(:,ni_idx), (nodeI(:,ni_BQ_image) == image))

            do ii = 1, size(pack_link_idx)
                !% cycle through link indexs in a partition 
                Lidx    => pack_link_idx(ii)
                NodeUp  => linkI(Lidx,li_Mnode_u)
                NodeDn  => linkI(Lidx,li_Mnode_d)
                NodeUpTyp => nodeI(NodeUp,ni_node_type)
                NodeDnTyp => nodeI(NodeDn,ni_node_type)

                !% THIS IS A HACK CODE:
                !% the central idea can be more organized
                !% First populate the elemI array
                !% only links and multiface junction will contribute to elem arrays

                !% condition of a specific node to be in that partition that is junction
                if (any(pack_node_idx .eq. NodeUp) .and. (NodeUpTyp .eq.  nJm)) then

                    call handle_multi_branch_node (NodeUp, image, firstIdx)

                    call subdivide_link (Lidx, image, firstIdx)

                else

                    call subdivide_link (Lidx, image, firstIdx)
                endif


                if (any(pack_node_idx .eq. NodeDn) .and. (NodeDnTyp .eq.  nJm)) then

                    call handle_multi_branch_node (NodeDn, image, firstIdx)

                endif

                !% at this point we have the information of which elemI array idxs
                !% are connected to which links/nodes. Now mapping out the faces
                !% by cycling through the nodes should be simple
           
            enddo

        enddo

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine network_data_create
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine subdivide_link (Lidx, image, firstIdx)
    !-------------------------------------------------------------------------- 
    !
    !% subdivides the links into elements
    !
    !--------------------------------------------------------------------------

        integer, pointer       :: LinkElem, lAssignStatus
        integer, intent(in)    :: Lidx, image
        integer, intent(inout) :: firstIdx
        integer                :: lastIdx

        character(64) :: subroutine_name = 'subdivide_link'

    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        LinkElem     => linkI(Lidx,li_N_element)
        !% if the link is already assigned or not
        lAssignStatus => linkI(Lidx,li_assigned)
        lastIdx = firstIdx + LinkElem - oneI

        ! print*, Lidx, 'Lidx'
        ! print*, LinkElem, 'LinkElem'
        ! print*, lAssignStatus, 'lAssignStatus'
        ! print*, firstIdx,'firstIdx'
        ! print*, lastIdx, 'lastIdx'
        ! print*, image, 'image'

        if (image .eq. oneI) then
            if (lAssignStatus .eq. lUnassigned) then
                elemI(firstIdx:lastIdx,ei_elementType)    = linkI(Lidx,li_link_type)   
                elemI(firstIdx:lastIdx,ei_geometryType)   = linkI(Lidx,li_geometry)
                elemI(firstIdx:lastIdx,ei_link_Gidx_SWMM) = Lidx
                elemI(firstIdx:lastIdx,ei_node_Gidx_SWMM) = nullvalueI

                lAssignStatus =  lAssigned

                !% set the first element index for next link
                firstIdx = firstIdx + LinkElem
            endif
        else
            if (lAssignStatus .eq. lUnassigned) then
                !% populating elemI
                elemI(firstIdx:lastIdx,ei_elementType)[image]    = linkI(Lidx,li_link_type)   
                elemI(firstIdx:lastIdx,ei_geometryType)[image]   = linkI(Lidx,li_geometry)
                elemI(firstIdx:lastIdx,ei_link_Gidx_SWMM)[image] = Lidx
                elemI(firstIdx:lastIdx,ei_node_Gidx_SWMM)[image] = nullvalueI

                lAssignStatus =  lAssigned

                !% set the first element index for next link
                firstIdx = firstIdx + LinkElem
            endif

        endif

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine subdivide_link
    ! !
    ! !==========================================================================
    ! !==========================================================================
    ! !
    subroutine handle_multi_branch_node (Nidx, image, firstIdx)
    !-------------------------------------------------------------------------- 
    !
    !% subdivides the links into elements
    !
    !--------------------------------------------------------------------------

        integer, pointer       :: NodeElem, nAssignStatus
        integer, intent(in)    :: Nidx, image
        integer, intent(inout) :: firstIdx
        integer                :: lastIdx

        character(64) :: subroutine_name = 'handle_multi_branch_node'

    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        NodeElem      => J_elem_add
        !% if the link is already assigned or not
        nAssignStatus => nodeI(Nidx,ni_assigned)
        lastIdx = firstIdx + NodeElem - oneI


        if (image .eq. oneI) then
            if (nAssignStatus .eq. nUnassigned) then
                elemI(firstIdx:lastIdx,ei_elementType)    = nodeI(Nidx,ni_node_type)   
                elemI(firstIdx:lastIdx,ei_link_Gidx_SWMM) = nullvalueI
                elemI(firstIdx:lastIdx,ei_node_Gidx_SWMM) = Nidx

                nAssignStatus =  nAssigned

                !% set the first element index for next link
                firstIdx = firstIdx + NodeElem
            endif
        else
            if (nAssignStatus .eq. nUnassigned) then
                !% populating elemI
                elemI(firstIdx:lastIdx,ei_elementType)[image]    = nodeI(Nidx,ni_node_type)   
                elemI(firstIdx:lastIdx,ei_link_Gidx_SWMM)[image] = nullvalueI
                elemI(firstIdx:lastIdx,ei_node_Gidx_SWMM)[image] = Nidx

                nAssignStatus =  nAssigned

                !% set the first element index for next link
                firstIdx = firstIdx + NodeElem
            endif

        endif

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine handle_multi_branch_node
 
    !
    !==========================================================================
    ! END OF MODULE
    !==========================================================================
    !
end module network_define
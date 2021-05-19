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
        
        integer :: ii, jj, image
        integer :: ElemGlobalIdx, FaceGlobalIdx 
        integer :: ElemLocallIdx, FacelocallIdx
        integer :: lastElem, lastFace 
        integer, pointer :: N_Images
        integer, dimension(:), allocatable :: pack_links

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
                        
        do image = 1,  N_Images

            !% initializing local element and face id
            ElemLocallIdx = first_elem_index
            FacelocallIdx = first_face_index
            lastElem      = zeroR
            lastFace      = zeroR

            !% pack all the link indexes in a partition
            pack_links = pack(linkI(:,li_idx), (linkI(:,li_BQ_image) == image)
            
            do ii = size(pack_links)
                Lidx    => pack_links(ii)
                
                ! NodeDn  => linkI(Lidx,li_Mnode_d)

                !% handle the upstream node of that link
                call handle_upstream_node &
                    (Nidx, image, ElemLocalIdx, ElemGlobalIdx, FaceLocalIdx, FaceGlobalIdx)

                !% subdivide the link according to number of elements
                call subdivide_link (Lidx, image, ElemLocalIdx, &
                        ElemGlobalIdx, FaceLocalIdx, FaceGlobalIdx)

                !% handle the downstream node of that link
                call handle_this_node (NodeDn, image, ElemLocalIdx, &
                        ElemGlobalIdx, FaceLocalIdx, FaceGlobalIdx)


            enddo
        enddo

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine network_data_create
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine subdivide_link(Lidx, image, ElemLocalIdx, &
        ElemGlobalIdx, FaceLocalIdx, FaceGlobalIdx)
    !-------------------------------------------------------------------------- 
    !
    !% subdivides the links into elements
    !
    !--------------------------------------------------------------------------
        integer :: ii, LastElemIdx, LastFaceIdx

        integer, pointer       :: AssignStatus
        integer, intent(in)    :: Lidx
        integer, intent(in)    :: image
        integer, intent(inout) :: ElemLocalIdx, ElemGlobalIdx
        integer, intent(inout) :: FaceLocalIdx, FaceGlobalIdx


        character(64) :: subroutine_name = 'network_data_create'

    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        AssignStatus    => linkI(Lidx, li_assigned)

        if (AssignStatus .eq. lUnassigned) then

            do ii = 1, linkI(Lidx,li_N_element) 

                !% populating elemI
                elem2I(ElemLocalIdx,ei_Lidx)[image]           = ElemLocalIdx
                elem2I(ElemLocalIdx,ei_Gidx)[image]           = ElemGlobalIdx
                elem2I(ElemLocalIdx,ei_elementType)[image]    = linkI(Lidx,li_link_type)
                elem2I(ElemLocalIdx,ei_geometryType)[image]   = linkI(Lidx,li_geometry)
                elem2I(ElemLocalIdx,ei_link_Gidx_SWMM)[image] = Lidx
                elem2I(ElemLocalIdx,ei_Mface_dL)[image]       = FaceLocalIdx
                elem2I(ElemLocalIdx,ei_Mface_uL)[image]       = FaceLocalIdx - oneI

                !% populating faceI
                !% upstream face of an element will always be populated
                !% thus advancing face counter
                FaceLocalIdx  = FaceLocalIdx + oneI
                FaceGlobalIdx = FaceGlobalIdx + oneI

                faceI(FaceLocalIdx,fi_Lidx)                   = FaceLocalIdx
                faceI(FaceLocalIdx,fi_Gidx)                   = ElemGlobalIdx
                faceI(FaceLocalIdx,fi_Melem_uL)               = ElemLocalIdx
                faceI(FaceLocalIdx,fi_Melem_dL)               = ElemLocalIdx + oneI

                ElemLocalIdx  = ElemLocalIdx + oneI
                ElemGlobalIdx = ElemGlobalIdx + oneI


            end do

            AssignStatus = setAssigned(Lidx, li_assigned, lUnassigned, lAssigned, linkI)

        endif

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine subdivide_link
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine handle_upstream_node &
        (Nidx, image, ElemLocalIdx, ElemGlobalIdx, FaceLocalIdx, FaceGlobalIdx)

    !-------------------------------------------------------------------------- 
    !
    !% handels an upstream node of a specific link
    !
    !--------------------------------------------------------------------------
        
        integer, intent(in)    :: Lidx 
        integer, intent(in)    :: current_image
        integer, intent(inout) :: ElemLocalIdx, ElemGlobalIdx
        integer, intent(inout) :: FaceLocalIdx, FaceGlobalIdx

        integer, pointer :: Nidx, NodeType, AssignStatus, NJelem
        integer, pointer :: LinkU1, LinkU2, LinkU3, LinkD1, LinkD2, LinkD3

        integer :: ii, JMLocalIdx, JMGlobalIdx

        character(64) :: subroutine_name = 'handle_upstream_node'

    !--------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name

        Nidx         => linkI(Lidx,li_Mnode_u)
        NodeType     => nodeI(Nidx,ni_node_type)
        AssignStatus => nodeI(Nidx,ni_assigned)
        LinkU1       => nodeI(Nidx,ni_Mlink_u1)
        LinkU2       => nodeI(Nidx,ni_Mlink_u2)
        LinkU3       => nodeI(Nidx,ni_Mlink_u3)
        LinkD1       => nodeI(Nidx,ni_Mlink_d1)
        LinkD2       => nodeI(Nidx,ni_Mlink_d2)
        LinkD3       => nodeI(Nidx,ni_Mlink_d3)

        select case (NodeType)

            case (nBCup)

                if (AssignStatus .eq. nUnassigned) then

                    faceI(FaceLocalIdx, fi_Lidx)[image]     = FaceLocalIdx  !% update local face index in faceI
                    faceI(FaceLocalIdx, fi_Gidx)[image]     = FaceGlobalIdx !% update global face index in faceI
                    faceI(FaceLocalIdx, fi_Melem_uL)[image] = nullvalueI    !% upstream element map of nBCup is nullvalue
                    faceI(FaceLocalIdx, fi_Melem_dL)[image] = ElemLocalIdx

                    AssignStatus = setAssigned(Nidx, ni_assigned, nUnassigned, nAssigned, nodeI)

                endif

            case (nJ2)

                if (AssignStatus .eq. nUnassigned) then

                    faceI(FaceLocalIdx, fi_Lidx)[image]     = FaceLocalIdx        !% update local face index in faceI
                    faceI(FaceLocalIdx, fi_Gidx)[image]     = FaceGlobalIdx       !% update global face index in faceI
                    faceI(FaceLocalIdx, fi_Melem_uL)[image] = ElemLocalIdx  
                    faceI(FaceLocalIdx, fi_Melem_dL)[image] = ElemLocalIdx + oneI !% downstream element map of nBCup is nullvalue

                    AssignStatus = setAssigned(Nidx, ni_assigned, nUnassigned, nAssigned, nodeI)
                endif

            case (nJm)

                if (AssignStatus .eq. nUnassigned) then

                    NJelem => Nelem_in_Junction

                    do ii = 1, NJelem

                        if (ii .eq. Junction_main) then

                            elem2I(ElemLocalIdx,ei_Lidx)[image]           = ElemLocalIdx
                            elem2I(ElemLocalIdx,ei_Gidx)[image]           = ElemGlobalIdx
                            elem2I(ElemLocalIdx,ei_elementType)[image]    = eJunctionMain
                            elem2I(ElemLocalIdx,ei_node_Gidx_SWMM)[image] = Nidx
                            elem2I(ElemLocalIdx,ei_Mface_dL)[image]       = nullvalueI
                            elem2I(ElemLocalIdx,ei_Mface_uL)[image]       = nullvalueI

                            ElemLocalIdx  = ElemLocalIdx + oneI
                            ElemGlobalIdx = ElemGlobalIdx + oneI

                        elseif (ii .eq. Junction_branch_1_in) then

                            elem2I(ElemLocalIdx,ei_Lidx)[image]           = ElemLocalIdx
                            elem2I(ElemLocalIdx,ei_Gidx)[image]           = ElemGlobalIdx
                            elem2I(ElemLocalIdx,ei_elementType)[image]    = eJunctionBranch
                            elem2I(ElemLocalIdx,ei_node_Gidx_SWMM)[image] = Nidx
                            !% an upstream branch doesnot have any donwstream face
                            elem2I(ElemLocalIdx,ei_Mface_dL)[image]       = nullvalueI
                            elem2I(ElemLocalIdx,ei_Mface_uL)[image]       = FaceLocalIdx - oneI



                            !% populating faceI
                            faceI(FaceLocalIdx,fi_Lidx)                   = FaceLocalIdx
                            faceI(FaceLocalIdx,fi_Gidx)                   = ElemGlobalIdx
                            faceI(FaceLocalIdx,fi_Melem_uL)               = ElemLocalIdx
                            faceI(FaceLocalIdx,fi_Melem_dL)               = ElemLocalIdx + oneI

                            ElemLocalIdx  = ElemLocalIdx + oneI
                            ElemGlobalIdx = ElemGlobalIdx + oneI

                            FaceLocalIdx  = FaceLocalIdx + oneI
                            ElemGlobalIdx = ElemGlobalIdx + oneI

                        elseif (ii .eq. Junction_branch_2_out) then

                            elem2I(ElemLocalIdx,ei_Lidx)[image]           = ElemLocalIdx
                            elem2I(ElemLocalIdx,ei_Gidx)[image]           = ElemGlobalIdx
                            elem2I(ElemLocalIdx,ei_elementType)[image]    = eJunctionBranch
                            elem2I(ElemLocalIdx,ei_node_Gidx_SWMM)[image] = Nidx
                            elem2I(ElemLocalIdx,ei_Mface_dL)[image]       = FaceLocalIdx - oneI
                            ! an downstream branch doesnot have any upstream face
                            elem2I(ElemLocalIdx,ei_Mface_uL)[image]       = nullvalueI

                            ElemLocalIdx  = ElemLocalIdx + oneI
                            ElemGlobalIdx = ElemGlobalIdx + oneI

                            FaceLocalIdx  = FaceLocalIdx + oneI
                            ElemGlobalIdx = ElemGlobalIdx + oneI

                        elseif (ii .eq. Junction_branch_3_in) then

                            elem2I(ElemLocalIdx,ei_Lidx)[image]           = ElemLocalIdx
                            elem2I(ElemLocalIdx,ei_Gidx)[image]           = ElemGlobalIdx
                            elem2I(ElemLocalIdx,ei_elementType)[image]    = eJunctionBranch
                            elem2I(ElemLocalIdx,ei_node_Gidx_SWMM)[image] = Nidx
                            !% an upstream branch doesnot have any donwstream face
                            elem2I(ElemLocalIdx,ei_Mface_dL)[image]       = nullvalueI
                            elem2I(ElemLocalIdx,ei_Mface_uL)[image]       = FaceLocalIdx - oneI

                            ElemLocalIdx  = ElemLocalIdx + oneI
                            ElemGlobalIdx = ElemGlobalIdx + oneI

                            FaceLocalIdx  = FaceLocalIdx + oneI
                            ElemGlobalIdx = ElemGlobalIdx + oneI

                        elseif (ii .eq. Junction_branch_4_out) then

                            elem2I(ElemLocalIdx,ei_Lidx)[image]           = ElemLocalIdx
                            elem2I(ElemLocalIdx,ei_Gidx)[image]           = ElemGlobalIdx
                            elem2I(ElemLocalIdx,ei_elementType)[image]    = eJunctionBranch
                            elem2I(ElemLocalIdx,ei_node_Gidx_SWMM)[image] = Nidx
                            elem2I(ElemLocalIdx,ei_Mface_dL)[image]       = FaceLocalIdx - oneI
                            ! an downstream branch doesnot have any upstream face
                            elem2I(ElemLocalIdx,ei_Mface_uL)[image]       = nullvalueI

                            ElemLocalIdx  = ElemLocalIdx + oneI
                            ElemGlobalIdx = ElemGlobalIdx + oneI

                            FaceLocalIdx  = FaceLocalIdx + oneI
                            ElemGlobalIdx = ElemGlobalIdx + oneI

                        elseif (ii .eq. Junction_branch_5_in) then

                            elem2I(ElemLocalIdx,ei_Lidx)[image]           = ElemLocalIdx
                            elem2I(ElemLocalIdx,ei_Gidx)[image]           = ElemGlobalIdx
                            elem2I(ElemLocalIdx,ei_elementType)[image]    = eJunctionBranch
                            elem2I(ElemLocalIdx,ei_node_Gidx_SWMM)[image] = Nidx
                            !% an upstream branch doesnot have any donwstream face
                            elem2I(ElemLocalIdx,ei_Mface_dL)[image]       = nullvalueI
                            elem2I(ElemLocalIdx,ei_Mface_uL)[image]       = FaceLocalIdx - oneI

                            ElemLocalIdx  = ElemLocalIdx + oneI
                            ElemGlobalIdx = ElemGlobalIdx + oneI

                            FaceLocalIdx  = FaceLocalIdx + oneI
                            ElemGlobalIdx = ElemGlobalIdx + oneI

                        elseif (ii .eq. Junction_branch_6_out) then

                            elem2I(ElemLocalIdx,ei_Lidx)[image]           = ElemLocalIdx
                            elem2I(ElemLocalIdx,ei_Gidx)[image]           = ElemGlobalIdx
                            elem2I(ElemLocalIdx,ei_elementType)[image]    = eJunctionBranch
                            elem2I(ElemLocalIdx,ei_node_Gidx_SWMM)[image] = Nidx
                            elem2I(ElemLocalIdx,ei_Mface_dL)[image]       = FaceLocalIdx - oneI
                            ! an downstream branch doesnot have any upstream face
                            elem2I(ElemLocalIdx,ei_Mface_uL)[image]       = nullvalueI

                            ElemLocalIdx  = ElemLocalIdx + oneI
                            ElemGlobalIdx = ElemGlobalIdx + oneI

                            FaceLocalIdx  = FaceLocalIdx + oneI
                            ElemGlobalIdx = ElemGlobalIdx + oneI

                        else
                            print*, 'error: unexpected numer of junction branch encountered'
                            stop
                        endif


                    end do

                    AssignStatus = setAssigned(Nidx, ni_assigned, nUnassigned, nAssigned, nodeI)

                endif

            case default
                print *, 'error: unknown node type ',NodeType,' in ',subroutine_name
            stop
        end select


        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine handle_upstream_node
    !
    !==========================================================================
    !==========================================================================
    !
    !
    !==========================================================================
    !==========================================================================
    !
    function setAssigned &
        (thisX, assigned_idx, unassignedValue, assignedValue, arrayI) &
        result(f_result)
        !
        ! Change the assigned status or stop on error if the node or link is
        ! already assigned.
        !
        character(64) :: subroutine_name = 'setAssigned'
        integer :: f_result
        integer, intent(in) :: thisX, assigned_idx, unassignedValue, assignedValue
        integer, intent(in) :: arrayI(:,:)
        !
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        if (arrayI(thisX,assigned_idx) .eq. unassignedValue) then
            f_result = assignedValue
        else
            print *, thisX, assigned_idx, unassignedValue, assignedValue
            print *, arrayI(thisX,assigned_idx)
            print *, 'error: expected an unassigned node or link in ',subroutine_name
            STOP
        endif

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end function setAssigned
    !
    !==========================================================================
    ! END OF MODULE
    !==========================================================================
    !
end module network_define
module network_define
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Handles relationship between coarse link-node network and high-resolution
    !% element-face network. This module defines all the indexes and mappings
    !%
    !%==========================================================================
    !
    use interface_
    use utility_allocate
    use discretization
    use define_indexes
    use define_keys
    use define_globals
    use define_settings
    use utility_profiler
    use utility_crash, only: util_crashpoint

    implicit none
    private

    public :: network_define_toplevel

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine network_define_toplevel ()
        !%------------------------------------------------------------------
        !% Description:
        !% Initializes a element-face network from a link-node network.
        !%   Requires network links and nodes before execution
        !%-------------------------------------------------------------------
        !% Declarations:
            integer :: ii, jj
            character(64) :: subroutine_name = 'init_network_define_toplevel'
        !%-------------------------------------------------------------------
        !% Preliminaries:
            !if (crashYN) return
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            if (setting%Profile%useYN) call util_profiler_start (pfc_init_network_define_toplevel)
        !%------------------------------------------------------------------

        !% --- get the slope of each link given the node Z values
        call init_network_linkslope ()

        !% --- divide the link node networks in elements and faces
        call init_network_datacreate ()

        !% --- replace ni_elemface_idx of nJ2 nodes for the upstream elem
        !%     of the face associated with the node
        call init_network_update_nj2_elem ()

        sync all

        !%------------------------------------------------------------------
        !% Closing

        if (setting%Debug%File%network_define) then
            print *
            print *, '===================================================================' //&
            '==================================================================='
            print *, 'image = ', this_image()
            print *, '.......................Elements...............................'
            ! print *
            ! print *, 'a)   ei_Lidx       ei_Gidx     link_BQ    link_SWMM  node_BQ   node_SWMM'
            ! do jj = 1,N_elem(this_image())
            !     print *, elemI(jj,ei_Lidx), elemI(jj,ei_Gidx), elemI(jj,ei_link_Gidx_BIPquick), &
            !     elemI(jj,ei_link_Gidx_SWMM),elemI(jj,ei_node_Gidx_BIPquick),elemI(jj,ei_node_Gidx_SWMM)
            ! end do
            print *
            print *, 'b)   ei_Lidx      ei_Type      Mface_uL    Mface_dL   '
            do jj = 1,N_elem(this_image())
                print*, elemI(jj,ei_Lidx), reverseKey(elemI(jj,ei_elementType)), elemI(jj,ei_Mface_uL), elemI(jj,ei_Mface_dL)
            end do
            print *, '.......................Faces.............................'
            print *, 'c)     fi_Lidx     fi_Gidx    elem_uL     elem_dL   C_image    Ifidx'
            do jj = 1,N_face(this_image())
                print*, faceI(jj,fi_Lidx),faceI(jj,fi_Gidx),faceI(jj,fi_Melem_uL), &
                faceI(jj,fi_Melem_dL),faceI(jj,fi_Connected_image), faceI(jj,fi_Identical_Lidx)
            end do
            print *
            print *, 'd)     fi_Lidx     GElem_up    GElem_dn    node_SWMM    link_SWMM'
            do jj = 1,N_face(this_image())
                print *, faceI(jj,fi_Lidx),faceI(jj,fi_GhostElem_uL),&
                faceI(jj,fi_GhostElem_dL),faceI(jj, fi_node_idx_SWMM),&
                faceI(jj,fi_link_idx_SWMM)
            end do

            ! print *
            ! print *, '.......................Faces..................................'
            ! print *, 'e)     fi_Lidx     fi_BCtype   fYN_isInteriorFace    fYN_isSharedFace'//&
            ! '    fYN_isnull    fYN_isUpGhost    fYN_isDnGhost'
            ! do jj = 1,N_face(this_image())
            !     print *, faceI(jj,fi_Lidx),' ',faceI(jj,fi_BCtype),&
            !     '           ',faceYN(jj,fYN_isInteriorFace), &
            !     '                    ',faceYN(jj,fYN_isSharedFace), &
            !     '              ',faceYN(jj,fYN_isnull), &
            !     '            ',faceYN(jj,fYN_isUpGhost), &
            !     '              ',faceYN(jj,fYN_isDnGhost)
            ! end do
            print *, '===================================================================' //&
            '==================================================================='
            print *, ' '
        end if


        if (setting%Profile%useYN) call util_profiler_stop (pfc_init_network_define_toplevel)

        if (setting%Debug%File%network_define) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine network_define_toplevel
!%
!%==========================================================================
!% PRIVATE -- 1st level
!%==========================================================================
!%
    subroutine init_network_linkslope()
        !%------------------------------------------------------------------
        !% Description:
        !% compute the slope across each link
        !%------------------------------------------------------------------
        !% Declarations
            character(64) :: subroutine_name = 'init_network_linkslope'
            integer, pointer :: NodeUp, NodeDn, lType
            real(8), pointer :: zUp, zDn, Slope, Length
            integer          :: mm
        !%------------------------------------------------------------------
        !% Preliminaries:
            !if (crashYN) return
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        do mm = 1, N_link
            !% --- Inputs
            NodeUp      => link%I(mm,li_Mnode_u)
            NodeDn      => link%I(mm,li_Mnode_d)
            Length      => link%R(mm,lr_Length)
            lType       => link%I(mm,li_link_type)
            !% --- Output
            zUp         => link%R(mm,lr_ZbottomUp)
            zDn         => link%R(mm,lr_ZbottomDn)
            Slope       => link%R(mm,lr_Slope)

            if (((lType == lChannel) .or.     &
                 (lType == lPipe  )) .and.    &
                 (node%I(NodeUp,ni_node_type) == nJm)) then
                !% HACK: for upstream nJm, consider the offsets
                zUp = node%R(NodeUp,nr_Zbottom) + link%R(mm,lr_InletOffset) 
            else
                !% HACK: for upstream other nodes, ignore the offsets
                zUp = node%R(NodeUp,nr_Zbottom)
            end if

            if (((lType == lChannel) .or.     &
                 (lType == lPipe  )) .and.    &
                 (node%I(NodeDn,ni_node_type) == nJm)) then
                !% HACK: for downstream nJm, consider the offsets
                zDn = node%R(NodeDn,nr_Zbottom) + link%R(mm,lr_OutletOffset)
            else
                !% HACK: for other downstream nodes, ignore the offsets
                zDn = node%R(NodeDn,nr_Zbottom)
            end if 
            
            !%-------------------------------------------------------------------------
            !% HACK: Original length is used for slope calculation instead of adjusted
            !% length. Using adjusted lengths will induce errors in slope calculations
            !% and will distort the original network.
            !%-------------------------------------------------------------------------
            if ( (lType == lChannel) .or. (lType == lPipe)) then
                Slope = (zUp - zDn) / Length
            else
                Slope = zeroR
            end if
        end do

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%network_define) then
                !% --- provide output for debugging
                print *, subroutine_name,'--------------------------------'
                print *, 'link ID,               Slope,             length'
                do mm=1, N_link
                    print *, mm, link%R(mm,lr_Slope), link%R(mm,lr_Length)
                end do
                print *, ' '
                print *, ' node assigments '
                print *, 'link ID,         nodeUp,         nodeDn'
                do mm=1, N_link 
                    print *, ' '
                    print *, mm, link%I(mm,li_Mnode_u), link%I(mm,li_Mnode_d)
                    print *, trim(link%Names(mm)%str), ' ; ', trim(node%Names(link%I(mm,li_Mnode_u))%str), ' ; ', trim(node%Names(link%I(mm,li_Mnode_d))%str)
                end do

                print *, ' '
                do mm=1, N_node
                    print *, ' '
                    print *, mm, node%I(mm,ni_Mlink_u1), node%I(mm,ni_Mlink_d1)
                    print *, trim(node%names(mm)%str), ' ; '
                
                end do
            end if
   
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_network_linkslope
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_network_datacreate()
        !%------------------------------------------------------------------
        !% Description:
        !% creates the network of elements and faces from nodes and link
        !%
        !%------------------------------------------------------------------
        !% Declarations
            integer :: ii, image
            integer :: ElemGlobalCounter, FaceGlobalCounter
            integer :: ElemLocalCounter, FacelocalCounter
            integer :: unique_faces
            integer, pointer :: Lidx
            integer, pointer :: NodeUp, NodeDn, NodeUpTyp, NodeDnTyp
            character(64) :: subroutine_name = 'init_network_datacreate'
        !%-------------------------------------------------------------------
        !% Preliminaries    
            !if (crashYN) return
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% --- initializing global element and face index counter
        !%     these are added to through the network definition to get the index
        !%     space required to define arrays
        ElemGlobalCounter = oneI
        FaceGlobalCounter = zeroI

        !% --- initializing local element and face index countier
        ElemLocalCounter = oneI
        FacelocalCounter = zeroI

        !% --- Setting the local image value
        image = this_image()

        !% --- initialize the global indexes of elements and faces
        call init_network_set_global_indexes &
            (image, ElemGlobalCounter, FaceGlobalCounter)

        !% --- set the dummy element
        call init_network_set_dummy_elem ()

        !% --- handle all the links and nodes in a partition
        call init_network_handle_partition &
            (image, ElemLocalCounter, FacelocalCounter, ElemGlobalCounter, FaceGlobalCounter)

        !% --- finish mapping all the junction branch and faces that were not
        !%    handled in handle_link_nodes subroutine
        call init_network_map_nodes (image)

        !% --- set interior face logical
        call init_network_set_interior_faceYN ()

        !% --- shared faces are mapped by copying data from different images
        !%     thus a sync all is needed
        sync all

        !% --- set the same global face idx for shared faces across images
        call init_network_map_shared_faces (image)

        !% --- identify the boundary element connected to a shared faces
        call init_network_identify_boundary_element

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_network_datacreate
!%
!%==========================================================================
!%==========================================================================
!
    subroutine init_network_update_nj2_elem()
        !%-----------------------------------------------------------------
        !% Description:
        !% For nj2 nodes, assigns the ni_elem_idx as the element
        !% index that is upstream of the face.
        !% CAUTION: the node%I() is a global (non-coarray) but it is storing
        !% values for ni_elemface_idx that are ONLY correct on the 
        !% image ni_P_image -- on any other image you get a location that is NOT
        !% a correct element or face!
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, allocatable :: nJ2_nodes(:)
            integer              :: N_nJ2_nodes, ii
            character(64) :: subroutine_name = 'init_network_update_nj2_elem'
        !%-----------------------------------------------------------------
        !% Preliminaries
            !if (crashYN) return
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------

        N_nJ2_nodes = count((node%I(:, ni_node_type) == nJ2) .and. (node%I(:, ni_P_image) == this_image()))

        if (N_nJ2_nodes > 0) then

            nJ2_nodes = pack(node%I(:, ni_idx), (node%I(:, ni_node_type) == nJ2) &
                        .and. (node%I(:, ni_P_image) == this_image()))

            !% --- assign upstream element to an nJ2 face
            node%I(nJ2_nodes, ni_elem_idx) = faceI(node%I(nJ2_nodes, ni_face_idx),fi_Melem_uL)

        end if

        !%-----------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_network_update_nj2_elem
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_network_CC_elem_length_adjust ()
        !%--------------------------------------------------------------------------
        !%  RETAIN FOR POSSIBLE FUTURE USE
        !%--------------------------------------------------------------------------
        integer          :: ii
        integer, pointer :: AdjustType, elementType(:), elementIdx(:)
        real(8), pointer :: NominalLength, MinLengthFactor, elementLength(:)
        real(8)          :: MinElemLength
        character(64)    :: subroutine_name = 'init_network_CC_elem_length_adjust'
        !--------------------------------------------------------------------------

        print *, 'OBSOLETE 20230507 brh'
        stop 2309743

        ! AdjustType      => setting%Discretization%MinElemLengthMethod
        ! NominalLength   => setting%Discretization%NominalElemLength
        ! MinLengthFactor => setting%Discretization%MinElemLengthFactor

        ! elementIdx    => elemI(:,ei_Lidx)
        ! elementType   => elemI(:,ei_elementType)
        ! elementLength => elemR(:,er_Length)

        ! select case (AdjustType)

        ! case(RawElemLength)
        !     !% do not do any adjustment and return the raw network
        !     return

        ! case (ElemLengthAdjust)

        !     MinElemLength = NominalLength * MinLengthFactor

        !     do ii = 1,N_elem(this_image())
        !         if ((elementType(ii) == CC) .and. (elementLength(ii) < zeroR)) then
        !             print *, 'CODE ERROR: negative element length at element ',elementIdx(ii)
        !             write(*,"(A,f12.3,A,f12.3)") '       element length = ', elementLength(ii)
        !             write(*,"(A,i8)")            '       This element is in SWMM link ',elemI(elementIdx(ii),ei_link_Gidx_SWMM)
        !             write(*,"(A,A)")             '       This link name is            ',trim(link%Names(elemI(elementIdx(ii),ei_link_Gidx_SWMM))%str)
        !             write(*,"(A,f12.3)")         '       The SWMM link length is ',link%R(elemI(elementIdx(ii),ei_link_Gidx_SWMM),lr_Length)
        !             call util_crashpoint(209837)
        !         end if

        !         if ((elementType(ii) == CC) .and. (elementLength(ii) < MinElemLength)) then
        !             if (setting%Output%Verbose) then
        !                 !print*, 'In, ', subroutine_name
        !                 print *, ' '
        !                 write(*,"(A,i8,A,i5)")       '... small element detected at ElemIdx = ', elementIdx(ii), ' in processor = ',this_image()
        !                 write(*,"(A,f12.3,A,f12.3)") '       element length = ', elementLength(ii), ' is adjusted to ', MinElemLength
        !                 write(*,"(A,i8)")            '       This element is in SWMM link ',elemI(elementIdx(ii),ei_link_Gidx_SWMM)
        !                 write(*,"(A,A)")             '       This link name is            ',trim(link%Names(elemI(elementIdx(ii),ei_link_Gidx_SWMM))%str)
        !                 write(*,"(A,f12.3)")         '       The SWMM link length is ',link%R(elemI(elementIdx(ii),ei_link_Gidx_SWMM),lr_Length)

        !                 print *, 'RECOMMEND SHORTER NOMINAL ELEMENT LENGTH.'
        !                 print *, 'STOP HERE (working on new algorithm)'
        !                 call util_crashpoint(6098723)


        !             end if
        !             elementLength(ii) = MinElemLength
        !         end if
        !     end do

        ! case default
        !     print*, 'In, ', subroutine_name
        !     print *, 'CODE ERROR: AdjustType unknown for # ',AdjustType 
        !     print *, 'which has key ',trim(reverseKey(AdjustType))
        !     !stop 
        !     call util_crashpoint(89537)
        !     !return
        ! end select

        ! if (setting%Debug%File%network_define) &
        ! write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_network_CC_elem_length_adjust
!%   
!%==========================================================================    
!% PRIVATE -- 2nd Level
!%==========================================================================
!%
    subroutine init_network_set_global_indexes &
        (image, ElemGlobalCounter, FaceGlobalCounter)
        !%------------------------------------------------------------------
        !% Description:
        !% Adds the size of the elements and unique faces on each image
        !% to the global counters.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)     :: image
            integer, intent(inout)  :: ElemGlobalCounter, FaceGlobalCounter
            integer                 :: ii
            character(64) :: subroutine_name = 'init_network_set_global_indexes'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        if (image /= 1) then
           do ii=1, image-1
              ElemGlobalCounter = ElemGlobalCounter + N_elem(ii)
              FaceGlobalCounter = FaceGlobalCounter + N_unique_face(ii) + oneI
           end do
        end if

        !%-----------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%network_define) &
              write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_network_set_global_indexes
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_network_set_dummy_elem ()
        !%------------------------------------------------------------------
        !% Description:
        !% Creates the indexes for the dummy elements in each of the
        !% elemXX arrays.
        !%-------------------------------------------------------------------
        !% Declarations:
            !integer       :: dummyIdx !% changed to global 20220726
            character(64) :: subroutine_name = 'init_network_set_dummy_elem'
        !%-------------------------------------------------------------------
        !% Preliminaries   
            !if (crashYN) return
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% --- indexes for the dummy elements
        dummyIdx = max_caf_elem_N + N_dummy_elem

        !% --- set the elem arrays with dummy values
        elemI(dummyIdx, ei_Lidx)        = dummyIdx
        elemI(dummyIdx,ei_elementType)  = dummy
        elemYN(dummyIdx,eYN_isDummy)    = .true.

        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%network_define) &
             write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_network_set_dummy_elem
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_network_handle_partition &
        (image, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, FaceGlobalCounter)
        !%-------------------------------------------------------------------
        !% Description
        !% Traverse through all the links and nodes in a partition and creates
        !% elements and faces. This subroutine assumes there will be at least
        !% one link in a partition.
        !%--------------------------------------------------------------------
        !% Declarations
            integer, intent(in)     :: image
            integer, intent(inout)  :: ElemLocalCounter, FaceLocalCounter
            integer, intent(inout)  :: ElemGlobalCounter, FaceGlobalCounter

            integer                 :: ii, pLink
            integer, pointer        :: thisLink, upNode, dnNode
            integer, dimension(:), allocatable, target :: packed_link_idx

            character(64) :: subroutine_name = 'init_network_handle_partition'
        !%--------------------------------------------------------------------
        !% Preliminaries
            !if (crashYN) return
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%--------------------------------------------------------------------
        !% --- pack all the link indexes in a partition
        packed_link_idx = pack(link%I(:,li_idx), (link%I(:,li_P_image) == image))

        !% --- find the number of links in a partition
        pLink = size(packed_link_idx)

        !% --- cycling through all the links in a partition
        do ii = 1,pLink
            !% --- necessary pointers to the link and its connected nodes
            thisLink => packed_link_idx(ii)
            upNode   => link%I(thisLink,li_Mnode_u)
            dnNode   => link%I(thisLink,li_Mnode_d)

            !% --- handle the upstream node of the link to create elements and faces
            call init_network_handle_upstreamnode &
                (image, thisLink, upNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                FaceGlobalCounter)

            !% --- handle the link to create elements and faces
            call init_network_handle_link &
                (image, thisLink, upNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                FaceGlobalCounter)
  
            !% --- handle the downstream node of the link to create elements and faces
            call init_network_handle_downstreamnode &
                (image, thisLink, dnNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                FaceGlobalCounter)
                
        end do
        
        !%--------------------------------------------------------------------
        !% Closing
            deallocate(packed_link_idx) !% deallocate temporary array
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_network_handle_partition
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_network_map_nodes (image)
        !%-----------------------------------------------------------------
        !% Description
        !% map all the interior junction faces in an image
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: image

            integer :: ii, pNodes
            integer, pointer :: thisJunctionNode, nodeType
            integer, dimension(:), allocatable, target :: packed_node_idx, JunctionElementIdx

            character(64) :: subroutine_name = 'init_network_map_nodes'
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        
        !% --- pack all the interior node indexes in a partition to find face maps
        packed_node_idx = pack(node%I(:,ni_idx),                       &
                              ((node%I(:,ni_P_image)   == image) .and. &
                              ((node%I(:,ni_node_type) == nJm ) .or.   &
                               (node%I(:,ni_node_type) == nJ2 ) ) )    )

        !% --- number of interior nodes in a partition
        pNodes = size(packed_node_idx)

        !% --- cycle through all the nJm nodes and set the face maps
        do ii = 1,pNodes
            thisJunctionNode => packed_node_idx(ii)
            nodeType         => node%I(thisJunctionNode,ni_node_type)

            select case (nodeType)
                case (nJm)
                    JunctionElementIdx = pack( elemI(:,ei_Lidx), &
                                                ( elemI(:,ei_node_Gidx_BIPquick) == thisJunctionNode) )

                    call init_network_map_nJm_branches (image, thisJunctionNode, JunctionElementIdx)

                    !% --- deallocate temporary array
                    deallocate(JunctionElementIdx)

                case (nJ2)
                    call init_network_map_nJ2 (image, thisJunctionNode)

                case default    
                    write(*,*) 'CODE ERROR: unexpected case default in ',trim(subroutine_name)
                    print *, 'Node Type # of ',nodeType
                    print *, 'which has key of ',trim(reverseKey(nodeType))
                    call util_crashpoint(39705)
            end select
        end do

        !%------------------------------------------------------------------
        !% Closing
            !% --- deallocate temporary array
            deallocate(packed_node_idx)

            if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_network_map_nodes
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_network_map_shared_faces (image)
        !%------------------------------------------------------------------
        !% Description:
        !% set the global indexes for shared faces across images
        !%-------------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: image

            integer :: ii, NsharedFaces
            integer, pointer :: fLidx, nIdx, eUp, eDn, targetImage, nodeType
            logical, pointer :: isUpGhost, isDnGhost
            integer, dimension(:), allocatable, target ::  sharedFaces

            character(64) :: subroutine_name = 'init_network_map_shared_faces'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------

        !% --- pack the shared faces in an image
        sharedFaces = pack(faceI(:,fi_Lidx), faceYN(:,fYN_isSharedFace))

        !% --- fid the size of the pack
        NsharedFaces = size(sharedFaces)

        do ii = 1,NsharedFaces
            fLidx       => sharedFaces(ii)
            nIdx        => faceI(fLidx,fi_node_idx_BIPquick)
            nodeType    => node%I(nIdx,ni_node_type)

            select case (nodeType)
                case (nJ2)
                    call init_network_map_shared_nJ2_nodes (image, fLidx, nIdx)

                case (nJm)
                    call init_network_map_shared_nJm_nodes (image, fLidx, nIdx)

                case default    
                    write(*,*) 'CODE ERROR: unexpected case default in ',trim(subroutine_name)
                    print *, 'Node Type # of ',nodeType
                    print *, 'which has key of ',trim(reverseKey(nodeType))
                    call util_crashpoint(55934)
            end select
        end do

        !%-------------------------------------------------------------------
        !% Closing
            !% --- deallocate temporary array
            deallocate(sharedFaces)

            if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_network_map_shared_faces
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_network_handle_upstreamnode &
        (image, thisLink, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
        FaceGlobalCounter)
        !%----------------------------------------------------------------
        !% Description:
        !% Handle the node upstream of a link to create elements and faces
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in)     :: image, thisLink, thisNode
            integer, intent(inout)  :: ElemLocalCounter, FaceLocalCounter
            integer, intent(inout)  :: ElemGlobalCounter, FaceGlobalCounter
            integer                 :: ii
            integer, pointer        :: nAssignStatus, nodeType, linkUp
            character(64) :: subroutine_name = 'init_network_handle_upstreamnode'
        !%-----------------------------------------------------------------
        !% Preliminaries
            !if (crashYN) return
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------
      
        !% --- Check 1: If the node is in the partition
        if (node%I(thisNode,ni_P_image) == image) then

            !% --- necessary pointers
            nAssignStatus => node%I(thisNode,ni_assigned)
            nodeType      => node%I(thisNode,ni_node_type)
            
            select case (nodeType)

                !% --- Handle upstream boundary nodes
                case(nBCup, nJ1)                                        
                    !% --- Check 2: If the node has already been assigned
                    if (nAssignStatus == nUnassigned) then

                        !% --- Advance the local and global counter for the upstream
                        !%     boundary node
                        FaceLocalCounter  = FaceLocalCounter  + oneI
                        FaceGlobalCounter = FaceGlobalCounter + oneI

                        !% --- integer data
                        faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                        faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                        !% --- an upstream boundary face does not have any local upstream element
                        !%     thus, it is mapped to the dummy element
                        faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem
                        !% --- the downstream element is stored
                        faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter
                        if (nodeType == nBCup) then
                            faceI(FaceLocalCounter,fi_BCtype)   = BCup
                        else 
                            faceI(FaceLocalCounter,fi_BCtype)   = BCnone
                        end if

                        !% --- map the elem_idx to the downstream element
                        node%I(thisNode,ni_elem_idx)        = faceI(FaceLocalCounter,fi_Melem_dL)
                        !% --- map the fnode face to the new face
                        node%I(thisNode,ni_face_idx)        = FaceLocalCounter
                        !% --- set zbottom
                        faceR(FaceLocalCounter,fr_Zbottom)  = node%R(thisNode,nr_Zbottom)
                        !% --- set the node the face has been originated from
                        faceI(FacelocalCounter,fi_node_idx_BIPquick) = thisNode
                        faceI(FacelocalCounter,fi_node_idx_SWMM)     = thisNode
                        !% --- change the node assignmebt value
                        nAssignStatus =  nAssigned
                    end if

                !% --- Handle 2 branch junction nodes
                case (nJ2)
                    !% --- Check 2: If the node has already been assigned
                    if (nAssignStatus == nUnassigned) then

                        !% --- Advance the local and global counter for the upstream
                        !%     nJ2 node
                        FaceLocalCounter  = FaceLocalCounter + oneI
                        FaceGlobalCounter = FaceGlobalCounter + oneI

                        !% --- integer data
                        faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                        faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                        faceI(FacelocalCounter,fi_Melem_uL) = ElemLocalCounter - oneI
                        faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter
                        faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
                        !% --- set zbottom
                        faceR(FaceLocalCounter,fr_Zbottom)  = node%R(thisNode,nr_Zbottom)
                        !% --- set the node the face has been originated from
                        faceI(FacelocalCounter,fi_node_idx_BIPquick) = thisNode
                        !% --- Assign the face index to the node
                        node%I(thisNode,ni_face_idx)         = FaceLocalCounter
                        !% --- Assign the upstream element index to the node
                        node%I(thisNode,ni_elem_idx)         = faceI(FacelocalCounter,fi_Melem_uL)


                        !% --- Check 3: If the node is an edge node (meaning this node is the
                        !%     connecting node across partitions)
                        if (node%I(thisNode,ni_P_is_boundary) == EdgeNode) then

                            !% --- An upstream edge node indicates there are no local
                            !%     elements upstream of that node. Thus it is mapped to
                            !%     the dummy element
                            faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem

                            !% --- logical data
                            faceYN(FaceLocalCounter,fYN_isSharedFace) = .true.
                            faceYN(FaceLocalCounter,fYN_isUpGhost)    = .true.

                            !% --- find the connecting image to this face
                            linkUp  => node%I(thisNode,ni_Mlink_u1)

                            faceI(FaceLocalCounter,fi_Connected_image) = link%I(linkUp,li_P_image)

                            if (image < faceI(FaceLocalCounter,fi_Connected_image)) then
                                !% --- we only set the global indexes where the connection
                                !%     is in higher order than the current image.
                                !%     (for example if current image = 1 and connection is 2,
                                !%     we set the global counter. But when the current image = 2 but
                                !%     the connection is 1, we set it from init_network_map_shared_faces
                                !%     subroutine)
                                faceI(FaceLocalCounter,fi_Gidx) = FaceGlobalCounter
                            else
                                !% --- set global index as nullvalue for shared faces.
                                !%     these global indexes will be set later
                                faceI(FaceLocalCounter,fi_Gidx)     = nullvalueI
                            end if

                            !% --- set the swmm idx. if the node is phantom, it will not have
                            !%     any SWMM idx
                            if (.not. node%YN(thisNode,nYN_is_phantom_node)) then
                                faceI(FacelocalCounter,fi_node_idx_SWMM) = thisNode
                            endif

                        else
                            !% --- the node is not a edge node thus, the node cannot be a phantom node
                            !%     and the bipquick and SWMM idx will be the same
                            faceI(FacelocalCounter,fi_node_idx_SWMM) = thisNode
                        end if

                        !% --- change the node assignmebt value
                        nAssignStatus =  nAssigned
                    end if

                !% --- Handle junction nodes with more than 2 branches (multi branch junction node).
                case (nJm)

                    !% --- Check 2: If the node has already been assigned
                    if (nAssignStatus == nUnassigned) then

                        !% --- multibranch junction nodes will have both elements and faces.
                        !%     thus, a seperate subroutine is required to handle these nodes
                        call init_network_handle_nJm &
                            (image, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                            FaceGlobalCounter, nAssignStatus)
                    end if

                case default    
                    write(*,*) 'CODE ERROR: unexpected case default in ',trim(subroutine_name)
                    print*, 'error: node ' // node%Names(thisNode)%str // &
                            ' has an unexpected nodeType', nodeType
                    print *, 'which has key of ',trim(reverseKey(nodeType))
                    call util_crashpoint(987034)

            end select

        !% --- handle the node if it is not in the partition
        else

            !% --- Advance the local and global counter for the upstream
            !%     node outside the partition
            FaceLocalCounter  = FaceLocalCounter + oneI
            FaceGlobalCounter = FaceGlobalCounter + oneI

            !% --- integer data
            faceI(FacelocalCounter,fi_Lidx)     = FacelocalCounter
            faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
            faceI(FacelocalCounter,fi_Connected_image)   = node%I(thisNode,ni_P_image)
            faceI(FacelocalCounter,fi_link_idx_BIPquick) = thisLink
            faceI(FacelocalCounter,fi_link_idx_SWMM)     = link%I(thisLink,li_parent_link)

            !% --- set the face from the node it has been originated from
            faceI(FacelocalCounter,fi_node_idx_BIPquick) = thisNode

            !% --- Set the swmm idx.
            !%     If the node is phantom, it will not have any SWMM idx
            if (.not. node%YN(thisNode,nYN_is_phantom_node)) then
                faceI(FacelocalCounter,fi_node_idx_SWMM) = thisNode
            endif

            !% --- real data
            faceR(FaceLocalCounter,fr_Zbottom) = node%R(thisNode,nr_Zbottom)

            !% --- if the upstream node is not in the partiton,
            !%     the face map to upstream is mapped to
            !%     the dummy element
            faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem

            !% --- since no upstream node indicates start of a partiton,
            !%     the downstream element will be initialized elem idx
            faceI(FacelocalCounter,fi_Melem_dL) = ElemLocalCounter

            !% --- special condition if a element will be immediate downsteam of a JB
            if (node%I(thisNode,ni_node_type) == nJm) then
                elemYN(ElemLocalCounter,eYN_isElementDownstreamOfJB) = .true.
            endif

            !% --- since this is a shared face, it will have a copy in other image and they will
            !%     both share same global index. so, the face immediately after this shared face
            !%     will have the global index set from the init_network_set_global_indexes subroutine.
            !%     However, since the init_network_handle_link subroutine will advance the global face
            !%     count anyway, the count  here is needed to be adjusted by substracting one from the
            !%     count.
            FaceGlobalCounter = FaceGlobalCounter - oneI

            !% --- logical data
            faceYN(FacelocalCounter,fYN_isSharedFace)   = .true.
            faceYN(FacelocalCounter,fYN_isUpGhost)      = .true.

            if (image < faceI(FaceLocalCounter,fi_Connected_image)) then
                !% --- we only set the global indexes where the connection
                !%     is in higher order than the current image.
                !%     (for example if current image = 1 and connection is 2,
                !%     we set the global counter. But when the current image = 2 but
                !%     the connection is 1, we set it from init_network_map_shared_faces
                !%     subroutine)
                faceI(FaceLocalCounter,fi_Gidx) = FaceGlobalCounter
            else
                !% --- set global index as nullvalue for shared faces.
                !%     these global indexes will be set later
                faceI(FaceLocalCounter,fi_Gidx)     = nullvalueI
            end if
        end if

        !%-----------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_network_handle_upstreamnode
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_network_handle_link &
        (image, thisLink, upNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
        FaceGlobalCounter)
        !%-----------------------------------------------------------------
        !% Description:
        !% handle links in a partition
        !%-----------------------------------------------------------------
        !% Declarations:
        integer, intent(in)     :: image, thisLink, upNode
            integer, intent(inout)  :: ElemLocalCounter, FaceLocalCounter
            integer, intent(inout)  :: ElemGlobalCounter, FaceGlobalCounter

            integer                 :: ii
            real(8)                 :: zCenter, zDownstream
            integer, pointer        :: lAssignStatus, NlinkElem
            real(8), pointer        :: zUpstream

            character(64) :: subroutine_name = 'init_network_handle_link'

        !%-----------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------
        !% Aliases
            lAssignStatus => link%I(thisLink,li_assigned)
        !%-----------------------------------------------------------------

        if (lAssignStatus == lUnAssigned) then
            NlinkElem     => link%I(thisLink,li_N_element)
            zUpstream     => link%R(thisLink,lr_ZbottomUp)

            !% --- store the ID of the first (upstream) element in this link
            link%I(thisLink,li_first_elem_idx)   = ElemLocalCounter
            !% --- reference elevations at cell center
            zCenter     = zUpstream - onehalfR * link%R(thisLink,lr_ElementLength) * link%R(thisLink,lr_Slope)
            zDownstream = zUpstream            - link%R(thisLink,lr_ElementLength) * link%R(thisLink,lr_Slope)
            do ii = 1, NlinkElem
                !%................................................................
                !% Element arrays update
                !%................................................................

                !% --- integer data
                elemI(ElemLocalCounter,ei_Lidx)                 = ElemLocalCounter
                elemI(ElemLocalCounter,ei_Gidx)                 = ElemGlobalCounter

                !% --- set the element type
                if ((link%I(thisLink,li_link_type) == lPipe) .or. &
                    (link%I(thisLink,li_link_type) == lChannel)) then
                    elemI(ElemLocalCounter,ei_elementType)      = CC
                elseif (link%I(thisLink,li_link_type) == lWeir) then
                    elemI(ElemLocalCounter,ei_elementType)      = weir
                elseif (link%I(thisLink,li_link_type) == lOrifice) then
                    elemI(ElemLocalCounter,ei_elementType)      = orifice
                elseif (link%I(thisLink,li_link_type) == lPump) then
                    elemI(ElemLocalCounter,ei_elementType)      = pump
                elseif (link%I(thisLink,li_link_type) == lOutlet) then
                    elemI(ElemLocalCounter,ei_elementType)      = outlet
                endif

                elemI(ElemLocalCounter,ei_Mface_uL)             = FaceLocalCounter
                elemI(ElemLocalCounter,ei_Mface_dL)             = FaceLocalCounter + oneI
                elemI(ElemLocalCounter,ei_link_pos)             = ii
                elemI(ElemLocalCounter,ei_link_Gidx_BIPquick)   = thisLink
                elemI(ElemLocalCounter,ei_link_Gidx_SWMM)       = link%I(thisLink,li_parent_link)

                !% --- real data
                elemR(ElemLocalCounter,er_Length)           = link%R(thisLink,lr_ElementLength)
                elemR(ElemLocalCounter,er_Zbottom)          = zCenter

                !%................................................................
                !% Face arrays update
                !%................................................................

                if (ii < NlinkElem) then
                !% --- advance only the downstream interior face counter of a link element
                    FaceLocalCounter  = FaceLocalCounter  + oneI
                    FaceGlobalCounter = FaceGlobalCounter + oneI

                    !% --- face integer data
                    faceI(FaceLocalCounter,fi_Lidx)              = FaceLocalCounter
                    faceI(FaceLocalCounter,fi_Gidx)              = FaceGlobalCounter
                    faceI(FaceLocalCounter,fi_Melem_dL)          = ElemLocalCounter + oneI
                    faceI(FaceLocalCounter,fi_Melem_uL)          = ElemLocalCounter
                    faceI(FaceLocalCounter,fi_BCtype)            = doesnotexist
                    faceR(FaceLocalCounter,fr_Zbottom)           = zDownstream
                    faceI(FaceLocalCounter,fi_link_idx_BIPquick) = thisLink
                    faceI(FaceLocalCounter,fi_link_idx_SWMM)     = link%I(thisLink,li_parent_link)
                end if

                !% --- counter for element z bottom calculation
                zCenter     = zCenter     - link%R(thisLink,lr_ElementLength) * link%R(thisLink,lr_Slope)
                zDownstream = zDownstream - link%R(thisLink,lr_ElementLength) * link%R(thisLink,lr_Slope)

                !% --- Advance the element counter
                ElemLocalCounter  = ElemLocalCounter  + oneI
                ElemGlobalCounter = ElemGlobalCounter + oneI
            end do

            lAssignStatus = lAssigned
            link%I(thisLink,li_last_elem_idx)    = ElemLocalCounter - oneI

        end if

        !%-----------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_network_handle_link
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_network_handle_downstreamnode &
        (image, thisLink, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
        FaceGlobalCounter)
        !%-----------------------------------------------------------------
        !% Description
        !% handle the node downstream of a link
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in)    :: image, thisLink, thisNode
            integer, intent(inout) :: ElemLocalCounter, FaceLocalCounter
            integer, intent(inout) :: ElemGlobalCounter, FaceGlobalCounter

            integer, pointer :: nAssignStatus, nodeType, linkDn

            character(64) :: subroutine_name = 'init_network_handle_downstreamnode'
        !%-----------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------
        
        !% --- Check 1: Is the node is in the partition
        if (node%I(thisNode,ni_P_image) == image) then

            !% --- necessary pointers
            nAssignStatus => node%I(thisNode,ni_assigned)
            nodeType      => node%I(thisNode,ni_node_type)

            select case (nodeType)

                case(nBCdn)
                    !% --- Check 2: If the node has already been assigned
                    if (nAssignStatus == nUnassigned) then

                        !% --- Advance face local and global counters for BCdn node
                        FaceLocalCounter  = FaceLocalCounter  + oneI
                        FaceGlobalCounter = FaceGlobalCounter + oneI

                        !% --- integer data
                        !%     an downstream boundary face does not have any local downstream element
                        !%     thus, it is mapped to the dummy element
                        faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                        faceI(FacelocalCounter,fi_Gidx)     = FaceGlobalCounter
                        faceI(FacelocalCounter,fi_Melem_uL) = ElemLocalCounter - oneI
                        faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem
                        faceI(FaceLocalCounter,fi_BCtype)   = BCdn
                        !% --- set zbottom
                        faceR(FaceLocalCounter,fr_Zbottom)  = node%R(thisNode,nr_Zbottom)

                        !% --- node face is this face value
                        node%I(thisNode,ni_face_idx)        = FaceLocalCounter
                        !% --- nBCdn node element is mapped to usptream of face
                        node%I(thisNode,ni_elem_idx)        = faceI(FacelocalCounter,fi_Melem_uL)

                        !% --- set the node the face has been originated from
                        faceI(FacelocalCounter,fi_node_idx_BIPquick) = thisNode
                        faceI(FacelocalCounter,fi_node_idx_SWMM)     = thisNode

                        !% --- change the node assignmebt value
                        nAssignStatus =  nAssigned
                    end if

                case (nJ2)
                    !% --- Check 2: If the node has already been assigned
                    if (nAssignStatus == nUnassigned) then

                        !% --- Advance face local and global counters for nJ2 node
                        FaceLocalCounter  = FaceLocalCounter  + oneI
                        FaceGlobalCounter = FaceGlobalCounter + oneI

                        !% --- integer data
                        faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                        faceI(FaceLocalCounter,fi_Gidx)     = FaceGlobalCounter
                        faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
                        faceI(FacelocalCounter,fi_Melem_uL) = ElemLocalCounter - oneI
                        faceI(FacelocalCounter,fi_Melem_dL) = ElemLocalCounter
                        !% --- set zbottom
                        faceR(FaceLocalCounter,fr_Zbottom)  = node%R(thisNode,nr_Zbottom)
                        !% --- node face value is this face
                        node%I(thisNode,ni_face_idx)        = FaceLocalCounter
                        !% --- node element value is the upstream element
                        node%I(thisNode,ni_elem_idx)        = faceI(FacelocalCounter,fi_Melem_uL) 

                        !% --- set the node the face has been originated from
                        faceI(FacelocalCounter,fi_node_idx_BIPquick) = thisNode

                        !% --- integer data
                        if (node%I(thisNode,ni_P_is_boundary) == EdgeNode) then

                            !% --- A downstream edge node indicates there are no local
                            !%     elements downstream of that node
                            faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem

                            !% --- logical data
                            faceYN(FaceLocalCounter,fYN_isSharedFace) = .true.
                            faceYN(FaceLocalCounter,fYN_isDnGhost)    = .true.

                            !% --- find the connecting image to this face
                            linkDn  => node%I(thisNode,ni_Mlink_d1)

                            faceI(FaceLocalCounter,fi_Connected_image)    = link%I(linkDn,li_P_image)

                            if (image < faceI(FaceLocalCounter,fi_Connected_image)) then
                                !% --- we only set the global indexes where the connection
                                !%     is in higher order than the current image.
                                !%     (for example if current image = 1 and connection is 2,
                                !%     we set the global counter. But when the current image = 2 but
                                !%     the connection is 1, we set it from init_network_map_shared_faces
                                !%     subroutine)
                                faceI(FaceLocalCounter,fi_Gidx) = FaceGlobalCounter
                            else
                                !% --- set global index as nullvalue for shared faces.
                                !%     these global indexes will be set later
                                faceI(FaceLocalCounter,fi_Gidx)     = nullvalueI
                            end if

                            !% --- set the swmm idx.
                            !%     if the node is phantom, it will not have any SWMM idx
                            if (.not. node%YN(thisNode,nYN_is_phantom_node)) then
                                faceI(FacelocalCounter,fi_node_idx_SWMM) = thisNode
                            endif
                        else
                            !% --- the node is not a edge node thus, the node cannot be a phantom node
                            !%     and the bipquick and SWMM idx will be the same
                            faceI(FacelocalCounter,fi_node_idx_SWMM) = thisNode
                        end if

                        !% --- change the node assignmebt value
                        nAssignStatus =  nAssigned
                    end if

                case (nJm)
                    !% --- Check 2: If the node has already been assigned
                    if (nAssignStatus == nUnassigned) then

                        call init_network_handle_nJm &
                            (image, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
                            FaceGlobalCounter, nAssignStatus)

                    end if

                case default

                    print *
                    print *, 'In ', subroutine_name
                    print *, 'CODE ERROR: node ' // node%Names(thisNode)%str // &
                            ' has an unexpected nodeType', nodeType
                    call util_crashpoint(398704)

            end select

        else
            !% --- Advance face local and global counters for nodes outside of the partition
            FaceLocalCounter  = FaceLocalCounter  + oneI
            FaceGlobalCounter = FaceGlobalCounter + oneI
            !% --- if the downstream node is not in the partiton.
            !%     through subdivide_link_going_downstream subroutine
            !%     upstream map to the element has alrady been set.
            !%     However, downstream map has set to wrong value.
            !%     Thus, setting the map elem ds to dummy elem
            !% --- integer data
            faceI(FacelocalCounter,fi_Lidx)     = FaceLocalCounter
            faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem
            faceI(FacelocalCounter,fi_Melem_uL) = ElemLocalCounter - oneI
            faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist
            faceI(FacelocalCounter,fi_Connected_image)   = node%I(thisNode,ni_P_image)
            !% --- set the node the face has been originated from
            faceI(FacelocalCounter,fi_node_idx_BIPquick) = thisNode
            faceI(FacelocalCounter,fi_link_idx_BIPquick) = thisLink
            faceI(FaceLocalCounter,fi_link_idx_SWMM)     = link%I(thisLink,li_parent_link)

            !% --- Set the swmm idx.
            !%     If the node is phantom, it will not have any SWMM idx
            if (.not. node%YN(thisNode,nYN_is_phantom_node)) then
                faceI(FacelocalCounter,fi_node_idx_SWMM) = thisNode
            endif

            !% --- real data
            faceR(FaceLocalCounter,fr_Zbottom) = node%R(thisNode,nr_Zbottom)

            !% --- special condition if a element will be immediate upsteam of a JB
            if (node%I(thisNode,ni_node_type) == nJm) then
                elemYN(ElemLocalCounter,eYN_isElementUpstreamOfJB) = .true.
            endif

            !% --- logical data
            faceYN(FacelocalCounter,fYN_isSharedFace) = .true.
            faceYN(FaceLocalCounter,fYN_isDnGhost)    = .true.

            if (image < faceI(FaceLocalCounter,fi_Connected_image)) then
                !% --- we only set the global indexes where the connection
                !%     is in higher order than the current image.
                !%     (for example if current image = 1 and connection is 2,
                !%     we set the global counter. But when the current image = 2 but
                !%     the connection is 1, we set it from init_network_map_shared_faces
                !%     subroutine)
                faceI(FaceLocalCounter,fi_Gidx) = FaceGlobalCounter
            else
                !% --- set global index as nullvalue for shared faces.
                !%     these global indexes will be set later
                faceI(FaceLocalCounter,fi_Gidx)     = nullvalueI
            end if

        end if

        !%-----------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_network_handle_downstreamnode
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_network_handle_nJm &
        (image, thisNode, ElemLocalCounter, FaceLocalCounter, ElemGlobalCounter, &
        FaceGlobalCounter, nAssignStatus)
        !%-----------------------------------------------------------------
        !% Description
        !% subdivides the multi branch junctions into elements and faces
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in)    :: image, thisNode
            integer, intent(inout) :: ElemLocalCounter, FaceLocalCounter
            integer, intent(inout) :: ElemGlobalCounter, FaceGlobalCounter
            integer, intent(inout) :: nAssignStatus

            integer, pointer :: upBranchIdx, dnBranchIdx

            integer :: ii, upBranchSelector, dnBranchSelector

            character(64) :: subroutine_name = 'init_network_handle_nJm'
        !%-----------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------

        !%................................................................
        !% Junction main
        !%................................................................

        !% --- Element Arrays
        !% --- integer data
        elemI(ElemLocalCounter,ei_Lidx)                 = ElemLocalCounter
        elemI(ElemLocalCounter,ei_Gidx)                 = ElemGlobalCounter
        elemI(ElemLocalCounter,ei_elementType)          = JM
        elemI(ElemLocalCounter,ei_node_Gidx_BIPquick)   = thisNode
        !% --- a JM node will never be a phantom node. Thus, the BQuick and SWMM idx will be the same
        elemI(ElemLocalCounter,ei_node_Gidx_SWMM)       = thisNode
        !% --- save the number total number of real branches connected to a junction main
        elemSI(ElemLocalCounter,esi_JunctionMain_Total_Branches) = node%I(thisNode,ni_N_link_u) + node%I(thisNode,ni_N_link_d) 
        !% --- Assign junction main element to node
        node%I(thisNode,ni_elem_idx)                    = ElemLocalCounter
        !% --- a JM node is connected to faces through multipol JB branches, so the node face index is null
        node%I(thisNode,ni_face_idx)                    = nullvalueI

        !% --- real data
        elemR(ElemLocalCounter,er_Zbottom) = node%R(thisNode,nr_zbottom)

        !% --- Advance the element counter to 1st upstream branch
        ElemLocalCounter  = ElemLocalCounter  + oneI
        ElemGlobalCounter = ElemGlobalCounter + oneI

        !%................................................................
        !% Handle Junction Branches
        !%................................................................

        !% --- initialize selectors for upstream and downstream branches
        upBranchSelector = zeroI
        dnBranchSelector = zeroI

        !% --- loop through all the branches
        do ii = 1,max_branch_per_node

            !% --- common junction branch data
            !%     element arrays
            !%     integer data
            elemI(ElemLocalCounter,ei_Lidx)           = ElemLocalCounter
            elemI(ElemLocalCounter,ei_Gidx)           = ElemGlobalCounter
            elemI(ElemLocalCounter,ei_elementType)    = JB
            elemI(ElemLocalCounter,ei_node_Gidx_BIPquick)   = thisNode
            !% --- A JB will never come from a phantom node.
            !%     Thus, the BQuick and SWMM idx will be the same
            elemI(ElemLocalCounter,ei_node_Gidx_SWMM)       = thisNode

            !% --- real data
            elemR(ElemLocalCounter,er_Zbottom) = node%R(thisNode,nr_zbottom)

            !% --- advance the face counters for the branch
            FaceLocalCounter  = FaceLocalCounter  + oneI
            FaceGlobalCounter = FaceGlobalCounter + oneI

            !%......................................................
            !% Upstream Branches
            !%......................................................
            select case (mod(ii,2))
                case (1)
                !% --- finds odd number branches
                !%     all the odd numbers are upstream branches
                    upBranchSelector = upBranchSelector + oneI
                    !% --- pointer to upstream branch
                    upBranchIdx => node%I(thisNode,ni_idx_base1 + upBranchSelector)

                    !% --- elem array
                    !%     map the upstream face of the branch element
                    elemI(ElemLocalCounter,ei_Mface_uL) = FaceLocalCounter

                    !% --- face array
                    !%     integer data
                    faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                    faceI(FacelocalCounter,fi_Gidx)     = FaceGlobalCounter
                    faceI(FaceLocalCounter,fi_Melem_dL) = ElemLocalCounter
                    faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist

                    !% --- set the node the face has been originated from
                    faceI(FacelocalCounter,fi_node_idx_BIPquick) = thisNode
                    faceI(FaceLocalCounter,fi_node_idx_SWMM)     = thisNode

                    !% --- real branches
                    if (upBranchIdx /= nullvalueI) then
                        !% --- integer data
                        elemSI(ElemLocalCounter,esi_JunctionBranch_Exists)           = oneI
                        elemSI(ElemLocalCounter,esi_JunctionBranch_Link_Connection)  = upBranchIdx
                        elemR(ElemLocalCounter,er_Length) = setting%Discretization%NominalElemLength
                        elemYN(ElemLocalCounter,eYN_isElementUpstreamOfJB) = .true.
                        faceI(FaceLocalCounter,fi_link_idx_BIPquick) = upBranchIdx
                        faceI(FaceLocalCounter,fi_link_idx_SWMM)     = link%I(upBranchIdx,li_parent_link)
                        !% --- set zbottom
                        elemR(ElemLocalCounter,er_Zbottom)  = link%R(upBranchIdx,lr_ZbottomDn)
                        faceR(FaceLocalCounter,fr_Zbottom)  = elemR(ElemLocalCounter,er_Zbottom)
                        !% --- identifier for downstream junction branch faces
                        faceYN(FaceLocalCounter,fYN_isUpstreamJBFace) = .true.

                        !% --- Check 4: this node is the connecting node across partitions
                        if ( (node%I(thisNode,ni_P_is_boundary) == EdgeNode) &
                             .and.                                           &
                             (link%I(upBranchIdx,li_P_image)    /= image   ) &
                           )  then

                            !% --- since this is a shared face in the upstream direction,
                            !%     the up_map is set to dummy element
                            faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem
                            faceI(FaceLocalCounter,fi_Connected_image) = link%I(upBranchIdx,li_P_image)

                            !% --- logical data
                            faceYN(FacelocalCounter,fYN_isSharedFace) = .true.
                            faceYN(FaceLocalCounter,fYN_isUpGhost)    = .true.

                            if (image < faceI(FaceLocalCounter,fi_Connected_image)) then
                                !% --- we only set the global indexes where the connection
                                !%     is in higher order than the current image.
                                !%     (for example if current image = 1 and connection is 2,
                                !%     we set the global counter. But when the current image = 2 but
                                !%     the connection is 1, we set it from init_network_map_shared_faces
                                !%     subroutine)
                                faceI(FaceLocalCounter,fi_Gidx) = FaceGlobalCounter
                            else
                                !% --- set global index as nullvalue for shared faces.
                                !%     these global indexes will be set later
                                faceI(FaceLocalCounter,fi_Gidx)     = nullvalueI
                            end if
                        end if
                    else

                        !% --- since this is a null face in the upstream direction,
                        !%     the up_map is set to dummy element
                        faceI(FaceLocalCounter,fi_Melem_uL) = max_caf_elem_N + N_dummy_elem

                        call init_network_nullify_nJm_branch &
                            (ElemLocalCounter, FaceLocalCounter)
                    end if

                !%......................................................
                !% Downstream Branches
                !%......................................................
                case (0)
                    !% --- even number branches
                    !%     all the even numbers are downstream branches
                    dnBranchSelector = dnBranchSelector + oneI
                    !% --- pointer to upstream branch
                    dnBranchIdx => node%I(thisNode,ni_idx_base2 + dnBranchSelector)

                    !% --- elem array
                    !%     integer data
                    elemI(ElemLocalCounter,ei_Mface_dL) = FaceLocalCounter

                    !% --- face array
                    !%     integer data
                    faceI(FaceLocalCounter,fi_Lidx)     = FaceLocalCounter
                    faceI(FacelocalCounter,fi_Gidx)     = FaceGlobalCounter
                    faceI(FaceLocalCounter,fi_Melem_uL) = ElemLocalCounter
                    faceI(FaceLocalCounter,fi_BCtype)   = doesnotexist

                    !% --- set the node the face has been originated from
                    faceI(FacelocalCounter,fi_node_idx_BIPquick) = thisNode
                    faceI(FaceLocalCounter,fi_node_idx_SWMM)     = thisNode

                    !% --- Check 3: if the branch is a valid branch
                    if (dnBranchIdx /= nullvalueI) then
                        !% --- integer data
                        elemSI(ElemLocalCounter,esi_JunctionBranch_Exists)          = oneI
                        elemSI(ElemLocalCounter,esi_JunctionBranch_Link_Connection) = dnBranchIdx
                        elemR(ElemLocalCounter,er_Length) = setting%Discretization%NominalElemLength
                        elemYN(ElemLocalCounter,eYN_isElementDownstreamOfJB) = .true.
                        faceI(FacelocalCounter,fi_link_idx_BIPquick) = dnBranchIdx
                        faceI(FaceLocalCounter,fi_link_idx_SWMM)     = link%I(dnBranchIdx,li_parent_link)
                        !% --- set zbottom
                        elemR(ElemLocalCounter,er_Zbottom)  = link%R(dnBranchIdx,lr_ZbottomUp)
                        faceR(FaceLocalCounter,fr_Zbottom)  = elemR(ElemLocalCounter,er_Zbottom)
                        !% --- identifier for downstream junction branch faces
                        faceYN(FaceLocalCounter,fYN_isDownstreamJBFace) = .true.

                        !% --- Check 4: if the link connecting this branch is a part of this partition and
                        !%     the node is not an edge node (meaning this node is the connecting node
                        !%     across partitions)
                        if ( (node%I(thisNode,ni_P_is_boundary) == EdgeNode)  .and. &
                            (link%I(dnBranchIdx,li_P_image)    /= image   ) )  then

                            !% --- since this is a shared face in the downstream direction,
                            !%     the dn_map is set to dummy element
                            faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem
                            faceI(FaceLocalCounter,fi_Connected_image) = link%I(dnBranchIdx,li_P_image)

                            !% --- logical data
                            faceYN(FacelocalCounter,fYN_isSharedFace) = .true.
                            faceYN(FaceLocalCounter,fYN_isDnGhost)    = .true.

                            if (image < faceI(FaceLocalCounter,fi_Connected_image)) then
                                !% --- we only set the global indexes where the connection
                                !%     is in higher order than the current image.
                                !%     (for example if current image = 1 and connection is 2,
                                !%     we set the global counter. But when the current image = 2 but
                                !%     the connection is 1, we set it from init_network_map_shared_faces
                                !%     subroutine)
                                faceI(FaceLocalCounter,fi_Gidx) = FaceGlobalCounter
                            else
                                !% --- set global index as nullvalue for shared faces.
                                !%     these global indexes will be set later
                                faceI(FaceLocalCounter,fi_Gidx) = nullvalueI
                            end if
                        end if

                    else
                        !% --- since this is a null face in the downstream direction,
                        !%     the dn_map is set to dummy element
                        faceI(FaceLocalCounter,fi_Melem_dL) = max_caf_elem_N + N_dummy_elem

                        call init_network_nullify_nJm_branch &
                            (ElemLocalCounter, FaceLocalCounter)
                    end if

                case default
                    print *, 'CODE ERROR: unexpected case default'
                    print *, 'unsupported value for mod(ii,2) of ',mod(ii,2)
                    call util_crashpoint(65874)
            end select

            !% --- Advance the element counter for next branch
            ElemLocalCounter  = ElemLocalCounter  + oneI
            ElemGlobalCounter = ElemGlobalCounter + oneI
        end do

        !% --- set status to assigned
        nAssignStatus = nAssigned

        !%-----------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_network_handle_nJm
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_network_map_nJm_branches (image, thisJNode, JelemIdx)
        !%-----------------------------------------------------------------
        !% Description:
        !% maps all the multi-branch junction elements
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in)                       :: image, thisJNode
            integer, dimension(:), target, intent(in) :: JelemIdx

            integer          :: ii, upBranchSelector, dnBranchSelector
            integer          :: LinkFirstElem, LinkLastElem
            integer, pointer :: upBranchIdx, dnBranchIdx
            integer, pointer :: eIdx, fLidx
            character(64) :: subroutine_name = 'init_network_map_nJm_branches'
        !%-----------------------------------------------------------------
        !% Preliminaries
        if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------

        !% --- initialize selectors for upstream and downstream branches
        upBranchSelector = zeroI
        dnBranchSelector = zeroI

        !% --- cycle through the junction elements of map faces
        do ii = 1, max_branch_per_node + 1

            !% --- Cycling through all the junction elements including
            !%     junction main. Since the JM is ii=1 the 
            !%     upstream elements are even and the downstream
            !%     are odd

            if (ii > 1) then !% ii=1 is junction element
                !% --- all the even numbers are upstream branch elements
                select case (mod(ii,2))
                    case(0)
                        upBranchSelector = upBranchSelector + oneI
                        !% --- pointer to upstream branch
                        upBranchIdx => node%I(thisJNode,ni_idx_base1 + upBranchSelector)

                        !% --- condition for a link connecting this branch is valid and
                        !%     included in this partition.

                        if (upBranchIdx /= nullvalueI) then
                            if (link%I(upBranchIdx,li_P_image) == image) then
                                !% --- find the last element index of the link
                                LinkLastElem = link%I(upBranchIdx,li_last_elem_idx)

                                !% --- pointer to the specific branch element
                                eIdx => JelemIdx(ii)

                                !% --- find the downstream face index of that last element
                                fLidx => elemI(eIdx,ei_Mface_uL)

                                !% --- if the face is a shared face across images,
                                !%     it will not have any upstream local element
                                if ( .not. faceYN(fLidx,fYN_isSharedFace)) then
                                    !% --- the upstream face of the upstream branch will be the
                                    !%     last downstream face of the connected link
                                    !%     here, one important thing to remember is that
                                    !%     the upstrem branch elements does not have any
                                    !%     downstream faces.

                                    !% --- local d/s face map to element u/s of the branch
                                    elemI(LinkLastElem,ei_Mface_dL) = fLidx
                                    !% --- set the first element as the immediate downstream element of a JB
                                    elemYN(LinkLastElem,eYN_isElementUpstreamOfJB) = .true.

                                    !% --- local u/s element of the face
                                    faceI(fLidx,fi_Melem_uL) = LinkLastElem
                                end if
                            end if
                        end if

                    !% --- all odd numbers starting from 3 are downstream branch elements
                    case (1)
                        dnBranchSelector = dnBranchSelector + oneI
                        !% --- pointer to upstream branch
                        dnBranchIdx => node%I(thisJNode,ni_idx_base2 + dnBranchSelector)

                        !% --- condition for a link connecting this branch is valid and
                        !%     included in this partition.
                        if (dnBranchIdx /= nullvalueI)  then
                            if (link%I(dnBranchIdx,li_P_image) == image) then

                                !% --- find the first element index of the link
                                LinkFirstElem = link%I(dnBranchIdx,li_first_elem_idx)

                                !% --- pointer to the specific branch element
                                eIdx => JelemIdx(ii)

                                !% --- find the downstream face index of that last element
                                fLidx => elemI(eIdx,ei_Mface_dL)

                                !% --- if the face is a shared face across images,
                                !%     it will not have any upstream local element
                                !%     (not sure if we need this condition)
                                if ( .not. faceYN(fLidx,fYN_isSharedFace)) then

                                    !% --- the downstream face of the downstream branch will be the
                                    !%     first upstream face of the connected link
                                    !%     here, one important thing to remember is that
                                    !%     the downstream branch elements does not have any
                                    !%     upstream faces.

                                    !% --- local map to upstream face for elemI
                                    elemI(LinkFirstElem,ei_Mface_uL) = fLidx
                                    !% --- set the first element as the immediate downstream element of a JB
                                    elemYN(LinkFirstElem,eYN_isElementDownstreamOfJB) = .true.

                                    !% local downstream element of the face
                                    faceI(fLidx,fi_Melem_dL) = LinkFirstElem
                                end if
                            end if
                        end if

                    case default
                        print *, 'CODE ERROR: unexpected case default'
                        print *, 'Unsupported value for mod(ii,2) of ',mod(ii,2)
                        call util_crashpoint(99374)    
                end select
            end if
        end do

        !%-----------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_network_map_nJm_branches
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_network_map_shared_nJm_nodes (image, fLidx, nIdx)
        !%-----------------------------------------------------------------
        !% Description
        !% set the global index, map, and ghost element for nJm nodes
        !%-----------------------------------------------------------------
        !% Declarations:
        integer, intent(in) :: image, fLidx, nIdx

        integer             :: ii
        integer, pointer    :: fGidx, eUp, eDn, targetImage, branchIdx

        character(64) :: subroutine_name = 'init_network_map_shared_nJm_nodes'
        !%-----------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------
        !% Aliases
            fGidx       => faceI(fLidx,fi_Gidx)
            eUp         => faceI(fLidx,fi_Melem_uL)
            eDn         => faceI(fLidx,fi_Melem_dL)
            targetImage => faceI(fLidx,fi_Connected_image)
            branchIdx   => faceI(fLidx,fi_link_idx_BIPquick)
        !%-----------------------------------------------------------------

        do ii = 1,N_face(targetImage)

            if ((faceI(ii,fi_Connected_image)[targetImage]   == image)   .and. &
                (faceI(ii,fi_node_idx_BIPquick)[targetImage] == nIdx )   .and. &
                (faceI(ii,fi_link_idx_BIPquick)[targetImage] == branchIdx)) then

                !% --- set the local index of the indetical face
                faceI(ii,fi_Identical_Lidx)[targetImage] = fLidx

                !% --- find the local ghost element index of the connected image
                if (faceYN(ii,fYN_isUpGhost)[targetImage]) then
                    faceI(ii,fi_GhostElem_uL)[targetImage] = eUp

                elseif (faceYN(ii,fYN_isDnGhost)[targetImage]) then
                    faceI(ii,fi_GhostElem_dL)[targetImage] = eDn
                end if

                !% --- find the global index and set to target image
                if (faceI(ii,fi_Gidx)[targetImage] == nullvalueI) then
                    faceI(ii,fi_Gidx)[targetImage] = fGidx
                end if
            end if
        end do

        !%-----------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%network_define) &

            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_network_map_shared_nJm_nodes
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_network_map_nJ2 (image, thisJNode)
        !%-----------------------------------------------------------------
        !% Description:
        !% map all the nJ2 nodes. All the nJ2 node maps are handeled in the partition
        !% this is for the special cases where a disconnected nJ2 has not been mapped
        !% properly
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: image, thisJNode

            integer             :: ii, upBranchSelector, dnBranchSelector
            integer             :: LinkFirstElem, LinkLastElem
            integer, pointer    :: upBranchIdx, dnBranchIdx
            integer, pointer    :: eIdx, fLidx

            character(64) :: subroutine_name = 'init_network_map_nJ2'
        !%-----------------------------------------------------------------
        !% Preliminaries
        if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------

        !% --- initialize selecteros for upstream and downstream branches
        upBranchSelector = zeroI
        dnBranchSelector = zeroI

        upBranchSelector = upBranchSelector + oneI
        !% --- pointer to upstream branch
        upBranchIdx => node%I(thisJNode,ni_idx_base1 + upBranchSelector)

        !% --- condition for a link included in this partition.
        if (link%I(upBranchIdx,li_P_image) == image) then

            !% --- find the last element index of the link
            LinkLastElem = link%I(upBranchIdx,li_last_elem_idx)

            !% --- find the downstream face index of that last element
            fLidx => node%I(thisJNode,ni_face_idx)

            !% --- if the face is a shared face across images,
            !%     it will not have any upstream local element
            if ( .not. faceYN(fLidx,fYN_isSharedFace)) then
                !% --- local d/s face map to element u/s of the branch
                elemI(LinkLastElem,ei_Mface_dL) = fLidx

                !% --- local u/s element of the face
                faceI(fLidx,fi_Melem_uL) = LinkLastElem
            end if
        end if

        dnBranchSelector = dnBranchSelector + oneI
        !% --- pointer to upstream branch
        dnBranchIdx => node%I(thisJNode,ni_idx_base2 + dnBranchSelector)

        !% --- condition for a link included in this partition.
        if (link%I(dnBranchIdx,li_P_image) == image) then

            !% --- find the first element index of the link
            LinkFirstElem = link%I(dnBranchIdx,li_first_elem_idx)

            !% --- find the downstream face index of that last element
            fLidx => node%I(thisJNode,ni_face_idx)

            !% --- if the face is a shared face across images,
            !%     it will not have any downstream local element
            !%     (not sure if we need this condition)
            if ( .not. faceYN(fLidx,fYN_isSharedFace)) then

                !% --- local map to upstream face for elemI
                elemI(LinkFirstElem,ei_Mface_uL) = fLidx

                !% --- local downstream element of the face
                faceI(fLidx,fi_Melem_dL) = LinkFirstElem
            end if
        end if

        !%-----------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_network_map_nJ2
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_network_map_shared_nJ2_nodes (image, fLidx, nIdx)
        !%-----------------------------------------------------------------
        !% Description:
        !% set the global index, map, and ghost element for nJ2 nodes
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: image, fLidx, nIdx

            integer             :: ii
            integer, pointer    :: fGidx, eUp, eDn, targetImage
            logical, pointer    :: isUpGhost, isDnGhost

            character(64) :: subroutine_name = 'init_network_map_shared_nJ2_nodes'
        !%-----------------------------------------------------------------
        !% Preliminaries
        if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------
        !% Aliases
            fGidx       => faceI(fLidx,fi_Gidx)
            eUp         => faceI(fLidx,fi_Melem_uL)
            eDn         => faceI(fLidx,fi_Melem_dL)
            targetImage => faceI(fLidx,fi_Connected_image)
        !%-----------------------------------------------------------------

        do ii = 1,N_face(targetImage)

            if ((faceI(ii,fi_Connected_image)[targetImage]   == image) .and. &
                (faceI(ii,fi_node_idx_BIPquick)[targetImage] == nIdx)) then
                
                !% --- set the local index of the indetical face
                faceI(ii,fi_Identical_Lidx)[targetImage] = fLidx

                !% --- find the local ghost element index of the connected image
                if (faceYN(ii,fYN_isUpGhost)[targetImage]) then
                    faceI(ii,fi_GhostElem_uL)[targetImage]   = eUp
                elseif (faceYN(ii,fYN_isDnGhost)[targetImage]) then
                    faceI(ii,fi_GhostElem_dL)[targetImage] = eDn
                end if

                !% --- find the global index and set to target image
                if (faceI(ii,fi_Gidx)[targetImage] == nullvalueI) then
                    faceI(ii,fi_Gidx)[targetImage] = fGidx
                end if
            end if
        end do

        !%-----------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_network_map_shared_nJ2_nodes
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_network_nullify_nJm_branch (ElemIdx, FaceIdx)
        !%-----------------------------------------------------------------
        !% Description
        !% set all the values to zero for a null junction
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in)  :: ElemIdx, FaceIdx
            character(64) :: subroutine_name = 'init_network_nullify_nJm_branch'
        !%-----------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------

        !% --- set everything to zero for a non-existant branch
        elemR(ElemIdx,:)                            = zeroR
        elemSR(ElemIdx,:)                           = zeroR
        elemSGR(ElemIdx,:)                          = zeroR
        elemSI(ElemIdx,esi_JunctionBranch_Exists)   = zeroI
        faceR(FaceIdx,:)                            = zeroR
        faceYN(FaceIdx,fYN_isnull)                  = .true.
        elemYN(ElemIdx,eYN_isDummy)                 = .true.

        !%-----------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_network_nullify_nJm_branch
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_network_set_interior_faceYN ()
        !%-----------------------------------------------------------------
        !% Description
        !% set the logicals of fYN_isInteriorFace
        !%-----------------------------------------------------------------
        !% Declarations:
            character(64) :: subroutine_name = 'init_network_set_interior_faceYN'
        !--------------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------

        where ( (faceI(:,fi_BCtype)         ==  doesnotexist) &
                .and. &
                (faceYN(:,fYN_isnull)       .eqv. .false.     ) &
                .and. &
                (faceYN(:,fYN_isSharedFace) .eqv. .false.     ) )

            faceYN(:,fYN_isInteriorFace) = .true.
        endwhere

        !%-----------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_network_set_interior_faceYN
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_network_identify_boundary_element()
        !%-----------------------------------------------------------------
        !% Description:
        !% Identify and mark the element those are connected to a shared face
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, pointer :: eUp(:), eDn(:)
            character(64)    :: subroutine_name = 'init_network_identify_boundary_element'
        !%-----------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------
        
        eUp => faceI(:,fi_Melem_uL)
        eDn => faceI(:,fi_Melem_dL)

        where (faceYN(:,fYN_isUpGhost))
            elemYN(eDn,eYN_isBoundary_up) = .true.
        endwhere

        where (faceYN(:,fYN_isDnGhost))
            elemYN(eUp,eYN_isBoundary_dn) = .true.
        endwhere 

        !%-----------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_network_identify_boundary_element
!%
!%==========================================================================
!%==========================================================================
!%
    ! function init_network_nJm_branch_length (LinkIdx) result (BranchLength)
    !%    ARCHIVE FOR FUTURE USE
    !     !--------------------------------------------------------------------------
    !     !
    !     !% compute the length of a junction branch
    !     !
    !     !--------------------------------------------------------------------------

    !     integer, intent(in)  :: LinkIdx
    !     real(8)              :: BranchLength
    !     real(8), pointer     :: elem_nominal_length, elem_shorten_cof

    !     character(64) :: subroutine_name = 'init_network_nJm_branch_length'
    !     !--------------------------------------------------------------------------
    !     if (setting%Debug%File%network_define) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !     elem_nominal_length => setting%Discretization%NominalElemLength
    !     elem_shorten_cof    => setting%Discretization%JunctionBranchLengthFactor

    !     !% find the length of the junction branch
    !     if (link%I(LinkIdx,li_length_adjusted) == OneSideAdjust) then
    !         BranchLength = link%R(LinkIdx,lr_Length) - link%R(LinkIdx,lr_AdjustedLength)
    !     elseif (link%I(LinkIdx,li_length_adjusted) == BothSideAdjust) then
    !         BranchLength = (link%R(LinkIdx,lr_Length) - link%R(LinkIdx,lr_AdjustedLength))/twoR
    !     elseif (link%I(LinkIdx,li_length_adjusted) == DiagAdjust) then
    !         BranchLength = elem_shorten_cof * elem_nominal_length
    !     end if

    !     if (setting%Debug%File%network_define) &
    !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    ! end function init_network_nJm_branch_length
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
!%
end module network_define

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
    use array_index
    use initialization
    use data_keys
    use globals
    use junction  

    implicit none
    
    private
    
    public :: network_initiation
    
    integer:: debuglevel = 0
    
 contains
!
!========================================================================== 
!==========================================================================
!
 subroutine network_initiation &
    (linkI, linkR, linkYN, linkName, &
     nodeI, nodeR, nodeYN, nodeName, &
     elem2R, elem2I, elem2YN, elem2Name, &
     elemMR, elemMI, elemMYN, elemMName, &
     faceR,  faceI,  faceYN,  faceName)
!
! Initializes a element-face network from a link-node network.
! Requires network links and nodes before execution (e.g. module test_cases)
! 
 character(64) :: subroutine_name = 'network_initiation'
 
 integer,   target,         intent(in out)      :: linkI(:,:)
 integer,   target,         intent(in out)      :: nodeI(:,:)
 real,                      intent(in out)      :: linkR(:,:)    
 real,      target,         intent(in out)      :: nodeR(:,:)  
 
 type(string), intent(in out)   :: linkName(:) 
 type(string), intent(in out)   :: nodeName(:)
 
 logical,   intent(in out)    :: linkYN(:,:)
 logical,   intent(out)       :: nodeYN(:,:)
 
!%   elem2# are the values for elements that have only 2 faces
 real,       dimension(:,:), allocatable, target    :: elem2R       ! real data for elements with 2 faces
 integer,    dimension(:,:), allocatable, target    :: elem2I       ! integer data for elements with 2 faces
 logical,    dimension(:,:), allocatable, target    :: elem2YN      ! logical data for elements with 2 faces

 type(string), dimension(:), allocatable, target    :: elem2Name    ! array of character strings

!%   elemM# are the values for elements that have more than 2 faces
 real,       dimension(:,:), allocatable, target    :: elemMR       ! real data for elements with multi faces
 integer,    dimension(:,:), allocatable, target    :: elemMI       ! integer data for elements with multi faces
 logical,    dimension(:,:), allocatable, target    :: elemMYN      ! logical data for elements with multi faces

 type(string), dimension(:), allocatable, target    :: elemMName    ! array of character strings

!%   face# are the values for faces (always bounded by 2 elements)
 real,       dimension(:,:), allocatable, target    :: faceR       ! real data for faces
 integer,    dimension(:,:), allocatable, target    :: faceI       ! integer data for faces
 logical,    dimension(:,:), allocatable, target    :: faceYN      ! logical data for face

 type(string), dimension(:), allocatable, target    :: faceName    ! array of character strings
 
 integer :: ii

!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

!%   assign links to nodes
 call network_node_assignment (nodeI, linkI) 
 
!%   confirm that the link-node network is valid   
 call network_check_node_link_match (linkI, nodeI)

!%   check that sufficient BC locations have been identified 
 call network_check_BC (nodeI, N_node)
     
!%   get the slope of each link given the node Z values
 call network_get_link_slope (linkR, nodeR, linkI, nodeI)

! HACK - need a check here to look for non-monotonic zbottom in the link/node
! system. Where found the system needs to be subdivided into monotonic links.
! Typically this should only be an issue where the links are representing a
! high-resolution natural channel. 

!%   add sections of links to the nodes to create junctions 
 call network_adjust_link_length (linkR, nodeR, linkI, nodeI)      

!%   determine the system size 
 call network_count_elements_and_faces (nodeR, linkI, nodeI)

!%   allocate all the data storage 
 call allocate_data_storage &
    (elem2R, elemMR, faceR, elem2I, elemMI, faceI, elem2YN, elemMYN, faceYN, &
     elem2Name, elemMname, faceName)
 

!%   ensure elemMR (junction) array has zero values in the multiple geometry storage
 call initialize_array_zerovalues (elemMR) 
 call initialize_dummy_values &
    (elem2R, elem2I, elem2YN, &
     elemMR, elemMI, elemMYN, &
     faceR,  faceI,  faceYN)

!%   store the data from the network in the element and face arrays
 call network_data_create &
    (elem2R, elemMR, faceR, linkR, nodeR, elem2I, elemMI, faceI, &
     linkI, nodeI, nodeYN, elem2Name, elemMName, faceName, linkName, nodeName)   

!%   setup the geometric relationships between the junction branches and the main values
 call junction_geometry_setup (elemMR, elemMI) 
 
!%   assign branch mappings for faces 
 call junction_branch_assigned_to_faces (faceI, elemMI)
  

!% Debug output
 if ((debuglevel > 0) .or. (debuglevelall > 0)) then
    print *, subroutine_name,'----------------------------------------'
    print *, first_elem2_index, first_elem2_index+N_elem2-1
    ! Printout network for debugging
    print *
    print *, 'a)       ii  ,   e2i_idx  ,  e2i_elem_type  '
    do ii=first_elem2_index, first_elem2_index+N_elem2-1
        print *, ii, elem2I(ii,e2i_idx), elem2I(ii,e2i_elem_type)
    enddo

    print *
    print *, 'b)       ii  ,e2i_geometry, e2i_link_ID, e2i_link_Pos'
    do ii=first_elem2_index, first_elem2_index+N_elem2-1
        print *, ii, elem2I(ii,e2i_geometry), elem2I(ii,e2i_link_ID), elem2I(ii,e2i_link_Pos)
    enddo
    print *
    print *, 'c)       ii  ,e2i_roughness_type, e2i_link_ID, e2i_link_Pos'
    do ii=first_elem2_index, first_elem2_index+N_elem2-1
        print *, ii, elem2I(ii,e2i_roughness_type), elem2I(ii,e2i_link_ID), elem2I(ii,e2i_link_Pos)
    enddo

    print *
    print *, 'd)       ii  ,e2i_Mface_u, e2i_Mface_d,  e2i_link_ID, e2i_link_Pos'
    do ii=first_elem2_index, first_elem2_index+N_elem2-1
        print *, ii, elem2I(ii,e2i_Mface_u), elem2I(ii,e2i_Mface_d), elem2I(ii,e2i_link_ID), elem2I(ii,e2i_link_Pos)
    enddo
    
    print *
    print *, 'e)       ii  ,e2r_Length, e2r_Topwidth, e2r_Zbottom'
    do ii=first_elem2_index, first_elem2_index+N_elem2-1
        print *, ii, elem2R(ii,e2r_Length), elem2R(ii,e2r_Topwidth), elem2R(ii,e2r_Zbottom)
    enddo
    
    print *
    print *, '------------- faces -----------------expecting ',N_face
    print *
    print *, 'h)       ii,         idx,     Melem_u,    Melem_d,    etype_u,    etype_d'
    do ii=first_face_index, first_face_index+N_face-1
        print *, ii, faceI(ii,fi_idx), faceI(ii,fi_Melem_u), faceI(ii,fi_Melem_d), &
                     faceI(ii,fi_etype_u), faceI(ii,fi_etype_d)
    enddo
    
    print * 
    print *, 'i)       ii,         Zbottom,     Topwidth'
    do ii=first_face_index, first_face_index+N_face-1
        print *, ii, faceR(ii,fr_Zbottom), faceR(ii,fr_Topwidth)
    enddo
    
    print *
    print *, '------------- junctions ---------------------'
    print *
    print *, 'j)       ii  , face maps u1, u2, d1 '
    do ii= first_elemM_index, first_elemM_index+N_elemM-1
        print *, ii, elemMI(ii,eMi_Mface_u1), elemMI(ii,eMi_Mface_u2), elemMI(ii,eMi_Mface_d1) 
    enddo
    
    print *
    print *, 'k)       ii  , Length, u1,  u2, d1, '
    do ii= first_elemM_index, first_elemM_index+N_elemM-1
        print *, ii, elemMR(ii,eMr_Length), elemMR(ii,eMr_Length_u1), elemMR(ii,eMr_Length_u2), elemMR(ii,eMr_Length_d1) 
    enddo

    print *
    print *, 'l)       ii  , Topwidth, u1, u2, d1 '
    do ii= first_elemM_index, first_elemM_index+N_elemM-1
        print *, ii, elemMR(ii,eMr_Topwidth), elemMR(ii,eMr_Topwidth_u1), elemMR(ii,eMr_Topwidth_u2), elemMR(ii,eMr_Topwidth_d1) 
    enddo
    
    print *
    print *, 'm)       ii  , Zbottom, u1, u2, d1 '
    do ii= first_elemM_index, first_elemM_index+N_elemM-1
        print *, ii, elemMR(ii,eMr_Zbottom), elemMR(ii,eMr_Zbottom_u1), elemMR(ii,eMr_Zbottom_u2), elemMR(ii,eMr_Zbottom_d1) 
    enddo
 endif 
 
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine network_initiation
!
!========================================================================== 
!
! PRIVATE BELOW HERE
!
!==========================================================================
!
 subroutine network_node_assignment &
    (nodeI, linkI)
!
! assign maps from nodes to links that are consistent with link maps
! This assumes that all links are assigned with maps to nodes correctly
!
 character(64) :: subroutine_name = 'network_node_assignment'
 
 integer,   intent(in out)    :: nodeI(:,:)
 integer,   intent(in)        :: linkI(:,:)
 
 integer :: ii, thisnode
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
! assign zeros for accumulators
 nodeI(:,ni_N_link_d) = 0    
 nodeI(:,ni_N_link_u) = 0    
 
 nodeI(:,ni_Mlink_d1) = nullvalueI
 nodeI(:,ni_Mlink_d2) = nullvalueI
 nodeI(:,ni_Mlink_d3) = nullvalueI 
 nodeI(:,ni_Mlink_u1) = nullvalueI
 nodeI(:,ni_Mlink_u2) = nullvalueI
 nodeI(:,ni_Mlink_u3) = nullvalueI  
  
!%  cycle through links to assign nodes
 do ii=1,N_link
    !%  get the upstream node of the link  
    thisnode = linkI(ii,li_Mnode_u)
    print *, 'this node=', thisnode
    if ((thisnode < 1) .or. (thisnode > N_node)) then
        print *, ii,'= this link'
        print *, thisnode,'= upstream node assigned'
        print *, 'error: problem in link-node definitions'
        stop
    else
        nodeI(thisnode,ni_N_link_d) = nodeI(thisnode,ni_N_link_d) +1     
    endif

    if (nodeI(thisnode,ni_N_link_d) > dnstream_face_per_elemM) then
        print *, dnstream_face_per_elemM,' = allowed downstream links per node'
        print *, ii,'= this link'
        print *, thisnode,'= this node'
        print *, 'error: attempt to assign too many downstream links to one node'
        stop
    endif

    select case (nodeI(thisnode,ni_N_link_d))
        case (1)
            nodeI(thisnode,ni_Mlink_d1) = ii
        case (2)
            nodeI(thisnode,ni_Mlink_d2) = ii
        case (3)
            nodeI(thisnode,ni_Mlink_d3) = ii
        case default
            print *, 'error - attempt to assign more than 3 downstream links to one node'
            stop
    end select
    
!    print *, '----------'
!    print *, ii
!    print *, nodeI(:,ni_Mlink_d1)
!    print *, nodeI(:,ni_Mlink_d2)
!    print *, nodeI(:,ni_Mlink_d3)
!    print *, nodeI(:,ni_N_link_d)
!
    !% get the downstream node of the link
    thisnode = linkI(ii,li_Mnode_d)
    if ((thisnode < 1) .or. (thisnode > N_node)) then
        print *, ii,'= this link'
        print *, thisnode,'= downstream node assigned'
        print *, 'error: problem in link-node definitions'
        stop
    else
        nodeI(thisnode,ni_N_link_u) = nodeI(thisnode,ni_N_link_u) +1     
    endif
    
    if (nodeI(thisnode,ni_N_link_u) > upstream_face_per_elemM) then
        print *,upstream_face_per_elemM,' = allowed upstream links per node'
        print *, ii,'= this link'
        print *, thisnode,'= this node'
        print *, 'error: attempt to assign too many upstream links to one node'
        stop
    endif

    select case (nodeI(thisnode,ni_N_link_u))
        case (1)
            nodeI(thisnode,ni_Mlink_u1) = ii
        case (2)
            nodeI(thisnode,ni_Mlink_u2) = ii
        case (3)
            nodeI(thisnode,ni_Mlink_u3) = ii
        case default
            print *, 'error - attempt to assign more than 3 uptream links to one node'
            stop
    end select   

!    print *, '----------'
!    print *, ii
!    print *, nodeI(:,ni_Mlink_u1)
!    print *, nodeI(:,ni_Mlink_u2)
!    print *, nodeI(:,ni_Mlink_u3)
!    print *, nodeI(:,ni_N_link_u)
!
    
 end do
  
!%  cycle through links to assign nodes
! do ii=1,N_link
! 
!    !%  get the upstream node of the link (if it exists)    
!    thisnode = linkI(ii,li_Mnode_u)
!    
!    print *, ii,'=Link; ', thisnode,'= upstream node assigned'
!    
!    !%  look for next available downstream link position
!    if (thisnode > 0) then
!        if (nodeI(thisnode,ni_Mlink_d1) == nullvalueI) then
!            nodeI(thisnode,ni_Mlink_d1) = ii
!        else
!            print *, 'node ',thisnode
!            print *, 'error: attempt to assign 2 downstream links to 1 node in ',subroutine_name
!            stop
!        endif
!        
!        !%  increment the downstream link counter
!        nodeI(thisnode,ni_N_link_d) = nodeI(thisnode,ni_N_link_d)+1
!    endif
!    
!    !%  get the downstream node of the link (if it exists)
!    thisnode = linkI(ii,li_Mnode_d)
!    
!    print *, ii,'=Link; ', thisnode,'= downstream node assigned'
!    
!    !%  look for the next available upstream link position
!    if (thisnode > 0) then
!
!        if (nodeI(thisnode,ni_Mlink_u1) == nullvalueI) then
!            nodeI(thisnode,ni_Mlink_u1) = ii
!        elseif (nodeI(thisnode,ni_Mlink_u2) == nullvalueI) then
!                nodeI(thisnode,ni_Mlink_u2) = ii
!        else
!            print *, 'node ',thisnode
!            print *, 'error: attempt to assign 3 upstream links to 1 node in ',subroutine_name
!            stop
!        endif
!        
!        !%  increment the upstream link counter
!        nodeI(thisnode,ni_N_link_u) = nodeI(thisnode,ni_N_link_u)+1
!    endif
! enddo 
! 


 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine network_node_assignment
!
!========================================================================== 
!==========================================================================
!
 subroutine network_check_node_link_match &
    (linkI, nodeI)
!
! check that the upstream and downstream linkages are consistent
! and complete
!    
 character(64) :: subroutine_name = 'network_check_node_link_match'

 integer, target,   intent(in) :: linkI(:,:)
 integer,           intent(in) :: nodeI(:,:)

 integer :: ii, mm
 integer, pointer :: Unode, Dnode

 logical :: ifound

!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name
 
 do ii = 1,N_link
 
    !%  upstream and downstream nodes to this link
    Unode => linkI(ii,li_Mnode_u)
    Dnode => linkI(ii,li_Mnode_d)

    !%  check that upstream node has the correct downstream link  
    ifound = .false.
    do mm=1,dnstream_face_per_elemM
        if (nodeI(Unode,ni_MlinkDn(mm)) == ii) then
            ifound = .true.
        endif
    enddo
    if (.not. ifound) then
        print *, 'downstream link ',ii,' not found for node ', Unode
        stop
    endif

    !%  check that downstream node has the correct upstream link
    ifound = .false.
    do mm=1,upstream_face_per_elemM
        if (nodeI(Dnode,ni_MlinkUp(mm)) == ii) then
            ifound = .true.
        endif
    enddo
    if (.not. ifound) then
        print *, 'upstream link ',ii,' not found for node ', Dnode
        stop
    endif

 enddo

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name   
 end subroutine network_check_node_link_match    
!
!==========================================================================
!==========================================================================
!
 subroutine network_check_BC &
    (nodeI, N_node)
!    
! check that BC nodes have correct upstream and downstream links
!
    character(64) :: subroutine_name = 'network_check_BC'
    integer, intent(in)     :: N_node
    integer, intent(in)     :: nodeI(:,:)
    
    integer :: ii
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name
    
 do ii = 1,N_node
    if ( nodeI(ii,ni_node_type) == nBCup ) then
        if ( nodeI(ii,ni_Mlink_d1) < 1 ) then
            print *, 'node = ',ii
            print *, 'error: A downstream link is required at an upstream BC node in ',subroutine_name
            stop
        endif
        
        if ( nodeI(ii,ni_N_link_d) > 1 ) then
            print *, 'node = ',ii  
            print *, 'error: An upstream BC node can only have one downstream link in ',subroutine_name
            STOP
        endif
        
        if ( nodeI(ii,ni_Mlink_u1) > 1 ) then
            print *, 'node = ',ii
            print *, 'link = ',nodeI(ii,ni_Mlink_u1)
            print *, 'error: Upstream BC node cannot have upstream link in ',subroutine_name
            stop
        endif
       
    endif
    
    if ( nodeI(ii,ni_node_type) == nBCdn ) then
        if ( nodeI(ii,ni_Mlink_u1) < 1 ) then
            print *, 'node = ',ii
            print *, 'error: An upstream link is required at downstream BC node in ',subroutine_name
            stop
        endif
        
        if ( nodeI(ii,ni_N_link_u) > 1 ) then
            print *, 'node = ',ii  
            print *, 'error: A downstream BC node can only have one upstream link in ',subroutine_name
            STOP
        endif

        if ( nodeI(ii,ni_Mlink_d1) > 1 ) then
            print *, 'node = ',ii
            print *, 'link = ',nodeI(ii,ni_Mlink_d1)
            print *, 'error: Downstream BC node cannot have downstream link in ',subroutine_name
            stop
        endif
        
    endif
 enddo
    
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name   
 end subroutine network_check_BC
!
!==========================================================================
!==========================================================================
!
 subroutine network_get_link_slope &
    (linkR, nodeR, linkI, nodeI)
!
! compute the slope across each link
! 
 character(64) :: subroutine_name = 'network_get_link_slope'
 
 integer,   target,         intent(in)      :: linkI(:,:)
 integer,   target,         intent(in)      :: nodeI(:,:)
 real,                      intent(in out)  :: linkR(:,:)    
 real,      target,         intent(in)      :: nodeR(:,:)

 integer,   pointer :: nUp, nDn
 real,      pointer :: zUp, zDn
 integer :: mm
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 

 do mm = 1,N_link  
    nUp => linkI(mm,li_Mnode_u)
    nDn => linkI(mm,li_Mnode_d)
    zUp => nodeR(nUp,nr_Zbottom)
    zDn => nodeR(nDn,nr_Zbottom)
    linkR(mm,lr_Slope) = (zUp - zDn) / linkR(mm,lr_Length)
 end do
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) then
    !% provide output for debugging
    print *, subroutine_name,'--------------------------------'
    print *, 'link ID,  Slope, length'
    do mm=1,N_link
        print *, mm, linkR(mm,lr_Slope), linkR(mm,lr_Length)
    end do 
 endif

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine network_get_link_slope
!
!========================================================================== 
!==========================================================================
!
 subroutine network_adjust_link_length &
    (linkR, nodeR, linkI, nodeI )
!    
! Adjust the link lengths so that a portion goes to each junction
! when the junctions need physical dimensions. This is necessary when
! a junction has more than an single upstream and downstream reach.
!
! We take away a portion of a nominal element length at a junction if
! there is more than one element in the upstream or downstream link.
! If there is only one element in the link, then subtract a portion of 
! the actual link length. 
!
! The "lost" link length is added to the junction (node)
!
 character(64) :: subroutine_name = 'network_adjust_link_length'

 integer,   target,         intent(in out)  :: linkI(:,:)
 integer,                   intent(in)      :: nodeI(:,:)
 real,      target,         intent(in out)  :: linkR(:,:)    
 real,      target,         intent(in out)  :: nodeR(:,:)

 integer :: ii, mm

 real ::  delta
 
 real,    pointer   :: linkLength(:)
 integer, pointer   :: linkNelem(:)
 
 real,    pointer   :: element_nominal_length(:)
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name  

 linkLength             => linkR(:,lr_Length)
 element_nominal_length => linkR(:,lr_ElementLength)
 linkNelem              => linkI(:,li_N_element)

 do ii = 1,N_node 
    !%   only applies to nJm junctions -- not to nJ2 junctions
    if (nodeI(ii,ni_node_type) == nJm) then    
        !%   shorten the upstream links
        if (nodeI(ii,ni_N_link_u) > 0) then
            call link_shortening &
                (linkLength, nodeR, nodeI, ii, &
                 ni_N_link_u, ni_MlinkUp, nr_ElementLengthUp, element_nominal_length)
        endif
        
        !%   shorten the downstream links
        if (nodeI(ii,ni_N_link_d) > 0) then
            call link_shortening &
                (linkLength, nodeR, nodeI, ii,  &
                 ni_N_link_d, ni_MlinkDn, nr_ElementLengthDn, element_nominal_length)
        endif
    endif
 enddo

 
!%  get the uniform elements that subdivide a link
 do ii = 1,N_link
        
    linkNelem(ii) = floor(linkLength(ii)/element_nominal_length(ii))
    
    if (linkNelem(ii) > 0) then
        delta = ( linkLength(ii) - element_nominal_length(ii) * linkNelem(ii)) &
              / linkNelem(ii)
        linkR(ii,lr_ElementLength) = element_nominal_length(ii) + delta
    else
        !%  for case when only 1 nominal element will fit in link
        linkNelem(ii) = 1
        linkR(ii,lr_ElementLength) = linkLength(ii)
    endif
 enddo
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) then
    !%  debugging output
    print *
    print *, subroutine_name,'----------------------------------'
    print *, 'link number, element_length, N_element'
    do ii=1,N_link
        print *, ii, linkR(ii, lr_ElementLength), linkI(ii, li_N_element)
    enddo
 endif
    
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine network_adjust_link_length
!
!==========================================================================
!==========================================================================
!
 subroutine network_count_elements_and_faces &
    (nodeR,  linkI, nodeI)
!
! computes global values for element, face, and BC counts
!
 character(64) :: subroutine_name = 'network_count_elements_and_faces'

 real,    intent(in)   :: nodeR(:,:)
 integer, intent(in)   :: linkI(:,:)
 integer, intent(in)   :: nodeI(:,:)      

 integer :: ii, NupstreamBC, NdownstreamJ, NdnstreamBC

!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!% sum the number of elements in the links
 N_elem2 = sum(linkI(:,li_N_element), MASK=(linkI(:,li_N_element)> 0))  

!% add junctions that have physical dimesions - these become elements
 do ii=1,N_node
    if ( (nodeR(ii,nr_ElementLength_u1) > 0.0) .or. &
         (nodeR(ii,nr_ElementLength_d1) > 0.0) ) then
         N_elemM = N_elemM +1
    endif
 enddo
 
!% count the number of upstream BC nodes
 NupstreamBC = count(nodeI(:,ni_node_type) == nBCup)
 
!% count the number of upstream BC nodes
 NdnstreamBC = count(nodeI(:,ni_node_type) == nBCdn)
 
!% Add ghost elements outside of BC
 N_elem2 = N_elem2 + NupstreamBC + NdnstreamBC  

!% count the numer of downstream faces at junction elements
 NdownstreamJ = sum(nodeI(:,ni_N_link_d), MASK = (nodeI(:,ni_node_type) == nJm))

!% Faces are between elem2, upstream BC, and all downstreamJ.
!% Note the NdnstreamBC are already counted in the initial link count.
 N_face = sum(linkI(:,li_N_element)) + NupstreamBC + NdownstreamJ
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) then
    !%  debugging output
    print *
    print *, subroutine_name,'----------------------------------' 
    print *, 'expecting # elem2 = ',N_elem2
    print *, 'expecting # elemM = ',N_elemM
    print *, 'expecting # face  = ',N_face
 endif

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name    
 end subroutine network_count_elements_and_faces   
!
!==========================================================================  
!==========================================================================
!
 subroutine network_data_create & 
    (elem2R, elemMR, faceR, linkR, nodeR, elem2I, elemMI, faceI, &
     linkI, nodeI, nodeYN, elem2Name, elemMName, faceName, linkName, nodeName)
 !
 ! creates the network of elements and junctions from nodes and links
 !
 character(64) :: subroutine_name = 'network_data_create'

 integer, intent(in out) :: elem2I(:,:), elemMI(:,:), faceI(:,:)
 real,    intent(in out) :: elem2R(:,:), elemMR(:,:), faceR(:,:)

 integer, intent(in out) :: linkI(:,:)
 
 integer, target, intent(in out) :: nodeI(:,:)
 
 real,    target, intent(in)     :: linkR(:,:), nodeR(:,:)
 
 logical, target, intent(in) :: nodeYN(:,:)

 integer :: ii, lastElem2, thisElem2, lastElemM, thisElemM, lastFace, thisFace

 integer :: N_BCnodes_d, thisNode, thisLink
 
 logical, pointer :: nodeYNmask(:)
 
 integer, pointer :: nodesDownstream(:)
 
 real, pointer  :: Zdownstream
 
 type(string),  intent(in out)  :: elem2Name(:), elemMName(:), faceName(:)
 type(string),  intent(in)      :: linkName(:), nodeName(:)

!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name
 
 lastElem2 = zeroI ! last elem2 that was assigned
 lastElemM = zeroI ! last elemM that was assigned
 thisElem2 = first_elem2_index ! this elem2 that is ready to be assigned
 thisElemM = first_elemM_index ! this elemM that is ready to be assigned
 lastFace  = zeroI
 thisFace  = first_face_index
    
 nodeYNmask => nodeYN(:,nYN_temp1)
 nodeYNmask =  (nodeI(:,ni_node_type) == nBCdn)
 
 nodesDownstream => nodeI(:,ni_temp1)
 N_BCnodes_d = count(nodeYNmask)
 
 if (N_BCnodes_d == 0) then
    print *, 'No downstream channel boundary nodes found - code development incomplete for dowstream pipe'
    STOP
 elseif (N_BCnodes_d > 0) then
    !NOTE - need to test/modify this for multiple downstream nodes.
    ! get a list of downstream nodes
    nodesDownstream(1:N_BCnodes_d) = pack(nodeI(:,ni_idx),nodeYNmask)
 else
    print *, 'N_BCnodes_d ',N_BCnodes_d
    print *, 'error: unexpected else in ',subroutine_name
    stop
 endif

 do ii=1,N_BCnodes_d

    if ((debuglevel > 0) .or. (debuglevelall > 0)) then
        print *, 'starting downstream node ', nodesDownstream(ii)
    endif
    
    thisNode = nodesDownstream(ii)
    if (nodeI(thisNode,ni_assigned) == nAssigned) then
        cycle ! node is already assigned!
    endif
    nodeI(thisNode,ni_assigned) = setAssigned(thisNode, ni_assigned, nUnassigned, nAssigned, nodeI)
    
    !% error checking for downstream BC
    if ( nodeI(thisNode,ni_N_link_d) /= 0 ) then
        print *, 'error: A downstream BC node should not have a downstream link in ',subroutine_name
        STOP
    endif        
    if ( nodeI(thisNode,ni_N_link_u) /= 1 ) then
        print *, 'error: A downstream BC node must have exactly one upstream link in ',subroutine_name
        STOP
    endif

    !% set the upstream link (one only) of the downstream BC node
    thisLink = nodeI(thisNode,ni_Mlink_u1)
    if (linkI(thisLink,li_assigned) == lAssigned) then
        print *, 'error: it is not clear how the node was not assigned but the link is assigned in ',subroutine_name
        stop
    endif
    linkI(thisLink,li_assigned) = setAssigned(thisLink, li_assigned, lUnassigned, lAssigned, linkI)

    !% Assign the ghost element - use values for upstream link
    elem2I(thisElem2,e2i_idx)            = thisElem2
    elem2I(thisElem2,e2i_elem_type)      = eBCdn
    elem2I(thisElem2,e2i_geometry)       = linkI(thisLink,li_geometry)
    elem2I(thisElem2,e2i_roughness_type) = linkI(thisLink,li_roughness_type)
    elem2I(thisElem2,e2i_link_ID)        = thisLink
    elem2I(thisElem2,e2i_link_Pos)       = nullvalueI
    elem2I(thisElem2,e2i_Mface_u)        = thisFace
    elem2I(thisElem2,e2i_Mface_d)        = nullvalueI
    
    elem2R(thisElem2,e2r_Length)         = linkR(thislink,lr_ElementLength)
    
    !% note the following has a minus as we are going downstream for the ghost
    elem2R(thisElem2,e2r_Zbottom)        = nodeR(thisNode,nr_Zbottom)          &
                                            - 0.5 * linkR(thisLink,lr_Slope)   &
                                            * linkR(thislink,lr_ElementLength)
    
    select case (linkI(thisLink,li_geometry))
        case (lRectangular)
            elem2R(thisElem2,e2r_BreadthScale) = linkR(thisLink,lr_BreadthScale)
            elem2R(thisElem2,e2r_Topwidth)     = linkR(thisLink,lr_BreadthScale)
            !faceR(thisFace,fr_Topwidth)    = linkR(thisLink,lr_Breadth)
        case (lParabolic)
            elem2R(thisElem2,e2r_BreadthScale) = zeroR
            elem2R(thisElem2,e2r_Topwidth) = twoR &
                * sqrt(linkR(thisLink,lr_InitialDepth)/linkR(thisLink,lr_ParabolaValue))
        case (lTrapezoidal)
            elem2R(thisElem2,e2r_BreadthScale) = linkR(thisLink,lr_BreadthScale)
            elem2R(thisElem2,e2r_Topwidth) = linkR(thisLink,lr_BreadthScale)   &
                    + linkR(thisLink,lr_InitialDepth)                          &
                    * (linkR(thisLink,lr_LeftSlope) + linkR(thisLink,lr_RightSlope))
        case (lTriangle)
            elem2R(thisElem2,e2r_BreadthScale) = zeroR
            elem2R(thisElem2,e2r_Topwidth)     = linkR(thisLink,lr_InitialDepth) &
                    * (linkR(thisLink,lr_LeftSlope) + linkR(thisLink,lr_RightSlope))
        case (lWidthDepth)
            elem2R(thisElem2,e2r_Topwidth)     = linkR(thisLink,lr_Topwidth)
            elem2R(thisElem2,e2r_BreadthScale) = linkR(thisLink,lr_BreadthScale)
            faceR(thisFace,fr_Topwidth)        = linkR(thisLink,lr_Topwidth)
        case default
            print *, 'error: case statement is incomplete in ',subroutine_name
            stop
    end select
 
    ! use the node name for the ghost element
    elem2Name(thisElem2) = nodeName(thisNode)
 
    ! store the face info for the downstream BC node
    faceI(thisFace,fi_idx)          = thisFace
    faceI(thisFace,fi_type)         = fBCdn
    faceI(thisFace,fi_Melem_u)      = thisElem2+1         
    faceI(thisFace,fi_Melem_d)      = thisElem2
    faceI(thisFace,fi_etype_u)      = linkI(thisLink,li_link_type)
    faceI(thisFace,fi_etype_d)      = elem2I(thisElem2,e2i_elem_type)
    faceI(thisFace,fi_jump_type)    = nullvalueI
    faceI(thisFace,fi_link_ID)      = nullvalueI
    faceI(thisFace,fi_link_Pos)     = nullvalueI
    faceI(thisFace,fi_node_ID)      = thisNode
    
    faceR(thisFace,fr_Zbottom)      = nodeR(thisNode,nr_Zbottom)
    
    faceName(thisFace)              = nodeName(thisNode)
    
    lastFace = thisFace
    thisFace = thisFace+1
    
    lastElem2 = thisElem2
    thisElem2 = thisElem2+1
    
    !% subdivide the first link upstream of the downstream BC (only one link allowed)
    Zdownstream => nodeR(thisNode,nr_Zbottom)
    call subdivide_link_going_upstream &
        (lastElem2, thisElem2, lastFace, thisFace, thisLink,   &
         elem2I, faceI, linkI, elem2R, faceR, linkR, nodeR, Zdownstream) 
    nullify(Zdownstream)     
           
    !% advance the node counter
    thisNode = linkI(thisLink,li_Mnode_u)
    if (nodeI(thisNode,ni_assigned) == nAssigned) then
        cycle ! node is already assigned!
    endif
    nodeI(thisNode,ni_assigned) = setAssigned(thisNode, ni_assigned, nUnassigned, nAssigned, nodeI)
     
    !% deal with the next node 
    call handle_thisnode &
         (lastElem2, thisElem2, lastElemM, thisElemM,          &
          lastFace,  thisFace,  thisNode,  thisLink,           &
          elem2I, elemMI, faceI, linkI, nodeI, elem2R, elemMR, faceR, linkR, nodeR, &
          elem2Name, elemMName, faceName, linkName, nodeName)        
 enddo
    
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine network_data_create

!==========================================================================
!==========================================================================
!
 recursive subroutine handle_thisnode &
     (lastElem2, thisElem2, lastElemM, thisElemM,          &
      lastFace,  thisFace,  thisNode,  thisLink,           &
      elem2I, elemMI, faceI, linkI, nodeI, elem2R, elemMR, faceR, linkR, nodeR, &
      elem2Name, elemMName, faceName, linkName, nodeName)   
!
! Subroutine to handle a node, identify its links, and create elements
! Called recursively in a crawl through the node/link system
!      
 character(64) :: subroutine_name = 'handle_thisnode'
 integer, intent(in out) :: lastElem2, thisElem2, lastElemM, thisElemM, lastFace, thisFace
 integer, intent(in out) :: thisNode,  thisLink

 integer, intent(in out) :: elem2I(:,:), elemMI(:,:), faceI(:,:)

 integer,               intent(in out) :: linkI(:,:)
 integer, target,       intent(in out) :: nodeI(:,:)
 
 real,    intent(in out) :: elem2R(:,:), elemMR(:,:), faceR(:,:)
 
 real,  target,    intent(in)     :: linkR(:,:),  nodeR(:,:)
 
 type(string),  intent(in out)  :: elem2Name(:), elemMName(:), faceName(:)
 type(string),  intent(in)      :: linkName(:), nodeName(:)

 integer ::  ii, mm, jElem
 
 integer, pointer :: dlink
 
 real :: zDownstream

 integer, dimension(upstream_face_per_elemM) :: linkSet
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 if (nodeI(thisNode,ni_node_type) == nBCup) then 

    if ((debuglevel > 0) .or. (debuglevelall > 0)) then
        print *
        print *, subroutine_name,'----------------------------------------'
        print *, 'this node is nBCup ', thisnode
    endif
    
    ! note that only one downstream link is allowed to an upstream BC
    if (nodeI(thisNode,ni_N_link_d) > 1) then
        print *, 'error: only 1 downstream link allowed on an upstream BC node in ',subroutine_name
        stop
    endif
 
    ! get the downstream link - note there can only be one at an upstream BC
    dlink => nodeI(thisNode,ni_Mlink_d1) 
    
    ! the face that was previously defined when the link was subdivided
    thisFace = linkI(dlink,li_Mface_u) 
 
    ! store the ghost element
    elem2I(thisElem2,e2i_idx)            = thisElem2
    elem2I(thisElem2,e2i_elem_type)      = eBCup
    elem2I(thisElem2,e2i_geometry)       = elem2I(lastElem2,e2i_geometry)
    elem2I(thisElem2,e2i_roughness_type) = elem2I(lastElem2,e2i_roughness_type)
    elem2I(thisElem2,e2i_link_ID)        = dLink
    elem2I(thisElem2,e2i_link_Pos)       = nullvalueI
    elem2I(thisElem2,e2i_Mface_u)        = nullvalueI
    elem2I(thisElem2,e2i_Mface_d)        = thisFace

    elem2R(thisElem2,e2r_Length)         = elem2R(lastElem2,e2r_Length)
    elem2R(thisElem2,e2r_Zbottom)        = nodeR(thisNode,nr_Zbottom)  &
                                         + 0.5 * linkR(dlink,lr_Slope)  &
                                         * elem2R(lastElem2,e2r_Length)
    
    select case (elem2I(thisElem2,e2i_geometry))
        case (eRectangular)
            elem2R(thisElem2,e2r_Topwidth) = elem2R(lastElem2,e2r_Topwidth)
            elem2R(thisElem2,e2r_BreadthScale) = elem2R(lastElem2,e2r_BreadthScale)
            !faceR(thisFace,fr_Topwidth)    = elem2R(lastElem2,e2r_Topwidth)
        case (eParabolic)
            elem2R(thisElem2,e2r_Topwidth) = elem2R(lastElem2,e2r_Topwidth)
            elem2R(thisElem2,e2r_BreadthScale) = elem2R(lastElem2,e2r_BreadthScale)
            faceR(thisFace,fr_Topwidth)    = elem2R(lastElem2,e2r_Topwidth)
        case (eTrapezoidal)
            elem2R(thisElem2,e2r_Topwidth) = elem2R(lastElem2,e2r_Topwidth)
            elem2R(thisElem2,e2r_BreadthScale) = elem2R(lastElem2,e2r_BreadthScale)
            faceR(thisFace,fr_Topwidth)    = elem2R(lastElem2,e2r_Topwidth)
        case (eTriangle)
            elem2R(thisElem2,e2r_Topwidth) = elem2R(lastElem2,e2r_Topwidth)
            elem2R(thisElem2,e2r_BreadthScale) = elem2R(lastElem2,e2r_BreadthScale)
            faceR(thisFace,fr_Topwidth)    = elem2R(lastElem2,e2r_Topwidth)
        case (eWidthDepth)
            elem2R(thisElem2,e2r_Topwidth) = elem2R(lastElem2,e2r_Topwidth)
            elem2R(thisElem2,e2r_BreadthScale) = elem2R(lastElem2,e2r_BreadthScale)
            faceR(thisFace,fr_Topwidth)    = elem2R(lastElem2,e2r_Topwidth)
        case default
            print *, 'error: case statement is incomplete in ',subroutine_name
            stop
    end select
    
    ! use the node name for the ghost element
    elem2Name(thisElem2) = nodeName(thisNode)
    
    faceR(thisFace,fr_Zbottom)      = nodeR(thisNode,nr_Zbottom)
            
    faceI(thisFace,fi_idx)          = thisFace
    faceI(thisFace,fi_type)         = fBCup
    faceI(thisFace,fi_Melem_u)      = thisElem2           
    faceI(thisFace,fi_Melem_d)      = linkI(dlink,li_Melem_u)
    faceI(thisFace,fi_etype_u)      = eBCup
    faceI(thisFace,fi_etype_d)      = elem2I(linkI(dlink,li_Melem_u),e2i_elem_type)
    faceI(thisFace,fi_jump_type)    = nullvalueI
    faceI(thisFace,fi_link_ID)      = nullvalueI
    faceI(thisFace,fi_link_Pos)     = nullvalueI
    faceI(thisFace,fi_node_ID)      = thisNode  
    
    faceName(thisFace)              = nodeName(thisNode)
    
    nullify(dlink)
    
    ! set these values to null as the next step will not be connected
    thisNode = nullvalueI
    lastFace = nullvalueI
    lastElem2 = nullvalueI
    
    ! increment counters
    thisFace = thisFace+1
    thisElem2 = thisElem2 + 1
     
 elseif (nodeI(thisNode,ni_node_type) == nBCdn) then  
 
    print *
    print *, 'this node is nBCdn ', thisnode
    print *, 'error: A downstream BC is not expected during the node cycling in ',subroutine_name
    stop
        
 elseif (nodeI(thisNode,ni_node_type) == nJ2) then

    if ((debuglevel > 0) .or. (debuglevelall > 0)) then
        print *
        print *, subroutine_name,'----------------------------------------'
        print *, 'this node is nJ2 ', thisnode
    endif
     
    thisLink = nodeI(thisNode,ni_Mlink_u1)
    linkI(thisLink,li_assigned) = setAssigned(thisLink, li_assigned, lUnassigned, lAssigned, linkI)

    !%  data that were unassigned for last face
    faceI(lastFace,fi_Melem_u)      = thisElem2  
    faceI(lastFace,fi_etype_u)      = linkI(thisLink,li_link_type)
    faceI(lastFace,fi_node_ID)      = thisNode  
    !faceI(lastFace,fi_idx)          = lastFace
    !faceI(lastFace,fi_Melem_d)      = lastElem2
    !faceI(lastFace,fi_etype_d)      = elem2I(lastElem2,e2i_elem_type)
    !faceI(lastFace,fi_jump_type)    = nullvalueI
    !faceI(lastFace,fi_link_ID)      = nullvalueI
    !faceI(lastFace,fi_link_Pos)     = nullvalueI

    faceR(lastface,fr_Zbottom)      = nodeR(thisNode,nr_Zbottom)
    
    faceI(lastFace,fi_type) = setFaceType &
        (faceI(lastFace,fi_etype_u), faceI(lastFace,fi_etype_d))
        
    faceName(lastFace) = nodeName(thisNode)    

    zDownstream = nodeR(thisNode,nr_Zbottom)    
    call subdivide_link_going_upstream &
        (lastElem2, thisElem2, lastFace, thisFace, thisLink,   &
         elem2I, faceI, linkI, elem2R, faceR, linkR, nodeR, zDownstream ) 
    
    !%  get the next upstream link     
    thisNode = linkI(thisLink,li_Mnode_u)     
         
    if (thisNode > 0) then
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, 'Recursion: handle_thisnode ', thisNode

        call handle_thisnode &
            (lastElem2, thisElem2, lastElemM, thisElemM,          &
             lastFace,  thisFace,  thisNode,  thisLink,           &
             elem2I, elemMI, faceI, linkI, nodeI, elem2R, elemMR, faceR, linkR, nodeR , &
             elem2Name, elemMName, faceName, linkName, nodeName)    
    endif
         
    
 elseif (nodeI(thisNode,ni_node_type) == nJm) then

    if ((debuglevel > 0) .or. (debuglevelall > 0)) then
        print *
        print *, subroutine_name,'----------------------------------------'
        print *, 'this node is nJm ', thisnode !, nodeI(thisNode,ni_N_link_u)
    endif

    linkSet(:) = nullvalueI
    
    !%  assign the link ID for the links connected to the node
    linkSet(1) = nodeI(thisNode,ni_Mlink_u1)
    do ii=2,3
        if (upstream_face_per_elemM > ii-1) then
            if (nodeI(thisNode,ni_N_link_u) > ii-1) then
                linkSet(ii) = nodeI(thisNode,ni_MlinkUp(ii))
            endif
        endif  
    enddo
    if (upstream_face_per_elemM  > 3) then
        print *, 'error: code not designed for upstream_face_per_elem > 3 in ',subroutine_name
        STOP
    endif            
    if (nodeI(thisNode,ni_N_link_u) > 3) then
        print *, 'error: attempt to use more than 3 upstream links in ',subroutine_name
        STOP
    endif   

    print *, 'Junction with downstream of thisLink ',thisLink
    !print *, 'linkSet ',linkSet(:)
    !print *, 'linkAss ',linkI(linkSet(1),li_assigned),linkI(linkSet(2),li_assigned)
    
    !%  Create the junction element we will use here and increment for the next(recursion)
    jElem = thisElemM
    lastElemM = thisElemM
    thisElemM = thisElemM+1
    elemMI(jElem,eMi_idx)             = jElem
    elemMI(jElem,eMi_elem_type)       = eJunctionChannel ! HACK - PIPE NOT HANDLED
    elemMI(jElem,eMi_geometry)        = eRectangular ! HACK - NEED AN APPROACH TO ASSIGN
    elemMI(jElem,eMi_nfaces)          = nodeI(thisNode,ni_N_link_u) + nodeI(thisNode,ni_N_link_d)
    elemMI(jElem,eMi_nfaces_d)        = nodeI(thisNode,ni_N_link_d)
    elemMI(jElem,eMi_nfaces_u)        = nodeI(thisNode,ni_N_link_u)
    elemMI(jElem,eMi_roughness_type)  = nullvalueI            
    elemMI(jElem,eMi_node_ID)         = thisNode
    !elemMI(jElem,eMi_link_Pos)        = nullvalueI 
    
    !%  the following assumes the topwidth and length are functions of branches
    elemMR(jElem,eMr_Topwidth)        = zeroR ! do not use nullvalueR with geometry
    elemMR(jElem,eMr_BreadthScale)    = zeroR ! do not use nullvalueR with geometry
    elemMR(jElem,eMr_Length)          = zeroR ! do not use nullvalueR with geometry  
    elemMR(jElem,eMr_Zbottom)         = nodeR(thisNode,nr_Zbottom)
    
    !%  assign values from adjacent downstream elements/links
    do mm=1,nodeI(thisNode,ni_N_link_d)
        dlink => nodeI(thisNode,ni_MlinkDn(mm))        
        elemMR(jElem,eMr_LengthDn(mm))   =  nodeR(thisNode,nr_ElementLengthDn(mm))
        elemMR(jElem,eMr_ZbottomDn(mm))  =  nodeR(thisNode,nr_Zbottom) &
                                            - 0.5* elemMR(jElem,eMr_LengthDn(mm)) &
                                            * linkR(dlink,lr_Slope)       
        select case (elemMI(jElem,eMi_geometry))
            case (eRectangular)
                elemMR(jElem,eMr_TopwidthDn(mm))     =  linkR(dlink,lr_BreadthScale)        
                elemMR(jElem,eMr_BreadthscaleDn(mm)) =  linkR(dlink,lr_BreadthScale)        
            case default
                print *, 'elemMI(jElem,eMi_geometry) ', elemMI(jElem,eMi_geometry)
                print *, 'error: case statement is incomplete in ',subroutine_name
                stop
        end select
    end do
    
    elemMName(JElem) = nodeName(thisNode) 
    
    nullify(dlink)
    
    !%  assign values from adjacent upstream elements/links
    do mm=1,nodeI(thisNode,ni_N_link_u)
        dlink => nodeI(thisNode,ni_MlinkUp(mm))
        elemMR(jElem,eMr_LengthUp(mm))   = nodeR(thisNode,nr_ElementLengthUp(mm))
        elemMR(jElem,eMr_ZbottomUp(mm))  = nodeR(thisNode,nr_Zbottom) &
                                            + 0.5* elemMR(jElem,eMr_LengthUp(mm)) &
                                            * linkR(dlink,lr_Slope)           
        select case (elemMI(jElem,eMi_geometry))
            case (eRectangular)
                elemMR(jElem,eMr_TopwidthUp(mm))     =  linkR(dlink,lr_BreadthScale)        
                elemMR(jElem,eMr_BreadthScaleUp(mm)) =  linkR(dlink,lr_BreadthScale)        
            case default
                print *, 'error: case statement is incomplete in ',subroutine_name
                stop
        end select
    end do
    nullify(dlink)    
        
    !%  Check for downstream elements that have been assigned and store mappings
    do mm=1,nodeI(thisNode,ni_N_link_d)
        dlink => nodeI(thisNode,ni_MlinkDn(mm))
        if (linkI(dlink,li_assigned) == lAssigned) then
            !print *, dlink
            elemMI(jElem,eMi_MfaceDn(mm)) = linkI(dlink,li_Mface_u)
            faceI(linkI(dlink,li_Mface_u), fi_Melem_u) = jElem
            faceI(linkI(dlink,li_Mface_u), fi_etype_u) = elemMI(jElem,eMi_elem_type)
            faceI(linkI(dlink,li_Mface_u), fi_type) = setFaceType &
                 (faceI(linkI(dlink,li_Mface_u),fi_etype_u), &
                  faceI(linkI(dlink,li_Mface_u),fi_etype_d) )
        endif
    enddo
    nullify(dlink)
    
    !%  remove already-assigned links from the set to be cycled
    do mm=1,upstream_face_per_elemM
        if (linkSet(mm) > 0) then
            if (linkI(linkSet(mm),li_assigned) == lAssigned) then
                linkSet(mm) = 0
            endif
        endif
    enddo
             
    !%  assign the faces around the junction for upstream faces
    do mm=1,upstream_face_per_elemM
        dlink => nodeI(thisNode,ni_MlinkUp(mm))
        if (linkSet(mm) > 0) then   
            elemMI(jElem,eMi_MfaceUp(mm))   = thisFace
            faceI(thisFace,fi_idx)          = thisFace
            faceI(thisFace,fi_type)         = nullvalueI ! assigned later
            faceI(thisFace,fi_Melem_u)      = nullvalueI ! upstream element assigned later, before subdivide_link
            faceI(thisFace,fi_Melem_d)      = jElem 
            faceI(thisFace,fi_etype_u)      = linkI(linkSet(mm),li_link_type)
            faceI(thisFace,fi_etype_d)      = eJunctionChannel
            faceI(thisFace,fi_jump_type)    = nullvalueI
            faceI(thisFace,fi_node_ID)      = thisNode
            faceI(thisFace,fi_link_ID)      = linkSet(mm)
            faceI(thisFace,fi_link_Pos)     = nullvalueI
            
            faceR(thisFace,fr_Zbottom) = nodeR(thisNode,nr_Zbottom) &
                                        + elemMR(jElem,eMr_LengthUp(mm)) &
                                        * linkR(dlink,lr_Slope)
                                        
            faceName(thisFace) = nodeName(thisNode)                            
                
            lastFace = thisFace
            thisFace = thisFace+1
        endif
    enddo
    nullify(dlink)
    
    !%  cycle through the links and subdivide
    do mm = 1,upstream_face_per_elemM
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, 'cycle through linkSet ', mm
        if (linkSet(mm) > 0) then
            if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, 'linkSet(mm) found ',mm
            
            lastElemM = jElem ! back to the element we created (otherise lastElemM has recursion)
            lastFace = elemMI(jElem,eMi_MfaceUp(mm))
            faceI(lastFace,fi_Melem_u) = thisElem2
            
            faceI(lastFace,fi_type) = setFaceType &
                (faceI(lastFace,fi_etype_u), faceI(lastFace,fi_etype_d))
            
            zDownstream = faceR(lastFace,fr_Zbottom)
            
            call subdivide_link_going_upstream &
                (lastElem2, thisElem2, lastFace, thisFace, linkSet(mm),   &
                 elem2I, faceI, linkI, elem2R, faceR, linkR, nodeR, zDownstream) 
                                  
            ! choose the next node to analyze
            thisNode = linkI(linkSet(mm),li_Mnode_u)
            
            ! note that the prior link has been assigned elements
            linkI(linkSet(mm),li_assigned) = setAssigned &
                (linkSet(mm), li_assigned, lUnassigned, lAssigned, linkI)
                
            linkSet(mm) = 0

            if (thisNode > 0) then
                if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, 'Recursion: handle_thisnode ', thisNode

                call handle_thisnode &
                    (lastElem2, thisElem2, lastElemM, thisElemM,          &
                     lastFace,  thisFace,  thisNode,  thisLink,           &
                     elem2I, elemMI, faceI, linkI, nodeI, elem2R, elemMR, faceR, linkR, nodeR, &
                     elem2Name, elemMName, faceName, linkName, nodeName)     
            endif

        endif
    enddo
 else
    print *, 'error: Unknown value for ni_node_type of ',ni_node_type,' in ',subroutine_name
    STOP
 endif

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name 
 end subroutine handle_thisnode   

!==========================================================================
!==========================================================================
!
 subroutine subdivide_link_going_upstream &
    (lastElem2, thisElem2, lastFace, thisFace, thisLink,   &
     elem2I, faceI, linkI, elem2R, faceR, linkR, nodeR, zDownstream)  
!
! Subdivide thisLink from the linkI array into a number of smaller 
! elements and store in the elemI array
!
! this assumes that the face data for thisFace has already been stored
! 
 character(64) :: subroutine_name = 'subdivide_link_going_upstream' 
 integer, intent(in out) :: lastElem2, thisElem2, lastFace, thisFace
 integer, intent(in)     :: thisLink

 integer, intent(in out) :: elem2I(:,:), faceI(:,:)

 integer, intent(in out)     :: linkI(:,:)
 
 real,    intent(in out) :: elem2R(:,:), faceR(:,:)
 real,    intent(in)     :: linkR(:,:), nodeR(:,:)
 
 real,    intent(in)     :: zDownstream
 
 real :: zcenter, zface

 integer :: mm, pp
!-------------------------------------------------------------------------- 

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

!%  store the ID of the first (downstream) element in this link
 linkI(thisLink,li_Melem_d) = thisElem2 ! the new element is the 1st (downstream)
 linkI(thisLink,li_Mface_d) = lastFace ! the old face is the 1st
 
!%  reference elevations at cell center and cell face
 zcenter = zDownstream - 0.5 * linkR(thislink,lr_ElementLength) * linkR(thislink,lr_Slope)
 zface   = zDownstream
 
 do mm = 1,linkI(thisLink,li_N_element)
    !%  store the elem info
    elem2I(thisElem2,e2i_idx)              = thisElem2
    elem2I(thisElem2,e2i_elem_type)        = linkI(thisLink,li_link_type)
    elem2I(thisElem2,e2i_geometry)         = linkI(thislink,li_geometry)
    elem2I(thisElem2,e2i_roughness_type)   = linkI(thisLink,li_roughness_type)
    elem2I(thisElem2,e2i_link_ID)          = thisLink
    elem2I(thisElem2,e2i_link_Pos)         = mm
    elem2I(thisElem2,e2i_Mface_d)          = lastFace
    elem2I(thisElem2,e2i_Mface_u)          = thisFace
    
    elem2R(thisElem2,e2r_Length)           = linkR(thislink,lr_ElementLength)
    
    zcenter = zcenter + linkR(thislink,lr_Slope) * linkR(thislink,lr_ElementLength)
    zface   = zface   + linkR(thislink,lr_Slope) * linkR(thislink,lr_ElementLength)
    
    elem2R(thisElem2,e2r_Zbottom)          = zcenter
        
    select case (linkI(thisLink,li_geometry))
        case (lRectangular)
            elem2R(thisElem2,e2r_BreadthScale) = linkR(thisLink,lr_BreadthScale)
            elem2R(thisElem2,e2r_Topwidth)     = linkR(thisLink,lr_BreadthScale)
            !faceR(thisFace,fr_Topwidth)    = linkR(thisLink,lr_Breadth)
        case (lParabolic)
            elem2R(thisElem2,e2r_BreadthScale) = zeroR
            elem2R(thisElem2,e2r_Topwidth) = twoR &
                * sqrt(linkR(thisLink,lr_InitialDepth)/linkR(thisLink,lr_ParabolaValue))
        case (lTrapezoidal)
            elem2R(thisElem2,e2r_BreadthScale) = linkR(thisLink,lr_BreadthScale)
            elem2R(thisElem2,e2r_Topwidth) = linkR(thisLink,lr_BreadthScale)   &
                    + linkR(thisLink,lr_InitialDepth)                          &
                    * (linkR(thisLink,lr_LeftSlope) + linkR(thisLink,lr_RightSlope))
        case (lTriangle)
            elem2R(thisElem2,e2r_BreadthScale) = zeroR
            elem2R(thisElem2,e2r_Topwidth)     = linkR(thisLink,lr_InitialDepth) &
                    * (linkR(thisLink,lr_LeftSlope) + linkR(thisLink,lr_RightSlope))
        case (lWidthDepth)
            elem2R(thisElem2,e2r_Topwidth)     = linkR(thisLink,lr_Topwidth)
            elem2R(thisElem2,e2r_BreadthScale) = linkR(thisLink,lr_BreadthScale)
            faceR(thisFace,fr_Topwidth)    = linkR(thisLink,lr_Topwidth)
        case default
            print *, 'error: case statement is incomplete in ',subroutine_name
            stop
    end select
        
    !%  store the face info on upstream face
    faceI(thisFace,fi_idx)          = thisFace
    faceI(thisFace,fi_Melem_u)      = thisElem2+1           
    faceI(thisFace,fi_Melem_d)      = thisElem2
    faceI(thisFace,fi_etype_u)      = linkI(thisLink,li_link_type)
    faceI(thisFace,fi_etype_d)      = linkI(thisLink,li_link_type)
    faceI(thisFace,fi_jump_type)    = nullvalueI
    faceI(thisFace,fi_node_ID)      = nullvalueI
    faceI(thisFace,fi_link_ID)      = thisLink
    faceI(thisFace,fi_link_Pos)     = mm
    
    faceR(thisFace,fr_Zbottom)      = zface

    faceI(thisFace,fi_type) = setFaceType &
        (faceI(thisFace,fi_etype_u), faceI(thisFace,fi_etype_d))
        
    lastElem2 = thisElem2
    thisElem2 = thisElem2+1 ! ready to be assigned in next loop
    lastFace = thisFace
    thisFace = thisFace+1 ! ready to be assigned in next loop
 enddo
 
!%  store the ID of the last (upstream) element in this link
 linkI(thisLink,li_Melem_u) = lastElem2 ! the last link stored
 linkI(thisLink,li_Mface_u) = lastFace ! the last face stored

!%  reset for the last face element mapping info - upstream is not yet known.
 faceI(lastFace,fi_Melem_u)  = nullvalueI
 faceI(lastFace,fi_etype_u)  = nullvalueI
 faceI(lastFace,fi_node_ID)  = nullvalueI
 faceI(lastFace,fi_link_ID)  = nullvalueI
 faceI(lastFace,fi_link_Pos) = nullvalueI

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name 
 end subroutine subdivide_link_going_upstream
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
    
 if (arrayI(thisX,assigned_idx) == unassignedValue) then
    f_result = assignedValue
 else
    print *, thisX, assigned_idx, unassignedValue, assignedValue
    print *, arrayI(thisX,assigned_idx)
    print *, 'error: expected an unassigned node or link in ',subroutine_name
    STOP
 endif
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end function setAssigned

!==========================================================================
!==========================================================================
!
 function setFaceType &
    (up_elem_type, dn_elem_type) &
    result(f_result)
!
! defines whether a face is a BC, channel, pipe or multiple connector
!    
 character(64) :: subroutine_name = 'setFaceType'
 
 integer :: f_result
 
 integer, intent(in)    ::  up_elem_type, dn_elem_type
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name
 
 if ( up_elem_type /= dn_elem_type ) then
    if (up_elem_type == eBCup) then
        f_result = fBCup
    elseif (dn_elem_type == eBCdn) then
        f_result = fBCdn
    else
        f_result = fMultiple
    endif
 else
    if (up_elem_type == fChannel) then
        f_result = fChannel
    elseif (up_elem_type == fPipe) then
        f_result = fPipe 
    else
        print *, 'upstream element: ',up_elem_type
        print *, 'dnstream element: ',dn_elem_type
        print *, 'Unexpected if resolution '
        STOP
    endif
 endif 
    
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end function setFaceType
!
!========================================================================== 
!==========================================================================
!
 subroutine link_shortening &
    (linkLength, nodeR, nodeI, thisNode, niNlinkX, niMlinkX, &
     nrElementLengthX, element_nominal_length)
!
! decrements the overall length of a link and adds a portion to
! a junction.
!    
! niNlinkX = ni_N_link_u or ni_N_link_d  (index to number of up/down links)
! niMlinkX = ni_MlinkUp or ni_MlinkDn  (index to map to up/down links)
! nrElementLengthX = nr_ElementLengthUp or nr_ElementLengthDn

 character(64) :: subroutine_name = 'link_shortening'
 
 real,              intent(in out)  ::  linkLength(:)
 real,    target,   intent(in out)  ::  nodeR(:,:)
 
 integer, target,   intent(in)      :: nodeI(:,:)

 integer,           intent(in)      :: thisNode, niNlinkX, nrElementLengthX(:), niMlinkX(:)
 
 real,              intent(in)      :: element_nominal_length(:)
 
 real                   :: delta ! local variable
 
 real,      pointer     :: elemLength
 
 integer,   pointer     :: tlink

 integer    :: mm
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
    

!%  iterate through the up/down (X) links at this node
 do mm=1,nodeI(thisNode,niNlinkX)
 
    !%  use pointers to simplify code
    elemLength => nodeR(thisNode,nrElementLengthX(mm)) ! elem length location
    tlink      => nodeI(thisNode,niMlinkX(mm)) ! link index up or down
      
    if (tlink > 0) then
        if ( linkLength(tlink) .LE. 1.5*element_nominal_length(tlink) ) then
            !%  where there is only one element in a link
            delta = 0.25*linkLength(tlink)
        else
            !%  where there are multiple elements in a link
            delta = 0.33*element_nominal_length(tlink)
        endif
        !%  decrement link
        linkLength(tlink) = linkLength(tlink) - delta
        !%  store junction length in nodeR (target)
        elemLength = delta
    endif
 enddo
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine link_shortening 
!
!========================================================================== 
! END OF MODULE network_define
!========================================================================== 
 end module network_define

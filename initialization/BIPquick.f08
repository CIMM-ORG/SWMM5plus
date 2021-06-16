module BIPquick

  !%***********************************************************************************
  !%
  !% Description:
  !%   This module holds the public BIPquick_partitioning subroutine that is called by the 
  !%   partitioning.f08 utility module.  The private subroutines are used by the BIPquick
  !%   main to do the fairly complicated stuff it does.
  !%
  !%***********************************************************************************
  
  use define_indexes
  use define_globals
  use define_settings
    
   implicit none
  
   private
  
   public :: init_partitioning_BIPquick
  
   real(8), parameter :: precision_matching_tolerance = 1.0D-5 ! a tolerance parameter for whether or not two real(8) numbers are equal
  
   !% The columns for B_nodeR
   integer, parameter :: directweight = oneI
   integer, parameter :: totalweight = twoI
  
   !% The columns for B_nodeI
   integer, parameter :: upstream1 = oneI
   integer, parameter :: upstream2 = twoI
   integer, parameter :: upstream3 = threeI
  
   !HACK - need to figure out how the number of processors is determined
   integer, parameter :: processors = 3

   !HACK - I'm not sure this is still needed
   integer, parameter :: lr_target_default = 1.0
  
   contains
  
  !--------------------------------------------------------------------------
  !% What subroutines do I need to recreate?
  !%    BIPquick main - public subroutine that does the following
  !%      Allocate B_nodeI (to contain the upstream nodes), B_nodeR (for nr_directweight_u, nr_totalweight_u)
  !%      Populate B_nodeI (bip_network_processing)
  !%      phantom_naming_convention()
  !%      calc_directweight()
  !%      calculate_partition_threshold()
  !%      Initialize flag arrays (only add the ones that I'm sure I need)
  !%      do loop for processors
  !%          total_weight_assigner() - includes calc_upstream_weight()
  !%          calc_effective_root() - function that returns effective_root
  !%          trav_subnetwork()    - assigns nodes to images, sets the visited_flags to .true.
  !%          calc_spanning_link()        - determines if any of the links span the partition_threshold
  !%          calc_phantom_node_loc()   - determines where in the link to put the phantom node
  !%          do while loop for case 3
  !%          create_phantom_node()   - create the phantom node in the nodeI, nodeR, linkI, linkR, B_nodeI, B_nodeR
  !%      assign_links()      - use the existing funny logic to assign links to the images
  !%      check_is_boundary() - looks at the adjacent links for every node and increments the ni_P_is_boundary column
  
  subroutine init_partitioning_BIPquick()
      ! ----------------------------------------------------------------------------------------------------------------
      !
      ! Description:
      !   This subroutine serves as the BIPquick main.  It is the only public subroutine from the BIPquick module and its
      !   output is populated columns ni_P_image, ni_P_is_boundary in nodeI and li_P_image in linkI.
      !
      ! -----------------------------------------------------------------------------------------------------------------
    character(64) :: subroutine_name = 'BIPquick_subroutine'
  
    integer       :: mp, ii, jj
    real(8)       :: partition_threshold, max_weight
    logical       :: ideal_exists = .false.
    integer       :: spanning_link = nullValueI
    integer       :: effective_root  
    integer       :: phantom_node_idx, phantom_link_idx
  
    real(8) :: start, intermediate, finish
    call cpu_time(start)
    ! -----------------------------------------------------------------------------------------------------------------
    if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
  
    !% Allocate the B_nodeI and B_nodeR temporary arrays needed for BIPquick
    call bip_allocate_arrays()
  
    !% This subroutine populates B_nodeI with the upstream neighbors for a given node
    call bip_network_processing()
  
    !% This subroutine determines what the phantom nodes will be named
    call phantom_naming_convention(phantom_node_idx, phantom_link_idx)
  
    !% This subroutine populates the directweight column of B_nodeR
    call calc_directweight()

    do ii = 1, size(B_nodeR, 1)
      print*, nodeI(ii, ni_idx), B_nodeR(ii, directweight)
    end do

    stop
  
    do mp = 1, processors
      print*, "Partition", mp
  
      B_nodeR(:, totalweight) = 0.0
  
      !% This subroutine populates the totalweight column of B_nodeR and calculates max_weight
      call calc_totalweight()
  
      partition_threshold = max_weight/real(processors - mp + 1, 8)
  
      !% This subroutine determines if there is an ideal partition possible and what the effective root is
      effective_root = calc_effective_root()
  
      if (ideal_exists) then
  
        !% This subroutine traverses the subnetwork upstream of the effective root
        !% and populates the partition column of nodeI
        call trav_subnetwork()
  
      else
  
        !% This subroutine checks if the partition threshold is spanned by any link
        call calc_spanning_link()
  
        do while (spanning_link == nullValueI)
  
            !% This subroutine houses the litany of steps that are required for a Case 3 partition
            call trav_casethree()
            exit !HACK until we can fill in the guts of this part
        end do
  
        !% This subroutine creates a phantom node/link and adds it to nodeI/linkI
        call phantom_node_generator()
  
        !% This subroutine does the same thing as the previous call to trav_subnetwork()
        call trav_subnetwork()
  
      endif 
  
    end do
  
    !% This subroutine assigns network links to images on the basis of their endpoint nodes
    call trav_assign_link()
  
    !% This subroutine deallocates all of the array variables used in BIPquick
    call bip_deallocate_arrays()
    if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
  end subroutine init_partitioning_BIPquick
  !
  !==========================================================================   
  !========================================================================== 
  !
  subroutine bip_allocate_arrays()
      ! ----------------------------------------------------------------------------------------------------------------
      !
      ! Description:
      !   This subroutine allocates the temporary arrays that are needed for the BIPquick algorithm
      !   These temporary arrays are initialized in Globals, so they're not needed as arguments to BIPquick subroutines
      !
      ! -----------------------------------------------------------------------------------------------------------------
      character(64) :: subroutine_name = 'bip_allocate_arrays'
      ! -----------------------------------------------------------------------------------------------------------------
      if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
  
      allocate(B_nodeI(size(nodeI,1), max_us_branch_per_node))
      B_nodeI(:,:) = nullValueI
  
      allocate(B_nodeR(size(nodeR,1), twoI))
      B_nodeR(:,:) = zeroR
  
      allocate(visited_flag_weight(size(nodeI, oneI)))
      visited_flag_weight(:) = .false.
  
      allocate(visit_network_mask(size(nodeI, oneI)))
      visit_network_mask(:) = .false.
  
      allocate(partition_boolean(size(linkI, oneI)))
      partition_boolean(:) = .false.
  
      allocate(weight_range(size(linkI, oneI), twoI))
      weight_range(:,:) = zeroR
  
      if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
  end subroutine bip_allocate_arrays
  !
  !==========================================================================   
  !========================================================================== 
  !
  subroutine bip_deallocate_arrays()
      ! ----------------------------------------------------------------------------------------------------------------
      !
      ! Description:
      !   This subroutine deallocates the temporary arrays that are needed for the BIPquick algorithm
      !
      ! -----------------------------------------------------------------------------------------------------------------
      character(64) :: subroutine_name = 'bip_deallocate_arrays'
      ! -----------------------------------------------------------------------------------------------------------------
      if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
  
      deallocate(B_nodeI)
      deallocate(B_nodeR)
      deallocate(visited_flag_weight)
      deallocate(visit_network_mask)
      deallocate(partition_boolean)
      deallocate(weight_range)
  
      if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
  end subroutine bip_deallocate_arrays
  !
  !==========================================================================   
  !========================================================================== 
  !
  subroutine bip_network_processing()
    ! ----------------------------------------------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine is a preprocessing step.
    !   The network is traversed once and the B_nodeI array is populated with the adjacent upstream nodes.  This
    !   step eliminates the need to continuously jump back and forth between the nodeI/linkI arrays during subsequent
    !   network traversals.
    !
    ! ----------------------------------------------------------------------------------------------------------------
    character(64) :: subroutine_name = 'bip_network_processing'
  
    integer :: upstream_link, upstream_node, uplink_counter
    integer ii, jj, uplinks
    ! ----------------------------------------------------------------------------------------------------------------
      if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name

      !% Iterate through the nodes array
      do ii= 1,size(nodeI,1)
        ! if ( mod(ii, 1) == 0 ) then
        !     print*, "processing upstream nodes of ni_idx:", ii
        ! endif
        
        !% The number of links upstream of a node
        uplink_counter = nodeI(ii, ni_N_link_u)
        print*, nodeI(ii, ni_idx), "has", nodeI(ii, ni_N_link_u), "upstream links"

        !% Iterate through the links upstream of a node
        do uplinks= 1, uplink_counter

            !% If the link entry is not nullValueI (i.e. if it exists)
            if ( nodeI(ii, ni_idx_base1 + uplinks) /= nullValueI ) then

                !% Go to the upstream link to find the next upstream node
                upstream_link = nodeI(ii, ni_idx_base1 + uplinks)
                upstream_node = linkI(upstream_link, li_Mnode_u)

                !% Add the adjacent upstream node to B_nodeI
                B_nodeI(ii, uplinks) = upstream_node

            endif
        enddo
     enddo
    
      if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
  end subroutine bip_network_processing
  !
  !============================================================================ 
  !============================================================================ 
  ! 
  subroutine BIPquick_Optimal_Hardcode()
      integer :: ii
  
      nodeI(:, ni_P_image) = (/1, 1, 1, 2, 3, 3, 1, 1, 2, 2, 2, 2, 3, 3, 3, 1, 3/)
      linkI(:, li_P_image) = (/1, 1, 1, 2, 3, 3, 3, 2, 2, 2, 2, 2, 3, 3, 3, 1, 1, 1/)
  
  end subroutine BIPquick_Optimal_Hardcode    
  !
  !============================================================================ 
  !============================================================================ 
  ! 
  subroutine BIPquick_YJunction_Hardcode()
    integer :: ii
  
    nodeI(:, ni_P_image) = (/1, 2, 1, 3, -998877, -998877/)
    nodeI(:, ni_P_is_boundary) = (/0, 0, 1, 0, -998877, -998877/)
    linkI(:, li_P_image) = (/1, 2, 3, -998877, -998877/)
  
  end subroutine BIPquick_YJunction_Hardcode
  !
  !============================================================================ 
  !============================================================================ 
  !
  function calc_link_weights(link_index) result(weight)
    ! ----------------------------------------------------------------------------
    !  
    ! Description:
    !   the weight attributed to each link (that will ultimately be assigned to the 
    !   downstream node) are normalized by lr_Target.  This gives an estimate of 
    !   computational complexity. In the future lr_Target can be customized for each 
    !   link.
    !  
    ! ----------------------------------------------------------------------------
    character(64)   :: function_name = 'calc_link_weights'
  
    integer, intent(in) :: link_index
    real(8)             :: weight
    ! --------------------------------------------------------------------------
    if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name

    !% The link weight is equal to the link length divided by the element length
    weight = linkR(link_index, lr_Length) / linkR(link_index, lr_ElementLength)
  
    if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
  end function calc_link_weights
  !
  !============================================================================
  !============================================================================
  !
  subroutine calc_directweight()
   !----------------------------------------------------------------------------
   !
   ! Description:
   !   This subroutine looks at the upstream links for each node, calculates their
   !   link weights using weighting_function(), and sums those weights to yield the
   !   B_nodeR(node, directweight)
   !
   !----------------------------------------------------------------------------
    character(64) :: subroutine_name = 'calc_directweight'
    real(8)       :: calc_link_weights_output ! HACK
  
    real(8)  :: lr_target
    integer :: rootnode_index, links_row, upstream_links
    integer :: ii, jj
  
   !--------------------------------------------------------------------------
    if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
    
    !% Calculates directweight for each node
    do ii= 1,size(nodeI,1)

      !% Need a loop bc multiple links might have a given node as its downstream endpoint
      do jj=1,size(linkI(:, li_Mnode_d))

        !% If the link has the current node as a downstream endpoint
        if (linkI(jj, li_Mnode_d) == nodeI(ii, ni_idx)) then
            
          !% The directweight for that node is the running total of link weights
          B_nodeR(ii, directweight) = B_nodeR(ii, directweight) &
            + calc_link_weights(linkI(jj, li_idx))
        endif
      enddo
    enddo
   
    if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
  end subroutine calc_directweight
   !
   !============================================================================
   !============================================================================
   !
  recursive subroutine calc_upstream_weight()
     !-----------------------------------------------------------------------------
     !
     ! Description: Recursive subroutine that visits each node upstream of some root
     !  and adds the directweight to the root's totalweight.  This recursive subroutine
     !  is called for each node remaining in the network.
     !
     !-----------------------------------------------------------------------------
  
     character(64) :: subroutine_name = 'upstream_weight_calculation'

    !  real(8), intent(in out) :: B_nodeR(:,:), nodeR(:,:), linkR(:,:)
    !  integer, intent(in out) :: B_nodeI(:,:), nodeI(:,:), linkI(:,:)
     integer :: upstream_node_list(3)
    !  integer, intent(in) :: weight_index
     integer :: root, node_row_contents, link_idx, new_root
     integer :: link_row_contents, node_upstream
    !  logical, intent(in out) :: visited_flag_weight(:)
    !  logical, intent(in) :: visit_network_mask(:)
     integer :: ii, jj, kk
     !--------------------------------------------------------------------------
     if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
  
  
     if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
  end subroutine calc_upstream_weight
   !
   !============================================================================
   !============================================================================
   !
  subroutine calc_totalweight()
    !-----------------------------------------------------------------------------
    !
    ! Description: This subroutine drives the calc_upstream_weight() recursive
    !  subroutine.  If a node remains in the network (i.e. hasn't been assigned to a
    !  partition yet), then it is passed as a root to the calc_upstream_weight().
    !
    !-----------------------------------------------------------------------------
  
    character(64) :: subroutine_name = 'calc_totalweight'
  
    ! real(8), intent(in out) :: B_nodeR(:,:), nodeR(:,:), linkR(:,:)
    ! integer, intent(in out) :: B_nodeI(:,:), nodeI(:,:), linkI(:,:), B_node_Partition(:,:), B_link_Partition(:,:)
    ! integer, intent(in out) :: weight_index
    ! logical, intent(in out) :: visited_flag_weight(:)
    ! logical, intent(in out) :: visit_network_mask(:)
    ! real(8), intent(in out) :: max_weight
    integer :: ii, root
    real(8) :: nullValue = nullValueI
    !--------------------------------------------------------------------------
    if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
    if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name

  end subroutine calc_totalweight
  !
  !============================================================================
  !============================================================================
  !
  recursive subroutine trav_subnetwork() 
    !-----------------------------------------------------------------------------
    !
    ! Description: This recursive subroutine visits every node upstream of a root node
    !  (where the root node is the "effective_root") and assigns that node to the current
    !  image.  It also updates the visit_network_mask(:) array to "remove" that node from
    !  the network.
    !
    !-----------------------------------------------------------------------------
  
    character(64) :: subroutine_name = 'trav_subnetwork'
  
    ! logical, intent(in out) :: visit_network_mask(:)
    integer :: upstream_node_list(3)
    integer :: node_row_contents, link_row_contents, new_root, node_upstream
    ! integer, intent(in) :: root, proc
    integer :: ii, jj, kk
    ! integer, intent (in out) :: print_counter
    !--------------------------------------------------------------------------
    if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
  
  
    if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
  end subroutine trav_subnetwork
   !
   !============================================================================
   !============================================================================
   !
  subroutine trav_assign_link()
   !-----------------------------------------------------------------------------
   !
   ! Description: This subroutine is used to assign links to images.  BIPquick first
   !  assigns nodes to images, and then uses that info to assign links to images.
   !
   !-----------------------------------------------------------------------------
  
    character(64) :: subroutine_name = 'trav_assign_link'
  
    !--------------------------------------------------------------------------
    if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
  
  
    if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
  end subroutine trav_assign_link
  !
  !============================================================================
  !============================================================================
  !
  function calc_effective_root() result (effective_root)
   !-----------------------------------------------------------------------------
   !
   ! Description: This function is used to search for the effective root.  The effective
   !  root can fit one of two descriptions: it can have a totalweight that is exactly
   !  equal to the partition_threshold (Case 1), or it can have a totalweight that is the 
   !  floor of the range (partition_threshold, max_weight] (nearest_overestimate_node, Case 3).
   !
   !-----------------------------------------------------------------------------

    character(64) :: subroutine_name = 'calc_effective_root'
    integer :: effective_root
  
    ! real(8), intent(in) :: max_weight, partition_threshold
    ! logical, intent(in out) :: ideal_exists
    real(8) :: nearest_overestimate
    ! real(8), intent(in) :: B_nodeR(:,:)
    ! integer, intent(in) :: nodeI(:,:)
    integer :: ii
    integer :: nullValue = nullValueI
    !--------------------------------------------------------------------------
    if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
  
    effective_root = oneI ! HACK
  
    if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
  end function calc_effective_root
  !
  !============================================================================
  !============================================================================
  !
  subroutine calc_spanning_link()
   !-----------------------------------------------------------------------------
   !
   ! Description: This subroutine is used to search for a link that spans the 
   !  partition_threshold (Case 2).  A link is considered to "span" if the 
   !  partition_threshold is in the range (totalweight_upstream_node, 
   !  totalweight_upstream_node + weight_of_link).
   !
   !-----------------------------------------------------------------------------
  
    character(64) :: subroutine_name = 'calc_spanning_link'
  
    ! integer, intent(in out) :: spanning_link
    ! integer, intent(out) :: spanning_node_upstream
    ! real(8), allocatable, intent(in out) :: weight_range(:,:)
    ! real(8), intent(in) :: lr_target, partition_threshold
    ! logical, intent(in out) :: partition_boolean(:)
    integer :: upstream_node
    integer :: ii, jj
    !--------------------------------------------------------------------------
    if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
  
    if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
  end subroutine calc_spanning_link
  !
  !==========================================================================
  !==========================================================================
  !
  function calc_phantom_node_loc() result(length_from_start)
    !-----------------------------------------------------------------------------
   !
   ! Description:  This function is used to calculate how far along the spanning_link 
   !  (from the upstream node) the phantom_node should be placed.
  
      ! HACK - might want to consider snapping this to the nearest element subdivision
   !
   !-----------------------------------------------------------------------------  
    character(64)   :: function_name = 'calc_phantom_node_loc'
  
    ! real(8), intent(in) :: weight_range(:,:)
    ! real(8), intent(in) :: partition_threshold, lr_target
    ! integer, intent(in) :: spanning_link
    real(8) :: length_from_start, total_length, start, weight_ratio
    integer :: ii
    !--------------------------------------------------------------------------
    if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name
   
    if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
  end function calc_phantom_node_loc
  !
  !==========================================================================
  !==========================================================================
  !
  subroutine phantom_naming_convention(phantom_node_idx, phantom_link_idx)
   !-----------------------------------------------------------------------------
   !
   ! Description: This subroutine is used to establish what the first phantom node/link
   !  will be named.  The phantom_node/link is simply named for the index of the nodeI/linkI
   !  array in which it is populated.
   !
   !-----------------------------------------------------------------------------
  
    character(64)   :: function_name = 'phantom_naming_convention'
  
    integer, intent(in out) :: phantom_node_idx, phantom_link_idx
    !--------------------------------------------------------------------------
    if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name
     
    !% The phantom_node_idx is equal to 1 more than the largest ni_idx (that is not nullValueI)
    phantom_node_idx = maxval(nodeI(:, ni_idx), MASK= ( nodeI(:, ni_idx) /= nullValueI )) + 1

    !% The phantom_link_idx is equal to 1 more than the largest li_idx (that is not nullValueI)
    phantom_link_idx = maxval(linkI(:, li_idx), MASK= ( linkI(:, ni_idx) /= nullValueI )) + 1
     
    if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
  end subroutine phantom_naming_convention
  !
  !==========================================================================
  !==========================================================================
  !
  subroutine phantom_node_generator()
   !-----------------------------------------------------------------------------
   !
   ! Description: This subroutine populates the nodeI/linkI/linkR arrays with the phantom
   !  node/link that has been generated by a Case 2 system.  The nodeI entries
   !  are set from the properties of phantom_nodes, the linkI/R phantom entries are copied
   !  from the linkI/R real entries then updated.
   !
   !-----------------------------------------------------------------------------
    character(64) :: subroutine_name = 'phantom_node_generator'
  
    ! real(8), intent(in) :: start_point
    ! integer, intent(in) :: spanning_link, phantom_node_idx, phantom_link_idx, spanning_node_upstream
    ! integer, intent(in) :: mp
    ! real(8), intent(in) :: partition_threshold
    ! real(8), intent(in)  :: lr_target
    ! integer, intent(in out) :: phantom_array_location
    integer :: phantom_counter = 0
    integer :: phantom_name, phantom_array_location_link
    integer :: downstream_node, upstream_spanning_node
    integer :: ii, jj, kk
    !--------------------------------------------------------------------------
    if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
  
  
    if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
  end subroutine phantom_node_generator
  !
  !==========================================================================
  !==========================================================================
  !
  subroutine trav_casethree()
   !-----------------------------------------------------------------------------
   !
   ! Description: This subroutine drives the steps required to reduce a Case 3 network
   !  into a Case 2 network (that has a spanning link).  This reduction involves calling
   !  trav_subnetwork() on the effective_root, updating the partition_threshold, and 
   !  searching the remaining network for Case 1 or 2.  This iterative reduction occurs
   !  in a do while loop until a Case 1 or 2 network is found.
   !
   !-----------------------------------------------------------------------------
  
    character(64) :: subroutine_name = 'trav_casethree'
      
    !--------------------------------------------------------------------------
    if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
   
   
    if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
   end subroutine trav_casethree
  !
  !============================================================================
  !============================================================================
  !
   end module BIPquick
  !==========================================================================
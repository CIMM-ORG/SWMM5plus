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

 contains

!--------------------------------------------------------------------------
!% What subroutines do I need to recreate?
!%    BIPquick main - public subroutine that does the following
!%      Allocate B_nodeI (to contain the upstream nodes), B_nodeR (for nr_directweight_u, nr_totalweight_u)
!%      Populate B_nodeI (network_node_preprocessing)
!%      phantom_naming_convention()
!%      local_node_weighting()
!%      calculate_partition_threshold()
!%      Initialize flag arrays (only add the ones that I'm sure I need)
!%      do loop for processors
!%          total_weight_assigner() - includes upstream_weight_calculation()
!%          ideal_partition_check() - function that returns effective_root
!%          subnetwork_carving()    - assigns nodes to images, sets the visited_flags to .true.
!%          spanning_check()        - determines if any of the links span the partition_threshold
!%          linear_interpolator()   - determines where in the link to put the phantom node
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

  integer       :: mp
  real(8)       :: partition_threshold, max_weight
  logical       :: ideal_exists = .false.
  integer       :: spanning_link = nullValueI
  integer       :: effective_root

  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name

  !% Allocate the B_nodeI and B_nodeR temporary arrays needed for BIPquick
  call allocate_BIPquick_arrays()

  !% This subroutine populates B_nodeI with the upstream neighbors for a given node
  call network_node_preprocessing()

  !% This subroutine determines what the phantom nodes will be named
  call phantom_naming_convention()

  !% This subroutine populates the directweight column of B_nodeR
  call local_node_weighting()

  do mp = 1, processors
    print*, "Partition", mp

    B_nodeR(:, totalweight) = 0.0

    !% This subroutine populates the totalweight column of B_nodeR and calculates max_weight
    call nr_totalweight_assigner()

    partition_threshold = max_weight/real(processors - mp + 1, 8)

    !% This subroutine determines if there is an ideal partition possible and what the effective root is
    effective_root = ideal_partition_check()

    if (ideal_exists) then

      !% This subroutine traverses the subnetwork upstream of the effective root
      !% and populates the partition column of nodeI
      call subnetwork_carving()

    else

      !% This subroutine checks if the partition threshold is spanned by any link
      call spanning_check()

      do while (spanning_link == nullValueI)

          !% This subroutine houses the litany of steps that are required for a Case 3 partition
          call non_ideal_partitioning()
          exit !HACK until we can fill in the guts of this part
      end do

      !% This subroutine creates a phantom node/link and adds it to nodeI/linkI
      call phantom_node_generator()

      !% This subroutine does the same thing as the previous call to subnetwork_carving()
      call subnetwork_carving()

    endif 

  end do

  !% This subroutine assigns network links to images on the basis of their endpoint nodes
  call assign_links()

  !% This subroutine deallocates all of the array variables used in BIPquick
  call deallocate_BIPquick_arrays()
  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end subroutine init_partitioning_BIPquick
!    
!==========================================================================   
!========================================================================== 
!
subroutine allocate_BIPquick_arrays()
    ! ----------------------------------------------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine allocates the temporary arrays that are needed for the BIPquick algorithm
    !   These temporary arrays are initialized in Globals, so they're not needed as arguments to BIPquick subroutines
    !
    ! -----------------------------------------------------------------------------------------------------------------
    character(64) :: subroutine_name = 'allocate_BIPquick_arrays'
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
end subroutine allocate_BIPquick_arrays
!    
!==========================================================================   
!========================================================================== 
!
subroutine deallocate_BIPquick_arrays()
    ! ----------------------------------------------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine deallocates the temporary arrays that are needed for the BIPquick algorithm
    !
    ! -----------------------------------------------------------------------------------------------------------------
    character(64) :: subroutine_name = 'deallocate_BIPquick_arrays'
    ! -----------------------------------------------------------------------------------------------------------------
    if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name

    deallocate(B_nodeI)
    deallocate(B_nodeR)
    deallocate(visited_flag_weight)
    deallocate(visit_network_mask)
    deallocate(partition_boolean)
    deallocate(weight_range)

    if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end subroutine deallocate_BIPquick_arrays
!    
!==========================================================================   
!========================================================================== 
!
  subroutine network_node_preprocessing()
  ! ----------------------------------------------------------------------------------------------------------------
  !
  ! Description:
  !   This subroutine is a preprocessing step.
  !   The network is traversed once and the B_nodeI array is populated with the adjacent upstream nodes.  This
  !   step eliminates the need to continuously jump back and forth between the nodeI/linkI arrays during subsequent
  !   network traversals.
  !
  ! ----------------------------------------------------------------------------------------------------------------
  character(64) :: subroutine_name = 'network_node_preprocessing'
  ! ----------------------------------------------------------------------------------------------------------------
    if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
  
  
    if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
  end subroutine network_node_preprocessing
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
function weighting_function() result(weight)
    ! ----------------------------------------------------------------------------
    !
    ! Description:
    !   The weight attributed to each link (that will ultimately be assigned to the 
    !   downstream node) are normalized by lr_Target.  This gives an estimate of 
    !   computational complexity. In the future lr_Target can be customized for each 
    !   link.
    !
    ! ----------------------------------------------------------------------------
  character(64)   :: function_name = 'weighting_function'

    ! --------------------------------------------------------------------------
  if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name

  if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
end function weighting_function
!
!============================================================================
!============================================================================
!
subroutine local_node_weighting()
 !----------------------------------------------------------------------------
 !
 ! Description:
 !   this function takes each node index (ni_idx) and finds that ni_idx in the 
 !   downstream node column of the links array (li_Mnode_d).  From this row in 
 !   links the link weight is grabbed and ascribed to the node-in-questions local 
 !   weight (nr_directweight_u).
 !
 !----------------------------------------------------------------------------
  character(64) :: subroutine_name = 'local_node_weighting'
  real(8)       :: weighting_function_output ! HACK
 !--------------------------------------------------------------------------
  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name

  call null_value_convert()

  weighting_function_output = weighting_function()

  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end subroutine local_node_weighting
 !
 !============================================================================
 !============================================================================
 !
recursive subroutine upstream_weight_calculation()
   !-----------------------------------------------------------------------------
   !
   ! Description: Recursive subroutine that visits each node upstream of some root
   !  and adds the directweight to the root's totalweight.  This recursive subroutine
   !  is called for each node remaining in the network.
   !
   !-----------------------------------------------------------------------------

   character(64) :: subroutine_name = 'upstream_weight_calculation'

   !--------------------------------------------------------------------------
   if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name


   if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end subroutine upstream_weight_calculation
 !
 !============================================================================
 !============================================================================
 !
subroutine nr_totalweight_assigner()
 !-----------------------------------------------------------------------------
 !
 ! Description: This subroutine drives the upstream_weight_calculation() recursive
 !  subroutine.  If a node remains in the network (i.e. hasn't been assigned to a
 !  partition yet), then it is passed as a root to the upstream_weight_calculation().
 !
 !-----------------------------------------------------------------------------

  character(64) :: subroutine_name = 'nr_totalweight_assigner'

 !--------------------------------------------------------------------------
  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name


  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end subroutine nr_totalweight_assigner
 !
 !============================================================================
 !============================================================================
 !
recursive subroutine subnetwork_carving() 
 !-----------------------------------------------------------------------------
 !
 ! Description: This recursive subroutine visits every node upstream of a root node
 !  (where the root node is the "effective_root") and assigns that node to the current
 !  image.  It also updates the visit_network_mask(:) array to "remove" that node from
 !  the network.
 !  
 !-----------------------------------------------------------------------------

  character(64) :: subroutine_name = 'subnetwork_carving'

 !--------------------------------------------------------------------------
  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name


  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end subroutine subnetwork_carving
 !
 !============================================================================
 !============================================================================
 !
subroutine subnetworks_links()
 !-----------------------------------------------------------------------------
 !
 ! Description: This subroutine is used to assign links to images.  BIPquick first
 !  assigns nodes to images, and then uses that info to assign links to images.
 !
 !-----------------------------------------------------------------------------

   character(64) :: subroutine_name = 'subnetworks_links'

 !--------------------------------------------------------------------------
  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name


  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end subroutine subnetworks_links
!
!============================================================================
!============================================================================
!
function ideal_partition_check() result (effective_root)
 !-----------------------------------------------------------------------------
 !
 ! Description: This function is used to search for the effective root.  The effective
 !  root can fit one of two descriptions: it can have a totalweight that is exactly
 !  equal to the partition_threshold (Case 1), or it can have a totalweight that is the 
 !  floor of the range (partition_threshold, max_weight] (nearest_overestimate_node, Case 3).
 !
 !-----------------------------------------------------------------------------

  character(64) :: subroutine_name = 'ideal_partition_check'
  integer :: effective_root
 !--------------------------------------------------------------------------
  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name

  effective_root = oneI ! HACK


  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end function ideal_partition_check
!
!============================================================================
!============================================================================
!
subroutine spanning_check()
 !-----------------------------------------------------------------------------
 !
 ! Description: This subroutine is used to search for a link that spans the 
 !  partition_threshold (Case 2).  A link is considered to "span" if the 
 !  partition_threshold is in the range (totalweight_upstream_node, 
 !  totalweight_upstream_node + weight_of_link).
 !
 !-----------------------------------------------------------------------------

  character(64) :: subroutine_name = 'spanning_check'
 !--------------------------------------------------------------------------
 if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name

 if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end subroutine spanning_check
!
!==========================================================================
!==========================================================================
!
function linear_interpolator() result(length_from_start)
 !-----------------------------------------------------------------------------
 !
 ! Description:  This function is used to calculate how far along the spanning_link 
 !  (from the upstream node) the phantom_node should be placed.

    ! HACK - might want to consider snapping this to the nearest element subdivision
 !
 !-----------------------------------------------------------------------------

  character(64)   :: function_name = 'linear_interpolator'

 !--------------------------------------------------------------------------
  if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name
 
  if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
end function linear_interpolator
!
!==========================================================================
!==========================================================================
!
subroutine phantom_naming_convention()
 !-----------------------------------------------------------------------------
 !
 ! Description: This subroutine is used to establish what the first phantom node/link
 !  will be named.  The phantom_node/link is simply named for the index of the nodeI/linkI
 !  array in which it is populated.
 !
 !-----------------------------------------------------------------------------

  character(64)   :: function_name = 'phantom_naming_convention'

 !--------------------------------------------------------------------------
  if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name


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

 !--------------------------------------------------------------------------
  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name



  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end subroutine phantom_node_generator
!
!==========================================================================
!==========================================================================
!
subroutine non_ideal_partitioning()
 !-----------------------------------------------------------------------------
 !
 ! Description: This subroutine drives the steps required to reduce a Case 3 network
 !  into a Case 2 network (that has a spanning link).  This reduction involves calling
 !  trav_subnetwork() on the effective_root, updating the partition_threshold, and 
 !  searching the remaining network for Case 1 or 2.  This iterative reduction occurs
 !  in a do while loop until a Case 1 or 2 network is found.
 !
 !-----------------------------------------------------------------------------

  character(64) :: subroutine_name = 'non_ideal_partitioning'
 
  !--------------------------------------------------------------------------
   if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
 
 
   if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
 end subroutine non_ideal_partitioning
!
!============================================================================
!============================================================================
!
 end module BIPquick
!==========================================================================

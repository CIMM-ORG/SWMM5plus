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

    print*, "Entered init_partitioning_BIPquick"
    call allocate_BIPquick_arrays()

   
    call deallocate_BIPquick_arrays()
    print*, "Leaving init_partitioning_BIPquick"
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
    B_nodeR(:,:) = nullValueR

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
     integer :: upstream_link, upstream_node, uplink_counter, root_node
     integer ii, jj, uplinks
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

     

  end subroutine BIPquick_Optimal_Hardcode    
!
!============================================================================ 
!============================================================================ 
! 
 subroutine BIPquick_YJunction_Hardcode()
     integer :: ii

 end subroutine BIPquick_YJunction_Hardcode
!
!============================================================================ 
!============================================================================ 
!
 subroutine BIPquick_subroutine()
     !-----------------------------------------------------------------------------
     ! Description:
     !   This subroutine is the main BIPquick subroutine
     !   It uses the link-node arrays initialized in $ call initialize_linknode_arrays() $ in the main.f08
     !   It also uses B_nodeI, B_linkI arrays that are initialized in allocate_storage.f08
     !   Two dummy arrays, B_nodeI and B_nodeR are used to contain some of the BIPquick specific parameters
     ! Method:
     !    
     !-----------------------------------------------------------------------------



 end subroutine BIPquick_subroutine
!
!============================================================================ 
!============================================================================ 
!
  !function sorted_metric() result(sorted_connectivity_metric)
      !
      ! This subroutine determines how many points of connectivity exist in the
      ! BIPquick-sorted network
      !character(64) :: function_name = 'sorted_metric'

      !real(8), intent(in) :: nodeMatrix(:, :), linkMatrix(:,:), subnetwork_container_links(:, :, :)
      !integer, intent(in) :: multiprocessors
      !real(8), allocatable :: sorted_nodes(:, :)
      !integer :: upstream_node, downstream_node
      !integer :: ii, jj, mp
      !integer :: node_basis_counter = 0
      !integer :: node_contained_counter = 0
      !integer :: sorted_connectivity_metric
      !integer :: nullValue = -998877
      !--------------------------------------------------------------------------
      !if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name



   !if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
  !end function sorted_metric
!
!============================================================================ 
!============================================================================ 
!
  !function unsorted_metric() &
      !result(unsorted_connectivity_metric)
      !
      ! This subroutine determines how many points of connectivity exist in the
      ! BIPquick-sorted network
      !character(64) :: function_name = 'sorted_metric'

      !--------------------------------------------------------------------------
      !if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name



      !if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
  !end function unsorted_metric
 !
 !============================================================================
 !============================================================================
 !
  !function weighting_function() !result(weight)
      !----------------------------------------------------------------------------
      !
      ! Description:
      !   the weight attributed to each link (that will ultimately be assigned to the 
      !   downstream node) are normalized by lr_Target.  This gives an estimate of 
      !   computational complexity. In the future lr_Target can be customized for each 
      !   link.
      !
      !----------------------------------------------------------------------------
      !character(64)   :: function_name = 'weighting_function'

      !--------------------------------------------------------------------------
      !if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name

      !if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
  !end function weighting_function
!
!============================================================================
!============================================================================
!
subroutine null_value_convert()
 !----------------------------------------------------------------------------
 !
 ! Description:
 !   this function is used to convert the null values that are defaulted in the 
 !   B_nr_directweight_u and B_nr_totalweight_u columns of the nodes array into float 
 !   zeros.
 !
 !----------------------------------------------------------------------------
  character(64) :: subroutine_name = 'null_value_convert'

  !real(8), intent(in out) :: array(:)
  integer :: ii
  !real(8) :: nullValue = nullvalueR
 !--------------------------------------------------------------------------
  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name


  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end subroutine null_value_convert
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

  !real(8)  :: lr_target
  !real(8), intent(in)  :: lr_target_default
  !integer, intent(in) :: nodeI(:,:), linkI(:,:)
  !real(8), intent(in out) :: nodeR(:,:), linkR(:,:), B_nodeR(:,:)
  !integer :: rootnode_index, links_row
  integer ii, jj

 !--------------------------------------------------------------------------
  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name


  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end subroutine local_node_weighting
 !
 !============================================================================
 !============================================================================
 !
recursive subroutine upstream_weight_calculation()
   !-----------------------------------------------------------------------------
   !
   ! Description:
   !
   !
   ! Method:
   !    
   !-----------------------------------------------------------------------------

   character(64) :: subroutine_name = 'upstream_weight_calculation'

   !  real(8), intent(in out) :: B_nodeR(:,:), nodeR(:,:), linkR(:,:)
   !  integer, intent(in out) :: B_nodeI(:,:), nodeI(:,:), linkI(:,:), B_node_Partition(:,:), B_link_Partition(:,:)
   !  integer :: upstream_node_list(3)
   !  integer, intent(in) :: weight_index
   !  integer :: root, node_row_contents, link_idx, new_root
   !  integer :: link_row_contents, node_upstream
   !  logical, intent(in out) :: visited_flag_weight(:)
   !  logical, intent(in) :: visit_network_mask(:)
   !  integer :: ii, jj, kk

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
 ! Description:
 !
 !
 ! Method:
 !    
 !-----------------------------------------------------------------------------

  character(64) :: subroutine_name = 'nr_totalweight_assigner'

  !real(8), intent(in out) :: B_nodeR(:,:), nodeR(:,:), linkR(:,:)
  !nteger, intent(in out) :: B_nodeI(:,:), nodeI(:,:), linkI(:,:), B_node_Partition(:,:), B_link_Partition(:,:)
  !integer, intent(in out) :: weight_index
  !logical, intent(in out) :: visited_flag_weight(:)
  !logical, intent(in out) :: visit_network_mask(:)
  !real(8), intent(in out) :: max_weight
  !integer :: ii, root
  !real(8) :: nullValue = nullValueI

 !--------------------------------------------------------------------------
  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name


  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end subroutine nr_totalweight_assigner
 !
 !============================================================================
 !============================================================================
 !
recursive subroutine subnetwork_carving() 

  character(64) :: subroutine_name = 'subnetwork_carving'

 !  logical, intent(in out) :: visit_network_mask(:)
 !  real(8), intent (in out) :: subnetwork_container_nodes (:,:,:)
 !  real(8), intent(in) :: nodeMatrix(:,:), linkMatrix(:,:)
 !  integer :: upstream_node_list(3)
 !  integer :: node_row_contents, link_row_contents, new_root, node_upstream
 !  integer, intent(in) :: root, proc
 !  integer :: ii, jj, kk
 !  integer, intent (in out) :: print_counter

 !--------------------------------------------------------------------------
  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name


  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end subroutine subnetwork_carving
 !
 !============================================================================
 !============================================================================
 !
subroutine subnetworks_links()

   character(64) :: subroutine_name = 'subnetworks_links'

 !  real(8) :: endpoint1, endpoint2
 !  real(8), intent(in) :: nodes_container(:,:)
 !  real(8), intent(in out) :: subnetwork_container_links(:, :, :)
 !  real(8), intent(in) :: linkMatrix(:,:)
 !  integer, allocatable :: potential_endpoints(:)
 !  integer, intent(in out) :: accounted_for_links(:)
 !  integer, intent(in out) :: link_counter
 !  integer, intent(in) :: proc
 !  integer :: ii, jj, linkCounter,mp
 !  real(8) :: accountingLink
 !  logical :: accountedLink = .false.
 ! !--------------------------------------------------------------------------
  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name



  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end subroutine subnetworks_links
!
!============================================================================
!============================================================================
!

!  function ideal_partition_check &
 !     (ideal_exists, max_weight, partition_threshold, B_nodeR, nodeI) result (effective_root)

 !  character(64) :: subroutine_name = 'ideal_partition_check'
 !  real(8), intent(in) :: max_weight, partition_threshold
 !  logical, intent(in out) :: ideal_exists
 !  integer :: effective_root
 !  real(8) :: nearest_overestimate
 !  real(8), intent(in) :: B_nodeR(:,:)
 !  integer, intent(in) :: nodeI(:,:)
 !  integer :: ii
 !  integer :: nullValue = nullValueI
 ! !--------------------------------------------------------------------------
 !  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name


 !  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
!  end function ideal_partition_check
 !
 !============================================================================
 !============================================================================
 !
subroutine spanning_check()

  character(64) :: subroutine_name = 'spanning_check'

 !  integer, intent(in out) :: spanning_link
 !  integer, intent(out) :: spanning_node_upstream
 !  real(8), allocatable, intent(in out) :: weight_range(:,:)
 !  real(8), intent(in) :: linkMatrix(:,:), nodeMatrix(:,:)
 !  real(8), intent(in) :: lr_target, partition_threshold
 !  logical, intent(in out) :: partition_boolean(:)
 !  integer :: upstream_node
 !  integer :: ii, jj

 ! !--------------------------------------------------------------------------
 if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name

 5568 if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end subroutine spanning_check
!
!==========================================================================
!==========================================================================
!
! !  function linear_interpolator(partition_threshold, spanning_link, linkMatrix, &
 !     weight_range, lr_target) result(length_from_start)

 !  character(64)   :: function_name = 'linear_interpolator'

 !  real(8), intent(in) :: linkMatrix(:,:), weight_range(:,:)
 !  real(8), intent(in) :: partition_threshold, lr_target
 !  integer, intent(in) :: spanning_link
 !  real(8) :: length_from_start, total_length, start, weight_ratio
 !  integer :: ii
 ! !--------------------------------------------------------------------------
 !  if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name
 !
 !  if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
! !  end function linear_interpolator
!
!==========================================================================
!==========================================================================
!
subroutine phantom_naming_convention2()

 character(64)   :: function_name = 'phantom_naming_convention2'

 !  real(8), intent(in) :: linkMatrix(:,:), nodeMatrix(:,:)
 !  integer, intent(out) :: phantom_node_idx, phantom_link_idx
 !  integer :: max_node_idx, max_link_idx
 !  integer, allocatable :: node_indices(:), link_indices(:)
 ! !--------------------------------------------------------------------------
 if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name



 if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
end subroutine phantom_naming_convention2
!
!==========================================================================
!==========================================================================
!
subroutine phantom_naming_convention()
 !-----------------------------------------------------------------------------
 !
 ! Description:
 !
 !
 ! Method:
 !    
 !-----------------------------------------------------------------------------

  character(64)   :: function_name = 'phantom_naming_convention'

 !  integer, intent(in) :: linkI(:,:), nodeI(:,:)
 !  integer, intent(out) :: phantom_node_idx, phantom_link_idx
 !  integer :: max_node_idx, max_link_idx
 !  integer, allocatable :: node_indices(:), link_indices(:)
 ! !--------------------------------------------------------------------------
  if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name


  if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
end subroutine phantom_naming_convention
!
!==========================================================================
!==========================================================================
!
subroutine phantom_node_generator()
 ! !
 ! !
  character(64) :: subroutine_name = 'phantom_node_generator'

 !  real(8), intent(in) :: start_point
 !  integer, intent(in) :: spanning_link, phantom_node_idx, phantom_link_idx, spanning_node_upstream
 !  integer, intent(in) :: mp
 !  real(8), intent(in) :: partition_threshold
 !  real(8), intent(in)  :: lr_target
 !  real(8), intent(in out) :: linkMatrix(:,:), nodeMatrix(:,:)
 !  real(8), intent (in out) :: subnetwork_container_nodes (:,:,:)
 !  real(8), intent (in out) :: subnetwork_container_links (:,:,:)
 !  integer, intent(in) :: n_rows_excluding_header_node
 !  integer, intent(in) :: n_rows_excluding_header_link
 !  integer, intent(in out) :: phantom_array_location
 !  integer :: phantom_counter = 0
 !  integer :: phantom_name, phantom_array_location_link
 !  integer :: downstream_node, upstream_spanning_node
 !  integer :: ii, jj, kk
 ! !--------------------------------------------------------------------------
  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name



  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end subroutine phantom_node_generator
!
!==========================================================================
!==========================================================================
!
!  function number_of_lines_in_file(iunit) result(n_lines)

 !  integer,intent(in)  :: iunit   ! the file unit number
 !  integer             :: n_lines ! the number of lines in the file

 !  character(len=1)    :: tmp
 !  integer             :: istat
 !  character(64) :: function_name = 'number_of_lines_in_file'
 ! !--------------------------------------------------------------------------
 !  if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name

 !  if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
!  end function number_of_lines_in_file
!
!==========================================================================
!==========================================================================
!
! function check_links(linkMatrix, accounted_for_links) result(missing_links)

 !  character(64) :: function_name = 'check_links'
 !  real(8), intent(in) :: linkMatrix(:,:)
 !  integer, intent(in) :: accounted_for_links(:)
 !  integer :: ii, jj, linkMatrix_count, accounted_count, missing_links

 ! !--------------------------------------------------------------------------
 !  if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name

 !  if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
!  end function check_links
!
!==========================================================================
!==========================================================================
!
subroutine reorganize_arrays()
 ! !
 ! ! this function is used to convert the null values that are defaulted in the
 ! ! B_nr_directweight_u and B_nr_totalweight_u columns of the nodes array into float
 ! ! zeros.
 character(64) :: subroutine_name = 'reorganize_arrays'
 !  real(8), intent (in) :: subnetwork_container_nodes (:,:,:)
 !  real(8), intent (in) :: subnetwork_container_links (:,:,:)
 !  real(8), intent(in out) :: nodeMatrix(:,:), linkMatrix(:,:)
 !  integer, intent(in) :: multiprocessors

 !  real(8), allocatable :: reorganizedNodes(:,:)
 !  real(8), allocatable :: reorganizedLinks(:,:)

 !  integer :: nodesRowCounter = 1, linksRowCounter = 1
 !  integer :: ii, jj, mp

 ! !--------------------------------------------------------------------------
  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name


  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
end subroutine reorganize_arrays

 end module BIPquick
!==========================================================================

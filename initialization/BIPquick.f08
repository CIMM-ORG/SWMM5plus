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
! !  subroutine network_node_preprocessing()
!     ! ----------------------------------------------------------------------------------------------------------------
!     !
!     ! Description:
!     !   This subroutine is a preprocessing step.
!     !   The network is traversed once and the B_nodeI array is populated with the adjacent upstream nodes.  This
!     !   step eliminates the need to continuously jump back and forth between the nodeI/linkI arrays during subsequent
!     !   network traversals.
!     !
!     ! ----------------------------------------------------------------------------------------------------------------
!     character(64) :: subroutine_name = 'network_node_preprocessing'
!     integer :: upstream_link, upstream_node, uplink_counter, root_node
!     integer ii, jj, uplinks
!     ! ----------------------------------------------------------------------------------------------------------------
!      if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
    
!      do ii= 1,size(nodeI,1)
!         if ( mod(ii, 1) == 0 ) then
!             print*, "processing upstream nodes of ni_idx:", ii
!         endif
!         root_node = nodeI(ii, ni_idx)
!         uplink_counter = nodeI(ii, ni_N_link_u)
!         print*, nodeI(ii, ni_idx), "has", nodeI(ii, ni_N_link_u), "upstream links"
!         do uplinks= 1, uplink_counter
!             if ( nodeI(ii, ni_idx_base1 + uplinks) /= nullValueI ) then
!                 upstream_link = nodeI(ii, ni_idx_base1 + uplinks)
!                 do jj= 1,size(linkI,1)
!                     if ( linkI(jj, li_idx) == upstream_link ) then
!                         upstream_node = linkI(jj, li_Mnode_u)
!                         B_nodeI(ii, uplinks) = upstream_node
!                     endif
!                 enddo
!             endif
!         enddo
!      enddo
    
!      if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
!      end subroutine network_node_preprocessing
    !
    !============================================================================ 
    !============================================================================ 
    ! 
    

 end module BIPquick
!==========================================================================

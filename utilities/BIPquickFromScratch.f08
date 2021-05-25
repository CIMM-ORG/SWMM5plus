 module BIPquickFromScratch

!%***********************************************************************************
!%
!% Description:
!%   This module holds the public BIPquick_partitioning subroutine that is called by the 
!%   partitioning.f08 utility module.  The private subroutines are used by the BIPquick
!%   main to do the fairly complicated stuff it does.
!%
!%***********************************************************************************

 use assign_index
 use globals
 use setting_definition, only: setting

 implicit none

 private

 public :: BIPquick_partitioning

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

 subroutine BIPquick_partitioning()
    print*, "Entered BIPquick_partitioning"
    call allocate_BIPquick_arrays()


    call deallocate_BIPquick_arrays()
    print*, "Leaving BIPquick_partitioning"
 end subroutine BIPquick_partitioning

 subroutine allocate_BIPquick_arrays()
    allocate(B_nodeI(size(nodeI,1), max_us_branch_per_node))
    B_nodeI(:,:) = nullValueI
    allocate(B_nodeR(size(nodeR,1), twoI))
    B_nodeR(:,:) = nullValueR
 end subroutine allocate_BIPquick_arrays

 subroutine deallocate_BIPquick_arrays()
    deallocate(B_nodeI)
    deallocate(B_nodeR)
 end subroutine deallocate_BIPquick_arrays

 end module BIPquickFromScratch
!==========================================================================

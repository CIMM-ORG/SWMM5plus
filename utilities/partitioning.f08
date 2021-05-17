module partitioning

    use array_index
    use data_keys
    use globals
    use setting_definition, only: setting
    use BIPquick
    use initialization

    implicit none

!***********************************************************************************
!
! Description:
!   This module controls the partitioning algorithm used.  Currently the options are
!       - BIPquick
!       - Default
!
!***********************************************************************************

    private

    public :: execute_partitioning
   
    integer, parameter :: B_ni_idx_Partition = 1 ! the node index number
    integer, parameter :: B_ni_Partition_No = 2 ! the Partition number to which that node index belongs
    integer, parameter :: B_ni_is_boundary = 3 ! a binary marker that is 1 when the node is shared between partitions in the link-node paradigm
   
    integer, parameter :: B_li_idx_Partition = 1 ! the link index number
    integer, parameter :: B_li_Partition_No = 2 ! the Partition number to which that link index belongs

    integer, pointer :: setP_N_images => setting%Partitioning%N_Image
    logical, pointer :: setP_UseBIPquick => setting%Partitioning%UseBIPquick
    logical, pointer :: setP_UseDefault => setting%Partitioning%UseDefault

contains


subroutine execute_partitioning()
    logical :: partition_correct
    ! --- Testing Partitioning Module
   call partitioning_algorithm_check()

   if ( setP_UseDefault .eqv. .true. ) then
        call default_partitioning()
   else if ( setP_UseBIPquick .eqv. .true. ) then
        call BIPquick_YJunction_Hardcode()
   end if

    ! This subroutine checks to see if the default partitioning is working correctly for the hard-coded case
   partition_correct = default_performance_check()

   print*, "*** partitioning is complete", partition_correct
end subroutine 


subroutine partitioning_algorithm_check()
    ! print*, setting%Partitioning%UseBIPquick, setting%Partitioning%UseDefault
    if ( (setP_UseBIPquick .eqv. .true.) .and. (setP_UseDefault .eqv. .true.) )  then
        print*, "There are two partitioning algorithms being used"
        stop
    else if ( (setP_UseBIPquick .eqv. .false.) .and. (setP_UseDefault .eqv. .false.) ) then
        print*, "No partitioning algorithms have been specified, default partitioning will be used"
        setP_UseDefault = .true.
    else
        if ( setP_UseBIPquick .eqv. .true. ) then
            print*, "Using BIPquick Partitioning"
        else if ( setP_UseDefault .eqv. .true. ) then
            print*, "Using Default Partitioning"
        end if
    end if
end subroutine partitioning_algorithm_check


subroutine default_partitioning()
    integer :: ii, jj, N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2
    integer :: total_num_elements, num_attributed_elements, assigning_image
    integer :: current_node_image, adjacent_link_image
    integer, allocatable, dimension(:) :: adjacent_links
    real(8) :: partition_threshold
    logical :: partition_correct

! ----------------------------------------------------------------------------------------------------------------
    ! This subroutine populates the P_nodeI, P_linkI arrays
    ! Rather than applying the BIPquick routine, the default partitioning is going to work by
    !   - Counting the total number of elements expected, using that to calculate the partition threshold
    !   - Iterating through the nodeI array until the number of elements expected exceeds the partition threshold
    !   - Iterating through the linkI array until the number of elements expected exceeds the partition threshold
    !   - Iterating again through the nodeI array to determine if the adjacent links are on different processors
! -----------------------------------------------------------------------------------------------------------------
    ! Determines the number of nodes of each type for the purpose of calculating partition threshold
    call count_node_types(N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2)
    ! print*, num_nJm_nodes, num_one_elem_nodes, num_zero_elem_nodes

    ! print*, sum(linkI(:, li_N_element))

    ! HACK The total number of elements is the sum of the elements from the links, plus the number of each node_type
    ! multiplied by how many elements are expected for that node_type
    total_num_elements = sum(linkI(:, li_N_element)) + (N_nBCup * N_elem_nBCup) + (N_nBCdn * N_elem_nBCdn) + &
        (N_nJm * N_elem_nJm) + (N_nStorage * N_elem_nStorage) + (N_nJ2 * N_elem_nJ2)
    partition_threshold = total_num_elements / real(setP_N_images)
    ! print*, total_num_elements, setP_N_images, partition_threshold

    ! This loop counts the elements attributed to each link, and assigns the link to an image
    num_attributed_elements = 0
    assigning_image = 1
    do ii = 1, size(linkI, 1)
        num_attributed_elements = num_attributed_elements + linkI(ii, li_N_element)
        ! If the number of elements is greater than the partition threshold, reset the number of elements and increment the image
        if ( num_attributed_elements > partition_threshold) then
            num_attributed_elements = 0
            ! This is a check to make sure that links aren't added to an image that doesn't exist
            if ( assigning_image /= setP_N_images ) then 
                assigning_image = assigning_image + 1
            end if
        end if

        P_linkI(ii, B_li_idx_Partition) = linkI(ii, ni_idx)
        P_linkI(ii, B_li_Partition_No) = assigning_image
    end do

    ! This loop counts the elements attributed to each node, and assigns the node to an image
    ! It also determines if that node has an adjacent link on a different image
    allocate(adjacent_links(6))
    do ii = 1, size(nodeI, 1)
        if ( (nodeI(ii, ni_node_type) == nBCup) &
            .or. (nodeI(ii, ni_node_type) == nBCdn) &
            .or. (nodeI(ii, ni_node_type) == nStorage) ) then
            num_attributed_elements = num_attributed_elements + 1
        else if ( nodeI(ii, ni_node_type) == nJm ) then
            num_attributed_elements = num_attributed_elements + 7
        end if

        ! If the number of attributed nodes exceeds the partition_threshold, then the remaining nodes are assigned to a new image
        if ( num_attributed_elements > partition_threshold) then
            num_attributed_elements = 0
            ! This is a check to make sure that nodes aren't added to an image that doesn't exist
            if ( assigning_image /= setP_N_images ) then 
                assigning_image = assigning_image + 1
            end if
        end if

        ! Fills in the P_nodeI array
        P_nodeI(ii, B_ni_idx_Partition) = nodeI(ii, ni_idx)
        P_nodeI(ii, B_ni_Partition_No) = assigning_image
        P_nodeI(ii, B_ni_is_boundary) = 0

        ! This bit of code checks the current node image, and compares it to the images of the adjacent links
        current_node_image = P_nodeI(ii, B_ni_Partition_No)
        adjacent_links = nodeI(ii, ni_Mlink_u1:ni_Mlink_d3)
        do jj = 1, size(adjacent_links)
            if ( adjacent_links(jj) == nullValueI ) then
                cycle
            end if
            adjacent_link_image = P_linkI(adjacent_links(jj), B_li_Partition_No)
            ! If the adjacent link and current node are on different images, then that node is a boundary
            if ( adjacent_link_image /= current_node_image ) then
                P_nodeI(ii, B_ni_is_boundary) = 1
            end if
        end do
    end do

end subroutine default_partitioning


function default_performance_check() result(partition_correct)
    integer :: ii
    logical :: partition_correct
    integer, allocatable, dimension(:,:) :: PartCheck_nodeI, PartCheck_linkI
    logical, allocatable, dimension(:,:) :: ArraySame_nodeI, ArraySame_linkI

    ! Allocate and initialize the correct partition arrays to be checked against
    allocate(PartCheck_nodeI(size(nodeI,1), B_ni_is_boundary))
    allocate(PartCheck_linkI(size(linkI,1), B_li_Partition_No))
    allocate(ArraySame_nodeI(size(nodeI,1), B_ni_is_boundary))
    allocate(ArraySame_linkI(size(linkI,1), B_li_Partition_No))

    PartCheck_nodeI(:, B_ni_idx_Partition) = (/1,2,3,4/) 
    PartCheck_nodeI(:, B_ni_Partition_No) = (/2,2,3,3/)
    PartCheck_nodeI(:, B_ni_is_boundary) = (/1,0,1,1/)

    PartCheck_linkI(:, B_li_idx_Partition) = (/1,2,3/)
    PartCheck_linkI(:, B_li_Partition_No) = (/1,2,2/)

    ! Assume that the partition arrays are going to be incorrect
    partition_correct = .false.

    ! Create a logical array of if the two arrays match
    ArraySame_nodeI(:,:) = ( PartCheck_nodeI == P_nodeI )
    ArraySame_linkI(:,:) = ( PartCheck_linkI == P_linkI )
    
    if ( all(ArraySame_nodeI) .eqv. .true. ) then 
        print*, "The node arrays are partitioned correctly"
    else
        print*, "There is a mistake in the P_nodeI"
        print*, ArraySame_nodeI(:,:)
        do ii = 1, size(P_nodeI, 1)
            print*, P_nodeI(ii, :)
            print*, PartCheck_nodeI(ii, :)
       end do
    end if

    if ( all(ArraySame_linkI) .eqv. .true. ) then 
        print*, "The link arrays are partitioned correctly"
    else
        print*, "There is a mistake in the P_linkI"
        print*, ArraySame_linkI(:,:)
        do ii = 1, size(P_linkI, 1)
            print*, P_linkI(ii, :)
            print*, PartCheck_linkI(ii, :)
       end do
    end if

    if ( (all(ArraySame_nodeI) .eqv. .true.) .and. (all(ArraySame_linkI) .eqv. .true.) ) then
        partition_correct = .true.
    end if

end function default_performance_check


end module partitioning
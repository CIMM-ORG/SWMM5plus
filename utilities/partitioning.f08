module partitioning

    use assign_index
    use data_keys
    use globals
    use setting_definition, only: setting
    use BIPquickFromScratch
    use initialization

    implicit none

!%***********************************************************************************
!%
!% Description:
!%   This module controls the partitioning algorithm used.  Currently the options are
!%       - BIPquick
!%       - Default
!%       - Random
!%
!%***********************************************************************************

    private

    public :: execute_partitioning

    integer, pointer :: setP_N_images => setting%Partitioning%N_Image

contains

!
!==========================================================================
!==========================================================================
!

subroutine execute_partitioning()
    ! --------------------------------------------------------
    !
    ! Description:
    !   The purpose of this subroutine is to check which partitioning
    !   algorithm should be used, then call that algorithm, then 
    !   check that the output is correct (if debug == true)
    !
    !---------------------------------------------------------
    logical :: partition_correct
    integer :: connectivity, ii
    real(8) :: part_size_balance
    character(64) :: subroutine_name = 'execute_partitioning'

    !% --------------------------------------------------------

    !% Determine which partitioning method is being used
    if (setting%Partitioning%PartitioningMethod == Default) then
        if (setting%Verbose) print*, "Using Default Partitioning"
        call default_partitioning()
    else if (setting%Partitioning%PartitioningMethod == BQuick) then
        if (setting%Verbose) print*, "Using BIPquick Partitioning Check"
        call BIPquick_partitioning()
    else if (setting%Partitioning%PartitioningMethod == Random) then
        if (setting%Verbose) print*, "Using Random Partitioning"
        call random_partitioning()
    end if

    if (setting%Debug%File%partitioning) then
        print *, '*** leave ', subroutine_name

        do ii = 1, size(P_nodeI, 1)
            print*, P_nodeI(ii, :)
        end do
        do ii = 1, size(P_linkI, 1)
            print*, P_linkI(ii, :)
        end do

        !% This subroutine checks to see if the default partitioning is working correctly for the hard-coded case
        ! partition_correct = default_performance_check()
        connectivity = partition_diagnostic_connectivity()
        part_size_balance = partition_diagnostic_partsizebalance()

        print*, "*** partitioning is complete", connectivity, part_size_balance
    end if
end subroutine

!
!==========================================================================
!==========================================================================
!

subroutine default_partitioning()
    ! ----------------------------------------------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine populates the P_nodeI, P_linkI arrays
    !   Rather than applying the BIPquick routine, the default partitioning is going to work by
    !       - Counting the total number of elements expected, using that to calculate the partition threshold
    !       - Iterating through the nodeI array until the number of elements expected exceeds the partition threshold
    !       - Iterating through the linkI array until the number of elements expected exceeds the partition threshold
    !       - Iterating again through the nodeI array to determine if the adjacent links are on different processors
    !
    ! -----------------------------------------------------------------------------------------------------------------
    integer :: ii, jj, N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2
    integer :: total_num_elements, num_attributed_elements, assigning_image
    integer :: current_node_image, adjacent_link_image
    integer, allocatable, dimension(:) :: adjacent_links
    real(8) :: partition_threshold
    logical :: partition_correct

!% ----------------------------------------------------------------------------------------------------------------
    !% This subroutine populates the P_nodeI, P_linkI arrays
    !% Rather than applying the BIPquick routine, the default partitioning works by
    !%   - Counting the total number of elements expected, using that to calculate the partition threshold
    !%   - Iterating through the nodeI array until the number of elements expected exceeds the partition threshold
    !%   - Iterating through the linkI array until the number of elements expected exceeds the partition threshold
    !%   - Iterating again through the nodeI array to determine if the adjacent links are on different processors
!% -----------------------------------------------------------------------------------------------------------------
    !% Determines the number of nodes of each type for the purpose of calculating partition threshold
    call count_node_types(N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2)


    !% HACK The total number of elements is the sum of the elements from the links, plus the number of each node_type
    !% multiplied by how many elements are expected for that node_type
    total_num_elements = sum(linkI(:, li_N_element)) + (N_nBCup * N_elem_nBCup) + (N_nBCdn * N_elem_nBCdn) + &
        (N_nJm * N_elem_nJm) + (N_nStorage * N_elem_nStorage) + (N_nJ2 * N_elem_nJ2)
    partition_threshold = total_num_elements / real(setP_N_images)

    !% This loop counts the elements attributed to each link, and assigns the link to an image
    num_attributed_elements = 0
    assigning_image = 1
    do ii = 1, size(linkI, 1)

        !% The num_attributed elements is incremented by the li_N_element for that link
        num_attributed_elements = num_attributed_elements + linkI(ii, li_N_element)

        !% The link is added to the P_linkI array
        P_linkI(ii, P_li_idx_Partition) = linkI(ii, ni_idx)
        P_linkI(ii, P_li_Partition_No) = assigning_image

        !% If the link is the last link, we need to not reset the num_attributed_elem going into the nodes loop
        if ( ii == size(linkI, 1) ) then 

            !% The last link is assigned to the current image and the link do-loop is exited
            P_linkI(ii, P_li_idx_Partition) = linkI(ii, ni_idx)
            P_linkI(ii, P_li_Partition_No) = assigning_image
            exit
        end if

        !% If the number of elements is greater than the partition threshold, reset the number of elements and increment the image
        if ( (num_attributed_elements > partition_threshold) ) then
            !% This is a check to make sure that links aren't added to an image that doesn't exist
            if ( assigning_image /= setP_N_images ) then
                assigning_image = assigning_image + 1
            end if
            num_attributed_elements = 0
        end if
    end do

    !% This loop counts the elements attributed to each node, and assigns the node to an image
    !% It also determines if that node has an adjacent link on a different image
    allocate(adjacent_links(6))
    do ii = 1, size(nodeI, 1)

        !% This if statement increments the num_attributed_elements by the number of elements associated with that node type
        if ( nodeI(ii, ni_node_type) == nBCup ) then
            num_attributed_elements = num_attributed_elements + N_elem_nBCup
        else if ( nodeI(ii, ni_node_type) == nBCdn ) then
            num_attributed_elements = num_attributed_elements + N_elem_nBCdn
        else if ( nodeI(ii, ni_node_type) == nStorage ) then
            num_attributed_elements = num_attributed_elements + N_elem_nStorage
        else if ( nodeI(ii, ni_node_type) == nJ2 ) then
            num_attributed_elements = num_attributed_elements + N_elem_nJ2
        else if ( nodeI(ii, ni_node_type) == nJm ) then
            num_attributed_elements = num_attributed_elements + N_elem_nJm
        end if

        !% If the number of attributed nodes exceeds the partition_threshold, then the remaining nodes are assigned to a new image
        if ( num_attributed_elements > partition_threshold) then
            num_attributed_elements = 0
            !% This is a check to make sure that nodes aren't added to an image that doesn't exist
            if ( assigning_image /= setP_N_images ) then
                assigning_image = assigning_image + 1
            end if
        end if

        !% Fills in the P_nodeI array
        P_nodeI(ii, P_ni_idx_Partition) = nodeI(ii, ni_idx)
        P_nodeI(ii, P_ni_Partition_No) = assigning_image
        P_nodeI(ii, P_ni_is_boundary) = 0

        !% This bit of code checks the current node image, and compares it to the images of the adjacent links
        current_node_image = P_nodeI(ii, P_ni_Partition_No)
        adjacent_links = nodeI(ii, ni_Mlink_u1:ni_Mlink_d3)
        do jj = 1, size(adjacent_links)
            if ( adjacent_links(jj) == nullValueI ) then
                cycle
            end if
            adjacent_link_image = P_linkI(adjacent_links(jj), P_li_Partition_No)
            !% If the adjacent link and current node are on different images, then that node is a boundary
            if ( adjacent_link_image /= current_node_image ) then
                P_nodeI(ii, P_ni_is_boundary) = P_nodeI(ii, P_ni_is_boundary) + 1
            end if
        end do
    end do
    deallocate(adjacent_links)

end subroutine default_partitioning

!
!==========================================================================
!==========================================================================
!

subroutine random_partitioning()
    integer :: ii, jj, N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2
    integer :: total_num_elements, num_attributed_elements, assigning_image
    integer :: current_node_image, adjacent_link_image
    integer, allocatable, dimension(:) :: adjacent_links, elem_per_image
    logical, allocatable, dimension(:) :: image_full
    real(8) :: partition_threshold, rand_num

    !% ----------------------------------------------------------------------------------------------------------------
    !% This subroutine populates the P_nodeI, P_linkI arrays
    !% Rather than applying the BIPquick routine, the default random partitioning works by
    !%   - Counting the total number of elements expected, using that to calculate the partition threshold
    !%   - Iterating through the linkI array, randomly assigning links to images and counting the number of elements ascribed to each image
    !%   - Iterating through the nodeI array, randomly assigning nodes to images and counting the number of elements ascribed to each image
    !%   - Iterating again through the nodeI array to determine if the adjacent links are on different processors
    !% -----------------------------------------------------------------------------------------------------------------

    !% Determines the number of nodes of each type for the purpose of calculating partition threshold
    call count_node_types(N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2)

    !% HACK The total number of elements is the sum of the elements from the links, plus the number of each node_type
    !% multiplied by how many elements are expected for that node_type
    total_num_elements = sum(linkI(:, li_N_element)) + (N_nBCup * N_elem_nBCup) + (N_nBCdn * N_elem_nBCdn) + &
        (N_nJm * N_elem_nJm) + (N_nStorage * N_elem_nStorage) + (N_nJ2 * N_elem_nJ2)
    partition_threshold = ( total_num_elements / real(setP_N_images) )

    !% Initialize the arrays that will hold the number of elements already on an image (and whether that image is full)
    allocate(elem_per_image(setP_N_images))
    elem_per_image(:) = 0
    allocate(image_full(setP_N_images))
    image_full(:) = .false.

    !% This loop counts the elements attributed to each link, and assigns the link to an image
    do ii = 1, size(linkI, 1)

        !% Calculates a random number and maps it onto an image number
        call random_number(rand_num)
        assigning_image = int(rand_num*setP_N_images) + 1

        !% If the image number selected is already full, pick a new number
        do while ( image_full(assigning_image) .eqv. .true. )
            call random_number(rand_num)
            assigning_image = int(rand_num*setP_N_images) + 1
        end do

        !% elem_per_image is incremented by li_N_element for the current link
        elem_per_image(assigning_image) = elem_per_image(assigning_image) + linkI(ii, li_N_element)

        !% If the number of elements is greater than the partition threshold, that image number is closed
        !% Note, this check after the assigning_image has been selected allows for images be over-filled
        if ( elem_per_image(assigning_image) > partition_threshold ) then
            image_full(assigning_image) = .true.
        end if

        !% Assign the link to the current image
        P_linkI(ii, P_li_idx_Partition) = linkI(ii, ni_idx)
        P_linkI(ii, P_li_Partition_No) = assigning_image
    end do

    !% This loop counts the elements attributed to each node, and assigns the node to an image
    do ii = 1, size(nodeI, 1)

        !% Calculates a random number and maps it onto an image number
        call random_number(rand_num)
        assigning_image = int(rand_num*setP_N_images) + 1

        !% If the image number selected is already full, pick a new number
        do while ( image_full(assigning_image) .eqv. .true. )
            call random_number(rand_num)
            assigning_image = int(rand_num*setP_N_images) + 1
        end do

        !% elem_per_image is incremented by the number of elements associated with each node type
        if ( nodeI(ii, ni_node_type) == nBCup ) then
            elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nBCup
        else if ( nodeI(ii, ni_node_type) == nBCdn ) then
            elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nBCdn
        else if ( nodeI(ii, ni_node_type) == nStorage ) then
            elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nStorage
        else if ( nodeI(ii, ni_node_type) == nJ2 ) then
            elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nJ2
        else if ( nodeI(ii, ni_node_type) == nJm ) then
            elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nJm
        end if
        
        !% If the number of elements is greater than the partition threshold, that image number is closed
        !% Note, this check after the assigning_image has been selected allows for images be over-filled 
        if ( elem_per_image(assigning_image) > partition_threshold ) then
            image_full(assigning_image) = .true.
        end if

        !% Fills in the P_nodeI array
        P_nodeI(ii, P_ni_idx_Partition) = nodeI(ii, ni_idx)
        P_nodeI(ii, P_ni_Partition_No) = assigning_image
        P_nodeI(ii, P_ni_is_boundary) = 0

        !% This bit of code checks the current node image, and compares it to the images of the adjacent links
        current_node_image = P_nodeI(ii, P_ni_Partition_No)
        adjacent_links = nodeI(ii, ni_Mlink_u1:ni_Mlink_d3)
        do jj = 1, size(adjacent_links)
            if ( adjacent_links(jj) == nullValueI ) then
                cycle
            end if
            adjacent_link_image = P_linkI(adjacent_links(jj), P_li_Partition_No)
            !% If the adjacent link and current node are on different images, then that node is a boundary
            if ( adjacent_link_image /= current_node_image ) then
                P_nodeI(ii, P_ni_is_boundary) = P_nodeI(ii, P_ni_is_boundary) + 1
            end if
        end do
    end do
    deallocate(adjacent_links)
    deallocate(elem_per_image)
    deallocate(image_full)

end subroutine random_partitioning

!    
!==========================================================================   
!========================================================================== 
!

function partition_diagnostic_partsizebalance() result(part_size_balance)
    !% This function is used to determine the discrepancy between the images with the largest and smallest number of elements
    integer :: part_size_balance
    integer :: ii, current_image, max_elem, min_elem
    integer, allocatable, dimension(:) :: elem_per_image

    !% Reallocate the elem_per_image array because not all algorithms use this
    allocate(elem_per_image(setP_N_images))
    elem_per_image(:) = 0

    !% Iterate through the P_linkI array
    do ii = 1, size(P_linkI, 1)

        !% The current image is the one to which the current link has been assigned
        current_image = P_linkI(ii, P_li_Partition_No)
        
        !% Iterate the number of elements for the current image by li_N_element for that link
        elem_per_image(current_image) = elem_per_image(current_image) + linkI(ii, li_N_element)
    end do

    !% Iterate through the P_nodeI array
    do ii = 1, size(P_nodeI, 1)

        !% The current image is the one to which the current link has been assigned
        current_image = P_nodeI(ii, P_ni_Partition_No)

        !% elem_per_image for the current image is incremented by the number of elements associated with each node type
        if ( nodeI(ii, ni_node_type) == nBCup ) then
            elem_per_image(current_image) = elem_per_image(current_image) + N_elem_nBCup
        else if ( nodeI(ii, ni_node_type) == nBCdn ) then
            elem_per_image(current_image) = elem_per_image(current_image) + N_elem_nBCdn
        else if ( nodeI(ii, ni_node_type) == nStorage ) then
            elem_per_image(current_image) = elem_per_image(current_image) + N_elem_nStorage
        else if ( nodeI(ii, ni_node_type) == nJ2 ) then
            elem_per_image(current_image) = elem_per_image(current_image) + N_elem_nJ2
        else if ( nodeI(ii, ni_node_type) == nJm ) then
            elem_per_image(current_image) = elem_per_image(current_image) + N_elem_nJm
        end if
    end do

    !% The maximum and minimum number of elements for an image are determined
    max_elem = maxval(elem_per_image(:))
    min_elem = minval(elem_per_image(:))
    print*, "Elem_per_image", elem_per_image(:)

    !% The difference between the max and min number of elements per image is the part_size_balance objective metric
    part_size_balance = max_elem - min_elem
        
end function partition_diagnostic_partsizebalance

!    
!==========================================================================   
!========================================================================== 
!

function partition_diagnostic_connectivity() result(connectivity)
    !% This function is used to determine how many boundary nodes exist in this partition set
    integer :: connectivity

    !% The number of boundary nodes in P_nodeI is the connectivity objective metric
    !% HACK - I might need to rethink this as some boundary nodes contributed more than 1 connection point
    connectivity = sum(P_nodeI(:, P_ni_is_boundary))
end function partition_diagnostic_connectivity

!    
!==========================================================================   
!========================================================================== 
!

function default_performance_check() result(partition_correct)
    integer :: ii
    logical :: partition_correct
    integer, allocatable, dimension(:,:) :: PartCheck_nodeI, PartCheck_linkI
    logical, allocatable, dimension(:,:) :: ArraySame_nodeI, ArraySame_linkI

    !% Allocate and initialize the correct partition arrays to be checked against
    allocate(PartCheck_nodeI(size(nodeI,1), P_ni_is_boundary))
    allocate(PartCheck_linkI(size(linkI,1), P_li_Partition_No))
    allocate(ArraySame_nodeI(size(nodeI,1), P_ni_is_boundary))
    allocate(ArraySame_linkI(size(linkI,1), P_li_Partition_No))

    PartCheck_nodeI(:, P_ni_idx_Partition) = (/1,2,3,4/)
    PartCheck_nodeI(:, P_ni_Partition_No) = (/2,2,3,3/)
    PartCheck_nodeI(:, P_ni_is_boundary) = (/1,0,1,1/)

    PartCheck_linkI(:, P_li_idx_Partition) = (/1,2,3/)
    PartCheck_linkI(:, P_li_Partition_No) = (/1,2,2/)

    !% Assume that the partition arrays are going to be incorrect
    partition_correct = .false.

    !% Create a logical array of if the two arrays match
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

    deallocate(PartCheck_nodeI)
    deallocate(PartCheck_linkI)
    deallocate(ArraySame_nodeI)
    deallocate(ArraySame_linkI)

end function default_performance_check


end module partitioning
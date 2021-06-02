module partitioning

    use assign_index
    use data_keys
    use globals
    use setting_definition, only: setting
    use utility
    use BIPquickFromScratch

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

    public :: execute_partitioning, count_node_types

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

    call allocate_partitioning_arrays()

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
    else if (setting%Partitioning%PartitioningMethod == BLink) then
        if (setting%Verbose) print*, "Using Balanced Link Partitioning"
        call balanced_link_partitioning()
    end if
    if (setting%Debug%File%partitioning) then
        print *, '*** leave ', subroutine_name

        do ii = 1, size(nodeI, 1)
            print*, nodeI(ii, ni_idx), nodeI(ii, ni_P_image:ni_P_is_boundary)
        end do
        do ii = 1, size(linkI, 1)
            print*, linkI(ii, li_idx), linkI(ii, li_P_image)
        end do

        !% This subroutine checks to see if the default partitioning is working correctly for the hard-coded case
        ! partition_correct = default_performance_check()
        connectivity = partition_diagnostic_connectivity()
        part_size_balance = partition_diagnostic_partsizebalance()

        print*, "*** partitioning is complete", connectivity, part_size_balance
    end if

    call deallocate_partitioning_arrays()
end subroutine
!
!==========================================================================
!==========================================================================
!
subroutine allocate_partitioning_arrays()
    allocate(adjacent_links(max_us_branch_per_node+max_ds_branch_per_node))
    allocate(elem_per_image(num_images()))
    allocate(image_full(num_images()))
end subroutine allocate_partitioning_arrays
!
!==========================================================================
!==========================================================================
!
subroutine deallocate_partitioning_arrays()
    deallocate(adjacent_links)
    deallocate(elem_per_image)
    deallocate(image_full)
end subroutine deallocate_partitioning_arrays
!
!==========================================================================
!==========================================================================
!
subroutine default_partitioning()
    ! ----------------------------------------------------------------------------------------------------------------
    !
    ! Description:
    !   The default partitioning algorithm populates the partitioning columns of the linkI-nodeI arrays.  Rather than
    !   assigning links-nodes to images topologically (as in BIPquick), the default partitioning algorithm assigns them
    !   in the order in which they appear in the link-node arrays.
    !
    ! -----------------------------------------------------------------------------------------------------------------
    integer :: ii, jj, N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2
    integer :: total_num_elements, num_attributed_elements, assigning_image
    integer :: current_node_image, adjacent_link_image
    real(8) :: partition_threshold
    logical :: partition_correct

!% ----------------------------------------------------------------------------------------------------------------

    !% Determines the number of nodes of each type for the purpose of calculating partition threshold
    call count_node_types(N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2)

    !% HACK The total number of elements is the sum of the elements from the links, plus the number of each node_type
    !% multiplied by how many elements are expected for that node_type
    total_num_elements = sum(linkI(:, li_N_element)) + (N_nBCup * N_elem_nBCup) + (N_nBCdn * N_elem_nBCdn) + &
        (N_nJm * N_elem_nJm) + (N_nStorage * N_elem_nStorage) + (N_nJ2 * N_elem_nJ2)
    partition_threshold = total_num_elements / real(num_images())

    !% This loop counts the elements attributed to each link, and assigns the link to an image
    num_attributed_elements = 0
    assigning_image = 1
    do ii = 1, size(linkI, 1)

        !% The num_attributed elements is incremented by the li_N_element for that link
        num_attributed_elements = num_attributed_elements + linkI(ii, li_N_element)

        !% The link's P column is assigned
        linkI(ii, li_P_image) = assigning_image

        !% If the link is the last link, we need to not reset the num_attributed_elem going into the nodes loop
        if ( ii == size(linkI, 1) ) then

            !% The last link is assigned to the current image and the link do-loop is exited
            linkI(ii, li_P_image) = assigning_image
            exit
        end if

        !% If the number of elements is greater than the partition threshold, reset the number of elements and increment the image
        if ( (num_attributed_elements > partition_threshold) ) then
            !% This is a check to make sure that links aren't added to an image that doesn't exist
            if ( assigning_image /= num_images() ) then
                assigning_image = assigning_image + 1
            end if
            num_attributed_elements = 0
        end if
    end do

    !% This loop counts the elements attributed to each node, and assigns the node to an image
    !% It also determines if that node has an adjacent link on a different image
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
            if ( assigning_image /= num_images() ) then
                assigning_image = assigning_image + 1
            end if
        end if

        !% Fills in the nodeI array P columns
        nodeI(ii, ni_P_image) = assigning_image
        nodeI(ii, ni_P_is_boundary) = 0

        !% This bit of code checks the current node image, and compares it to the images of the adjacent links
        current_node_image = nodeI(ii, ni_P_image)
        adjacent_links = nodeI(ii, ni_Mlink_u1:ni_Mlink_d3)
        do jj = 1, size(adjacent_links)
            if ( adjacent_links(jj) == nullValueI ) then
                cycle
            end if
            adjacent_link_image = linkI(adjacent_links(jj), li_P_image)
            !% If the adjacent link and current node are on different images, then that node is a boundary
            if ( adjacent_link_image /= current_node_image ) then
                nodeI(ii, ni_P_is_boundary) = nodeI(ii, ni_P_is_boundary) + 1
            end if
        end do
    end do

end subroutine default_partitioning
!
!==========================================================================
!==========================================================================
!
subroutine random_partitioning()
    ! ----------------------------------------------------------------------------------------------------------------
    !
    ! Description:
    !   The random partitioning algorithm populates the partitioning columns of the linkI-nodeI arrays.  An alternative
    !   to the default partitioning algorithm, the random partitioning algorithm looks at each link and node and assigns
    !   it to a random image (after checking to ensure that image is not full).
    !
    ! -----------------------------------------------------------------------------------------------------------------
    integer :: ii, jj, N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2
    integer :: total_num_elements, num_attributed_elements, assigning_image
    integer :: current_node_image, adjacent_link_image
    real(8) :: partition_threshold, rand_num
    !% ----------------------------------------------------------------------------------------------------------------

    !% Determines the number of nodes of each type for the purpose of calculating partition threshold
    call count_node_types(N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2)

    !% HACK The total number of elements is the sum of the elements from the links, plus the number of each node_type
    !% multiplied by how many elements are expected for that node_type
    total_num_elements = sum(linkI(:, li_N_element)) + (N_nBCup * N_elem_nBCup) + (N_nBCdn * N_elem_nBCdn) + &
        (N_nJm * N_elem_nJm) + (N_nStorage * N_elem_nStorage) + (N_nJ2 * N_elem_nJ2)
    partition_threshold = ( total_num_elements / real(num_images()) )

    !% Initialize the arrays that will hold the number of elements already on an image (and whether that image is full)
    elem_per_image(:) = 0
    image_full(:) = .false.

    !% This loop counts the elements attributed to each link, and assigns the link to an image
    do ii = 1, size(linkI, 1)

        !% Calculates a random number and maps it onto an image number
        call random_number(rand_num)
        assigning_image = int(rand_num*num_images()) + 1

        !% If the image number selected is already full, pick a new number
        do while ( image_full(assigning_image) .eqv. .true. )
            call random_number(rand_num)
            assigning_image = int(rand_num*num_images()) + 1
        end do

        !% elem_per_image is incremented by li_N_element for the current link
        elem_per_image(assigning_image) = elem_per_image(assigning_image) + linkI(ii, li_N_element)

        !% If the number of elements is greater than the partition threshold, that image number is closed
        !% Note, this check after the assigning_image has been selected allows for images be over-filled
        if ( elem_per_image(assigning_image) > partition_threshold ) then
            image_full(assigning_image) = .true.
        end if

        !% Assign the link to the current image
        linkI(ii, li_P_image) = assigning_image
    end do

    !% This loop counts the elements attributed to each node, and assigns the node to an image
    do ii = 1, size(nodeI, 1)

        !% Calculates a random number and maps it onto an image number
        call random_number(rand_num)
        assigning_image = int(rand_num*num_images()) + 1

        !% If the image number selected is already full, pick a new number
        do while ( image_full(assigning_image) .eqv. .true. )
            call random_number(rand_num)
            assigning_image = int(rand_num*num_images()) + 1
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

        !% Assigns the nodes to an image, initializes the is_boundary check to 0
        nodeI(ii, ni_P_image) = assigning_image
        nodeI(ii, ni_P_is_boundary) = 0

        !% This bit of code checks the current node image, and compares it to the images of the adjacent links
        current_node_image = nodeI(ii, ni_P_image)
        adjacent_links = nodeI(ii, ni_Mlink_u1:ni_Mlink_d3)
        do jj = 1, size(adjacent_links)
            if ( adjacent_links(jj) == nullValueI ) then
                cycle
            end if
            adjacent_link_image = linkI(adjacent_links(jj), li_P_image)
            !% If the adjacent link and current node are on different images, then that node is a boundary
            if ( adjacent_link_image /= current_node_image ) then
                nodeI(ii, ni_P_is_boundary) = nodeI(ii, ni_P_is_boundary) + 1
            end if
        end do
    end do

end subroutine random_partitioning
!
!==========================================================================
!==========================================================================
!
subroutine balanced_link_partitioning()
    integer :: ii, jj
    integer :: clink, clink_image, assigned_image
    integer :: start_id, end_id
    integer :: count, remainder, rank

    if (N_link < num_images()) then
        call default_partitioning()
    else
        do rank = 0, num_images()-1
            count = N_link / num_images()
            remainder = mod(N_link, num_images())

            if (rank < remainder) then
                ! The first 'remainder' ranks get 'count + 1' tasks each
                start_id = rank * (count + 1)
                end_id = start_id + count
            else
                ! The remaining 'size - remainder' ranks get 'count' task each
                start_id = rank * count + remainder
                end_id = start_id + (count - 1)
            end if

            linkI(start_id+1:end_id+1, li_P_image) = rank+1
        end do
        do ii = 1, N_node
            assigned_image = nullvalueI
            do jj = 1, (max_us_branch_per_node + max_ds_branch_per_node)
                clink = nodeI(ii, ni_idx_base1+jj)
                if (clink /= nullvalueI) then
                    clink_image = linkI(clink, li_P_image)
                    if (clink_image < assigned_image) then
                        assigned_image = clink_image
                    end if
                end if
            end do
            nodeI(ii, ni_P_image) = assigned_image
        end do
    end if
end subroutine balanced_link_partitioning
!
!==========================================================================
!==========================================================================
!
function partition_diagnostic_partsizebalance() result(part_size_balance)
    ! ----------------------------------------------------------------------------------------------------------------
    !
    ! Description:
    !   This function is used to calculate the part size balance metric for the given partition set (i.e. the output
    !   from any of the partitioning algorithms).  The part size balance metric is equal to the greatest number of
    !   elements assigned to a single processor minus the smallest number of elements assigned to a single processor.
    !
    ! -----------------------------------------------------------------------------------------------------------------
    integer :: part_size_balance
    integer :: ii, current_image, max_elem, min_elem
    ! -----------------------------------------------------------------------------------------------------------------

    !% Reset the elem_per_image array to all zeros
    elem_per_image(:) = 0

    !% Iterate through the linkI array
    do ii = 1, size(linkI, 1)

        !% The current image is the one to which the current link has been assigned
        current_image = linkI(ii, li_P_image)

        !% Iterate the number of elements for the current image by li_N_element for that link
        elem_per_image(current_image) = elem_per_image(current_image) + linkI(ii, li_N_element)
    end do

    !% Iterate through the nodeI array
    do ii = 1, size(nodeI, 1)

        !% The current image is the one to which the current link has been assigned
        current_image = nodeI(ii, ni_P_image)

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
    ! ----------------------------------------------------------------------------------------------------------------
    !
    ! Description:
    !   This function is used to calculate the connectivity metric for the given partition set (i.e. the output
    !   from any of the partitioning algorithms).  The connectivity metric is equal to the sum of the ni_is_boundary
    !   column of the nodeI array after the partition has been completed.  Note: the ni_is_boundary column is given a
    !   value of 0 when the node is an internal partition node, and is incremented by 1 for each additional partition
    !   the node touches.
    !
    ! -----------------------------------------------------------------------------------------------------------------
    integer :: connectivity
    ! -----------------------------------------------------------------------------------------------------------------

    !% The sum of the ni_is_boundary column is the connectivity
    connectivity = sum(nodeI(:, ni_P_is_boundary))
end function partition_diagnostic_connectivity
!
!==========================================================================
!==========================================================================
!
! function default_performance_check() result(partition_correct)
!     integer :: ii
!     logical :: partition_correct
!     integer, allocatable, dimension(:,:) :: PartCheck_nodeI, PartCheck_linkI
!     logical, allocatable, dimension(:,:) :: ArraySame_nodeI, ArraySame_linkI

!     !% Allocate and initialize the correct partition arrays to be checked against
!     allocate(PartCheck_nodeI(size(nodeI,1), P_ni_is_boundary))
!     allocate(PartCheck_linkI(size(linkI,1), P_li_Partition_No))
!     allocate(ArraySame_nodeI(size(nodeI,1), P_ni_is_boundary))
!     allocate(ArraySame_linkI(size(linkI,1), P_li_Partition_No))

!     PartCheck_nodeI(:, P_ni_idx_Partition) = (/1,2,3,4/)
!     PartCheck_nodeI(:, P_ni_Partition_No) = (/2,2,3,3/)
!     PartCheck_nodeI(:, P_ni_is_boundary) = (/1,0,1,1/)

!     PartCheck_linkI(:, P_li_idx_Partition) = (/1,2,3/)
!     PartCheck_linkI(:, P_li_Partition_No) = (/1,2,2/)

!     !% Assume that the partition arrays are going to be incorrect
!     partition_correct = .false.

!     !% Create a logical array of if the two arrays match
!     ArraySame_nodeI(:,:) = ( PartCheck_nodeI == P_nodeI )
!     ArraySame_linkI(:,:) = ( PartCheck_linkI == P_linkI )

!     if ( all(ArraySame_nodeI) .eqv. .true. ) then
!         print*, "The node arrays are partitioned correctly"
!     else
!         print*, "There is a mistake in the P_nodeI"
!         print*, ArraySame_nodeI(:,:)
!         do ii = 1, size(P_nodeI, 1)
!             print*, P_nodeI(ii, :)
!             print*, PartCheck_nodeI(ii, :)
!        end do
!     end if

!     if ( all(ArraySame_linkI) .eqv. .true. ) then
!         print*, "The link arrays are partitioned correctly"
!     else
!         print*, "There is a mistake in the P_linkI"
!         print*, ArraySame_linkI(:,:)
!         do ii = 1, size(P_linkI, 1)
!             print*, P_linkI(ii, :)
!             print*, PartCheck_linkI(ii, :)
!        end do
!     end if

!     if ( (all(ArraySame_nodeI) .eqv. .true.) .and. (all(ArraySame_linkI) .eqv. .true.) ) then
!         partition_correct = .true.
!     end if

!     deallocate(PartCheck_nodeI)
!     deallocate(PartCheck_linkI)
!     deallocate(ArraySame_nodeI)
!     deallocate(ArraySame_linkI)

! end function default_performance_check


end module partitioning
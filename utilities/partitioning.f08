module partitioning

    use array_index
    use data_keys
    use globals
    use setting_definition, only: setting
    use BIPquick

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

    public :: default_partitioning, partitioning_algorithm_check
   
    integer, parameter :: B_ni_idx_Partition = 1 ! the node index number
    integer, parameter :: B_ni_Partition_No = 2 ! the Partition number to which that node index belongs
    integer, parameter :: B_ni_is_boundary = 3 ! a binary marker that is 1 when the node is shared between partitions in the link-node paradigm
   
    integer, parameter :: B_li_idx_Partition = 1 ! the link index number
    integer, parameter :: B_li_Partition_No = 2 ! the Partition number to which that link index belongs


contains

subroutine partitioning_algorithm_check()
    print*, setting%Partitioning%UseBIPquick, setting%Partitioning%UseDefault
    if ( (setting%Partitioning%UseBIPquick .eqv. .true.) .and. (setting%Partitioning%UseDefault .eqv. .true.) )  then
        print*, "There are two partitioning algorithms being used"
        stop
    else if ( (setting%Partitioning%UseBIPquick .eqv. .false.) .and. (setting%Partitioning%UseDefault .eqv. .false.) ) then
        print*, "No partitioning algorithms have been specified, default partitioning will be used"
        setting%Partitioning%UseDefault = .true.
    else
        if ( setting%Partitioning%UseBIPquick .eqv. .true. ) then
            print*, "Using BIPquick Partitioning"
        else if ( setting%Partitioning%UseDefault .eqv. .true. ) then
            print*, "Using Default Partitioning"
        end if
    end if
end subroutine partitioning_algorithm_check

subroutine default_partitioning()
    integer :: ii, jj, num_nJm_nodes, num_one_elem_nodes, num_zero_elem_nodes
    integer :: total_num_elements, num_attributed_elements, assigning_image
    integer :: current_node_image, adjacent_link_image
    integer, allocatable, dimension(:) :: adjacent_links
    real(8) :: partition_threshold

! ----------------------------------------------------------------------------------------------------------------
    ! This subroutine populates the P_nodeI, P_linkI arrays
    ! Rather than applying the BIPquick routine, the default partitioning is going to work by
    !   - Counting the total number of elements expected, using that to calculate the partition threshold
    !   - Iterating through the nodeI array until the number of elements expected exceeds the partition threshold
    !   - Iterating through the linkI array until the number of elements expected exceeds the partition threshold
    !   - Iterating again through the nodeI array to determine if the adjacent links are on different processors
! -----------------------------------------------------------------------------------------------------------------
    ! Determines the number of nodes of each type for the purpose of calculating partition threshold
    call count_node_types(num_nJm_nodes, num_one_elem_nodes, num_zero_elem_nodes)
    ! print*, num_nJm_nodes, num_one_elem_nodes, num_zero_elem_nodes

    ! HACK This is a temporary hardcode until Gerardo can populate this column from the CFL condition
    linkI(:, li_N_element) = 10
    ! print*, sum(linkI(:, li_N_element))

    ! HACK The total number of elements is the sum of the elements from the links, plus the number of each node_type
    ! multiplied by how many elements are expected for that node_type
    total_num_elements = sum(linkI(:, li_N_element)) + num_nJm_nodes*7 + num_one_elem_nodes*1 + num_zero_elem_nodes*0
    partition_threshold = total_num_elements / real(setting%Partitioning%Num_Images_Setting)
    ! print*, total_num_elements, setting%Partitioning%Num_Images_Setting, partition_threshold

    ! This loop counts the elements attributed to each link, and assigns the link to an image
    num_attributed_elements = 0
    assigning_image = 1
    do ii = 1, size(linkI, 1)
        num_attributed_elements = num_attributed_elements + linkI(ii, li_N_element)
        ! If the number of elements is greater than the partition threshold, reset the number of elements and increment the image
        if ( num_attributed_elements > partition_threshold) then
            num_attributed_elements = 0
            ! This is a check to make sure that links aren't added to an image that doesn't exist
            if ( assigning_image /= setting%Partitioning%Num_Images_Setting ) then 
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
            if ( assigning_image /= setting%Partitioning%Num_Images_Setting ) then 
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

subroutine count_node_types(num_nJm_nodes, num_one_elem_nodes, num_zero_elem_nodes)
    integer, intent(in out) :: num_nJm_nodes, num_one_elem_nodes, num_zero_elem_nodes
    integer :: ii, N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2

    ! This subroutine uses the vectorized count() function to search the array for number of instances of each node type
    N_nBCup = count(nodeI(:, ni_node_type) == nBCup)
    N_nBCdn = count(nodeI(:, ni_node_type) == nBCdn)
    N_nJm = count(nodeI(:, ni_node_type) == nJM)
    N_nStorage = count(nodeI(:, ni_node_type) == nStorage)
    N_nJ2 = count(nodeI(:, ni_node_type) == nJ2)

    ! The nodes that correspond to having 7, 1, and 0 attributed elements are summed together
    num_nJm_nodes = N_nJm
    num_one_elem_nodes = N_nBCup + N_nBCdn + N_nStorage
    num_zero_elem_nodes = N_nJ2

end subroutine count_node_types


end module partitioning
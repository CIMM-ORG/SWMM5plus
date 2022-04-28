module partitioning

    use define_keys
    use define_globals
    use define_indexes
    use define_settings, only: setting
    use discretization, only: init_discretization_nominal
    use utility
    use utility_allocate
    use BIPquick
    use utility_deallocate
    use utility_crash
    use utility_profiler

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

    public :: partitioning_toplevel

contains
!
!==========================================================================
!==========================================================================
!
subroutine partitioning_toplevel()
    ! --------------------------------------------------------
    ! Description:
    !   The purpose of this subroutine is to check which partitioning
    !   algorithm should be used, then call that algorithm, then
    !   check that the output is correct (if debug == true)
    !---------------------------------------------------------
        logical :: partition_correct
        integer :: connectivity, ii, nn
        real(8) :: part_size_balance
        character(64) :: subroutine_name = 'partitioning_toplevel'
    !% --------------------------------------------------------
        if (crashYN) return
    !% --------------------------------------------------------
    call util_count_node_types(N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2, N_nJ1)

    !% the allcoation should probably should only be for image=1
    !% but we need to then be careful with deallocation
    call util_allocate_partitioning_arrays() 

    !% --- check for using a single processor or multiprocessor
    if ( num_images() == 1 ) then
        node%I(:, ni_P_image) = oneI
        node%I(:, ni_P_is_boundary) = zeroI
        link%I(:, li_P_image) = oneI
        if (setting%Output%Verbose) print*, "... Using one processor, bypassing partitioning"
    else
        !% --- confine the partitioning computations to a single image
        if (this_image() == 1) then
            !% Determine which partitioning method is being used
            !print *   !% this is needed because SWMM-C doesn't have a newline after their last printout
            if (setting%Partitioning%PartitioningMethod == Default) then
                if (setting%Output%Verbose) write(*,"(A)")  "... using Default Partitioning..."
                call init_partitioning_default()
            else if (setting%Partitioning%PartitioningMethod == Random) then
                if (setting%Output%Verbose) write(*,"(A)")  "... using Random Partitioning..."
                call init_partitioning_random()
            else if (setting%Partitioning%PartitioningMethod == BLink) then
                if (setting%Output%Verbose) write(*,"(A)")  "... using Balanced Link Partitioning..."
                call init_partitioning_linkbalance()
            else if (setting%Partitioning%PartitioningMethod == BQuick) then
                if (setting%Output%Verbose) write(*,"(A)") "... using BIPquick Partitioning..."
                call BIPquick_toplevel()
                ! call init_partitioning_bquick_diagnostic ()
            else
                print *, "Error, partitioning method not supported"
                call util_crashpoint(87095)
                return
                !stop 
            end if
        end if

        call util_crashstop(44873)

        !% broadcast partitioning results to all images
        call co_broadcast(node%I, source_image=1)
        call co_broadcast(node%R, source_image=1)
        call co_broadcast(node%YN, source_image=1)
        call co_broadcast(link%I, source_image=1)
        call co_broadcast(link%R, source_image=1)
        call co_broadcast(link%YN, source_image=1)
        sync all
    end if

    if (setting%Debug%File%partitioning) then
        print *, "Node Partitioning"
        print *
        do ii = 1, size(node%I, 1)
            if ( ii <= N_node ) then
                print*, node%Names(ii)%str, node%I(ii, ni_idx), node%I(ii, ni_P_image:ni_P_is_boundary)
            else
                print*, node%I(ii, ni_idx), node%I(ii, ni_P_image:ni_P_is_boundary)
            endif
        end do

        print *, "Link Partitioning"
        print *
        do ii = 1, size(link%I, 1)
            if ( ii <= N_link ) then
                print*, link%Names(ii)%str, link%I(ii, li_idx), link%I(ii, li_P_image), link%I(ii, li_parent_link), &
                    link%I(ii, li_Mnode_u:li_Mnode_d)
            else
                print*, link%I(ii, li_idx), link%I(ii, li_P_image), link%I(ii, li_parent_link), &
                    link%I(ii, li_Mnode_u:li_Mnode_d)
            endif

        end do

        !% This subroutine checks to see if the default partitioning is working correctly for the hard-coded case
        ! partition_correct = default_performance_check()
        connectivity = init_partitioning_metric_connectivity()
        ! part_size_balance = init_partitioning_metric_partsizebalance()
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        write(*,"(2(A,i5),A)") &
        "completed partitioning (", connectivity, ") | [Processor ", this_image(), "]" ! part_size_balance
    end if

    N_node = count(node%I(:,ni_idx) /= nullvalueI)
    N_link = count(link%I(:,li_idx) /= nullvalueI)

    call util_deallocate_partitioning_arrays()

    call init_timer_stop()

end subroutine partitioning_toplevel
!
!==========================================================================
!==========================================================================
!
subroutine init_partitioning_bquick_diagnostic ()
    ! --------------------------------------------------------
    ! Description:
    !   The purpose of this subroutine is to delete phantom links
    !   and nodes so that the network remain intact
    !---------------------------------------------------------
    integer :: ii
    integer, pointer :: pNode, pLink, sLink, dnNidx, upNidx
    integer, pointer :: plinkUp, pLinkDn
    real(8), pointer :: nominalElemLength
    integer, dimension(:), allocatable, target :: nodeIndexes
    character(64) :: subroutine_name = 'init_partitioning_bquick_diagnostic'
!% --------------------------------------------------------
    if (crashYN) return
!% --------------------------------------------------------
    !% pointers
    nominalElemLength => setting%Discretization%NominalElemLength

    !% pack all the node indexes excluding the nullvalues
    nodeIndexes = pack(node%I(:,ni_idx), (node%I(:,ni_idx) /= nullvalueI))

    do ii = 1, size(nodeIndexes, 1)
        pNode => nodeIndexes(ii)
        if (node%YN(pNode,nYN_is_phantom_node)) then
            !% if both the upstream and downstream links are phantom
            if (link%YN(node%I(pNode,ni_Mlink_u1),lYN_isPhantomLink) .and. &
                link%YN(node%I(pNode,ni_Mlink_d1),lYN_isPhantomLink)) then
                plinkUp => node%I(pNode,ni_Mlink_u1)
                pLinkDn => node%I(pNode,ni_Mlink_d1)
                upNidx  => link%I(plinkUp,li_Mnode_u)
                dnNidx  => link%I(pLinkDn,li_Mnode_d)
            !% if the upstream link is a phantom link
            else if (link%YN(node%I(pNode,ni_Mlink_u1),lYN_isPhantomLink)) then
                pLink  => node%I(pNode,ni_Mlink_u1)
                sLink  => node%I(pNode,ni_Mlink_d1)
                upNidx => link%I(pLink,li_Mnode_u)
                dnNidx => link%I(sLink,li_Mnode_d)
            !% if the downstream link is a phantom link
            else if (link%YN(node%I(pNode,ni_Mlink_d1),lYN_isPhantomLink)) then
                sLink  => node%I(pNode,ni_Mlink_u1)
                pLink  => node%I(pNode,ni_Mlink_d1)
                upNidx => link%I(sLink,li_Mnode_u)
                dnNidx => link%I(pLink,li_Mnode_d)
            !% should not reach this error condition
            else
                print*, 'In subroutine', subroutine_name
                print*, 'Error: phantom node', pNode, 'doesnot have any up or dn phantom link'
                !stop 
                call util_crashpoint(147856)
                return
            end if

            !% print diagnistic of the spanning and phantom links
            if (this_image() == 1) then
                if (link%YN(node%I(pNode,ni_Mlink_u1),lYN_isPhantomLink) .and. &
                    link%YN(node%I(pNode,ni_Mlink_d1),lYN_isPhantomLink)) then
                    print*, 'Phantom link detected both upstream and downstream'
                    print*, pNode,                     ' = phantom node index'
                    print*, node%I(pNode,ni_P_image),  ' = phantom node image'
                    print*, plinkUp,                   ' = up phantom link index'
                    print*, link%R(plinkUp,lr_length), ' = up phantom link length'
                    print*, link%I(plinkUp,li_P_image),' = up phantom link in image'
                    print*, upNidx,                    ' = node up idx of the phantom link'
                    print*, node%I(upNidx,ni_P_image), ' = node up of phantom link in image'
                    print*, pLinkDn,                   ' = dn phantom link index'
                    print*, link%R(pLinkDn,lr_length), ' = dn phantom link length'
                    print*, link%I(pLinkDn,li_P_image),' = dn phantom link in image'
                    print*, dnNidx,                    ' = node dn idx of the spanning link'
                    print*, node%I(dnNidx,ni_P_image), ' = node dn of spanning link in image'
                    print*
                else if (link%YN(node%I(pNode,ni_Mlink_u1),lYN_isPhantomLink)) then
                    print*, 'Upstream phantom link detected'
                    print*, pNode,                     ' = phantom node index'
                    print*, node%I(pNode,ni_P_image),  ' = phantom node image'
                    print*, pLink,                     ' = up phantom link index'
                    print*, link%R(pLink,lr_length),   ' = up phantom link length'
                    print*, link%I(pLink,li_P_image),  ' = up phantom link in image'
                    print*, upNidx,                    ' = node up idx of the phantom link'
                    print*, node%I(upNidx,ni_P_image), ' = node up of phantom link in image'
                    print*, sLink,                     ' = dn spanning link index'
                    print*, link%R(sLink,lr_length),   ' = dn spanning link length'
                    print*, link%I(sLink,li_P_image),  ' = dn spanning link in image'
                    print*, dnNidx,                    ' = node dn idx of the spanning link'
                    print*, node%I(dnNidx,ni_P_image), ' = node dn of spanning link in image'
                    print*
                else if (link%YN(node%I(pNode,ni_Mlink_d1),lYN_isPhantomLink)) then
                    print*, 'Downstream phantom link detected'
                    print*, pNode,                     ' = phantom node index'
                    print*, node%I(pNode,ni_P_image),  ' = phantom node image'
                    print*, sLink,                     ' = up spanning link index'
                    print*, link%R(sLink,lr_length),   ' = up spanning link length'
                    print*, link%I(sLink,li_P_image),  ' = up spanning link in image'
                    print*, upNidx,                    ' = node idx up of the spanning link'
                    print*, node%I(upNidx,ni_P_image), ' = node up of spanning link in image'
                    print*, pLink,                     ' = dn phantom link index'
                    print*, link%R(pLink,lr_length),   ' = dn phantom link length'
                    print*, link%I(pLink,li_P_image),  ' = dn phantom link in image'
                    print*, dnNidx,                    ' = node idx dn of the phantom link'
                    print*, node%I(dnNidx,ni_P_image), ' = node dn of phantom link in image'
                    print*
                end if
            end if 

            !% adjust only when the spanning or phantom link is smaller 
            !% than the nominal element length and the phantom link and node
            !% are on different images

            !% HACK: if the phantom link and nodes are on a same image, 
            !% that implys
            ! if ( (link%R(pLink,lr_length) .le. nominalElemLength) .or. &
            !      (link%R(sLink,lr_length) .le. nominalElemLength)) then
            !     !% print out for the adjusted links    
            !     if (this_image() == 1) then
            !         print*, 'Adjusting Link', pLink, ' and ', sLink
            !     end if

            ! end if    
        end if
    end do 

    !% deacclocate the temporary array
    deallocate(nodeIndexes)

end subroutine init_partitioning_bquick_diagnostic
!
!==========================================================================
!==========================================================================
!
subroutine init_partitioning_default()
    ! ----------------------------------------------------------------------------------------------------------------
    !
    ! Description:
    !   The default partitioning algorithm populates the partitioning columns of the link%I-node%I arrays.  Rather than
    !   assigning links-nodes to images topologically (as in BIPquick), the default partitioning algorithm assigns them
    !   in the order in which they appear in the link-node arrays.
    !
    ! -----------------------------------------------------------------------------------------------------------------
    integer :: ii, jj, N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2, N_nJ1
    integer :: total_num_elements, num_attributed_elements, assigning_image
    integer :: current_node_image, adjacent_link_image
    real(8) :: partition_threshold
    logical :: partition_correct

    if (crashYN) return
    !% Determines the number of nodes of each type for the purpose of calculating partition threshold
    call util_count_node_types(N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2, N_nJ1)

    !% HACK The total number of elements is the sum of the elements from the links, plus the number of each node_type
    !% multiplied by how many elements are expected for that node_type
    total_num_elements = sum(link%I(:, li_N_element))         &
                             + (N_nBCup    * N_elem_nBCup)    &
                             + (N_nBCdn    * N_elem_nBCdn)    &
                             + (N_nJm      * N_elem_nJm)      &
                             + (N_nStorage * N_elem_nStorage) &
                             + (N_nJ2      * N_elem_nJ2)      &
                             + (N_nJ1      * N_elem_nJ1)                   !% brh 20211217
    partition_threshold = total_num_elements / real(num_images())

    !% This loop counts the elements attributed to each link, and assigns the link to an image
    num_attributed_elements = 0
    assigning_image = 1
    do ii = 1, size(link%I, 1)

        !% The num_attributed elements is incremented by the li_N_element for that link
        num_attributed_elements = num_attributed_elements + link%I(ii, li_N_element)

        !% The link's P column is assigned
        link%I(ii, li_P_image) = assigning_image

        !% If the link is the last link, we need to not reset the num_attributed_elem going into the nodes loop
        if ( ii == size(link%I, 1) ) then

            !% The last link is assigned to the current image and the link do-loop is exited
            link%I(ii, li_P_image) = assigning_image
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
    do ii = 1, size(node%I, 1)

        !% This if statement increments the num_attributed_elements by the number of elements associated with that node type

        ! if ( node%I(ii, ni_node_type) == nBCup ) then
        !     num_attributed_elements = num_attributed_elements + N_elem_nBCup
        ! else if ( node%I(ii, ni_node_type) == nBCdn ) then
        !     num_attributed_elements = num_attributed_elements + N_elem_nBCdn
        ! else if ( node%I(ii, ni_node_type) == nStorage ) then
        !     num_attributed_elements = num_attributed_elements + N_elem_nStorage
        ! else if ( node%I(ii, ni_node_type) == nJ2 ) then
        !     num_attributed_elements = num_attributed_elements + N_elem_nJ2
        ! else if ( node%I(ii, ni_node_type) == nJm ) then
        !     num_attributed_elements = num_attributed_elements + N_elem_nJm
        ! end if

        !% brh20211217 -- revised to add nJ1
        select case (node%I(ii, ni_node_type))
        case (nBCup)
            num_attributed_elements = num_attributed_elements + N_elem_nBCup
        case (nBCdn)
            num_attributed_elements = num_attributed_elements + N_elem_nBCdn
        case (nStorage)
            num_attributed_elements = num_attributed_elements + N_elem_nStorage
        case (nJ1)
            num_attributed_elements = num_attributed_elements + N_elem_nJ1
        case (nJ2)
            num_attributed_elements = num_attributed_elements + N_elem_nJ2
        case (nJM)
            num_attributed_elements = num_attributed_elements + N_elem_nJm
        case default 
            print *, 'CODE ERROR: unknown node type # of ',node%I(ii, ni_node_type)
            print *, 'which has key of ',trim(reverseKey(node%I(ii, ni_node_type)))
            !stop 
            call util_crashpoint(1098226)
            return
        end select


        !% If the number of attributed nodes exceeds the partition_threshold, then the remaining nodes are assigned to a new image
        if ( num_attributed_elements > partition_threshold) then
            num_attributed_elements = 0
            !% This is a check to make sure that nodes aren't added to an image that doesn't exist
            if ( assigning_image /= num_images() ) then
                assigning_image = assigning_image + 1
            end if
        end if
        !% Fills in the node%I array P columns
        node%I(ii, ni_P_image) = assigning_image
        node%I(ii, ni_P_is_boundary) = 0

        !% This bit of code checks the current node image, and compares it to the images of the adjacent links
        current_node_image = node%I(ii, ni_P_image)
        !adjacent_links = node%I(ii, ni_Mlink_u1:ni_Mlink_d3)
        adjacent_links = node%I(ii, ni_MlinkStart:ni_MlinkEnd)        !% brh20211219
        do jj = 1, size(adjacent_links)
            if ( adjacent_links(jj) == nullValueI ) then
                cycle
            end if
            adjacent_link_image = link%I(adjacent_links(jj), li_P_image)
            !% If the adjacent link and current node are on different images, then that node is a boundary
            if ( adjacent_link_image /= current_node_image ) then
                node%I(ii, ni_P_is_boundary) = node%I(ii, ni_P_is_boundary) + 1
            end if
        end do
    end do

end subroutine init_partitioning_default
!
!==========================================================================
!==========================================================================
!
subroutine init_partitioning_random()
    ! ----------------------------------------------------------------------------------------------------------------
    !
    ! Description:
    !   The random partitioning algorithm populates the partitioning columns of the link%I-node%I arrays.  An alternative
    !   to the default partitioning algorithm, the random partitioning algorithm looks at each link and node and assigns
    !   it to a random image (after checking to ensure that image is not full).
    !
    ! -----------------------------------------------------------------------------------------------------------------
    integer :: ii, jj, N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2, N_nJ1
    integer :: total_num_elements, num_attributed_elements, assigning_image
    integer :: current_node_image, adjacent_link_image
    real(8) :: partition_threshold, rand_num
    !% ----------------------------------------------------------------------------------------------------------------
    if (crashYN) return
    !% Determines the number of nodes of each type for the purpose of calculating partition threshold
    call util_count_node_types(N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2, N_nJ1)

    !% HACK The total number of elements is the sum of the elements from the links, plus the number of each node_type
    !% multiplied by how many elements are expected for that node_type
    total_num_elements = sum(link%I(:, li_N_element))        &
                            + (N_nBCup * N_elem_nBCup)       &
                            + (N_nBCdn * N_elem_nBCdn)       &
                            + (N_nJm * N_elem_nJm)           &
                            + (N_nStorage * N_elem_nStorage) &
                            + (N_nJ2 * N_elem_nJ2)           &
                            + (N_nJ1 * N_elem_nJ1)                            !% brh20211217
    partition_threshold = ( total_num_elements / real(num_images()) )

    !% Initialize the arrays that will hold the number of elements already on an image (and whether that image is full)
    elem_per_image(:) = 0
    image_full(:) = .false.

    !% This loop counts the elements attributed to each link, and assigns the link to an image
    do ii = 1, size(link%I, 1)

        !% Calculates a random number and maps it onto an image number
        call random_number(rand_num)
        assigning_image = int(rand_num*num_images()) + 1

        !% If the image number selected is already full, pick a new number
        do while ( image_full(assigning_image) .eqv. .true. )
            call random_number(rand_num)
            assigning_image = int(rand_num*num_images()) + 1
        end do

        !% elem_per_image is incremented by li_N_element for the current link
        elem_per_image(assigning_image) = elem_per_image(assigning_image) + link%I(ii, li_N_element)

        !% If the number of elements is greater than the partition threshold, that image number is closed
        !% Note, this check after the assigning_image has been selected allows for images be over-filled
        if ( elem_per_image(assigning_image) > partition_threshold ) then
            image_full(assigning_image) = .true.
        end if

        !% Assign the link to the current image
        link%I(ii, li_P_image) = assigning_image
    end do

    !% This loop counts the elements attributed to each node, and assigns the node to an image
    do ii = 1, size(node%I, 1)

        !% Calculates a random number and maps it onto an image number
        call random_number(rand_num)
        assigning_image = int(rand_num*num_images()) + 1

        !% If the image number selected is already full, pick a new number
        do while ( image_full(assigning_image) .eqv. .true. )
            call random_number(rand_num)
            assigning_image = int(rand_num*num_images()) + 1
        end do

        !% elem_per_image is incremented by the number of elements associated with each node type
        ! if ( node%I(ii, ni_node_type) == nBCup ) then
        !     elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nBCup
        ! else if ( node%I(ii, ni_node_type) == nBCdn ) then
        !     elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nBCdn
        ! else if ( node%I(ii, ni_node_type) == nStorage ) then
        !     elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nStorage
        ! else if ( node%I(ii, ni_node_type) == nJ2 ) then
        !     elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nJ2
        ! else if ( node%I(ii, ni_node_type) == nJm ) then
        !     elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nJm
        ! end if

        !% brh20211217 -- revised to add nJ1
        select case (node%I(ii, ni_node_type))
        case (nBCup)
            elem_per_image(assigning_image) = elem_per_image(assigning_image)+ N_elem_nBCup
        case (nBCdn)
            elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nBCdn
        case (nStorage)
            elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nStorage
        case (nJ1)
            elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nJ1
        case (nJ2)
            elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nJ2
        case (nJM)
            elem_per_image(assigning_image) = elem_per_image(assigning_image) + N_elem_nJm
        case default 
            print *, 'CODE ERROR: unknown node type # of ',node%I(ii, ni_node_type)
            print *, 'which has key of ',trim(reverseKey(node%I(ii, ni_node_type)))
            !stop 
            call util_crashpoint(1098226)
            return
        end select

        !% If the number of elements is greater than the partition threshold, that image number is closed
        !% Note, this check after the assigning_image has been selected allows for images be over-filled
        if ( elem_per_image(assigning_image) > partition_threshold ) then
            image_full(assigning_image) = .true.
        end if

        !% Assigns the nodes to an image, initializes the is_boundary check to 0
        node%I(ii, ni_P_image) = assigning_image
        node%I(ii, ni_P_is_boundary) = 0

        !% This bit of code checks the current node image, and compares it to the images of the adjacent links
        current_node_image = node%I(ii, ni_P_image)
        !adjacent_links = node%I(ii, ni_Mlink_u1:ni_Mlink_d3)
        adjacent_links = node%I(ii, ni_MlinkStart:ni_MlinkEnd)        !% brh20211219
        do jj = 1, size(adjacent_links)
            if ( adjacent_links(jj) == nullValueI ) then
                cycle
            end if
            adjacent_link_image = link%I(adjacent_links(jj), li_P_image)
            !% If the adjacent link and current node are on different images, then that node is a boundary
            if ( adjacent_link_image /= current_node_image ) then
                node%I(ii, ni_P_is_boundary) = node%I(ii, ni_P_is_boundary) + 1
            end if
        end do
    end do

end subroutine init_partitioning_random
!
!==========================================================================
!==========================================================================
!
subroutine init_partitioning_linkbalance()
!-----------------------------------------------------------------------------
!
! Description:
!   a balanced partitioning algorithm which distriutes all the links equally
!   to the available number of processors
!
!-----------------------------------------------------------------------------
    integer :: ii, jj
    integer :: clink, clink_image, assigned_image
    integer :: start_id, end_id
    integer :: count, remainder, rank

    character(64) :: subroutine_name = 'init_partitioning_linkbalance'

!-----------------------------------------------------------------------------
    if (crashYN) return
    if (setting%Debug%File%partitioning) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    if (SWMM_N_link < num_images()) then
        call init_partitioning_default()
    else
        do rank = 0, num_images()-1
            count = SWMM_N_link / num_images()
            remainder = mod(SWMM_N_link, num_images())

            if (rank < remainder) then
                ! The first 'remainder' ranks get 'count + 1' tasks each
                start_id = rank * (count + 1)
                end_id = start_id + count
            else
                ! The remaining 'size - remainder' ranks get 'count' task each
                start_id = rank * count + remainder
                end_id = start_id + (count - 1)
            end if

            link%I(start_id+1:end_id+1, li_P_image) = rank+1
        end do
        node%I(:, ni_P_is_boundary) = 0
        do ii = 1, N_node
            assigned_image = nullvalueI
            do jj = 1,max_branch_per_node
                clink = node%I(ii, ni_idx_base1+jj)
                if (clink /= nullvalueI) then
                    clink_image = link%I(clink, li_P_image)
                    if ( (assigned_image /= nullValueI) .and. &
                        (assigned_image /= clink_image) ) then
                        node%I(ii, ni_P_is_boundary) = node%I(ii, ni_P_is_boundary) + 1
                    end if
                    if (clink_image < assigned_image) then
                        assigned_image = clink_image
                    end if
                end if
            end do
            node%I(ii, ni_P_image) = assigned_image
        end do
    end if
    if (setting%Debug%File%partitioning) then
        print *, link%I(:, li_P_image), "node%I(:, li_P_image)"
        print *, node%I(:, ni_P_image), "node%I(:, ni_P_image)"
        print *, node%I(:, ni_P_is_boundary), "node%I(:, ni_P_is_boundary)"
    end if

    if (setting%Debug%File%partitioning)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
end subroutine init_partitioning_linkbalance
!
!==========================================================================
!==========================================================================
!
function init_partitioning_metric_partsizebalance() result(part_size_balance)
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
    if (crashYN) return
    !% Reset the elem_per_image array to all zeros
    elem_per_image(:) = 0

    !% Iterate through the link%I array
    do ii = 1, size(link%I, 1)

        !% The current image is the one to which the current link has been assigned
        current_image = link%I(ii, li_P_image)

        !% Iterate the number of elements for the current image by li_N_element for that link
        elem_per_image(current_image) = elem_per_image(current_image) + link%I(ii, li_N_element)
    end do

    !% Iterate through the node%I array
    do ii = 1, size(node%I, 1)

        !% The current image is the one to which the current link has been assigned
        current_image = node%I(ii, ni_P_image)

        ! !% elem_per_image for the current image is incremented by the number of elements associated with each node type
        ! if ( node%I(ii, ni_node_type) == nBCup ) then
        !     elem_per_image(current_image) = elem_per_image(current_image) + N_elem_nBCup
        ! else if ( node%I(ii, ni_node_type) == nBCdn ) then
        !     elem_per_image(current_image) = elem_per_image(current_image) + N_elem_nBCdn
        ! else if ( node%I(ii, ni_node_type) == nStorage ) then
        !     elem_per_image(current_image) = elem_per_image(current_image) + N_elem_nStorage
        ! else if ( node%I(ii, ni_node_type) == nJ2 ) then
        !     elem_per_image(current_image) = elem_per_image(current_image) + N_elem_nJ2
        ! else if ( node%I(ii, ni_node_type) == nJm ) then
        !     elem_per_image(current_image) = elem_per_image(current_image) + N_elem_nJm
        ! end if

        !% brh20211217 -- revised to add nJ1
        select case (node%I(ii, ni_node_type))
        case (nBCup)
            elem_per_image(current_image) = elem_per_image(current_image) + N_elem_nBCup
        case (nBCdn)
            elem_per_image(current_image) = elem_per_image(current_image)+ N_elem_nBCdn
        case (nStorage)
            elem_per_image(current_image) = elem_per_image(current_image)+ N_elem_nStorage
        case (nJ1)
            elem_per_image(current_image) = elem_per_image(current_image) + N_elem_nJ1
        case (nJ2)
            elem_per_image(current_image) = elem_per_image(current_image) + N_elem_nJ2
        case (nJM)
            elem_per_image(current_image) = elem_per_image(current_image) + N_elem_nJm
        case default 
            print *, 'CODE ERROR: unknown node type # of ',node%I(ii, ni_node_type)
            print *, 'which has key of ',trim(reverseKey(node%I(ii, ni_node_type)))
            !stop 
            call util_crashpoint(73875)
            return
        end select
    end do

     

    !% The maximum and minimum number of elements for an image are determined
    max_elem = maxval(elem_per_image(:))
    min_elem = minval(elem_per_image(:))
    print*, "Elem_per_image", elem_per_image(:)

    !% The difference between the max and min number of elements per image is the part_size_balance objective metric
    part_size_balance = max_elem - min_elem

end function init_partitioning_metric_partsizebalance
!
!==========================================================================
!==========================================================================
!
function init_partitioning_metric_connectivity() result(connectivity)
    ! ----------------------------------------------------------------------------------------------------------------
    !
    ! Description:
    !   This function is used to calculate the connectivity metric for the given partition set (i.e. the output
    !   from any of the partitioning algorithms).  The connectivity metric is equal to the sum of the ni_is_boundary
    !   column of the node%I array after the partition has been completed.  Note: the ni_is_boundary column is given a
    !   value of 0 when the node is an internal partition node, and is incremented by 1 for each additional partition
    !   the node touches.
    !
    ! -----------------------------------------------------------------------------------------------------------------
    integer :: connectivity
    ! -----------------------------------------------------------------------------------------------------------------

    !% The sum of the ni_is_boundary column is the connectivity
    connectivity = sum(node%I(:, ni_P_is_boundary))
end function init_partitioning_metric_connectivity
!
!==========================================================================
!==========================================================================
!
subroutine init_timer_stop ()

    integer(kind=8) :: crate, cmax, cval
    
    if (this_image() == 1) then
        call system_clock(count=cval,count_rate=crate,count_max=cmax)
        setting%Time%WallClock%PartitionEnd = cval
    end if

end subroutine init_timer_stop
!
!==========================================================================
!   End module
!==========================================================================
!
end module partitioning
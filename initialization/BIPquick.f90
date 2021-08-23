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
    use discretization, only: init_discretization_nominal
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

contains

    ! -----------------------------------------------------------------------------------------------------------------

    subroutine init_partitioning_BIPquick()
    ! ----------------------------------------------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine serves as the BIPquick main.  It is the only public subroutine from the BIPquick module and its
    !   output is populated columns ni_P_image, ni_P_is_boundary in node%I and li_P_image in link%I.
    !
    ! -----------------------------------------------------------------------------------------------------------------
        character(64) :: subroutine_name = 'BIPquick_subroutine'

        integer       :: mp, ii, jj, image
        real(8)       :: partition_threshold, max_weight, phantom_node_start
        logical       :: ideal_exists = .false.
        integer       :: spanning_link = nullValueI
        integer       :: effective_root
        integer       :: phantom_node_idx, phantom_link_idx
        integer       :: connectivity

        real(8) :: start, intermediate, finish
        call cpu_time(start)
        ! -----------------------------------------------------------------------------------------------------------------
        if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),subroutine_name

        !% One processor bypass for BIPquick
        if ( num_images() == 1 ) then
            node%I(:, ni_P_image) = oneI
            node%I(:, ni_P_is_boundary) = zeroI
            link%I(:, li_P_image) = oneI
            return
        end if

        !% Initialize the temporary arrays needed for BIPquick
        call bip_initialize_arrays()

        !% This subroutine populates B_nodeI with the upstream neighbors for a given node
        call bip_network_processing()

        !% Determine what the phantom nodes will be named
        phantom_node_idx = N_node + 1
        phantom_link_idx = N_link + 1

        !% This subroutine populates the directweight column of B_nodeR
        call calc_directweight()

        !% BIPquick sweeps through the network a finite number of times
        do mp = 1, num_images()

            !% Save the current processor as image (used as input to trav_subnetwork)
            image = mp
            ! print*, "Partition", mp

            !% Reset the node totalweight column, the ideal_exists boolean, and spanning_link integer
            B_nodeR(:, totalweight) = 0.0
            ideal_exists = .false.
            spanning_link = nullValueI

            !% This subroutine populates the totalweight column of B_nodeR and calculates max_weight
            call calc_totalweight(max_weight)

            ! do ii = 1, size(B_nodeR, 1)
            !     print*, node%I(ii, ni_idx), B_nodeR(ii, :)
            ! end do

            !% The partition_threshold is the current max_weight divided by the number of processors remaining (including current mp)
            partition_threshold = max_weight/real(num_images() - mp + 1, 8)

            !% This subroutine determines if there is an ideal partition possible and what the effective root is
            effective_root = calc_effective_root(ideal_exists, max_weight, partition_threshold)

            if (ideal_exists) then

                !% This subroutine traverses the subnetwork upstream of the effective root
                !% and populates the partition column of node%I
                call trav_subnetwork(effective_root, image)

            else

                !% This subroutine checks if the partition threshold is spanned by any link (Case 2)
                call calc_spanning_link(spanning_link, partition_threshold)

                !% While the spanning link doesn't exist (i.e. when the system is still Case 3)
                do while ( (spanning_link == nullValueI) .and. ( ideal_exists .eqv. .false. ) )

                    !% This subroutine houses the litany of steps that are required for a Case 3 partition
                    call trav_casethree(effective_root, spanning_link, image, partition_threshold, max_weight, ideal_exists)

                    !% HACK - I'm fairly sure that this do-loop will work for repeated instances of Case 3

                end do

                !% The outcomes of the Case 3 while loop are a Case 1
                if ( ideal_exists .eqv. .true. ) then

                !% In which case traverse the subnetwork from the effective_root
                call trav_subnetwork(effective_root, image)

                !% Or Case 2
                else if ( spanning_link /= nullValueI ) then

                !% In which case the distance along the spanning_link to the phantom node is calculated
                phantom_node_start = calc_phantom_node_loc(spanning_link, partition_threshold)

                !% This subroutine creates a phantom node/link and adds it to node%I/link%I
                call phantom_node_generator&
                (spanning_link, partition_threshold, phantom_node_start, phantom_node_idx, phantom_link_idx)

                !% This subroutine does the same thing as the previous call to trav_subnetwork()
                call trav_subnetwork(phantom_node_idx, image)

                phantom_node_idx = phantom_node_idx + 1
                phantom_link_idx = phantom_link_idx + 1

                else
                print*, "Something has gone wrong in BIPquick Case 3, there is no ideal exists or spanning link"
                stop

                end if

            endif

        end do

        !% This subroutine assigns network links to images on the basis of their endpoint nodes
        call trav_assign_link()

        !% This subroutine calculates the ni_P_is_boundary column of the node%I array
        call calc_is_boundary()

        connectivity = connectivity_metric()

        if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),subroutine_name
    end subroutine init_partitioning_BIPquick
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine bip_initialize_arrays()
    ! ----------------------------------------------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine allocates the temporary arrays that are needed for the BIPquick algorithm
    !   These temporary arrays are initialized in Globals, so they're not needed as arguments to BIPquick subroutines
    !
    ! -----------------------------------------------------------------------------------------------------------------
        character(64) :: subroutine_name = 'bip_allocate_arrays'
    ! -----------------------------------------------------------------------------------------------------------------
        if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),subroutine_name

        B_nodeI(:,:) = nullValueI
        B_nodeR(:,:) = zeroR
        totalweight_visited_nodes(:) = .false.
        partitioned_nodes(:) = .false.
        partitioned_links(:) = .false.
        weight_range(:,:) = zeroR
        accounted_for_links = .false.

        if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),subroutine_name
    end subroutine bip_initialize_arrays
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
    !   step eliminates the need to continuously jump back and forth between the node%I/link%I arrays during subsequent
    !   network traversals.
    !
    ! ----------------------------------------------------------------------------------------------------------------
        character(64) :: subroutine_name = 'bip_network_processing'

        integer :: upstream_link, upstream_node, uplink_counter
        integer ii, jj, uplinks
    ! ----------------------------------------------------------------------------------------------------------------
        if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),subroutine_name

        !% Iterate through the nodes array
        do ii= 1,size(node%I,1)

            if ( node%I(ii, ni_idx) == nullValueI ) then
                cycle
            end if

            !% The number of links upstream of a node
            uplink_counter = node%I(ii, ni_N_link_u)

            if (setting%Debug%File%BIPquick) then
                write (*,"(A,I2,A)") "Node " // node%Names(node%I(ii, ni_idx))%str &
                    // " has", node%I(ii, ni_N_link_u), " upstream links"
            end if

            !% Iterate through the links upstream of a node
            do uplinks= 1, uplink_counter
                !% If the link entry is not nullValueI (i.e. if it exists)
                if ( node%I(ii, ni_idx_base1 + uplinks) /= nullValueI ) then
                    !% Go to the upstream link to find the next upstream node
                    upstream_link = node%I(ii, ni_idx_base1 + uplinks)
                    upstream_node = link%I(upstream_link, li_Mnode_u)

                    !% Add the adjacent upstream node to B_nodeI
                    B_nodeI(ii, uplinks) = upstream_node
                endif
            enddo
        enddo

        if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),subroutine_name
    end subroutine bip_network_processing
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
        real(8)             :: weight, length, element_length
    ! --------------------------------------------------------------------------
        if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),function_name

        !% Sometimes the Interface gives garbage for these values so I need to adjust
        length = link%R(link_index, lr_Length)
        if ( (length < 0.0) .or. (length > nullValueI) ) then
        length = 1.0
        end if

        !% In particular sometimes the lr_ElementLength can be Infinity
        element_length = link%R(link_index, lr_ElementLength)
        if ( (element_length < 0.0) .or. (element_length > length) ) then
        element_length = 1.0
        end if

        !% The link weight is equal to the link length divided by the element length
        weight = length / element_length

        if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),function_name
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
        !   link weights using calc_link_weights(), and sums those weights to yield the
        !   B_nodeR(node, directweight)
        !
        !----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'calc_directweight'

        real(8)  :: lr_target
        integer :: rootnode_index, links_row, upstream_links
        integer :: ii, jj

        !--------------------------------------------------------------------------
        if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),subroutine_name

        !% Calculates directweight for each node
        do ii = 1, size(node%I,1)

        if ( node%I(ii, ni_idx) == nullValueI ) then
        cycle
        endif

        !% Need a loop bc multiple links might have a given node as its downstream endpoint
        do jj=1,size(link%I(:, li_Mnode_d))

        !% If the link has the current node as a downstream endpoint
        if (link%I(jj, li_Mnode_d) == node%I(ii, ni_idx)) then

        !% The directweight for that node is the running total of link weights
        B_nodeR(ii, directweight) = B_nodeR(ii, directweight) &
        + calc_link_weights(link%I(jj, li_idx))
        endif
        enddo
        enddo

        if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),subroutine_name
    end subroutine calc_directweight
    !
    !============================================================================
    !============================================================================
    !
    recursive subroutine calc_upstream_weight(weight_index, root)
        !-----------------------------------------------------------------------------
        !
        ! Description: Recursive subroutine that visits each node upstream of some root
        !  and adds the directweight to the root's totalweight.  This recursive subroutine
        !  is called for each node remaining in the network.
        !
        !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'calc_upstream_weight'

        integer :: upstream_node_list(3)
        integer, intent(in out) :: weight_index, root
        integer :: jj
        !--------------------------------------------------------------------------
        if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),subroutine_name

        !% If the node has not been visited this traversal (protective against cross connection bugs)
        !% and the node has not already been partitioned
        if ( (totalweight_visited_nodes(root) .eqv. .false.) .and. (partitioned_nodes(root) .eqv. .false.) ) then

        !% Mark the current root node as having been visited
        totalweight_visited_nodes(root) = .true.

        !% The totalweight of the weight_index node is increased by the root node's directweight
        B_nodeR(weight_index, totalweight) = B_nodeR(weight_index, totalweight) + B_nodeR(root, directweight)

        !% The adjacent upstream nodes are saved
        upstream_node_list = B_nodeI(root,:)

        !% Iterate through the adjacent upstream nodes
        do jj= 1, size(upstream_node_list)

        !% If the upstream node exists
        if( upstream_node_list(jj) /= nullValueI) then

        !% The call the recursive calc_upstream_weight on the weight_index node
        !% and the adjacent upstream node as the new root
        call calc_upstream_weight(weight_index, upstream_node_list(jj))
        endif
        enddo
        endif

        if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),subroutine_name
    end subroutine calc_upstream_weight
    !
    !============================================================================
    !============================================================================
    !
    subroutine calc_totalweight(max_weight)
        !-----------------------------------------------------------------------------
        !
        ! Description: This subroutine drives the calc_upstream_weight() recursive
        !  subroutine.  If a node remains in the network (i.e. hasn't been assigned to a
        !  partition yet), then it is passed as a root to the calc_upstream_weight().
        !
        !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'calc_totalweight'

        real(8), intent(in out) :: max_weight
        integer :: ii, weight_index

        !--------------------------------------------------------------------------
        if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),subroutine_name

        !% Calculates the totalweight for all nodes
        do ii=1, size(node%I,1)

        if ( node%I(ii, ni_idx) == nullValueI ) then
        cycle
        end if

        !% Provided that the node has not already been assigned to a partition
        if ( node%I(ii, ni_P_image) == nullValueI ) then

        !% The boolean for visited nodes during the upstream traversal is reset
        totalweight_visited_nodes(:) = .false.

        !% The weight_index is saved so that the recursive calc_upstream_weight knows where to add the totalweight updates
        weight_index = ii

        !% The weight_index (i.e. the current node) is passed to the first iteration
        !% of calc_upstream_weight as both the node being updated and the root node for traversal
        call calc_upstream_weight(weight_index, weight_index)
        endif
        enddo

        !% The max_weight is the largest totalweight value for this partition
        max_weight = (maxval(B_nodeR(:, totalweight)))

        if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),subroutine_name

    end subroutine calc_totalweight
    !
    !============================================================================
    !============================================================================
    !
    recursive subroutine trav_subnetwork(root, image)
        !-----------------------------------------------------------------------------
        !
        ! Description: This recursive subroutine visits every node upstream of a root node
        !  (where the root node is the "effective_root") and assigns that node to the current
        !  image.  It also updates the partitioned_nodes(:) array to "remove" that node from
        !  the network.
        !
        !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'trav_subnetwork'

        integer, intent(in out) :: root, image
        integer :: upstream_node_list(3)
        integer :: ii, jj, kk
        !--------------------------------------------------------------------------
        if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),subroutine_name

        !% If the root node has not been added to a partition
        if  ( partitioned_nodes(root) .eqv. .false. ) then

        !% Mark it as having been added to a partition
        partitioned_nodes(root) = .true.

        !% Add that node to the current image
        node%I(root, ni_P_image) = image

        !% Save the adjacent upstream nodes
        upstream_node_list = B_nodeI(root, :)

        !% Find the links that are in the subnetwork and mark them as being added to a partition
        do jj = 1, size(link%I, 1)
        if ( link%I(jj, li_Mnode_d) == root ) then
        partitioned_links(jj) = .true.
        end if
        end do

        !% Iterate through the upstream nodes
        do jj = 1, size(upstream_node_list)

        !% If the upstream node exists
        if ( upstream_node_list(jj) /= nullValueI ) then

        !% call the recursive subroutine on the new root node
        call trav_subnetwork(upstream_node_list(jj), image)
        endif
        enddo

        endif


        if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),subroutine_name
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

        integer :: potential_endpoints(size(node%I,1))
        integer :: endpoint_up, endpoint_dn, dn_image
        integer :: jj

        !--------------------------------------------------------------------------
        if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),subroutine_name

        !% Save the system nodes as potential endpoints
        potential_endpoints(:) = node%I(:, ni_idx)

        !% HACK - I'm not sure if this logic is air tight

        !% For each link, if the link is not nullValueI
        do jj=1, size(link%I,1)
        if ( link%I(jj, li_idx) /= nullValueI ) then

        !% Save the endpoints of that link
        endpoint_up = link%I(jj, li_Mnode_u)
        endpoint_dn = link%I(jj, li_Mnode_d)

        !% If the endpoints are in the nodes array and the link has not been accounted for
        if ( any(potential_endpoints(:) == endpoint_up) .and. &
        any(potential_endpoints(:) == endpoint_dn) .and. &
        ( accounted_for_links(jj) .eqv. .false.) ) then

        !% Assign the link to the downstream node's image
        dn_image = node%I(endpoint_dn, ni_P_image)
        link%I(jj, li_P_image) = dn_image

        !% Set the accounted for boolean to true
        accounted_for_links(jj) = .true.
        endif
        endif
        enddo

        !% This do loop just checks any links that somehow slipped through
        do jj=1, size(accounted_for_links, 1)
        if ( ( accounted_for_links(jj) .eqv. .false. ) &
        .and. ( link%I(jj, li_idx) /= nullValueI ) ) then
        endpoint_dn = link%I(jj, li_Mnode_d)
        dn_image = node%I(endpoint_dn, ni_P_image)
        link%I(jj, li_P_image) = dn_image
        accounted_for_links(jj) = .true.
        end if
        end do

        if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),subroutine_name
    end subroutine trav_assign_link
    !
    !============================================================================
    !============================================================================
    !
    function calc_effective_root(ideal_exists, max_weight, partition_threshold) result (effective_root)
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

        real(8), intent(in) :: max_weight, partition_threshold
        logical, intent(in out) :: ideal_exists
        real(8) :: nearest_overestimate
        integer :: ii
        !--------------------------------------------------------------------------
        if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),subroutine_name

        !% The nearest overestimate is set above the max_weight as a buffer
        nearest_overestimate = max_weight*1.1

        !% The effective_root is initialized as a nullValueI
        effective_root = nullValueI

        !% Searching through each node
        do ii=1, size(node%I,1)

        !% If the node has already been partitioned then go to the next one
        if (partitioned_nodes(ii) .eqv. .true. ) then
        cycle
        end if

        !% If the node's totalweight matches the partition_threshold to within a tolerance
        if ( abs ((B_nodeR(ii, totalweight) - partition_threshold)/partition_threshold) &
        < precision_matching_tolerance )  then

        !% Then the effective root is set and the ideal (Case 1) boolean is set to true
        effective_root = node%I(ii, ni_idx)
        ideal_exists = .true.
        exit
        endif

        !% Alternatively, if the totalweight is greater than the partition threshold and
        !% less than the nearest overestimate
        if (&
        (B_nodeR(ii, totalweight) > partition_threshold) .and. &
        (B_nodeR(ii, totalweight) < nearest_overestimate) &
        ) then

        !% Then update the nearest overestimate and set the effective root
        nearest_overestimate = B_nodeR(ii, totalweight)
        effective_root = node%I(ii, ni_idx)
        endif
        enddo

        !% The effective root is the one that most nearly overestimates the partition threshold
        !% reverse The Price Is Right style
        ! print*, "Effective root:", effective_root

        if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),subroutine_name
    end function calc_effective_root
    !
    !============================================================================
    !============================================================================
    !
    subroutine calc_spanning_link(spanning_link, partition_threshold)
        !-----------------------------------------------------------------------------
        !
        ! Description: This subroutine is used to search for a link that spans the
        !  partition_threshold (Case 2).  A link is considered to "span" if the
        !  partition_threshold is in the range (totalweight_upstream_node,
        !  totalweight_upstream_node + weight_of_link).
        !
        !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'calc_spanning_link'

        integer, intent(in out) :: spanning_link
        real(8), intent(in)     :: partition_threshold
        integer :: weight_range(2)
        integer :: upstream_node
        integer :: ii, jj
        !--------------------------------------------------------------------------
        if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),subroutine_name

        !% Check each link for spanning the partition threshold
        do jj=1, size(link%I,1)

        if ( link%I(jj, li_idx) == nullValueI ) then
        cycle
        end if

        !% Save the upstream node of the current link
        upstream_node = link%I(jj, li_Mnode_u)

        !% The first entry of the weight_range is the upstream node's totalweight
        weight_range(oneI) = B_nodeR(upstream_node, totalweight)

        !% The second entry is the first entry + the weight of the link
        weight_range(twoI) = weight_range(oneI) + calc_link_weights(link%I(jj, li_idx))

        ! print*, "Link", jj, "Weight Range", weight_range(:), "Partition Threshold", partition_threshold, &
        !   "Partitioned Links", partitioned_links(jj)

        !% If the partition threshold is between the weight_range entries
        !% and that link has not yet been partitioned
        if( (weight_range(oneI) < partition_threshold) .and. &
        (partition_threshold < weight_range(twoI)) .and. &
        (partitioned_links(jj) .eqv. .false.) ) then

        !% The current link is the spanning link
        spanning_link = link%I(jj, li_idx)
        ! print*, "The new spanning link is:  ", spanning_link

        !% Mark this link as being partitioned
        partitioned_links(jj) = .true.

        !% Only need one spanning link - if found, exit
        exit
        endif
        enddo

        if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),subroutine_name
    end subroutine calc_spanning_link
    !
    !==========================================================================
    !==========================================================================
    !
    function calc_phantom_node_loc(spanning_link, partition_threshold) result(length_from_start)
        !-----------------------------------------------------------------------------
        !
        ! Description:  This function is used to calculate how far along the spanning_link
        !  (from the upstream node) the phantom_node should be placed.

        ! HACK - might want to consider snapping this to the nearest element subdivision
        !
        !-----------------------------------------------------------------------------
        character(64)   :: function_name = 'calc_phantom_node_loc'

        real(8), intent(in) :: partition_threshold
        integer, intent(in) :: spanning_link
        real(8) :: length_from_start, total_length, start_weight, weight_ratio, link_weight
        integer :: upstream_node
        !--------------------------------------------------------------------------
        if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),function_name

        !% The length of the spanning_link
        total_length = link%R(spanning_link, lr_Length)

        !% The weight of the spanning_link
        link_weight = calc_link_weights(spanning_link)

        !% The upstream node from the spanning_link
        upstream_node = link%I(spanning_link, li_Mnode_u)

        !% The totalweight of the upstream node
        start_weight = B_nodeR(upstream_node, totalweight)

        !% The weight_ratio yields the factor that the link would have to be
        !% to have a downstream weight equal to the partition_threshold
        weight_ratio = (partition_threshold - start_weight) / link_weight

        !% Multiply the total_length by the weight_ratio to get the distance from
        !% the upstream node that the phantom_node should be generated
        length_from_start = weight_ratio * total_length

        if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),function_name
    end function calc_phantom_node_loc
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine phantom_node_generator &
        (spanning_link, partition_threshold, phantom_node_start, phantom_node_idx, phantom_link_idx)
        !-----------------------------------------------------------------------------
        !
        ! Description: This subroutine populates the node%I/link%I/link%R arrays with the phantom
        !  node/link that has been generated by a Case 2 system.  The node%I entries
        !  are set from the properties of phantom_nodes, the link%I/R phantom entries are copied
        !  from the link%I/R real entries then updated.
        !
        !-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'phantom_node_generator'

        integer, intent(in out)   :: phantom_node_idx, phantom_link_idx
        real(8), intent(in)       :: phantom_node_start, partition_threshold
        integer, intent(in)       :: spanning_link
        integer :: upstream_node_list(3)
        integer :: downstream_node, upstream_node
        integer :: kk
        real    :: l1, l2, y1, y2
        !--------------------------------------------------------------------------
        if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),subroutine_name

        !% The phantom node index is given
        node%I(phantom_node_idx, ni_idx) = phantom_node_idx

        !% The phantom node type is guaranteed to be a simple 2 link junction
        !% HACK - need to figure out if I make this more robust than the integer
        node%I(phantom_node_idx, ni_node_type) = nJ2

        !% The phantom node is guaranteed to have one upstream and one downstream link
        node%I(phantom_node_idx, ni_N_link_u) = oneI
        node%I(phantom_node_idx, ni_N_link_d) = oneI

        !% Initialize the upstream and downstream links as null values
        node%I(phantom_node_idx, ni_Mlink_u1:ni_Mlink_d3) = nullValueI

        !% The upstream link for the phantom node is the spanning link
        node%I(phantom_node_idx, ni_Mlink_u1) = spanning_link

        !% The downstream link for the phantom node is the phantom link
        node%I(phantom_node_idx, ni_Mlink_d1) = phantom_link_idx

        !% Reset the phantom node directweight to 0.0 (for cleanliness, this won't matter)
        B_nodeR(phantom_node_idx, directweight) = 0.0

        !% By definition, the phantom node totalweight will be the partition_threshold
        B_nodeR(phantom_node_idx, totalweight) = partition_threshold

        !% Find the adjacent upstream nodes for the B_nodeI traversal array
        upstream_node = link%I(spanning_link, li_Mnode_u)
        upstream_node_list(:) = (/upstream_node, nullValueI, nullValueI/)
        B_nodeI(phantom_node_idx, :) = upstream_node_list

        !% Copy the link row entries from the spanning link to the phantom link
        link%I(phantom_link_idx, :) = link%I(spanning_link, :)
        link%R(phantom_link_idx, :) = link%R(spanning_link, :)

        !% The phantom link length is the spanning_link length - the phantom node location
        link%R(phantom_link_idx, lr_Length) = link%R(spanning_link, lr_Length) - phantom_node_start
        link%R(spanning_link, lr_Length) = phantom_node_start
        call init_discretization_nominal(phantom_link_idx)
        call init_discretization_nominal(spanning_link)

        !% Save the original downstream node for the spanning link
        downstream_node = link%I(spanning_link, li_Mnode_d)

        !% The downstream node for the spanning link is set as the phantom node
        link%I(spanning_link, li_Mnode_d) = phantom_node_idx

        !% Reduce the downstream node directweight by the spanning link's new length
        B_nodeR(downstream_node, directweight) = B_nodeR(downstream_node, directweight) &
        - calc_link_weights(spanning_link)

        y1 = node%R(upstream_node, nr_Zbottom)
        y2 = node%R(downstream_node, nr_Zbottom)
        l1 = phantom_node_start
        l2 = link%R(phantom_link_idx, lr_Length)
        !% Interpolate zBottom
        node%R(phantom_node_idx, nr_Zbottom) = y2 + l2*(y1 - y2)/(l1 + l2)
        !% Interpolate InitialDepth
        y1 = node%R(upstream_node, nr_InitialDepth)
        y2 = node%R(downstream_node, nr_InitialDepth)
        node%R(phantom_node_idx, nr_InitialDepth) = y2 + l2*(y1 - y2)/(l1 + l2)

        !% Checks the adjacent nodes that were originally upstream of the downstream node
        upstream_node_list(:) = B_nodeI(downstream_node, :)
        do kk = 1, size(upstream_node_list)

            !% If the adjacent upstream node is the upstream node from the spanning link
            if ( upstream_node_list(kk) == upstream_node ) then

            !% Then replace it with the phantom node in B_nodeI
            B_nodeI(downstream_node, kk) = phantom_node_idx

            !% Also replace the downstream node's upstream link with the phantom link
            node%I(downstream_node, ni_idx_base1 + kk) = phantom_link_idx

            end if
        end do

        !% The resets the phantom index to having the phantom link and phantom node (as upstream node)
        link%I(phantom_link_idx, li_idx) = phantom_link_idx
        link%I(phantom_link_idx, li_Mnode_u) = phantom_node_idx

        !% HACK - need to check with Saz/Gerardo in Network_Define to see if I missed anything here

        if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),subroutine_name
    end subroutine phantom_node_generator
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine trav_casethree(effective_root, spanning_link, image, partition_threshold, max_weight, ideal_exists)
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

        integer, intent(in out)   :: effective_root, spanning_link, image
        real(8), intent(in out)   :: partition_threshold, max_weight
        logical, intent(in out)   :: ideal_exists
        integer   :: upstream_node_list(3), upstream_node
        real(8)   :: upstream_link_length, upstream_weight, total_clipped_weight
        integer   :: jj

    !--------------------------------------------------------------------------
        if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),subroutine_name

        print*, "Case 3 Found, effective_root is: ", effective_root

        !% The upstream node list is used to chose a branch for removal
        !% The effective root is guaranteed to have at least one upstream node
        upstream_node_list(:) = B_nodeI(effective_root, :)
        upstream_node = upstream_node_list(oneI)

        !% Find the link who's upstream node is the upstream_node
        do jj=1,size(link%I,1)

        !% If the link index is nullValue then cycle
        if ( link%I(jj, li_idx) == nullValueI ) then
        cycle
        end if

        !% If the link has upstream_node as its upstream node
        if (link%I(jj, li_Mnode_u) == upstream_node) then

        !% The clipped weight is the upstream totalweight + the link weight
        total_clipped_weight = B_nodeR(upstream_node, totalweight) &
        + calc_link_weights(jj)
        ! print*, "The weight being clipped from the effective node is", total_clipped_weight

        !% Tags the link as being partitioned (removes it from the spanning link potential)
        partitioned_links(jj) = .true.

        !% Exit saves time and records jj for later use
        exit
        endif
        enddo

        !% Reduce the effective_root directweight by the link length
        B_nodeR(effective_root, directweight) = &
        B_nodeR(effective_root, directweight) - calc_link_weights(jj)

        !% Reduce the effective_root totalweight by the total_clipped_weight
        B_nodeR(effective_root, totalweight) = &
        B_nodeR(effective_root, totalweight) - total_clipped_weight

        if ( total_clipped_weight <= zeroR ) then
            print*, "BIPquick Case 3: Haven't removed any weight"
            stop
        end if

        !% Reduce the partition_threshold by the total_clipped_weight too
        partition_threshold = partition_threshold - total_clipped_weight

        print*, "Calling subnetwork carving on ", upstream_node

        !% Assigns the link to the current image and removes it from future assignment
        link%I(jj, li_P_image) = image
        accounted_for_links(jj) = .true.

        !% Assign the subnetwork induced on the upstream node to the current image
        call trav_subnetwork(upstream_node, image)

        !% Checks the remaining network for a spanning_link
        call calc_spanning_link(spanning_link, partition_threshold)

        !% Resets the effective root to reflect updated system
        effective_root = calc_effective_root(ideal_exists, max_weight, partition_threshold)

        if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),subroutine_name
    end subroutine trav_casethree
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine calc_is_boundary()
    !-----------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the number of nodes that exist as
    !   a boundary between 1 or more partitions.  If a node has adjacent links that
    !   have been assigned to other partitions than its own, the ni_P_is_boundary
    !   column is incremented.
    !
    !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'calc_is_boundary'

        integer    :: adjacent_links(6), link_image
        integer    :: ii, kk

    !--------------------------------------------------------------------------
        if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),subroutine_name

        !% Initialize the ni_P_is_boundary column to 0
        node%I(:, ni_P_is_boundary) = zeroI

        !% Check each node in the network
        do ii = 1, size(node%I, 1)

        !% Create a list of links that are adjacent to the node
        adjacent_links = node%I(ii, ni_Mlink_u1:ni_Mlink_d3)

        !% Iterate through that list
        do kk = 1, size(adjacent_links)

        !% If the adjacent link doesn't exist, skip it
        if ( adjacent_links(kk) == nullValueI ) then
        cycle
        end if

        !% Check the image that the link has been assigned to (from trav_assign_link)
        link_image = link%I(adjacent_links(kk), li_P_image)

        !% If the link and the image are on separate images, increment the ni_P_is_boundary
        if ( link_image /= node%I(ii, ni_P_image) ) then
        node%I(ii, ni_P_is_boundary) = node%I(ii, ni_P_is_boundary) + 1
        end if
        end do
        end do

        if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),subroutine_name

    end subroutine calc_is_boundary
    !
    !==========================================================================
    !==========================================================================
    !

    function connectivity_metric() result(connectivity)
    !-----------------------------------------------------------------------------
    !
    ! Description: This subroutine is used to calculate the number of nodes that belong
    !   to multiple partitions
    !
    !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'connectivity_metric'

        integer  :: connectivity, ii

    !--------------------------------------------------------------------------
        if (setting%Debug%File%BIPquick) print *, '*** enter ', this_image(),subroutine_name

        connectivity = 0

        !% The sum of the ni_P_is_boundary column is the connectivity value
        do ii = 1, size(node%I, 1)
        connectivity = connectivity + node%I(ii, ni_P_is_boundary)
        end do

        if (setting%Debug%File%BIPquick) print *, '*** leave ', this_image(),subroutine_name

    end function connectivity_metric
    !
    !==========================================================================
    !==========================================================================
    !
end module BIPquick
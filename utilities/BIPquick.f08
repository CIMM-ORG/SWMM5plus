!==========================================================================
!2019-11-11 ==> Contributed by Eddie Tiernan
 module BIPquick ! the module name that is referenced in the main.f08

! the modules that need to precede BIPquick
 use assign_index
 use globals
 use setting_definition, only: setting
 
! this designation just means that every variable has to be given an explicit type/size
 implicit none

! the module's subroutines are private by default, allows some generically named subroutines to be used here and elsewhere
 private

! the public subroutines are ones that can be called from any module that has $ use BIPquick $ at the top
 public :: BIPquick_subroutine, BIPquick_Optimal_Hardcode, BIPquick_YJunction_Hardcode
 
 real(8), parameter :: precision_matching_tolerance = 1.0D-5 ! a tolerance parameter for whether or not two real(8) numbers are equal
 
 integer, parameter :: B_nr_directweight_u = 1 ! the cumulative weight of the links directly upstream of a node
 integer, parameter :: B_nr_totalweight_u  = 2 ! the cumulative weight of all links upstream of a node

 integer, parameter :: B_ni_idx_Partition = 1 ! the node index number
 integer, parameter :: B_ni_Partition_No = 2 ! the Partition number to which that node index belongs
 integer, parameter :: B_ni_is_boundary = 3 ! a binary marker that is 1 when the node is shared between partitions in the link-node paradigm

 integer, parameter :: B_li_idx_Partition = 1 ! the link index number
 integer, parameter :: B_li_Partition_No = 2 ! the Partition number to which that link index belongs
 
 contains
 
!-------------------------------------------------------------------------- 

! These two subroutines just populate the B_nodeI, B_linkI arrays with hardcoded values that *would* be the output of BIPquick_subroutine (if it worked)
! OPTIMAL.inp is the name of the system input file this hardcode subroutine emulates
 subroutine BIPquick_Optimal_Hardcode() 
    ! integer, dimension(:,:), intent(in out)  :: P_nodeI
    ! integer, dimension(:,:), intent(in out)   :: P_linkI
    ! integer,  dimension(:,:), intent(in out)  :: linkI
    ! integer,  dimension(:,:), intent(in out)  :: nodeI
    integer :: ii

    P_nodeI(:, B_ni_Partition_No) = (/1, 1, 1, 2, 3, 3, 1, 1, 2, 2, 2, 2, 3, 3, 3, 1, 3/)
    P_linkI(:, B_li_Partition_No) = (/1, 1, 1, 2, 3, 3, 3, 2, 2, 2, 2, 2, 3, 3, 3, 1, 1, 1/)

    do ii = 1, size(nodeI, 1)
        P_nodeI(ii, B_ni_idx_Partition) = nodeI(ii, ni_idx)
    enddo

    do ii = 1, size(linkI, 1)
        P_linkI(ii, B_li_idx_Partition) = linkI(ii, ni_idx)
    enddo

end subroutine BIPquick_Optimal_Hardcode

! Y_Junction_NetworkDefineTest.inp is the name of the system input file this hardcode subroutine emulates
subroutine BIPquick_YJunction_Hardcode() 
    integer :: ii

    print*, "The YJunction Arrays have been initialized and are size", size(P_nodeI,1), size(P_linkI,1)

    !B_nodeI(:, B_ni_Partition_No) = (/1, 2, 1, 3/)
    !B_nodeI(:, B_ni_is_boundary) = (/0, 0, 1, 0/)
    !B_linkI(:, B_li_Partition_No) = (/1, 2, 3/)


    P_nodeI(:, B_ni_Partition_No) = (/1, 2, 1, 3/)
    P_nodeI(:, B_ni_is_boundary) = (/0, 0, 1, 0/)
    P_linkI(:, B_li_Partition_No) = (/1, 2, 3/)

        
    
    linkI(:,li_BQ_image) = (/1, 2, 3/)
    nodeI(:,ni_BQ_image) = (/1, 2, 1, 3/)
    nodeI(:,ni_BQ_edge) = (/0, 0, 1, 0/)

    print*, "The P_nodeI array looks like"

    do ii = 1, size(nodeI, 1)
        P_nodeI(ii, B_ni_idx_Partition) = nodeI(ii, ni_idx)
        print*, P_nodeI(ii, :)
    enddo

    print*, "The P_linkI array looks like"

    do ii = 1, size(linkI, 1)
        P_linkI(ii, B_li_idx_Partition) = linkI(ii, ni_idx)
        print*, P_linkI(ii, :)
    enddo


end subroutine BIPquick_YJunction_Hardcode

! This subroutine is the main BIPquick subroutine
! It uses the link-node arrays initialized in $ call initialize_linknode_arrays() $ in the main.f08
! It also uses B_nodeI, B_linkI arrays that are initialized in allocate_storage.f08
! Two dummy arrays, B_nodeI and B_nodeR are used to contain some of the BIPquick specific parameters
 subroutine BIPquick_subroutine(linkI, nodeI, linkR, nodeR)
     real(8)    :: lr_target_default = 1.0                         ! for the time being, the target length of an element is a hardcoded parameter
     integer :: n_rows_in_file_node, n_rows_in_file_link    ! counter for the number of rows in the node/link .csv files
     integer :: n_rows_excluding_header_node, n_rows_excluding_header_link  ! number of rows in the node/link .csv files excluding the header, used to determine the size of the arrays
     integer :: n_rows_plus_processors_node, n_rows_plus_processors_link    ! the NodeMatrix and LinkMatrix arrays are of the size of the .csv file plus the number of processors, for phantom nodes/links
     integer :: multiprocessors = 3                         ! for the OPTIMAL example, the number of processors is 3.  This is a project dependent parameter
     integer :: phantom_node_idx, phantom_link_idx
     integer :: spanning_node_upstream
     integer :: weight_index = -998877                      ! the node_idx that is called in the nr_totalweight_assigner() that keeps track of which node's totalweight is being updated
     real(8) :: partition_threshold                        ! the collective weight (i.e. length) being searched for
     real(8) :: max_weight = 0.0                           ! the cumulative weight of the graph
     integer :: effective_root = -998877                    ! a dummy node that either has the partition_threshold as it's weight or is a near overestimate of the partition_threshold
     integer :: effective_root2 = -998877
     logical :: ideal_exists = .false.                      ! boolean to check whether the partition_threshold is found exactly at a node
     integer :: spanning_link = -998877                     ! link_idk for the link that "spans" (or abstractly contains the partition_threshold along its length)
     real(8) :: start_point = 0.0                          ! distance from the upstream node that the phantom_node should be placed
     integer :: root                                        ! node_idx for whichever node is currently being updated
     real(8) :: upstream_link_length = 0.0                 ! used to determine the weight being added to the node in nr_totalweight_assigner()
     integer :: upstream_node = -998877                     ! dummy node_idx used to track which node is being jumped to next in nr_totalweight_assigner()
     real(8) :: upstream_weight = 0.0                      ! sum of the upstream_link_lengths that are feeding into the current root node
     real(8) :: total_clipped_weight = 0.0                 ! amount of weight that has been removed due to a pseudo-partition.  this only occurs in case 3 
     integer :: phantom_array_location                      ! the row number in the nodeMatrix that will be updated with the phantom_node information
     real(8) :: local_max_weight = 0.0

     real(8), dimension(:,:), allocatable :: nodeMatrix        ! the nodeMatrix is the array from the nodes.csv file that is being reorganized by BIPquick
     real(8), dimension(:,:), allocatable :: linkMatrix        ! the linkMatrix is the array from the links.csv file that is being reorganized by BIPquick
     integer,  dimension(:,:),              intent(in out)  :: linkI
     integer,  dimension(:,:),              intent(in out)  :: nodeI
     integer,  dimension(:,:), allocatable  :: B_node_Partition
     integer,  dimension(:,:), allocatable  :: B_link_Partition
     real(8),  dimension(:,:),                 intent(in out)  :: linkR
     real(8),  dimension(:,:),                 intent(in out)  :: nodeR
     real(8),  dimension(:,:), allocatable  :: B_nodeR
     integer, dimension(:,:), allocatable   :: B_nodeI

     real(8), dimension(:,:), allocatable :: weight_range      ! this is an array of tuples that contains the max and min totalweights on each link, used to locate a spanning_link
     real(8), dimension(:,:), allocatable :: nodes_container   ! dummy array used in the BIPquick post-processing step 
     
     real(8), dimension(:,:,:), allocatable :: subnetwork_container_nodes      ! 3D array that stores the node info for each identified 'part', or partition component
     real(8), dimension(:,:,:), allocatable :: subnetwork_container_links      ! 3D array that stores the link info for each identified 'part'
     
     logical, allocatable, dimension(:) :: visited_flag_weight                  ! list of booleans that determines whether that node's directweight has already contributed to the weight_index node's totalweight
     logical, allocatable, dimension(:) :: visit_network_mask                   ! boolean list that tracks whether a node already belongs to a part
     logical, allocatable, dimension(:) :: partition_boolean                    ! boolean list that tracks whether a link already belongs to a part
     logical, allocatable, dimension(:) :: link_accounting
     
     integer, allocatable, dimension(:) :: accounted_for_links                  ! post-processing list of links for reordering the linkMatrix
     
     integer:: ii, jj, kk, mp                                   ! counters: ii - row in nodeMatrix, jj - row in linkMatrix, kk - secondary row counter for node/linkMatrix, mp - for each multiprocessor
     integer :: print_counter = 0
     integer :: sorted_connectivity_metric, unsorted_connectivity_metric, link_counter, missed_counter, while_counter, missing_links
     

     real(8) :: start, intermediate, finish
     call cpu_time(start)

     if ( setting%Partitioning%PartitioningMethod == P02 ) then

        ! Allocate and set the temporary arrays. (Multiprocessors - 1) represents the maximum number of phantom nodes
        ! B_node Partition will hold [ni_idx, Partition_No]
        allocate(B_node_Partition(size(nodeI,1) + multiprocessors - 1, 2))
        B_node_Partition(:,:) = nullValueI
        do ii = 1, size(nodeI, 1)
            B_node_Partition(ii, B_ni_idx_Partition) = nodeI(ii, ni_idx)
        enddo
        
        ! B_link Partition will hold [li_idx, Partition_No]
        allocate(B_link_Partition(size(linkI,1) + multiprocessors - 1, 2))
        B_link_Partition(:,:) = nullValueI
        do ii = 1, size(linkI, 1)
            B_link_Partition(ii, B_ni_idx_Partition) = linkI(ii, ni_idx)
        enddo

        ! B_nodeR will hold [B_nr_directweight_u, B_nr_totalweight_u]
        allocate(B_nodeR(size(nodeR, 1) + multiprocessors - 1, 2))
        B_nodeR(:,:) = nullValueR

        ! B_nodeI will hold [upstream_node1, 2, 3]
        allocate(B_nodeI(size(nodeI, 1) + multiprocessors - 1, max_us_branch_per_node))
        B_nodeI(:,:) = nullValueI

        do ii = 1, size(nodeI,1)
            print*, nodeI(ii,ni_idx), nodeI(ii, ni_node_type), nodeI(ii, ni_Mlink_u1:ni_Mlink_d3)
        enddo
        print*, '_______'
        do ii = 1, size(linkI,1)
            print*, linkI(ii,li_idx), linkI(ii, li_Mnode_u:li_Mnode_d), linkR(ii, lr_Length)
        enddo
        print*, '_______'

        call network_node_preprocessing(nodeI, linkI, B_nodeI)
        print*, "printing the upstream nodes"
        do ii = 1, size(B_nodeI,1)
            print*, B_nodeI(ii,:)
        enddo        
                  
         ! the idx of the phantom nodes are based on how many nodes exist, so this function determines the number of digits in the last node/link.  Then the phantom_index starts at 10^digits
         call phantom_naming_convention(nodeI, linkI, phantom_node_idx, phantom_link_idx)
                  
         ! the boolean lists need to be allocated and initialized as containing all .false. values
         allocate(visited_flag_weight(size(nodeI,1)))
         visited_flag_weight(:) = .false.
         
         allocate(visit_network_mask(size(nodeI,1)))
         visit_network_mask(:) = .false.
         
         allocate(partition_boolean(size(linkI,1)))
         partition_boolean(:) = .false.

         print*, "----------------------------------------" 

        !  determine the weight directly upstream of each node
         call local_node_weighting(nodeI, linkI, nodeR, linkR, B_nodeR, lr_target_default)

         do ii = 1, size(B_nodeR,1)
             print*, B_nodeR(ii, B_nr_directweight_u)
         enddo

         print*, "----------------------------------------"  
         
         ! allocate the tuple list to be the size of the linkMatrix
         allocate (weight_range(size(linkI,1),2))
          
         ! determine the weight directly upstream of each node
         local_max_weight = maxval(B_nodeR(:, B_nr_directweight_u))
         print*, "Local node max weight: ", local_max_weight

            
         ! this do loop is the BIPquick main process.  For each multiprocessor a Partition of nodes is populated.  The contents of this do loop determine which nodes are included.
         do mp = 1, multiprocessors
            print*, "Partition: ", mp
            ! for the identification of each part, the totalweight of all nodes must be reset and reassigned
            B_nodeR(:, B_nr_totalweight_u) = 0.0
            
            call nr_totalweight_assigner(B_nodeR, nodeR, linkR, B_nodeI, nodeI, linkI, B_node_Partition, & 
                B_link_Partition, weight_index, max_weight, visited_flag_weight, visit_network_mask)
            print*, "Finished totalweight assigner"
            print*, "The max weight is: ", max_weight
                        
            ! the partition_threshold for the next part is always equal to the size of the remaining graph divided by the number of processors remaining
            partition_threshold = max_weight/real(multiprocessors - mp + 1, 8)
            print*, "The partition threshold is: ", partition_threshold
                        
        	! the effective_root is determined either as the node that exactly equals the partition_threshold, or rather the node that is the nearest overestimate of the partition_threshold
            print*, "The effective root and totalweight are"
            effective_root = ideal_partition_check &
            (ideal_exists, max_weight, partition_threshold, B_nodeR, nodeI)
            stop
         enddo
            
!         	! within the ideal_partition_check, a boolean ideal_exists becomes .true. if the effective_root exactly equals the partition_threshold
!             if(ideal_exists .eqv. .true.) then
!         	   print*, "Case 1 Found"
!                 print*, "Print counter is ", print_counter

!         		! if an ideal node is found, the subnetwork_carving() function is called on that node to extract the upstream sub-graph as a complete part
!                 print*, nodeMatrix(effective_root+1, B_ni_idx)
!                 call subnetwork_carving &
!                     (effective_root + 1, mp, subnetwork_container_nodes, visit_network_mask, & 
!                      nodeMatrix, linkMatrix, print_counter)
!             else
!         		! if an ideal node is NOT found, the network needs to be checked for a link that "spans" the partition_threshold.  
!         		! "Span" means that the partition_threshold exists in the range [upstream_node_weight, downstream_node_weight]
!         		! The weight_range tuple array is initialized as null rather than zeroes because zero is a potentially real(8) value
!                 weight_range(:,:) = -998877
!         		! The spanning_check() function populates the weight_range array and then searches the array for any row (and 			! corresponding link) that bounds the partition_threshold
!                 call spanning_check &
!                     (spanning_link, spanning_node_upstream, weight_range, linkMatrix, nodeMatrix, lr_target, &
!                      partition_threshold, partition_boolean)
!                 print*, spanning_link, " is the spanning link"
!         		! If a spanning link is not found, that indicates no precise partition is attainable (Case 3)

!             	if ( spanning_link /= -998877 ) then
!             		print*, "Case 2 Found"
!                     print*, "The new spanning link is:  ", spanning_link
!             	endif 

!                 while_counter = 0
!                 do while (spanning_link == -998877)
!                     if ( effective_root2 /= effective_root ) then
!                         print*, "The effective root has changed ", effective_root2, effective_root
!                         while_counter = 0
!                     endif
!         	    	print*, "Case 3 Found, while counter is: ", while_counter
!                     subnetwork_container_nodes(mp, effective_root+1, :) = nodeMatrix(effective_root+1, :)
!                     upstream_node = nodeMatrix(effective_root+1, B_ni_u1_idx + while_counter)
!                     print*, upstream_node
!                     if ( upstream_node == -998877 ) then
!                         print*, "The upstream_node is empty but the spanning_link is still empty too"
!                         effective_root2 = ideal_partition_check &
!                         (ideal_exists, max_weight, partition_threshold, nodeMatrix)
!                         if ( effective_root2 == effective_root ) then
!                             print*, "The effective root is still ", effective_root, " but there are no upstream paths"
!                         endif
!                         stop
!                     endif
!                     do jj=1,size(linkMatrix,1)
!                         if (linkMatrix(jj, B_li_Mnode_u) == upstream_node) then
!                             upstream_link_length = linkMatrix(jj,B_lr_Length)
!                             upstream_weight = nodeMatrix(upstream_node+1, B_nr_totalweight_u)
!                             total_clipped_weight = upstream_weight + &
!                                         weighting_function(lr_target, upstream_link_length)
!                             print*, "The weight being clipped from the effective node is", total_clipped_weight
!                             partition_boolean(jj) = .true.
!                             !exit
!                         endif
!                     enddo
                        
!                     nodeMatrix(effective_root+1, B_nr_directweight_u) = &
!                         nodeMatrix(effective_root+1, B_nr_directweight_u) - weighting_function(lr_target, upstream_link_length)
!                     nodeMatrix(effective_root+1, B_nr_totalweight_u) = &
!                         nodeMatrix(effective_root+1, B_nr_totalweight_u) - total_clipped_weight
        		
!                     partition_threshold = partition_threshold - total_clipped_weight

!                     print*, "Calling subnetwork carving on ", upstream_node

!                     while_counter = while_counter + 1
                    
!                     call subnetwork_carving &
!                         (upstream_node + 1, mp, subnetwork_container_nodes, &
!                         visit_network_mask, nodeMatrix, linkMatrix, print_counter)


!                     call spanning_check &
!                         (spanning_link, spanning_node_upstream, weight_range, linkMatrix, nodeMatrix, &
!                         lr_target, partition_threshold, partition_boolean)

!                     effective_root2 = effective_root

!                     effective_root = ideal_partition_check &
!                         (ideal_exists, max_weight, partition_threshold, nodeMatrix)

!                     print*, "Print counter is ", print_counter
!                 enddo
                
!                 start_point = linear_interpolator(partition_threshold, &
!                 spanning_link, linkMatrix, weight_range, lr_target)
!                 call phantom_node_generator(spanning_link, start_point, mp, &
!                     partition_threshold, linkMatrix, nodeMatrix, &
!                     subnetwork_container_nodes, subnetwork_container_links, & 
!                     n_rows_excluding_header_node, n_rows_excluding_header_link, &
!                     phantom_array_location, phantom_node_idx, phantom_link_idx, lr_target, spanning_node_upstream)
                
!                 call subnetwork_carving &
!                     (int(nodeMatrix(phantom_array_location, B_ni_idx)) + 1, mp, &
!                     subnetwork_container_nodes, visit_network_mask, nodeMatrix, &
!                     linkMatrix, print_counter)
!             endif
            
!             do ii=1, size(visit_network_mask,1)
!                 if (visit_network_mask(ii) .eqv. .true.) then
!                     nodeMatrix(ii, B_nr_directweight_u) = 0.0
!                 endif
!             enddo
!         	print*, "End of a partition!"
!                 print*, "******************************************************"
!          enddo

!          call cpu_time(intermediate)
!          print*, "Time for partitioning - ", intermediate-start

!          !open(unit=29, file='Subcontainer2NetworkNode.csv', status='unknown')

!          !do ii=1, n_rows_excluding_header_node
!          !   write(29, '(*(I0 : ", "))') int(subnetwork_container_nodes(2, ii,:))
!          !enddo

!          ! Allocates and initializes two arrays that are used in postprocessing
!          allocate(accounted_for_links(size(linkMatrix,1)))
!          accounted_for_links(:) = -998877
!          link_counter = 1

!          allocate(link_accounting(size(linkMatrix,1)))
!          link_accounting(:) = .false.
         
!          allocate(nodes_container(size(nodeMatrix,1),size(nodeMatrix,2)))
!          nodes_container(:,:) = -998877
         
!          ! For each pre-defined part, populate a dummy variable with the nodes information for that part.
!          ! Call the subnetwork_links() function which determines, on the basis of a parts nodes, which links are also on that part
!          do mp = 1, size(subnetwork_container_nodes,1)
!         	print*, "Filling subnetwork links: ", mp
!             nodes_container(:,:) = subnetwork_container_nodes(mp,:,:)
!             call subnetworks_links (mp, nodes_container, subnetwork_container_links, &
!                 linkMatrix, accounted_for_links, link_counter)
!          enddo


!         ! print*, "Subcontainer Nodes"
!         ! do ii=1, size(subnetwork_container_nodes,1)
!         !	do jj=1, size(nodeMatrix,1)
!         !    		print*, subnetwork_container_nodes(ii,jj, B_ni_idx: 4)
!         !	enddo
!         !	print*, "-------------------------------------"
!         ! enddo
!         ! print*, "Subcontainer Links"
!         ! do ii=1, size(subnetwork_container_links,1)
!         !    do jj=1, size(linkMatrix,1)
!         !            print*, subnetwork_container_links(ii,jj, B_li_idx: 4)
!         !    enddo
!         !    print*, "-------------------------------------"
!         ! enddo

!          missing_links = check_links(linkMatrix, accounted_for_links) 
!          print*, "The number of links missing is: ", missing_links

!          sorted_connectivity_metric = sorted_metric(nodeMatrix, linkMatrix, subnetwork_container_links, multiprocessors)
!          print*, "The connectivity of the sorted matrix is: ", sorted_connectivity_metric
!          unsorted_connectivity_metric = unsorted_metric(nodeMatrix, linkMatrix, multiprocessors, &
!          partition_threshold, n_rows_excluding_header_node)
!          print*, "The connectivity of the unsorted matrix is: ", unsorted_connectivity_metric

!          call cpu_time(intermediate)
!          print*, "Time for processing links and metrics - ", intermediate-start
         
!          ! Calls a preprocessing function that collapses the separated graph parts back into a single array
!          call reorganize_arrays(nodeMatrix, linkMatrix, multiprocessors, &
!                     subnetwork_container_nodes, subnetwork_container_links)
         
!         ! open(unit=21, file='IdealNetworkNode.csv', status='unknown')
!         ! open(unit=22, file='IdealNetworkLink.csv', status='unknown')

!         ! do ii=1, n_rows_excluding_header_node
!         !    write(21, '(*(I0 : ", "))') int(nodeMatrix(ii,:))
!         ! enddo
         
!         ! do ii=1, n_rows_excluding_header_link
!         !    write(22, '(*(f10.5 : ", "))') linkMatrix(ii,:)
!         ! enddo

!          !print*, "Reorganized Nodes"
!          !do ii=1, size(nodeMatrix,1)
!          !   print*, nodeMatrix(ii,1:4)
!          !enddo
         
!          !print*, "Reorganized Links"
!          !do ii=1, size(linkMatrix,1)
!          !   print*, linkMatrix(ii,1:4)
!          !enddo

!          call cpu_time(finish)
!          print*, "Total Time - ", finish-start
     endif
end subroutine BIPquick_subroutine
!
!============================================================================ 
!============================================================================ 

!  function sorted_metric(nodeMatrix, linkMatrix, subnetwork_container_links, multiprocessors) result(sorted_connectivity_metric)
! !
! ! This subroutine determines how many points of connectivity exist in the
! ! BIPquick-sorted network
!  character(64) :: function_name = 'sorted_metric'
 
!  real(8), intent(in) :: nodeMatrix(:, :), linkMatrix(:,:), subnetwork_container_links(:, :, :)
!  integer, intent(in) :: multiprocessors
!  real(8), allocatable :: sorted_nodes(:, :)
!  integer :: upstream_node, downstream_node
!  integer :: ii, jj, mp
!  integer :: node_basis_counter = 0
!  integer :: node_contained_counter = 0
!  integer :: sorted_connectivity_metric
!  integer :: nullValue = -998877
! !-------------------------------------------------------------------------- 
!  if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name

!  allocate(sorted_nodes(multiprocessors, size(nodeMatrix,1)))
!  sorted_nodes(:,:) = nullValue
 
!  do ii = 1, size(nodeMatrix, 1)
!     if ( nodeMatrix(ii, B_ni_idx) /= nullValue ) then
!         node_basis_counter = node_basis_counter + 1
!     endif
!  enddo

!  do mp = 1, multiprocessors
!      do jj = 1, size(linkMatrix, 1)
!         if ( subnetwork_container_links(mp, jj, B_li_idx) /= nullValue) then
!             upstream_node = linkMatrix(jj, B_li_Mnode_u)
!             downstream_node = linkMatrix(jj, B_li_Mnode_d)
!             if ( all(sorted_nodes(mp, :) /= upstream_node) ) then
!                 sorted_nodes(mp, upstream_node + 1) = upstream_node
!             endif
!             if ( all(sorted_nodes(mp, :) /= downstream_node) ) then
!                 sorted_nodes(mp, downstream_node + 1) = downstream_node
!             endif
!         endif
!     enddo
!  enddo

!  do ii = 1, size(sorted_nodes,1)
!     if ( any(sorted_nodes(2,:) == sorted_nodes(1,ii)) .eqv. .true. ) then
!         print*, sorted_nodes(1,ii), " is the node_idx of the boundary"
!     endif
!  enddo


!  do mp = 1, multiprocessors
!     do ii = 1, size(nodeMatrix, 1)
!         if ( sorted_nodes(mp, ii) /= nullValue ) then
!             node_contained_counter = node_contained_counter + 1
!         endif
!     enddo
!  enddo

!  sorted_connectivity_metric = node_contained_counter - node_basis_counter
 
!  if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
!  end function sorted_metric
! !
! !============================================================================ 
! !============================================================================ 
! ! 
!  function unsorted_metric(nodeMatrix, linkMatrix, multiprocessors, partition_threshold, n_rows_excluding_header_node) & 
!     result(unsorted_connectivity_metric)
! !
! ! This subroutine determines how many points of connectivity exist in the
! ! BIPquick-sorted network
!  character(64) :: function_name = 'sorted_metric'
 
!  real(8), intent(in) :: nodeMatrix(:, :), linkMatrix(:, :)
!  real(8), intent(in) :: partition_threshold
!  real(8), allocatable :: default_nodes(:, :, :), default_links(:, :, :)
!  integer, intent(in) :: multiprocessors, n_rows_excluding_header_node
!  integer :: ii, jj, mp
!  integer :: node_basis_counter = 0
!  real(8) :: running_length = 0.0
!  integer :: upstream_node, downstream_node
!  integer :: node_contained_counter = 0
!  integer :: unsorted_connectivity_metric
!  integer :: nullValue = -998877
! !-------------------------------------------------------------------------- 
!  if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name

!  allocate(default_nodes(multiprocessors, size(nodeMatrix,1), size(nodeMatrix,2)))
!  allocate(default_links(multiprocessors, size(linkMatrix,1), size(linkMatrix,2)))

!  default_nodes(:,:,:) = nullValue
!  default_links(:,:,:) = nullValue
 
!  do ii = 1, size(nodeMatrix, 1)
!     if ( nodeMatrix(ii, B_ni_idx) /= nullValue ) then
!         node_basis_counter = node_basis_counter + 1
!     endif
!  enddo

!  mp = 1
!  do jj = 1, size(linkMatrix, 1)
!     if ( linkMatrix(jj, B_li_idx) /= nullValue) then
!         running_length = running_length + linkMatrix(jj, B_lr_Length)
!         default_links(mp, jj, B_li_idx) = linkMatrix(jj, B_li_idx)
!         upstream_node = linkMatrix(jj, B_li_Mnode_u)
!         downstream_node = linkMatrix(jj, B_li_Mnode_d)
!         if ( all(default_nodes(mp, :, B_ni_idx) /= upstream_node) ) then
!             default_nodes(mp, upstream_node + 1, :) = nodeMatrix(upstream_node + 1, :)
!         endif
!         if ( all(default_nodes(mp, :, B_ni_idx) /= downstream_node) ) then
!             default_nodes(mp, downstream_node + 1, :) = nodeMatrix(downstream_node + 1, :)
!         endif
!         if ( running_length >= partition_threshold ) then
!             running_length = 0
!             mp = mp + 1
!         endif
!     endif
!  enddo

!  open(unit=101, file='DefaultNodes_1.csv', status='unknown')
!  open(unit=102, file='DefaultNodes_2.csv', status='unknown')
!  open(unit=103, file='DefaultNodes_3.csv', status='unknown')


!  do ii=1, n_rows_excluding_header_node
!     write(101, '(*(I0 : ", "))') int(default_nodes(1, ii,:))
!  enddo
 
!  do ii=1, n_rows_excluding_header_node
!     write(102, '(*(I0 : ", "))') int(default_nodes(2, ii,:))
!  enddo

!  do ii=1, n_rows_excluding_header_node
!     write(103, '(*(I0 : ", "))') int(default_nodes(3, ii,:))
!  enddo

!  do mp = 1, multiprocessors
!     do ii = 1, size(nodeMatrix, 1)
!         if ( default_nodes(mp, ii, B_ni_idx) /= nullValue ) then
!             node_contained_counter = node_contained_counter + 1
!         endif
!     enddo
!  enddo

!  unsorted_connectivity_metric = node_contained_counter - node_basis_counter
 
!  if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
!  end function unsorted_metric
! !
! !============================================================================ 
! !============================================================================ 
! ! 
 function weighting_function(lr_target, link_length) result(weight)
!
! the weight attributed to each link (that will ultimately be assigned to the 
! downstream node) are normalized by lr_Target.  This gives an estimate of 
! computational complexity. In the future lr_Target can be customized for each 
! link.
 character(64)   :: function_name = 'weighting_function'
 
 real(8), intent(in)  :: lr_target
 real(8), intent(in)  :: link_length
 real(8) :: weight
!-------------------------------------------------------------------------- 
 if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name
 
 weight = link_length/lr_target
 
 if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
 end function weighting_function
!
!============================================================================ 
!============================================================================ 
! 
 subroutine null_value_convert(array)
!
! this function is used to convert the null values that are defaulted in the 
! B_nr_directweight_u and B_nr_totalweight_u columns of the nodes array into float 
! zeros.
 character(64) :: subroutine_name = 'null_value_convert'
 
 real(8), intent(in out) :: array(:)
 integer :: ii
 real(8) :: nullValue = nullvalueR
!-------------------------------------------------------------------------- 
 if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
 
 where (array(:) == nullValue)
     array(:) = 0.0
 endwhere
 
 if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
 end subroutine null_value_convert
!
!============================================================================ 
!============================================================================ 
! 
 subroutine local_node_weighting(nodeI, linkI, nodeR, linkR, B_nodeR, lr_target_default)
!
! this function takes each node index (ni_idx) and finds that ni_idx in the 
! downstream node column of the links array (li_Mnode_d).  From this row in 
! links the link weight is grabbed and ascribed to the node-in-questions local 
! weight (nr_directweight_u).
 character(64) :: subroutine_name = 'local_node_weighting'
 
 real(8)  :: lr_target
 real(8), intent(in)  :: lr_target_default
 integer, intent(in) :: nodeI(:,:), linkI(:,:)
 real(8), intent(in out) :: nodeR(:,:), linkR(:,:), B_nodeR(:,:)
 integer :: rootnode_index, links_row
 integer ii, jj
 
!-------------------------------------------------------------------------- 
 if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
 
 call null_value_convert(B_nodeR(:,B_nr_directweight_u))
 call null_value_convert(B_nodeR(:,B_nr_totalweight_u))
 
 do ii= 1,size(nodeI,1) ! This weighting function occurs for each node
    rootnode_index = nodeI(ii, ni_idx) ! Assign to variable the node index
    links_row = 1 ! Initialize the links_row which points to the index of the link
    do jj=1,size(linkI(:, li_Mnode_d))
        if (linkI(jj, li_Mnode_d) == rootnode_index) then
            lr_target = linkR(jj, lr_ElementLength)
            if (lr_target == nullValueR .or. lr_target < 0) then
                lr_target = lr_target_default
            endif
            print*, lr_target, "*"
            B_nodeR(ii, B_nr_directweight_u) &
                = B_nodeR(ii, B_nr_directweight_u) &
                + weighting_function(lr_target, linkR(jj, lr_Length))
        endif
        links_row = links_row + 1
    enddo
 enddo
 
 if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
 end subroutine local_node_weighting
!
!============================================================================ 
!============================================================================ 
!  
 subroutine network_node_preprocessing(nodeI, linkI, B_nodeI)
!
! This subroutine will be called immediately following the population of the node/linkMatrices
! The upstream adjacent nodes will be added as additional columns for each node.  This will enhance code speedup,
! as it will no longer be necessary to jump back and forth between nodeMatrix and linkMatrix
 character(64) :: subroutine_name = 'network_node_preprocessing'
 
 integer, intent(in) :: nodeI(:,:), linkI(:,:)
 integer, intent(in out) :: B_nodeI(:,:)
 integer :: upstream_link, upstream_node, uplink_counter, root_node
 integer ii, jj, uplinks
 
!-------------------------------------------------------------------------- 
 if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
 
 do ii= 1,size(nodeI,1)
    if ( mod(ii, 1) == 0 ) then
        print*, "processing upstream nodes of ni_idx:", ii
    endif
    root_node = nodeI(ii, ni_idx)
    uplink_counter = nodeI(ii, ni_N_link_u)
    print*, nodeI(ii, ni_idx), "has", nodeI(ii, ni_N_link_u), "upstream links"
    do uplinks= 1, uplink_counter
        if ( nodeI(ii, ni_idx_base1 + uplinks) /= nullValueI ) then
            upstream_link = nodeI(ii, ni_idx_base1 + uplinks)
            do jj= 1,size(linkI,1)
                if ( linkI(jj, li_idx) == upstream_link ) then
                    upstream_node = linkI(jj, li_Mnode_u)
                    B_nodeI(ii, uplinks) = upstream_node
                endif
            enddo
        endif  
    enddo
 enddo

 if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
 end subroutine network_node_preprocessing
!
!============================================================================ 
!============================================================================ 
! 
 recursive subroutine upstream_weight_calculation(weight_index, root, B_nodeR, nodeR, linkR, B_nodeI, nodeI, linkI, & 
                        B_node_Partition, B_link_Partition, visited_flag_weight, visit_network_mask)

 character(64) :: subroutine_name = 'upstream_weight_calculation'
 
 real(8), intent(in out) :: B_nodeR(:,:), nodeR(:,:), linkR(:,:)
 integer, intent(in out) :: B_nodeI(:,:), nodeI(:,:), linkI(:,:), B_node_Partition(:,:), B_link_Partition(:,:)
 integer :: upstream_node_list(3)
 integer, intent(in) :: weight_index
 integer :: root, node_row_contents, link_idx, new_root
 integer :: link_row_contents, node_upstream
 logical, intent(in out) :: visited_flag_weight(:)
 logical, intent(in) :: visit_network_mask(:)
 integer :: ii, jj, kk
 
!-------------------------------------------------------------------------- 
 if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
	 if ( (visited_flag_weight(root) .eqv. .false.) .and. (visit_network_mask(root) .eqv. .false.) ) then
		visited_flag_weight(root) = .true.
		B_nodeR(weight_index, B_nr_totalweight_u) &
		    = B_nodeR(weight_index, B_nr_totalweight_u) &
		    + B_nodeR(root, B_nr_directweight_u)
        print*, nodeI(weight_index, ni_idx), "was increased by", nodeI(root, ni_idx)

		! upstream_node_list(:) = nodeI(root, B_ni_u1_idx:B_ni_u3_idx)
        upstream_node_list = B_nodeI(root,:)		
        do jj= 1, size(upstream_node_list)
		    if( upstream_node_list(jj) /= -998877) then
                print*, upstream_node_list
                call upstream_weight_calculation(weight_index, upstream_node_list(jj)+1, & 
                    B_nodeR, nodeR, linkR, B_nodeI, nodeI, linkI, B_node_Partition, & 
                    B_link_Partition, visited_flag_weight, visit_network_mask)
		        ! call upstream_weight_calculation(weight_index, upstream_node_list(jj)+1, &
		        ! nodeMatrix, linkMatrix, visited_flag_weight, &
		        ! visit_network_mask)
		    endif
		enddo  
    endif
 
 if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
 end subroutine upstream_weight_calculation
!
!============================================================================ 
!============================================================================ 
! 
 subroutine nr_totalweight_assigner(B_nodeR, nodeR, linkR, B_nodeI, nodeI, linkI, B_node_Partition, & 
                B_link_Partition, weight_index, max_weight, visited_flag_weight, visit_network_mask)

 character(64) :: subroutine_name = 'nr_totalweight_assigner'
 
 real(8), intent(in out) :: B_nodeR(:,:), nodeR(:,:), linkR(:,:)
 integer, intent(in out) :: B_nodeI(:,:), nodeI(:,:), linkI(:,:), B_node_Partition(:,:), B_link_Partition(:,:)
 integer, intent(in out) :: weight_index 
 logical, intent(in out) :: visited_flag_weight(:)
 logical, intent(in out) :: visit_network_mask(:)
 real(8), intent(in out) :: max_weight
 integer :: ii, root
 real(8) :: nullValue = nullValueI
 
!-------------------------------------------------------------------------- 
 if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
 
 do ii=1, size(B_node_Partition,1)
    if ( mod(ii, 1) == 0 )then
        print*, "Currently weighting node: ", ii
    endif

    if (B_node_Partition(ii, B_ni_idx_Partition) /= nullValue ) then
        visited_flag_weight(:) = .false.
        weight_index = ii
        print*, "The weight_index is", ii
        call upstream_weight_calculation(weight_index, & 
                ii, B_nodeR, nodeR, linkR, B_nodeI, nodeI, linkI, & 
                B_node_Partition, B_link_Partition, & 
                visited_flag_weight, visit_network_mask)
    endif
 enddo
 
 max_weight = (maxval(B_nodeR(:, B_nr_totalweight_u)))
 
 if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
 end subroutine nr_totalweight_assigner
! !
! !============================================================================ 
! !============================================================================ 
! ! 
!  recursive subroutine subnetwork_carving &
!     (root, proc, subnetwork_container_nodes, visit_network_mask, &
!      nodeMatrix, linkMatrix, print_counter)

!  character(64) :: subroutine_name = 'subnetwork_carving'
 
!  logical, intent(in out) :: visit_network_mask(:)
!  real(8), intent (in out) :: subnetwork_container_nodes (:,:,:)
!  real(8), intent(in) :: nodeMatrix(:,:), linkMatrix(:,:)
!  integer :: upstream_node_list(3)
!  integer :: node_row_contents, link_row_contents, new_root, node_upstream
!  integer, intent(in) :: root, proc
!  integer :: ii, jj, kk
!  integer, intent (in out) :: print_counter
 
! !-------------------------------------------------------------------------- 
!  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
 
!  if  ( visit_network_mask(root) .eqv. .false. ) then
!     visit_network_mask(root) = .true.
!     subnetwork_container_nodes(proc, root, :) = nodeMatrix(root, :)
!     print_counter = print_counter + 1

!     upstream_node_list(:) = nodeMatrix(root, B_ni_u1_idx:B_ni_u3_idx)
!     !print*, root, upstream_node_list(:)
!     do jj= 1, size(upstream_node_list)
!         if ( upstream_node_list(jj) /= -998877 ) then
!             call subnetwork_carving(upstream_node_list(jj)+1, proc, &
!             subnetwork_container_nodes, visit_network_mask, nodeMatrix, linkMatrix, print_counter)
!         endif
!     enddo

!  elseif( visit_network_mask(root) .eqv. .true.) then
!         subnetwork_container_nodes(proc, root, :) = nodeMatrix(root, :)
!         print_counter = print_counter + 1
!  endif
 
!  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
!  end subroutine subnetwork_carving
! !
! !============================================================================ 
! !============================================================================ 
! ! 
!  subroutine subnetworks_links &
!     (proc,nodes_container, subnetwork_container_links, linkMatrix, &
!      accounted_for_links, link_counter)

!  character(64) :: subroutine_name = 'subnetworks_links'
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
!  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
 
!  allocate (potential_endpoints(size(nodes_container,1)))
!  potential_endpoints(:) = nodes_container(:, B_ni_idx)
 
!  do ii=1, size(linkMatrix,1)
!     if ( linkMatrix(ii, B_li_idx) /= -998877 ) then
!         endpoint1 = linkMatrix(ii, B_li_Mnode_u)
!         endpoint2 = linkMatrix(ii, B_li_Mnode_d)
!         if ( any(potential_endpoints(:) == endpoint1) .and. &
!              any(potential_endpoints(:) == endpoint2) .and. &
!              all(accounted_for_links(:) /= linkMatrix(ii, B_li_idx)) ) then
!             subnetwork_container_links(proc, ii,:) = linkMatrix(ii,:)
!             accounted_for_links(link_counter) = linkMatrix(ii, B_li_idx)
!             link_counter = link_counter + 1
!         endif
!     endif
!  enddo
 
!  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
!  end subroutine subnetworks_links
! !
! !============================================================================ 
! !============================================================================ 
! ! 
 function ideal_partition_check &
    (ideal_exists, max_weight, partition_threshold, B_nodeR, nodeI) result (effective_root)

 character(64) :: subroutine_name = 'ideal_partition_check'
 real(8), intent(in) :: max_weight, partition_threshold
 logical, intent(in out) :: ideal_exists
 integer :: effective_root
 real(8) :: nearest_overestimate
 real(8), intent(in) :: B_nodeR(:,:)
 integer, intent(in) :: nodeI(:,:)
 integer :: ii
 integer :: nullValue = nullValueI
!-------------------------------------------------------------------------- 
 if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name

 nearest_overestimate = max_weight*1.1
 effective_root = nullValue

 do ii=1, size(nodeI,1)
    if ( abs ((B_nodeR(ii, B_nr_totalweight_u) - partition_threshold)/partition_threshold) &
            < precision_matching_tolerance )  then
        effective_root = nodeI(ii, ni_idx)
        ideal_exists = .true.
        exit
    endif
    if (&
        (B_nodeR(ii, B_nr_totalweight_u) > partition_threshold) .and. &
        (B_nodeR(ii, B_nr_totalweight_u) < nearest_overestimate) &
       ) then
       nearest_overestimate = B_nodeR(ii, B_nr_totalweight_u)
       effective_root = nodeI(ii, ni_idx)
       print*, effective_root, B_nodeR(ii, B_nr_totalweight_u)
    endif
 enddo
 
 if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
 end function ideal_partition_check
! !
! !============================================================================ 
! !============================================================================ 
! ! 
!  subroutine spanning_check &
!     (spanning_link, spanning_node_upstream, weight_range, linkMatrix, nodeMatrix, lr_target, &
!      partition_threshold, partition_boolean)
! !
! !
! !
!  character(64) :: subroutine_name = 'spanning_check'
 
!  integer, intent(in out) :: spanning_link
!  integer, intent(out) :: spanning_node_upstream
!  real(8), allocatable, intent(in out) :: weight_range(:,:)
!  real(8), intent(in) :: linkMatrix(:,:), nodeMatrix(:,:)
!  real(8), intent(in) :: lr_target, partition_threshold
!  logical, intent(in out) :: partition_boolean(:)
!  integer :: upstream_node
!  integer :: ii, jj
 
! !-------------------------------------------------------------------------- 
!  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
 
!  spanning_link = -998877

!  do jj=1, size(linkMatrix,1)
!     upstream_node = linkMatrix(jj, B_li_Mnode_u)
!     do ii=1, size(nodeMatrix,1)
!         if(nodeMatrix(ii, B_ni_idx) == upstream_node) then
!             weight_range(jj,1) = nodeMatrix(ii, B_nr_totalweight_u)
!         endif
!     enddo
!     weight_range (jj, 2) = weight_range (jj, 1) + &
!         weighting_function(lr_target, linkMatrix(jj, B_lr_Length))
!  enddo

!  do jj=1, size(weight_range,1)
!     if(&
!         (weight_range(jj,1) < partition_threshold) .and. &
!         (partition_threshold < weight_range(jj,2)) .and. &
!         (partition_boolean(jj) .eqv. .false.) &
!       ) then
!         spanning_link = linkMatrix(jj,B_li_idx)
!         print*, "The new spanning link is:  ", spanning_link
!         spanning_node_upstream = linkMatrix(jj, B_li_Mnode_u)
!         partition_boolean(jj) = .true.
!         GOTO 5568
!     endif
!  enddo
 
!  5568 if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name  
!  end subroutine spanning_check
! !
! !========================================================================== 
! !==========================================================================
! ! 
!  function linear_interpolator(partition_threshold, spanning_link, linkMatrix, &
!     weight_range, lr_target) result(length_from_start)

!  character(64)   :: function_name = 'linear_interpolator'
 
!  real(8), intent(in) :: linkMatrix(:,:), weight_range(:,:)
!  real(8), intent(in) :: partition_threshold, lr_target
!  integer, intent(in) :: spanning_link
!  real(8) :: length_from_start, total_length, start, weight_ratio
!  integer :: ii
! !-------------------------------------------------------------------------- 
!  if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name
!  do ii= 1,size(linkMatrix,1)
!     if (int(linkMatrix(ii, B_li_idx)) == spanning_link) then
!         total_length = linkMatrix(ii, B_lr_Length)
!         start = weight_range(ii,1)
!     endif
!  enddo
!  weight_ratio = (partition_threshold - start)&
!     /weighting_function(lr_target, total_length)
!  length_from_start = weight_ratio*total_length
!  if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
!  end function linear_interpolator
! !
! !========================================================================== 
! !==========================================================================
! ! 
!  subroutine phantom_naming_convention2(nodeMatrix, linkMatrix, phantom_node_idx, phantom_link_idx)

!  character(64)   :: function_name = 'phantom_naming_convention2'
 
!  real(8), intent(in) :: linkMatrix(:,:), nodeMatrix(:,:)
!  integer, intent(out) :: phantom_node_idx, phantom_link_idx
!  integer :: max_node_idx, max_link_idx
!  integer, allocatable :: node_indices(:), link_indices(:)
! !-------------------------------------------------------------------------- 
!  if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name
 
!  allocate(node_indices(size(nodeMatrix,1)))
!  allocate(link_indices(size(linkMatrix,1)))
 
!  node_indices(:) = int(nodeMatrix(:, B_ni_idx))
!  max_node_idx = maxval(node_indices(:))
 
!  link_indices(:) = int(linkMatrix(:, B_ni_idx))
!  max_link_idx = maxval(link_indices(:))
  
!  phantom_node_idx = max_node_idx
!  phantom_link_idx = max_link_idx
 
!  if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
!  end subroutine phantom_naming_convention2
! !
! !========================================================================== 
! !==========================================================================
! ! 
 subroutine phantom_naming_convention(nodeI, linkI, phantom_node_idx, phantom_link_idx)

 character(64)   :: function_name = 'phantom_naming_convention'

 integer, intent(in) :: linkI(:,:), nodeI(:,:)
 integer, intent(out) :: phantom_node_idx, phantom_link_idx
 integer :: max_node_idx, max_link_idx
 integer, allocatable :: node_indices(:), link_indices(:)
!-------------------------------------------------------------------------- 
 if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name

 allocate(node_indices(size(nodeI,1)))
 allocate(link_indices(size(linkI,1)))

 node_indices(:) = int(nodeI(:, ni_idx))
 max_node_idx = maxval(node_indices(:))

 link_indices(:) = int(linkI(:, li_idx))
 max_link_idx = maxval(link_indices(:))
    
 phantom_node_idx = max_node_idx + 1
 phantom_link_idx = max_link_idx + 1

 if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
 end subroutine phantom_naming_convention
! !
! !========================================================================== 
! !==========================================================================
! ! 
!  subroutine phantom_node_generator(spanning_link, start_point, mp, &
!     partition_threshold, linkMatrix, nodeMatrix, subnetwork_container_nodes, &
!     subnetwork_container_links, n_rows_excluding_header_node, &
!     n_rows_excluding_header_link, phantom_array_location, phantom_node_idx, phantom_link_idx, &
!     lr_target, spanning_node_upstream)
! !
! !
!  character(64) :: subroutine_name = 'phantom_node_generator'
 
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
!  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
 
!  phantom_name = int(phantom_node_idx) + phantom_counter
 
!  phantom_array_location = n_rows_excluding_header_node + mp
 
!  phantom_array_location_link = n_rows_excluding_header_link + mp
 
!  nodeMatrix(phantom_array_location, B_node_id) = phantom_name
 
!  nodeMatrix(phantom_array_location, B_ni_idx) = phantom_name
 
!  nodeMatrix(phantom_array_location, B_ni_node_type) = int(0)
 
!  nodeMatrix(phantom_array_location, B_ni_N_link_u) = int(1)

!  nodeMatrix(phantom_array_location, B_ni_N_link_d) = int(1)
 
!  nodeMatrix(phantom_array_location, B_ni_Mlink_u2:B_ni_Mlink_d3) = -998877
 
!  nodeMatrix(phantom_array_location, B_ni_Mlink_u1) = spanning_link
 
!  nodeMatrix(phantom_array_location, B_ni_Mlink_d1) = phantom_link_idx + phantom_counter
 
!  nodeMatrix(phantom_array_location, B_nr_directweight_u) = 0.0
 
!  nodeMatrix(phantom_array_location, B_nr_totalweight_u) = partition_threshold

!  nodeMatrix(phantom_array_location, B_ni_u1_idx) = spanning_node_upstream
 
!  do jj=1, size(linkMatrix,1)
!     if(linkMatrix(jj,B_li_idx) == spanning_link) then
!         linkMatrix(phantom_array_location_link,:) = linkMatrix(jj,:)
        
!         linkMatrix(phantom_array_location_link,B_lr_Length) &
!             = linkMatrix(jj,B_lr_Length) - start_point
            
!         downstream_node = linkMatrix(jj,B_li_Mnode_d)
	
!         do ii=1, size(nodeMatrix,1)
!             if (nodeMatrix(ii, B_ni_idx) == downstream_node) then
!                 nodeMatrix(ii,B_nr_directweight_u) = &
!                     nodeMatrix(ii,B_nr_directweight_u) - &
!                     weighting_function(lr_target, start_point)
! 				do kk=B_ni_u1_idx, B_ni_u3_idx
! 					if ( nodeMatrix(ii, kk) == spanning_node_upstream ) then
! 						nodeMatrix(ii, kk) = phantom_name
! 					endif
! 				enddo
!             endif
!         enddo
        
!         linkMatrix(jj, B_li_Mnode_d) = phantom_name
        
!         linkMatrix(jj,B_lr_Length) = start_point
!         exit
!     endif
!  enddo
 
!  linkMatrix(phantom_array_location_link, B_link_id) = phantom_link_idx + phantom_counter
 
!  linkMatrix(phantom_array_location_link, B_li_idx) = phantom_link_idx + phantom_counter
 
!  linkMatrix(phantom_array_location_link, B_li_Mnode_u) = phantom_name
 
!  phantom_counter = phantom_counter + 1 

!  !do ii=1, size(nodeMatrix, 1)
!  !   print*, nodeMatrix(ii,:)
!  !enddo
 
!  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name  
!  end subroutine phantom_node_generator
! !
! !============================================================================ 
! !============================================================================ 
! !
!  function number_of_lines_in_file(iunit) result(n_lines)
 
!  integer,intent(in)  :: iunit   ! the file unit number
!  integer             :: n_lines ! the number of lines in the file

!  character(len=1)    :: tmp
!  integer             :: istat
!  character(64) :: function_name = 'number_of_lines_in_file'
! !-------------------------------------------------------------------------- 
!  if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name
!  rewind(iunit)
!  n_lines = 0
!  do
!      read(iunit,fmt='(A1)',iostat=istat) tmp
!      if (is_iostat_end(istat)) exit
!      n_lines = n_lines + 1
!  end do
!  rewind(iunit)
!  if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
!  end function number_of_lines_in_file
! !
! !============================================================================ 
! !============================================================================ 
! ! 

! function check_links(linkMatrix, accounted_for_links) result(missing_links)
 
!  character(64) :: function_name = 'check_links'
!  real(8), intent(in) :: linkMatrix(:,:)
!  integer, intent(in) :: accounted_for_links(:)
!  integer :: ii, jj, linkMatrix_count, accounted_count, missing_links

! !-------------------------------------------------------------------------- 
!  if (setting%Debug%File%BIPquick) print *, '*** enter ',function_name

!  linkMatrix_count = 0
!  accounted_count = 0

!  do ii = 1, size(linkMatrix, 1)
!     if ( linkMatrix(ii, B_li_idx) /= -998877 ) then
!         linkMatrix_count = linkMatrix_count + 1
!     endif
!  enddo

!  do jj = 1, size(accounted_for_links, 1)
!     if ( accounted_for_links(jj) /= -998877 ) then
!         accounted_count = accounted_count + 1
!     endif
!  enddo

!  missing_links = linkMatrix_count - accounted_count

!  if (setting%Debug%File%BIPquick) print *, '*** leave ',function_name
!  end function check_links
! !
! !============================================================================ 
! !============================================================================ 
! ! 
!  subroutine reorganize_arrays(nodeMatrix, linkMatrix, multiprocessors, &
!             subnetwork_container_nodes, subnetwork_container_links)
! !
! ! this function is used to convert the null values that are defaulted in the 
! ! B_nr_directweight_u and B_nr_totalweight_u columns of the nodes array into float 
! ! zeros.
!  character(64) :: subroutine_name = 'reorganize_arrays'
!  real(8), intent (in) :: subnetwork_container_nodes (:,:,:)
!  real(8), intent (in) :: subnetwork_container_links (:,:,:)
!  real(8), intent(in out) :: nodeMatrix(:,:), linkMatrix(:,:)
!  integer, intent(in) :: multiprocessors
 
!  real(8), allocatable :: reorganizedNodes(:,:)
!  real(8), allocatable :: reorganizedLinks(:,:)
 
!  integer :: nodesRowCounter = 1, linksRowCounter = 1
!  integer :: ii, jj, mp
 
! !-------------------------------------------------------------------------- 
!  if (setting%Debug%File%BIPquick) print *, '*** enter ',subroutine_name
 
!  allocate(reorganizedNodes(size(nodeMatrix,1),size(nodeMatrix,2)))
!  reorganizedNodes = -998877
 
!  do mp=1, multiprocessors
!     do ii=1, size(subnetwork_container_nodes,2)
!         do jj=1, size(reorganizedNodes,1)
!             if (subnetwork_container_nodes(mp,ii,B_ni_idx) == reorganizedNodes(jj,B_ni_idx)) then
!                 GOTO 5566
!             endif
!         enddo
!         if (subnetwork_container_nodes(mp,ii,B_ni_idx) == -998877) then
!             GOTO 5566
!         endif
        
!         reorganizedNodes(nodesRowCounter,:) = subnetwork_container_nodes(mp,ii,:)
!         nodesRowCounter = nodesRowCounter + 1
!  5566   continue
!     enddo
!  enddo
 
!  nodeMatrix(:,:) = reorganizedNodes(:,:)
 
!  allocate(reorganizedLinks(size(linkMatrix,1),size(linkMatrix,2)))
!  reorganizedLinks = -998877
 
!  do mp=1, multiprocessors
!     do ii=1, size(subnetwork_container_links,2)
!         do jj=1, size(reorganizedLinks,1)
!             if (subnetwork_container_links(mp,ii,B_ni_idx) == reorganizedLinks(jj,B_ni_idx)) then
!                 GOTO 5567
!             endif
!         enddo
!         if (subnetwork_container_links(mp,ii,B_ni_idx) == -998877) then
!             GOTO 5567
!         endif
        
!         reorganizedLinks(linksRowCounter,:) = subnetwork_container_links(mp,ii,:)
!         linksRowCounter = linksRowCounter + 1
!  5567   continue
!     enddo
!  enddo
 
!  linkMatrix(:,:) = reorganizedLinks(:,:)
 
!  if (setting%Debug%File%BIPquick) print *, '*** leave ',subroutine_name
!  end subroutine reorganize_arrays
! !
! !========================================================================== 
! !==========================================================================
! !

 

 end module BIPquick
!==========================================================================

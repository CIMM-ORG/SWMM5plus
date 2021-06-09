module initialization
    use interface
    use partitioning
    use discretization
    use utility_allocate
    use utility_array
    use define_indexes
    use define_keys
    use define_globals
    use initial_condition
    use network_define
    use utility, only: util_export_linknode_csv
    use utility_array
    use define_settings, only: setting

    implicit none

!-----------------------------------------------------------------------------
!
! Description:
!    General initialization of data structures (not including network)
!
! Method:
!    Creates the arrays index structures that are used for accessing data.
!    Arguably, this could be done more simply but we want the fundamental
!    column indexes in array_index to be parameters rather than variables. By
!    using parameters we reduce the possibility of accidentally changing a
!    column definition.
!
! Note on naming:
!    The driver subroutine is named after the driver module (in this case,
!    initialization).  Subsequent subroutines are name such that the subroutine
!    name is essentially a path "init_<module>_<subroutine_name>"
!-----------------------------------------------------------------------------

    private

    public :: initialize_all

contains
    !
    !==========================================================================
    ! PUBLIC
    !==========================================================================
    !
    subroutine initialize_all()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   a public subroutine that calls all the private initialization subroutines
    !
    !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'initialize_all'

    !-----------------------------------------------------------------------------

        if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name

        !% def_load_settings is one of the few subroutines in the Definition modules
        call def_load_settings(setting%Paths%setting)
        if (this_image() == 1) then
            call execute_command_line("if [ -d debug ]; then rm -r debug; fi && mkdir debug")
        end if

        call init_read_arguments()

        if (setting%Verbose) print *, "Simulation Starts"

        call interface_init()

        call init_linknode_arrays()

        call init_partitioning()
        
        sync all 
        !% HACK: this sync call is probably not needed

        call init_network()

        call initial_condition_setup ()

        if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name
    end subroutine initialize_all
    !
    !==========================================================================
    ! PRIVATE
    !==========================================================================
    !
    subroutine init_linknode_arrays()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   Retrieves data from SWMM C interface and populates link and node tables
    !
    !-----------------------------------------------------------------------------

        integer       :: ii, total_n_links
        logical       :: l1, l2
        
        character(64) :: subroutine_name = 'init_linknode_arrays'

    !-----------------------------------------------------------------------------

        if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name

        if (.not. api_is_initialized) then
            print *, "ERROR: API is not initialized"
            stop
        end if

        ! Allocate storage for link & node tables
        call util_allocate_linknode()

        nodeI(:,ni_N_link_u) = 0
        nodeI(:,ni_N_link_d) = 0

        do ii = 1, N_link
            linkI(ii,li_idx) = ii
            linkI(ii,li_link_type) = interface_get_link_attribute(ii, link_type)
            linkI(ii,li_geometry) = interface_get_link_attribute(ii, link_geometry)
            linkI(ii,li_Mnode_u) = interface_get_link_attribute(ii, link_node1) + 1 ! node1 in C starts from 0
            linkI(ii,li_Mnode_d) = interface_get_link_attribute(ii, link_node2) + 1 ! node2 in C starts from 0

            ! HACK This is a temporary hardcode until Gerardo can populate this column from the CFL condition
            linkI(ii, li_N_element) = 10

            nodeI(linkI(ii,li_Mnode_d), ni_N_link_u) = nodeI(linkI(ii,li_Mnode_d), ni_N_link_u) + 1
            nodeI(linkI(ii,li_Mnode_d), ni_idx_base1 + nodeI(linkI(ii,li_Mnode_d), ni_N_link_u)) = ii
            nodeI(linkI(ii,li_Mnode_u), ni_N_link_d) = nodeI(linkI(ii,li_Mnode_u), ni_N_link_d) + 1
            nodeI(linkI(ii,li_Mnode_u), ni_idx_base2 + nodeI(linkI(ii,li_Mnode_u), ni_N_link_d)) = ii

            linkI(ii,li_InitialDepthType) = 1 ! TODO - get from params file
            linkR(ii,lr_Length) = interface_get_link_attribute(ii, conduit_length)

            ! linkR(ii,lr_TopWidth): defined in network_define.f08
            linkR(ii,lr_BreadthScale) = interface_get_link_attribute(ii, link_xsect_wMax)
            ! linkR(ii,lr_Slope): defined in network_define.f08
            linkR(ii,lr_LeftSlope) = interface_get_link_attribute(ii, link_left_slope)
            linkR(ii,lr_RightSlope) = interface_get_link_attribute(ii, link_right_slope)
            linkR(ii,lr_Roughness) = interface_get_link_attribute(ii, conduit_roughness)
            linkR(ii,lr_InitialFlowrate) = interface_get_link_attribute(ii, link_q0)
            linkR(ii,lr_InitialUpstreamDepth) = interface_get_node_attribute(linkI(ii,li_Mnode_u), node_initDepth)
            linkR(ii,lr_InitialDnstreamDepth) = interface_get_node_attribute(linkI(ii,li_Mnode_d), node_initDepth)
            linkR(ii,lr_InitialDepth) = (linkR(ii,lr_InitialDnstreamDepth) + linkR(ii,lr_InitialUpstreamDepth)) / 2.0
        end do

        do ii = 1, N_node
            total_n_links = nodeI(ii,ni_N_link_u) + nodeI(ii,ni_N_link_d)
            nodeI(ii, ni_idx) = ii
            if (interface_get_node_attribute(ii, node_type) == oneI) then ! OUTFALL
                nodeI(ii, ni_node_type) = nBCdn
            else if (total_n_links == twoI) then
                nodeI(ii, ni_node_type) = nJ2
            else if (total_n_links > twoI) then
                nodeI(ii, ni_node_type) = nJm
            end if
            ! Nodes with nBCup are defined in inflow.f08 -> (inflow_load_inflows)
            l1 = interface_get_node_attribute(ii, node_has_extInflow) == 1
            l2 = interface_get_node_attribute(ii, node_has_dwfInflow) == 1
            if (l1 .or. l2) then
                nodeYN(ii, nYN_has_inflow) = .true.
                ! if (total_n_links == 1) then
                !     nodeI(ii, ni_node_type) = nBCup
                ! end if
                if (nodeI(ii,ni_N_link_u) == zeroI) then
                    nodeI(ii, ni_node_type) = nBCup
                end if   
            end if

            nodeR(ii,nr_InitialDepth) = interface_get_node_attribute(ii, node_initDepth)
            nodeR(ii,nr_Zbottom) = interface_get_node_attribute(ii, node_invertElev)
        end do

        if (setting%Debug%File%initialization) then
            if (this_image() == 1) then
                call util_export_linknode_csv()
            end if
        end if

        if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name
    end subroutine init_linknode_arrays
    !
    !==========================================================================
    !==========================================================================
    !

    subroutine init_partitioning()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine calls the public subroutine from the utility module,
    !   partitioning.f08. It also calls a public subroutine from the temporary
    !   coarray_partition.f08 utility module that defines how big the coarrays
    !   must be.
    !
    !-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'init_partitioning'
    !-----------------------------------------------------------------------------

        if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name
        
        !% find the number of elements in a link based on nominal element length   
        call init_discretization_nominal()

        !% in order to keep the main() clean, move the following two subroutines here, BIPquick can be removed
        call init_partitioning_method()

        !% adjust the link lenghts by cutting off a certain portion for the junction branch
        !% this subroutine is called here to correctly estimate the number of elements and faces
        !% to allocate the coarrays. HACK: it can be moved someplace more suitable
        call init_discretization_adjustlinklength()

        !% calculate the largest number of elements and faces to allocate the coarrays
        call init_coarray_length()

        !% allocate elem and face coarrays
        call util_allocate_elemX_faceX()

        !% allocate colum idxs of elem and face arrays for pointer operation
        call util_allocate_columns()

        if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name

    end subroutine init_partitioning
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_coarray_length()
        ! for coarray length determination
        integer :: nimgs_assign
        integer, allocatable :: unique_imagenum(:)
        integer :: ii, jj, kk, idx, counter, elem_counter=0, face_counter=0, junction_counter=0

        integer :: duplicated_face_counter=0
        integer, allocatable :: node_index(:), link_index(:), temp_arr(:)
        character(64) :: subroutine_name = 'init_coarray_length'

        if (setting%Debug%File%utility_array) print *, '*** enter ',subroutine_name

        call util_image_number_calculation(nimgs_assign, unique_imagenum)

        allocate(N_elem(size(unique_imagenum,1)))
        allocate(N_face(size(unique_imagenum,1)))
        allocate(N_unique_face(size(unique_imagenum,1)))

        do ii=1, size(unique_imagenum,1)
            node_index = PACK([(counter, counter=1,size(nodeI,1))], nodeI(:, ni_P_image) .eq. unique_imagenum(ii))
            link_index = PACK([(counter, counter=1,size(linkI,1))], linkI(:, li_P_image) .eq. unique_imagenum(ii))
            ! create corresponding indices for node and link in this image

            ! The number of elements and faces is actually decided by the junctions
            ! So we will calculate the number of junction and base on different scenarios to decided
            ! how many elem/face are assigned to each image
            junction_counter = count(nodeI(node_index, ni_node_type) == nJm)

            !% first calculate the number of nodes in each partition, assign elems/faces for junctions
            elem_counter = elem_counter + J_elem_add * junction_counter
            face_counter = face_counter + J_face_add * junction_counter

            !% loop through the links and calculate the internal faces between elements
            do jj = 1, size(link_index,1)
                idx = link_index(jj)
                face_counter = face_counter + linkI(idx, li_N_element) - 1 !% internal faces between elems, e.g. 5 elements have 4 internal faces
                elem_counter = elem_counter + linkI(idx, li_N_element) ! number of elements
            enddo

            !% now we loop through the nodes and count the node faces
            do jj = 1, size(node_index,1)
                idx = node_index(jj)
                if (nodeI(idx, ni_node_type) .eq. nJ2) then
                    face_counter = face_counter + 1 !% add the face of 1-to-1 junction between 2 links
                elseif (nodeI(idx, ni_node_type) .eq. nBCup) then
                    face_counter = face_counter +1 !% add the upstream faces
                elseif (nodeI(idx, ni_node_type) .eq. nBCdn) then
                    face_counter = face_counter +1 !% add the downstream faces
                endif !% multiple junction faces already counted
            enddo

            !% Now we count the space for duplicated faces
            do jj = 1, size(link_index,1)
                idx = link_index(jj)
                !% check upstream node first
                if ( ( nodeI(linkI(idx, li_Mnode_u), ni_P_is_boundary) .eq. 1) .and. &
                    ( nodeI(linkI(idx, li_Mnode_u), ni_P_image) .ne. ii) ) then
                    face_counter = face_counter +1
                    duplicated_face_counter = duplicated_face_counter + 1
                endif
                ! then downstream node
                if ( ( nodeI(linkI(idx, li_Mnode_d), ni_P_is_boundary) .eq. 1) .and. &
                    ( nodeI(linkI(idx, li_Mnode_d), ni_P_image) .ne. ii) ) then
                    face_counter = face_counter +1
                    duplicated_face_counter = duplicated_face_counter + 1
                endif

            enddo

            N_elem(ii) = elem_counter
            N_face(ii) = face_counter
            N_unique_face(ii) = face_counter - duplicated_face_counter

            elem_counter = zeroI ! reset the counter
            face_counter = zeroI
            junction_counter = zeroI
            duplicated_face_counter = zeroI
        enddo

        max_caf_elem_N = maxval(N_elem)
        max_caf_face_N = maxval(N_face) ! assign the max value

        if (setting%Debug%File%utility_array) then
            do ii = 1, size(unique_imagenum,1)
                print*, 'Image => ', ii
                print*, 'Elements expected ', N_elem(ii)
                print*, 'Faces expected    ', N_face(ii)
            end do
        endif

        if (setting%Debug%File%utility_array)  print *, '*** leave ',subroutine_name

    end subroutine init_coarray_length
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine init_read_arguments()
        integer :: ii
        logical :: arg_param = .false.
        character(len=8) :: param
        character(len=256) :: arg

        do ii = 1, iargc()
            call getarg(ii, arg)
            if (.not. arg_param) then
                param = arg
                if (ii == 1) then
                    if (arg(:1) == '-') then
                        print *, "ERROR: it is necessary to define the path to the .inp file"
                        stop
                    end if
                    setting%Paths%inp = arg
                elseif ((trim(arg) == "-s") .or. & ! user provides settings file
                        ((trim(arg) == "--test") .or. (trim(arg) == "-t"))) then  ! hard coded test case
                    arg_param = .true.
                elseif ((trim(arg) == '--verbose') .or. (trim(arg) == "-v")) then
                    setting%Verbose = .true.
                elseif ((trim(arg) == '--warnings-off') .or. (trim(arg) == "-woff")) then
                    setting%Warning = .false.
                else
                    write(*, *) 'The argument ' // trim(arg) // ' is unsupported'
                    stop
                end if
            else
                arg_param = .false.
                if (trim(param) == '-s') then
                    setting%Paths%setting = arg
                elseif ((trim(arg) == "--test") .or. (trim(arg) == "-t")) then
                    setting%TestCase%UseTestCase = .true.
                    setting%TestCase%TestName = trim(arg)
                    if (trim(arg) == 'simple_channel') then
                    else if (trim(arg) == 'simple_orifice') then
                    else
                        write(*, *) 'The test case ' // trim(arg) // ' is unsupported. Please use one of the following:'
                        print *, new_line('')
                        print *, "simple_channel, simple_orifice, simple_pipe"
                        print *, "simple_weir, swashes, waller_creek"
                        print *, "y_channel, y_storage_channel"
                        stop
                    end if
                elseif (trim(param) == '--run-tests') then
                end if
            end if
        end do
    end subroutine init_read_arguments
    !
    !==========================================================================
    ! END OF MODULE
    !==========================================================================
    !
end module initialization

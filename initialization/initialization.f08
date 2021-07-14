module initialization
    use define_keys
    use define_globals
    use define_settings
    use define_indexes
    use interface
    use partitioning
    use discretization
    use utility_allocate
    use utility_array
    use initial_condition
    use network_define
    use utility, only: util_export_linknode_csv, util_create_warning_file
    use utility_array


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
        if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name
    !-----------------------------------------------------------------------------
        !% set the branchsign global -- this is used for junction branches (JB)
        !% for upstream (+1) and downstream (-1)
        !% HACK: for clarity and consistency, this probably should be moved into
        !% the init_network. Placed here for the time being in case we need it
        !% for translating link/node from SWMM-C or partitioning.
        branchsign(1:2:max_branch_per_node-1) = +oneR
        branchsign(2:2:max_branch_per_node)   = -oneR

        !% load the settings.json file with the default setting% model control structure
        !% def_load_settings is one of the few subroutines in the Definition modules
        call def_load_settings(setting%Paths%setting)

        !% execute the command line options provided when the code is run
        if (this_image() == 1) then
            call execute_command_line ("if [ -d debug ]; then rm -r debug; fi && mkdir debug")
        end if

        !% Create or replace warning files

        if(setting%Warning) call util_create_warning_file()
        
        
        !% read and store the command-line options
        call init_read_arguments ()

        if (setting%Verbose) print *, "Simulation Starts"

        !% initialize the API with the SWMM-C code
        call interface_init ()

        !% set up and store the SWMM-C link-node arrays in equivalent Fortran arrays
        call init_linknode_arrays ()

        !% partition the network for multi-processor parallel computation
        call init_partitioning ()

        call init_network ()

        call init_IC_setup ()

        call init_bc()

        !% wait for all the processors to reach this stage before starting the time loop
        sync all

        !% wait for all the processors to reach this stage before starting the time loop
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
    !   Retrieves data from EPA-SWMM interface and populates link and node tables
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

        node%I(:,ni_N_link_u) = 0
        node%I(:,ni_N_link_d) = 0

        do ii = 1, N_link
            link%I(ii,li_idx) = ii
            link%I(ii,li_link_type) = interface_get_link_attribute(ii, api_link_type)
            link%I(ii,li_geometry) = interface_get_link_attribute(ii, api_link_geometry)
            link%I(ii,li_Mnode_u) = interface_get_link_attribute(ii, api_link_node1) + 1 ! node1 in C starts from 0
            link%I(ii,li_Mnode_d) = interface_get_link_attribute(ii, api_link_node2) + 1 ! node2 in C starts from 0

            ! HACK This is a temporary hardcode until Gerardo can populate this column from the CFL condition
            link%I(ii, li_N_element) = 10

            node%I(link%I(ii,li_Mnode_d), ni_N_link_u) = node%I(link%I(ii,li_Mnode_d), ni_N_link_u) + 1
            node%I(link%I(ii,li_Mnode_d), ni_idx_base1 + node%I(link%I(ii,li_Mnode_d), ni_N_link_u)) = ii
            node%I(link%I(ii,li_Mnode_u), ni_N_link_d) = node%I(link%I(ii,li_Mnode_u), ni_N_link_d) + 1
            node%I(link%I(ii,li_Mnode_u), ni_idx_base2 + node%I(link%I(ii,li_Mnode_u), ni_N_link_d)) = ii

            !% HACK All links have the same initial depth type which is the default one
            !% a better approach would be to allow specific links to have specific depth
            !% types via an external JSON file for links whose path can be specified in
            !% setting%Link%PropertiesFile
            link%I(ii,li_InitialDepthType) = setting%Link%DefaultInitDepthType
            link%R(ii,lr_Length) = interface_get_link_attribute(ii, api_conduit_length)

            ! link%R(ii,lr_TopWidth): defined in network_define.f08
            link%R(ii,lr_BreadthScale) = interface_get_link_attribute(ii, api_link_xsect_wMax)
            ! link%R(ii,lr_Slope): defined in network_define.f08
            link%R(ii,lr_LeftSlope) = interface_get_link_attribute(ii, api_link_left_slope)
            link%R(ii,lr_RightSlope) = interface_get_link_attribute(ii, api_link_right_slope)
            link%R(ii,lr_Roughness) = interface_get_link_attribute(ii, api_conduit_roughness)
            link%R(ii,lr_InitialFlowrate) = interface_get_link_attribute(ii, api_link_q0)
            link%R(ii,lr_InitialUpstreamDepth) = interface_get_node_attribute(link%I(ii,li_Mnode_u), api_node_initDepth)
            link%R(ii,lr_InitialDnstreamDepth) = interface_get_node_attribute(link%I(ii,li_Mnode_d), api_node_initDepth)
            link%R(ii,lr_InitialDepth) = (link%R(ii,lr_InitialDnstreamDepth) + link%R(ii,lr_InitialUpstreamDepth)) / 2.0
        end do

        do ii = 1, N_node
            total_n_links = node%I(ii,ni_N_link_u) + node%I(ii,ni_N_link_d)
            node%I(ii, ni_idx) = ii
            if (interface_get_node_attribute(ii, api_node_type) == API_OUTFALL) then
                node%I(ii, ni_node_type) = nBCdn
            else if ((total_n_links == twoI)         .and. &
                     (node%I(ii,ni_N_link_u) == oneI) .and. &
                     (node%I(ii,ni_N_link_d) == oneI) )then
                node%I(ii, ni_node_type) = nJ2
            else if (total_n_links >= twoI) then
                node%I(ii, ni_node_type) = nJm
            end if
            l1 = interface_get_node_attribute(ii, api_node_has_extInflow) == 1
            l2 = interface_get_node_attribute(ii, api_node_has_dwfInflow) == 1
            if (l1 .or. l2) then
                node%YN(ii, nYN_has_inflow) = .true.
                if (node%I(ii,ni_N_link_u) == zeroI) then ! No upstream links
                    node%I(ii, ni_node_type) = nBCup
                end if
            end if

            node%R(ii,nr_InitialDepth) = interface_get_node_attribute(ii, api_node_initDepth)
            node%R(ii,nr_Zbottom) = interface_get_node_attribute(ii, api_node_invertElev)
        end do

        !% Update Link/Node names
        call interface_update_linknode_names()

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
    subroutine init_bc()
    !%-----------------------------------------------------------------------------
    !%
    !% Description:
    !%    Initializes boundary connditions
    !%
    !% Notes:
    !%    The structures are general enough to support 3 types of BCs:
    !%
    !%    BCup: updstream boundary condition which can be inflow or head BC
    !%    BCdn: downstream boundary condition which can be inflow or head BC
    !%    BClat: lateral inflow coming into and nJ2 or nJm node.
    !%
    !%    However, the code only supports inflow BCs for BCup and BClat,
    !%    and head BCs for BCdn, mimimcking EPA-SWMM 5.13 functionalities.
    !%    Further developments allowing other types of inflow and head BCs,
    !%    should store the respective BC in either the BC%inflowX or the
    !%    BC%headX arrays defining the corresponding type of BC (i.e., BCup,
    !%    BCdn, and BClat) in the BC%xI(:,bi_category) column.
    !%
    !%-----------------------------------------------------------------------------
        integer :: ii, nidx, ntype, counter_bc_er
        integer :: ntseries, nbasepat
        logical :: has_dwf_inflow, has_ext_inflow
        character(64) :: subroutine_name = "init_bc"
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%initialization)  print *, '*** enter ', subroutine_name

        call util_allocate_bc()

        !% Convention to denote that xR_timeseries arrays haven't been fetched
        if (N_flowBC > 0) then
            BC%flowI(:,bi_fetch) = 1
            BC%flowIdx(:) = 0
            !% Convention to denote association between nodes and face/elements
            !% BCup and BCdn BCs are associated with faces, thus bi_elem_idx is null
            !% BClat BCs are associated with elements, thus bi_face_idx is null
            BC%flowI(:, bi_face_idx) = nullvalueI
            BC%flowI(:, bi_elem_idx) = nullvalueI
        end if
        if (N_headBC > 0) then
            BC%headI(:,bi_fetch) = 1
            BC%headIdx(:) = 0
        end if

        !% Initialize Inflow BCs
        if (N_flowBC > 0) then
            do ii = 1, N_flowBC
                nidx = node%P%have_flowBC(ii)
                ntype = node%I(nidx, ni_node_type)

                !% Check if node has inflow BC
                has_dwf_inflow = &
                    (interface_get_node_attribute(nidx, api_node_has_dwfInflow) == 1)
                has_ext_inflow = &
                    (interface_get_node_attribute(nidx, api_node_has_extInflow) == 1)

                print *, nidx, api_node_has_extInflow, has_dwf_inflow, has_ext_inflow
                !% Handle Inflow BCs (BCup and BClat only)
                if (has_dwf_inflow .or. has_ext_inflow) then
                    if ((ntype == nJm) .or. (ntype == nJ2)) then
                        BC%flowI(ii, bi_category) = BClat
                        BC%flowI(ii, bi_elem_idx) = node%I(nidx, ni_elemface_idx) !% elem idx
                    else if (ntype == nBCup) then
                        BC%flowI(ii, bi_category) = BCup
                        BC%flowI(ii, bi_face_idx) = node%I(nidx, ni_elemface_idx) !% face idx
                    else
                        print *, "Error, BC type can't be an inflow BC for node " // node%Names(nidx)%str
                        stop
                    end if

                    BC%flowI(ii, bi_node_idx) = nidx
                    BC%flowI(ii, bi_idx) = ii
                    nbasepat = &
                        interface_get_node_attribute(nidx, api_node_extInflow_basePat)
                    ntseries = &
                        interface_get_node_attribute(nidx, api_node_extInflow_tSeries)

                    !% BC does not have fixed value if its associated with dwfInflow
                    !% or if extInflow has tseries or pattern
                    BC%flowI(ii, bi_subcategory) = BCQ_tseries
                    if (.not. has_dwf_inflow) then !% extInflow only
                        if ((ntseries == -1) .and. (nbasepat /= -1)) then
                            BC%flowI(ii, bi_subcategory) = BCQ_fixed
                        end if
                    end if
                else
                    print *, "There is an error, only nodes with extInflow or dwfInflow can have inflow BC"
                    stop
                end if
            end do
        end if

        !% Initialize Head BCs
        if (N_headBC > 0) then
            do ii = 1, N_headBC
                nidx = node%P%have_headBC(ii)
                ntype = node%I(nidx, ni_node_type)

                if (ntype == nBCdn) then
                    BC%headI(ii, bi_category) = BCdn
                    BC%headI(ii, bi_face_idx) = node%I(ii, ni_elemface_idx) !% face idx
                else
                    print *, "Error, BC type can't be a head BC for node " // node%Names(nidx)%str
                    stop
                end if

                BC%headI(ii, bi_idx) = ii
                BC%headI(ii, bi_node_idx) = nidx
                if (interface_get_node_attribute(nidx, ni_node_subtype) == API_FREE_OUTFALL) then
                    BC%headI(ii, bi_subcategory) = BCH_free
                else if (interface_get_node_attribute(nidx, ni_node_subtype) == API_NORMAL_OUTFALL) then
                    BC%headI(ii, bi_subcategory) = BCH_normal
                else if (interface_get_node_attribute(nidx, ni_node_subtype) == API_FIXED_OUTFALL) then
                    BC%headI(ii, bi_subcategory) = BCH_fixed
                else if (interface_get_node_attribute(nidx, ni_node_subtype) == API_TIDAL_OUTFALL) then
                    BC%headI(ii, bi_subcategory) = BCH_tidal
                else if (interface_get_node_attribute(nidx, ni_node_subtype) == API_TIMESERIES_OUTFALL) then
                    BC%headI(ii, bi_subcategory) = BCH_tseries
                end if
            end do
        end if
        if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name
    end subroutine init_bc
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

        !% Set the network partitioning method used for multi-processor parallel computation
        call init_partitioning_method()

        !% adjust the link lengths by cutting off a certain portion for the junction branch
        !% this subroutine is called here to correctly estimate the number of elements and faces
        !% to allocate the coarrays.
        !% HACK: This might be moved someplace more suitable?
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
            node_index = PACK([(counter, counter=1,size(node%I,1))], node%I(:, ni_P_image) == unique_imagenum(ii))
            link_index = PACK([(counter, counter=1,size(link%I,1))], link%I(:, li_P_image) == unique_imagenum(ii))
            ! create corresponding indices for node and link in this image

            ! The number of elements and faces is actually decided by the junctions
            ! So we will calculate the number of junction and base on different scenarios to decided
            ! how many elem/face are assigned to each image
            junction_counter = count(node%I(node_index, ni_node_type) == nJm)

            !% first calculate the number of nodes in each partition, assign elems/faces for junctions
            elem_counter = elem_counter + J_elem_add * junction_counter
            face_counter = face_counter + J_face_add * junction_counter

            !% loop through the links and calculate the internal faces between elements
            do jj = 1, size(link_index,1)
                idx = link_index(jj)
                face_counter = face_counter + link%I(idx, li_N_element) - 1 !% internal faces between elems, e.g. 5 elements have 4 internal faces
                elem_counter = elem_counter + link%I(idx, li_N_element) ! number of elements
            enddo

            !% now we loop through the nodes and count the node faces
            do jj = 1, size(node_index,1)
                idx = node_index(jj)
                if (node%I(idx, ni_node_type) == nJ2) then
                    face_counter = face_counter + 1 !% add the face of 1-to-1 junction between 2 links
                elseif (node%I(idx, ni_node_type) == nBCup) then
                    face_counter = face_counter +1 !% add the upstream faces
                elseif (node%I(idx, ni_node_type) == nBCdn) then
                    face_counter = face_counter +1 !% add the downstream faces
                endif !% multiple junction faces already counted
            enddo

            !% Now we count the space for duplicated faces
            do jj = 1, size(link_index,1)
                idx = link_index(jj)
                !% check upstream node first
                if ( ( node%I(link%I(idx, li_Mnode_u), ni_P_is_boundary) == 1) .and. &
                    ( node%I(link%I(idx, li_Mnode_u), ni_P_image) .ne. ii) ) then
                    face_counter = face_counter +1
                    duplicated_face_counter = duplicated_face_counter + 1
                endif
                ! then downstream node
                if ( ( node%I(link%I(idx, li_Mnode_d), ni_P_is_boundary) == 1) .and. &
                    ( node%I(link%I(idx, li_Mnode_d), ni_P_image) .ne. ii) ) then
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

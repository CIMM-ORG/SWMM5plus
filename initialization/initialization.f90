module initialization
    use boundary_conditions
    use define_keys
    use define_globals
    use define_settings
    use define_indexes
    use discretization
    use initial_condition
    use interface
    use network_define
    use partitioning
    use pack_mask_arrays, only: pack_nodes
    use utility_allocate
    use utility_array
    use utility_output
    use utility_array
    use utility_profiler
    !use utility_prof_jobcount
    use pack_mask_arrays
    use output

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
    !%
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    subroutine initialize_all()
    !%-----------------------------------------------------------------------------
    !%
    !% Description:
    !%   a public subroutine that calls all the private initialization subroutines
    !%
    !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'initialize_all'
    !%-----------------------------------------------------------------------------
        if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        !% ---  Define project & settings paths
        call getcwd(setting%Paths%project)
        setting%paths%setting = trim(setting%Paths%project) // "/definitions/settings.json"

        !% read and store the command-line options
        call init_read_arguments ()

        !% set the branchsign global -- this is used for junction branches (JB)
        !% for upstream (+1) and downstream (-1)
        !% HACK: for clarity and consistency, this probably should be moved into
        !% the init_network. Placed here for the time being in case we need it
        !% for translating link/node from SWMM-C or partitioning.
        branchsign(1:max_branch_per_node-1:2) = +oneR
        branchsign(2:max_branch_per_node:2)   = -oneR

        !% load the settings.json file with the default setting% model control structure
        !% def_load_settings is one of the few subroutines in the Definition modules
        call def_load_settings()

        !% read and store the command-line options
        call init_read_arguments ()

        if (setting%Verbose) print *, "Simulation Starts"

        !% set up the profiler
        if (setting%Profile%YN) then
            call util_allocate_profiler ()   
            call util_profiler_start (pfc_initialize_all)
        end if

        !% initialize the API with the SWMM-C code
        call interface_init ()

        !if (setting%Verbose) print *, "begin link-node processing"

        !% set up and store the SWMM-C link-node arrays in equivalent Fortran arrays
        call init_linknode_arrays ()

        !if (setting%Verbose) print *, "begin partitioning"

        call init_partitioning()

        !% HACK -- to this point the above could all be done on image(1) and then
        !% distributed to the other images. This might create problems in ensuring
        !% that all the data gets copied over when new stuff is added. Probably OK
        !% to keep the above as computing on all images until the code is near complete

        !% HACK: this sync call is probably not needed
        sync all

        !if (setting%Verbose) print *, "begin network definition"

        call init_network_define_toplevel ()

        !if (setting%Verbose) print *, "begin reading csv"

        !% read in link names for output
        call output_read_csv_link_names()
        call output_read_csv_node_names()

        !if (setting%Verbose) print *, "begin initializing boundary conditions"

        !% initialize boundary conditions
        call init_bc()

        call init_time()

        if (setting%Verbose) then
            if (this_image() == 1) then
            if ((N_link > 5000) .or. (N_node > 5000)) then
                print *, "begin setting initial conditions (this takes several minutes for big systems)"
                print *, "This system has ", SWMM_N_link, " links and ", SWMM_N_node, " nodes"
                print *, "The finite-volume system is ", sum(N_elem(:)), " elements"
            endif
            endif
        endif
        call init_IC_setup ()

        !if (setting%Verbose) print *, "begin setup of output files"

        !% creating output_folders and files
        call util_output_clean_folders()
        call util_output_create_folders()

        if ((this_image() == 1) .and. setting%Debug%Input) call util_output_export_linknode_input()
        if (setting%Debug%Output) then
            call util_output_create_elemR_files()
            call util_output_create_faceR_files()
            call util_output_create_summary_files()
        end if
        if (setting%Debug%Output .or. setting%Output%report) then
            call output_create_link_files()
            call output_create_node_files()
        end if

        if (setting%Profile%YN) call util_profiler_stop (pfc_initialize_all)

        !% wait for all the processors to reach this stage before starting the time loop
        sync all
        
        if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine initialize_all
    !%
    !%==========================================================================
    !% PRIVATE
    !%==========================================================================
    !%
    subroutine init_linknode_arrays()
    !%-----------------------------------------------------------------------------
    !%
    !% Description:
    !%   Retrieves data from EPA-SWMM interface and populates link and node tables
    !% Note:
    !%   The order in which link and nodes are populated coincides with the
    !%   order in which links and nodes are allocated in EPA-SWMM data structures
    !%   Keeping the same order is important to be able to locate node/link data
    !%   by label and not by index, reusing EPA-SWMM functionalities.
    !%-----------------------------------------------------------------------------

        integer       :: ii, total_n_links

        character(64) :: subroutine_name = 'init_linknode_arrays'

    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        if (.not. api_is_initialized) then
            print *, "ERROR: API is not initialized"
            stop "in " // subroutine_name
        end if

        !% Allocate storage for link & node tables
        call util_allocate_linknode()

        link%I(:,li_num_phantom_links) = 0
        node%I(:,ni_N_link_u) = 0
        node%I(:,ni_N_link_d) = 0

        do ii = 1, SWMM_N_link
            link%I(ii,li_idx) = ii
            link%I(ii,li_link_type) = interface_get_link_attribute(ii, api_link_type)
            link%I(ii,li_weir_type) = interface_get_link_attribute(ii, api_weir_type)
            link%I(ii,li_orif_type) = interface_get_link_attribute(ii, api_orifice_type)
            link%I(ii,li_pump_type) = interface_get_link_attribute(ii, api_pump_type)
            link%I(ii,li_geometry) = interface_get_link_attribute(ii, api_link_geometry)
            link%I(ii,li_Mnode_u) = interface_get_link_attribute(ii, api_link_node1) + 1 ! node1 in C starts from 0
            link%I(ii,li_Mnode_d) = interface_get_link_attribute(ii, api_link_node2) + 1 ! node2 in C starts from 0
            link%I(ii,li_parent_link) = ii

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

            !% link%R(ii,lr_TopWidth): defined in network_define.f08
            link%R(ii,lr_BreadthScale) = interface_get_link_attribute(ii, api_link_xsect_wMax)
            !% link%R(ii,lr_Slope): defined in network_define.f08
            link%R(ii,lr_LeftSlope) = interface_get_link_attribute(ii, api_link_left_slope)
            link%R(ii,lr_RightSlope) = interface_get_link_attribute(ii, api_link_right_slope)
            link%R(ii,lr_Roughness) = interface_get_link_attribute(ii, api_conduit_roughness)
            link%R(ii,lr_InitialFlowrate) = interface_get_link_attribute(ii, api_link_q0)
            link%R(ii,lr_InitialUpstreamDepth) = interface_get_node_attribute(link%I(ii,li_Mnode_u), api_node_initDepth)
            link%R(ii,lr_InitialDnstreamDepth) = interface_get_node_attribute(link%I(ii,li_Mnode_d), api_node_initDepth)
            link%R(ii,lr_InitialDepth) = (link%R(ii,lr_InitialDnstreamDepth) + link%R(ii,lr_InitialUpstreamDepth)) / 2.0
            link%R(ii,lr_FullDepth) = interface_get_link_attribute(ii, api_link_xsect_yFull)
            link%R(ii,lr_InletOffset) = interface_get_link_attribute(ii,api_link_offset1)
            link%R(ii,lr_OutletOffset) = interface_get_link_attribute(ii,api_link_offset2)

            !% special element attributes
            link%I(ii,li_weir_EndContrations) = interface_get_link_attribute(ii, api_weir_end_contractions)
            link%R(ii,lr_DischargeCoeff1) = interface_get_link_attribute(ii, api_discharge_coeff1)
            link%R(ii,lr_DischargeCoeff2) = interface_get_link_attribute(ii, api_discharge_coeff2)

            !% SWMM5 doesnot distinct between channel and conduit
            !% however we need that distinction to set up the init condition
            if ( (link%I(ii,li_link_type) == lPipe)          .and. &
                 ( &
                 (link%I(ii,li_geometry) == lRectangular)    .or. &
                 (link%I(ii,li_geometry) == lTrapezoidal)    .or. &
                 (link%I(ii,li_geometry) == lPower_function) .or. &
                 (link%I(ii,li_geometry) == lRect_triang)    .or. &
                 (link%I(ii,li_geometry) == lRect_round)     .or. &
                 (link%I(ii,li_geometry) == lMod_basket)     .or. &   
                 (link%I(ii,li_geometry) == lIrregular)) ) then

                link%I(ii,li_link_type) = lChannel
            end if
        end do

        do ii = 1, N_node
            total_n_links = node%I(ii,ni_N_link_u) + node%I(ii,ni_N_link_d)
            node%I(ii, ni_idx) = ii
            if (interface_get_node_attribute(ii, api_node_type) == API_OUTFALL) then
                node%I(ii, ni_node_type) = nBCdn
            else if ((total_n_links == twoI)          .and. &
                     (node%I(ii,ni_N_link_u) == oneI) .and. &
                     (node%I(ii,ni_N_link_d) == oneI) )then
                node%I(ii, ni_node_type) = nJ2
            else if (total_n_links >= twoI) then
                node%I(ii, ni_node_type) = nJm
            end if

            node%YN(ii, nYN_has_extInflow) = interface_get_node_attribute(ii, api_node_has_extInflow) == 1
            node%YN(ii, nYN_has_dwfInflow) = interface_get_node_attribute(ii, api_node_has_dwfInflow) == 1

            if (node%YN(ii, nYN_has_extInflow) .or. node%YN(ii, nYN_has_dwfInflow)) then
                node%YN(ii, nYN_has_inflow) = .true.
                if ((node%I(ii,ni_N_link_u) == zeroI) .and. (total_n_links == oneI)) then
                    node%I(ii, ni_node_type) = nBCup
                end if
            end if

            node%R(ii,nr_InitialDepth) = interface_get_node_attribute(ii, api_node_initDepth)
            node%R(ii,nr_Zbottom) = interface_get_node_attribute(ii, api_node_invertElev)
            node%I(ii, ni_pattern_resolution) = interface_get_BC_resolution(ii)
        end do

        !% Update Link/Node names
        call interface_update_linknode_names()

    end subroutine init_linknode_arrays
    !%
    !%==========================================================================
    !%==========================================================================
    !%
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
        character(64) :: subroutine_name = "init_bc"
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        if (setting%Profile%YN) call util_profiler_start (pfc_init_bc)

        call pack_nodes()
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
            BC%flowR_timeseries = nullValueR
        end if
        if (N_headBC > 0) then
            BC%headI = nullvalueI
            BC%headI(:,bi_fetch) = 1
            BC%headIdx(:) = 0
            BC%headR_timeseries = nullValueR
        end if

        !% Initialize Inflow BCs
        if (N_flowBC > 0) then
            do ii = 1, N_flowBC
                nidx = node%P%have_flowBC(ii)
                ntype = node%I(nidx, ni_node_type)

                !% Handle Inflow BCs (BCup and BClat only)
                if (node%YN(nidx, nYN_has_extInflow) .or. node%YN(nidx, nYN_has_dwfInflow)) then
                    if ((ntype == nJm) .or. (ntype == nJ2)) then
                        BC%flowI(ii, bi_category) = BClat
                        BC%flowI(ii, bi_elem_idx) = node%I(nidx, ni_elemface_idx) !% elem idx
                    else if (ntype == nBCup) then
                        BC%flowI(ii, bi_category) = BCup
                        BC%flowI(ii, bi_face_idx) = node%I(nidx, ni_elemface_idx) !% face idx
                    else
                        print *, "Error, BC type can't be an inflow BC for node " // node%Names(nidx)%str
                        stop "in " // subroutine_name
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
                    if (.not. node%YN(nidx, nYN_has_dwfInflow)) then !% extInflow only
                        if ((ntseries == -1) .and. (nbasepat /= -1)) then
                            BC%flowI(ii, bi_subcategory) = BCQ_fixed
                        end if
                    end if
                else
                    print *, "There is an error, only nodes with extInflow or dwfInflow can have inflow BC"
                    stop "in " // subroutine_name
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
                    BC%headI(ii, bi_face_idx) = node%I(nidx, ni_elemface_idx) !% face idx
                else
                    print *, "Error, BC type can't be a head BC for node " // node%Names(nidx)%str
                    stop "in " // subroutine_name
                end if

                BC%headI(ii, bi_idx) = ii
                BC%headI(ii, bi_node_idx) = nidx

                if (interface_get_node_attribute(nidx, api_node_outfall_type) == API_FREE_OUTFALL) then
                    BC%headI(ii, bi_subcategory) = BCH_free
                else if (interface_get_node_attribute(nidx, api_node_outfall_type) == API_NORMAL_OUTFALL) then
                    BC%headI(ii, bi_subcategory) = BCH_normal
                else if (interface_get_node_attribute(nidx, api_node_outfall_type) == API_FIXED_OUTFALL) then
                    BC%headI(ii, bi_subcategory) = BCH_fixed
                else if (interface_get_node_attribute(nidx, api_node_outfall_type) == API_TIDAL_OUTFALL) then
                    BC%headI(ii, bi_subcategory) = BCH_tidal
                else if (interface_get_node_attribute(nidx, api_node_outfall_type) == API_TIMESERIES_OUTFALL) then
                    BC%headI(ii, bi_subcategory) = BCH_tseries
                end if
            end do
        end if

        call bc_step()
        call pack_bc()

        if (setting%Profile%YN) call util_profiler_stop (pfc_init_bc)

        if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine init_bc
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine init_partitioning()
    !%-----------------------------------------------------------------------------
    !%
    !% Description:
    !%   This subroutine calls the public subroutine from the utility module,
    !%   partitioning.f08. It also calls a public subroutine from the temporary
    !%   coarray_partition.f08 utility module that defines how big the coarrays
    !%   must be.
    !%
    !%-----------------------------------------------------------------------------
        integer       :: ii
        character(64) :: subroutine_name = 'init_partitioning'
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        if (setting%Profile%YN) call util_profiler_start (pfc_init_partitioning) 

        !% find the number of elements in a link based on nominal element length
        do ii = 1, SWMM_N_link
            call init_discretization_nominal(ii)
        end do

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

        if (setting%Profile%YN) call util_profiler_stop (pfc_init_partitioning)

        if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

    end subroutine init_partitioning
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine init_coarray_length()
        !% for coarray length determination
        integer :: nimgs_assign
        integer, allocatable :: unique_imagenum(:)
        integer :: ii, jj, kk, idx, counter, elem_counter=0, face_counter=0, junction_counter=0

        integer :: duplicated_face_counter=0
        integer, allocatable :: node_index(:), link_index(:), temp_arr(:)
        character(64) :: subroutine_name = 'init_coarray_length'

        if (setting%Debug%File%utility_array) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        call util_image_number_calculation(nimgs_assign, unique_imagenum)

        allocate(N_elem(num_images()))
        allocate(N_face(num_images()))
        allocate(N_unique_face(num_images()))

        do ii=1, num_images()

            node_index = PACK([(counter, counter=1,size(node%I,1))], node%I(:, ni_P_image) == unique_imagenum(ii))
            link_index = PACK([(counter, counter=1,size(link%I,1))], link%I(:, li_P_image) == unique_imagenum(ii))
            !% create corresponding indices for node and link in this image

            !% The number of elements and faces is actually decided by the junctions
            !% So we will calculate the number of junction and base on different scenarios to decided
            !% how many elem/face are assigned to each image
            junction_counter = count(node%I(node_index, ni_node_type) == nJm)

            !% first calculate the number of nodes in each partition, assign elems/faces for junctions
            elem_counter = elem_counter + J_elem_add * junction_counter
            face_counter = face_counter + J_face_add * junction_counter

            !% loop through the links and calculate the internal faces between elements
            do jj = 1, size(link_index,1)
                idx = link_index(jj)
                face_counter = face_counter + link%I(idx, li_N_element) - 1 !% internal faces between elems, e.g. 5 elements have 4 internal faces
                elem_counter = elem_counter + link%I(idx, li_N_element) ! number of elements
            end do

            !% now we loop through the nodes and count the node faces
            do jj = 1, size(node_index,1)
                idx = node_index(jj)
                if (node%I(idx, ni_node_type) == nJ2) then
                    face_counter = face_counter + 1 !% add the face of 1-to-1 junction between 2 links
                elseif (node%I(idx, ni_node_type) == nBCup) then
                    face_counter = face_counter +1 !% add the upstream faces
                elseif (node%I(idx, ni_node_type) == nBCdn) then
                    face_counter = face_counter +1 !% add the downstream faces
                end if !% multiple junction faces already counted
            end do

            !% Now we count the space for duplicated faces
            do jj = 1, size(link_index,1)
                idx = link_index(jj)
                !% check upstream node first
                if ( ( node%I(link%I(idx, li_Mnode_u), ni_P_is_boundary) == 1) .and. &
                    ( node%I(link%I(idx, li_Mnode_u), ni_P_image) .ne. ii) ) then
                    face_counter = face_counter +1
                    duplicated_face_counter = duplicated_face_counter + 1
                end if
                !% then downstream node
                if ( ( node%I(link%I(idx, li_Mnode_d), ni_P_is_boundary) == 1) .and. &
                    ( node%I(link%I(idx, li_Mnode_d), ni_P_image) .ne. ii) ) then
                    face_counter = face_counter +1
                    duplicated_face_counter = duplicated_face_counter + 1
                end if
            end do

            N_elem(ii) = elem_counter
            N_face(ii) = face_counter
            N_unique_face(ii) = face_counter - duplicated_face_counter

            elem_counter = zeroI ! reset the counter
            face_counter = zeroI
            junction_counter = zeroI
            duplicated_face_counter = zeroI

        end do

        max_caf_elem_N = maxval(N_elem)
        max_caf_face_N = maxval(N_face) ! assign the max value

        if (setting%Debug%File%utility_array) then
            do ii = 1, size(unique_imagenum,1)
                print*, 'Processor => ', ii
                print*, 'Elements expected ', N_elem(ii)
                print*, 'Faces expected    ', N_face(ii)
            end do
        end if

        if (setting%Debug%File%utility_array)  &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

    end subroutine init_coarray_length
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine init_read_arguments()
        integer :: ii
        logical :: arg_param = .false.
        character(len=8) :: param
        character(len=256) :: arg
        character(64) :: subroutine_name = "init_read_arguments"

        if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        do ii = 1, iargc()
            call getarg(ii, arg)
            if (.not. arg_param) then
                param = arg
                if (ii == 1) then
                    if (arg(:1) == '-') then
                        print *, "ERROR: it is necessary to define the path to the .inp file"
                        stop "in " // subroutine_name
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
                    stop "in " // subroutine_name
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
                        stop "in " // subroutine_name
                    end if
                elseif (trim(param) == '--run-tests') then
                end if
            end if
        end do

        if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine init_read_arguments
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine init_time()
        logical :: doHydraulics

        setting%Time%Dt = setting%Time%Hydraulics%Dt
        setting%Time%Now = 0
        setting%Time%Step = 0
        setting%Time%Hydraulics%Step = 0
        setting%Time%Hydrology%Step = 0
        if (.not. setting%Simulation%useHydrology) setting%Time%Hydrology%Dt = nullValueR
        !% Initialize report step
        setting%Output%reportStep = int(setting%Output%reportStartTime / setting%Output%reportDt)

        if (setting%Time%Hydrology%Dt < setting%Time%Hydraulics%Dt) then
            stop "Error: Hydrology time step can't be smaller than hydraulics time step"
        end if
    end subroutine init_time
    !%
    !%==========================================================================
    !% END OF MODULE
    !%==========================================================================
    !%
end module initialization

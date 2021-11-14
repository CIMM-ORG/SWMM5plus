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
    use utility_allocate
    use utility_array
    use utility_datetime
    use utility_output
    use utility_profiler
    use utility_files
    use pack_mask_arrays
    use output

    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !%    General initialization of data structures (not including network)
    !%
    !% Method:
    !%    Creates the arrays index structures that are used for accessing data.
    !%    Arguably, this could be done more simply but we want the fundamental
    !%    column indexes in array_index to be parameters rather than variables. By
    !%    using parameters we reduce the possibility of accidentally changing a
    !%    column definition.
    !%
    !% Note on naming:
    !%    The driver subroutine is named after the driver module (in this case,
    !%    initialization).  Subsequent subroutines are name such that the subroutine
    !%    name is essentially a path "init_<module>_<subroutine_name>"
    !%-----------------------------------------------------------------------------

    private
    public :: initialize_toplevel

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine initialize_toplevel ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !%   a public subroutine that calls all the private initialization subroutines
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'initialize_toplevel'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        !% --- CPU and wall-clock timers
        call init_model_timer()

        !% --- assign and store unit numbers for files
        call util_file_assign_unitnumber ()

        !% --- get command line assignments and store
        call util_file_get_commandline ()

        !% --- input project paths and filenames from command line arguments
        !%     note that all files and folders must exist
        ! call util_file_setup_input_paths_and_files()

        if (setting%Output%Verbose) &
            write(*,"(2A,i5,A)") new_line(" "), 'begin initialization [Processor ', this_image(), "] ..."

        !% --- set the branchsign global -- this is used for junction branches (JB)
        !%     for upstream (+1) and downstream (-1)
        !%     HACK: for clarity and consistency, this probably should be moved into
        !%     the init_network. Placed here for the time being in case we need it
        !%     for translating link/node from SWMM-C or partitioning.
        branchsign(1:max_branch_per_node-1:2) = +oneR
        branchsign(2:max_branch_per_node:2)   = -oneR

        !% --- load the settings.json file with the default setting% model control structure
        !%     def_load_settings is one of the few subroutines in the Definition modules
        call def_load_settings()

        !% --- initialize the time stamp used for output (must be after json is read)
        call init_timestamp ()


        !% HACK
        !% --- Read and process the command-line options a second time to prevent overwrite
        !%     from json file and reprocess.
        !%     Possibly replace this later with a set of settings that are saved
        !%     and simply overwrite after the settings.json is loaded.
        call util_file_assign_unitnumber ()
        call util_file_get_commandline ()
        call util_file_setup_input_paths_and_files()

        !% --- output file directories
        call util_file_setup_output_folders()

        !%  --- finish setting all the paths to output folders
        sync all

        if (setting%Output%Verbose) then
            write(*,"(A)") "Simulation Starts..."
            write(*,"(A)") 'Using the following files:'
            write(*,"(A)") 'Input file  : '//trim(setting%File%inp_file)
            write(*,"(A)") 'Report file : '//trim(setting%File%rpt_file)
            write(*,"(A)") 'Output file : '//trim(setting%File%out_file)
            write(*,"(A)") 'Settings file: '//trim(setting%File%setting_file)
        end if

        !% --- set up the profiler
        if (setting%Profile%YN) then
            call util_allocate_profiler ()
            call util_profiler_start (pfc_initialize_all)
        end if

        !%  --- finish setting all the file paths before initialing the interface
        sync all
        !% --- initialize the API with the SWMM-C code
        call interface_init ()

        !% --- set up and store the SWMM-C link-node arrays in equivalent Fortran arrays
        if (setting%Output%Verbose) print *, "begin link-node processing"
        call init_linknode_arrays ()

        !% --- store the SWMM-C curves in equivalent Fortran arrays
        if (setting%Output%Verbose) print *, "begin SWMM5 curve processing"
        call init_curves()

        if (setting%Output%Verbose) print *, "begin partitioning"
        call init_partitioning()

        !% --- HACK -- to this point the above could all be done on image(1) and then
        !% distributed to the other images. This might create problems in ensuring
        !% that all the data gets copied over when new stuff is added. Probably OK
        !% to keep the above as computing on all images until the code is near complete

        sync all

        if (setting%Output%Verbose) print *, "begin network definition"
        call init_network_define_toplevel ()

        !% --- HACK --- NEEDS FIXING FOR MULTILEVEL OUTPUT WITH MULTIPLE PROCESSOR
        ! if (num_images() > 1) then
        !     print *, "ERROR (code) -- NEED TO FIX HACK"
        !     stop 9347044
        ! else
        !     !% temporarily assign the node index to the Global swmm index
        !     !% only good for one processor operation
        !     faceI(:,fi_node_Gidx_SWMM) = faceI(:,fi_node_idx)
        ! end if
        !% ---- END HACK

        !% --- setup for csv output of links and nodes
        !brh20211006 call outputD_read_csv_link_names()
        !brh20211006 call outputD_read_csv_node_names()

        !% --- initialize boundary condition
        if (setting%Output%Verbose) print *, "begin initializing boundary conditions"
        call init_bc()
        call init_time()

        if (setting%Output%Verbose) then
            if (this_image() == 1) then
            if ((N_link > 5000) .or. (N_node > 5000)) then
                print *, "begin setting initial conditions (this takes several minutes for big systems)"
                print *, "This system has ", SWMM_N_link, " links and ", SWMM_N_node, " nodes"
                print *, "The finite-volume system is ", sum(N_elem(:)), " elements"
            endif
            endif
        endif
        call init_IC_toplevel ()

        print *
        print *, 'WORK NEEDED HERE ',73794
        print *
        !% --- designate/select the nodes/links for output
        call output_COMMON_nodes_selection ()
        call output_COMMON_links_selection ()
        !% --- designate the corresponding elements for output
        call outputML_element_selection ()
        !% --- deisgnate the corresponding face to output
        call outputML_face_selection ()
        !% --- create packed arrays of elem row numbers that are output
        call pack_element_outputML ()
        !% --- create packed arrays of face row numbers that are output
        call pack_face_outputML ()

        !print *, elemP(1:npack_elemP(ep_Output_Elements),ep_Output_Elements)
        !print *
        !print *, faceP(1:npack_faceP(fp_Output_Faces),fp_Output_Faces)
        !stop 34780

        !% --- compute the N_OutElem for each image
        call outputML_size_OutElem_by_image ()
        !% --- compute the N_OutFace for each imaige
        call outputML_size_OutFace_by_image ()
        !% --- setup the output element data types
        call outputML_element_outtype_selection ()
        !% -- setup the output face data types
        call outputML_face_outtype_selection ()
        !% --- create storage space for multi-level output data
        call util_allocate_outputML_storage ()
        !% --- create storage for output times
        call util_allocate_outputML_times ()
        !% --- create storage for output binary filenames
        call util_allocate_outputML_filenames ()

        !if (setting%Output%Verbose) print *, "begin setup of output files"

        !% creating output_folders and files
        !% brh 20211004 -- moved this functionality into util_file_setup_output_files
        !call util_output_clean_folders()
        !call util_output_create_folders()

        !brh20211006 COMMENTING OUT THE OUTPUT BY CSV
        !brh20211006 if ((this_image() == 1) .and. setting%Debug%Setup) call util_output_export_linknode_input()
        !brh20211006 if (setting%Debug%Output) then
        !brh20211006     call util_output_create_elemR_files()
        !brh20211006     call util_output_create_faceR_files()
        !brh20211006     call util_output_create_summary_files()
        !brh20211006 end if
        !brh20211006 if (setting%Debug%Output .or. setting%Output%report) then
        !brh20211006     call output_create_link_files()
        !brh20211006     call output_create_node_files()
        !brh20211006 end if

        ! if (setting%Profile%YN) call util_profiler_stop (pfc_initialize_all)

        !% wait for all the processors to reach this stage before starting the time loop
        sync all

        if (icrash) then  !% if crash in initialization, write the output and exit
            call outputML_store_data (.true.)
            return
        end if

        if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine initialize_toplevel
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine init_model_timer()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% starts and stores the CPU clock and wall clock time
        !%-----------------------------------------------------------------------------
        !% store CPU clock start time
        call cpu_time(setting%Time%CPU%EpochStartSeconds)

        !% store Real time
        call system_clock(count=setting%Time%Real%EpochStartSeconds)

    end subroutine init_model_timer
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_timestamp ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% initializes time stamp used for output files
        !%-----------------------------------------------------------------------------
        integer :: thisTime(8), ii, thisunit, ios
        character(len=4) :: cyear
        character(len=2) :: cmonth, cday, chour, cmin
        character(len=13) :: datetimestamp
        character(64) :: subroutine_name = 'init_timestamp'
        !%-----------------------------------------------------------------------------
        if (icrash) return
        call date_and_time(values = thisTime)
        write(cyear, "(i4)") thisTime(1)
        write(cmonth,"(i2)") thisTime(2)
        write(cday,  "(i2)") thisTime(3)
        write(chour, "(i2)") thisTime(5)
        write(cmin,  "(i2)") thisTime(6)
        if (thisTime(2) < 10) then
            cmonth = '0'//adjustl(cmonth)
        end if
        if (thisTime(3) < 10) then
            cday = '0'//adjustl(cday)
        end if
        if (thisTime(5) < 10) then
            chour = '0'//adjustl(chour)
        end if
        if (thisTime(6) < 10) then
            cmin = '0'//adjustl(cmin)
        end if

        if (this_image() == 1) then
            datetimestamp = cyear//cmonth//cday//'_'//chour//cmin
        endif

        call co_broadcast (datetimestamp, source_image=1)

        setting%Time%DateTimeStamp = datetimestamp

        ! print*, 'image', this_image()
        ! print *, setting%Time%DateTimeStamp

        ! !% --- distribute to all processors
        ! !% --- HACK using a write/read file as the setting varialble is not a coarray
        ! if (this_image() == 1) then
        !     open(newunit = thisunit, &
        !         file = 'temp_fortran.txt',    &
        !         action = 'write', &
        !         iostat = ios)
        !     print*, 'ios', ios
        !     if (ios /= 0) then
        !         write(*,"(A)") 'ERROR (CODE) file temp_fortran.txt could not be opened for writing.'
        !         write(*,"(A)") 'File purpose is write/reading for syncing non-coarrays across images'
        !         stop
        !     end if
        !     write(thisunit,"(A)") setting%Time%DateTimeStamp
        !     close(thisunit)
        ! end if
        ! !% testing
        ! !open(newunit = thisunit, &
        ! !    file = 'temp_fortran.txt',    &
        ! !    action = 'read', &
        ! !    iostat = ios)
        ! !read(thisunit,"(A)")  datetimestamp
        ! !print *, datetimestamp

        ! !% read sequentially into other images
        ! do ii = 2,num_images()
        !     open(newunit = thisunit, &
        !         file = 'temp_fortran.txt',    &
        !         action = 'read', &
        !         iostat = ios)
        !     if (ios /= 0) then
        !         write(*,"(A)") 'ERROR (CODE) temp_fortran.txt file could not be opened for reading.'
        !         write(*,"(A)") 'File purpose is write/reading for syncing non-coarrays across images'
        !         stop
        !     end if
        !     read(thisunit,"(A)") setting%Time%DateTimeStamp
        ! end do

        sync all

    end subroutine init_timestamp
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_linknode_arrays()
        !%-----------------------------------------------------------------------------
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
        if (icrash) return
        if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        if (.not. api_is_initialized) then
            print *, "ERROR: API is not initialized"
            stop
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
            else if (interface_get_node_attribute(ii, api_node_type) == API_STORAGE) then
                node%I(ii, ni_node_type) = nJm
                node%YN(ii, nYN_has_storage) = .true.
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

            node%R(ii,nr_InitialDepth)      = interface_get_node_attribute(ii, api_node_initDepth)
            node%R(ii,nr_Zbottom)           = interface_get_node_attribute(ii, api_node_invertElev)
            node%R(ii,nr_FullDepth)         = interface_get_node_attribute(ii, api_node_fullDepth)
            node%R(ii,nr_StorageConstant)   = interface_get_node_attribute(ii, api_node_StorageConstant)
            node%R(ii,nr_StorageCoeff)      = interface_get_node_attribute(ii, api_node_StorageCoeff)
            node%R(ii,nr_StorageExponent)   = interface_get_node_attribute(ii, api_node_StorageExponent)
            node%I(ii,ni_curve_ID)          = interface_get_node_attribute(ii, api_node_StorageCurveID) + 1
            node%I(ii,ni_pattern_resolution) = interface_get_BC_resolution(ii)
        end do

        !% Update Link/Node names
        call interface_update_linknode_names()

        if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

    end subroutine init_linknode_arrays
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_curves()
        !%-----------------------------------------------------------------------------
        !% Description:
        !%   Retrieves data from EPA-SWMM interface and populates curve curves
        !%-----------------------------------------------------------------------------

        integer       :: ii, jj, additional_storage_curves, Total_curves

        character(64) :: subroutine_name = 'init_curves'

        !%-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        if (.not. api_is_initialized) then
            print *, "ERROR: API is not initialized"
            stop
        end if

        !% we create additional curves for functional storage as well
        !% this allocates the space for functional storage curve
        additional_storage_curves = count((node%YN(:, nYN_has_storage)) .and. &
                                          (node%I(:,ni_curve_ID) == 0))

        Total_Curves = additional_storage_curves + SWMM_N_Curve
        if (Total_Curves > SWMM_N_Curve) N_Curve = SWMM_N_Curve + additional_storage_curves

        !% allocate the number of curve objets from SWMM5
        call util_allocate_curves()

        do ii = 1, SWMM_N_Curve
            curve(ii)%ID = ii
            curve(ii)%Type = interface_get_table_attribute(ii, api_table_type)
            !% get the number of entries in a curve
            curve(ii)%NumRows = interface_get_num_table_entries(ii)
            !% allocate the value space
            call util_allocate_curve_entries (ii,curve(ii)%NumRows)
            !% get the first entry of the curve
            curve(ii)%ValueArray(1,:) = interface_get_first_entry_table(ii)
            !% populate the rest of the curves
            do jj = 2,curve(ii)%NumRows
                curve(ii)%ValueArray(jj,:) = interface_get_next_entry_table(ii, API_CURVE)
            end do
        end do

        if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

    end subroutine init_curves
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
        if (icrash) return
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
                    if (.not. node%YN(nidx, nYN_has_dwfInflow)) then !% extInflow only
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
                    BC%headI(ii, bi_face_idx) = node%I(nidx, ni_elemface_idx) !% face idx
                else
                    print *, "Error, BC type can't be a head BC for node " // node%Names(nidx)%str
                    stop
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
        if (icrash) return
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
        !%-----------------------------------------------------------------------------
        !%
        !% Description:
        !% Determines the overall length of the common coarray to handle the different
        !% number of elements on each processor
        !%
        !%-----------------------------------------------------------------------------
        integer :: nimgs_assign
        integer, allocatable :: unique_imagenum(:)
        integer :: ii, jj, kk, idx, counter, elem_counter=0, face_counter=0, junction_counter=0

        integer :: duplicated_face_counter=0
        integer, allocatable :: node_index(:), link_index(:), temp_arr(:)
        character(64) :: subroutine_name = 'init_coarray_length'
        !%-----------------------------------------------------------------------------
        if (icrash) return
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
    subroutine init_time()
        !%-----------------------------------------------------------------------------
        !%
        !% Description:
        !%
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------

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

!%==========================================================================
!% END OF MODULE
!%==========================================================================
!%
end module initialization

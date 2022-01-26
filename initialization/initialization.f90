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
    use utility
    use utility_allocate
    use utility_array
    use utility_datetime
    use utility_output
    use utility_profiler
    use utility_files
    use output
    use pack_mask_arrays

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
        !%-------------------------------------------------------------------
        !% Description:
        !% Calls all the initialization subroutines
        !%-------------------------------------------------------------------
        !% Declarations
            integer :: ii,jj
            integer, pointer :: Npack, thisP(:)
            integer, allocatable :: tempP(:)
            character(64) :: subroutine_name = 'initialize_toplevel'
            !% temporary debugging
            integer    :: elemInLink(100), nEleminLink, iset(5)
            integer    :: thislink_idx, thislink_image
            integer    :: thisnode_idx, thisnode_image, elemJM_idx
            integer    :: iUpSet(max_up_branch_per_node,4)
            integer    :: iDnSet(max_dn_branch_per_node,4)
            integer    :: nUpBranch, nDnBranch
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (icrash) return
            if (setting%Debug%File%initialization) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------    
        !% --- set the CPU and wall-clock timers
        call init_model_timer()

        !% --- define the reverse keys (used mainly for debugging)
        call define_keys_reverse()
        !% NOTES:
        !%    reverseKey(ii) gives you the text name of the ii key
        !%    there are also two useful subroutines:
        !%       call define_keys_printByNumber() !% command-line writes a full list of keys by number
        !%       call define_keys_printByName()   !% command-line writes a full list of keys in alphabetical order

        !% --- assign and store unit numbers for input/output files
        call util_file_assign_unitnumber ()

        !% --- get command line assignments and store
        call util_file_get_commandline ()

        !% --- setup the input project paths and filenames from command line arguments.
        !%        Note that all files and folders must exist or you get error condition.
        !%        This is needed here so that -p command line option works
        call util_file_setup_input_paths_and_files()
        
        !% --- load the settings.json file with the default setting% model control structure
        !%         define_settings_load is one of the few subroutines in the Definition modules
        !%         If the file is not found, the defaults in define_settings.f90 are used 
        call define_settings_load()
        
        !% --- if the settings.json file was read we need to re-process the command-line 
        !%        options a second time to prevent overwrite from json file.
        !%        That is, settings on the command line take precedence over the json file
        if (setting%JSON_FoundFileYN) then
            call util_file_assign_unitnumber ()
            call util_file_get_commandline ()
            call util_file_setup_input_paths_and_files()
        end if

        !% --- initialize the time stamp used for output (must be after json is read)
        call init_timestamp ()

        !% --- setup the output file directories. 
        !%        This will create a new directory with a timestamp for output
        call util_file_setup_output_folders()

        sync all

        !% -- program header
        if ((setting%Output%Verbose) .and. (this_image() == 1)) &
             call util_print_programheader ()    

        !% --- set up the profiler
        if (setting%Profile%useYN) then
            call util_allocate_profiler ()
            call util_profiler_start (pfc_initialize_all)
        else 
            !% continue without profiler    
        end if

        !% --- initialize the coarrays that depend on number of images
        !%     and not on number of links/nodes, elements or faces.
        !if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, "begin initialize secondary coarrays"
        call util_allocate_secondary_coarrays ()

        !% --- initialize the API with the SWMM-C code
        !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin interface between SWMM-C and 5+"
        call interface_init ()

        !% --- set up and store the SWMM-C link-node arrays in equivalent Fortran arrays
        !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin link-node processing"
        call init_linknode_arrays ()
        
        !% --- initialize globals that are run-time dependent
        !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin initialize globals"
        call init_globals()

        !% --- allocate storage for subcatchment arrays
        if (setting%Simulation%useHydrology) then 
            !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin subcatchment allocation"
            call util_allocate_subcatch()
        else    
            !% continue without hydrology    
        end if

        !% --- store the SWMM-C curves in equivalent Fortran arrays
        !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin SWMM5 curve processing"
        call init_curves()

        !% --- break the link-node system into partitions for multi-processor operation
        !if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, "begin link-node partitioning"
        call init_partitioning()

        sync all 

        !% --- default keys  brh20220103
        elemI(:,ei_elementType) = undefinedKey
        elemI(:,ei_geometryType) = undefinedKey
        elemI(:,ei_HeqType) = undefinedKey
        elemI(:,ei_QeqType) = undefinedKey
        elemI(:,ei_specificType) = undefinedKey
        
        !% --- translate the link-node system into a finite-volume network
        if (setting%Simulation%useHydraulics) then 
            !if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, "begin network definition"
            call network_define_toplevel ()
        else 
            if (this_image() == 1) then
                write(*,*) 'USER ERROR: setting.Simulation.useHydraulics == .false.'
                write(*,*) '...this presently is not supported in SWMM5+'
            end if
            stop 879815
        end if                                   

        sync all 

        !% --- iinitialize boundary and ghost elem arrays for inter image data transfer
        call init_boundary_ghost_elem_array ()

        sync all

        !% --- initialize the time variables
        !if (setting%Output%Verbose) print *, "begin initializing time"
        call init_time()

        !% --- initialize the subcatchments connecting to SWMM-C
        if (setting%Simulation%useHydrology) then 
            if (SWMM_N_subcatch > 0) then
                if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin subcatchment initialization"
                call init_subcatchment()
            else 
                if (this_image() == 1) then
                    write(*,'(A)') 'setting.Simulation.useHydrology requested, but no subcatchments found.'
                    write(*,'(A)') '...skipping hydrology in this simulation.'
                end if
                setting%Simulation%useHydrology = .false.
            end if
        else 
            !% continue without hydrology    
        end if

        sync all

        !% --- initialize the output reports
        !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin initializing output report"
        call init_report()

        !% --- set up initial conditions in the FV network
        if (setting%Simulation%useHydraulics) then !% brh20211208 -- only if N_link > 0
            if ((setting%Output%Verbose) .and. (this_image() == 1)) then 
                !print *, "begin initial conditions"
                if (this_image() == 1) then
                    if ((N_link > 5000) .or. (N_node > 5000)) then
                        write(*,"(A)") " ... setting initial conditions -- may take several minutes for big systems ..."
                        write(*,"(A,i8,A,i8,A)") "      SWMM system has ", SWMM_N_link, " links and ", SWMM_N_node, " nodes"
                        write(*,"(A,i8,A)")      "      FV system has   ", sum(N_elem(:)), " elements"
                    else 
                        !% no need to write for small systems
                    end if
                else 
                    !% only print for 1 processor
                end if
            else 
                !% be silent    
            end if    
            !% --- initial conditions all are setup here
            call init_IC_toplevel ()
        else 
            if (this_image() == 1) then
                write(*,*) 'USER ERROR: setting.Simulation.useHydraulics == .false.'
                write(*,*) '...this presently is not supported in SWMM5+'
            end if
            stop 9378975
        end if          

        !% initialize volume conservation storage for debugging
        elemR(:,er_VolumeConservation) = zeroR    

        !% --- setup the multi-level finite-volume output
        !%        Ideally, this should be a procedure accessed in the output module, 
        !%        but that caused linking problems due to pack/mask calls
        if (setting%Output%Report%provideYN) then 
            if (setting%Simulation%useHydraulics) then !% 
                !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin setup of output files"
                !% --- Get the output element and face locations
                !print *, 'starting outputML_selection'
                call outputML_selection ()                
                !% --- Create packed arrays of elem row numbers that are output
                !print *, 'starting pack_element_outputML'
                call pack_element_outputML ()
                !% --- Create packed arrays of face row numbers that are output
                !print *, 'starting pack_face_outputML'
                call pack_face_outputML ()
                !% --  Setup of the outputML data structures
                !print *, 'starting outputML_setup', this_image()
                call outputML_setup ()
            else 
                if (this_image() == 1) then
                    write(*,*) 'USER ERROR: setting.Simulation.useHydraulics == .false.'
                    write(*,*) '...this presently is not supported in SWMM5+'
                end if
                stop 487587  
            end if  
        else 
            !% continue without any output files                                      
        end if
        !% --- wait for all processors before exiting to the time loop
        sync all
    
        !%------------------------------------------------------------------- 
        !% Closing
            call init_check_setup_conditions()
            call init_timer_stop ()
            if (setting%Simulation%stopAfterInitializationYN) then
                if (this_image() == 1) then
                    write(*,*) ' '
                    write(*,*) '********************************************************'
                    write(*,*) '** Stopping after initialization for review due to -R **'
                    write(*,*) '** as command-line argument or due to setting the     **' 
                    write(*,*) '** stopAfterInitializationYN = true in json file.     **'
                    write(*,*) '** Remove the -R from the command line and/or change  **'
                    write(*,*) '** the json file to run a full simulation.            **'
                    write(*,*) '********************************************************'
                end if
                stop 333578
            end if
            if (icrash) then  !% if crash in initialization, write the output and exit
                if (setting%Output%Report%provideYN) then 
                    call outputML_store_data (.true.)
                end if
                return
            end if
            if (setting%Debug%File%initialization)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine initialize_toplevel
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine init_model_timer()
        !%------------------------------------------------------------------
        !% Description:
        !% starts and stores the CPU clock and wall clock time
        !%------------------------------------------------------------------
        !% Declarations:
        integer(kind=8) :: crate, cmax, cval
        !%-------------------------------------------------------------------
        !% store CPU clock start time
        call cpu_time(setting%Time%CPU%EpochStartSeconds)

        !% get the Wall clock time
        if (this_image() == 1) then
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%CountRate = crate
            setting%Time%WallClock%Start = cval
        end if

        


    end subroutine init_model_timer
!%
!%==========================================================================
!%==========================================================================
!   
    subroutine init_timer_stop ()

        integer(kind=8) :: crate, cmax, cval
        
        if (this_image() == 1) then
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%InitializationEnd = cval
        end if

    end subroutine init_timer_stop
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
        !% Declarations
            integer       :: ii, total_n_links
            character(64) :: subroutine_name = 'init_linknode_arrays'
        !%-----------------------------------------------------------------------------
        !% Preliminaries
            if (icrash) return
            if (setting%Debug%File%initialization) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

            if (.not. api_is_initialized) then
                print *, "ERROR: API is not initialized"
                stop
            end if
        !%-----------------------------------------------------------------------------
        !% Allocate storage for link & node tables
        call util_allocate_linknode()

        link%I(:,li_num_phantom_links) = 0
        node%I(:,ni_N_link_u) = 0
        node%I(:,ni_N_link_d) = 0

        do ii = 1, SWMM_N_link
            link%I(ii,li_idx) = ii
            link%I(ii,li_link_type) = interface_get_linkf_attribute(ii, api_linkf_type)
            link%I(ii,li_weir_type) = interface_get_linkf_attribute(ii, api_linkf_weir_type)
            link%I(ii,li_orif_type) = interface_get_linkf_attribute(ii, api_linkf_orifice_type)
            link%I(ii,li_outlet_type) = interface_get_linkf_attribute(ii, api_linkf_outlet_type)
            link%I(ii,li_pump_type) = interface_get_linkf_attribute(ii, api_linkf_pump_type)
            link%I(ii,li_geometry) = interface_get_linkf_attribute(ii, api_linkf_geometry)
            link%I(ii,li_Mnode_u) = interface_get_linkf_attribute(ii, api_linkf_node1) + 1 ! node1 in C starts from 0
            link%I(ii,li_Mnode_d) = interface_get_linkf_attribute(ii, api_linkf_node2) + 1 ! node2 in C starts from 0
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
            link%R(ii,lr_Length) = interface_get_linkf_attribute(ii, api_linkf_conduit_length)
            !% link%R(ii,lr_TopWidth): defined in network_define.f08
            link%R(ii,lr_BreadthScale) = interface_get_linkf_attribute(ii, api_linkf_xsect_wMax)
            !% link%R(ii,lr_Slope): defined in network_define.f08
            link%R(ii,lr_LeftSlope) = interface_get_linkf_attribute(ii, api_linkf_left_slope)
            link%R(ii,lr_RightSlope) = interface_get_linkf_attribute(ii, api_linkf_right_slope)
            link%R(ii,lr_Roughness) = interface_get_linkf_attribute(ii, api_linkf_conduit_roughness)
            link%R(ii,lr_InitialFlowrate) = interface_get_linkf_attribute(ii, api_linkf_q0)
            !write(*,*) 'api_nodef_initDepth 1'
            link%R(ii,lr_InitialUpstreamDepth) = interface_get_nodef_attribute(link%I(ii,li_Mnode_u), api_nodef_initDepth)
            !write(*,*) 'api_nodef_initDepth 2'
            link%R(ii,lr_InitialDnstreamDepth) = interface_get_nodef_attribute(link%I(ii,li_Mnode_d), api_nodef_initDepth)
            
            link%R(ii,lr_InitialDepth) = (link%R(ii,lr_InitialDnstreamDepth) + link%R(ii,lr_InitialUpstreamDepth)) / 2.0
            link%R(ii,lr_FullDepth) = interface_get_linkf_attribute(ii, api_linkf_xsect_yFull)
            link%R(ii,lr_InletOffset) = interface_get_linkf_attribute(ii,api_linkf_offset1)
            link%R(ii,lr_OutletOffset) = interface_get_linkf_attribute(ii,api_linkf_offset2)

            !% special element attributes
            link%I(ii,li_weir_EndContrations) = interface_get_linkf_attribute(ii, api_linkf_weir_end_contractions)
            link%I(ii,li_curve_id) = interface_get_linkf_attribute(ii, api_linkf_curveid)
            link%R(ii,lr_DischargeCoeff1) = interface_get_linkf_attribute(ii, api_linkf_discharge_coeff1)
            link%R(ii,lr_DischargeCoeff2) = interface_get_linkf_attribute(ii, api_linkf_discharge_coeff2)
            link%R(ii,lr_SideSlope) = interface_get_linkf_attribute(ii, api_linkf_weir_side_slope)
            !% SWMM5 doesnot distinguish between channel and conduit
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
            !% brh20211207s
            link%YN(ii,lYN_isOutput) = (interface_get_linkf_attribute(ii,api_linkf_rptFlag) == 1)
            !% brh20211207e
        end do

        !% check for errors in number of connections
        ! print *, ' '
        ! print *, '=========================================================='
        ! print *, '  in ',trim(subroutine_name)
        ! print *, 'printing number of connections up and down on each node'
        do ii = 1,N_node

            !print *, ii, node%I(ii, ni_N_link_u),  node%I(ii, ni_N_link_d)

            if (node%I(ii, ni_N_link_u) > max_up_branch_per_node) then
                if (this_image() == 1) then
                    write(*,*) 'FATAL ERROR IN INPUT FILE'
                    write(*,"(A,i4,A)") 'One or more nodes have more than ',max_up_branch_per_node,' upstream connections'
                    write(*,*) 'Unfortunately, this connection limit is a hard-coded limit of SWMM5+ an cannot be exceeded.'
                end if
                stop 387666
            end if
            if (node%I(ii, ni_N_link_u) > max_dn_branch_per_node) then
                if (this_image() == 1) then
                    write(*,*) 'FATAL ERROR IN INPUT FILE'
                    write(*,"(A,i4,A)") 'One or more nodes have more than ',max_dn_branch_per_node,' downstream connections'
                    write(*,*) 'Unfortunately, this connection limit is a hard-coded limit of SWMM5+ an cannot be exceeded.'
                end if
                stop 86752
            end if

        end do

        !stop 398706
        !write(*,*) 
        !write(*,*) 'FINISHED WITH LINKS ---------------------------------------------------------'
        !write(*,*) 'N_node = ',N_node
        !write(*,*)

        
        do ii = 1, N_node
            !write(*,*) '======= starting node ',ii
            !write(*,*)
            total_n_links = node%I(ii,ni_N_link_u) + node%I(ii,ni_N_link_d)
            node%I(ii, ni_idx) = ii
            !%
            !% --- handle special types of nodes
            if (interface_get_nodef_attribute(ii, api_nodef_type) == API_OUTFALL) then
                !write(*,*) '... is outfall type'
                node%I(ii, ni_node_type) = nBCdn
            else if (interface_get_nodef_attribute(ii, api_nodef_type) == API_STORAGE) then
                !write(*,*) '... is storage type'
                node%I(ii, ni_node_type) = nJm
                node%YN(ii, nYN_has_storage) = .true.
            else 
                !% --- classify by number of links connected
                select case (total_n_links)
                    case (oneI)
                        !write(*,*) '... is 1 junction is an upstream BC
                        node%I(ii, ni_node_type) = nJ1
                    case (twoI)
                        !write(*,*) '... is 2 junction type'        
                        node%I(ii, ni_node_type) = nJ2
                    case default 
                        !write(*,*) '... is 3+ junction type'
                        node%I(ii, ni_node_type) = nJm
                end select
            end if 
        
            ! if (interface_get_nodef_attribute(ii, api_nodef_type) == API_OUTFALL) then
            !     !write(*,*) '... is outfall type'
            !     node%I(ii, ni_node_type) = nBCdn
            ! else if (interface_get_nodef_attribute(ii, api_nodef_type) == API_STORAGE) then
            !     !write(*,*) '... is storage type'
            !     node%I(ii, ni_node_type) = nJm
            !     node%YN(ii, nYN_has_storage) = .true.
            ! else if ((total_n_links == twoI)          .and. &
            !          (node%I(ii,ni_N_link_u) == oneI) .and. &
            !          (node%I(ii,ni_N_link_d) == oneI) )then
            !     !write(*,*) '... is 2 junction type'        
            !     node%I(ii, ni_node_type) = nJ2
            ! else if (total_n_links >= twoI) then
            !     !write(*,*) '... is 3+ junction type'
            !     node%I(ii, ni_node_type) = nJm
            ! else if (total_n_links == oneI) then  !% brh 20211217
            !     !write(*,*) '... is 1 junction is an upstream BC
            !     node%I(ii, ni_node_type) = nJ1
            ! else 
            !     write(*,*) 'CODE or INP FILE ERROR, unexpected else condition '
            !     write(*,*) 'Node type is undefined for node',ii
            !     stop 98075               
            ! end if
            !write(*,*)

            !write(*,*) 'call api_nodef_has_extInflow'
            node%YN(ii, nYN_has_extInflow) = (interface_get_nodef_attribute(ii, api_nodef_has_extInflow) == 1)
            !write(*,*) '... nYN_has_extInflow = ',node%YN(ii, nYN_has_extInflow)
            !write(*,*)

            !write(*,*) 'call api_nodef_has_dwfInflow'
            node%YN(ii, nYN_has_dwfInflow) = (interface_get_nodef_attribute(ii, api_nodef_has_dwfInflow) == 1)
            !write(*,*) '... nYN_has_dwfInflow =', node%YN(ii,nYN_has_dwfInflow)
            !write(*,*)

            !% --- set up the inflows
            if (node%YN(ii, nYN_has_extInflow) .or. node%YN(ii, nYN_has_dwfInflow)) then
                !% set inflow to true for any node type
                node%YN(ii, nYN_has_inflow) = .true.
                !% change the node type of an nJ1 with inflow to nBCup
                if (node%I(ii, ni_node_type) == nJ1) node%I(ii, ni_node_type) = nBCup
                !if ((node%I(ii,ni_N_link_u) == zeroI) .and. (total_n_links == oneI)) then
                !    node%I(ii, ni_node_type) = nBCup
                !end if
            end if

            !write(*,*) 'call api_nodef_initDepth'
            node%R(ii,nr_InitialDepth)      = interface_get_nodef_attribute(ii, api_nodef_initDepth)
            !write(*,*) '... nr_InitialDepth = ',node%R(ii,nr_InitialDepth)
            !write(*,*)

            !write(*,*) 'call api_nodef_invertElev'
            node%R(ii,nr_Zbottom)           = interface_get_nodef_attribute(ii, api_nodef_invertElev)
            !write(*,*) '... nr_Zbottom = ', node%R(ii,nr_Zbottom) 
            !write(*,*)

            !write(*,*) 'call api_nodef_fullDepth -- may be zero!'
            node%R(ii,nr_FullDepth)         = interface_get_nodef_attribute(ii, api_nodef_fullDepth)
            !write(*,*) '... nr_FullDepth = ',node%R(ii,nr_FullDepth) 
            !write(*,*)

            !write(*,*) 'call api_nodef_StorageConstant'
            node%R(ii,nr_StorageConstant)   = interface_get_nodef_attribute(ii, api_nodef_StorageConstant)
            !write(*,*) '... nr_StorageConstant = ',node%R(ii,nr_StorageConstant)
            !write(*,*)

            !write(*,*) 'call api_nodef_StorageCoeff'
            node%R(ii,nr_StorageCoeff)      = interface_get_nodef_attribute(ii, api_nodef_StorageCoeff)
            !write(*,*) '... nr_StorageCoeff = ',node%R(ii,nr_StorageCoeff)
            !write(*,*)
            
            !write(*,*) 'call api_nodef_StorageExponent'
            node%R(ii,nr_StorageExponent)   = interface_get_nodef_attribute(ii, api_nodef_StorageExponent)
            !write(*,*) '... nr_StorageExponent = ',node%R(ii,nr_StorageExponent)
            !write(*,*)

            !write(*,*) 'call api_nodef_StorageCurveID'
            node%I(ii,ni_curve_ID)          = interface_get_nodef_attribute(ii, api_nodef_StorageCurveID)
            !write(*,*) '... ni_curve_ID = ',node%I(ii,ni_curve_ID)
            !write(*,*)

            !write(*,*) 'call interface_get_BC_resolution'
            node%I(ii,ni_pattern_resolution) = interface_get_BC_resolution(ii)
            !write(*,*) '... ni_pattern_resolution = ',node%I(ii,ni_pattern_resolution)
            !write(*,*)  

            !% brh20211207s
            !write(*,*) 'call api_nodef_rptFlag'
            node%YN(ii,nYN_isOutput)          = (interface_get_nodef_attribute(ii, api_nodef_rptFlag) == 1)
            !write(*,*) '... nYN_isOutput = ',node%YN(ii,nYN_isOutput)
            !write(*,*)
            !% brh20211207e
        end do

        !% Update Link/Node names
        call interface_update_linknode_names()

        !%-----------------------------------------------------------------------------
        !% closing
            ! print *,  'at end of ',trim(subroutine_name)
            ! print *, 'node data'
            
            ! print *, 'idx,    nodeType,    linkU,    linkD,   curveID, patternRes,    assigned'
            ! do ii=1,N_node
            !     write(*,"(10i8)") node%I(ii,ni_idx), node%I(ii,ni_node_type), node%I(ii,ni_N_link_u), node%I(ii,ni_N_link_d) &
            !     , node%I(ii,ni_curve_ID), node%I(ii,ni_pattern_resolution), node%I(ii,ni_assigned)
            ! end do

            if (setting%Debug%File%initialization)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_linknode_arrays
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_globals()
        !%------------------------------------------------------------------
        !% Description:
        !% Initializes globals values that are run-time based, i.e.
        !% they depend on a size
        !%------------------------------------------------------------------
        !% branchsign global is used for junction branches (JB)
        !%     for upstream (+1) and downstream (-1)
        branchsign(1:max_branch_per_node-1:2) = +oneR
        branchsign(2:max_branch_per_node:2)   = -oneR
        !%------------------------------------------------------------------
        !% Closing
    end subroutine init_globals
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
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if (.not. api_is_initialized) then
            print *, "ERROR: API is not initialized"
            stop 84785
        end if

        !% we create additional curves for functional storage as well
        !% this allocates the space for functional storage curve
        additional_storage_curves = count((node%YN(:,nYN_has_storage)) .and. &
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
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_curves
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
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% brh20211208s
        !% if there are no links, the system cannot be partitioned
        if (N_link == 0) then
            if (this_image() == 1) then
                write(*,*) '******************************************************'
                write(*,*) '*                    FATAL ERROR                     *'
                write(*,*) '* The SWMM input file does not include any links.    *'
                write(*,*) '* The SWMM5+ code requires at least one link to run. *'
                write(*,*) '* This run was stopped without any output.           *'
                write(*,*) '******************************************************'
            end if
            stop 9705322
            !write(*,*) '*** WARNING no conduit/channel links found, using SWMM hydrology only'
            setting%Simulation%useHydraulics = .false.
            return
        end if
        !% brh20211208e    

        if (setting%Profile%useYN) call util_profiler_start (pfc_init_partitioning)

        !% find the number of elements in a link based on nominal element length
        do ii = 1, SWMM_N_link
            call init_discretization_nominal(ii)
        end do

        !% Set the network partitioning method used for multi-processor parallel computation
        call init_partitioning_method()
        sync all

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

        if (setting%Profile%useYN) call util_profiler_stop (pfc_init_partitioning)

        if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_partitioning
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_subcatchment()
        !%------------------------------------------------------------------
        !% Description:
        !% sets the connection between the subcatchments of SWMM-c and the
        !% elements that will get runoff
        !%
        !% This populates subcatchI() for the node index and element index
        !% that are connected to the subcatchment. Note that the node index
        !% stored is the SWMM5+ node index and the SWMM-C index is obtained
        !% by subctracting one.
        !%
        !% HACK -- not sure how this functions of a node appears on more
        !% than one image
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: ii
            integer, pointer :: nodeIdx(:), elemIdx(:), nodeType(:)
            integer, pointer :: tface
            logical, pointer :: isToNode(:)
            character(64) :: subroutine_name = 'init_subcatchment'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (icrash) return
            if (setting%Debug%File%initialization) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
            nodeIdx  => subcatchI(:,si_runoff_nodeIdx) 
            elemIdx  => subcatchI(:,si_runoff_elemIdx)
            isToNode => subcatchYN(:,sYN_hasRunoff)
            nodeType => node%I(:,ni_node_type)
        !%------------------------------------------------------------------

        !% set the counter for the number of subcatchments to an element to zero
        !% THIS IS ONLY NEEDED IF WE NEED TO GO FROM elem -> subcatch, which we should avoid
        !elemI(:, ei_Nsubcatch) = zeroI

        !% cycle through subcatchments to set connections to runoff nodes in SWMM-C
        do ii=1,SWMM_N_subcatch
            !%Add 1 to the SWMM-C node to get the SWMM5+ node
            nodeIdx(ii) = interface_get_subcatch_runoff_nodeIdx(ii-1)+oneI
            if (nodeIdx(ii) < 1) then !% not a runoff node (SWMM-C flag)
                isToNode(ii) = .false.
            else
                isToNode(ii) = .true.
            end if         

            !% only handle the subcatchment on images with its connected node
            if (this_image() .eq. node%I(nodeIdx(ii), ni_P_image)) then
                subcatchI(ii,si_runoff_P_image) = this_image()
                select case (nodeType(nodeIdx(ii)))
                    case (nJ2)
                        !% for a node that is a face, the subcatch connects to the element
                        !% upstream of the face, which is defined by ni_elemface_idx
                        elemIdx(ii) = node%I(nodeIdx(ii), ni_elemface_idx)
                    case (nJm)
                        !% for a node that is a multi-branch junction, subcatch connects to 
                        !% the element itself, which is defined by ni_elemface_idx
                        elemIdx(ii) = node%I(nodeIdx(ii), ni_elemface_idx)
                    case (nBCup,nJ1)
                        !% for a node that is an upstream BC or dead end the subcatch connects 
                        !% into the first element downstream of the face
                        !% Here ni_elemface_idx holds the face index
                        tface => node%I(nodeIdx(ii),ni_elemface_idx) 
                        if (tface .ne. nullvalueI) then 
                            elemIdx(ii) = faceI(tface,fi_Melem_dL)
                        else
                            elemIdx(ii) = nullvalueI
                        end if
                    case (nBCdn)
                        !% for a node that is an downstreamstream BC, the subcatch connects 
                        !% first element upstreamstream of the face
                        !% into the Here ni_elemface_idx holds the face index
                        tface => node%I(nodeIdx(ii),ni_elemface_idx)
                        if (tface .ne. nullvalueI) then 
                            elemIdx(ii) = faceI(tface,fi_Melem_uL)
                        else
                            elemIdx(ii) = nullvalueI
                        end if
                    case default 
                        write(*,*) 'ERROR CODE: unexpected case default in '//trim(subroutine_name)
                end select
                !% store logical for elem
                if (elemIdx(ii) .ne. nullvalueI) then 
                    elemYN(elemIdx(ii), eYN_hasSubcatchRunOff) = .true.
                    !elemI(elemIdx(ii), ei_Nsubcatch) = elemI(elemIdx(ii), ei_Nsubcatch) + oneI
                end if
            end if
        end do

        ! do ii = 1,SWMM_N_subcatch
        !     print *, ii
        !     print *,  ii, subcatchI(:,si_runoff_nodeIdx) ,  subcatchI(:,si_runoff_elemIdx) 
        ! end do

        !do ii = 1,size(node%I,DIM=1)
        !    print *, node%I(ii,ni_idx), node%I(ii,ni_node_type), reverseKey(node%I(ii,ni_node_type))
        !end do

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%initialization)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_subcatchment
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_time()
        !%------------------------------------------------------------------
        !% Description:
        !% initializes the time either using the SWMM input file or the
        !% json file data (selected by setting.Time.useSWMMinpYN)
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        !% HACK start time is always measured from zero (need to fix for hotstart)
        setting%Time%Start = zeroR 
        setting%Time%Now   = zeroR
        setting%Time%Step  = zeroI

        if (setting%Time%useSWMMinpYN) then 
            !% set the start/stop times and time steps from SWMM *.inp file
            setting%Time%StartEpoch    = setting%SWMMinput%StartEpoch
            setting%Time%EndEpoch      = setting%SWMMinput%EndEpoch
            setting%Time%Hydraulics%Dt = setting%SWMMinput%RouteStep
            setting%Time%Hydrology%Dt  = setting%SWMMinput%WetStep
            ! HACK ??                  = setting%SWMMinput%DryStep
            ! HACK ??                  = setting%SWMMinput%TotalDuration
        else 
            !% use values from json file
        end if

        !% Translate epoc endtime to seconds from a zero start time
        !% use floor() to match approachin SWMM-C
        setting%Time%End = real(floor(                            &
                (setting%Time%EndEpoch - setting%Time%StartEpoch) &
                 * real(secsperday)),KIND=8)

        !% null out the wet step if not using hydrology
        if (.not. setting%Simulation%useHydrology) setting%Time%Hydrology%Dt = nullValueR

        !%------------------------------------------------------------------
        !% Closing
    end subroutine init_time
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine init_report()

        !% if setting require the SWMM input file values, then overwrite setting values
        if (setting%Output%Report%useSWMMinpYN) then 
            setting%Output%Report%StartTime    = util_datetime_epoch_to_secs(setting%SWMMinput%ReportStartTimeEpoch)
            setting%Output%Report%TimeInterval = setting%SWMMinput%ReportTimeInterval
            if ((setting%Output%Verbose) .and. (this_image() == 1)) then
                write(*,*) ' '
                write(*,"(A)") ' ... using report start time and time interval from SWMM input file (*.inp)'
                write(*,*) ' '
            end if
        else 
            if ((setting%Output%Verbose) .and. (this_image() == 1)) then
                write(*,*) ' '
                write(*,"(A)") '... using  report start time and time interval from *.json file'
                write(*,*) ' '
            end if
        end if

        !% if selected report time is before the start time
        if (setting%Output%Report%StartTime < setting%Time%Start) then 
            setting%Output%Report%StartTime = setting%Time%Start
        else 
            !% continue
        end if

        if (setting%Output%Report%TimeInterval < zeroR) then 
            if (this_image() == 1) then
                write(*,*) '***************************************************************'
                write(*,*) '** WARNING -- selected report time interval is zero or less, **'
                write(*,*) '**          so all reporting will be suppressed              **'
                write(*,*) '***************************************************************'
            end if
            setting%Output%Report%provideYN = .false.
            setting%Output%Report%suppress_MultiLevel_Output = .true.
            setting%Output%Report%ThisStep = 1
        else 
            !% Initialize report step -- 
            !% Determine how many report steps have already been missed before
            !% the output reports are actually written
            setting%Output%Report%ThisStep = int( &
                        ( setting%Output%Report%StartTime - setting%Time%Start ) &
                        / setting%Output%Report%TimeInterval )
        end if

    end subroutine init_report
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
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

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
                elseif ((node%I(idx, ni_node_type) == nBCup) .or. (node%I(idx, ni_node_type) == nJ1)) then  !% brh20211217
                    face_counter = face_counter +1 !% add the upstream faces
                elseif (node%I(idx, ni_node_type) == nBCdn) then
                    face_counter = face_counter +1 !% add the downstream faces
                elseif (node%I(idx, ni_node_type) == nJm) then   
                    !% skip -- faces are counted elsewhere
                else 
                    print *, jj, node%I(idx, ni_node_type), reverseKey(node%I(idx, ni_node_type))
                    print *, reverseKey(nJ1), reverseKey(nJ2), reverseKey(nBCup), reverseKey(nBCdn)
                    print *, 'CODE ERROR, unexpected else'
                    stop 390715
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
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_coarray_length
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine init_boundary_ghost_elem_array()
        !
        !--------------------------------------------------------------------------
        ! initialize ghost and boundary element arrays for inter image data transfer
        !--------------------------------------------------------------------------
        integer          :: ii, jj, NSfaces, eset_local(4), eBGset(4)
        integer, pointer :: Nfaces, fIdx, fGidx, eUp, eDn, ci, BeUp, BeDn
        integer, dimension(:), allocatable, target :: packed_shared_face_idx
        character(64)    :: subroutine_name = 'init_boundary_ghost_elem_array'
        !--------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% only initialize inter-image data transfer array for ore than one processor
        if (num_images() > 1) then

            !% allocate elemB and elemG data structure
            call util_allocate_boundary_ghost_elem_array()

            !% number of faces in this image
            Nfaces => N_face((this_image()))
            !% count the number of shared faces in this image
            NSfaces = count(faceYN(1:Nfaces,fYN_isSharedFace))
            !% packed indexes of shared faces in this image
            packed_shared_face_idx = pack(faceI(1:Nfaces,fi_Lidx),faceYN(1:Nfaces,fYN_isSharedFace))
            !% column indexes of elemI local data needed to be transferred 
            eset_local = [ei_Lidx, ei_Gidx, ei_Mface_uL, ei_Mface_dL]
            !% column indexes of the spot for receiving local elemI data
            eBGset     = [ebgi_elem_Lidx, ebgi_elem_Gidx, ebgi_Mface_uL, ebgi_Mface_dL]

            print*
            do ii = 1,NSfaces
                fIdx => packed_shared_face_idx(ii) 
                eUp  => faceI(fIdx,fi_Melem_uL)
                eDn  => faceI(fIdx,fi_Melem_dL)

                !% local integer data transfer
                elemB%I(ii,ebgi_idx) = ii

                if (faceYN(fIdx,fYN_isUpGhost)) then
                    !% copy the local integer data over
                    elemB%I(ii,eBGset) = elemI(eDn,eset_local)
                    !% save the position of the boundary array index in the elemI array
                    elemI(eDn,ei_BoundaryArray_idx) = ii
                    !% save the position of the boundary array index in the faceI array
                    faceI(fIdx,fi_BoundaryElem_dL) = ii
                else if (faceYN(fIdx,fYN_isDnGhost)) then
                    !% copy the local integer data over
                    elemB%I(ii,eBGset) = elemI(eUp,eset_local)
                    !% save the position of the boundary array index in the elemI array
                    elemI(eUp,ei_BoundaryArray_idx) = ii
                    !% save the position of the boundary array index in the faceI array
                    faceI(fIdx,fi_BoundaryElem_uL) = ii
                else if (faceYN(fIdx,fYN_isUpGhost) .and. faceYN(fIdx,fYN_isDnGhost)) then
                    print*, 'error: condition not handeled single element with two shared faces'
                    stop 914744
                else
                    print*, 'should not reach this condition'
                    stop 54673
                end if
            end do

            !% wait and sync untill all the images update their elemB arrays with local data
            sync all

            !% now loop through all the shared faces again to find the boundary array location of the ghost 
            !% element in a remote image
            do ii = 1,NSfaces
                fIdx  => packed_shared_face_idx(ii)
                fGidx => faceI(fIdx,fi_Gidx)
                eUp   => faceI(fIdx,fi_Melem_uL)
                eDn   => faceI(fIdx,fi_Melem_dL)
                ci    => faceI(fIdx,fi_Connected_image)
                BeUp  => faceI(fIdx,fi_BoundaryElem_uL)
                BeDn  => faceI(fIdx,fi_BoundaryElem_dL)
                do jj = 1,N_face(ci)
                    if ((faceI(jj,fi_Connected_image)[ci] == this_image()) .and. &
                        (faceI(jj,fi_Gidx)[ci] == fGidx)) then
                        !% find the local ghost element index of the connected image
                        if (faceYN(jj,fYN_isUpGhost)[ci]) then
                            faceI(jj,fi_BoundaryElem_uL)[ci] = BeUp
                        else if (faceYN(jj,fYN_isDnGhost)[ci]) then
                            faceI(jj,fi_BoundaryElem_dL)[ci] = BeDn
                        end if
                    end if
                end do
            end do

            !% sync all the processors
            sync all

            deallocate(packed_shared_face_idx) !% deallocate temporary arrays
        end if

        if (setting%Debug%File%network_define) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_boundary_ghost_elem_array
!%    
!%==========================================================================
!%==========================================================================
!%
    subroutine init_check_setup_conditions ()
        !%------------------------------------------------------------------
        !% Description:
        !% provides warnings of SWMM5+ settings for features that are
        !% still in development and/or testing.
        !%------------------------------------------------------------------
            logical :: ifound = .false.
        !%------------------------------------------------------------------  
        !% Preliminaries: 
            if (this_image() .ne. 1) return
        !%------------------------------------------------------------------ 
        write(*,*) ' '
        write(*,'(A)') '*******************************************************************'
        
        if (.not. setting%ZeroValue%UseZeroValues) then
            write(*,'(A)') '** setting.ZeroValue.UseZeroValues = false, which can cause infinity'
            write(*,'(A)') '** when zero depths encountered. Strongly recommend using true.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%VariableDT%ApplyYN) then
            write(*,'(A)') '** setting.VariableDT%ApplyYN = false, which is not fully tested'
            write(*,'(A)') '** and is likely to cause unknown problems. '
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%Time%useSWMMinpYN) then
            write(*,'(A)') '** setting.Time.useSWMMinYN = false, which means the simulation'
            write(*,'(A)') '** will use the timing data from the json file rather than from'
            write(*,'(A)') '** the SWMM *.inp file. Use this option with caution.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (setting%TestCase%UseTestCaseYN) then
            write(*,'(A)') '** setting.TestCase.UseTestCaseYN = true, which means the simulation'
            write(*,'(A)') '** will entirely ignore the SWMM input file and the json file.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%Solver%PreissmannSlot) then
            write(*,'(A)') '** setting.Solver.PreissmannSlot = false, which has not been fully tested.'
            write(*,'(A)') '** Unknown problems will occur in surcharged conditions. '
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (setting%Solver%SubtractReferenceHead) then
            write(*,'(A)') '** setting.Solver.SubtractReferenceHead = true, which is not fully tested with hydrology'
            write(*,'(A)') '** if there are problems with head values, return setting to false. '
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%Simulation%useHydraulics) then
            write(*,'(A)') '** setting.Simulation.useHydraulics = false, which implies the user'
            write(*,'(A)') '** wants a SWMM simulation without SWMM5+ hydraulics. This option '
            write(*,'(A)') '** is incomplete.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (setting%Simulation%useHydrology) then
            write(*,'(A)') '** setting.Simulation.useHydrology = true. '
            write(*,'(A)') '** the interface between SWMM/SWMM5+ for hydrology is still in '
            write(*,'(A)') '** development and the user may experience bugs. Please report'
            write(*,'(A)') '** these to the code development team'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (setting%Profile%useYN) then
            write(*,'(A)') '** setting.Profile.useYN = true. '
            write(*,'(A)') '** This code profiler is still in development, so results should '
            write(*,'(A)') '** be used with caution.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%Output%Verbose) then
            write(*,'(A)') '** setting.Output.Verbose = false. '
            write(*,'(A)') '** While SWMM5+ is in development we recommend using true. '
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%Output%Warning) then
            write(*,'(A)') '** setting.Output.Warning = false. '
            write(*,'(A)') '** While SWMM5+ is in development we recommend using true. '
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (setting%Junction%isDynamicYN) then
            write(*,'(A)') '** setting.Junction.isDynamic = true. We recommend false. '
            write(*,'(A)') '** Dynamic junctions are in development and provide inconsistent results. '
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%File%force_folder_creationYN) then
            write(*,'(A)') '** setting.File.force_folder_creationYN = false. We recommend true. '
            write(*,'(A)') '** If this is false, SWMM5+ will have write errors if output folders '
            write(*,'(A)') '** do not exist. This can cause the code to crash without warning.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%File%UseCommandLineFoldersYN) then
            write(*,'(A)') '** setting.File.UseCommandLineFoldersYN = false. We recommend true. '
            write(*,'(A)') '** Other approaches for foldres have not been fully tested.  '
            write(*,'(A)') '** We are not sure what will happen on your computer.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (setting%BC%disableInterpolationYN) then
            write(*,'(A)') '** setting.BC.disableInterpolationYN = true. We recommend false. '
            write(*,'(A)') '** This disables interpolation of time series and simply takes the'
            write(*,'(A)') '** next highest value.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%Output%Report%useSWMMinpYN) then
            write(*,'(A)') '** setting.Output.Report.useSWMMinYN = false, which means the simulation'
            write(*,'(A)') '** will use the report data from the json file rather than from'
            write(*,'(A)') '** the SWMM *.inp file. Use this option with caution.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%Output%Report%provideYN) then
            write(*,'(A)') '** setting.Output.Report.provideYN = false, which means the simulation'
            write(*,'(A)') '** will suppress the output reporting. Use this option with caution.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (setting%Output%Report%suppress_MultiLevel_Output) then
            write(*,'(A)') '** setting.Report.suppress_MultiLevel_Output = true, which means '
            write(*,'(A)') '** the simulation will suppress all the output reporting from hydraulics.'
            write(*,'(A)') '** Presently, the only output from hydraulics is through the ML scheme.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%Limiter%Velocity%UseLimitMaxYN) then
            write(*,'(A)') '** setting.Limiter.Velocity.UseLimitMaxYN = false. We recommend true. '
            write(*,'(A)') '** Without a velocity limiter, a minor instability can become'
            write(*,'(A)') '** catastrophic.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%Limiter%Dt%UseLimitMinYN) then
            write(*,'(A)') '** setting.Limiter.Dt.UseLimitMinYN = false. We recommend true. '
            write(*,'(A)') '** Set to false the variable dt could become very small as the '
            write(*,'(A)') '** when the simulation has gone unstable, resulting in the code '
            write(*,'(A)') '** appearing to hang up.'
            write(*,'(A)') '** '
            ifound = .true.
        end if

        if (ifound) then
            write(*,'(A)') '**                          WARNING'
            write(*,'(A)') '** The above are possible issues with user settings for this version of SWMM5+.  '
            write(*,'(A)') '** '
            write(*,'(A)') '*******************************************************************'
        end if
        !%------------------------------------------------------------------
    end subroutine init_check_setup_conditions    
!%    
!%==========================================================================
!% END OF MODULE
!%==========================================================================
!%
end module initialization

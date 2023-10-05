module initialization
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Provides the top-level initialization of SWMM5+
    !%
    !% Method:
    !% Creates and populates the data storage used for SWMM5+. This includes 
    !% typed data structures for the link/node data from EPA SWMM and 
    !% 2D arrays of data for SWMM5+ hydraulics. Interface calls are used
    !% to access EPA SWMM and read data from the *.inp file. The toplevel
    !% subroutine calls the partitioning algorithm to subdivide the link/node
    !% network for parallel processing. Each processor is assigned a portion
    !% of the network and has 2D arrays of local element/face data allocated.
    !%==========================================================================

    use IFPORT
    use boundary_conditions
    use define_keys
    use define_globals
    use define_settings
    use define_indexes
    use discretization
    use initial_condition
    use interface_
    use network_define
    use partitioning
    use culvert_elements, only: culvert_parameter_values
    use utility
    use utility_allocate
    use utility_array
    use utility_datetime
    use utility_output
    use utility_profiler
    use utility_files
    use utility_key_default
    use output
    use pack_mask_arrays
    use utility_crash
    use xsect_tables

    implicit none

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
            real(8)              :: arbitraryreal = 0.d0
            character(64)        :: subroutine_name = 'initialize_toplevel'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%initialization) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------  
        !% --- Set a small real based on machine precision
        !%     This produces a number that is significantly larger than machine
        !%     precision so that it can be a usable number, but small enough
        !%     to be irrelevant.
        setting%Eps%Machine = tenR**(-floor(sqrt(-log10(tiny(arbitraryreal)))))

        !% --- set the CPU and wall-clock timers
        call init_model_timer()

        !% --- define the reverse keys (used mainly for debugging)
        call define_keys_reverse()
        call define_apikeys_reverse()
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

        !% --- create duplicate input files
        !%     this is required because each image needs its own copy of the input files
        call util_file_duplicate_input ()
        call util_crashstop(2983)
        sync all

        !% --- initialize the time stamp used for output (must be after json is read)
        call init_timestamp ()
        sync all

        !% --- setup the output file directories. 
        !%     This will create a new directory with a timestamp for output
        call util_file_setup_output_folders()
        sync all

        !% --- print program header
        if ((setting%Output%Verbose) .and. (this_image() == 1)) &
             call util_print_programheader ()    

        !% --- set up the profiler
        if (setting%Profile%useYN) then
            call util_allocate_profiler ()
            call util_profiler_start (pfc_initialize_all)
        else 
            !% continue without profiler    
        end if

        !% --- set constraints on the setting% structure
        call util_setting_constraints ()

        !% --- initialize the coarrays that depend on number of images
        !%     and not on number of links/nodes, elements or faces.
        !if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, "begin initialize secondary coarrays"
        call util_allocate_secondary_coarrays ()

        !% --- initialize the API with the SWMM-C code
        !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin interface between SWMM-C and 5+"
        call interface_init ()
        call util_crashstop(43974)

        !% --- set up and store the SWMM-C link-node arrays in equivalent Fortran arrays
        !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin link-node processing"
        call init_linknode_arrays ()
        call util_crashstop(31973)

        !% --- initialize ForceMain settings (determines if FM is used)
        !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin Forcemain setting"
        call init_ForceMain_setting ()

        !% --- initialize Adjustments from EPA SWMM input file
       ! if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin get adjustments"
        call interface_get_adjustments ()

        !% --- setup the irregular transect arrays associated with SWMM-C input links
        !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin transect_arrays"
        call init_link_transect_array()
        call util_crashstop(42873)

        !% --- initialize globals that are run-time dependent
        !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin initialize globals"
        call init_globals()
        
        !% --- store the SWMM-C curves in equivalent Fortran arrays
        !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin SWMM5 curve processing"
        call init_curves()
        call util_crashstop(53454)

        !% --- read in profiles from .inp file and create 
        !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin SWMM5 profile processing"
        if (this_image() .eq. 1) then 
            call init_profiles()
        end if

        !% --- initialize culverts
        !if (setting%Output%Verbose) print *, "begin initializing culverts"
        call init_culvert()

        !% --- kinematic viscosity for water
        call init_viscosity()

        !%==========================================================================
        !%                      BEGIN PARTITIONING FOR PARALLEL                            
        !%      AFTER THIS POINT WE HAVE INSERTED NEW NODES AND SPLINT LINKS    
        !%==========================================================================
    
        !% --- break the link-node system into partitions for multi-processor operation
        !if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, "begin link-node partitioning"
        call init_partitioning()
        call util_crashstop(5297)

        !% --- HACK WORK NEEDED: we need to ensure that any phantom link defined in partitioning
        !% is NOT a culvert.i.e., the original portion of the link from SWMM must be defined as
        !% the culvert and the phantom must be upstream/downstream. 

        !% --- initialize types to the undefined key number
        elemI(:,ei_elementType)  = undefinedKey
        elemI(:,ei_geometryType) = undefinedKey
        elemI(:,ei_HeqType)      = undefinedKey
        elemI(:,ei_QeqType)      = undefinedKey
      
        !% --- error checking
        if (.not. setting%Simulation%useHydraulics) then 
            if (this_image() == 1) then
                write(*,*) 'USER CONFIGURATION ERROR setting.Simulation.useHydraulics == .false.'
                write(*,*) '...this presently is not supported in SWMM5+'
            end if
            call util_crashpoint(8815)
        end if  
        call util_crashstop(1973)

        !%==========================================================================
        !%                NETWORK DEFINITION ON EACH PROCESSOR IMAGE
        !%==========================================================================
        !% --- translate the link-node system into a finite-volume network
        !if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,"begin network define"
        call network_define_toplevel ()
        call util_crashstop(3293)

        !% --- LINK-ELEM DATA BROADCAST
        !%     ensures that all images have the unique data they need from other images after
        !%     partitioning and network definition
        !if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,"begin init linkarray broadcast"
        call init_linkarray_broadcast()
        call util_crashstop(550987)

        !% --- initialize boundary and ghost elem arrays for inter image data transfer
        !if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, "begin init boundary ghost"
        call init_boundary_ghost_elem_array ()
        call util_crashstop(2293)
        
        !% --- initialize the time variables
        !if (setting%Output%Verbose) print *, "begin initializing time"
        call init_time()

        !% --- initialize simple controls from json file
        !if (setting%Output%Verbose) print *, "begin initializing simulation controls"
        call init_simulation_controls() 

        !% --- HYDROLOGY
        if (setting%Simulation%useHydrology) then 
            if (setting%SWMMinput%N_subcatch > 0) then
                !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin subcatchment initialization"
                call init_subcatchment()
            else 
                if (this_image() == 1) then
                    write(*,'(A)') ' ... setting.Simulation.useHydrology requested, but no subcatchments found.'
                    write(*,'(A)') ' ... skipping hydrology in this simulation.'
                end if
                setting%Simulation%useHydrology = .false.
            end if
        else 
            !% continue without hydrology    
        end if
        call util_crashstop(320983)

        !%==========================================================================
        !%                                   OUTPUT SETUP
        !%==========================================================================
        !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin initializing output report"
        call init_report()

        !%==========================================================================
        !%                     SETUP INITIAL CONDITIONS ON ELEMENTS
        !%==========================================================================
        !% --- command line warning for long wait times
        if ((setting%Output%Verbose) .and. (this_image() == 1)) then 
            if ((N_link > 5000) .or. (N_node > 5000)) then
                    write(*,"(A)") " ... setting initial conditions --this may take several minutes for big systems ..."
                    write(*,"(A,i8,A,i8,A)") "      SWMM system has ", setting%SWMMinput%N_link, " links and ", setting%SWMMinput%N_node, " nodes"
                    write(*,"(A,i8,A)")      "      FV system has   ", sum(N_elem(:)), " elements"
            else 
                    !% --- no need to warn for small systems
            end if
        else 
            !% be silent    
        end if    

        !% --- initial conditions
        !if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, "begin init IC_toplevel"
        call init_IC_toplevel ()       
        call util_crashstop(4429873)

        !%-----------------------------------------------------------------------
        !%        NOTE: WE CAN CALL PACKED MAPS ep_... and fp_... AFTER THIS POINT
        !%-----------------------------------------------------------------------

        !% --- initialize blowup limits
        call util_crash_initialize

        !% --- allocate other temporary arrays (initialized to null)
        call util_allocate_temporary_arrays()

        !% --- initialize volume conservation storage for debugging
        elemR(:,er_VolumeConservation) = zeroR    

        !% --- setup the multi-level finite-volume output
        !%        HACK -- Ideally, this should be a procedure accessed in the output module, 
        !%        but that caused linking problems due to use of pack/mask calls
        if (setting%Output%Report%provideYN) then 
            if (setting%Simulation%useHydraulics) then !% 
                !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin setup of output files"
                !% --- Get the output element and face locations
                call outputML_selection ()                
                !% --- Create packed arrays of elem row numbers that are output
                call pack_element_outputML ()
                !% --- Create packed arrays of face row numbers that are output
                call pack_face_outputML ()
                !% --  Setup of the outputML data structures
                call outputML_setup ()
            else 
                if (this_image() == 1) then
                    write(*,*) 'USER CONFIGURATION ERROR setting.Simulation.useHydraulics == .false.'
                    write(*,*) '...this presently is not supported in SWMM5+'
                end if
                call util_crashpoint(487587)  
            end if  
        else 
            !% continue without any output files                                      
        end if
        call util_crashstop(103897)

        !% --- wait for all processors before exiting to the time loop
        sync all
 
        !%------------------------------------------------------------------- 
        !% Closing
            !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin init_check_setup_conditions"
            call init_check_setup_conditions()

            !if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin init_timer_stop"
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
                call util_crashpoint(333578)
            end if
            call util_crashstop(440987)

            if ((setting%Output%Verbose) .and. (this_image() == 1)) then 
                 print *, 'finished initialization'
                 print *, ' '
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
        !% --- store CPU clock start time
        call cpu_time(setting%Time%CPU%EpochStartSeconds)

        !% --- get/store the Wall clock time
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
        !%------------------------------------------------------------------
        !% Description:
        !% stops the timer at end of initialization
        !%------------------------------------------------------------------
        !% Declarations:
            integer(kind=8) :: crate, cmax, cval
        !%------------------------------------------------------------------
        
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
        !%------------------------------------------------------------------
        !% Description:
        !% initializes time stamp used for output files
        !% This is a character array in the form YYYYMMDD_hhmm
        !%-------------------------------------------------------------------
        !% Declarations:
            integer           :: thisTime(8)
            character(len=4)  :: cyear
            character(len=2)  :: cmonth, cday, chour, cmin
            character(len=13) :: datetimestamp
            character(64)     :: subroutine_name = 'init_timestamp'
        !%-------------------------------------------------------------------
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

        !% --- compute only on the first image
        if (this_image() == 1) then
            datetimestamp = cyear//cmonth//cday//'_'//chour//cmin
        endif

        !% --- ensure all processors use the same datetime stamp
        call co_broadcast (datetimestamp, source_image=1)

        setting%Time%DateTimeStamp = datetimestamp

        sync all

    end subroutine init_timestamp
!!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_linknode_arrays()
        !%------------------------------------------------------------------
        !% Description:
        !%   Retrieves data from EPA-SWMM interface and populates link and 
        !%   and node tables
        !% Note:
        !%   The order in which link and nodes are populated coincides with
        !%   the order in which links and nodes are allocated in EPA-SWMM 
        !%   data structures. Keeping the same order is important to be able 
        !%   to locate node/link data by label and not by index, reusing 
        !%   EPA-SWMM functionalities.
        !%------------------------------------------------------------------
        !% Declarations   
            integer          :: ii, jj, total_n_links, link_idx
            integer, pointer :: linkUp, linkDn
            logical          :: noerrorfound
            character(64)    :: subroutine_name = 'init_linknode_arrays'
        !%--------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%initialization) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

            if (.not. api_is_initialized) then
                print *, "CODE ERROR API is not initialized"
                call util_crashpoint(39873)
            end if
        !%-----------------------------------------------------------------------------
        !% --- Allocate storage for link & node tables
        call util_allocate_linknode()

        !% --- Allocate subcatchment storage
        if (setting%Simulation%useHydrology) then 
            call util_allocate_subcatch()
        else    
            !% --- continue without hydrology    
        end if

        !% --- Set default for all link and node keys
        call util_key_default_linknode()

        !% --- initialize number of links for each node to zero
        node%I(:,ni_N_link_u) = zeroI
        node%I(:,ni_N_link_d) = zeroI
        !% --- set default for extra surcharge depth
        !%     The InfiniteExtraDepthValue is used for SWMM5+ to identify junctions that
        !%     never overflow (i.e., without a manhole) and thus can be represented
        !%     as nJ2 faces between elements when there are only two connections.
        !%     We default ALL junction overflow height to the Infinite value. This
        !%     default is changed when link/node data is read in (further below)
        node%R(:,nr_OverflowHeightAboveCrown) = setting%Junction%InfiniteExtraDepthValue 

        !% -----------------------
        !% --- LINK DATA
        !% -----------------------
        !% --- cycle through links using do loop as the API between EPA SWMM and SWMM5+ 
        !%     is designed for individual links, not for array processing.
        do ii = 1, setting%SWMMinput%N_link

            !% --- store the basic link data
            link%I(ii,li_idx) = ii
            link%I(ii,li_link_direction) = interface_get_linkf_attribute(ii, api_linkf_direction,.true.)
            link%I(ii,li_link_type)      = interface_get_linkf_attribute(ii, api_linkf_type,     .true.)
            link%I(ii,li_link_sub_type)  = interface_get_linkf_attribute(ii, api_linkf_sub_type, .true.)
            link%I(ii,li_geometry)       = interface_get_linkf_attribute(ii, api_linkf_geometry, .true.)
            link%I(ii,li_barrels)        = interface_get_linkf_attribute(ii, api_linkf_conduit_barrels, .true.)
            link%I(ii,li_culvertCode)    = interface_get_linkf_attribute(ii, api_linkf_culvertCode, .true.)
                
            !% --- identify the upstream and downstream node indexes
            if (link%I(ii,li_link_direction) == oneI) then
                !% --- for standard channel/conduits where upstream is 
                !%     higher bottom elevation than downstream
                link%I(ii,li_Mnode_u) = interface_get_linkf_attribute(ii, api_linkf_node1,.true.) + oneI ! node1 in C starts from 0
                link%I(ii,li_Mnode_d) = interface_get_linkf_attribute(ii, api_linkf_node2,.true.) + oneI ! node2 in C starts from 0

                !% offset 1 is upstream for link with positive slope
                link%R(ii,lr_InletOffset)        = interface_get_linkf_attribute(ii, api_linkf_offset1,.false.)

                !% offset 2 is downstream for link with positive slope
                link%R(ii,lr_OutletOffset)       = interface_get_linkf_attribute(ii, api_linkf_offset2,.false.)

            else if (link%I(ii,li_link_direction) == -oneI) then
                !% --- when upstream has lower bottom elevation than downstream,
                !%     EPA-SWMM swaps the connection to prevent a negative slope. 
                !%     SWMM5+ can handle a negative slope in FV method, so switch the 
                !%     nodes back to their correct orientation
                link%I(ii,li_Mnode_u) = interface_get_linkf_attribute(ii, api_linkf_node2,.true.) + oneI ! node2 in C starts from 0
                link%I(ii,li_Mnode_d) = interface_get_linkf_attribute(ii, api_linkf_node1,.true.) + oneI ! node1 in C starts from 0

                !% offset 2 is upstream for link with negative/zero slope
                link%R(ii,lr_InletOffset)        = interface_get_linkf_attribute(ii, api_linkf_offset2, .false.)

                !% offset 1 is downstream for link with positive slope
                link%R(ii,lr_OutletOffset)       = interface_get_linkf_attribute(ii, api_linkf_offset1, .false.)

            else
                write(*,*) 'Fatal error: link direction should be 1 or -1'
                call util_crashpoint(7298734)
            end if 

            !% --- the "parent" link is the input EPA-SWMM link, 
            !%     which may be broken up (later) into 2 SWMM5+ links by partitioning
            link%I(ii,li_parent_link) = ii

            !% --- increment the connection counter for the node downstream
            node%I(link%I(ii,li_Mnode_d), ni_N_link_u) = node%I(link%I(ii,li_Mnode_d), ni_N_link_u) + oneI

            !% --- set the maps for the downstream node to upstream link (ni_Mlink_u#)
            !%     NOTE: this makes use of the ordering of the ni_Mlink_u# in define_indexes
            node%I(link%I(ii,li_Mnode_d), ni_idx_base1 + node%I(link%I(ii,li_Mnode_d), ni_N_link_u)) = ii

            !% --- increment the connection counter for the node upstream
            node%I(link%I(ii,li_Mnode_u), ni_N_link_d) = node%I(link%I(ii,li_Mnode_u), ni_N_link_d) + oneI

            !% --- set the maps for the upstream node to the downstream link (ni_Mlink_d#)
            !%     NOTE: this makes use of the ordering of the ni_Mlink_d# in define_indexes
            node%I(link%I(ii,li_Mnode_u), ni_idx_base2 + node%I(link%I(ii,li_Mnode_u), ni_N_link_d)) = ii

            !% --- set the approach used for computing the initial depth
            !%     HACK -- At this time, all links have the same initial depth type which is 
            !%     the default stored in setting. A better approach would be to allow different links
            !%     to have different depth types. However, this requires changes to the *.inp file or
            !%     development of a new auxiliary input file.
            link%I(ii,li_InitialDepthType)   = setting%Link%DefaultInitDepthType
            link%R(ii,lr_Length)             = interface_get_linkf_attribute(ii, api_linkf_conduit_length,   .false.)
            link%R(ii,lr_BreadthScale)       = interface_get_linkf_attribute(ii, api_linkf_xsect_wMax,       .false.)
            link%R(ii,lr_LeftSlope)          = interface_get_linkf_attribute(ii, api_linkf_left_slope,       .false.)
            link%R(ii,lr_RightSlope)         = interface_get_linkf_attribute(ii, api_linkf_right_slope,      .false.)
            link%R(ii,lr_Roughness)          = interface_get_linkf_attribute(ii, api_linkf_conduit_roughness,.false.)
            link%R(ii,lr_FullDepth)          = interface_get_linkf_attribute(ii, api_linkf_xsect_yFull,      .false.)
            link%R(ii,lr_FullArea)           = interface_get_linkf_attribute(ii, api_linkf_xsect_aFull,      .false.)
            link%R(ii,lr_FullHydRadius)      = interface_get_linkf_attribute(ii, api_linkf_xsect_rFull,      .false.)
            link%R(ii,lr_BottomDepth)        = interface_get_linkf_attribute(ii, api_linkf_xsect_yBot,       .false.)
            link%R(ii,lr_BottomRadius)       = interface_get_linkf_attribute(ii, api_linkf_xsect_rBot,      .false.)
            link%R(ii,lr_FlowrateInitial)    = interface_get_linkf_attribute(ii, api_linkf_q0,               .false.)
            link%R(ii,lr_FlowrateLimit)      = interface_get_linkf_attribute(ii, api_linkf_qlimit,           .false.)
            link%R(ii,lr_Kconduit_MinorLoss) = interface_get_linkf_attribute(ii, api_linkf_cLossAvg,         .false.)
            link%R(ii,lr_Kentry_MinorLoss)   = interface_get_linkf_attribute(ii, api_linkf_cLossInlet,       .false.)
            link%R(ii,lr_Kexit_MinorLoss)    = interface_get_linkf_attribute(ii, api_linkf_cLossOutlet,      .false.)
            link%R(ii,lr_SeepRate)           = interface_get_linkf_attribute(ii, api_linkf_seepRate,         .false.)
            link%R(ii,lr_ForceMain_Coef)     = interface_get_linkf_attribute(ii, api_linkf_forcemain_coef,   .false.)
            link%R(ii,lr_Setting)            = interface_get_linkf_attribute(ii, api_linkf_setting,          .false.)
            link%R(ii,lr_TimeLastSet)        = interface_get_linkf_attribute(ii, api_linkf_timelastset,     .false.)

            !% --- Note that link%R(ii,lr_Slope) and link%R(ii,lr_TopWidth) are defined in network_define.f08 
            !%     because SWMM5 reverses negative slope

            !% --- Note, link depths are NOT initialized here because these are determined by node attributes
            
    
            !% --- special element attributes
            link%I(ii,li_weir_EndContractions) = interface_get_linkf_attribute(ii, api_linkf_weir_end_contractions,.true.)
            link%I(ii,li_RoadSurface)          = interface_get_linkf_attribute(ii, api_linkf_weir_road_surface,    .true.)
            link%I(ii,li_curve_id)             = interface_get_linkf_attribute(ii, api_linkf_curveid,              .true.)
            link%R(ii,lr_DischargeCoeff1)      = interface_get_linkf_attribute(ii, api_linkf_discharge_coeff1,     .false.)
            link%R(ii,lr_DischargeCoeff2)      = interface_get_linkf_attribute(ii, api_linkf_discharge_coeff2,     .false.)
            link%R(ii,lr_initSetting)          = interface_get_linkf_attribute(ii, api_linkf_initSetting,          .false.)
            link%R(ii,lr_yOn)                  = interface_get_linkf_attribute(ii, api_linkf_yOn,                  .false.)
            link%R(ii,lr_yOff)                 = interface_get_linkf_attribute(ii, api_linkf_yOff,                 .false.)
            link%R(ii,lr_SideSlope)            = interface_get_linkf_attribute(ii, api_linkf_weir_side_slope,      .false.)
            link%R(ii,lr_RoadWidth)            = interface_get_linkf_attribute(ii, api_linkf_weir_road_width,      .false.)

            if (interface_get_linkf_attribute(ii, api_linkf_hasFlapGate,.true.) == oneI) then
                link%YN(ii,lYN_hasFlapGate)   = .true.
            else
                link%YN(ii,lYN_hasFlapGate)   = .false.
            end if

            !% --- EPA SWMM5 does not distinguish between open channels and closed conduits
            !%     however SWMM5+ needs that distinction to set up the initial conditions
            !%     So we look for EPA SWMMM conduits (given lPipe, above) and convert to lChannel
            if ( (link%I(ii,li_link_type) == lPipe)          .and. &
                 ( &
                    (link%I(ii,li_geometry) == lRectangular)    .or. &
                    (link%I(ii,li_geometry) == lTrapezoidal)    .or. &
                    (link%I(ii,li_geometry) == lTriangular)     .or. &
                    (link%I(ii,li_geometry) == lParabolic)      .or. &
                    (link%I(ii,li_geometry) == lPower_function) .or. &
                    (link%I(ii,li_geometry) == lIrregular)           &
                  )  &
                ) then

                link%I(ii,li_link_type) = lChannel
            end if

            !% --- check if a weir can surcharge
            !%     returns false if not a weir
            link%YN(ii,lYN_weir_CanSurcharge) = (interface_get_linkf_attribute(ii,api_linkf_weir_can_surcharge,.true.) == oneI)

            if (link%I(ii,li_link_type) == lWeir) then
                !% --- set road surface types
                if (link%I(ii,li_RoadSurface) == API_NOSURFACE) then
                    link%I(ii,li_RoadSurface) = NoRoadSurface
                else if (link%I(ii,li_RoadSurface) == API_PAVED) then
                    link%I(ii,li_RoadSurface) = Paved
                else if (link%I(ii,li_RoadSurface) == API_GRAVEL) then
                    link%I(ii,li_RoadSurface) = Gravel
                else
                    if (this_image() == 1) then
                        write(*,*) 'USER CONFIGURATION ERROR for roadway weir'
                        write(*,"(A,i4,A)") 'A roadway weir does not have an allowed road surface type'
                        write(*,"(A)")      'Allowed types are NOSURFACE, PAVED, GRAVEL'
                        write(*,"(A,i4)")   'Failure at link ',link%I(ii,li_idx)
                        write(*,"(A)")      'link name '//trim(link%Names(ii)%str)
                    end if
                    call util_crashpoint(548976)
                end if
            else
                link%I(ii,li_RoadSurface) = nullValueI
            end if
            
            ! !% HACK CODE FOR TESTING:
            ! !% for filled circular cross-sections, swmm always sets inlet and outlet offsets
            ! !% for the bottom filled elevation. For now, I am removing those for testing
            ! if (link%I(ii,li_geometry) ==  lFilled_circular) then
            !     link%R(ii,lr_InletOffset)  = link%R(ii,lr_InletOffset)  - link%R(ii,lr_BottomDepth)
            !     link%R(ii,lr_OutletOffset) = link%R(ii,lr_OutletOffset) - link%R(ii,lr_BottomDepth) 
            ! end if

            !% --- Irregular cross-sections (TRANSECTS in SWMM input file)
            if (link%I(ii,li_geometry) == lIrregular) then
               link%I(ii,li_transect_idx) = interface_get_linkf_attribute(ii, api_linkf_transectidx,.true.)
            end if

            !% --- set output links
            link%YN(ii,lYN_isOutput) = (interface_get_linkf_attribute(ii,api_linkf_rptFlag,.true.) == 1)

            !% HACK: experimental code
            !% --- replace smaller links using an equivalent orifice
            if (setting%Discretization%UseEquivalentOrifice) then

                if ( (link%R(ii,lr_Length) < setting%Discretization%MinLinkLength)                   .and. &
                     ((link%I(ii,li_link_type) == lpipe) .or. (link%I(ii,li_link_type) == lchannel)) .and. &
                     (.not. link%YN(ii,lYN_isPhantomLink))                                            ) then
                    
                    if (this_image() == 1) then
                        print*, 'WARNING: Converting link to equivalent orifice'
                        print*, 'Link index is ',ii,' link name is ',  trim(link%Names(ii)%str)
                        print*, 'has lenght of ', link%R(ii,lr_Length) 
                        print*, 'Which is smaller than user defined minimum link lenght of ', setting%Discretization%MinLinkLength
                        print*
                    end if
                    
                    !% set the link type type as Orifice
                    link%I(ii,li_link_type) = lOrifice
                    !% set the sub orifice type as side orifice
                    link%I(ii,li_link_sub_type) = lSideOrifice 
                    !% set a default discharge coefficient from settings
                    link%R(ii,lr_DischargeCoeff1) = setting%Discretization%EquivalentOrificeDischargeCoeff
                    !% set orifice time to operate to zero
                    link%R(ii,lr_DischargeCoeff2) = zeroR
                    !% set the default geometry of the equivalent orifice as circular
                    link%I(ii,li_geometry) = lCircular
                    !% reset the orifice opening from the original link full area
                    link%R(ii,lr_FullDepth) = sqrt(fourR * link%R(ii,lr_FullArea) / setting%Constant%pi)
                    !% reset the length of the element as minimum link lenght (will be reset later)
                    link%R(ii,lr_Length) = setting%Discretization%MinLinkLength
                end if

            end if

        end do

        !% --- ERROR CHECK for number of connections
        do ii = 1,N_node
            if (node%I(ii, ni_N_link_u) > max_up_branch_per_node) then
                if (this_image() == 1) then
                    write(*,*) 'USER CONFIGURATION ERROR for node connections'
                    write(*,"(A,i4,A)") 'One or more nodes have more than ',max_up_branch_per_node,' upstream connections'
                    write(*,*) 'Unfortunately, this connection limit is a hard-coded limit of SWMM5+ an cannot be exceeded.'
                    write(*,*) 'First error found at node ',ii
                    write(*,*) 'Node name ',trim(node%Names(ii)%str)
                end if
                call util_crashpoint(387666)
            end if

            if (node%I(ii, ni_N_link_d) > max_dn_branch_per_node) then
                if (this_image() == 1) then
                    write(*,*) 'USER CONFIGURATION ERROR for node connections'
                    write(*,"(A,i4,A)") 'One or more nodes have more than ',max_dn_branch_per_node,' downstream connections'
                    write(*,*) 'Unfortunately, this connection limit is a hard-coded limit of SWMM5+ an cannot be exceeded.'
                    write(*,*) 'First error found at at node ',ii
                    write(*,*) 'Node name ',trim(node%Names(ii)%str)
                end if
                call util_crashpoint(86752)
            end if
        end do

        !% --- Store the Link/Node names (moved here 20221216)
        call interface_update_linknode_names()

        !% -----------------------
        !% --- SUBCATCHMENT -- see also init_subcatchment
        !% -----------------------
        node%I(:,ni_routeFrom) = nullvalueI !% initialization 
        if (setting%Simulation%useHydrology) then
            do ii=1,setting%SWMMinput%N_subcatch
                !% --- Subtract 1 from the SWMM5+ subcatchment index to get the 
                !%     EPA SWMM subcatchment index; Add 1 to the EPA SWMM node index
                !%     for the SWMM5+ node index
                subcatchI(ii,si_runoff_nodeIdx) = interface_get_subcatch_runoff_nodeIdx(ii-1)+oneI
                if (subcatchI(ii,si_runoff_nodeIdx) < oneI) then !% not a runoff node (EPA SWMM flag)
                    subcatchYN(ii,sYN_hasRunoff) = .false.
                    subcatchI(ii,si_runoff_nodeIdx) = nullvalueI
                else
                    subcatchYN(ii,sYN_hasRunoff) = .true.
                    node%I(subcatchI(ii,si_runoff_nodeIdx),ni_routeFrom) = ii
                end if  
            end do
        end if
    
        !% -----------------------
        !% --- NODE DATA
        !% -----------------------
        do ii = 1, N_node

            !% --- get the total number of links connected to this node by
            !%     adding upstream and downstream
            total_n_links = node%I(ii,ni_N_link_u) + node%I(ii,ni_N_link_d)
            node%I(ii, ni_idx) = ii

            !% --- check for node inflows
            node%YN(ii, nYN_has_extInflow) = (interface_get_nodef_attribute(ii, api_nodef_has_extInflow) == 1)
            node%YN(ii, nYN_has_dwfInflow) = (interface_get_nodef_attribute(ii, api_nodef_has_dwfInflow) == 1)

            !% --- get initial depths
            node%R(ii,nr_InitialDepth)      = interface_get_nodef_attribute(ii, api_nodef_initDepth)

            !% error check
            if (node%R(ii,nr_InitialDepth) < zeroR) then
                print *, 'USER CONFIGURATION ERROR OR CODE ERROR for node depth and invert'
                print *, 'The initial depth at a node is less than the invert elevation'
                print *, 'This may occur if the user provided an input of depth at an '
                print *, 'outfall rather than stage elevation.  Otherwise, this is'
                print *, 'possibly a bug in the code'
                print *, 'Node # ',ii
                print *, 'Node Name ', trim(node%Names(ii)%str)
                print *, 'initial depth ',node%R(ii,nr_InitialDepth) 
                call util_crashpoint(77987232)
            end if

            !% --- node geometry
            node%R(ii,nr_Zbottom)           = interface_get_nodef_attribute(ii, api_nodef_invertElev)
            node%R(ii,nr_FullDepth)         = interface_get_nodef_attribute(ii, api_nodef_fullDepth)

            !% --- Total pressure head above max depth allowed for surcharge
            !%     If 0 then node cannot surcharge, so exceeding depth means water either ponds or is lost
            node%R(ii,nr_OverflowHeightAboveCrown) = interface_get_nodef_attribute(ii, api_nodef_surDepth)

            !% --- storage equations
            !%     Note that Storage is converted to SI units in api.c/api_get_nodef_attribute
            node%R(ii,nr_StorageConstant)   = interface_get_nodef_attribute(ii, api_nodef_StorageConstant)
            node%R(ii,nr_StorageExponent)   = interface_get_nodef_attribute(ii, api_nodef_StorageExponent)
            node%R(ii,nr_StorageCoeff)      = interface_get_nodef_attribute(ii, api_nodef_StorageCoeff)
            node%I(ii,ni_curve_ID)          = interface_get_nodef_attribute(ii, api_nodef_StorageCurveID)
            node%R(ii,nr_StorageFevap)      = interface_get_nodef_attribute(ii, api_nodef_StorageFevap)

            !% --- ponded area
            if (setting%SWMMinput%AllowPonding) then
                node%R(ii,nr_PondedArea) = interface_get_nodef_attribute(ii, api_nodef_PondedArea)
            else
                node%R(ii,nr_PondedArea) = zeroR
            end if

            !% --- set if node is designiated for output
            node%YN(ii,nYN_isOutput)  = (interface_get_nodef_attribute(ii, api_nodef_rptFlag) == 1)

            !%
            !% --- Assign required node types nJm, nJ1, nJ2, nBCdn,
            !%     Note that defined storage is ALWAYS nJM
            !%     The goal is to identify nodes that have only two connections and could
            !%     be represented by a face (nJ2) between two elements in the SWMM5+ FV system rather 
            !%     than requiring the more complicated junction (nJm) solution.
            if (interface_get_nodef_attribute(ii, api_nodef_type) == API_OUTFALL) then
                !% --- OUTFALL NODES
                node%I(ii, ni_node_type) = nBCdn

                !% --- get the SWMMoutfall index (i.e. the 'k' in Outfall[k].vRouted in EPASWMM)
                node%I(ii,ni_SWMMoutfallIdx) = interface_get_nodef_attribute(ii,api_nodef_outfall_idx)

                !% --- check for a flap gate on an outfall
                if (interface_get_nodef_attribute(ii, api_nodef_hasFlapGate) == oneI) then
                    node%YN(ii, nYN_hasFlapGate) = .true.
                else
                    node%YN(ii, nYN_hasFlapGate) = .false.
                end if

                !% --- check for routeTo subcatchment
                node%I(ii,ni_routeTo) = interface_get_nodef_attribute(ii,api_nodef_RouteTo)
                if (node%I(ii,ni_routeTo) == -oneI) then
                    node%I(ii,ni_routeTo) = nullvalueI
                elseif ( (node%I(ii,ni_routeTo) > zeroI)                           &
                     .and.                                                         &
                         (node%I(ii,ni_routeTo) .le. setting%SWMMinput%N_subcatch) &
                     ) then
                    !% correct value found
                else
                    print *, 'CODE ERROR unexpected value for ni_routeTo'
                    print *, 'value obtained is ',node%I(ii,ni_routeTo)
                    print *, 'allowable values are -1 or > 0 but less than number of subcatchments'
                    print *, 'number of subcatchments is ',setting%SWMMinput%N_subcatch
                    print *, 'Node index ',ii
                    print *, 'Node name  ',trim(node%Names(ii)%str)
                    call util_crashpoint(69873)
                end if

            else if (interface_get_nodef_attribute(ii, api_nodef_type) == API_STORAGE) then
                !% --- STORAGE NODES (always nJm)
                node%I(ii, ni_node_type)     = nJm
                node%YN(ii, nYN_has_storage) = .true.

            else 
                !% --- OTHER NODES
                !%     classify by number of links connected
                select case (total_n_links)
                    case (oneI)
                        node%I(ii, ni_node_type) = nJ1
                    case (twoI)      
                        node%I(ii, ni_node_type) = nJ2
                    case default 
                        node%I(ii, ni_node_type) = nJm
                end select
            end if 
            
            !% --- nJ2 strictly has one upstream and one downstream link and
            !%     cannot be a subcatchment outlet  
            !%     Other cases where an nJ2 has only (i) two upstream and no downstream links,
            !%     or (ii) two downstream and no upstream links, will be set to nJm
            if ( (node%I(ii,ni_node_type)  ==  nJ2)             &
                .and.                                           &
                    ((node%I(ii,ni_N_link_u)   >   oneI)  .or.  &
                     (node%I(ii,ni_N_link_d)   >   oneI)  .or.  &
                     (node%I(ii,ni_routeFrom) .ne. nullvalueI)  &
                    )                                           &
                ) then
                    !% ... switching to a 2 link nJm junction type'
                    node%I(ii, ni_node_type) = nJm
            end if

            !% ==========================================================================
            !% --- Further discrimination between 2-element junctions that are nJ2
            !%     and those that are classed nJm. Note that all defined STORAGE 
            !%     junctions are already set to nJm, so this only applies to junctions 
            !%     defined in SWMM input file without explicit storage
            !%  
            !%     The following "or" conditions must be met for an nJ2:
            !%     1. at least one connected element is open-channel AND ponding_Area = 0 AND
            !%        the OverflowDepth = 0
            !%     2. both elements are NOT open channel AND the junction extra
            !%        surcharge depth == Junction.InfiniteExtraDepthValue 
            !%        (i.e., no possible overflow or ponding)
            !%     3. Downstream link may NOT be a Type1 Pump
            !%     In addition, the offsets of connected links must be zero, unless
            !%     the connected link is a weir or orifice (their offset has a different
            !%     meaning.)
            !%     Key point is that nJ2 cannot have overflow or ponding, so if
            !%     at least one side is open channel and ponding area = 0 and the
            !%     surcharge extra depth = 0 it is treated as open channel 
            !%     (i.e., overflow occurs in the adjacent channel element) so
            !%     it can be nJ2.  If both sides are closed types (conduit, weir,
            !%     i.e., not open channel) then the junction must also be
            !%     be closed; thus if a value (other than InfiniteExtraDepthValue)
            !%     is provided for the extra surcharge, then the junction must be treated
            !%     as an nJm rather than nJ2. That is, setting the Surcharge Extra Depth
            !%     to the InfiniteExtraDepthValue implies a non-vented connection that
            !%     can be treated as a face.
            !%     In general, existence of non-zero offsets require an nJm unless 
            !%     the offset is associated with a weir or orifice

            if (node%I(ii, ni_node_type) == nJ2) then
                !% --- local aliases for the upstream and downstream links. These should
                !%     be guaranteed to be in the u1 and d1 positions
                linkUp => node%I(ii,ni_Mlink_u1)
                linkDn => node%I(ii,ni_Mlink_d1)
                    
                !% --- phantom nodes will always be a nJ2
                if  (node%YN(ii,nYN_is_phantom_node)) then
                    !% --- no action: retain nJ2

                !% --- special channels and conduits that allow nJ2
                !%     note the "meters_per_ft" conditional is to allow an input file
                !%     in CFS to use the infinite depth value in the settings to be
                !%     interpreted as feet.
                elseif  ( ( (link%I(linkUp,li_link_type) .eq. lChannel)                     &
                            .or.                                                            &
                            (link%I(linkDn,li_link_type) .eq. lChannel)                     &
                          )                                                                 &
                          .and.                                                             &
                          (node%R(ii,nr_PondedArea) == zeroR)                               &
                          .and.                                                             &
                          (  (node%R(ii,nr_OverflowHeightAboveCrown) == zeroR)              & 
                              .or.                                                          &
                             (node%R(ii,nr_OverflowHeightAboveCrown)                        &
                               == setting%Junction%InfiniteExtraDepthValue)                 &
                             .or.                                                           &
                             ( (node%R(ii,nr_OverflowHeightAboveCrown)                              &
                                  < setting%Junction%InfiniteExtraDepthValue*meters_per_ft + 0.01d0) &
                              .and.                                                                 &
                              (node%R(ii,nr_OverflowHeightAboveCrown)                               &
                                  > setting%Junction%InfiniteExtraDepthValue*meters_per_ft - 0.01d0) &
                             )                                                                      & 
                          )                                                                         &
                        ) then
                        !% retain nJ2 OPEN CHANNEL FACE
                        !% --- if either link is an open channel AND the ponded area
                        !%     is zero then the junction is an nJ2 face where any
                        !%     overflow is handled by adjacent channel (i.e. lost). Otherwise 
                        !%     reverts to nJm element with its own overflow/ponding. 
                        !%     Note that if ponding is OFF but the ponded area
                        !%     is defined, then the element is treated as nJm with
                        !%     overflow above the Surcharge Extra Depth
                        !% --- no action: retain nJ2

                elseif ( (link%I(linkUp,li_link_type) .ne. lChannel)         &
                         .and.                                               &
                         (link%I(linkDn,li_link_type) .ne. lChannel)         &
                         .and.                                               &
                         (( (node%R(ii,nr_OverflowHeightAboveCrown)                              &
                                < setting%Junction%InfiniteExtraDepthValue + 0.01d0)             &
                            .and.                                                                &
                            (node%R(ii,nr_OverflowHeightAboveCrown)                              &
                                > setting%Junction%InfiniteExtraDepthValue - 0.01d0)             &
                            )                                                                    &   
                           .or.                                              &
                           ( (node%R(ii,nr_OverflowHeightAboveCrown)                               &
                                < setting%Junction%InfiniteExtraDepthValue*meters_per_ft + 0.01d0) &
                            .and.                                                                  &
                            (node%R(ii,nr_OverflowHeightAboveCrown)                                &
                                > setting%Junction%InfiniteExtraDepthValue*meters_per_ft - 0.01d0) &
                            )                                                                      & 
                          )                                                                        &
                          ) then
                        !% nJ2 CLOSED CONDUIT FACE
                        !% --- if both links are NOT open channel AND the OverflowDepth
                        !%     is equal to the InfiniteExtraDepthValue, then this is retained 
                        !%     as an nJ2 (unvented)  face. Otherwise switched to a vented nJM element.
                        !%     HACK -- if 1000 ft is put in a CFS input file or 1000 m in an SI
                        !%     input file this is treated as infinite depth.
                        !%  
                        
                        !% --- no action: retain nJ2

                elseif  (( (link%I(linkUp,li_link_type) .eq. lWeir)         &
                          .or.                                              &
                           (link%I(linkDn,li_link_type) .eq. lWeir)         &
                         )                                                  &
                        .and. (.not. setting%Weir%ForceWeirNodesToJM)       &
                        ) then          
                        !% nJ2 Weir face 
                        !% --- an nJ2 weir face is retained without as long as
                        !%     the force is not in place. Note that nodes with
                        !%     storage already have nJM, so this does not affect
                        !%     them 
                        
                        !% no action: retain nJ2

                elseif  (( (link%I(linkUp,li_link_type) .eq. lOrifice)         &
                            .or.                                               &
                           (link%I(linkDn,li_link_type) .eq. lOrifice)         &
                           )                                                   &
                          .and. (.not. setting%Orifice%ForceOrificeNodesToJM)  &
                         ) then    
                        !% nJ2 Orifice face 
                        !% --- an nJ2 orifice face is retained without as long as
                        !%     the force is not in place. Note that nodes with
                        !%     storage already have nJM, so this does not affect
                        !%     them 
                        
                        !% no action: retain nJ2
 
                else
                    !% --- switch to nJm

                    node%I(ii, ni_node_type) = nJm

                    ! write(*,*) '...NOTE: ',trim(node%Names(ii)%str),' is held as nJM node rather than nJ2 (faces).'
                    ! write(*,*) '   This occurs because the input SurchargeDepth is less than the InfiniteDepthValue'
                    ! write(*,*) '   Input SurcharegeDepth is ',node%R(ii,nr_OverflowHeightAboveCrown)
                    ! write(*,*) '   InfiniteDepthValue is    ',setting%Junction%InfiniteExtraDepthValue
                end if

                !% --- regardless of the above, if either link up or down is a
                !%     multi-barrel link, then the node must be an nJm node
                if  (  (link%I(linkUp,li_barrels) > oneI)     &
                        .or.                                  &
                        (link%I(linkDn,li_barrels) > oneI)    &
                    ) then       
                    node%I(ii, ni_node_type) = nJm   
                    
                end if

                !% --- regardless of the above, if the downstream link is
                !%     at type 1 pump, then the node must be nJm
                if (  (link%I(linkDn,li_link_type) .eq. lPump)      &
                       .and.                                        &
                      (link%I(linkDn,li_link_sub_type)  .eq. lType1Pump) &
                    ) then
                    node%I(ii, ni_node_type) = nJm 
                end if

                !% --- regardless of the above, if either link up or down
                !%     is a culvert then the node must be nJm
                if (    (link%I(linkUp,li_culvertCode) > 0)      &
                        .or.                                     &
                        (link%I(linkDn,li_culvertCode) > 0)      &
                    ) then
                    node%I(ii, ni_node_type) = nJm             
                end if   

            
                !% --- further check on offsets for any nJ2 that passed the prior
                !%     restrictions. In general, we have an nJm if there are 
                !%     any offsets, except if the offset is a weir or orifice.
                if (link%R(linkUp,lr_OutletOffset) .ne. zeroR) then
                    !% --- offsets are OK for upstream weir or orifice links  
                    if (  (link%I(linkUp,li_link_type) .eq. lWeir)       &
                        .or.                                          &
                        (link%I(linkUp,li_link_type) .eq. lOrifice)    &
                        ) then    
                        !% --- retain nJ2
                    else
                        !% --- switch to nJm
                        node%I(ii, ni_node_type) = nJm
                    end if
                end if


                if (link%R(linkDn,lr_InletOffset) .ne. zeroR) then
                    !% --- offsets are OK for downstream weir or orifice links  
                    if (  (link%I(linkDn,li_link_type) .eq. lWeir)       &
                        .or.                                          &
                        (link%I(linkDn,li_link_type) .eq. lOrifice)    &
                        ) then    
                        !% --- retain nJ2
                    else
                        !% --- switch to nJm
                        node%I(ii, ni_node_type) = nJm
                    end if
                end if 

                !% --- force some or all of the nJ2 to nJm (used for debugging)
                !% --- global forcing of all nodes
                if (setting%Junction%ForceNodesJM ) then
                    !% --- switch to nJm 
                    node%I(ii, ni_node_type) = nJm
                else
                    !% -- forcing of weir-adjacent nodes, only
                    if ( ((link%I(linkDn,li_link_type) .eq. lWeir) &
                          .or.                                     &
                          (link%I(linkUp,li_link_type) .eq. lWeir) &
                         ) .and.                                   &
                         (setting%Weir%ForceWeirNodesToJM)         &
                        ) then 
                        node%I(ii,ni_node_type) = nJm
                    end if
                    !% -- forcing of orifice-adjacent nodes, only
                    if ( ((link%I(linkDn,li_link_type) .eq. lOrifice)  &
                          .or.                                         &
                          (link%I(linkUp,li_link_type) .eq. lOrifice)  &
                         ) .and.                                       &
                         (setting%Orifice%ForceOrificeNodesToJM)      &
                        ) then 
                        node%I(ii,ni_node_type) = nJm
                    end if
                end if
            end if

            !% --- end nJ2, nJm processing
            !% ==========================================================================

            !% --- inflows (must be initialized after nJ1, nJ2 are set)
            if (node%YN(ii, nYN_has_extInflow) .or. node%YN(ii, nYN_has_dwfInflow)) then
                !% set inflow to true for any node type
                node%YN(ii, nYN_has_inflow) = .true.
                !% change the node type of an nJ1 with inflow to nBCup
                if (node%I(ii, ni_node_type) == nJ1) node%I(ii, ni_node_type) = nBCup
            end if


            !% --- calculate the full volumes of the links upstream of a node
            !%     this volume will be used later the distribute lateral inflows across links
            node%R(ii,nr_UpLinksFullVolume) = zeroR

            do jj = 1,node%I(ii,ni_N_link_u)
                !% pointer to the upstream links of that node
                link_idx = node%I(ii,ni_idx_base1 + jj)

                !% only calculate the volumes of open channels 
                if (link%I(link_idx,li_link_type) == lchannel) then
                    !% add the the volume of the upstream links
                    node%R(ii,nr_UpLinksFullVolume) = node%R(ii,nr_UpLinksFullVolume) &
                                    + link%R(link_idx,lr_FullArea) * link%R(link_idx,lr_Length)    
                end if

            end do

            !% --- note pattern initialization MUST be called after inflows are set
            node%I(ii,ni_pattern_resolution) = interface_get_BC_resolution(ii)
        end do

        !% --- find the volume fraction metric for channel links to distribute lateral inflows
        do ii = 1, setting%SWMMinput%N_link
            if (link%I(ii,li_link_type) == lchannel) then
                link%R(ii,lr_VolumeFractionMetric) = link%R(ii,lr_FullArea) * link%R(ii,lr_Length) &
                                                / node%R(link%I(ii,li_Mnode_d),nr_UpLinksFullVolume)
            else
                !% --- lateral inflows will not be distributed for non channel links.
                link%R(ii,lr_VolumeFractionMetric) = zeroR
            end if

            !% save the node from which lateral inflows will be distributed
            link%I(ii,li_lateralInflowNode)    = link%I(ii,li_Mnode_d)

        end do

        !% --- error checking for disconnected nodes
        noerrorfound = .true.
        do ii = 1,N_node
            if (node%I(ii,ni_N_link_u) + node%I(ii,ni_N_link_d) == zeroI) then
                noerrorfound = .false.
                print *, 'USER CONFIGURATION ERROR disconnected node.'
                print *, 'A node has been found that does not connect to any links.'
                print *, 'It must be commented out in the input file.'
                print *, 'Node name is ',trim(node%Names(ii)%str)
            end if
        end do
        if (.not. noerrorfound) then
            call util_crashpoint(72813)
        endif
     
        !% --- ERROR CHECK connections for outfalls
        do ii = 1,N_node
            if (interface_get_nodef_attribute(ii, api_nodef_type) == API_OUTFALL) then
                !% --- check if outfall has more than one upstream connection
                if( node%I(ii, ni_N_link_u) > 1) then
                    write(*,*) 'USER CONFIGURATION ERROR for outfall'
                    write(*,*) 'Outfall has more than 1 upstream link, which is not supported by SWMM5+'
                    write(*,"(A,A)") 'Node name in file is ',trim(node%Names(ii)%str)
                    call util_crashpoint(99374)
                endif
                !% --- check if outfall has a downstream connection
                if ( node%I(ii, ni_N_link_d) > 0) then
                    write(*,*) 'USER CONFIGURATION ERROR for outfall'
                    write(*,*) 'Outfall has a downstream connection, which is not supported by SWMM5+'
                    write(*,"(A,A)") 'Node name iin file is ',trim(node%Names(ii)%str)
                    call util_crashpoint(847822)
                end if
            end if
        end do

        !% Check for small links while using nominal element length discretization 
        if (setting%Discretization%UseNominalElemLength) then
            do ii = 1, N_link
                if ( (link%I(ii,li_link_type) == lChannel) .or. (link%I(ii,li_link_type) == lPipe) ) then
                    if (link%R(ii,lr_Length) < (real(setting%Discretization%MinElementPerLink,8) * setting%Discretization%NominalElemLength)) then
                        print *, 'SWMM input file links too small for selected NominalElemLength and MinElementPerLink'
                        print *, 'Found link length of ',link%R(ii,lr_Length)
                        print *, 'Link index is ',ii,' link name is ',  trim(link%Names(ii)%str)
                        print *, 'setting.Discretization.NominalElemLength is ',setting%Discretization%NominalElemLength
                        print *, 'setting.Discretization.MinElementPerLink is ',setting%Discretization%MinElementPerLink
                        print *, 'Decrease nominal element length to less than', link%R(ii,lr_Length)/ real(setting%Discretization%MinElementPerLink,8)
                        print *, 'or modify link length in SWMM input file'
                        print *, 'NOTE: SWMM5+ requires MinElementPerLink of 3 or greater'
                        call util_crashpoint(447298)
                    end if
                end if
            end do
        end if

        if (setting%Debug%File%initialization) then
            !%-----------------------------------------------------------------------------
            print *, 'idx,    nodeType,    linkU,    linkD,   curveID, patternRes'
            do ii=1,N_node
                write(*,"(10i8)") node%I(ii,ni_idx), node%I(ii,ni_node_type), node%I(ii,ni_N_link_u), node%I(ii,ni_N_link_d) &
                , node%I(ii,ni_curve_ID), node%I(ii,ni_pattern_resolution)
            end do

            print *, 'idx,    LinkType,  nodeU,   nodeD, curveId'
            do ii=1,N_Link
                write(*,"(10i8)") link%I(ii,li_idx), link%I(ii,li_link_type), link%I(ii,li_Mnode_u), link%I(ii,li_Mnode_d) &
                , link%I(ii,li_curve_id) 
            end do
        end if 

        if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"


    end subroutine init_linknode_arrays
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine init_link_transect_array ()
        !%------------------------------------------------------------------
        !% Description:
        !% Initializes the transect storage for EPA-SWMM irregular cross-sections
        !% Note that the transect index number DOES NOT correspond to the link
        !% index number. There are generally more links than transects, and a
        !% transect may be assigned to more than one link.
        !% Note that setting%SWMMinput%N_transect is set from EPA-SWMM in interface_init()
        !%------------------------------------------------------------------
        !% Declarations     
            integer :: ii, jj, kk
            integer :: nflip(1)
            real(8), allocatable :: Darray(:), Aarray(:), Warray2(:,:)
            real(8), pointer     :: depthU(:,:), areaForDepthU(:,:)
            real(8), pointer     :: hydradForDepthU(:,:),widthForDepthU(:,:)
            real(8), pointer     :: areaU(:,:), depthForAreaU(:,:)
            real(8), pointer     :: hydradForAreaU(:,:),widthForAreaU(:,:)
            real(8), pointer     :: depthFull(:), areaFull(:), hydradFull(:)
            real(8), pointer     :: widthFull(:), widthMax(:)
            real(8), pointer     :: areaBelowBreadthMax(:), depthAtBreadthMax(:)
            real(8)              :: depthIncrement, areaIncrement
        !%------------------------------------------------------------------
        !% Preliminaries
            
            if (setting%SWMMinput%N_transect < oneI) return

            !% --- get the number of depth levels in the SWMM-C transect tables from EPA-SWMM
            setting%SWMMInput%N_transect_depth_items = interface_get_N_TRANSECT_TBL()

            !% --- match the transect # in EPA-SWM in SWMM5+
            !%     HACK, we may want make these independent in the future
            N_transect_depth_items  = setting%SWMMInput%N_transect_depth_items
            N_transect_area_items   = setting%SWMMInput%N_transect_depth_items

            !% local arrays
            allocate(Darray(N_transect_depth_items))
            Darray(:) = (/ (jj,jj=zeroI,N_transect_depth_items-oneI)/)

            allocate(Aarray(N_transect_area_items))
            Aarray(:) = (/ (jj,jj=zeroI,N_transect_area_items-oneI)/)

            !% --- a temporary width array
            allocate(Warray2(setting%SWMMinput%N_transect, N_transect_depth_items))

            !% --- allocation transect storage
            call util_allocate_link_transect()

        !%------------------------------------------------------------------
        !% Alias:
            !% --- tables based on uniform discretization of depth
            depthU              => link%transectTableDepthR(:,:,tt_depth)
            areaForDepthU       => link%transectTableDepthR(:,:,tt_area)
            hydradForDepthU     => link%transectTableDepthR(:,:,tt_hydradius)
            widthForDepthU      => link%transectTableDepthR(:,:,tt_width)

            !% tables based on uniform discretization of area
            areaU               => link%transectTableAreaR(:,:,tt_area)
            depthForAreaU       => link%transectTableAreaR(:,:,tt_depth)
            hydradForAreaU      => link%transectTableAreaR(:,:,tt_hydradius)
            widthForAreaU       => link%transectTableAreaR(:,:,tt_width)

            areaFull            => link%transectR(:,tr_areaFull)
            depthFull           => link%transectR(:,tr_depthFull)
            widthMax            => link%transectR(:,tr_widthMax)
            widthFull           => link%transectR(:,tr_widthFull)
            hydradFull          => link%transectR(:,tr_hydRadiusFull)
            areaBelowBreadthMax => link%transectR(:,tr_areaBelowBreadthMax)
            depthAtBreadthMax   => link%transectR(:,tr_depthAtBreadthMax)
        !%------------------------------------------------------------------

        !% --- get transect scalar date from EPA-SWMM
        do ii=1,setting%SWMMinput%N_transect
            link%transectI(ii,ti_idx) = ii
            link%transectR(ii,tr_depthFull)         = interface_get_transectf_attribute(ii,api_transectf_yFull)
            link%transectR(ii,tr_areaFull)          = interface_get_transectf_attribute(ii,api_transectf_aFull)
            link%transectR(ii,tr_hydRadiusFull)     = interface_get_transectf_attribute(ii,api_transectf_rFull)
            link%transectR(ii,tr_widthMax)          = interface_get_transectf_attribute(ii,api_transectf_wMax)
            link%transectR(ii,tr_depthAtBreadthMax) = interface_get_transectf_attribute(ii,api_transectf_ywMax)
            link%transectR(ii,tr_sectionFactor)     = interface_get_transectf_attribute(ii,api_transectf_sMax)
            link%transectR(ii,tr_areaAtMaxFlow)     = interface_get_transectf_attribute(ii,api_transectf_aMax)
            link%transectR(ii,tr_lengthFactor)      = interface_get_transectf_attribute(ii,api_transectf_lengthFactor)
            link%transectR(ii,tr_ManningsN)         = interface_get_transectf_attribute(ii,api_transectf_roughness)
        end do

        !% --- get the transect ID names from EPA-SWMM
        call interface_update_transectID_names ()

        !% --- get the transect table info from EPA-SWMM (normalized)
        !%     this gets Transect.aTbl, Transect.hradTbl Transect.widthTbl Transect.areaTbl
        !%     and stores in the transectTableDepthR() array
        call interface_get_transect_table ()

        !% --- convert EPA-SWMM normalized tables to dimensional values for processing here
        areaForDepthU   = areaForDepthU   * spread(areaFull,  twoI,N_transect_depth_items)
        widthForDepthU  = widthForDepthU  * spread(widthMax,  twoI,N_transect_depth_items)
        hydradForDepthU = hydradForDepthU * spread(hydradFull,twoI,N_transect_depth_items)

        !% --- compute additional transect data
        !%     EPA-SWMM only stores width, area, and hydraulic radius for irregular
        !%     cross-sections with uniform depth discretization. SWMM5+ stores
        !%     the actual depth discretization as well.
        do ii=1,setting%SWMMinput%N_transect
            !% --- EPA-SWMM assigns linear depth increment 
            !%     see function transect_validate in transect.c
            !%     Compute and store the uniformly-discretized depth
            depthIncrement = depthFull(ii) / ( real(N_transect_depth_items-1,8) )
            depthU(ii,:) = depthIncrement * Darray

            !% --- full width is not stored by EPA-SWMM, but is needed for
            !%     SWMM5+. This value is TTable width data at largest depth value
            widthFull(ii) = widthForDepthU(ii,N_transect_depth_items)               
        end do

        !% ------------------------------------------------------------------------
        !% --- BUGFIX -- REMOVE THIS IF EPA-SWMM UPDATED FOR Transect.ywMax
        !%     EPA-SWMM does not store the ywMax in Transect[index].ywMax
        !%     See xsect_setIrregXsectParams() in xsect.c
        !%     The value of ywMax is computed for each cross-section, but
        !%     is not stored as part of the Transect object.
        !%     Possible fix in EPA-SWMM is: Transect[index].ywMax = xsect->ywMax;
        !%     Here we provide a separate computation of yMax, which gets stored
        !%     in transectR(:,tr_depthAtBreadthMax)

            !% --- compute the change in width for each depth increment
            Warray2 = widthForDepthU(:,2:N_transect_depth_items) - widthForDepthU(:,1:N_transect_depth_items-1)

            !% --- identify negative width increments (decreasing width)
            where (Warray2 < zeroR)
                Warray2 = -oneR
            elsewhere
                Warray2 = oneR
            endwhere
            !% --- find the depth at max breadth
            do ii=1,setting%SWMMinput%N_transect
                !% --- get the index for the flip from increasing to decreasing width
                nflip = findloc(Warray2(ii,:),-oneR)
                !% --- if width never decreases, use max depth
                if (nflip(1) .eq. zeroI) nflip = N_transect_depth_items
                !% --- set the depth from the depth table
                depthAtBreadthMax(ii) = depthU(ii,nflip(1)) 
            end do
        !% END BUGFIX
        !% ------------------------------------------------------------------------

        !% --- update the area below the maximum width
        do ii=1,setting%SWMMinput%N_transect
            !% --- lookup the area below the max breadth
            areaBelowBreadthMax(ii) = xsect_table_lookup_singular ( &
                                    depthAtBreadthMax(ii), areaForDepthU(ii,:))
            !print *, 'areabelowbreadthmax ',areaBelowBreadthMax(ii)
        end do

        

        !% --- compute the uniform area transect tables by using the
        !%     nonuniform area as a lookup table for the uniform depth
        do ii=1,setting%SWMMinput%N_transect
            !% --- get the uniform increments of area
            areaIncrement = areaFull(ii) / ( dble(N_transect_area_items-1) )
            areaU(ii,:) = areaIncrement * Aarray

            !% --- call the lookups for the first level with error checking
            depthForAreaU(ii,1) = xsect_nonuniform_lookup_singular &
                (areaU(ii,1), areaForDepthU(ii,:), depthU(ii,:), .true.)

            hydradForAreaU(ii,1) = xsect_nonuniform_lookup_singular &
                (depthForAreaU(ii,1),depthU(ii,:), hydradForDepthU(ii,:),.true. )

            widthForAreaU(ii,1) = xsect_nonuniform_lookup_singular &
                    (depthForAreaU(ii,1), depthU(ii,:), widthForDepthU(ii,:),.true. )

            !% --- cycle through the uniform area values in this transect to 
            !%     find the associated depth     
            do kk=2,N_transect_area_items
                depthForAreaU(ii,kk) = xsect_nonuniform_lookup_singular &
                    (areaU(ii,kk), areaForDepthU(ii,:), depthU(ii,:), .false.)
            end do

            !% --- using the new depthForAreaU, find the hydradius at those depths
            do kk=2,N_transect_area_items
                hydradForAreaU(ii,kk) = xsect_nonuniform_lookup_singular &
                    (depthForAreaU(ii,kk), depthU(ii,:), hydradForDepthU(ii,:),.false. )
            end do

            !% --- using the new depthForAreaU, find the width at those depths
            do kk=2,N_transect_area_items
                widthForAreaU(ii,kk) = xsect_nonuniform_lookup_singular &
                    (depthForAreaU(ii,kk), depthU(ii,:), widthForDepthU(ii,:),.false. )
            end do
   
        end do

        !% --- Normalize all tables
            depthU          = depthU          / spread(depthFull  ,2,N_transect_depth_items)
            areaForDepthU   = areaForDepthU   / spread(areaFull   ,2,N_transect_depth_items)
            hydradForDepthU = hydradForDepthU / spread(hydradFull ,2,N_transect_depth_items)
            widthForDepthU  = widthForDepthU  / spread(widthMax   ,2,N_transect_depth_items)

            areaU           = areaU           / spread(areaFull   ,2,N_transect_area_items)
            depthForAreaU   = depthForAreaU   / spread(depthFull  ,2,N_transect_area_items)
            hydradForAreaU  = hydradForAreaU  / spread(hydradFull ,2,N_transect_area_items)
            widthForAreaU   = widthForAreaU   / spread(widthMax   ,2,N_transect_area_items)

        !%------------------------------------------------------------------
        !% Closing:
            deallocate(Darray)
            deallocate(Aarray)
            deallocate(Warray2)

    end subroutine init_link_transect_array
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
        !% Declarations:
            integer :: ii
            integer, dimension (setting%SWMMinput%N_subcatch) :: runon_count
            integer, pointer :: subRunon
        !%------------------------------------------------------------------
        !% branchsign global is used for junction branches (JB)
        !%     for upstream (+1) and downstream (-1)
            branchsign(1:max_branch_per_node-1:2) = +oneR
            branchsign(2:max_branch_per_node:2)   = -oneR
        !%------------------------------------------------------------------
        !% Subcatchments need additional indexes if RunOns exist
            if (setting%Simulation%useHydrology) then
                runon_count(:) = zeroI
                do ii=1,N_node
                    !% --- only BCdn nodes are outfalls and  canhave RunOn to catchment
                    if (node%I(ii,ni_node_type) .ne. nBCdn ) cycle
                    !print *, 'routeTo for ',ii,' is ', node%I(ii,ni_routeTo)
                    !% --- identify the subcatchment that this node runs onto
                    subRunon =>  node%I(ii,ni_routeTo)

                    !% --- a null value means the outfall is not to a subcatchment
                    if (node%I(ii,ni_routeTo) .eq. nullvalueI) cycle 

                    !% --- count the number of routed elements to each catchment
                    if (subRunon > setting%SWMMinput%N_subcatch) then
                        print *, 'CODE ERROR mismatch in subcatchment count'
                        call util_crashpoint(609873)
                    else
                        !% --- increment the counter of runons to each subcatchment
                        runon_count(subRunon) = runon_count(subRunon)+1
                    end if
                end do

                !% --- get the maximum number of runon to a single subcatchment
                N_subcatch_runon = maxval(runon_count)
                if (N_subcatch_runon > 0) then
                    !% --- initialize the subcatchment runon
                    call util_allocate_subcatch_runon ()

                else 
                    !% --- set size of arrays to 1 and set values to nullvalueI
                    !%     If any of these are called as column index, we expect
                    !%     to get a segmentation fault
                    allocate(si_RunOn_nodeIdx(1))
                    allocate(si_RunOn_SWMMoutfallIdx(1))
                    allocate(si_RunOn_faceIdx(1))
                    allocate(si_RunOn_P_image(1))
                    allocate(sr_RunOn_Volume(1))
                    si_RunOn_SWMMoutfallIdx = nullvalueI
                    si_RunOn_nodeIdx        = nullvalueI
                    si_RunOn_faceIdx        = nullvalueI
                    si_RunOn_P_image        = nullvalueI
                    sr_RunOn_Volume         = nullvalueI
                end if
            end if

        !% Closing
    end subroutine init_globals
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_curves()
        !%------------------------------------------------------------------
        !% Description:
        !%   Retrieves data from EPA-SWMM interface and populates curve curves
        !%------------------------------------------------------------------
        !% Declarations
            integer       :: ii, jj, additional_storage_curves
            character(64) :: subroutine_name = 'init_curves'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

            if (.not. api_is_initialized) then
                print *, "CODE ERROR API is not initialized"
                call util_crashpoint(84785)
            end if
        !%------------------------------------------------------------------
        !% --- we create additional curves for functional storage as well
        !%     this allocates the space for functional storage curve
        additional_storage_curves = count((node%YN(:,nYN_has_storage)) .and. &
                                          (node%I(:,ni_curve_ID) == 0))

        !% --- assign the global total number of storage curves
        N_Total_Curves = additional_storage_curves + setting%SWMMinput%N_curve

        !% --- allocate the number of curve objects from SWMM5
        call util_allocate_curves()
   
        do ii = 1, setting%SWMMinput%N_curve
            curve(ii)%ID = ii
            curve(ii)%Type = interface_get_table_attribute(ii, api_table_type)
            
            !% --- get the number of entries in a curve
            curve(ii)%NumRows = interface_get_num_table_entries(ii)
            
            !% --- allocate the value space
            call util_allocate_curve_entries (ii,curve(ii)%NumRows)

            !% --- get the first entry of the curve and store them at the first two columns
            curve(ii)%ValueArray(1,curve_read_col_1:curve_read_col_2) = interface_get_first_entry_table(ii)

            !% --- populate the rest of the curves
            do jj = 2,curve(ii)%NumRows
                curve(ii)%ValueArray(jj,curve_read_col_1:curve_read_col_2) = interface_get_next_entry_table(ii, API_CURVE)
            end do
            
        end do

        if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_curves
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_simulation_controls()
        !%------------------------------------------------------------------
        !% Description:
        !%   initialize simulation controls from SWMM5+ settings file
        !%------------------------------------------------------------------
        !% Declarations
            integer              :: ii, jj, lIdx
            integer, pointer     :: nControls
            integer, allocatable :: packedElemArray(:)
            character(64)        :: subroutine_name = 'init_simulation_controls'
        !%------------------------------------------------------------------
        !% Aliases:        
            !% pointer to the number of control in the settings file
            nControls => setting%Control%NumControl
        
        !%------------------------------------------------------------------
        !% Preliminaries
            if (nControls < oneI) return !% no controls

            if (setting%Debug%File%initialization) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------

        allocate(setting%Control%LinkIdx(nControls))
        setting%Control%LinkIdx = nullvalueI

        allocate(setting%Control%ElemIdx(nControls))
        setting%Control%ElemIdx = nullvalueI

        !% --- cycle through the controls to connect controlled EPA SWMM link to
        !%     SWMM5+ elements.
        do ii = 1,nControls
            do jj = 1,size(link%I,dim=1)
                if (link%Names(jj)%str == trim(setting%Control%Links(ii))) then
                    lIdx = jj
                    setting%Control%LinkIdx(ii) = lIdx
                    packedElemArray = pack(elemI(:,ei_Lidx),(elemI(:,ei_link_Gidx_SWMM) == lIdx))
                    if (size(packedElemArray)<oneI) then
                        print*, packedElemArray, 'packedElemArray'
                        print*, "Error - json file - setting " // 'Could not find the link ', trim(setting%Control%Links(ii)), ' in the network'
                    end if
                    !% Take the first element of the subsequent link as controlled element.
                    !% For conduits the first element will be controlled
                    !% for orifice and weirs it will not matter
                    setting%Control%ElemIdx(ii) = packedElemArray(1)
                    deallocate(packedElemArray)
                end if
            end do
            !% check for valid links in the control settings
            if (setting%Control%LinkIdx(ii) == nullvalueI) then
                print*, "Error - json file - setting " // 'Could not find then link ', trim(setting%Control%Links(ii)), ' in the network'
                call util_crashpoint(6298733)
            end if
        end do 

        !%------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_simulation_controls
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_profiles
        !%------------------------------------------------------------------
        !% Description:
        !% sets up SWMM5+ for output of link profiles that are designated
        !% in the EPA SWMM *.inp file
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: thisUnit, file_loc
            logical :: doesExist, isOpen
        !%------------------------------------------------------------------
        !% --- alias for the unit number
        thisUnit = setting%File%UnitNumber%inp_file
        file_loc = 0

        !% --- check that file exists
        inquire(unit=thisUnit,EXIST=doesExist)
        if (.not. doesExist) then 
            write(*,"(A)")   'USER CONFIGURATION ERROR OR CODE ERROR trying to open *.inp file that does not exist'
            write(*,"(A,A)") 'Filename: ',trim(setting%File%inp_file)
            call util_crashpoint(5509877)
        end if

        !% --- file should NOT be open at this time, if not there is a code error
        inquire(FILE=trim(setting%File%inp_file),OPENED=isOpen)
        if (isOpen) then 
            write(*,"(A)") 'CODE ERROR input file should not be open when processing profiles'
            write(*,"(A,A)") 'Filename: ',trim(setting%File%inp_file)
            call util_crashpoint(3385782)
        end if
        
        open(thisUnit, file = trim(setting%File%inp_file), STATUS = 'OLD', ACTION = 'READ')

        call init_profiles_count (thisUnit, file_loc)

        !% --- set the maximum number of items in a profile
        !%     assumes every link is separated by one node
        N_items_in_profile = (N_links_in_profile*twoI) + oneI 

        !% --- exit if no profiles were found
        if (N_Profiles == 0) then
            print *, "... no profiles found"
            return
        endif

        !% --- allocate storage for the profiles
        call util_allocate_output_profiles()

        !% --- store the output_profile_link_names and output_profile_names
        call init_profiles_read (thisUnit, file_loc)

        close (thisUnit)

        !% --- store the output profile link idexes
        call init_profiles_get_link_idx ()

        !% --- add nodes between the profile links and check continuity
        call init_profiles_add_nodes ()
   
    end subroutine init_profiles
!%
!%==========================================================================
!%==========================================================================
!%   
    subroutine init_profiles_count (thisUnit, file_loc)
        !%------------------------------------------------------------------
        !% Description
        !% Counts the number of profiles and the maximum number of
        !% links in any profile (which are globals)
        !%------------------------------------------------------------------
            integer, intent(in)    :: thisUnit !% open file for reading
            integer, intent(inout) :: file_loc !% loc in file for profile header
            character(max_names_string_length*N_links_on_profile_line+N_char_of_profile_name) &
                :: lineRead, nameRead, profile_name, profile_name_last, link_names

            integer :: delimiter_loc, nLinks, mmLine, temploc
            integer :: read_status, index_of_start
            logical :: EndOfProfile  = .false.
            logical :: first_profile = .true.
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        delimiter_loc = 2

        !% --- initialize global counters
        N_links_in_profile = 0
        N_profiles = 0

        !% --- initialize decision variables
        first_profile     = .true.
        profile_name_last = "NULL"

        !% --- initialize local coutners
        nLinks = 0

        !% --- Read through the input file looking for the profiles 
        !%     This first read simply counts the profiles and finds the
        !%     maximum number of links in the largest profile
        mmLine=0 !% line counter, used for debugging only

        EndOfProfile = .false.
        do while (.not. EndOfProfile) !% --- outer loop
            !% --- cycle through the input file
            mmLine=mmLine+1 

            !% --- store start of this line
            temploc = FTELL(thisUnit)
            !% --- read line of input file
            read(thisUnit, "(A)", iostat = read_status) lineRead
            if (read_status /= 0) then 
                !% -- end of file or read error
                exit
            endif
            lineRead = ADJUSTL(lineRead)

            if(lineRead .eq. "[PROFILES]") then
                !% --- Set file_loc at the [PROFILES] header
                !%     Allows later rewind of the file to this offset location 
                if(first_profile) then 
                    file_loc = temploc
                    first_profile = .false.
                end if

                !% --- read the next line below the header       
                read(thisUnit, "(A)", iostat = read_status) lineRead
                EndOfProfile = util_read_blankline_or_EOF (read_status, lineRead)
                if (EndOfProfile) exit !% --- break outer do loop without any profiles defined

                !% --- remove any leading blanks
                lineRead = ADJUSTL(lineRead)

                !% --- cycle for comment lines
                do while (index(lineRead,';') == 1)
                    !% --- read the next line        
                    read(thisUnit, "(A)", iostat = read_status) lineRead
                    EndOfProfile = util_read_blankline_or_EOF (read_status, lineRead)
                    if (EndOfProfile) exit !% --- break outer do loop without any profiles defined
                    !% --- remove any leading blanks
                    lineRead = ADJUSTL(lineRead)
                end do

                do while (.not. EndOfProfile) !% --- inner loop

                    !% --- first item on the line is the name of the profile
                    !%     which is in quotations, so the location of the
                    !%     start of the string is 2
                    nameRead      = lineRead 
                    delimiter_loc = 2
    
                    !% --- subdividing nameRead for profile name and link names 
                    !%     find end of the profile name (last quote mark)
                    !%     true = backwards (from end of string)
                    index_of_start = index(nameRead,"""",.true.)
                    !% --- store string with the profile name (includes quotation marks)
                    profile_name   = trim(ADJUSTL(nameRead(1:index_of_start+1)))
                    !% --- store string of only the link names
                    link_names     = trim(ADJUSTL(nameRead(index_of_start+1:len(nameRead))))

                    if (profile_name_last .eq. "NULL") then 
                        !% --- this is the first profile
                        N_Profiles = 1
                        nLinks = 0
                    else
                        !% --- second or later profile line
                        if (profile_name .eq. profile_name_last) then 
                            !% --- continuing the same profile
                            !%     nLinks will accumulate
                        else 
                            !% --- adding a new profile (prior profile done)
                            !% --- check if prior profile had max links
                            N_links_in_profile = max(N_links_in_profile,nLinks)
                            N_Profiles = N_Profiles+1
                            !% --- reset the link counter
                            nLinks = 0
                        end if
                    end if

                    !% --- count the number of links in this profile by iteratively looking for
                    !%     the blank delimiter and removing the data before it
                    do while (delimiter_loc > 1)
                        !% --- find next blank character between link names
                        delimiter_loc = index(link_names," ")
                        if (delimiter_loc .le. 1) exit
                        !% --- remove the preceding string data, leading blanks, and trailing blanks
                        link_names = trim(ADJUSTL(link_names(delimiter_loc+1:)))
                        nLinks = nLinks + 1
                    end do
                    !% --- check for new maximum
                    N_links_in_profile = max(N_links_in_profile,nLinks)

                    !% --- store this profile name as the last one read
                    profile_name_last = profile_name

                    !% --- read the next line
                    read(thisUnit, "(A)", iostat = read_status) lineRead
                    EndOfProfile = util_read_blankline_or_EOF (read_status, lineRead)
                    if (EndOfProfile) exit !% --- break inner and outer loops on blank or EOF
                    !% --- remove any leading blanks
                    lineRead = ADJUSTL(lineRead)

                    !% --- cycle for comment lines
                    do while (index(lineRead,';') == 1)
                        !% --- read the next line        
                        read(thisUnit, "(A)", iostat = read_status) lineRead
                        EndOfProfile = util_read_blankline_or_EOF (read_status, lineRead)
                        if (EndOfProfile) exit !% --- break loops
                        !% --- remove any leading blanks
                        lineRead = ADJUSTL(lineRead)
                    end do

                end do !% --- inner loop
                if (EndOfProfile) exit  !% --- break outer loop
            else
                !% no profiles on this line
            endif
            
        end do !% --- outer loop

        !% allocate the number of links per profile
        allocate(Number_of_links_per_profile(N_Profiles),stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'Number_of_links_per_profile')
        Number_of_links_per_profile(:) = nullValueI

    end subroutine init_profiles_count
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_profiles_read (thisUnit, file_loc)
        !%------------------------------------------------------------------
        !% Description
        !% Reads in the links for the profiles and stores the 
        !% output_profile_link_names and output_profile_names
        !%------------------------------------------------------------------
            integer, intent(in) :: thisUnit, file_loc

            character(max_names_string_length * N_links_on_profile_line + N_char_of_profile_name) &
                    :: lineRead, this_name, profile_name_last, link_names, profile_name, &
                       this_link
            integer :: ierror, read_status, iprof, jlink, delimiter_loc
            integer :: index_of_start
            logical :: EndOfProfile = .false.
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        
        !% --- initialize decision variables
        profile_name_last = "NULL"
        
        !% --- rewind this unit and reset to [PROFILE] header line
        rewind(thisUnit)
        ierror = fseek(thisUnit, file_loc, 0)

        !% --- if position not found, then exit
        if (ierror /= 0) then 
            print *, '... [PROFILE] location not found in *.inp file'
            return
        end if
        
        !% --- begin reading above [PROFILE] header
        read(thisUnit, "(A)", iostat = read_status) lineRead
        do while (lineRead == ' ')
            read(thisUnit, "(A)", iostat = read_status) lineRead
        end do
        if(lineRead .eq. "[PROFILES]") then

            !% --- read the next line below the header       
            read(thisUnit, "(A)", iostat = read_status) lineRead
            EndOfProfile = util_read_blankline_or_EOF (read_status, lineRead)
            if (EndOfProfile) return !% --- break outer do loop without any profiles defined

            !% --- remove any leading blanks
            lineRead = ADJUSTL(lineRead)

             !% --- cycle for comment lines
            do while (index(lineRead,';') == 1)
                !% --- read the next line        
                read(thisUnit, "(A)", iostat = read_status) lineRead
                EndOfProfile = util_read_blankline_or_EOF (read_status, lineRead)
                if (EndOfProfile) exit !% --- break outer do loop without any profiles defined
                !% --- remove any leading blanks
                lineRead = ADJUSTL(lineRead)
            end do

            iprof = 0
            jlink = 0
            do while (.not. EndOfProfile) 

                !% --- first item on the line is the name of the profile
                !%     which is in quotations, so the location of the
                !%     start of the string is 2
                delimiter_loc = 2

                !% --- subdividing nameRead for profile name and link names 
                !%     find end of the profile name (last quote mark)
                !%     true = backwards (from end of string)
                index_of_start = index(lineRead,"""",.true.)
                !% --- store string with the profile name (includes quotation marks)
                profile_name   = trim(ADJUSTL(lineRead(1:index_of_start+1)))
                !% --- store string of only the link names
                link_names     = trim(ADJUSTL(lineRead(index_of_start+1:len(lineRead))))

                !% --- store the trimmed profile name
                if (profile_name .eq. profile_name_last) then 
                    !% --- this is a continuing profile
                    !%     name is already stored
                    !%     do not change iprof
                else
                    !% --- starting a new profile
                    !%     get the profile name without quotation marks
                    iprof = iprof+1
                    !% reset the jlink counter
                    jlink = 0
                    !% --- check profile counter
                    if (iprof > size(output_profile_link_names,1)) then 
                        print *, 'CODE ERROR profile counter is too small'
                        print *, 'storage is ',size(output_profile_link_names,1)
                        print *, 'iprof = ',iprof
                        call util_crashpoint(71197)
                        return
                    end if

                    !% --- trim the name of quotation marks
                    this_name = trim(profile_name(2:index(profile_name,'"',.true.)-1))

                    !% --- check length before storing
                    if (len_trim(this_name) .le. stringLength_HDF5) then
                        output_profile_names(iprof) = trim(this_name)
                    else
                        print *, 'USER CONFIGURATION ERROR profile name is too long'
                        print *, 'SWMM5+ value for stringLength_HDF5 is ',stringLength_HDF5 
                        print *, 'profile has name ',trim(this_name)
                        print *, 'which is of length ',len_trim(this_name)
                        call util_crashpoint(610874)
                        return
                    end if
                    !% --- set last name for checking continuing profile
                    profile_name_last = profile_name
                end if

                !% --- cycle through the link names to identify links
                !% --- find next blank character between link names
                delimiter_loc = index(link_names," ")
                do while (delimiter_loc > 1)
                    !% --- increment link counter
                    jlink = jlink + 1
                    !% --- update the number of links per profile 
                    Number_of_links_per_profile(iprof) = jlink
                    !% --- check link counter
                    if (jlink > size(output_profile_link_names,2)) then 
                        print *, 'CODE ERROR profile link counter is too small'
                        print *, 'storage is ',size(output_profile_link_names,2)
                        print *, 'jlink = ',jlink 
                        call util_crashpoint(71198)
                        return
                    end if
                    !% --- get this link name
                    this_link = trim(ADJUSTL(link_names(1:delimiter_loc)))
                    !% --- check length before storing
                    if (len_trim(this_link) .le. stringLength_HDF5) then 
                        output_profile_link_names(iprof,jlink) = trim(this_link)
                    else
                        print *, 'USER CONFIGURATION ERROR link name is too long'
                        print *, 'link name in profile exceeded the SWMM5+ '
                        print *, 'value for stringLength_HDF5 of ',stringLength_HDF5 
                        print *, 'link with name ',trim(this_name)
                        print *, 'is of length ',len_trim(this_name)
                        call util_crashpoint(2209874)
                        return
                    end if
                    !% --- remove this link name from the link name list
                    link_names = trim(ADJUSTL(link_names(delimiter_loc+1:)))
                    !% --- get the next delimiter
                    delimiter_loc = index(link_names," ")    
                    !print *, 'link names ',link_names                
                end do

                !% --- read the next line
                read(thisUnit, "(A)", iostat = read_status) lineRead
                EndOfProfile = util_read_blankline_or_EOF (read_status, lineRead)

            end do

        else
            print *, 'CODE ERROR misalignment of PROFILE reading'
            call util_crashpoint(7209873)
            return
        end if

    end subroutine init_profiles_read
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_profiles_get_link_idx ()
        !%------------------------------------------------------------------
        !% Descriptions
        !% adds nodes between links in the output_profile
        !% checks for link-node-link continuity
        !%------------------------------------------------------------------
        !% Declarations
            integer :: pp, mm, kk
            logical :: foundLast
        !%------------------------------------------------------------------
        !% Preliminaries
            if (N_Profiles < 1) return
        !%------------------------------------------------------------------
        !% 
        !% --- cycle through the profiles
        do pp =1,N_profiles
            !% --- cycle through links of the system
            do kk = 1,N_link
                !% --- look for links in the profile
                if (any(output_profile_link_names .eq. trim(link%names(kk)%str))) then 
                    do mm = 1,N_links_in_profile
                        if (output_profile_link_names(pp,mm) .eq.  trim(link%names(kk)%str)) then 
                            !% --- store the link index
                            output_profile_link_idx(pp,mm) = link%I(kk,li_idx)
                            ! print *, 'this link in profile'
                            ! print *, pp,mm, output_profile_link_idx(pp,mm), trim(output_profile_link_names(pp,mm))
                        else
                            !% --- kk link is not a profile link
                        end if
                    end do !% profile links
                end if
            end do !% links
        end do !% pp profiles

        !% --- check for link names that do not exist
        do pp = 1, N_profiles
            do mm = 1,Number_of_links_per_profile(pp) 
                if (( trim(output_profile_link_names(pp,mm)) .ne. "NULL") &
                    .and.                                                        &
                    (output_profile_link_idx(pp,mm) .eq. nullvalueI)) then 
                    !% --- profile name not found
                    print *, 'USER CONFIGURATION ERROR Profile link name not found'
                    print *, 'The profile includes link name ', trim(output_profile_link_names(pp,mm))
                    print *, 'which was not found in the system links of the *.inp file '
                    call util_crashpoint(2297445)
                    return
                else
                    !% --- no problem
                end if
            end do
        end do

        !% --- check for continuity in storage and that there are at least 2 links
        !%     in a profile
        do pp = 1,N_Profiles
            foundLast = .false.
            do mm = 1,N_links_in_profile
                if (output_profile_link_idx(pp,mm) .ne. nullvalueI) then
                    !% --- profile has index
                    if (foundLast) then 
                        !% --- an index should not occur after foundLast is true
                        print *, 'CODE ERROR in output profiles'
                        print *, 'The output_profile_link_idx(:,:) is not contiguous'
                        print *, 'for profile ',trim(output_profile_names(pp))
                        print *, 'at link ', trim(output_profile_link_names(pp,mm))
                        call util_crashpoint(2209874)
                        return
                    else 
                        !% --- no problems
                    end if
                else
                    foundLast = .true.
                    if (mm < 2) then 
                        print *, 'USER CONFIGURATION ERROR in output profiles'
                        print *, 'The output profile does not have a least two valid links'
                        print *, 'for profile ',trim(output_profile_names(pp))
                        call util_crashpoint(5198733)
                        return
                    else
                        !% --- no problem 
                    end if
                end if
            end do
        end do

    end subroutine init_profiles_get_link_idx   
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_profiles_add_nodes ()
        !%------------------------------------------------------------------
        !% Descriptions
        !% adds nodes between links in the output_profile_ids
        !% checks for link-node-link continuity
        !%------------------------------------------------------------------
        !% Declarations
            integer :: pp, mm
            integer, pointer :: thisNode, nextNode, testNode, thisLink, nextLink 
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        do pp = 1,N_Profiles
            !% --- store the starting node (upstream of 1st link)
            thisNode => link%I(output_profile_link_idx(pp,1),li_Mnode_u)
            output_profile_node_names(pp,1) = trim(node%Names(thisNode)%str)
            output_profile_ids(pp,1) = thisNode
            !% --- store the starting link (item 2)
            output_profile_ids(pp,2) = output_profile_link_idx(pp,1)
            !% --- cycle through the links
            do mm = 1,Number_of_links_per_profile(pp)-1
                thisLink => output_profile_link_idx(pp,mm)
                nextLink => output_profile_link_idx(pp,mm+1)
                !% --- get the downstream node of thisLink
                thisNode => link%I(thisLink,li_Mnode_d)
                !% --- store the node data
                output_profile_ids(pp,2*mm+1) = thisNode
                output_profile_node_names(pp,mm+1) = trim(node%Names(thisNode)%str)
                !% --- check if next link exists
                if (nextLink .ne. nullvalueI) then
                    !% --- get the upstream node of nextLink
                    testNode => link%I(nextLink,li_Mnode_u)
                    !% --- check that the next link has the same upstream node
                    if (thisNode .ne. testNode) then 
                        print *, 'USER CONFIGURATION ERROR in output profile'
                        print *, 'Profile is not contiguous'
                        print *, 'Upstream link          ',trim(output_profile_link_names(pp,mm))
                        print *, '...has downstream node ',trim(node%Names(thisNode)%str)
                        print *, 'Downstream link        ',trim(output_profile_link_names(pp,mm+1))
                        print *, '... has upstream node  ',trim(node%Names(testNode)%str)
                        print *, 'The correct downstream link should be one that connects to'
                        print *, '... the downstream node.'
                        call util_crashpoint(7220987)
                    else
                        !% --- store next link
                        output_profile_ids(pp,2*mm+2) = nextLink
                        !% --- store next node
                        nextNode => link%I(nextLink,li_Mnode_d)
                        output_profile_ids(pp,2*mm+3) = nextNode     
                        output_profile_node_names(pp,mm+2) = trim(node%Names(nextNode)%str)
                    end if
                else
                    !% --- null value found for next link, profile is done
                    cycle !% mm
                endif
            end do !% mm profile links
        end do !% pp profiles   

        !% --- test output SAVE FOR DEBUGGING
        ! do pp = 1,N_profiles
        !     print *, 'Profile Name:',trim(output_profile_names(pp))
        !     !print *, output_profile_ids(pp,:)
        !     do mm=1,N_items_in_profile
        !         !print *, '      ', mm,  output_profile_ids(pp,mm)
        !         if (output_profile_ids(pp,mm) .ne. nullvalueI) then
        !             if (mod(mm,2) .eq. 0) then 
        !                 print *, 'link: ',mm, output_profile_ids(pp,mm), trim(output_profile_link_names(pp,1+((mm-1)/2)))
        !             else
        !                 print *, 'node: ', mm, output_profile_ids(pp,mm), trim(output_profile_node_names(pp,1+((mm-1)/2)))
        !             end if
        !         end if
        !     end do
        ! end do

        !% --- check for nodes and links that have the same names
        do pp = 1,N_profiles
            do mm = 1,Number_of_links_per_profile(pp)
                if (any(output_profile_node_names(pp,:) .eq. trim(output_profile_link_names(pp,mm)))) then 
                    print *, 'USER CONFIGURATION ERROR found link and node with identical names'
                    print *, 'To prevent confusion, SWMM5+ requires links and nodes to have '
                    print *, 'unique names in profiles. Problem found with '
                    print *, 'link and node named:', trim(output_profile_link_names(pp,mm))
                    call util_crashpoint(5109873)
                    return
                end if
            end do
        end do

    end subroutine init_profiles_add_nodes   
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine init_partitioning()
        !%------------------------------------------------------------------
        !% Description:
        !%   This subroutine calls the public subroutine from the utility module,
        !%   partitioning.f08. It also calls a public subroutine from the temporary
        !%   coarray_partition.f08 utility module that defines how big the coarrays
        !%   must be.
        !%
        !%------------------------------------------------------------------
            integer       :: ii
            character(64) :: subroutine_name = 'init_partitioning'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%initialization) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

            !% if there are no links, the system cannot be partitioned
            if (N_link == 0) then
                if (this_image() == 1) then
                    write(*,*) '******************************************************'
                    write(*,*) '*          USER CONFIGURATION ERROR                  *'
                    write(*,*) '* The SWMM input file does not include any links.    *'
                    write(*,*) '* The SWMM5+ code requires at least one link to run. *'
                    write(*,*) '* This run was stopped without any output.           *'
                    write(*,*) '******************************************************'
                end if
                call util_crashpoint(970532)
                return
            end if 
        
            if (setting%Profile%useYN) call util_profiler_start (pfc_init_partitioning)
        !%------------------------------------------------------------------

        !% --- find the number of elements in a link based on nominal element length
        do ii = 1, setting%SWMMinput%N_link
            call init_discretization_nominal(ii)
        end do

        !% --- Set the network partitioning method used for multi-processor parallel computation
        call partitioning_toplevel()
        sync all

        !% HOLD FOR FUTURE
        !% --- Compute the amount of a conduit length that is added to a connected junction.
        !%     This modifies the conduit length itself if setting%Discretization%AdjustLinkLengthForJunctionBranchYN
        !%     is true. The junction itself is setup in init_network_nJm_branch_length()
        !call init_discretization_adjustlinklength()

        !% --- calculate the largest number of elements and faces to allocate the coarrays
        call init_coarray_length()

        !% --- allocate elem and face coarrays
        call util_allocate_elemX_faceX()
        call util_key_default_elemX()
        call util_key_default_face()

        !% --- allocate column indexes of elem and face arrays for pointer operation
        call util_allocate_columns()

        !%------------------------------------------------------------------
        !% Closing
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
        !% HACK -- not sure how this functions if a node appears on more
        !% than one image
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: ii
            integer, pointer :: nodeIdx(:), elemIdx(:), nodeType(:)
            integer, pointer :: tface, thisSub, thisRunon
            character(64) :: subroutine_name = 'init_subcatchment'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%initialization) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
            nodeIdx  => subcatchI(:,si_runoff_nodeIdx) 
            elemIdx  => subcatchI(:,si_runoff_elemIdx)
            nodeType => node%I(:,ni_node_type)
        !%------------------------------------------------------------------

        !% --- cycle through subcatchments to set connections to runoff nodes in EPA SWMM
        do ii=1,setting%SWMMinput%N_subcatch
                
            if (nodeIdx(ii) .eq. nullvalueI) cycle
          
            !% --- the runoff image (partition) is the node image (partition)
            subcatchI(ii,si_runoff_P_image) = node%I(nodeIdx(ii), ni_P_image)

            !% --- only set the elements and faces for subcatchment on images with 
            !%     its (single) connected runoff
            if (this_image() .eq. node%I(nodeIdx(ii), ni_P_image)) then
                select case (nodeType(nodeIdx(ii)))

                    case (nJ2)
                        print *, 'CODE ERROR'
                        print *, 'AT THIS TIME, A SUBCATCHMENT SHOULD NEVER CONNECT WITH nJ2 JUNCTION'
                        print *, 'Check that all indexes are set correctly'
                        call util_crashpoint(448723)

                        !% for a node that is a face, the subcatch connects to the element
                        !% upstream of the face if that element is CC.
                        tface => node%I(nodeIdx(ii),ni_face_idx) 
                        elemIdx(ii) = faceI(tface,fi_Melem_uL)

                        select case (elemI(elemIdx(ii),ei_elementType))
                            case (CC)
                                !% --- use the elemIdx already computed
                            case (JB)
                                !% --- use the upstream JM
                                elemIdx(ii) = elemSI(elemIdx(ii),esi_JB_Main_Index)
                            case default
                                !% --- switch to the downstream element from the face
                                elemIdx(ii) = faceI(tface,fi_Melem_uL)

                                select case (elemI(elemIdx(ii),ei_elementType))
                                    case (CC)
                                        !% --- use the elemIdx already computed
                                    case (JB)
                                        !% --- use the downstream JM
                                        elemIdx(ii) = elemSI(elemIdx(ii),esi_JB_Main_Index)
                                    case default
                                        print *, 'USER CONFIGURATION ERROR a subcatchment input was defined'
                                        print *, 'to a node without storage, which is an nJ2 face.'
                                        print *, 'The upstream and downstream elements from the face'
                                        print *, 'are diagnostic elements, which cannot take a subcatchment'
                                        print *, 'inflow. Reconfigure so that the subcatchment outflow is into a'
                                        print *, 'storage node or into a node adjacent to a conduit link or' 
                                        print *, 'a junction with storage.'
                                        print *, 'node ID ',nodeIdx(ii)
                                        print *, 'node name ',trim(node%Names(nodeIdx(ii))%str)
                                        call util_crashpoint(4023987)
                                end select
                        end select

                    case (nJm)
                        !% --- for a node that is a multi-branch junction, subcatch connects to 
                        !%     the element itself
                        !%     Note this is only allowed for runoff. Runon requires an outfall
                        !%     which (by limitation of EPA-SWMM) must have only a single link
                        !%     connecte to it.
                        elemIdx(ii) = node%I(nodeIdx(ii), ni_elem_idx)

                    case (nBCup,nJ1)
                        !% --- for a node that is an upstream BC or dead end the subcatch connects 
                        !%     into the first element downstream of the face
                        tface => node%I(nodeIdx(ii),ni_face_idx) 
                        if (tface .ne. nullvalueI) then 
                            elemIdx(ii) = faceI(tface,fi_Melem_dL)

                            select case (elemI(elemIdx(ii),ei_elementType))
                                case (CC)
                                    !% --- use the elemIdx already computed
                                case (JB)
                                    !% --- use the associated JM element
                                    elemIdx(ii) = elemSI(elemIdx(ii),esi_JB_Main_Index)
                                case default
                                    print *, 'USER CONFIGURATION ERROR a subcatchment input was defined'
                                    print *, 'to an upstream (inflow) BC node but the downstream'
                                    print *, 'elements from the node are diagnostic elements, which' 
                                    print *, 'cannot take a subcatchment inflow.'
                                    print *, 'Reconfigure so that the subcatchment outflow is into a'
                                    print *, 'storage node or into a node adjacent to a conduit link or' 
                                    print *, 'a junction with storage.'
                                    print *, 'node ID ',nodeIdx(ii)
                                    print *, 'node name ',trim(node%Names(nodeIdx(ii))%str)
                                    call util_crashpoint(4023987)
                            end select

                        else
                            elemIdx(ii) = nullvalueI
                            print *, 'CODE ERROR, unexpected null downstream element for upstream BC'
                            print *, 'node ID ',nodeIdx(ii)
                            print *, 'node name ',trim(node%Names(nodeIdx(ii))%str)
                            call util_crashpoint(1098223)
                        end if

                    case (nBCdn)
                        !% --- for a node that is an downstreamstream BC, the subcatch connects 
                        !% first element upstreamstream of the face
                        tface => node%I(nodeIdx(ii),ni_face_idx)
                        if (tface .ne. nullvalueI) then 
                            elemIdx(ii) = faceI(tface,fi_Melem_uL)
                        
                            select case (elemI(elemIdx(ii),ei_elementType))
                                case (CC)
                                    !% --- use the elemIdx already computed
                                case (JB)
                                    !% --- use the associated JM element
                                    elemIdx(ii) = elemSI(elemIdx(ii),esi_JB_Main_Index)
                                case default
                                    print *, 'USER CONFIGURATION ERROR a subcatchment input was defined'
                                    print *, 'to an downstream (head) BC node but the upstream'
                                    print *, 'elements from the node are diagnostic elements, which' 
                                    print *, 'cannot take a subcatchment inflow.'
                                    print *, 'Reconfigure so that the subcatchment outflow is into a'
                                    print *, 'storage node or into a node adjacent to a conduit link or' 
                                    print *, 'a junction with storage.'
                                    print *, 'node ID ',nodeIdx(ii)
                                    print *, 'node name ',trim(node%Names(nodeIdx(ii))%str)
                                    call util_crashpoint(4023987)
                            end select

                        else
                            elemIdx(ii) = nullvalueI
                            print *, 'CODE ERROR, unexpected null upstream element for downstream BC'
                            print *, 'node ID ',nodeIdx(ii)
                            print *, 'node name ',trim(node%Names(nodeIdx(ii))%str)
                            call util_crashpoint(1098245)
                        end if

                    case default 
                        print *, ii, nodeIdx(ii)
                        write(*,*) 'CODE ERROR unexpected case default in '//trim(subroutine_name)
                        print *, 'node index ',nodeIdx(ii)
                        print *, 'ID ', trim(node%Names(nodeIdx(ii))%str)
                        print *, 'Node Type # of ',nodeType(nodeIdx(ii))
                        print *, 'which has key of ',trim(reverseKey(nodeType(nodeIdx(ii))))
                        call util_crashpoint(3758974)
                end select

                !% store logical to designated runoff element
                if (elemIdx(ii) .ne. nullvalueI) then 
                    elemYN(elemIdx(ii), eYN_hasSubcatchRunOff) = .true.
                end if
            end if
        end do

        !% --- initialize the runon counter for each subcatchment
        subcatchI (:,si_RunOn_count) = zeroI
        subcatchYN(:,sYN_hasRunOn)   = .false.  

        !% --- cycle through nodes if there are runons to subcatchments.
        !%     Outfalls store the runon connections from nodes to subcatchments 
        if (N_subcatch_runon > 0) then
            !% --- cycle through the nodes
            do ii=1,N_node

                !% --- look for nodes with valid routeTo data
                if  ( (node%I(ii,ni_routeTo) > 0)                               &
                    .and.                                                       &
                    (node%I(ii,ni_routeTo) .le. setting%SWMMinput%N_subcatch) &
                    ) then

                    !% --- set an alias for this subcatchment being routed to
                    thisSub => node%I(ii,ni_routeTo)
                    !% --- set the logical control
                    subcatchYN(thisSub,sYN_hasRunOn) = .true.
                    !% --- increment counter for number of runon connections to this
                    !%     subcatchment by 1
                    subcatchI (thisSub,si_RunOn_count) = subcatchI (thisSub,si_RunOn_count) + 1
                    !% --- alias for the column in the si_RunOn_...() arrays.
                    thisRunon => subcatchI (thisSub,si_RunOn_count)
                    !% --- store the node and face indexes for this runon
                    subcatchI(thisSub,si_RunOn_nodeIdx(thisRunon)) = ii
                    subcatchI(thisSub,si_RunOn_faceIdx(thisRunon)) = node%I(ii,ni_face_idx)
                    !% --- store the EPA SWMM outfall index used to export runon back to EPA SWMM
                    subcatchI(thisSub,si_RunOn_SWMMoutfallIdx(thisRunOn)) = node%I(ii,ni_SWMMoutfallIdx)
                    !% --- store the image (partition) associated with the node
                    subcatchI(thisSub,si_RunOn_P_image(thisRunon)) = node%I(ii,ni_P_image)
                    !% --- initialize the volume storage accumulators
                    subcatchR(thisSub,sr_RunOn_Volume(:)) = zeroR
                else 
                    !% no action
                end if

            end do
        else
            !% --- no action, initialization already done.
        end if
    
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
        !% Declarations
            real(8)  :: ttime
        !%------------------------------------------------------------------

        !% HACK start time is always measured from zero (need to fix for hotstart)
        setting%Time%Start = zeroR 
        setting%Time%Now   = zeroR
        setting%Time%Step  = zeroR

        setting%Time%Hydraulics%Dt = setting%VariableDT%InitialDt

        if (setting%Time%useSWMMinpYN) then 
            !% set the start/stop times and time steps from SWMM *.inp file
            setting%Time%StartEpoch    = setting%SWMMinput%StartEpoch
            setting%Time%EndEpoch      = setting%SWMMinput%EndEpoch
            !setting%Time%Hydraulics%Dt = setting%SWMMinput%RouteStep -- do not use!
            setting%Time%Hydrology%Dt  = setting%SWMMinput%Hydrology_WetStep
            ! HACK (not implemented)                  = setting%SWMMinput%DryStep
            ! HACK (not implemented)                  = setting%SWMMinput%TotalDuration
        else 
            !% --- use values from json file
        end if

        !% --- set the next time for updating climate if hydrology
        !%     is not used
        if ((.not. setting%Simulation%useHydrology)            &
            .and.                                              &
            (setting%Climate%useHydraulicsEvaporationTF)) then

            setting%Climate%LastTimeUpdate = setting%Time%Now

            setting%Climate%NextTimeUpdate = setting%Climate%LastTimeUpdate  &
                + setting%Climate%HydraulicsOnlyIntervalHours * seconds_per_hour
        end if

        !% --- SWMM uses epoch, in days, which has a precision of ~1e-4 seconds
        !%     To prevent precision from having an influence in microseconds, we
        !%     take the epoch difference, convert to seconds, then multiply by 10^4
        !%     and then round to the nearest integer.  We then divide by 10^4 to
        !%     get the number of seconds.  A final application of floor() and conversion
        !%     back to real ensures we only have whole seconds

        ttime = (setting%Time%EndEpoch - setting%Time%StartEpoch) * real(secsperday,KIND=8)
        ttime = util_datetime_seconds_precision (ttime)
        setting%Time%End = ttime

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
        !%------------------------------------------------------------------
        !% Description:
        !% initializes the output report time interval
        !%------------------------------------------------------------------

        !% --- if setting requires the SWMM input file values, then overwrite setting values
        if (setting%Output%Report%useSWMMinpYN) then 
            setting%Output%Report%StartTime    = util_datetime_epoch_to_secs(setting%SWMMinput%ReportStartTimeEpoch)
            setting%Output%Report%TimeInterval = setting%SWMMinput%ReportTimeInterval
            if ((setting%Output%Verbose) .and. (this_image() == 1)) then
                write(*,"(A)") ' ... using report start time and time interval from SWMM input file (*.inp)'
            end if
        else 
            if ((setting%Output%Verbose) .and. (this_image() == 1)) then
                write(*,"(A)") '... using  report start time and time interval from *.json file'
            end if
        end if

        !% --- if selected report time is before the start time use the start time
        if (setting%Output%Report%StartTime < setting%Time%Start) then 
            setting%Output%Report%StartTime = setting%Time%Start
        else 
            !% continue
        end if

        if (setting%Output%Report%TimeInterval < zeroR) then 
            if (this_image() == 1) then
                write(*,*) '***************************************************************'
                write(*,*) '** WARNING -- selected report time interval is zero or less, **'
                write(*,*) '**          so all output will be suppressed              **'
                write(*,*) '***************************************************************'
            end if
            setting%Output%Report%provideYN = .false.
            setting%Output%Report%suppress_MultiLevel_Output = .true.
            setting%Output%Report%ThisStep = 1
        else 
            !% --- Initialize the first report step
            !%     Determine how many report steps have already been missed before
            !%     the output reports are actually written
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
        !%------------------------------------------------------------------
        !% Description:
        !% Determines the overall length of the common coarray to handle the different
        !% number of elements on each processor
        !%-----------------------------------------------------------------
        !% Declarations
            integer :: nimgs_assign
            integer, allocatable :: unique_imagenum(:)
            integer :: ii, jj, idx, counter, elem_counter=0, face_counter=0, junction_counter=0

            integer :: duplicated_face_counter=0
            integer, allocatable :: node_index(:), link_index(:)
            character(64) :: subroutine_name = 'init_coarray_length'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%utility_array) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------

        call util_image_number_calculation(nimgs_assign, unique_imagenum)

        call util_allocate_scalar_for_images ()

        do ii=1, num_images()

            !% --- create corresponding indices for node and link in this image
            node_index = PACK([(counter, counter=1,size(node%I,1))], node%I(:, ni_P_image) == unique_imagenum(ii))
            link_index = PACK([(counter, counter=1,size(link%I,1))], link%I(:, li_P_image) == unique_imagenum(ii))
            
            !% --- The number of elements and faces is affected by the number of nJm junctions
            !%     So we will calculate the number of junction and use this in the setup
            junction_counter = count(node%I(node_index, ni_node_type) == nJm)

            !% --- first calculate the number of nodes in each partition, assign elems/faces for junctions
            !%     J_elem_add is the total number of elements reserved for each junction
            !%     J_face_add is the total number of faces reserved for each junction
            elem_counter = elem_counter + J_elem_add * junction_counter
            face_counter = face_counter + J_face_add * junction_counter

            !% --- loop through the links and count the internal faces between elements
            do jj = 1, size(link_index,1)
                idx = link_index(jj)
                face_counter = face_counter + link%I(idx, li_N_element)-1 !% internal faces between elems, e.g. 5 elements have 4 internal faces
                elem_counter = elem_counter + link%I(idx, li_N_element)   !% number of elements
            end do

            !% --- loop through the nodes and count the node faces
            do jj = 1, size(node_index,1)
                idx = node_index(jj)
                if (node%I(idx, ni_node_type) == nJ2) then
                    face_counter = face_counter + 1 !% add the face of 1-to-1 junction between 2 links
                elseif ((node%I(idx, ni_node_type) == nBCup) .or. (node%I(idx, ni_node_type) == nJ1)) then  
                    face_counter = face_counter +1 !% add the upstream faces
                elseif (node%I(idx, ni_node_type) == nBCdn) then
                    face_counter = face_counter +1 !% add the downstream faces
                elseif (node%I(idx, ni_node_type) == nJm) then   
                    !% skip -- faces are counted elsewhere
                else 
                    print *, jj, node%I(idx, ni_node_type), reverseKey(node%I(idx, ni_node_type))
                    print *, reverseKey(nJ1), reverseKey(nJ2), reverseKey(nBCup), reverseKey(nBCdn)
                    print *, 'CODE ERROR, unexpected else'
                    call util_crashpoint(390715)
                end if !% multiple junction faces already counted
            end do

            !% --- count the space for duplicated faces
            do jj = 1, size(link_index,1)
                idx = link_index(jj)
                !% --- check upstream node first
                if ( ( node%I(link%I(idx, li_Mnode_u), ni_P_is_boundary) == 1)  &
                     .and.                                                      &
                     ( node%I(link%I(idx, li_Mnode_u), ni_P_image) .ne. ii)     &
                    ) then
                    face_counter = face_counter +1
                    duplicated_face_counter = duplicated_face_counter + 1
                end if
                !% --- check  downstream node
                if ( ( node%I(link%I(idx, li_Mnode_d), ni_P_is_boundary) == 1)  & 
                     .and.                                                      &
                     ( node%I(link%I(idx, li_Mnode_d), ni_P_image) .ne. ii)     &
                    ) then
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

        !% --- assign the max size across all the coarrays
        max_caf_elem_N = maxval(N_elem)
        max_caf_face_N = maxval(N_face) 

        if (setting%Debug%File%utility_array) then
            do ii = 1, size(unique_imagenum,1)
                print *, 'Processor =>      ', ii
                print *, 'Elements expected ', N_elem(ii)
                print *, 'Faces expected    ', N_face(ii)
            end do
        end if

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%utility_array)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_coarray_length
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine init_linkarray_broadcast()
        !%------------------------------------------------------------------
        !% Description
        !% Broadcasts the link data unique to a particular link image
        !% to the other images (e.g., the local indexes of elements
        !% in the link). This is clunky because you cannot broadcast
        !% using a packed array of indices -- you have to do the entire array.
        !% This approach works because every image knows which links belong
        !% to which images through link%I(:,li_P_image)
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: ii, kk, cset(3)
            integer, allocatable :: linkThisImage(:), tempI(:)[:]
        !%------------------------------------------------------------------

        sync all
        do ii=1,num_images()
            !% --- get all the links for each image
            !%     Note that this_image() conducts this for every image, but
            !%     only when this_image()==ii will the values be the ones 
            !%     that are broadcast
            linkThisImage = pack(link%I(:,li_idx), link%I(:,li_P_image) == ii)

            !% --- check that this image has links, if so allocate 
            !%     the temporary integer space as a coarray
            if (size(linkThisImage) > 0) then
                !% --- allocate every loop because we have different # of links on each image
                allocate(tempI(size(linkThisImage))[*])
                tempI(:) = zeroI
            else
                cycle
            end if

            !% --- columns to be broadcast
            cset(1) = li_N_element
            cset(2) = li_first_elem_idx
            cset(3) = li_last_elem_idx

            !% --- cycle through the columns to broadcast
            do kk=1,3  !% increase this if more data needs to be broadcast in cset
                !% --- store the link data in the single array
                !%     note that if this data is all nullvalueI unless ii==this_image()
                tempI = link%I(linkThisImage,cset(kk))[ii]
                !% --- broadcast from the source image to all the others
                !%     which overwrites the nullvalueI
                call co_broadcast (tempI, source_image=ii)
                !%
                if (ii .ne. this_image()) then
                    !% store the broadcast data back in the link array
                    link%I(linkThisImage,cset(kk)) = tempI
                end if
            end do
            
            !% --- note we have to deallocate and allocate in every
            !%     loop because the different images have different numbers of elements
            deallocate(tempI)
        end do    

        !%------------------------------------------------------------------
        !% Closing
            deallocate(linkThisImage)
        
    end subroutine init_linkarray_broadcast
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine init_boundary_ghost_elem_array()
        !%------------------------------------------------------------------
        !% Description:
        !% initialize ghost and boundary element arrays for inter image data transfer
        !%-------------------------------------------------------------------
        !% Declarations
            integer          :: ii, jj, NSfaces, eset_local(4), eBGset(4)
            integer, pointer :: Nfaces, fIdx, fGidx, eUp, eDn, ci, BeUp, BeDn
            integer, dimension(:), allocatable, target :: packed_shared_face_idx
            character(64)    :: subroutine_name = 'init_boundary_ghost_elem_array'
        !-------------------------------------------------------------------
        !% Preliminaries
            !% only initialize inter-image data transfer array for more than one processor
            if (num_images() .le. 1) return

            if (setting%Debug%File%network_define) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !-------------------------------------------------------------------
        
        !% --- allocate elemB and elemG data structure
        call util_allocate_boundary_ghost_elem_array()

        !% --- number of faces in this image
        Nfaces => N_face((this_image()))
        !% --- count the number of shared faces in this image
        NSfaces = count(faceYN(1:Nfaces,fYN_isSharedFace))
        !% --- packed indexes of shared faces in this image
        packed_shared_face_idx = pack(faceI(1:Nfaces,fi_Lidx),faceYN(1:Nfaces,fYN_isSharedFace))
        !% --- column indexes of elemI local data needed to be transferred 
        eset_local = [ei_Lidx, ei_Gidx, ei_Mface_uL, ei_Mface_dL]
        !% --- column indexes of the spot for receiving local elemI data
        eBGset     = [ebgi_elem_Lidx, ebgi_elem_Gidx, ebgi_Mface_uL, ebgi_Mface_dL]

        !% HACK 
        !% The two do-loops below provide diagnostic for shared faces in the elemB array. 
        !% This could be  simplified, by replacing BeUp and BeDn pointers in the second
        !% do-loope with ii. However, the current state of mapping is more self explanatory
        !% and easier to debug

        !% --- cycle through faces to identify elements for the elemB set and store their
        !%     corresponding boundary index in te elemI and faceI arrays 
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
                print *, 'CODE ERROR unexpected else' 
                print *, 'Unexpected ELSE IF condition for single element with two shared faces'
                call util_crashpoint(914744)
            else
                print *, 'CODE ERROR unexpected else' 
                print *, 'Unexpected ELSE condition for shared faces'
                call util_crashpoint(54673)
            end if
        end do

        !% --- wait and sync until all the images update their elemB arrays with local data
        sync all

        !% --- loop through all the shared faces again to find the boundary array location 
        !%     of the ghost element in a remote image
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

                    !% --- find the local ghost element index of the connected image [ci]
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

        !-------------------------------------------------------------------
        !% Closing
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
            logical :: ifail  = .false.
        !%------------------------------------------------------------------  
        !% Preliminaries: 
            if (this_image() .ne. 1) return
        !%------------------------------------------------------------------ 
        !write(*,*) ' '
        !write(*,'(A)') '*******************************************************************'
 
        
        if (.not. setting%ZeroValue%UseZeroValues) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.ZeroValue.UseZeroValues = false, which is not tested and not presently allowed'
            write(*,'(A)') '** '
            ifound = .true.
            ifail  = .true.
        end if
        if (.not. setting%VariableDT%ApplyYN) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.VariableDT%ApplyYN = false'
            write(*,'(A)') '** This option is not tested and not presently allowed.'
            write(*,'(A)') '** '
            ifound = .true.
            ifail  = .true.
        end if
        if (.not. setting%Time%useSWMMinpYN) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.Time.useSWMMinYN = false'
            write(*,'(A)') '** The simulation will use the timing data from the json file '
            write(*,'(A)') '** rather than fromthe SWMM *.inp file.' 
            write(*,'(A)') '** Use this option with caution.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (setting%TestCase%UseTestCaseYN) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.TestCase.UseTestCaseYN = true'
            write(*,'(A)') '** This option is not tested and not presently allowed.'
            write(*,'(A)') '** '
            ifound = .true.
            ifail  = .true.
        end if
        if (.not. setting%Solver%PreissmannSlot%useSlotTF) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.Solver.PreissmannSlot.useSlotTF = false'
            write(*,'(A)') '** This option is not tested and not presently allowed.'
            write(*,'(A)') '** '
            ifound = .true.
            ifail  = .true.
        end if
        if (setting%Solver%SubtractReferenceHead) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.Solver.SubtractReferenceHead = true'
            write(*,'(A)') '** This option is not tested and not presently allowed.'
            write(*,'(A)') '** '
            ifound = .true.
            ifail  = .true.
        end if
        if (.not. setting%Simulation%useHydraulics) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.Simulation.useHydraulics = false'
            write(*,'(A)') '** This option is not tested and not presently allowed.'
            write(*,'(A)') '** '
            ifound = .true.
            ifail  = .true.
        end if
        if (setting%Profile%useYN) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.Profile.useYN = true. '
            write(*,'(A)') '** This code timing profiler is still in development, so results should '
            write(*,'(A)') '** be used with caution.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%Output%Verbose) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.Output.Verbose = false. '
            write(*,'(A)') '** While SWMM5+ is in development we recommend using true. '
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%Output%Warning) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.Output.Warning = false. '
            write(*,'(A)') '** While SWMM5+ is in development we recommend using true. '
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%File%force_folder_creationYN) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.File.force_folder_creationYN = false.'
            write(*,'(A)') '** This option is not tested and not presently allowed.'
            write(*,'(A)') '** '
            ifound = .true.
            ifail  = .true.
        end if
        if (.not. setting%File%UseCommandLineFoldersYN) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.File.UseCommandLineFoldersYN = false. '
            write(*,'(A)') '** This option is not tested and not presently allowed.'
            write(*,'(A)') '** '
            ifound = .true.
            ifail  = .true.
        end if
        if (setting%BC%disableInterpolationYN) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.BC.disableInterpolationYN = true.'
            write(*,'(A)') '** While SWMM5+ is in development we recommend using false. '
            write(*,'(A)') '** This disables interpolation of time series and simply takes the'
            write(*,'(A)') '** next highest value.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%Output%Report%useSWMMinpYN) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.Output.Report.useSWMMinYN = false.'
            write(*,'(A)') '** The simulation will use the report data from the json file '
            write(*,'(A)') '** rather than from the SWMM *.inp file.'
            write(*,'(A)') '** Use this option with caution.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%Output%Report%provideYN) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.Output.Report.provideYN = false.'
            write(*,'(A)') '** The simulation will suppress the output reporting.'
            write(*,'(A)') '** Use this option with caution.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (setting%Output%Report%suppress_MultiLevel_Output) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.Report.suppress_MultiLevel_Output = true.'
            write(*,'(A)') '** The simulation will suppress all the output reporting from hydraulics.'
            write(*,'(A)') '** Presently, the only output from hydraulics is through the ML scheme.'
            write(*,'(A)') '** Use this option with caution.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%Limiter%Velocity%UseLimitMaxYN) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.Limiter.Velocity.UseLimitMaxYN = false.'
            write(*,'(A)') '** We recommend using true'
            write(*,'(A)') '** Without a velocity limiter, a minor instability can become'
            write(*,'(A)') '** catastrophic.'
            write(*,'(A)') '** '
            ifound = .true.
        end if
        if (.not. setting%Limiter%Dt%UseLimitMinYN) then
            if (.not. ifound) then 
                write(*,*) ' '
                write(*,'(A)') '*******************************************************************'
            end if
            write(*,'(A)') '** setting.Limiter.Dt.UseLimitMinYN = false.  '
            write(*,'(A)') '** We recommend using true.'
            write(*,'(A)') '** The variable dt will become very small if the simulation'
            write(*,'(A)') '** is going unstable, resulting in the code seeming to hang.'
            write(*,'(A)') '** Using the minimum value ensures the code stops if it is'
            write(*,'(A)') '** setting a small time step due to instability.'
            write(*,'(A)') '** '
            ifound = .true.
        end if

        if (ifound) then
            write(*,'(A)') '**                          WARNING'
            write(*,'(A)') '** The above are possible issues with user settings for this version of SWMM5+.  '
            write(*,'(A)') '** '
            write(*,'(A)') '*******************************************************************'
        end if

        if (ifail) then 
            write(*,'(A)') '** '
            write(*,'(A)') '** USER CONFIGURATION ERROR FOUND'
            write(*,'(A)') '*******************************************************************'
            call util_crashpoint(798723)
        end if

        !%------------------------------------------------------------------
    end subroutine init_check_setup_conditions    
!% 
!%==========================================================================
!%==========================================================================
!%
    subroutine init_viscosity ()
        !%------------------------------------------------------------------
        !% Description:
        !% Initializes the kinematic viscosity of water
        !% Note that temperature limits are hard-coded based on
        !% data used in util_kinematic_viscosity_from_temperature
        !%------------------------------------------------------------------
        !% Declarations
            real(8), pointer :: temperature, viscosity
        !%------------------------------------------------------------------
        !% Aliases
            temperature => setting%Constant%water_temperature
            viscosity   => setting%Constant%water_kinematic_viscosity
        !%------------------------------------------------------------------
        if (      (temperature .ge. -8.d0) &
            .and. (temperature .le. 70.d0)   ) then
                viscosity = util_kinematic_viscosity_from_temperature(temperature)
        else
            print *, 'USER CONFIGURATION ERROR water temperature'
            print *, 'Water temperature outside valid range of -8C < T < 70 C,'
            print *, 'which is the valid range for computation of viscosity'
            print *, 'Value in setting.Constant.water_temperature is ',temperature
            call util_crashpoint(408734)
        end if

    end subroutine init_viscosity
!% 
!%==========================================================================
!%==========================================================================
!%
    subroutine init_ForceMain_setting ()
        !%------------------------------------------------------------------
        !% Description:
        !% Sets the AllowForceMainTF to false if no FM are found in the
        !% SWMMinput file links. 
        !%------------------------------------------------------------------
        !%
        !% --- check to see if any Force Main links exists
        if (setting%Solver%ForceMain%AllowForceMainTF) then
            if (any(link%I(:,li_geometry) .eq. lForce_main)) then
                !% --- maintain "true" for AllowForceMainTF
            else 
                !% --- no force main found
                setting%Solver%ForceMain%AllowForceMainTF = .false.
            end if
        end if

        !% --- check to see if user requests all closed conduits to be
        !%     treated as force mains. 
        if (setting%Solver%ForceMain%FMallClosedConduitsTF) then 
            if (.not. setting%Solver%ForceMain%AllowForceMainTF) then
                !% --- inconsistent between using force main and
                !%     setting all to closed conduits
                if (this_image() == 1) then
                    print *, 'USER CONFIGURATION ERROR Inconsistency in settings.'
                    print *, 'setting.Solver.ForceMain.AllowForceMainTF is false '
                    print *, 'setting.Solver.ForceMain.FMallClosedConduitsTF is true.'
                    print *, 'If you want to use Force main with all closed conduits then set '
                    print *, 'setting.Solver.ForceMain.AllowForceMainTF = true '
                    print *, 'otherwise set setting.Solver.ForceMain.FMallClosedConduits = false '
                end if
                call util_crashpoint(4298732)
            end if
        else 
            !% --- FMallClosedConduits = false is compatible with any value of
            !%     AllowForceMainTF
        end if

    end subroutine init_ForceMain_setting
!% 
!%==========================================================================
!%==========================================================================
!%   
    subroutine init_culvert ()
        !%------------------------------------------------------------------
        !% Description:
        !% Inialize culvert parameters and numbers.
        !%------------------------------------------------------------------

        !% --- set up the culvert parameters array
        call culvert_parameter_values ()

    end subroutine init_culvert
!%  
!%==========================================================================
!% END OF MODULE
!%==========================================================================
!%
end module initialization

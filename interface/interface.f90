module interface

    use iso_c_binding
    use c_library
    use utility
    use utility_datetime
    use define_indexes
    use define_keys
    use define_api_keys
    use define_globals
    use define_settings, only: setting

    implicit none

    private

!%==========================================================================
!% PUBLIC
!%==========================================================================

    type(c_lib_type), public :: c_lib
    logical, public :: api_is_initialized = .false.

    !% Public subroutines/functions
    public :: interface_init
    public :: interface_finalize
    ! public :: interface_run_step
    public :: interface_get_nodef_attribute
    public :: interface_get_linkf_attribute
    public :: interface_get_table_attribute
    public :: interface_get_num_table_entries
    public :: interface_get_first_entry_table
    public :: interface_get_next_entry_table
    public :: interface_get_obj_name_len
    public :: interface_update_linknode_names
    public :: interface_get_BC_resolution
    public :: interface_get_next_inflow_time
    public :: interface_get_next_head_time
    public :: interface_get_flowBC
    public :: interface_get_headBC
    public :: interface_find_object
    public :: interface_update_nodeResult
    public :: interface_update_linkResult
    public :: interface_write_output_line
    public :: interface_export_link_results
    public :: interface_call_runoff_execute
    public :: interface_get_subcatch_runoff
    public :: interface_get_subcatch_runoff_nodeIdx
    public :: interface_get_NewRunoffTime

!%==========================================================================
!% PRIVATE
!%==========================================================================

    ! Interface with SWMM shared library
    interface
        !% -------------------------------------------------------------------------------
        !% Simulation
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_initialize(inp_file, report_file, out_file, run_routing) &
            BIND(C, name="api_initialize")
            use, intrinsic :: iso_c_binding
            implicit none
            character(kind=c_char), intent(in) :: inp_file
            character(kind=c_char), intent(in) :: report_file
            character(kind=c_char), intent(in) :: out_file
            integer(c_int),  value, intent(in) :: run_routing
        end function api_initialize
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_finalize() &
            BIND(C, name="api_finalize")
            use, intrinsic :: iso_c_binding
            implicit none
        end function api_finalize
        !% -------------------------------------------------------------------------------
        real(c_double) function api_run_step() &
            BIND(C, name="api_run_step")
            use, intrinsic :: iso_c_binding
            implicit none
        end function api_run_step
        !% -------------------------------------------------------------------------------
        ! --- Property-extraction
        !% -------------------------------------------------------------------------------
        !% * During Simulation
        integer(c_double) function api_get_node_results(node_name, inflow, overflow, depth, volume) &
            BIND(C, name="api_get_node_results")
            use, intrinsic :: iso_c_binding
            implicit none
            character(kind=c_char), intent(in   ) :: node_name
            real(c_float),          intent(inout) :: inflow
            real(c_float),          intent(inout) :: overflow
            real(c_float),          intent(inout) :: depth
            real(c_float),          intent(inout) :: volume
        end function api_get_node_results
        !% -------------------------------------------------------------------------------
        integer(c_double) function api_get_link_results(link_name, flow, depth, volume) &
            BIND(C, name="api_get_link_results")
            use, intrinsic :: iso_c_binding
            implicit none
            character(kind=c_char), intent(in   ) :: link_name
            real(c_float),          intent(inout) :: flow
            real(c_float),          intent(inout) :: depth
            real(c_float),          intent(inout) :: volume
        end function api_get_link_results
        !% -------------------------------------------------------------------------------
        !% * After Initialization
        real(c_double) function api_get_start_datetime() &
            BIND(C, name="api_get_start_datetime")
            use, intrinsic :: iso_c_binding
            implicit none
        end function api_get_start_datetime
        !% -------------------------------------------------------------------------------
        real(c_double) function api_get_end_datetime() &
            BIND(C, name="api_get_end_datetime")
            use, intrinsic :: iso_c_binding
            implicit none
        end function api_get_end_datetime
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_get_flowBC(node_idx, current_datetime, flowBC) &
            BIND(C, name="api_get_flowBC")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: node_idx
            real(c_double), value, intent(in   ) :: current_datetime
            real(c_double),        intent(inout) :: flowBC
        end function api_get_flowBC
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_get_headBC(node_idx, current_datetime, headBC) &
            BIND(C, name="api_get_headBC")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: node_idx
            real(c_double), value, intent(in   ) :: current_datetime
            real(c_double),        intent(inout) :: headBC
        end function api_get_headBC
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_get_SWMM_times &
            (starttime_epoch, endtime_epoch, report_start_datetime, report_step, &
             hydrology_step, hydrology_dry_step, hydraulic_step, total_duration) &
            BIND(C, name="api_get_SWMM_times")
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double), intent(inout) :: starttime_epoch
            real(c_double), intent(inout) :: endtime_epoch
            real(c_double), intent(inout) :: report_start_datetime
            integer(c_int), intent(inout) :: report_step
            integer(c_int), intent(inout) :: hydrology_step
            integer(c_int), intent(inout) :: hydrology_dry_step
            real(c_double), intent(inout) :: hydraulic_step
            real(c_double), intent(inout) :: total_duration
        end function api_get_SWMM_times
        !% -------------------------------------------------------------------------------
        real(c_double) function api_get_NewRunoffTime() &
            BIND(C, name="api_get_NewRunoffTime")
            use, intrinsic :: iso_c_binding
            implicit none 
        end function api_get_NewRunoffTime
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_get_nodef_attribute(node_idx, attr, value) &
            BIND(C, name="api_get_nodef_attribute")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: node_idx
            integer(c_int), value, intent(in   ) :: attr
            real(c_double),        intent(inout) :: value
        end function api_get_nodef_attribute
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_get_linkf_attribute(link_idx, attr, value) &
            BIND(C, name="api_get_linkf_attribute")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: link_idx
            integer(c_int), value, intent(in   ) :: attr
            real(c_double),        intent(inout) :: value
        end function api_get_linkf_attribute
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_get_num_objects(obj_type) &
            BIND(C, name="api_get_num_objects")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: obj_type
        end function api_get_num_objects
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_get_object_name_len(object_idx, object_type, len_value) &
            BIND(C, name="api_get_object_name_len")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: object_idx
            integer(c_int), value, intent(in   ) :: object_type
            integer(c_int),        intent(inout) :: len_value
        end function api_get_object_name_len
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_get_object_name(object_idx, object_name, object_type) &
            BIND(C, name="api_get_object_name")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value,  intent(in   ) :: object_idx
            character(kind=c_char), intent(inout) :: object_name
            integer(c_int), value,  intent(in   ) :: object_type
        end function api_get_object_name
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_get_num_table_entries(table_idx, table_type, num_entries) &
            BIND(C, name="api_get_num_table_entries")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: table_idx
            integer(c_int), value, intent(in   ) :: table_type
            integer(c_int),        intent(inout) :: num_entries
        end function api_get_num_table_entries
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_get_table_attribute(table_idx, attr, value) &
            BIND(C, name="api_get_table_attribute")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: table_idx
            integer(c_int), value, intent(in   ) :: attr
            real(c_double),        intent(inout) :: value
        end function api_get_table_attribute
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_get_first_entry_table(table_idx, table_type, x, y) &
            BIND(C, name="api_get_first_entry_table")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: table_idx
            integer(c_int), value, intent(in   ) :: table_type
            real(c_double),        intent(inout) :: x
            real(c_double),        intent(inout) :: y
        end function api_get_first_entry_table
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_get_next_entry_table(table_idx, table_type, x, y) &
            BIND(C, name="api_get_next_entry_table")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: table_idx
            integer(c_int), value, intent(in   ) :: table_type
            real(c_double),        intent(inout) :: x
            real(c_double),        intent(inout) :: y
        end function api_get_next_entry_table
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_get_next_entry_tseries(tseries_idx) &
            BIND(C, name="api_get_next_entry_tseries")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in   ) :: tseries_idx
        end function api_get_next_entry_tseries
        !% -------------------------------------------------------------------------------
        ! --- Output Writing (Post Processing)
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_write_output_line(t) &
            BIND(C, name="api_write_output_line")
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double), value, intent(in) :: t
        end function api_write_output_line
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_update_nodeResult(node_idx, resultType, newNodeResult) &
            BIND(C, name="api_update_nodeResult")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: node_idx
            integer(c_int), value, intent(in) :: resultType
            real(c_double), value, intent(in) :: newNodeResult
        end function api_update_nodeResult
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_update_linkResult(link_idx, resultType, newLinkResult) &
            BIND(C, name="api_update_linkResult")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: link_idx
            integer(c_int), value, intent(in) :: resultType
            real(c_double), value, intent(in) :: newLinkResult
        end function api_update_linkResult
        !% -------------------------------------------------------------------------------
        ! --- Print-out
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_export_linknode_properties(units) &
            BIND(C, name="api_export_linknode_properties")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: units
        end function api_export_linknode_properties
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_export_link_results(link_idx) &
            BIND(C, name="api_export_link_results")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: link_idx
        end function api_export_link_results
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_export_node_results(node_idx) &
            BIND(C, name="api_export_node_results")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: node_idx
        end function api_export_node_results
        !% -------------------------------------------------------------------------------
        ! --- Hydrology
        !% -------------------------------------------------------------------------------   
        integer(c_int) function api_call_runoff_execute() &
            BIND(C, name='api_call_runoff_execte')
            use, intrinsic :: iso_c_binding
            implicit none 
        end function api_call_runoff_execute
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_get_subcatch_runoff(sc_idx,runoff) &
            BIND(C, name="api_get_subcatch_runoff")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in)    :: sc_idx
            real(c_double),        intent(inout) :: runoff
        end function api_get_subcatch_runoff
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_get_subcatch_runoff_nodeIdx(sc_idx, node_idx) &
            BIND(C, name="api_get_subcatch_runoff_nodeIdx")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in)    :: sc_idx
            integer(c_int),        intent(inout) :: node_idx
        end function api_get_subcatch_runoff_nodeIdx
        !% -------------------------------------------------------------------------------
        ! --- Utils
        !% -------------------------------------------------------------------------------
        integer(c_int) function api_find_object(object_type, object_name) &
            BIND(C, name="api_find_object")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value,  intent(in) :: object_type
            character(kind=c_char), intent(in) :: object_name
        end function api_find_object

    end interface
!%
!%==========================================================================
!% Procedures
!%==========================================================================
!%

    procedure(api_initialize),                 pointer :: ptr_api_initialize
    procedure(api_finalize),                   pointer :: ptr_api_finalize
    procedure(api_run_step),                   pointer :: ptr_api_run_step
    procedure(api_get_node_results),           pointer :: ptr_api_get_node_results
    procedure(api_get_link_results),           pointer :: ptr_api_get_link_results
    procedure(api_get_start_datetime),         pointer :: ptr_api_get_start_datetime
    procedure(api_get_end_datetime),           pointer :: ptr_api_get_end_datetime
    procedure(api_get_flowBC),                 pointer :: ptr_api_get_flowBC
    procedure(api_get_headBC),                 pointer :: ptr_api_get_headBC
    procedure(api_get_SWMM_times),             pointer :: ptr_api_get_SWMM_times
    procedure(api_get_NewRunoffTime),          pointer :: ptr_api_get_NewRunoffTime
    procedure(api_get_nodef_attribute),        pointer :: ptr_api_get_nodef_attribute
    procedure(api_get_linkf_attribute),        pointer :: ptr_api_get_linkf_attribute
    procedure(api_get_num_objects),            pointer :: ptr_api_get_num_objects
    procedure(api_get_object_name_len),        pointer :: ptr_api_get_object_name_len
    procedure(api_get_object_name),            pointer :: ptr_api_get_object_name
    procedure(api_get_num_table_entries),      pointer :: ptr_api_get_num_table_entries
    procedure(api_get_table_attribute),        pointer :: ptr_api_get_table_attribute
    procedure(api_get_first_entry_table),      pointer :: ptr_api_get_first_entry_table
    procedure(api_get_next_entry_table),       pointer :: ptr_api_get_next_entry_table
    procedure(api_get_next_entry_tseries),     pointer :: ptr_api_get_next_entry_tseries
    procedure(api_write_output_line),          pointer :: ptr_api_write_output_line
    procedure(api_update_nodeResult),          pointer :: ptr_api_update_nodeResult
    procedure(api_update_linkResult),          pointer :: ptr_api_update_linkResult
    procedure(api_export_linknode_properties), pointer :: ptr_api_export_linknode_properties
    procedure(api_export_link_results),        pointer :: ptr_api_export_link_results
    procedure(api_export_node_results),        pointer :: ptr_api_export_node_results
    procedure(api_find_object),                pointer :: ptr_api_find_object
    procedure(api_call_runoff_execute),        pointer :: ptr_api_call_runoff_execute
    procedure(api_get_subcatch_runoff),        pointer :: ptr_api_get_subcatch_runoff 
    procedure(api_get_subcatch_runoff_nodeIdx), pointer :: ptr_api_get_subcatch_runoff_nodeIdx 
    
    !% Error handling
    character(len = 1024) :: errmsg
    integer :: errstat

contains

!%=============================================================================
!% PUBLIC
!%=============================================================================

!%=============================================================================
!%  Simulation subroutines/functions
!%=============================================================================

    subroutine interface_init()
        !%---------------------------------------------------------------------
        !% Description:
        !%    initializes the EPA-SWMM shared library, creating input (.inp),
        !%    report (.rpt), and output (.out) files, necessary to run simulation with
        !%    EPA-SWMM. It also updates the number of objects in the SWMM model, i.e.,
        !%    number of links, nodes, and tables, and defines the start and end
        !%    simulation times.
        !%----------------------------------------------------------------------
            integer :: ppos, num_args, error
            character(64) :: subroutine_name = 'interface_init'
        
        !% Preliminaries:
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        setting%File%inp_file = trim(setting%File%inp_file) // c_null_char
        setting%File%rpt_file = trim(setting%File%rpt_file) // c_null_char
        setting%File%out_file = trim(setting%File%out_file) // c_null_char
        c_lib%filename = trim(setting%File%library_folder) // "/libswmm5.so"

        !% --- initialize the api between SWMM-C and SWMM5+
        !%     This returns and stores the SWMM-C input and output filenames
        call load_api_procedure("api_initialize")
        error = ptr_api_initialize( &
            setting%File%inp_file, &
            setting%File%rpt_file, &
            setting%File%out_file, &
            dummyI)
        call print_api_error(error, subroutine_name)
        api_is_initialized = .true.

        !% --- Get number of objects in SWMM-C
        SWMM_N_link = get_num_objects(API_LINK)
        N_link = SWMM_N_link

        SWMM_N_node = get_num_objects(API_NODE)
        N_node = SWMM_N_node

        SWMM_N_Curve = get_num_objects(API_CURVE)
        N_curve = SWMM_N_Curve

        SWMM_N_subcatch = get_num_objects(API_SUBCATCH)     

        SWMM_N_pollutant = get_num_objects(API_POLLUT)
        if (SWMM_N_pollutant > 0) then
            write(*,*) ' '
            write(*,*) '*************************************************************'
            write(*,*) '*                                                           *'
            write(*,*) '*  WARNING WARNING WARNING WARNING WARNING WARNING WARNING  *'
            write(*,*) '*                                                           *'
            write(*,*) '* The SWMM input file contains one or more pollutants.      *'
            write(*,*) '* The current version of SWMM5+ does not support pollutants *'
            write(*,*) '* The simulation will be run, but pollutants will not       *'
            write(*,*) '* be included in the hydraulic computations.                *'
            write(*,*) '*************************************************************'
        end if

               
        if ((N_link == 200) .AND. (N_node == 200)) then
            print *, ''
            print *, 'BUG WARNING location ',980879,' in ',subroutine_name
            print *, '...if the SWMM code detects a parse error for the *.inp file then the ...'
            print *, '...get_num_objects function returns an error code 200 (SWMM parse error)... '
            print *, '...that is stored instead of the names of nodes and links...'
            print *, '...this has unexpected errors.'
            print *, ''
            print *, SWMM_N_link, SWMM_N_node
            print *, 'ERROR (input file): Appears to be parse error in the input file...'
            print *, '...where some links/nodes are either not connected or not identified...'
            print *, '...This can happen if nodes are renamed and some of the conduit connection did not get changed...'
            print *, '...This error might have been tripped accidently if the system has exactly...'
            print *, '...200 nodes and 200 links.'
            stop 2398760
        end if

        !% --- get the time start, end, and interval data from SWMM-C input file
        call interface_get_SWMM_times()

        !%----------------------------------------------------------------------
        !% closing
            if (setting%Debug%File%interface) then
                print *, new_line("")
                print *, "SWMM_N_link", SWMM_N_link
                print *, "SWMM_N_node", SWMM_N_node
                print *, new_line("")
                print *, "SWMM start time", setting%Time%StartEpoch
                print *, "SWMM end time", setting%Time%EndEpoch
                print *, "setting%time%start", setting%Time%Start
                print *, "setting%time%end", setting%Time%End
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            end if
    end subroutine interface_init
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_finalize()
        !%---------------------------------------------------------------------
        !% Description:
        !% stops EPA SWMMC, deletes pointers and closese the shared library
        !%---------------------------------------------------------------------
        !% Declarations:
            integer :: errstat
            character(64) :: subroutine_name = 'interface_finalize'
        !%---------------------------------------------------------------------
        !% Preliminaries:
           if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%---------------------------------------------------------------------

        call load_api_procedure("api_finalize")
        errstat = ptr_api_finalize()
        call print_api_error(errstat, subroutine_name)

        !%---------------------------------------------------------------------
        !% Closing        
           if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine interface_finalize
!%
!%=============================================================================
!%=============================================================================
!%
    ! subroutine interface_run_step()
    ! !%-----------------------------------------------------------------------------
    ! !% Description:
    ! !%    runs steps of EPA-SWMM model. If setting%Simulation% was defined
    ! !%    true, steps include routing model. If false, steps are for hydrology only
    ! !%-----------------------------------------------------------------------------
    !     real(8), pointer :: timeNow
    !     character(64)    :: subroutine_name = 'interface_run_step'
    ! !%-----------------------------------------------------------------------------

    !     if (setting%Debug%File%interface)  &
    !     write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !     timeNow => setting%Time%Now

    !     c_lib%procname = "api_run_step"
    !     call c_lib_load(c_lib, errstat, errmsg)
    !     if (errstat /= 0) then
    !         print *, "ERROR: " // trim(errmsg)
    !         stop
    !     end if
    !     call c_f_procpointer(c_lib%procaddr, ptr_api_run_step)

    !     timeNow = ptr_api_run_step(api)
    !     if (setting%Debug%File%interface)  &
    !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    ! end subroutine interface_run_step
!%
!%=============================================================================
!%   Property-extraction functions (only run after initialization)
!%=============================================================================
!%
    subroutine interface_update_linknode_names()
        !%-----------------------------------------------------------------------------
        !% Description:
        !%    Updates the link%Names and node%Names arrays which should've been
        !%    allocated by the time the subroutine is executed. Names are copied
        !%    from EPA-SWMM
        !%-----------------------------------------------------------------------------
            integer :: ii
            character(64) :: subroutine_name = "interface_update_linknode_names"
        !%-----------------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        do ii = 1, SWMM_N_link
            call load_api_procedure("api_get_object_name")
            errstat = ptr_api_get_object_name(ii-1, link%Names(ii)%str, API_LINK)

            if (errstat /= 0) then
                write(*, "(A,i2,A)") "API ERROR : ", errstat, " [" // subroutine_name // "]"
                stop 8673489
            end if
        end do

        do ii = 1, SWMM_N_node
            call load_api_procedure("api_get_object_name")
            errstat = ptr_api_get_object_name(ii-1, node%Names(ii)%str, API_NODE)

            if (errstat /= 0) then
                write(*, "(A,i2,A)") "API ERROR : ", errstat, " [" // subroutine_name // "]"
                stop 48705
            end if
        end do

        if (setting%Debug%File%interface) then
            print *, new_line("")
            print *, "List of Links"
            do ii = 1, SWMM_N_link
                print *, "- ", link%Names(ii)%str
            end do
            print *, new_line("")
            print *, "List of Nodes"
            do ii = 1, SWMM_N_node
                print *, "- ", node%Names(ii)%str
            end do
            print *, new_line("")
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        end if

    end subroutine interface_update_linknode_names
!%
!%=============================================================================
!%=============================================================================
!%
    integer(c_int) function interface_get_obj_name_len(obj_idx, obj_type) result(len_value)
        !%-----------------------------------------------------------------------------
        !% Description:
        !%    Returns the length of the name string associated to the EPA-SWMM object.
        !%    This function is necessary to allocate the entries of the link%Names and
        !%    node%Names arraysr. The function is currently compatible with NODE and
        !%    LINK types.
        !%-----------------------------------------------------------------------------
            integer, intent(in) :: obj_idx  ! index of the EPA-SWMM object
            integer, intent(in) :: obj_type ! type of EPA-SWMM object (API_NODE, API_LINK)
            integer                :: error
            character(64)          :: subroutine_name = "interface_get_obj_name_len"
        !%-----------------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
        write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        call load_api_procedure("api_get_object_name_len")
        error = ptr_api_get_object_name_len(obj_idx-1, obj_type, len_value)
        call print_api_error(error, subroutine_name)
        
        if (setting%Debug%File%interface) then
            print *, obj_idx, obj_type, len_value
        end if

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end function interface_get_obj_name_len
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_nodef_attribute(node_idx, attr) result(node_value)
        !%-----------------------------------------------------------------------------
        !% Description:
        !%    Retrieves node attributes from EPA-SWMM. API node attributes are
        !%    defined in define_api_keys.f08.
        !% Notes:
        !%    Fortran indexes are translated to C indexes and viceversa when
        !%    necessary. Fortran indexes always start from 1, whereas C indexes
        !%    start from 0.
        !%-----------------------------------------------------------------------------
            integer :: node_idx, attr, error
            real(c_double), target :: node_value
            character(64) :: subroutine_name = 'interface_get_nodef_attribute'
        !%-----------------------------------------------------------------------------

        !write(*,*) '--- interface_get_nodef_attribute', attr

        if (setting%Debug%File%interface)  &
            write(*,"(3(A,i5),A)") '*** enter ' // trim(subroutine_name) // &
            "(node_idx=", node_idx, ", attr=", attr, ")" // " [Processor ", this_image(), "]"

        if ((attr > N_api_nodef_attributes) .or. (attr < 1)) then
            print *, "error: unexpected node attribute value", attr
            stop 948705
        end if

        if ((node_idx > N_node) .or. (node_idx < 1)) then
            print *, "error: unexpected node index value", node_idx
            stop 397904
        end if

        !% Substracts 1 to every Fortran index (it becomes a C index)
        call load_api_procedure("api_get_nodef_attribute")
        error = ptr_api_get_nodef_attribute(node_idx-1, attr, node_value)
        !print *, '   node value ',node_value
        call print_api_error(error, subroutine_name)

        !% Adds 1 to every C index extracted from EPA-SWMM (it becomes a Fortran index)
        if ((attr == api_nodef_extInflow_tSeries) .or. (attr == api_nodef_extInflow_basePat)) then
            if (node_value /= -1) node_value = node_value + 1
        end if

        !write(*,*) '.................'
        if (setting%Debug%File%interface) &
            write(*,"(3(A,i5),A)") '*** leave ' // trim(subroutine_name) // &
            "(node_idx=", node_idx, ", attr=", attr, ")" // " [Processor ", this_image(), "]"
    end function interface_get_nodef_attribute
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_linkf_attribute(link_idx, attr) result(link_value)
        !%-----------------------------------------------------------------------------
        !% Description:
        !%    Retrieves link attributes from EPA-SWMM. API link attributes are
        !%    defined in define_api_keys.f08.
        !% Notes:
        !%    Fortran indexes are translated to C indexes and viceversa when
        !%    necessary. Fortran indexes always start from 1, whereas C indexes
        !%    start from 0.
        !%-----------------------------------------------------------------------------
            integer :: link_idx, attr, error

            real(c_double), target :: link_value
            character(64) :: thisposition
            character(64) :: subroutine_name = 'interface_get_linkf_attribute'
        !%-----------------------------------------------------------------------------

        if (setting%Debug%File%interface) &
            write(*,"(3(A,i5),A)") '*** enter ' // trim(subroutine_name) // &
            "(link_idx=", link_idx, ", attr=", attr, ")" // " [Processor ", this_image(), "]"

        if ((attr > N_api_total_linkf_attributes) .or. (attr < 1)) then
            print *, "error: unexpected link attribute value", attr
            stop 49870555
        end if

        if ((link_idx > SWMM_N_link) .or. (link_idx < 1)) then
            print *, "error: unexpected link index value", link_idx
            stop 9987355
        end if

        if (attr <= N_api_linkf_attributes) then
            ! Fortran index starts in 1, whereas in C starts in 0
            call load_api_procedure("api_get_linkf_attribute")
            error = ptr_api_get_linkf_attribute(link_idx-1, attr, link_value)
            thisposition = trim(subroutine_name)//'_A01'
            call print_api_error(error, thisposition)
        else if ( (attr > N_api_linkf_attributes) .and. &
                  (attr <= (N_api_linkf_attributes + N_api_linkf_type_attributes)) )then
            call load_api_procedure("api_get_linkf_attribute")
            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_type, link_value)
            thisposition = trim(subroutine_name)//'_B02'
            call print_api_error(error, thisposition)
            
            if (link_value == API_CONDUIT) then
                if (attr == api_linkf_type) then
                    link_value = lPipe
                else if (attr == api_linkf_weir_type) then
                    !write(*,*) '****** Not clear why this is null brh20211207 at 837894'
                    link_value = nullvalueI
                else if (attr == api_linkf_orifice_type) then
                    !write(*,*) '****** Not clear why this is null brh20211207 at 938765'
                    link_value = nullvalueI
                else if (attr == api_linkf_pump_type) then
                    !write(*,*) '****** Not clear why this is null brh20211207 at 875323'
                    link_value = nullvalueI
                else if (attr == api_linkf_outlet_type) then
                    !write(*,*) '****** Not clear why this is null brh20211207 at 3897895'
                    link_value = nullvalueI
                else
                    link_value = nullvalueI
                    write(*,*)
                    write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at 45782987'
                    write(*,*) '   attr = ',attr
                    write(*,*) '   allowable are ',api_linkf_type, api_linkf_weir_type, api_linkf_orifice_type, api_linkf_pump_type, api_linkf_outlet_type
                    write(*,*) '   skipping error condition!'
                    write(*,*) '******'
                end if

            else if (link_value == API_PUMP) then
                if (attr == api_linkf_type) then
                    link_value = lPump
                else if (attr == api_linkf_weir_type) then
                    write(*,*) '****** Not clear why this is null brh20211207 at 478974'
                    link_value = nullvalueI
                else if (attr == api_linkf_orifice_type) then
                    write(*,*) '****** Not clear why this is null brh20211207 at 3695'
                    link_value = nullvalueI
                else if (attr == api_linkf_pump_type) then
                    call load_api_procedure("api_get_linkf_attribute")
                    error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_pump_type, link_value)
                    thisposition = trim(subroutine_name)//'_C03'
                    call print_api_error(error, thisposition)
                    if (link_value == API_TYPE1_PUMP) then
                        link_value = lType1Pump
                    else if (link_value == API_TYPE2_PUMP) then
                        link_value = lType2Pump
                    else if (link_value == API_TYPE3_PUMP) then
                        link_value = lType3Pump
                    else if (link_value == API_TYPE4_PUMP) then
                        link_value = lType4Pump
                    else if (link_value == API_IDEAL_PUMP) then
                        link_value = lTypeIdealPump   
                    else
                        link_value = nullvalueI
                        write(*,*)
                        write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at  442789'
                        write(*,*) '   link_value = ',link_value
                        write(*,*) '   allowable are ',API_TYPE1_PUMP, API_TYPE2_PUMP, API_TYPE3_PUMP, API_TYPE4_PUMP, API_IDEAL_PUMP
                        write(*,*) '   skipping error condition!'
                        write(*,*) '******'                     
                    endif
                else if (attr == api_linkf_outlet_type) then
                    write(*,*) '****** Not clear why this is null brh20211207 at 773682'
                    link_value = nullvalueI
                !% brh20211207s    
                else
                    link_value = nullvalueI
                    write(*,*)
                    write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at  8836785'
                    write(*,*) '   attr = ',attr
                    write(*,*) '   allowable are ',api_linkf_type, api_linkf_weir_type, api_linkf_orifice_type, api_linkf_pump_type, api_linkf_outlet_type
                    write(*,*) '   skipping error condition!'
                    write(*,*) '******'          
                end if

            else if (link_value == API_ORIFICE) then
                if (attr == api_linkf_type) then
                    link_value = lOrifice
                else if (attr == api_linkf_weir_type) then
                    write(*,*) '****** Not clear why this is null brh20211207 at 82674'
                    link_value = nullvalueI
                else if (attr == api_linkf_orifice_type) then
                    call load_api_procedure("api_get_linkf_attribute")
                    error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_orifice_type, link_value)
                    thisposition = trim(subroutine_name)//'_D04'
                    call print_api_error(error, thisposition)
                    if (link_value == API_SIDE_ORIFICE) then
                        link_value = lSideOrifice
                    else if (link_value == API_BOTTOM_ORIFICE) then
                        link_value = lBottomOrifice  
                    else
                        link_value = nullvalueI
                        write(*,*)
                        write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at  442789'
                        write(*,*) '   link_value = ',link_value
                        write(*,*) '   allowable are ',API_SIDE_ORIFICE, API_BOTTOM_ORIFICE
                        write(*,*) '   skipping error condition!'
                        write(*,*) '******'     
                    endif
                else if (attr == api_linkf_pump_type) then
                    write(*,*) '****** Not clear why this is null brh20211207 at  8378904'
                    link_value = nullvalueI
                else if (attr == api_linkf_outlet_type) then
                    write(*,*) '****** Not clear why this is null brh20211207 at 3446789'
                    link_value = nullvalueI  
                else
                    link_value = nullvalueI
                    write(*,*)
                    write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at  9937685'
                    write(*,*) '   attr = ',attr
                    write(*,*) '   allowable are ',api_linkf_type, api_linkf_weir_type, api_linkf_orifice_type, api_linkf_pump_type, api_linkf_outlet_type
                    write(*,*) '   skipping error condition!'
                    write(*,*) '******'              
                end if

            else if (link_value == API_WEIR) then
                if (attr == api_linkf_type) then
                    link_value = lWeir
                else if (attr == api_linkf_weir_type) then
                    call load_api_procedure("api_get_linkf_attribute")
                    error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_weir_type, link_value)
                    thisposition = trim(subroutine_name)//'_E05'
                    call print_api_error(error, thisposition)
                    if (link_value == API_TRANSVERSE_WEIR) then
                        link_value = lTransverseWeir
                    else if (link_value == API_SIDEFLOW_WEIR) then
                        link_value = lSideFlowWeir
                    else if (link_value == API_VNOTCH_WEIR) then
                        link_value = lVnotchWeir
                    else if (link_value == API_TRAPEZOIDAL_WEIR) then
                        link_value = lTrapezoidalWeir
                    else if (link_value == API_ROADWAY_WEIR) then
                        link_value = lRoadWayWeir  
                    else
                        link_value = nullvalueI
                        write(*,*)
                        write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at 620873'
                        write(*,*) '   link_value = ',link_value
                        write(*,*) '   allowable are ',API_TRANSVERSE_WEIR, API_SIDEFLOW_WEIR, API_VNOTCH_WEIR, API_TRAPEZOIDAL_WEIR, API_ROADWAY_WEIR
                        write(*,*) '   skipping error condition!'
                        write(*,*) '******'  
                    endif
                else if (attr == api_linkf_orifice_type) then
                    write(*,*) '****** Not clear why this is null brh20211207 at 836637'
                    link_value = nullvalueI
                else if (attr == api_linkf_pump_type) then
                    write(*,*) '****** Not clear why this is null brh20211207 at 6637892'
                    link_value = nullvalueI
                else if (attr == api_linkf_outlet_type) then
                    write(*,*) '****** Not clear why this is null brh20211207 at 826784'
                    link_value = nullvalueI    
                else
                    link_value = nullvalueI
                    write(*,*)
                    write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at  9937685'
                    write(*,*) '   attr = ',attr
                    write(*,*) '   allowable are ',api_linkf_type, api_linkf_weir_type, api_linkf_orifice_type, api_linkf_pump_type, api_linkf_outlet_type
                    write(*,*) '   skipping error condition!'
                    write(*,*) '******'   
                end if

            else if (link_value == API_OUTLET) then
                if (attr == api_linkf_type) then
                    link_value = lOutlet
                else if (attr == api_linkf_weir_type) then
                    write(*,*) '****** Not clear why this is null brh20211207 at 226784'
                    link_value = nullvalueI
                else if (attr == api_linkf_orifice_type) then
                    write(*,*) '****** Not clear why this is null brh20211207 at 737890'
                    link_value = nullvalueI
                else if (attr == api_linkf_pump_type) then
                    write(*,*) '****** Not clear why this is null brh20211207 at 827894'
                    link_value = nullvalueI
                else if (attr == api_linkf_outlet_type) then
                    call load_api_procedure("api_get_linkf_attribute")
                    error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_outlet_type, link_value)
                    thisposition = trim(subroutine_name)//'_F05'
                    call print_api_error(error, thisposition)
                    if (link_value == API_NODE_DEPTH) then
                        link_value = lNodeDepth
                    else if (link_value == API_NODE_HEAD) then
                        link_value = lNodeHead 
                    else
                        link_value = nullvalueI
                        write(*,*)
                        write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at 220455'
                        write(*,*) '   link_value = ',link_value
                        write(*,*) '   allowable are ',API_NODE_DEPTH, API_NODE_HEAD
                        write(*,*) '   skipping error condition!'
                        write(*,*) '******'                       
                    endif
                else
                    link_value = nullvalueI
                    write(*,*)
                    write(*,*) '****** Unexpected else in ',trim(subroutine_name),' at  11947'
                    write(*,*) '   attr = ',attr
                    write(*,*) '   allowable are ',api_linkf_type, api_linkf_weir_type, api_linkf_orifice_type, api_linkf_pump_type, api_linkf_outlet_type
                    write(*,*) '   skipping error condition!'
                    write(*,*) '******'     
                end if                

            endif

        else
            call load_api_procedure("api_get_linkf_attribute")
            error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_type, link_value)
            thisposition = trim(subroutine_name)//'_E05'
            call print_api_error(error, thisposition)
            if (link_value == API_RECT_CLOSED) then
                if (attr == api_linkf_geometry) then
                    link_value = lRectangular_closed
                else if (attr == api_linkf_xsect_wMax) then
                    call load_api_procedure("api_get_linkf_attribute")
                    error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                    thisposition = trim(subroutine_name)//'_F06'
                    call print_api_error(error, thisposition)
                else if (attr == api_linkf_xsect_yFull) then
                    call load_api_procedure("api_get_linkf_attribute")
                    error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                    thisposition = trim(subroutine_name)//'_G07'
                    call print_api_error(error, thisposition)
                else
                    write(*,*) '****** Not clear why this is null brh20211207 at 87438687'
                    link_value = nullvalueR
                end if
            else if (link_value == API_RECT_OPEN) then
                if (attr == api_linkf_geometry) then
                    link_value = lRectangular
                else if (attr == api_linkf_xsect_wMax) then
                    call load_api_procedure("api_get_linkf_attribute")
                    error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                    thisposition = trim(subroutine_name)//'_H08'
                    call print_api_error(error, thisposition)
                else if (attr == api_linkf_xsect_yFull) then
                    call load_api_procedure("api_get_linkf_attribute")
                    error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                    thisposition = trim(subroutine_name)//'_I09'
                    call print_api_error(error, thisposition)
                else
                    write(*,*) '****** Not clear why this is null brh20211207 at 72678439'
                    link_value = nullvalueR
                end if
            else if (link_value == API_TRAPEZOIDAL) then
                if (attr == api_linkf_geometry) then
                    link_value = lTrapezoidal
                else if (attr == api_linkf_xsect_wMax) then
                    call load_api_procedure("api_get_linkf_attribute")
                    error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yBot, link_value)
                    thisposition = trim(subroutine_name)//'_J10'
                    call print_api_error(error, thisposition)
                else if (attr == api_linkf_xsect_yFull) then
                    call load_api_procedure("api_get_linkf_attribute")
                    error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                    thisposition = trim(subroutine_name)//'_K11'
                    call print_api_error(error, thisposition)
                else
                    write(*,*) '****** Not clear why this is null brh20211207 at 55278093'
                    link_value = nullvalueR
                end if
            else if (link_value == API_TRIANGULAR) then
                if (attr == api_linkf_geometry) then
                    link_value = lTriangular
                else if (attr == api_linkf_xsect_wMax) then
                    call load_api_procedure("api_get_linkf_attribute")
                    error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                    thisposition = trim(subroutine_name)//'_M12'
                    call print_api_error(error, thisposition)
                else if (attr == api_linkf_xsect_yFull) then
                    call load_api_procedure("api_get_linkf_attribute")
                    error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                    thisposition = trim(subroutine_name)//'_N13'
                    call print_api_error(error, thisposition)
                else
                    write(*,*) '****** Not clear why this is null brh20211207 at 267849'
                    link_value = nullvalueR
                end if
            else if (link_value == API_PARABOLIC) then
                if (attr == api_linkf_geometry) then
                    link_value = lParabolic
                else if (attr == api_linkf_xsect_wMax) then
                    call load_api_procedure("api_get_linkf_attribute")
                    error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                    thisposition = trim(subroutine_name)//'_O14'
                    call print_api_error(error, thisposition)
                else if (attr == api_linkf_xsect_yFull) then
                    call load_api_procedure("api_get_linkf_attribute")
                    error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                    thisposition = trim(subroutine_name)//'_P15'
                    call print_api_error(error, thisposition)
                else
                    write(*,*) '****** Not clear why this is null brh20211207 at 8270894987'
                    link_value = nullvalueR
                end if
            else if (link_value == API_CIRCULAR) then
                if (attr == api_linkf_geometry) then
                    link_value = lCircular
                else if (attr == api_linkf_xsect_wMax) then
                    call load_api_procedure("api_get_linkf_attribute")
                    error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_wMax, link_value)
                    thisposition = trim(subroutine_name)//'_Q16'
                    call print_api_error(error, thisposition)
                else if (attr == api_linkf_xsect_yFull) then
                    call load_api_procedure("api_get_linkf_attribute")
                    error = ptr_api_get_linkf_attribute(link_idx-1, api_linkf_xsect_yFull, link_value)
                    thisposition = trim(subroutine_name)//'_R17'
                    call print_api_error(error, thisposition)
                else
                    write(*,*) '****** Not clear why this is null brh20211207 at 8276438'
                    link_value = nullvalueR
                end if
            else
                write(*,*) '****** Not clear why this is null brh20211207 at 227589'
                link_value = nullvalueR
            end if
        end if
        if (setting%Debug%File%interface)  then
            write(*,"(3(A,i5),A)") '*** leave ' // trim(subroutine_name) // &
            "(link_idx=", link_idx, ", attr=", attr, ")" // " [Processor ", this_image(), "]"
        end if
    end function interface_get_linkf_attribute
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_table_attribute(table_idx, attr)
        !%-----------------------------------------------------------------------------
        !% Description:
        !%    Retrieves table attributes from EPA-SWMM. API table attributes are
        !%    defined in define_api_keys.f08.
        !% Notes:
        !%    Fortran indexes are translated to C indexes and viceversa when
        !%    necessary. Fortran indexes always start from 1, whereas C indexes
        !%    start from 0.
        !%-----------------------------------------------------------------------------
            integer :: table_idx, attr, error
            real(8) :: interface_get_table_attribute

            real(c_double), target :: table_value
            character(64) :: thisposition
            character(64) :: subroutine_name = 'interface_get_table_attribute'
        !%-----------------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if ((attr > N_api_total_table_attributes) .or. (attr < 1)) then
            print *, "error: unexpected table attribute value", attr
            stop 28704
        end if

        if ((table_idx > SWMM_N_Curve) .or. (table_idx < 1)) then
            print *, "error: unexpected table index value", table_idx
            stop 498734
        end if

        !% Substracts 1 to every Fortran index (it becomes a C index)
        if (attr == api_table_type) then
            call load_api_procedure("api_get_table_attribute")
            error = ptr_api_get_table_attribute(table_idx-1, attr, table_value)
            if (table_value == API_STORAGE_CURVE) then
                interface_get_table_attribute = StorageCurve
            else if (table_value == API_DIVERSION_CURVE) then
                interface_get_table_attribute = DiversionCurve
            else if (table_value == API_TIDAL_CURVE) then
                interface_get_table_attribute = TidalCurve
            else if (table_value == API_RATING_CURVE) then
                interface_get_table_attribute = RatingCurve
            else if (table_value == API_CONTROL_CURVE) then
                interface_get_table_attribute = ControlCurve
            else if (table_value == API_SHAPE_CURVE) then
                interface_get_table_attribute = ShapeCurve
            else if (table_value == API_WEIR_CURVE) then
                interface_get_table_attribute = WeirCurve
            else if (table_value == API_PUMP1_CURVE) then
                interface_get_table_attribute = Pump1Curve
            else if (table_value == API_PUMP2_CURVE) then
                interface_get_table_attribute = Pump2Curve
            else if (table_value == API_PUMP3_CURVE) then
                interface_get_table_attribute = Pump3Curve
            else if (table_value == API_PUMP4_CURVE) then
                interface_get_table_attribute = Pump4Curve
            else
                write(*,*) '****** Not clear why this is null brh20211207 at 873984'
                interface_get_table_attribute = nullvalueI
            end if
        else
            call load_api_procedure("api_get_table_attribute")
            error = ptr_api_get_table_attribute(table_idx-1, attr, table_value)
            call print_api_error(error, subroutine_name)
            interface_get_table_attribute = table_value
        end if

        if (setting%Debug%File%interface)  then
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            print *, "table", table_value, attr
        end if
    end function interface_get_table_attribute
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_num_table_entries(table_idx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !%    Retrieves table attributes from EPA-SWMM. API table attributes are
        !%    defined in define_api_keys.f08.
        !% Notes:
        !%    Fortran indexes are translated to C indexes and viceversa when
        !%    necessary. Fortran indexes always start from 1, whereas C indexes
        !%    start from 0.
        !%-----------------------------------------------------------------------------
            integer :: table_idx, table_type, error
            integer :: interface_get_num_table_entries
            integer(c_int), target :: table_entries
            character(64) :: thisposition
            character(64) :: subroutine_name = 'interface_get_num_table_entries'
        !%-----------------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if ((table_idx > SWMM_N_Curve) .or. (table_idx < 1)) then
            print *, "error: unexpected table index value", table_idx
            stop 835551
        end if

        !% Substracts 1 to every Fortran index (it becomes a C index)
        call load_api_procedure("api_get_num_table_entries")
        error = ptr_api_get_num_table_entries(table_idx-1, API_CURVE, table_entries)
        call print_api_error(error, subroutine_name)
        interface_get_num_table_entries = table_entries
        if (setting%Debug%File%interface)  then
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            print *, "table", table_entries
        end if
        
    end function interface_get_num_table_entries
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_first_entry_table(table_idx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !%    Retrieves first table entries from EPA-SWMM. API table attributes are
        !%    defined in define_api_keys.f08.
        !% Notes:
        !%    Fortran indexes are translated to C indexes and viceversa when
        !%    necessary. Fortran indexes always start from 1, whereas C indexes
        !%    start from 0.
        !%-----------------------------------------------------------------------------
            integer :: table_idx, error, success
            real(8) :: interface_get_first_entry_table(2)

            real(c_double), target :: x_entry, y_entry
            character(64) :: subroutine_name = 'interface_get_first_entry_table'
        !%-----------------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if ((table_idx > SWMM_N_Curve) .or. (table_idx < 1)) then
            print *, "error: unexpected table index value", table_idx
            stop 837014
        end if

        !% Substracts 1 to every Fortran index (it becomes a C index)
        call load_api_procedure("api_get_first_entry_table")
        success = ptr_api_get_first_entry_table(table_idx-1, API_CURVE, x_entry, y_entry)

        if (success == 0) then
            error = -1
        else
            error = 0
        end if

        call print_api_error(error, subroutine_name)
        interface_get_first_entry_table(1) = x_entry
        interface_get_first_entry_table(2) = y_entry

        if (setting%Debug%File%interface)  then
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        end if
    end function interface_get_first_entry_table
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_next_entry_table(table_idx, table_type)
        !%-----------------------------------------------------------------------------
        !% Description:
        !%    Retrieves next table entries from EPA-SWMM. API table attributes are
        !%    defined in define_api_keys.f08.
        !% Notes:
        !%    Fortran indexes are translated to C indexes and viceversa when
        !%    necessary. Fortran indexes always start from 1, whereas C indexes
        !%    start from 0.
        !%-----------------------------------------------------------------------------
        integer :: table_idx, table_type, error, success
        real(8) :: interface_get_next_entry_table(2)

        real(c_double), target :: x_entry, y_entry
        character(64) :: subroutine_name
        !%-----------------------------------------------------------------------------
        subroutine_name = 'interface_get_next_entry_table'
        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if ((table_idx > SWMM_N_Curve) .or. (table_idx < 1)) then
            print *, "error: unexpected table index value", table_idx
            stop 8367894
        end if

        !% Substracts 1 to every Fortran index (it becomes a C index)
        call load_api_procedure("api_get_next_entry_table")
        success = ptr_api_get_next_entry_table(table_idx-1, table_type, x_entry, y_entry)

        if (success == 0) then
            error = -1
        else
            error = 0
        end if

        call print_api_error(error, subroutine_name)
        interface_get_next_entry_table(1) = x_entry
        interface_get_next_entry_table(2) = y_entry

        if (setting%Debug%File%interface)  then
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        end if
    end function interface_get_next_entry_table
!%
!%=============================================================================
!%   Boundary Conditions (execute after initialization only)
!%=============================================================================
!%
    function interface_get_BC_resolution(node_idx) result(resolution)
        !%-----------------------------------------------------------------------------
        !% Description:
        !%    Computes the finest pattern resolution associated with a node's BC.
        !%    The resulting pattern type is stored in node%I(:,ni_pattern_resolution)
        !%    and is reused to fetch inflow BCs such that the amount of inflow points
        !%    is minimized.
        !% Notes:
        !%    * Currently patterns are only associated with inflow BCs. It is necessary
        !%      to preserve the order of api_monthly, api_weekend, api_daily, and
        !%      api_hourly in define_api_index.f08 for the function to work.
        !%    * The function is called during the intialization of the node%I table
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: node_idx
        integer             :: p0, p1, p2, p3, p4
        integer             :: resolution
        real(8)             :: baseline
        character(64)       :: subroutine_name
        !%-----------------------------------------------------------------------------
        subroutine_name = 'interface_get_BC_resolution'
        if (setting%Debug%File%interface)  &
        write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name)  //  &
            " [Processor ", this_image(), "]"

        resolution = nullvalueI

        if (node%YN(node_idx, nYN_has_inflow)) then ! Upstream/Lateral BC

            resolution = 0

            if (node%YN(node_idx, nYN_has_extInflow)) then

                !write(*,*) 'api_nodef_extInflow_basePat_type',api_nodef_extInflow_basePat_type
                p0 = interface_get_nodef_attribute(node_idx, api_nodef_extInflow_basePat_type)
                !write(*,*), p0
                !write(*,*), api_hourly_pattern, api_weekend_pattern, api_daily_pattern, api_monthly_pattern

                if (p0 == api_hourly_pattern) then
                    resolution = api_hourly
                else if (p0 == api_weekend_pattern) then
                    resolution = api_weekend
                else if (p0 == api_daily_pattern) then
                    resolution = api_daily
                else if (p0 == api_monthly_pattern) then
                    !write(*,*) 'api_nodef_extInflow_baseline',api_nodef_extInflow_baseline
                    baseline = interface_get_nodef_attribute(node_idx, api_nodef_extInflow_baseline)
                    if (baseline > 0) resolution = api_monthly
                else
                    !write(*,*) '--- no external inflow timeseries pattern for node ',node_idx
                    !write(*,*) '****** unexpected else in ',trim(subroutine_name), ' at 9873094'
                    !write(*,*), '   p0 read is ',p0
                    !write(*,*), '   p0 allowed are ',api_hourly_pattern, api_weekend_pattern, api_daily_pattern, api_monthly_pattern
                    !write(*,*), '   skipping error condition!'
                    !write(*,*) '******'  
                end if
            else
                !write(*,*) '--- no external inflows to node ',node_idx
            end if

            if (node%YN(node_idx, nYN_has_dwfInflow)) then

                write(*,*) 'call api_nodef_dwfInflow_hourly_pattern'
                p1 = interface_get_nodef_attribute(node_idx, api_nodef_dwfInflow_hourly_pattern)
                write(*,*) '   p1 = ',p1

                write(*,*) 'api_nodef_dwfInflow_weekend_pattern'
                p2 = interface_get_nodef_attribute(node_idx, api_nodef_dwfInflow_weekend_pattern)
                write(*,*) '   p2 = ',p2

                write(*,*) 'api_nodef_dwfInflow_daily_pattern'
                p3 = interface_get_nodef_attribute(node_idx, api_nodef_dwfInflow_daily_pattern)
                write(*,*) '   p3 = ',p3


                write(*,*) 'api_nodef_dwfInflow_monthly_pattern'
                p4 = interface_get_nodef_attribute(node_idx, api_nodef_dwfInflow_monthly_pattern)
                write(*,*) '   p4 = ',p4

                if (p1 > 0) then
                    resolution = max(api_hourly, resolution)
                else if (p2 > 0) then
                    resolution = max(api_weekend, resolution)
                else if (p3 > 0) then
                    resolution = max(api_daily, resolution)
                else if (p4 > 0) then
                    resolution = max(api_monthly, resolution)
                else
                    write(*,*) '***** unexpected else in ',trim(subroutine_name), 'at 3987044'
                    write(*,*), '   p0 read is ',p0
                    write(*,*), '   p0 allowed are ',api_hourly_pattern, api_weekend_pattern, api_daily_pattern, api_monthly_pattern
                    write(*,*), '   skipping error condition!'
                    write(*,*) '******'
                    !stop 3987044
                !% brh20211207e    
                end if
            !% brh20211208s
            else
                !write(*,*) '   no dwfInflow to this node'
            !% brh20211208e
            end if
        !% brh20211208s
        else
            !write(*,*) '   no dwf or ext inflows to this node'
        !% brh20211208e
        end if

    end function interface_get_BC_resolution
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_next_inflow_time(bc_idx, tnow) result(tnext)
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
        integer, intent(in) :: bc_idx
        real(8), intent(in) :: tnow
        real(8)             :: tnext, t1, t2, tnextp
        integer             :: nidx, nres, tseries, success
        character(64) :: subroutine_name 
        !%---------------------------------------------------------------------
        subroutine_name = 'interface_get_next_inflow_time'
        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        nidx = BC%flowI(bc_idx, bi_node_idx)
        if (.not. node%YN(nidx, nYN_has_inflow)) then
            print *, "Error, node " // node%Names(nidx)%str // " does not have an inflow"
        end if
        nres = node%I(nidx, ni_pattern_resolution)
        if (nres >= 0) then
            tnextp = util_datetime_get_next_time(tnow, nres)
            if (node%YN(nidx, nYN_has_extInflow)) then
                tseries = interface_get_nodef_attribute(nidx, api_nodef_extInflow_tSeries)
                if (tseries >= 0) then
                    success = get_next_entry_tseries(tseries)
                    tnext = interface_get_nodef_attribute(nidx, api_nodef_extInflow_tSeries_x1)
                    tnext = util_datetime_epoch_to_secs(tnext)
                    if (success == 0) then ! unsuccessful
                        tnext = interface_get_nodef_attribute(nidx, api_nodef_extInflow_tSeries_x2)
                        tnext = util_datetime_epoch_to_secs(tnext)
                        if (tnext == tnow) then
                            tnext = setting%Time%End
                            setting%BC%disableInterpolationYN = .true.
                        !% brh20211207s
                        else
                            ! write(*,*)
                            ! write(*,*) '****** Unexpected else in ',subroutine_name,' at 4479823'
                            ! write(*,*) '   tnext = ', tnext
                            ! write(*,*) '   tnow  = ', tnow
                            ! write(*,*) '   skipping error condition'
                            ! write(*,*) '******'
                            ! stop 4479823
                        !% brh20211207e    
                        end if
                    !% brh20211207s
                    else
                        ! write(*,*)
                        ! write(*,*) '****** Unexpected else in ',subroutine_name,' at 4589709'
                        ! write(*,*) '   success = ', success
                        ! write(*,*) '   skipping error condition'
                        ! write(*,*) '******'
                        ! stop 4589709
                    !% brh20211207e                          
                    end if
                else
                    tnext = setting%Time%End
                end if
            end if
            tnext = min(tnext, tnextp)
        else
            tnext = setting%Time%End
        end if

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end function interface_get_next_inflow_time
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_next_head_time(bc_idx, tnow) result(tnext)
        !%---------------------------------------------------------------------
        !% Description
        !%---------------------------------------------------------------------
        integer, intent(in) :: bc_idx
        real(8), intent(in) :: tnow
        real(8)             :: tnext, tnextp
        integer             :: nidx, nres, tseries
        character(64) :: subroutine_name
        !%---------------------------------------------------------------------
        subroutine_name = 'interface_get_next_head_time'
        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        nidx = BC%headI(bc_idx, bi_node_idx)
        if (BC%headI(bc_idx, bi_subcategory) == BCH_fixed) then
            tnext = setting%Time%End
        else
            print *, "Error, unsupported head boundary condition for node " // node%Names(nidx)%str
            stop
        end if

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end function interface_get_next_head_time
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_flowBC(bc_idx, tnow) result(bc_value)
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
        integer, intent(in) :: bc_idx
        real(8), intent(in) :: tnow
        integer             :: error, nidx
        real(8)             :: epochNow, bc_value
        character(64) :: subroutine_name
        !%---------------------------------------------------------------------
        subroutine_name = 'interface_get_flowBC'
        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        nidx = BC%flowI(bc_idx, bi_node_idx)
        epochNow = util_datetime_secs_to_epoch(tnow)
        call load_api_procedure("api_get_flowBC")
        error = ptr_api_get_flowBC(nidx-1, epochNow, bc_value)
        call print_api_error(error, subroutine_name)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end function interface_get_flowBC
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_headBC(bc_idx, tnow) result(bc_value)
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
        integer, intent(in) :: bc_idx
        real(8), intent(in) :: tnow
        integer             :: error, nidx
        real(8)             :: epochNow, bc_value
        character(64) :: subroutine_name
        !%---------------------------------------------------------------------
        subroutine_name = 'interface_get_headBC'
        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        nidx = BC%headI(bc_idx, bi_node_idx)
        epochNow = util_datetime_secs_to_epoch(tnow)
        call load_api_procedure("api_get_headBC")
        error = ptr_api_get_headBC(nidx-1, epochNow, bc_value)
        call print_api_error(error, subroutine_name)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end function interface_get_headBC
!%
!%=============================================================================
!%    Write Outputs (execute after initialization only)
!%=============================================================================
!%
    subroutine interface_export_link_results(link_idx)
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
        integer, intent(in) :: link_idx
        integer :: error
        character(64) :: subroutine_name = 'interface_export_link_results'
        !%---------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        call load_api_procedure("api_export_link_results")
        error = ptr_api_export_link_results(link_idx-1)
        call print_api_error(error, subroutine_name)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine interface_export_link_results
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_update_nodeResult(node_idx, result_type, node_result)
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
               integer, intent(in) :: node_idx, result_type
        real(8), intent(in) :: node_result
        integer             :: error
        character(64)       :: subroutine_name = "interface_update_nodeResult"
        !%----------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        call load_api_procedure("api_update_nodeResult")
        error = ptr_api_update_nodeResult(node_idx-1, result_type, node_result)
        call print_api_error(error, subroutine_name)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine interface_update_nodeResult
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_update_linkResult(link_idx, result_type, link_result)
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
        integer, intent(in) :: link_idx, result_type
        real(8), intent(in) :: link_result
        integer             :: error
        character(64)       :: subroutine_name = "interface_update_linkResult"
        !%----------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        call load_api_procedure("api_update_linkResult")
        error = ptr_api_update_linkResult(link_idx-1, result_type, link_result)
        call print_api_error(error, subroutine_name)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine interface_update_linkResult
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_write_output_line(reportTime)
        !%--------------------------------------------------------------------------
        !% Description:
        !%    Writes .out file with SWMM5+ data
        !%--------------------------------------------------------------------------
            real(c_double), intent(in) :: reportTime ! time in seconds
            integer                    :: error
            character(64)              :: subroutine_name = "interface_write_output_line"
        !%--------------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        call load_api_procedure("api_write_output_line")
        error = ptr_api_write_output_line(reportTime)
        call print_api_error(error, subroutine_name)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine interface_write_output_line
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_get_SWMM_times()
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
        integer                :: error
        real(c_double), target :: reportStart_epoch, RouteStep 
        real(c_double), target :: starttime_epoch, endtime_epoch
        real(c_double), target :: TotalDuration
        integer(c_int), target :: reportStep, WetStep, DryStep
        character(64)          :: subroutine_name = 'interface_get_SWMM_times'
        !%----------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% --- reportStart is given in epoch datetime
        !% --- reportStep and hydroStep are given in integer seconds
        !%
        !% --- Note that if reportStart is earlier than StartEpoch from get_start_datetime()
        !% ... then SWMM will make reportStart = start datetime.
        !%
        !% --- This generates EPA-SWMM Error Code: -1 if the reportStart is after the end time.
        call load_api_procedure("api_get_SWMM_times")
        error = ptr_api_get_SWMM_times( &
            starttime_epoch,  &
            endtime_epoch,    &
            reportStart_epoch, &
            reportStep,  &
            WetStep,     &
            DryStep,     &
            RouteStep,   &
            TotalDuration)
        call print_api_error(error, subroutine_name)

        setting%SWMMinput%StartEpoch            = starttime_epoch
        setting%SWMMinput%EndEpoch              = endtime_epoch                      
        setting%SWMMinput%ReportStartTimeEpoch  = reportStart_epoch
        setting%SWMMinput%ReportTimeInterval    = reportStep
        setting%SWMMinput%WetStep               = WetStep
        setting%SWMMinput%RouteStep             = RouteStep
        setting%SWMMinput%DryStep               = DryStep
        setting%SWMMinput%TotalDuration         = TotalDuration

        !print *, '............'
        !print *, 'in ',trim(subroutine_name), starttime_epoch, endtime_epoch, setting%SWMMinput%EndEpoch
        !print *, '............'

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine interface_get_SWMM_times
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_find_object(object_type, object_name) result(object_idx)
        !%---------------------------------------------------------------------
        !% Description:
        !% Returns the index of the object, or 0 if the object couldn't be found
        !%---------------------------------------------------------------------
        character(*), intent(in) :: object_name
        integer, intent(in) :: object_type
        integer :: object_idx
        character(64) :: subroutine_name = 'interface_find_object'
        !%---------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        call load_api_procedure("api_find_object")
        object_idx = ptr_api_find_object(object_type, trim(object_name)//c_null_char) + 1

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end function interface_find_object
!%
!%=============================================================================
!%=============================================================================
!%   
    subroutine interface_call_runoff_execute()
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
        integer             :: error
        character(64) :: subroutine_name = 'interface_call_runoff_execute'     
        !%---------------------------------------------------------------------   
        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"
        
        ! print *, '...just before loading runoff_execute in interface_call_runoff_execute'

        call load_api_procedure("api_call_runoff_execute")
        error = ptr_api_call_runoff_execute()   
        call print_api_error(error, subroutine_name)
    
        ! print *, "...finishing interface_call_runoff_execute"
        
        if (setting%Debug%File%interface)  &
             write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
        
    end subroutine interface_call_runoff_execute
!%
!%=============================================================================
!%=============================================================================
!%   
    real(8) function interface_get_subcatch_runoff (sc_idx) result(runoff)
        !%---------------------------------------------------------------------
        !% Description:
        !% Gets the SWMM-C "newRunoff" for j=sc_idx for the data
        !% Subcatch[j].newRunoff in subcatch.c/subcatch_getRunoff
        !%---------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: sc_idx  !% index of this subcatchment
            integer             :: error
            character(64) :: subroutine_name = 'interface_get_subcatch_runoff'     
        !%---------------------------------------------------------------------   
        !% Preliminaries
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"
        !%---------------------------------------------------------------------

        call load_api_procedure("api_get_subcatch_runoff")
        error = ptr_api_get_subcatch_runoff(sc_idx, runoff)   
        call print_api_error(error, subroutine_name)

        !%---------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%interface)  &
             write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
        
    end function interface_get_subcatch_runoff
!%
!%=============================================================================
!%=============================================================================
!%  
    integer function interface_get_subcatch_runoff_nodeIdx (sc_idx) result(node_idx)    
        !%---------------------------------------------------------------------
        !% Description:
        !% Gets the SWMM-C node index that a subcatchment runs off to
        !%---------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: sc_idx  !% index of this subcatchment
            integer             :: error
            character(64) :: subroutine_name = 'interface_get_subcatch_runoff_nodeIdx'     
        !%---------------------------------------------------------------------   
        !% Preliminaries
            if (setting%Debug%File%interface)  &
                write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"
        !%---------------------------------------------------------------------

        call load_api_procedure("api_get_subcatch_runoff_nodeIdx")
        error = ptr_api_get_subcatch_runoff_nodeIdx(sc_idx, node_idx)   
        call print_api_error(error, subroutine_name)
        
        !%---------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end function interface_get_subcatch_runoff_nodeIdx
!%
!%=============================================================================
!% PRIVATE
!%=============================================================================
!%
    subroutine load_api_procedure(api_procedure_name)
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
        character(kind=c_char, len=*) :: api_procedure_name
        character(64) :: subroutine_name = 'load_api_procedure'
        !%---------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name)  // " (" // api_procedure_name // ")" // &
                " [Processor ", this_image(), "]"

        !% Open the shared library
        call c_lib_open(c_lib, errstat, errmsg)
        c_lib%procname = api_procedure_name
        call c_lib_load(c_lib, errstat, errmsg)

        !% Loads shared library funcitonalities
        select case (api_procedure_name)
            case ("api_initialize")
                call c_f_procpointer(c_lib%procaddr, ptr_api_initialize)
            case ("api_finalize")
                print *, '9781053 calling api_finalize'
                call c_f_procpointer(c_lib%procaddr, ptr_api_finalize)
            case ("api_get_nodef_attribute")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_nodef_attribute)
            case ("api_get_linkf_attribute")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_linkf_attribute)
            case ("api_get_table_attribute")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_table_attribute)
            case ("api_get_num_objects")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_num_objects)
            case ("api_get_object_name_len")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_object_name_len)
            case ("api_get_object_name")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_object_name)
            case ("api_get_num_table_entries")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_num_table_entries)
            case ("api_get_first_entry_table")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_first_entry_table)
            case ("api_get_next_entry_table")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_next_entry_table)
            case ("api_get_start_datetime")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_start_datetime)
            case ("api_get_end_datetime")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_end_datetime)
            case ("api_get_flowBC")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_flowBC)
            case ("api_get_headBC")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_headBC)
            case ("api_get_SWMM_times")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_SWMM_times)
            case ("api_get_next_entry_tseries")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_next_entry_tseries)
            case ("api_find_object")
                call c_f_procpointer(c_lib%procaddr, ptr_api_find_object)
            case ("api_run_step")
                call c_f_procpointer(c_lib%procaddr, ptr_api_run_step)
            case ("api_export_link_results")
                call c_f_procpointer(c_lib%procaddr, ptr_api_export_link_results)
            case ("api_write_output_line")
                call c_f_procpointer(c_lib%procaddr, ptr_api_write_output_line)
            case ("api_update_nodeResult")
                call c_f_procpointer(c_lib%procaddr, ptr_api_update_nodeResult)
            case ("api_update_linkResult")
                call c_f_procpointer(c_lib%procaddr, ptr_api_update_linkResult)
            case ("api_call_runoff_execute")    
                call c_f_procpointer(c_lib%procaddr, ptr_api_call_runoff_execute)
            case ("api_get_subcatch_runoff")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_subcatch_runoff)    
            case ("api_get_subcatch_runoff_nodeIdx")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_subcatch_runoff_nodeIdx)        
            case ("api_get_NewRunoffTime")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_NewRunoffTime)
            case default
                write(*,"(A,A)") "Error, procedure " // api_procedure_name // &
                 " has not been handled in load_api_procedure"
                stop
        end select

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name)  // " (" // api_procedure_name // ")" // &
                " [Processor ", this_image(), "]"
    end subroutine load_api_procedure
!%
!%=============================================================================
!%=============================================================================
!%
    function get_next_entry_tseries(k) result(success)
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
        integer, intent(in   ) :: k
        integer                :: success
        character(64)          :: subroutine_name
        !%---------------------------------------------------------------------
        subroutine_name = 'get_next_entry_tseries'
        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        call load_api_procedure("api_get_next_entry_tseries")
        success = ptr_api_get_next_entry_tseries(k-1) ! Fortran to C convention

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end function get_next_entry_tseries
!%
!%=============================================================================
!%=============================================================================
!%
    function get_num_objects(obj_type)
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
        integer :: obj_type
        integer :: get_num_objects
        character(64) :: subroutine_name
        !%---------------------------------------------------------------------
        subroutine_name = 'get_num_objects'
        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        call load_api_procedure("api_get_num_objects")
        get_num_objects = ptr_api_get_num_objects(obj_type)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end function get_num_objects
!%
!%=============================================================================
!%=============================================================================
!%
    function get_start_datetime()
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
        real(c_double) :: get_start_datetime
        character(64) :: subroutine_name
        !%---------------------------------------------------------------------
        subroutine_name = 'get_start_datetime'
        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        call load_api_procedure("api_get_start_datetime")
        get_start_datetime = ptr_api_get_start_datetime()

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end function get_start_datetime
!%
!%=============================================================================
!%=============================================================================
!%
    function get_end_datetime()
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
        real(8) :: get_end_datetime
        character(64) ::  subroutine_name
        !%---------------------------------------------------------------------
        subroutine_name = 'get_end_datetime'
        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        call load_api_procedure("api_get_end_datetime")
        get_end_datetime = ptr_api_get_end_datetime()

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end function get_end_datetime
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_NewRunoffTime()
        !%---------------------------------------------------------------------
        !% Description:
        !% gets the latest NewRunoffTime from SWMM-C (in milliseconds)
        !% and converts to seconds
        !%---------------------------------------------------------------------
        real(8)       :: interface_get_NewRunoffTime
        character(64) :: subroutine_name
        !%---------------------------------------------------------------------
        subroutine_name = 'interface_get_NewRunoffTime'
        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        call load_api_procedure("api_get_NewRunoffTime")
        interface_get_NewRunoffTime = ptr_api_get_NewRunoffTime() / onethousandR

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end function interface_get_NewRunoffTime
!%
!%=============================================================================
!%=============================================================================
!%
        
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine print_api_error(error, subroutine_name)
        !%---------------------------------------------------------------------
        !% Description:
        !%---------------------------------------------------------------------
        integer, intent(in) :: error
        character(64), intent(in) :: subroutine_name
        !%---------------------------------------------------------------------

        if (error /= 0) then
            write(*, "(A,i5,A)") new_line("") // "EPA-SWMM Error Code: ", error, " in "// subroutine_name &
                // " (see SWMM report file for more details)"
            write(*, "(A)") "Report File Path: " // trim(setting%File%rpt_file)
            stop 
        end if
    end subroutine print_api_error
!%
!%=============================================================================
!% END MODULE interface
!%=============================================================================
!%
end module interface
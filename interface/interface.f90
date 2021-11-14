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
    public :: interface_get_node_attribute
    public :: interface_get_link_attribute
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

!%==========================================================================
!% PRIVATE
!%==========================================================================

    ! Interface with SWMM shared library
    abstract interface
        !% -------------------------------------------------------------------------------
        !% Simulation
        !% -------------------------------------------------------------------------------
        function api_initialize(inp_file, report_file, out_file, run_routing)
            use, intrinsic :: iso_c_binding
            implicit none
            character(len=*,kind=c_char), intent(in) :: inp_file
            character(len=*,kind=c_char), intent(in) :: report_file
            character(len=*,kind=c_char), intent(in) :: out_file
            integer(c_int),        value, intent(in) :: run_routing
            integer(c_int) :: api_initialize
        end function api_initialize
        !% -------------------------------------------------------------------------------
        subroutine api_finalize()
            use, intrinsic :: iso_c_binding
            implicit none
        end subroutine api_finalize
        !% -------------------------------------------------------------------------------
        function api_run_step()
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double) :: api_run_step
        end function api_run_step
        !% -------------------------------------------------------------------------------
        ! --- Property-extraction
        !% -------------------------------------------------------------------------------
        ! * After Initialization
        function api_get_start_datetime()
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double) :: api_get_start_datetime
        end function api_get_start_datetime
        !% -------------------------------------------------------------------------------
        function api_get_end_datetime()
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double) :: api_get_end_datetime
        end function api_get_end_datetime
        !% -------------------------------------------------------------------------------
        function api_get_flowBC(k, current_datetime)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: k
            real(c_double),        intent(in) :: current_datetime
            real(c_double)                    :: api_get_flowBC
        end function api_get_flowBC
        !% -------------------------------------------------------------------------------
        function api_get_headBC(k, current_datetime)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: k
            real(c_double),        intent(in) :: current_datetime
            real(c_double)                    :: api_get_headBC
        end function api_get_headBC
        !% -------------------------------------------------------------------------------
        function api_get_report_times &
            (report_start_datetime, report_step, hydrology_step)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: report_start_datetime
            type(c_ptr), value, intent(in) :: report_step
            type(c_ptr), value, intent(in) :: hydrology_step
            integer(c_int) :: api_get_report_times
        end function api_get_report_times
        !% -------------------------------------------------------------------------------
        function api_get_node_attribute(k, attr, value)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: k
            integer(c_int), value, intent(in) :: attr
            type(c_ptr),    value, intent(in) :: value
            integer(c_int) :: api_get_node_attribute
        end function api_get_node_attribute
        !% -------------------------------------------------------------------------------
        function api_get_link_attribute(k, attr, value)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: k
            integer(c_int), value, intent(in) :: attr
            type(c_ptr),    value, intent(in) :: value
            integer(c_int) :: api_get_link_attribute
        end function api_get_link_attribute
        !% -------------------------------------------------------------------------------
        function api_get_table_attribute(k, attr, value)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: k
            integer(c_int), value, intent(in) :: attr
            type(c_ptr),    value, intent(in) :: value
            integer(c_int) :: api_get_table_attribute
        end function api_get_table_attribute
        !% -------------------------------------------------------------------------------
        integer function api_get_num_objects(obj_type)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: obj_type
        end function api_get_num_objects
        !% -------------------------------------------------------------------------------
        function api_get_object_name_len(k, object_type, len_value)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: k
            integer(c_int), value, intent(in) :: object_type
            type(c_ptr),    value, intent(in) :: len_value
            integer(c_int)                    :: api_get_object_name_len
        end function api_get_object_name_len
        !% -------------------------------------------------------------------------------
        function api_get_object_name(k, object_name, object_type)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int),           value, intent(in) :: k
            character(c_char), dimension(*), intent(in) :: object_name
            integer(c_int),           value, intent(in) :: object_type
            integer(c_int) :: api_get_object_name
        end function api_get_object_name
        !% -------------------------------------------------------------------------------
        function api_get_num_table_entries(k, table_type, num_entries)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: k
            integer(c_int), value, intent(in) :: table_type
            type(c_ptr),    value, intent(in) :: num_entries
            integer(c_int) :: api_get_num_table_entries
        end function api_get_num_table_entries
        !% -------------------------------------------------------------------------------
        function api_get_first_entry_table(k, table_type, x, y)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: k
            integer(c_int), value, intent(in) :: table_type
            type(c_ptr),    value, intent(in) :: x
            type(c_ptr),    value, intent(in) :: y
            integer(c_int) :: api_get_first_entry_table
        end function api_get_first_entry_table
        !% -------------------------------------------------------------------------------
        function api_get_next_entry_table(k, table_type, x, y)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: k
            integer(c_int), value, intent(in) :: table_type
            type(c_ptr),    value, intent(in) :: x
            type(c_ptr),    value, intent(in) :: y
            integer(c_int) :: api_get_next_entry_table
        end function api_get_next_entry_table
        !% -------------------------------------------------------------------------------
        function api_get_next_entry_tseries(k)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: k
            integer(c_int)                    :: api_get_next_entry_tseries
        end function api_get_next_entry_tseries
        !% -------------------------------------------------------------------------------
        function api_find_object(object_type, object_name)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int),           value, intent(in) :: object_type
            character(c_char), dimension(*), intent(in) :: object_name
            integer(c_int) :: api_find_object
        end function api_find_object
        !% -------------------------------------------------------------------------------
        ! --- Write Output
        !% -------------------------------------------------------------------------------
        function api_export_link_results(link_idx)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: link_idx
            integer(c_int)                    :: api_export_link_results
        end function api_export_link_results
        !% -------------------------------------------------------------------------------
        function api_write_output_line(t)
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double), intent(in) :: t
            integer(c_int)             :: api_write_output_line
        end function api_write_output_line
        !% -------------------------------------------------------------------------------
        function api_update_nodeResult(node_idx, resultType, newNodeResult)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: node_idx
            integer(c_int), value, intent(in) :: resultType
            real(c_double),        intent(in) :: newNodeResult
            integer(c_int)                    :: api_update_nodeResult
        end function api_update_nodeResult
        !% -------------------------------------------------------------------------------
        function api_update_linkResult(link_idx, resultType, newLinkResult)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value, intent(in) :: link_idx
            integer(c_int), value, intent(in) :: resultType
            real(c_double), value, intent(in) :: newLinkResult
            integer(c_int)                    :: api_update_linkResult
        end function api_update_linkResult
        !% -------------------------------------------------------------------------------
    end interface
!%
!%==========================================================================
!% Procedures
!%==========================================================================
!%
    procedure(api_initialize),             pointer :: ptr_api_initialize
    procedure(api_finalize),               pointer :: ptr_api_finalize
    procedure(api_get_node_attribute),     pointer :: ptr_api_get_node_attribute
    procedure(api_get_link_attribute),     pointer :: ptr_api_get_link_attribute
    procedure(api_get_table_attribute),    pointer :: ptr_api_get_table_attribute
    procedure(api_get_num_objects),        pointer :: ptr_api_get_num_objects
    procedure(api_get_object_name_len),    pointer :: ptr_api_get_object_name_len
    procedure(api_get_object_name),        pointer :: ptr_api_get_object_name
    procedure(api_get_num_table_entries),  pointer :: ptr_api_get_num_table_entries
    procedure(api_get_first_entry_table),  pointer :: ptr_api_get_first_entry_table
    procedure(api_get_next_entry_table),   pointer :: ptr_api_get_next_entry_table
    procedure(api_get_start_datetime),     pointer :: ptr_api_get_start_datetime
    procedure(api_get_end_datetime),       pointer :: ptr_api_get_end_datetime
    procedure(api_get_flowBC),             pointer :: ptr_api_get_flowBC
    procedure(api_get_headBC),             pointer :: ptr_api_get_headBC
    procedure(api_get_report_times),       pointer :: ptr_api_get_report_times
    procedure(api_get_next_entry_tseries), pointer :: ptr_api_get_next_entry_tseries
    procedure(api_find_object),            pointer :: ptr_api_find_object
    procedure(api_run_step),               pointer :: ptr_api_run_step
    procedure(api_export_link_results),    pointer :: ptr_api_export_link_results
    procedure(api_write_output_line),      pointer :: ptr_api_write_output_line
    procedure(api_update_nodeResult),      pointer :: ptr_api_update_nodeResult
    procedure(api_update_linkResult),      pointer :: ptr_api_update_linkResult

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
    !%-----------------------------------------------------------------------------
    !% Description:
    !%    initializes the EPA-SWMM shared library, creating input (.inp),
    !%    report (.rpt), and output (.out) files, necessary to run simulation with
    !%    EPA-SWMM. It also updates the number of objects in the SWMM model, i.e.,
    !%    number of links, nodes, and tables, and defines the start and end
    !%    simulation times.
    !%-----------------------------------------------------------------------------
        integer :: ppos, num_args, error
        character(64) :: subroutine_name = 'interface_init'
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        setting%File%inp_file = trim(setting%File%inp_file) // c_null_char
        setting%File%rpt_file = trim(setting%File%rpt_file) // c_null_char
        setting%File%out_file = trim(setting%File%out_file) // c_null_char
        c_lib%filename = trim(setting%File%library_folder) // "/libswmm5.so"

        call load_api_procedure("api_initialize")
        error = ptr_api_initialize( &
            setting%File%inp_file, &
            setting%File%rpt_file, &
            setting%File%out_file, &
            setting%Simulation%useSWMMC)
        call print_api_error(error, subroutine_name)
        api_is_initialized = .true.

        !% Get number of objects

        SWMM_N_link = get_num_objects(API_LINK)
        N_link = SWMM_N_link
        SWMM_N_node = get_num_objects(API_NODE)
        N_node = SWMM_N_node
        SWMM_N_Curve = get_num_objects(API_CURVE)
        N_curve = SWMM_N_Curve

        print *
        print *, 'BUG WARNING location ',980879,' in ',subroutine_name
        print *, '...if the SWMM code detects a parse error for the *.inp file then the ...'
        print *, '...get_num_objects function returns an error code 200 (SWMM parse error)... '
        print *, '...that is stored instead of the names of nodes and links...'
        print *, '...this has unexpected errors.'
        if ((N_link == 200) .AND. (N_node == 200)) then
             print *, SWMM_N_link, SWMM_N_node
             print *, 'ERROR (input file): Appears to be parse error in the input file...'
             print *, '...where some links/nodes are either not connected or not identified...'
             print *, '...This can happen if nodes are renamed and some of the conduit connection did not get changed...'
             print *, '...This error might have been tripped accidently if the system has exactly...'
             print *, '...200 nodes and 200 links.'
             stop
        end if
        print *

        !% Defines start and end simulation times
        !% SWMM defines start and end dates as epoch times in days
        !% we need to transform those values to durations in seconds,
        !% such that our start time is zero and our end time is
        !% the total simulation duration in seconds.

        setting%Time%StartEpoch = get_start_datetime()
        setting%Time%EndEpoch = get_end_datetime()
        setting%Time%Start = 0
        setting%Time%End = (setting%Time%EndEpoch - setting%Time%StartEpoch) * real(secsperday)

        call interface_get_report_times()

        if (setting%Debug%File%interface) then
            print *, new_line("")
            print *, "SWMM_N_link", SWMM_N_link
            print *, "SWMM_N_node", SWMM_N_node
            print *, new_line("")
            print *, "SWMM start time", setting%Time%StartEpoch
            print *, "SWMM end time", setting%Time%EndEpoch
            print *, "setting%time%start", setting%Time%Start
            print *, "setting%time%end", setting%Time%End
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
        end if
    end subroutine interface_init
    !%
    !%=============================================================================
    !%=============================================================================
    !%
    subroutine interface_finalize()
    !%-----------------------------------------------------------------------------
    !% Description:
    !%    finalizes the EPA-SWMM shared library
    !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'interface_finalize'
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        call ptr_api_finalize()

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

    end subroutine interface_finalize
!%
!%=============================================================================
!%=============================================================================
!%
    ! subroutine interface_run_step()
    ! !%-----------------------------------------------------------------------------
    ! !% Description:
    ! !%    runs steps of EPA-SWMM model. If setting%Simulation%useSWMMC was defined
    ! !%    true, steps include routing model. If false, steps are for hydrology only
    ! !%-----------------------------------------------------------------------------
    !     real(8), pointer :: timeNow
    !     character(64)    :: subroutine_name = 'interface_run_step'
    ! !%-----------------------------------------------------------------------------

    !     if (setting%Debug%File%interface)  &
    !     write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

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
    !     write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

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
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        do ii = 1, SWMM_N_link
            call load_api_procedure("api_get_object_name")
            errstat = ptr_api_get_object_name(ii-1, link%Names(ii)%str, API_LINK)

            if (errstat /= 0) then
                write(*, "(A,i2,A)") "API ERROR : ", errstat, " [" // subroutine_name // "]"
                stop
            end if
        end do

        do ii = 1, SWMM_N_node
            call load_api_procedure("api_get_object_name")
            errstat = ptr_api_get_object_name(ii-1, node%Names(ii)%str, API_NODE)

            if (errstat /= 0) then
                write(*, "(A,i2,A)") "API ERROR : ", errstat, " [" // subroutine_name // "]"
                stop
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
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
        end if

    end subroutine interface_update_linknode_names
!%
!%=============================================================================
!%=============================================================================
!%
    integer function interface_get_obj_name_len(obj_idx, obj_type) result(obj_name_len)
    !%-----------------------------------------------------------------------------
    !% Description:
    !%    Returns the lenght of the name string associated to the EPA-SWMM object.
    !%    This function is necessary to allocate the entries of the link%Names and
    !%    node%Names arraysr. The function is currently compatible with NODE and
    !%    LINK types.
    !%-----------------------------------------------------------------------------
        integer, intent(in)    :: obj_idx  ! index of the EPA-SWMM object
        integer, intent(in)    :: obj_type ! type of EPA-SWMM object (API_NODE, API_LINK)
        integer                :: error
        type(c_ptr)            :: cptr_value
        real(c_double), target :: len_value
        character(64)          :: subroutine_name = "interface_get_obj_name_len"
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
        write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        cptr_value = c_loc(len_value)
        call load_api_procedure("api_get_object_name_len")
        error = ptr_api_get_object_name_len(obj_idx-1, obj_type, cptr_value)
        call print_api_error(error, subroutine_name)
        obj_name_len = len_value

        if (setting%Debug%File%interface) then
            print *, obj_idx, obj_type, obj_name_len
        end if

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end function interface_get_obj_name_len
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_node_attribute(node_idx, attr) result(attr_value)
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
        real(8) :: attr_value
        type(c_ptr) :: cptr_value
        real(c_double), target :: node_value
        character(64) :: subroutine_name = 'interface_get_node_attribute'
    !%-----------------------------------------------------------------------------

        cptr_value = c_loc(node_value)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        if ((attr > N_api_node_attributes) .or. (attr < 1)) then
            print *, "error: unexpected node attribute value", attr
            stop
        end if

        if ((node_idx > N_node) .or. (node_idx < 1)) then
            print *, "error: unexpected node index value", node_idx
            stop
        end if

        !% Substracts 1 to every Fortran index (it becomes a C index)
        call load_api_procedure("api_get_node_attribute")
        error = ptr_api_get_node_attribute(node_idx-1, attr, cptr_value)
        call print_api_error(error, subroutine_name)

        attr_value = node_value

        !% Adds 1 to every C index extracted from EPA-SWMM (it becomes a Fortran index)
        if ((attr == api_node_extInflow_tSeries) .or. (attr == api_node_extInflow_basePat)) then
            if (node_value /= -1) attr_value = attr_value + 1
        end if

        if (setting%Debug%File%interface) &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end function interface_get_node_attribute
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_link_attribute(link_idx, attr) result(attr_value)
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
        real(8) :: attr_value

        type(c_ptr) :: cptr_value
        real(c_double), target :: link_value
        character(64) :: thisposition
        character(64) :: subroutine_name = 'interface_get_link_attribute'
    !%-----------------------------------------------------------------------------

        cptr_value = c_loc(link_value)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        if ((attr > N_api_total_link_attributes) .or. (attr < 1)) then
            print *, "error: unexpected link attribute value", attr
            stop
        end if

        if ((link_idx > SWMM_N_link) .or. (link_idx < 1)) then
            print *, "error: unexpected link index value", link_idx
            stop
        end if

        if (attr <= N_api_link_attributes) then
            ! Fortran index starts in 1, whereas in C starts in 0
            call load_api_procedure("api_get_link_attribute")
            error = ptr_api_get_link_attribute(link_idx-1, attr, cptr_value)
            thisposition = trim(subroutine_name)//'_A01'
            call print_api_error(error, thisposition)
            attr_value = link_value
        else if ( (attr > N_api_link_attributes) .and. &
                  (attr <= (N_api_link_attributes + N_api_link_type_attributes)) )then
            call load_api_procedure("api_get_link_attribute")
            error = ptr_api_get_link_attribute(link_idx-1, api_link_type, cptr_value)
            thisposition = trim(subroutine_name)//'_B02'
            call print_api_error(error, thisposition)
            attr_value = link_value
            if (link_value == API_CONDUIT) then
                if (attr == api_link_type) then
                    attr_value = lPipe
                else if (attr == api_weir_type) then
                    attr_value = nullvalueI
                else if (attr == api_orifice_type) then
                    attr_value = nullvalueI
                else if (attr == api_pump_type) then
                    attr_value = nullvalueI
                end if

            else if (link_value == API_PUMP) then

                if (attr == api_link_type) then
                    attr_value = lPump
                else if (attr == api_weir_type) then
                    attr_value = nullvalueI
                else if (attr == api_orifice_type) then
                    attr_value = nullvalueI
                else if (attr == api_pump_type) then
                    call load_api_procedure("api_get_link_attribute")
                    error = ptr_api_get_link_attribute(link_idx-1, api_pump_type, cptr_value)
                    thisposition = trim(subroutine_name)//'_C03'
                    call print_api_error(error, thisposition)
                    attr_value = link_value
                    if (link_value == API_TYPE1_PUMP) then
                        attr_value = lType1Pump
                    else if (link_value == API_TYPE2_PUMP) then
                        attr_value = lType2Pump
                    else if (link_value == API_TYPE3_PUMP) then
                        attr_value = lType3Pump
                    else if (link_value == API_TYPE4_PUMP) then
                        attr_value = lType4Pump
                    else if (link_value == API_IDEAL_PUMP) then
                        attr_value = lTypeIdealPump
                    endif
                end if

            else if (link_value == API_ORIFICE) then
                if (attr == api_link_type) then
                    attr_value = lOrifice
                else if (attr == api_weir_type) then
                    attr_value = nullvalueI
                else if (attr == api_orifice_type) then
                    call load_api_procedure("api_get_link_attribute")
                    error = ptr_api_get_link_attribute(link_idx-1, api_orifice_type, cptr_value)
                    thisposition = trim(subroutine_name)//'_D04'
                    call print_api_error(error, thisposition)
                    attr_value = link_value
                    if (link_value == API_SIDE_ORIFICE) then
                        attr_value = lSideOrifice
                    else if (link_value == API_BOTTOM_ORIFICE) then
                        attr_value = lBottomOrifice
                    endif
                else if (attr == api_pump_type) then
                    attr_value = nullvalueI
                end if

            else if (link_value == API_WEIR) then
                if (attr == api_link_type) then
                    attr_value = lWeir
                else if (attr == api_weir_type) then
                    call load_api_procedure("api_get_link_attribute")
                    error = ptr_api_get_link_attribute(link_idx-1, api_weir_type, cptr_value)
                    thisposition = trim(subroutine_name)//'_E05'
                    call print_api_error(error, thisposition)
                    attr_value = link_value
                    if (link_value == API_TRANSVERSE_WEIR) then
                        attr_value = lType1Pump
                    else if (link_value == API_SIDEFLOW_WEIR) then
                        attr_value = lType2Pump
                    else if (link_value == API_VNOTCH_WEIR) then
                        attr_value = lType3Pump
                    else if (link_value == API_TRAPEZOIDAL_WEIR) then
                        attr_value = lType4Pump
                    else if (link_value == API_ROADWAY_WEIR) then
                        attr_value = lTypeIdealPump
                    endif
                else if (attr == api_orifice_type) then
                    attr_value = nullvalueI
                else if (attr == api_pump_type) then
                    attr_value = nullvalueI
                end if
            endif

        else
            call load_api_procedure("api_get_link_attribute")
            error = ptr_api_get_link_attribute(link_idx-1, api_link_xsect_type, cptr_value)
            thisposition = trim(subroutine_name)//'_E05'
            call print_api_error(error, thisposition)
            attr_value = link_value
            if (link_value == API_RECT_CLOSED) then
                if (attr == api_link_geometry) then
                    attr_value = lRectangular_closed
                else if (attr == api_link_xsect_wMax) then
                    call load_api_procedure("api_get_link_attribute")
                    error = ptr_api_get_link_attribute(link_idx-1, api_link_xsect_wMax, cptr_value)
                    thisposition = trim(subroutine_name)//'_F06'
                    call print_api_error(error, thisposition)
                    attr_value = link_value
                else if (attr == api_link_xsect_yFull) then
                    call load_api_procedure("api_get_link_attribute")
                    error = ptr_api_get_link_attribute(link_idx-1, api_link_xsect_yFull, cptr_value)
                    thisposition = trim(subroutine_name)//'_G07'
                    call print_api_error(error, thisposition)
                    attr_value = link_value
                else
                    attr_value = nullvalueR
                end if
            else if (link_value == API_RECT_OPEN) then
                if (attr == api_link_geometry) then
                    attr_value = lRectangular
                else if (attr == api_link_xsect_wMax) then
                    call load_api_procedure("api_get_link_attribute")
                    error = ptr_api_get_link_attribute(link_idx-1, api_link_xsect_wMax, cptr_value)
                    thisposition = trim(subroutine_name)//'_H08'
                    call print_api_error(error, thisposition)
                    attr_value = link_value
                else if (attr == api_link_xsect_yFull) then
                    call load_api_procedure("api_get_link_attribute")
                    error = ptr_api_get_link_attribute(link_idx-1, api_link_xsect_yFull, cptr_value)
                    thisposition = trim(subroutine_name)//'_I09'
                    call print_api_error(error, thisposition)
                    attr_value = link_value
                else
                    attr_value = nullvalueR
                end if
            else if (link_value == API_TRAPEZOIDAL) then
                if (attr == api_link_geometry) then
                    attr_value = lTrapezoidal
                else if (attr == api_link_xsect_wMax) then
                    call load_api_procedure("api_get_link_attribute")
                    error = ptr_api_get_link_attribute(link_idx-1, api_link_xsect_yBot, cptr_value)
                    thisposition = trim(subroutine_name)//'_J10'
                    call print_api_error(error, thisposition)
                    attr_value = link_value
                else if (attr == api_link_xsect_yFull) then
                    call load_api_procedure("api_get_link_attribute")
                    error = ptr_api_get_link_attribute(link_idx-1, api_link_xsect_yFull, cptr_value)
                    thisposition = trim(subroutine_name)//'_K11'
                    call print_api_error(error, thisposition)
                    attr_value = link_value
                else
                    attr_value = nullvalueR
                end if
            else if (link_value == API_TRIANGULAR) then
                if (attr == api_link_geometry) then
                    attr_value = lTriangular
                else if (attr == api_link_xsect_wMax) then
                    call load_api_procedure("api_get_link_attribute")
                    error = ptr_api_get_link_attribute(link_idx-1, api_link_xsect_wMax, cptr_value)
                    thisposition = trim(subroutine_name)//'_M12'
                    call print_api_error(error, thisposition)
                    attr_value = link_value
                else if (attr == api_link_xsect_yFull) then
                    call load_api_procedure("api_get_link_attribute")
                    error = ptr_api_get_link_attribute(link_idx-1, api_link_xsect_yFull, cptr_value)
                    thisposition = trim(subroutine_name)//'_N13'
                    call print_api_error(error, thisposition)
                    attr_value = link_value
                else
                    attr_value = nullvalueR
                end if
            else if (link_value == API_PARABOLIC) then
                if (attr == api_link_geometry) then
                    attr_value = lParabolic
                else if (attr == api_link_xsect_wMax) then
                    call load_api_procedure("api_get_link_attribute")
                    error = ptr_api_get_link_attribute(link_idx-1, api_link_xsect_wMax, cptr_value)
                    thisposition = trim(subroutine_name)//'_O14'
                    call print_api_error(error, thisposition)
                    attr_value = link_value
                else if (attr == api_link_xsect_yFull) then
                    call load_api_procedure("api_get_link_attribute")
                    error = ptr_api_get_link_attribute(link_idx-1, api_link_xsect_yFull, cptr_value)
                    thisposition = trim(subroutine_name)//'_P15'
                    call print_api_error(error, thisposition)
                    attr_value = link_value
                else
                    attr_value = nullvalueR
                end if
            else if (link_value == API_CIRCULAR) then
                if (attr == api_link_geometry) then
                    attr_value = lCircular
                else if (attr == api_link_xsect_wMax) then
                    call load_api_procedure("api_get_link_attribute")
                    error = ptr_api_get_link_attribute(link_idx-1, api_link_xsect_wMax, cptr_value)
                    thisposition = trim(subroutine_name)//'_Q16'
                    call print_api_error(error, thisposition)
                    attr_value = link_value
                else if (attr == api_link_xsect_yFull) then
                    call load_api_procedure("api_get_link_attribute")
                    error = ptr_api_get_link_attribute(link_idx-1, api_link_xsect_yFull, cptr_value)
                    thisposition = trim(subroutine_name)//'_R17'
                    call print_api_error(error, thisposition)
                    attr_value = link_value
                else
                    attr_value = nullvalueR
                end if
            else
                attr_value = nullvalueR
            end if
        end if
        if (setting%Debug%File%interface)  then
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
            ! print *, "LINK", link_value, attr
        end if
    end function interface_get_link_attribute
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

        type(c_ptr) :: cptr_value
        real(c_double), target :: table_value
        character(64) :: thisposition
        character(64) :: subroutine_name = 'interface_get_table_attribute'
    !%-----------------------------------------------------------------------------

        cptr_value = c_loc(table_value)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        if ((attr > N_api_total_table_attributes) .or. (attr < 1)) then
            print *, "error: unexpected table attribute value", attr
            stop
        end if

        if ((table_idx > SWMM_N_Curve) .or. (table_idx < 1)) then
            print *, "error: unexpected table index value", table_idx
            stop
        end if

        !% Substracts 1 to every Fortran index (it becomes a C index)
        if (attr == api_table_type) then
            call load_api_procedure("api_get_table_attribute")
            error = ptr_api_get_table_attribute(table_idx-1, attr, cptr_value)
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
                interface_get_table_attribute = nullvalueI
            end if
        else
            call load_api_procedure("api_get_table_attribute")
            error = ptr_api_get_table_attribute(table_idx-1, attr, cptr_value)
            call print_api_error(error, subroutine_name)
            interface_get_table_attribute = table_value
        end if

        if (setting%Debug%File%interface)  then
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
            print *, "table", table_value, attr
        end if
    end function interface_get_table_attribute
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_num_table_entries(table_idx) result(num_entries)
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
        integer :: num_entries

        type(c_ptr) :: cptr_value
        integer(c_int), target :: table_entries
        character(64) :: thisposition
        character(64) :: subroutine_name = 'interface_get_num_table_entries'
    !%-----------------------------------------------------------------------------

        cptr_value = c_loc(table_entries)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        if ((table_idx > SWMM_N_Curve) .or. (table_idx < 1)) then
            print *, "error: unexpected table index value", table_idx
            stop
        end if

        !% Substracts 1 to every Fortran index (it becomes a C index)
        call load_api_procedure("api_get_num_table_entries")
        error = ptr_api_get_num_table_entries(table_idx-1, API_CURVE, cptr_value)
        call print_api_error(error, subroutine_name)
        num_entries = table_entries

        if (setting%Debug%File%interface)  then
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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

        type(c_ptr) :: cptr_x_entry, cptr_y_entry
        real(c_double), target :: x_entry, y_entry
        character(64) :: subroutine_name = 'interface_get_first_entry_table'
    !%-----------------------------------------------------------------------------

        cptr_x_entry = c_loc(x_entry)
        cptr_y_entry = c_loc(y_entry)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        if ((table_idx > SWMM_N_Curve) .or. (table_idx < 1)) then
            print *, "error: unexpected table index value", table_idx
            stop
        end if

        !% Substracts 1 to every Fortran index (it becomes a C index)
        call load_api_procedure("api_get_first_entry_table")
        success = ptr_api_get_first_entry_table(table_idx-1, API_CURVE, cptr_x_entry, cptr_y_entry)

        if (success == 0) then
            error = -1
        else
            error = 0
        end if

        call print_api_error(error, subroutine_name)
        interface_get_first_entry_table(1) = x_entry
        interface_get_first_entry_table(2) = y_entry

        if (setting%Debug%File%interface)  then
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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

        type(c_ptr) :: cptr_x_entry, cptr_y_entry
        real(c_double), target :: x_entry, y_entry
        character(64) :: subroutine_name = 'interface_get_next_entry_table'
    !%-----------------------------------------------------------------------------

        cptr_x_entry = c_loc(x_entry)
        cptr_y_entry = c_loc(y_entry)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        if ((table_idx > SWMM_N_Curve) .or. (table_idx < 1)) then
            print *, "error: unexpected table index value", table_idx
            stop
        end if

        !% Substracts 1 to every Fortran index (it becomes a C index)
        call load_api_procedure("api_get_next_entry_table")
        success = ptr_api_get_next_entry_table(table_idx-1, table_type, cptr_x_entry, cptr_y_entry)

        if (success == 0) then
            error = -1
        else
            error = 0
        end if

        call print_api_error(error, subroutine_name)
        interface_get_next_entry_table(1) = x_entry
        interface_get_next_entry_table(2) = y_entry

        if (setting%Debug%File%interface)  then
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
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
        !%-----------------------------------------------------------------------------

        resolution = nullvalueI

        if (node%YN(node_idx, nYN_has_inflow)) then ! Upstream/Lateral BC

            resolution = 0

            if (node%YN(node_idx, nYN_has_extInflow)) then

                p0 = interface_get_node_attribute(node_idx, api_node_extInflow_basePat_type)

                if (p0 == api_hourly_pattern) then
                    resolution = api_hourly
                else if (p0 == api_weekend_pattern) then
                    resolution = api_weekend
                else if (p0 == api_daily_pattern) then
                    resolution = api_daily
                else if (p0 == api_monthly_pattern) then
                    baseline = interface_get_node_attribute(node_idx, api_node_extInflow_baseline)
                    if (baseline > 0) resolution = api_monthly
                end if

            end if

            if (node%YN(node_idx, nYN_has_dwfInflow)) then

                p1 = interface_get_node_attribute(node_idx, api_node_dwfInflow_hourly_pattern)
                p2 = interface_get_node_attribute(node_idx, api_node_dwfInflow_weekend_pattern)
                p3 = interface_get_node_attribute(node_idx, api_node_dwfInflow_daily_pattern)
                p4 = interface_get_node_attribute(node_idx, api_node_dwfInflow_monthly_pattern)

                if (p1 > 0) then
                    resolution = max(api_hourly, resolution)
                else if (p2 > 0) then
                    resolution = max(api_weekend, resolution)
                else if (p3 > 0) then
                    resolution = max(api_daily, resolution)
                else if (p4 > 0) then
                    resolution = max(api_monthly, resolution)
                end if

            end if
        end if

    end function interface_get_BC_resolution
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_next_inflow_time(bc_idx, tnow) result(tnext)
        integer, intent(in) :: bc_idx
        real(8), intent(in) :: tnow
        real(8)             :: tnext, t1, t2, tnextp
        integer             :: nidx, nres, tseries, success
        character(64) :: subroutine_name

        subroutine_name = 'interface_get_next_inflow_time'

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        nidx = BC%flowI(bc_idx, bi_node_idx)
        if (.not. node%YN(nidx, nYN_has_inflow)) then
            print *, "Error, node " // node%Names(nidx)%str // " does not have an inflow"
        end if
        nres = node%I(nidx, ni_pattern_resolution)
        if (nres >= 0) then
            tnextp = util_datetime_get_next_time(tnow, nres)
            if (node%YN(nidx, nYN_has_extInflow)) then
                tseries = interface_get_node_attribute(nidx, api_node_extInflow_tSeries)
                if (tseries >= 0) then
                    success = get_next_entry_tseries(tseries)
                    tnext = interface_get_node_attribute(nidx, api_node_extInflow_tSeries_x1)
                    tnext = util_datetime_epoch_to_secs(tnext)
                    if (success == 0) then ! unsuccessful
                        tnext = interface_get_node_attribute(nidx, api_node_extInflow_tSeries_x2)
                        tnext = util_datetime_epoch_to_secs(tnext)
                        if (tnext == tnow) then
                            tnext = setting%Time%End
                            setting%BC%disableInterpolation = .true.
                        end if
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
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end function interface_get_next_inflow_time
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_next_head_time(bc_idx, tnow) result(tnext)
        integer, intent(in) :: bc_idx
        real(8), intent(in) :: tnow
        real(8)             :: tnext, tnextp
        integer             :: nidx, nres, tseries
        character(64) :: subroutine_name

        subroutine_name = 'interface_get_next_head_time'

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        nidx = BC%headI(bc_idx, bi_node_idx)
        if (BC%headI(bc_idx, bi_subcategory) == BCH_fixed) then
            tnext = setting%Time%End
        else
            print *, "Error, unsupported head boundary condition for node " // node%Names(nidx)%str
            stop
        end if

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

    end function interface_get_next_head_time
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_flowBC(bc_idx, tnow) result(bc_value)
        integer, intent(in) :: bc_idx
        real(8), intent(in) :: tnow
        integer             :: nidx
        real(8)             :: epochNow, bc_value
        character(64) :: subroutine_name

        subroutine_name = 'interface_get_flowBC'

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        nidx = BC%flowI(bc_idx, bi_node_idx)
        epochNow = util_datetime_secs_to_epoch(tnow)
        call load_api_procedure("api_get_flowBC")
        bc_value = ptr_api_get_flowBC(nidx-1, epochNow)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

    end function interface_get_flowBC
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_get_headBC(bc_idx, tnow) result(bc_value)
        integer, intent(in) :: bc_idx
        real(8), intent(in) :: tnow
        integer             :: nidx
        real(8)             :: epochNow, bc_value
        character(64) :: subroutine_name

        subroutine_name = 'interface_get_headBC'

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        nidx = BC%headI(bc_idx, bi_node_idx)
        epochNow = util_datetime_secs_to_epoch(tnow)
        call load_api_procedure("api_get_headBC")
        bc_value = ptr_api_get_headBC(nidx-1, epochNow)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

    end function interface_get_headBC
!%
!%=============================================================================
!%    Write Outputs (execute after initialization only)
!%=============================================================================
!%
    subroutine interface_export_link_results(link_idx)
        integer, intent(in) :: link_idx
        integer :: error
        character(64) :: subroutine_name = 'interface_export_link_results'

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        call load_api_procedure("api_export_link_results")
        error = ptr_api_export_link_results(link_idx-1)
        call print_api_error(error, subroutine_name)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine interface_export_link_results
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_update_nodeResult(node_idx, result_type, node_result)
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: node_idx, result_type
        real(8), intent(in) :: node_result
        integer             :: error
        character(64)       :: subroutine_name = "interface_update_nodeResult"
        !%-----------------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        call load_api_procedure("api_update_nodeResult")
        error = ptr_api_update_nodeResult(node_idx-1, result_type, node_result)
        call print_api_error(error, subroutine_name)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine interface_update_nodeResult
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_update_linkResult(link_idx, result_type, link_result)
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: link_idx, result_type
        real(8), intent(in) :: link_result
        integer             :: error
        character(64)       :: subroutine_name = "interface_update_linkResult"
        !%-----------------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        call load_api_procedure("api_update_linkResult")
        error = ptr_api_update_linkResult(link_idx-1, result_type, link_result)
        call print_api_error(error, subroutine_name)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine interface_update_linkResult
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_write_output_line(reportTime)
    !%-----------------------------------------------------------------------------
    !% Description:
    !%    Writes .out file with SWMM5+ data
    !%-----------------------------------------------------------------------------
        real(c_double),intent(in) :: reportTime ! time in seconds
        integer                   :: error
        character(64)             :: subroutine_name = "interface_write_output_line"
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        call load_api_procedure("api_write_output_line")
        error = ptr_api_write_output_line(reportTime)
        call print_api_error(error, subroutine_name)

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine interface_write_output_line
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine interface_get_report_times()
    !%-----------------------------------------------------------------------------
        integer                :: error
        real(c_double), target :: reportStart
        integer(c_int), target :: reportStep, hydroStep
        type(c_ptr)            :: cptr_reportStart, cptr_reportStep, cptr_hydroStep
        character(64)          :: subroutine_name = 'interface_get_report_times'
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        cptr_reportStart = c_loc(reportStart)
        cptr_reportStep = c_loc(reportStep)
        cptr_hydroStep = c_loc(hydroStep)

        !% --- reportStart is given in epoch datetime
        !% --- reportStep and hydroStep are given in integer seconds
        !%
        !% --- Note that if reportStart is earlier than StartEpoch from get_start_datetime()
        !% ... then SWMM will make reportStart = start datetime.
        !%
        !% --- This generates EPA-SWMM Error Code: -1 if the reportStart is after the end time.
        call load_api_procedure("api_get_report_times")
        error = ptr_api_get_report_times(cptr_reportStart, cptr_reportStep, cptr_hydroStep)
        call print_api_error(error, subroutine_name)

        reportStart = util_datetime_epoch_to_secs(reportStart)

        setting%Output%reportStartTime = reportStart
        setting%Output%reportDt = reportStep
        setting%Time%Hydrology%Dt = hydroStep

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine interface_get_report_times
!%
!%=============================================================================
!%=============================================================================
!%
    function interface_find_object(object_type, object_name) result(object_idx)
        !% Returns the index of the object, or 0 if the object couldn't be found
        character(*), intent(in) :: object_name
        integer, intent(in) :: object_type
        integer :: object_idx
        character(64) :: subroutine_name

        subroutine_name = 'interface_find_object'

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        call load_api_procedure("api_find_object")
        object_idx = ptr_api_find_object(object_type, trim(object_name)//c_null_char) + 1

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

    end function interface_find_object
!%
!%=============================================================================
!% PRIVATE
!%=============================================================================
!%
    subroutine load_api_procedure(api_procedure_name)
        character(len=*) :: api_procedure_name
        character(64) :: subroutine_name = 'load_api_procedure'

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        !% Open the shared library
        call c_lib_open(c_lib, errstat, errmsg)
        c_lib%procname = api_procedure_name
        call c_lib_load(c_lib, errstat, errmsg)

        !% Loads shared library funcitonalities
        select case (api_procedure_name)
            case ("api_initialize")
                call c_f_procpointer(c_lib%procaddr, ptr_api_initialize)
            case ("api_finalize")
                call c_f_procpointer(c_lib%procaddr, ptr_api_finalize)
            case ("api_get_node_attribute")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_node_attribute)
            case ("api_get_link_attribute")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_link_attribute)
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
            case ("api_get_report_times")
                call c_f_procpointer(c_lib%procaddr, ptr_api_get_report_times)
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
            case default
                write(*,"(A,A)") "Error, procedure " // api_procedure_name // " cannot been handled"
                stop
        end select

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine load_api_procedure

    function get_next_entry_tseries(k) result(success)
        integer, intent(in   ) :: k
        integer                :: success
        character(64)          :: subroutine_name

        subroutine_name = 'get_next_entry_tseries'

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        call load_api_procedure("api_get_next_entry_tseries")
        success = ptr_api_get_next_entry_tseries(k-1) ! Fortran to C convention

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end function get_next_entry_tseries
!%
!%=============================================================================
!%=============================================================================
!%
    function get_num_objects(obj_type)

        integer :: obj_type
        integer :: get_num_objects
        character(64) :: subroutine_name

        subroutine_name = 'get_num_objects'

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        call load_api_procedure("api_get_num_objects")
        print *, "HERE", obj_type
        get_num_objects = ptr_api_get_num_objects(obj_type)
        print *, "HERE", 2

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

    end function get_num_objects
!%
!%=============================================================================
!%=============================================================================
!%
    function get_start_datetime()
        real(8) :: get_start_datetime
        character(64) :: subroutine_name = 'get_start_datetime'

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        call load_api_procedure("api_get_start_datetime")
        get_start_datetime = ptr_api_get_start_datetime()

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end function get_start_datetime
!%
!%=============================================================================
!%=============================================================================
!%
    function get_end_datetime()
        real(8) :: get_end_datetime
        character(64) :: subroutine_name

        subroutine_name = 'get_end_datetime'

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        call load_api_procedure("api_get_end_datetime")
        get_end_datetime = ptr_api_get_end_datetime()

        if (setting%Debug%File%interface)  &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end function get_end_datetime
!%
!%=============================================================================
!%=============================================================================
!%
    subroutine print_api_error(error, subroutine_name)
        integer, intent(in) :: error
        character(64), intent(in) :: subroutine_name

        if (error /= 0) then
            write(*, "(A,i5,A)") new_line("") // "EPA-SWMM Error Code: ", error, " in "// subroutine_name
            stop
        end if
    end subroutine print_api_error
!%
!%=============================================================================
!% END MODULE interface
!%=============================================================================
!%
end module interface
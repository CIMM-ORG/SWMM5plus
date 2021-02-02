module interface

    use errors
    use iso_c_binding
    use dll_mod
    use objects

    use data_keys ! (comment if debugging)

    implicit none

    public

    ! interface to C DLL
    abstract interface

        ! --- Simulation

        function api_initialize(inp_file, report_file, out_file)
            use, intrinsic :: iso_c_binding
            implicit none
            character(c_char), dimension(*) :: inp_file
            character(c_char), dimension(*) :: report_file
            character(c_char), dimension(*) :: out_file
            type(c_ptr) :: api_initialize
        end function api_initialize

        subroutine api_finalize(api)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
        end subroutine api_finalize

        ! --- Property-extraction

        ! * After Initialization

        function api_get_node_attribute(api, k, attr, value)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
            integer(c_int), value :: k
            integer(c_int), value :: attr
            type(c_ptr), value, intent(in) :: value
            integer(c_int) :: api_get_node_attribute
        end function api_get_node_attribute

        function api_get_link_attribute(api, k, attr, value)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
            integer(c_int), value :: k
            integer(c_int), value :: attr
            type(c_ptr), value, intent(in) :: value
            integer(c_int) :: api_get_link_attribute
        end function api_get_link_attribute

        function api_get_num_objects(api, obj_type)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
            integer(c_int), value :: obj_type
            integer(c_int) :: api_get_num_objects
        end function api_get_num_objects

        function api_get_first_table_entry(k, table_type, x, y)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value :: k
            integer(c_int), value :: table_type
            type(c_ptr), value, intent(in) :: x
            type(c_ptr), value, intent(in) :: y
            integer(c_int) :: api_get_first_table_entry
        end function api_get_first_table_entry

        function api_get_next_table_entry(k, table_type, x, y)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value :: k
            integer(c_int), value :: table_type
            type(c_ptr), value, intent(in) :: x
            type(c_ptr), value, intent(in) :: y
            integer(c_int) :: api_get_next_table_entry
        end function api_get_next_table_entry

        function api_get_pattern_factors(k, factors)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value :: k
            type(c_ptr), intent(inout) :: factors
            integer(c_int) :: api_get_pattern_factors
        end function api_get_pattern_factors

        function api_get_pattern_type(k)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value :: k
            integer(c_int) :: api_get_pattern_type
        end function api_get_pattern_type
    end interface

    character(len = 1024), private :: errmsg
    integer, private :: errstat
    integer, private :: debuglevel = 0
    type(dll_type), private :: dll
    type(os_type) :: os
    type(c_ptr) :: api
    logical :: api_is_initialized = .false.

    ! Error codes - Uncomment if debugging (also defined in globals.f08)
    ! integer, parameter :: nullvalueI = -998877
    ! real, parameter :: nullvalueR = -9.98877e16

    ! Number of objects
    integer :: num_nodes
    integer :: num_links
    integer :: num_curves
    integer :: num_tseries

    ! SWMM objects
    integer, parameter :: SWMM_NODE = 2
    integer, parameter :: SWMM_LINK = 3
    integer, parameter :: SWMM_CURVES = 7
    integer, parameter :: SWMM_TSERIES = 8
    integer, parameter :: API_NODES_WITH_EXTINFLOW = 1000
    integer, parameter :: API_NODES_WITH_DWFINFLOW = 1001

    ! SWMM XSECT_TYPES
    integer, parameter :: SWMM_RECT_CLOSED = 3
    integer, parameter :: SWMM_RECT_OPEN = 4
    integer, parameter :: SWMM_TRAPEZOIDAL = 5
    integer, parameter :: SWMM_TRIANGULAR = 6
    integer, parameter :: SWMM_PARABOLIC = 7

    ! SWMM+ XSECT_TYPES - Uncomment if debugging (also defined in data_keys.f08)
    ! integer, parameter :: lchannel = 1
    ! integer, parameter :: lpipe = 2
    ! integer, parameter :: lRectangular = 1
    ! integer, parameter :: lParabolic = 2
    ! integer, parameter :: lTrapezoidal = 3
    ! integer, parameter :: lTriangular = 4

    ! api_node_attributes
    integer, parameter :: node_ID = 1
    integer, parameter :: node_type = 2
    integer, parameter :: node_invertElev = 3
    integer, parameter :: node_initDepth = 4
    integer, parameter :: node_extInflow_tSeries = 5
    integer, parameter :: node_extInflow_basePat = 6
    integer, parameter :: node_extInflow_baseline = 7
    integer, parameter :: node_extInflow_sFactor = 8
    integer, parameter :: node_has_extInflow = 9
    integer, parameter :: node_dwfInflow_monthly_pattern = 10
    integer, parameter :: node_dwfInflow_daily_pattern = 11
    integer, parameter :: node_dwfInflow_hourly_pattern = 12
    integer, parameter :: node_dwfInflow_weekly_pattern = 13
    integer, parameter :: node_dwfInflow_avgvalue = 14
    integer, parameter :: node_has_dwfInflow = 15
    integer, parameter :: node_inflow = 16
    integer, parameter :: node_volume = 17
    integer, parameter :: node_overflow = 18
    integer, parameter :: num_node_attributes = 18

    ! api_link_attributes
    integer, parameter :: link_ID = 1
    integer, parameter :: link_subIndex = 2
    integer, parameter :: link_node1 = 3
    integer, parameter :: link_node2 = 4
    integer, parameter :: link_q0 = 5
    integer, parameter :: link_flow = 6
    integer, parameter :: link_depth = 7
    integer, parameter :: link_volume = 8
    integer, parameter :: link_froude = 9
    integer, parameter :: link_setting = 10
    integer, parameter :: link_left_slope = 11
    integer, parameter :: link_right_slope = 12
    integer, parameter :: conduit_roughness = 13
    integer, parameter :: conduit_length = 14
    integer, parameter :: num_link_attributes = 14
    ! --- xsect attributes
    integer, parameter :: link_type = 15
    integer, parameter :: link_xsect_type = 16
    integer, parameter :: link_geometry = 17
    integer, parameter :: link_xsect_wMax = 18
    integer, parameter :: link_xsect_yBot = 19
    integer, parameter :: num_link_xsect_attributes = 19 - num_link_attributes
    integer, parameter :: num_total_link_attributes = num_link_attributes + num_link_xsect_attributes

    procedure(api_initialize), pointer, private :: ptr_api_initialize
    procedure(api_finalize), pointer, private :: ptr_api_finalize
    procedure(api_get_node_attribute), pointer, private :: ptr_api_get_node_attribute
    procedure(api_get_link_attribute), pointer, private :: ptr_api_get_link_attribute
    procedure(api_get_num_objects), pointer, private :: ptr_api_get_num_objects
    procedure(api_get_first_table_entry), pointer, private :: ptr_api_get_first_table_entry
    procedure(api_get_next_table_entry), pointer, private :: ptr_api_get_next_table_entry
    procedure(api_get_pattern_factors), pointer, private :: ptr_api_get_pattern_factors
    procedure(api_get_pattern_type), pointer, private :: ptr_api_get_pattern_type
contains

    ! --- Simulation

    subroutine initialize_api()

        integer :: ppos, num_args
        character(len=256) :: inp_file ! absolute path to .inp
        character(len=256) :: rpt_file ! absolute path to .rpt
        character(len=256) :: out_file ! absolute path to .out
        character(len = 256) :: cwd
        character(64) :: subroutine_name

        subroutine_name = 'initialize_api'

        if (debuglevel > 0) print *, '*** enter ', subroutine_name

        ! Initialize C API
        api = c_null_ptr
        call init_os_type(1, os)
        call init_dll(dll)

        ! Get current working directory
        call getcwd(cwd)

        ! Retrieve .inp file path from args
        num_args = command_argument_count()
        if (num_args < 1) then
            print *, "error: path to .inp file was not defined"
            stop
        end if

        call get_command_argument(1, inp_file)

        ppos = scan(trim(inp_file), '.', back = .true.)
        if (ppos > 0) then
            rpt_file = inp_file(1:ppos) // "rpt"
            out_file = inp_file(1:ppos) // "out"
        end if

        inp_file = trim(inp_file) // c_null_char
        rpt_file = trim(rpt_file) // c_null_char
        out_file = trim(out_file) // c_null_char

        dll%filename = "libswmm5.so"

        ! Initialize API
        dll%procname = "api_initialize"
        call load_dll(os, dll, errstat, errmsg)
        call print_error(errstat, 'error: loading api_initialize')
        call c_f_procpointer(dll%procaddr, ptr_api_initialize)
        api = ptr_api_initialize(inp_file, rpt_file, out_file)

        num_links = get_num_objects(SWMM_LINK)
        num_nodes = get_num_objects(SWMM_NODE)
        num_curves = get_num_objects(SWMM_CURVES)
        num_tseries = get_num_objects(SWMM_TSERIES)

        if (debuglevel > 0)  print *, '*** leave ', subroutine_name

        api_is_initialized = .true.
    end subroutine initialize_api

    subroutine finalize_api()
        character(64) :: subroutine_name

        subroutine_name = 'finalize_api'

        if (debuglevel > 0) print *, '*** enter ', subroutine_name

        dll%procname = "api_finalize"
        call load_dll(os, dll, errstat, errmsg )
        call print_error(errstat, 'error: loading api_finalize')
        call c_f_procpointer(dll%procaddr, ptr_api_finalize)
        call ptr_api_finalize(api)
        if (errstat /= 0) then
            call print_error(errstat, dll%procname)
            stop
        end if
        if (debuglevel > 0)  print *, '*** leave ', subroutine_name

    end subroutine finalize_api

    ! --- Property-extraction

    ! * After Initialization

    function get_node_attribute(node_idx, attr)

        integer :: node_idx, attr, error
        real :: get_node_attribute
        character(64) :: subroutine_name
        type(c_ptr) :: cptr_value
        real (c_double), target :: node_value

        cptr_value = c_loc(node_value)

        subroutine_name = 'get_node_attr'

        if (debuglevel > 0) print *, '*** enter ', subroutine_name

        if ((attr > num_node_attributes) .or. (attr < 1)) then
            print *, "error: unexpected node attribute value", attr
            stop
        end if

        if ((node_idx > num_nodes) .or. (node_idx < 1)) then
            print *, "error: unexpected node index value", node_idx
            stop
        end if

        dll%procname = "api_get_node_attribute"
        call load_dll(os, dll, errstat, errmsg)
        call print_error(errstat, 'error: loading api_get_node_attribute')
        call c_f_procpointer(dll%procaddr, ptr_api_get_node_attribute)
        ! Fortran index starts in 1, whereas in C starts in 0
        error = ptr_api_get_node_attribute(api, node_idx-1, attr, cptr_value)
        call print_swmm_error_code(error)

        get_node_attribute = node_value

        if (debuglevel > 0)  then
            print *, '*** leave ', subroutine_name
            ! print *, "NODE", node_value, attr
        end if
    end function get_node_attribute

    function get_link_attribute(link_idx, attr)

        integer :: link_idx, attr, error
        real :: get_link_attribute
        character(64) :: subroutine_name
        type(c_ptr) :: cptr_value
        real (c_double), target :: link_value

        cptr_value = c_loc(link_value)

        subroutine_name = 'get_link_attribute'

        if (debuglevel > 0) print *, '*** enter ', subroutine_name

        if ((attr > num_total_link_attributes) .or. (attr < 1)) then
            print *, "error: unexpected link attribute value", attr
            stop
        end if

        if ((link_idx > num_links) .or. (link_idx < 1)) then
            print *, "error: unexpected link index value", link_idx
            stop
        end if

        dll%procname = "api_get_link_attribute"
        call load_dll(os, dll, errstat, errmsg )
        call print_error(errstat, 'error: loading api_get_link_attribute')
        call c_f_procpointer(dll%procaddr, ptr_api_get_link_attribute)

        if (attr <= num_link_attributes) then
            ! Fortran index starts in 1, whereas in C starts in 0
            error = ptr_api_get_link_attribute(api, link_idx-1, attr, cptr_value)
            call print_swmm_error_code(error)
            get_link_attribute = link_value
        else
            error = ptr_api_get_link_attribute(api, link_idx-1, link_xsect_type, cptr_value)
            call print_swmm_error_code(error)
            get_link_attribute = link_value
            if (link_value == SWMM_RECT_CLOSED) then
                if (attr == link_geometry) then
                    get_link_attribute = lRectangular
                else if (attr == link_type) then
                    get_link_attribute = lpipe
                else if (attr == link_xsect_wMax) then
                    error = ptr_api_get_link_attribute(api, link_idx-1, link_xsect_wMax, cptr_value)
                    call print_swmm_error_code(error)
                    get_link_attribute = link_value
                else
                    get_link_attribute = nullvalueR
                end if
            else if (link_value == SWMM_RECT_OPEN) then
                if (attr == link_geometry) then
                    get_link_attribute = lRectangular
                else if (attr == link_type) then
                    get_link_attribute = lchannel
                else if (attr == link_xsect_wMax) then
                    error = ptr_api_get_link_attribute(api, link_idx-1, link_xsect_wMax, cptr_value)
                    call print_swmm_error_code(error)
                    get_link_attribute = link_value
                else
                    get_link_attribute = nullvalueR
                end if
            else if (link_value == SWMM_TRAPEZOIDAL) then
                if (attr == link_geometry) then
                    get_link_attribute = lTrapezoidal
                else if (attr == link_type) then
                    get_link_attribute = lchannel
                else if (attr == link_xsect_wMax) then
                    error = ptr_api_get_link_attribute(api, link_idx-1, link_xsect_yBot, cptr_value)
                    call print_swmm_error_code(error)
                    get_link_attribute = link_value
                else
                    get_link_attribute = nullvalueR
                end if
            else if (link_value == SWMM_TRIANGULAR) then
                if (attr == link_geometry) then
                    get_link_attribute = lTriangular
                else if (attr == link_type) then
                    get_link_attribute = lchannel
                else if (attr == link_xsect_wMax) then
                    error = ptr_api_get_link_attribute(api, link_idx-1, link_xsect_wMax, cptr_value)
                    call print_swmm_error_code(error)
                    get_link_attribute = link_value
                else
                    get_link_attribute = nullvalueR
                end if
            else if (link_value == SWMM_PARABOLIC) then
                if (attr == link_geometry) then
                    get_link_attribute = lParabolic
                else if (attr == link_type) then
                    get_link_attribute = lchannel
                else if (attr == link_xsect_wMax) then
                    error = ptr_api_get_link_attribute(api, link_idx-1, link_xsect_wMax, cptr_value)
                    call print_swmm_error_code(error)
                    get_link_attribute = link_value
                else
                    get_link_attribute = nullvalueR
                end if
            else
                get_link_attribute = nullvalueR
            end if
        end if
        if (debuglevel > 0)  then
            print *, '*** leave ', subroutine_name
            ! print *, "LINK", link_value, attr
        end if
    end function get_link_attribute

    function get_num_objects(obj_type)

        integer :: obj_type
        integer :: get_num_objects
        character(64) :: subroutine_name

        subroutine_name = 'get_num_objects'

        if (debuglevel > 0) print *, '*** enter ', subroutine_name

        dll%procname = "api_get_num_objects"
        call load_dll(os, dll, errstat, errmsg )
        call print_error(errstat, 'error: loading api_get_num_objects')
        call c_f_procpointer(dll%procaddr, ptr_api_get_num_objects)
        get_num_objects = ptr_api_get_num_objects(api, obj_type)
        if (debuglevel > 0)  print *, '*** leave ', subroutine_name

    end function get_num_objects

    function get_first_table_entry(k, table_type, entries)
        integer, intent(in) :: k ! table id
        integer, intent(in) :: table_type
        real, dimension(2), intent(inout) :: entries
        integer :: get_first_table_entry
        type(c_ptr) :: cptr_x, cptr_y
        real (c_double), target :: x, y

        cptr_x = c_loc(x)
        cptr_y = c_loc(y)

        dll%procname = "api_get_first_table_entry"
        call load_dll(os, dll, errstat, errmsg )
        call print_error(errstat, 'error: loading api_get_first_table_entry')
        call c_f_procpointer(dll%procaddr, ptr_api_get_first_table_entry)
        get_first_table_entry = ptr_api_get_first_table_entry(k, table_type, cptr_x, cptr_y)
        print *, x, y
    end function get_first_table_entry

    function get_next_table_entry(k, table_type, entries)
        integer, intent(in) :: k ! table id
        integer, intent(in) :: table_type
        real, dimension(2), intent(inout) :: entries
        integer :: get_next_table_entry
        type(c_ptr) :: cptr_x, cptr_y
        real (c_double), target :: x, y

        cptr_x = c_loc(x)
        cptr_y = c_loc(y)

        dll%procname = "api_get_next_table_entry"
        call load_dll(os, dll, errstat, errmsg )
        call print_error(errstat, 'error: loading api_get_next_table_entry')
        call c_f_procpointer(dll%procaddr, ptr_api_get_next_table_entry)
        get_next_table_entry = ptr_api_get_next_table_entry(k, table_type, cptr_x, cptr_y)

        entries(1) = x
        entries(2) = y
    end function get_next_table_entry

    function get_inflow_tseries(k)
        integer, intent(in) :: k
        type(tseries) :: get_inflow_tseries

        integer :: success
        real, dimension(2) :: entries

        success = get_first_table_entry(k, SWMM_TSERIES, entries)
        call tables_add_entry(get_inflow_tseries%table, entries(1), entries(2))
        do while (.true.)
            success = get_next_table_entry(k, SWMM_TSERIES, entries)
            if (success == 0) exit
            call tables_add_entry(get_inflow_tseries%table, entries(1), entries(2))
        end do
    end function get_inflow_tseries

    function get_pattern_factors(k)
        integer, intent(in) :: k
        type(pattern) :: pfactors
        type(pattern) :: get_pattern_factors
        integer :: i, count
        real(c_double), dimension(24), target :: factors
        type(c_ptr) :: cptr_factors

        cptr_factors = c_loc(factors)

        if (k .ne. -1) then
            dll%procname = "api_get_pattern_factors"
            call load_dll(os, dll, errstat, errmsg )
            call print_error(errstat, 'error: loading api_get_pattern_factors')
            call c_f_procpointer(dll%procaddr, ptr_api_get_pattern_factors)
            count = ptr_api_get_pattern_factors(k, cptr_factors)

            dll%procname = "api_get_pattern_factors"
            call load_dll(os, dll, errstat, errmsg )
            call print_error(errstat, 'error: loading api_get_pattern_factors')
            call c_f_procpointer(dll%procaddr, ptr_api_get_pattern_factors)
            get_pattern_factors%count = ptr_api_get_pattern_factors(k, cptr_factors)
            get_pattern_factors%type = get_pattern_type(k)
            get_pattern_factors%factor = factors
        end if
    end function get_pattern_factors

    function get_pattern_type(k)
        integer, intent(in) :: k
        integer :: get_pattern_type

        dll%procname = "api_get_pattern_type"
        call load_dll(os, dll, errstat, errmsg )
        call print_error(errstat, 'error: loading api_get_pattern_type')
        call c_f_procpointer(dll%procaddr, ptr_api_get_pattern_type)
        get_pattern_type = ptr_api_get_pattern_type(k)
    end function get_pattern_type

    ! --- Utils
    subroutine print_swmm_error_code(error)
        integer, intent(in) :: error
        if (error .ne. 0) then
            print *, "SWMM Error Code: " , error
            stop
        end if
    end subroutine print_swmm_error_code

end module interface

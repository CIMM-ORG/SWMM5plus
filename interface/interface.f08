module interface

    use iso_c_binding
    use c_library
    use utility
    use utility_datetime
    use define_keys
    use define_api_keys
    use define_globals
    use define_settings, only: setting

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
        function api_get_object_name_len(api, k, object_type)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
            integer(c_int), value :: k
            integer(c_int), value :: object_type
            integer(c_int) :: api_get_object_name_len
        end function api_get_object_name_len

        function api_get_object_name(api, k, object_name, object_type)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
            integer(c_int), value :: k
            character(c_char), dimension(*) :: object_name
            integer(c_int), value :: object_type
            integer(c_int) :: api_get_object_name
        end function api_get_object_name

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

        function api_get_pattern_count(k)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value :: k
            integer(c_int) :: api_get_pattern_count
        end function api_get_pattern_count

        function api_get_pattern_factor(k, j)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value :: k, j
            real(c_double) :: api_get_pattern_factor
        end function api_get_pattern_factor

        function api_get_pattern_type(k)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value :: k
            integer(c_int) :: api_get_pattern_type
        end function api_get_pattern_type

        function api_get_start_datetime()
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double) :: api_get_start_datetime
        end function api_get_start_datetime

        function api_get_end_datetime()
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double) :: api_get_end_datetime
        end function api_get_end_datetime
    end interface

    type(c_lib_type) :: c_lib
    type(c_ptr) :: api
    character(len = 1024), private :: errmsg
    integer, private :: errstat
    logical :: api_is_initialized = .false.

    ! Time constants
    real(8) :: swmm_start_time ! in days
    real(8) :: swmm_end_time ! in days

    procedure(api_initialize), pointer, private :: ptr_api_initialize
    procedure(api_finalize), pointer, private :: ptr_api_finalize
    procedure(api_get_node_attribute), pointer, private :: ptr_api_get_node_attribute
    procedure(api_get_link_attribute), pointer, private :: ptr_api_get_link_attribute
    procedure(api_get_num_objects), pointer, private :: ptr_api_get_num_objects
    procedure(api_get_first_table_entry), pointer, private :: ptr_api_get_first_table_entry
    procedure(api_get_next_table_entry), pointer, private :: ptr_api_get_next_table_entry
    procedure(api_get_pattern_count), pointer, private :: ptr_api_get_pattern_count
    procedure(api_get_pattern_factor), pointer, private :: ptr_api_get_pattern_factor
    procedure(api_get_pattern_type), pointer, private :: ptr_api_get_pattern_type
    procedure(api_get_start_datetime), pointer, private :: ptr_api_get_start_datetime
    procedure(api_get_end_datetime), pointer, private :: ptr_api_get_end_datetime
    procedure(api_get_object_name_len), pointer, private :: ptr_api_get_object_name_len
    procedure(api_get_object_name), pointer, private :: ptr_api_get_object_name

contains

    ! --- Simulation

    subroutine interface_init()

        integer :: ppos, num_args
        character(64) :: subroutine_name

        subroutine_name = 'interface_init'

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        ! Initialize C API
        api = c_null_ptr

        if (setting%Paths%inp == "") then
            print *, "ERROR: it is necessary to define the path to the .inp file"
            stop
        end if
        ppos = scan(trim(setting%Paths%inp), '.', back = .true.)
        if (ppos > 0) then
            setting%Paths%rpt = setting%Paths%inp(1:ppos) // "rpt"
            setting%Paths%out = setting%Paths%inp(1:ppos) // "out"
        end if

        setting%Paths%inp = trim(setting%Paths%inp) // c_null_char
        setting%Paths%rpt = trim(setting%Paths%rpt) // c_null_char
        setting%Paths%out = trim(setting%Paths%out) // c_null_char

        c_lib%filename = trim(setting%Paths%project) // "/libswmm5.so"

        ! Initialize API
        c_lib%procname = "api_initialize"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_initialize)
        api = ptr_api_initialize(setting%Paths%inp, setting%Paths%rpt, setting%Paths%out)

        N_link = interface_get_num_objects(API_LINK)
        N_node = interface_get_num_objects(API_NODE)
        N_curve = interface_get_num_objects(API_CURVES)
        N_tseries = interface_get_num_objects(API_TSERIES)
        N_pattern = interface_get_num_objects(API_TIMEPATTERN)

        api_is_initialized = .true.

        swmm_start_time = interface_get_start_datetime()
        swmm_end_time = interface_get_end_datetime()

        setting%time%starttime = 0
        setting%time%endtime = (swmm_end_time - swmm_start_time) * real(secsperday)

        print *, new_line("")
        if (setting%Debug%File%interface) then
            print *, "N_link", N_link
            print *, "N_node", N_node
            print *, "N_curve", N_curve
            print *, "N_tseries", N_tseries
            print *, "N_pattern", N_pattern
            print *, '*** leave ', subroutine_name
        end if
    end subroutine interface_init

    subroutine interface_finalize()
        character(64) :: subroutine_name = 'interface_finalize'

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        c_lib%procname = "api_finalize"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_finalize)
        call ptr_api_finalize(api)
        if (setting%Debug%File%interface)  print *, '*** leave ', subroutine_name

    end subroutine interface_finalize

    ! --- Property-extraction

    ! * After Initialization

    function interface_get_start_datetime()
        real(8) :: interface_get_start_datetime
        character(64) :: subroutine_name = 'interface_get_start_datetime'

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        c_lib%procname = "api_get_start_datetime"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_start_datetime)
        interface_get_start_datetime = ptr_api_get_start_datetime()
        if (setting%Debug%File%interface)  print *, '*** leave ', subroutine_name
    end function interface_get_start_datetime

    function interface_get_end_datetime()
        real(8) :: interface_get_end_datetime
        character(64) :: subroutine_name

        subroutine_name = 'interface_get_end_datetime'

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        c_lib%procname = "api_get_end_datetime"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_end_datetime)
        interface_get_end_datetime = ptr_api_get_end_datetime()
        if (setting%Debug%File%interface)  print *, '*** leave ', subroutine_name
    end function interface_get_end_datetime

    subroutine interface_update_linknode_names()
        integer :: ii
        character(64) :: subroutine_name = "interface_update_linknode_names"

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        c_lib%procname = "api_get_object_name"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_object_name)

        do ii = 1, N_node
            call ptr_api_get_object_name(api, ii-1, nodeName(ii)%str, API_NODE)
        end do

        do ii = 1, N_link
            call ptr_api_get_object_name(api, ii-1, linkName(ii)%str, API_LINK)
        end do

        if (setting%Debug%File%interface) then
            print *, new_line("")
            print *, "List of Links"
            do ii = 1, N_link
                print *, "- ", linkName(ii)%str
            end do
            print *, new_line("")
            print *, "List of Nodes"
            do ii = 1, N_node
                print *, "- ", nodeName(ii)%str
            end do
            print *, new_line("")
            print *, '*** leave ', subroutine_name
        end if

    end subroutine interface_update_linknode_names

    function interface_get_obj_name_len(obj_idx, obj_type) result(obj_name_len)
        integer, intent(in) :: obj_idx
        integer, intent(in) :: obj_type
        integer :: obj_name_len
        character(64) :: subroutine_name = "interface_get_obj_name_len"

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        c_lib%procname = "api_get_object_name_len"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_object_name_len)
        obj_name_len = ptr_api_get_object_name_len(api, obj_idx-1, obj_type)

        if (setting%Debug%File%interface) then
            print *, obj_idx, obj_type, obj_name_len
        end if

        if (setting%Debug%File%interface)  print *, '*** leave ', subroutine_name
    end function interface_get_obj_name_len

    function interface_get_node_attribute(node_idx, attr)

        integer :: node_idx, attr, error
        real(8) :: interface_get_node_attribute
        type(c_ptr) :: cptr_value
        real(c_double), target :: node_value
        character(64) :: subroutine_name = 'interface_get_node_attr'

        cptr_value = c_loc(node_value)


        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        if ((attr > N_api_node_attributes) .or. (attr < 1)) then
            print *, "error: unexpected node attribute value", attr
            stop
        end if

        if ((node_idx > N_node) .or. (node_idx < 1)) then
            print *, "error: unexpected node index value", node_idx
            stop
        end if

        c_lib%procname = "api_get_node_attribute"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_node_attribute)
        ! Fortran index starts in 1, whereas in C starts in 0
        error = ptr_api_get_node_attribute(api, node_idx-1, attr, cptr_value)
        call interface_print_error(error)

        interface_get_node_attribute = node_value

        ! Fortran index correction
        if ((attr == api_node_extInflow_tSeries) .or. (attr == api_node_extInflow_basePat)) then
            if (node_value /= -1) interface_get_node_attribute = interface_get_node_attribute + 1
        end if

        if (setting%Debug%File%interface)  then
            print *, '*** leave ', subroutine_name
            ! print *, "NODE", node_value, attr
        end if
    end function interface_get_node_attribute

    function interface_get_link_attribute(link_idx, attr)

        integer :: link_idx, attr, error
        real(8) :: interface_get_link_attribute
        character(64) :: subroutine_name
        type(c_ptr) :: cptr_value
        real(c_double), target :: link_value

        cptr_value = c_loc(link_value)

        subroutine_name = 'interface_get_link_attribute'

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        if ((attr > N_api_total_link_attributes) .or. (attr < 1)) then
            print *, "error: unexpected link attribute value", attr
            stop
        end if

        if ((link_idx > N_link) .or. (link_idx < 1)) then
            print *, "error: unexpected link index value", link_idx
            stop
        end if

        c_lib%procname = "api_get_link_attribute"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_link_attribute)

        if (attr <= N_api_link_attributes) then
            ! Fortran index starts in 1, whereas in C starts in 0
            error = ptr_api_get_link_attribute(api, link_idx-1, attr, cptr_value)
            call interface_print_error(error)
            interface_get_link_attribute = link_value
        else
            error = ptr_api_get_link_attribute(api, link_idx-1, api_link_xsect_type, cptr_value)
            call interface_print_error(error)
            interface_get_link_attribute = link_value
            if (link_value == API_RECT_CLOSED) then
                if (attr == api_link_geometry) then
                    interface_get_link_attribute = lRectangular
                else if (attr == api_link_type) then
                    interface_get_link_attribute = lpipe
                else if (attr == api_link_xsect_wMax) then
                    error = ptr_api_get_link_attribute(api, link_idx-1, api_link_xsect_wMax, cptr_value)
                    call interface_print_error(error)
                    interface_get_link_attribute = link_value
                else
                    interface_get_link_attribute = nullvalueR
                end if
            else if (link_value == API_RECT_OPEN) then
                if (attr == api_link_geometry) then
                    interface_get_link_attribute = lRectangular
                else if (attr == api_link_type) then
                    interface_get_link_attribute = lchannel
                else if (attr == api_link_xsect_wMax) then
                    error = ptr_api_get_link_attribute(api, link_idx-1, api_link_xsect_wMax, cptr_value)
                    call interface_print_error(error)
                    interface_get_link_attribute = link_value
                else
                    interface_get_link_attribute = nullvalueR
                end if
            else if (link_value == API_TRAPEZOIDAL) then
                if (attr == api_link_geometry) then
                    interface_get_link_attribute = lTrapezoidal
                else if (attr == api_link_type) then
                    interface_get_link_attribute = lchannel
                else if (attr == api_link_xsect_wMax) then
                    error = ptr_api_get_link_attribute(api, link_idx-1, api_link_xsect_yBot, cptr_value)
                    call interface_print_error(error)
                    interface_get_link_attribute = link_value
                else
                    interface_get_link_attribute = nullvalueR
                end if
            else if (link_value == API_TRIANGULAR) then
                if (attr == api_link_geometry) then
                    interface_get_link_attribute = lTriangular
                else if (attr == api_link_type) then
                    interface_get_link_attribute = lchannel
                else if (attr == api_link_xsect_wMax) then
                    error = ptr_api_get_link_attribute(api, link_idx-1, api_link_xsect_wMax, cptr_value)
                    call interface_print_error(error)
                    interface_get_link_attribute = link_value
                else
                    interface_get_link_attribute = nullvalueR
                end if
            else if (link_value == API_PARABOLIC) then
                if (attr == api_link_geometry) then
                    interface_get_link_attribute = lParabolic
                else if (attr == api_link_type) then
                    interface_get_link_attribute = lchannel
                else if (attr == api_link_xsect_wMax) then
                    error = ptr_api_get_link_attribute(api, link_idx-1, api_link_xsect_wMax, cptr_value)
                    call interface_print_error(error)
                    interface_get_link_attribute = link_value
                else
                    interface_get_link_attribute = nullvalueR
                end if
            else
                interface_get_link_attribute = nullvalueR
            end if
        end if
        if (setting%Debug%File%interface)  then
            print *, '*** leave ', subroutine_name
            ! print *, "LINK", link_value, attr
        end if
    end function interface_get_link_attribute

    function interface_get_num_objects(obj_type)

        integer :: obj_type
        integer :: interface_get_num_objects
        character(64) :: subroutine_name

        subroutine_name = 'interface_get_num_objects'

        if (setting%Debug%File%interface)  print *, '*** enter ', subroutine_name

        c_lib%procname = "api_get_num_objects"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_num_objects)
        interface_get_num_objects = ptr_api_get_num_objects(api, obj_type)
        if (setting%Debug%File%interface)  print *, '*** leave ', subroutine_name

    end function interface_get_num_objects

    function interface_get_first_table_entry(k, table_type, entries)
        integer, intent(in) :: k ! table id
        integer, intent(in) :: table_type
        real(8), dimension(2), intent(inout) :: entries
        integer :: interface_get_first_table_entry
        type(c_ptr) :: cptr_x, cptr_y
        real(c_double), target :: x, y

        cptr_x = c_loc(x)
        cptr_y = c_loc(y)

        c_lib%procname = "api_get_first_table_entry"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_first_table_entry)
        interface_get_first_table_entry = ptr_api_get_first_table_entry(k-1, table_type, cptr_x, cptr_y) ! index starts at 0 in C

        entries(1) = x
        entries(2) = y
    end function interface_get_first_table_entry

    function interface_get_next_table_entry(k, table_type, entries)
        integer, intent(in) :: k ! table id
        integer, intent(in) :: table_type
        real(8), dimension(2), intent(inout) :: entries
        integer :: interface_get_next_table_entry
        type(c_ptr) :: cptr_x, cptr_y
        real(c_double), target :: x, y

        cptr_x = c_loc(x)
        cptr_y = c_loc(y)

        c_lib%procname = "api_get_next_table_entry"
        call c_lib_load(c_lib, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        end if
        call c_f_procpointer(c_lib%procaddr, ptr_api_get_next_table_entry)
        interface_get_next_table_entry = ptr_api_get_next_table_entry(k-1, table_type, cptr_x, cptr_y) ! index starts at 0 in C

        entries(1) = x
        entries(2) = y
    end function interface_get_next_table_entry

    function interface_get_pattern(k)
        integer, intent(in) :: k
        type(pattern) :: pfactors
        type(pattern) :: interface_get_pattern
        integer :: i, count

        if (k /= -1) then
            c_lib%procname = "api_get_pattern_count"
            call c_lib_load(c_lib, errstat, errmsg)
            if (errstat /= 0) then
                print *, "ERROR: " // trim(errmsg)
                stop
            end if
            call c_f_procpointer(c_lib%procaddr, ptr_api_get_pattern_count)
            interface_get_pattern%count = ptr_api_get_pattern_count(k-1)

            c_lib%procname = "api_get_pattern_factor"
            call c_lib_load(c_lib, errstat, errmsg)
            if (errstat /= 0) then
                print *, "ERROR: " // trim(errmsg)
                stop
            end if
            call c_f_procpointer(c_lib%procaddr, ptr_api_get_pattern_factor)
            do i = 1, 24
                interface_get_pattern%factor(i) = ptr_api_get_pattern_factor(k-1, i-1) ! index starts at 0 in C
            end do

            c_lib%procname = "api_get_pattern_type"
            call c_lib_load(c_lib, errstat, errmsg)
            if (errstat /= 0) then
                print *, "ERROR: " // trim(errmsg)
                stop
            end if
            call c_f_procpointer(c_lib%procaddr, ptr_api_get_pattern_type) ! index starts at 0 in C
            interface_get_pattern%ptype = ptr_api_get_pattern_type(k-1)
        end if
    end function interface_get_pattern

    subroutine interface_print_error(error)
        integer, intent(in) :: error
        if (error /= 0) then
            print *, "SWMM Error Code: " , error
            stop
        end if
    end subroutine interface_print_error
end module interface
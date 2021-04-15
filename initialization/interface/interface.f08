module interface

    use iso_c_binding
    use dll
    use setting_definition, only: setting
    use tables
    use datetime
    use data_keys
    use globals

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

        subroutine api_print_object_name(k, object_type)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value :: k
            integer(c_int), value :: object_type
        end subroutine api_print_object_name
    end interface

    character(len = 1024), private :: errmsg
    integer, private :: errstat
    integer, private :: debuglevel = 0
    type(dll_type), private :: dll
    type(c_ptr) :: api
    logical :: api_is_initialized = .false.

    ! Time constants
    real(8) :: swmm_start_time ! in days
    real(8) :: swmm_end_time ! in days

    integer, parameter :: num_node_attributes = node_overflow

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
    procedure(api_print_object_name), pointer, private :: ptr_api_print_object_name

contains

    ! --- Simulation

    subroutine initialize_api()

        integer :: ppos, num_args
        character(64) :: subroutine_name

        subroutine_name = 'initialize_api'

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** enter ', subroutine_name

        ! Initialize C API
        api = c_null_ptr

        if (setting%Paths%inp == "") then
            print *, "ERROR: it is necessary to define the path to the .inp file"
            stop
        endif
        ppos = scan(trim(setting%Paths%inp), '.', back = .true.)
        if (ppos > 0) then
            setting%Paths%rpt = setting%Paths%inp(1:ppos) // "rpt"
            setting%Paths%out = setting%Paths%inp(1:ppos) // "out"
        end if

        setting%Paths%inp = trim(setting%Paths%inp) // c_null_char
        setting%Paths%rpt = trim(setting%Paths%rpt) // c_null_char
        setting%Paths%out = trim(setting%Paths%out) // c_null_char

        dll%filename = trim(setting%Paths%project) // "/libswmm5.so"

        ! Initialize API
        dll%procname = "api_initialize"
        call load_dll(dll, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        endif
        call c_f_procpointer(dll%procaddr, ptr_api_initialize)
        api = ptr_api_initialize(setting%Paths%inp, setting%Paths%rpt, setting%Paths%out)

        N_link = get_num_objects(SWMM_LINK)
        N_node = get_num_objects(SWMM_NODE)
        N_curve = get_num_objects(SWMM_CURVES)
        N_tseries = get_num_objects(SWMM_TSERIES)
        N_pattern = get_num_objects(SWMM_TIMEPATTERN)

        api_is_initialized = .true.

        swmm_start_time = get_start_datetime()
        swmm_end_time = get_end_datetime()

        setting%time%starttime = 0
        setting%time%endtime = (swmm_end_time - swmm_start_time) * real(secsperday)

        if (N_tseries > 0) call load_all_tseries()
        if (N_pattern > 0) call load_all_patterns()

        print *, new_line("")
        if ((debuglevel > 0) .or. (debuglevelall > 0)) then
            print *, "N_link", N_link
            print *, "N_node", N_node
            print *, "N_curve", N_curve
            print *, "N_tseries", N_tseries
            print *, "N_pattern", N_pattern
            print *, '*** leave ', subroutine_name
        endif
    end subroutine initialize_api

    subroutine finalize_api()
        character(64) :: subroutine_name = 'finalize_api'

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** enter ', subroutine_name

        dll%procname = "api_finalize"
        call load_dll(dll, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        endif
        call c_f_procpointer(dll%procaddr, ptr_api_finalize)
        call ptr_api_finalize(api)
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name

    end subroutine finalize_api

    subroutine free_api()
        integer :: i

        if (allocated(all_tseries)) then
            do i = 1, N_tseries
                call free_table(all_tseries(i))
            enddo
            deallocate(all_tseries)
        endif

        if (allocated(all_patterns)) then
            deallocate(all_patterns)
        endif
    end subroutine free_api

    ! --- Property-extraction

    ! * After Initialization

    function get_start_datetime()
        real(8) :: get_start_datetime
        character(64) :: subroutine_name = 'get_start_datetime'

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** enter ', subroutine_name

        dll%procname = "api_get_start_datetime"
        call load_dll(dll, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        endif
        call c_f_procpointer(dll%procaddr, ptr_api_get_start_datetime)
        get_start_datetime = ptr_api_get_start_datetime()
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name
    end function get_start_datetime

    function get_end_datetime()
        real(8) :: get_end_datetime
        character(64) :: subroutine_name

        subroutine_name = 'get_end_datetime'

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** enter ', subroutine_name

        dll%procname = "api_get_end_datetime"
        call load_dll(dll, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        endif
        call c_f_procpointer(dll%procaddr, ptr_api_get_end_datetime)
        get_end_datetime = ptr_api_get_end_datetime()
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name
    end function get_end_datetime

    subroutine load_all_tseries()
        integer :: i
        integer :: success
        real(8), dimension(2) :: entries
        character(64) :: subroutine_name = 'load_all_tseries'

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** enter ', subroutine_name

        if (N_tseries == 0) return

        allocate(all_tseries(N_tseries))
        do i = 1, N_tseries
            all_tseries(i) = new_real_table(SWMM_TSERIES, 2)
            success = get_first_table_entry(i, SWMM_TSERIES, entries)
            if (success == 0) then
                print *, new_line("")
                print *, "ERROR: API can't get Tseries"
                stop
            endif
            call tables_add_entry(all_tseries(i), entries)
            do while (.true.)
                success = get_next_table_entry(i, SWMM_TSERIES, entries)
                if (success == 0) exit
                if (entries(1) < swmm_start_time) cycle
                if (entries(1) > swmm_end_time) exit
                call tables_add_entry(all_tseries(i), entries)
            end do
        end do
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name
    end subroutine load_all_tseries

    subroutine load_all_patterns()
        integer :: i = 1
        integer :: success
        real(8), dimension(2) :: entries

        if (N_pattern == 0) return

        allocate(all_patterns(N_pattern))
        do i = 1, N_pattern
            all_patterns(i) = get_pattern(i)
        end do
    end subroutine load_all_patterns

    function get_node_attribute(node_idx, attr)

        integer :: node_idx, attr, error
        real(8) :: get_node_attribute
        type(c_ptr) :: cptr_value
        real(c_double), target :: node_value
        character(64) :: subroutine_name = 'get_node_attr'

        cptr_value = c_loc(node_value)


        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** enter ', subroutine_name

        if ((attr > num_node_attributes) .or. (attr < 1)) then
            print *, "error: unexpected node attribute value", attr
            stop
        end if

        if ((node_idx > N_node) .or. (node_idx < 1)) then
            print *, "error: unexpected node index value", node_idx
            stop
        end if

        dll%procname = "api_get_node_attribute"
        call load_dll(dll, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        endif
        call c_f_procpointer(dll%procaddr, ptr_api_get_node_attribute)
        ! Fortran index starts in 1, whereas in C starts in 0
        error = ptr_api_get_node_attribute(api, node_idx-1, attr, cptr_value)
        call print_swmm_error_code(error)

        get_node_attribute = node_value

        ! Fortran index correction
        if ((attr == node_extInflow_tSeries) .or. (attr == node_extInflow_basePat)) then
            if (node_value /= -1) get_node_attribute = get_node_attribute + 1
        endif

        if (debuglevel > 0)  then
            print *, '*** leave ', subroutine_name
            ! print *, "NODE", node_value, attr
        end if
    end function get_node_attribute

    function get_link_attribute(link_idx, attr)

        integer :: link_idx, attr, error
        real(8) :: get_link_attribute
        character(64) :: subroutine_name
        type(c_ptr) :: cptr_value
        real(c_double), target :: link_value

        cptr_value = c_loc(link_value)

        subroutine_name = 'get_link_attribute'

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** enter ', subroutine_name

        if ((attr > num_total_link_attributes) .or. (attr < 1)) then
            print *, "error: unexpected link attribute value", attr
            stop
        end if

        if ((link_idx > N_link) .or. (link_idx < 1)) then
            print *, "error: unexpected link index value", link_idx
            stop
        end if

        dll%procname = "api_get_link_attribute"
        call load_dll(dll, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        endif
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

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** enter ', subroutine_name

        dll%procname = "api_get_num_objects"
        call load_dll(dll, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        endif
        call c_f_procpointer(dll%procaddr, ptr_api_get_num_objects)
        get_num_objects = ptr_api_get_num_objects(api, obj_type)
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name

    end function get_num_objects

    function get_first_table_entry(k, table_type, entries)
        integer, intent(in) :: k ! table id
        integer, intent(in) :: table_type
        real(8), dimension(2), intent(inout) :: entries
        integer :: get_first_table_entry
        type(c_ptr) :: cptr_x, cptr_y
        real(c_double), target :: x, y

        cptr_x = c_loc(x)
        cptr_y = c_loc(y)

        dll%procname = "api_get_first_table_entry"
        call load_dll(dll, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        endif
        call c_f_procpointer(dll%procaddr, ptr_api_get_first_table_entry)
        get_first_table_entry = ptr_api_get_first_table_entry(k-1, table_type, cptr_x, cptr_y) ! index starts at 0 in C

        entries(1) = x
        entries(2) = y
    end function get_first_table_entry

    function get_next_table_entry(k, table_type, entries)
        integer, intent(in) :: k ! table id
        integer, intent(in) :: table_type
        real(8), dimension(2), intent(inout) :: entries
        integer :: get_next_table_entry
        type(c_ptr) :: cptr_x, cptr_y
        real(c_double), target :: x, y

        cptr_x = c_loc(x)
        cptr_y = c_loc(y)

        dll%procname = "api_get_next_table_entry"
        call load_dll(dll, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        endif
        call c_f_procpointer(dll%procaddr, ptr_api_get_next_table_entry)
        get_next_table_entry = ptr_api_get_next_table_entry(k-1, table_type, cptr_x, cptr_y) ! index starts at 0 in C

        entries(1) = x
        entries(2) = y
    end function get_next_table_entry

    function get_pattern(k)
        integer, intent(in) :: k
        type(pattern) :: pfactors
        type(pattern) :: get_pattern
        integer :: i, count

        if (k .ne. -1) then
            dll%procname = "api_get_pattern_count"
            call load_dll(dll, errstat, errmsg)
            if (errstat /= 0) then
                print *, "ERROR: " // trim(errmsg)
                stop
            endif
            call c_f_procpointer(dll%procaddr, ptr_api_get_pattern_count)
            get_pattern%count = ptr_api_get_pattern_count(k-1)

            dll%procname = "api_get_pattern_factor"
            call load_dll(dll, errstat, errmsg)
            if (errstat /= 0) then
                print *, "ERROR: " // trim(errmsg)
                stop
            endif
            call c_f_procpointer(dll%procaddr, ptr_api_get_pattern_factor)
            do i = 1, 24
                get_pattern%factor(i) = ptr_api_get_pattern_factor(k-1, i-1) ! index starts at 0 in C
            end do

            dll%procname = "api_get_pattern_type"
            call load_dll(dll, errstat, errmsg)
            if (errstat /= 0) then
                print *, "ERROR: " // trim(errmsg)
                stop
            endif
            call c_f_procpointer(dll%procaddr, ptr_api_get_pattern_type) ! index starts at 0 in C
            get_pattern%ptype = ptr_api_get_pattern_type(k-1)
        end if
    end function get_pattern

    ! --- Utils

    subroutine print_object_name(k, object_type)
        integer, intent(in) :: k
        integer, intent(in) :: object_type

        dll%procname = "api_print_object_name"
        call load_dll(dll, errstat, errmsg)
        if (errstat /= 0) then
            print *, "ERROR: " // trim(errmsg)
            stop
        endif
        call c_f_procpointer(dll%procaddr, ptr_api_print_object_name) ! index starts at 0 in C
        call ptr_api_print_object_name(k-1, object_type)
    end subroutine print_object_name

    subroutine print_swmm_error_code(error)
        integer, intent(in) :: error
        if (error .ne. 0) then
            print *, "SWMM Error Code: " , error
            stop
        end if
    end subroutine print_swmm_error_code

end module interface
module interface
    use iso_c_binding
    use dll_mod
    use data_keys
    use globals

    implicit none

    public

    ! interface to C DLL
    abstract interface
        function api_initialize(f1, f2, f3)
            use, intrinsic :: iso_c_binding
            implicit none
            character(c_char), dimension(*) :: f1, f2, f3
            type(c_ptr) :: api_initialize
        end function api_initialize

        subroutine api_finalize(api)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
        end subroutine api_finalize

        function api_get_node_attribute(api, k, attr)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
            integer(c_int), value :: k
            integer(c_int), value :: attr
            real(c_double) :: api_get_node_attribute
        end function api_get_node_attribute

        function api_get_link_attribute(api, k, attr)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
            integer(c_int), value :: k
            integer(c_int), value :: attr
            real(c_double) :: api_get_link_attribute
        end function api_get_link_attribute

        function api_num_links(api)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
            integer(c_int) :: api_num_links
        end function api_num_links

        function api_num_nodes(api)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
            integer(c_int) :: api_num_nodes
        end function api_num_nodes

    end interface

    character(len = 1024), private :: errmsg
    integer, private :: errstat
    integer, private :: debuglevel = 0
    type(dll_type), private :: dll
    type(os_type) :: os
    type(c_ptr) :: api

    ! api_node_attributes
    integer, parameter :: node_ID = 1
    integer, parameter :: node_type = 2
    integer, parameter :: node_invertElev = 3
    integer, parameter :: node_initDepth = 4
    integer, parameter :: node_extInflow_tSeries = 5
    integer, parameter :: node_extInflow_basePat = 6
    integer, parameter :: node_extInflow_baseline = 7
    integer, parameter :: node_depth = 8
    integer, parameter :: node_inflow = 9
    integer, parameter :: node_volume = 10
    integer, parameter :: node_overflow = 11
    integer, parameter :: num_node_attributes = 11

    ! api_link_attributes
    integer, parameter :: link_ID = 1
    integer, parameter :: link_subIndex = 2
    integer, parameter :: link_type = 3
    integer, parameter :: link_node1 = 4
    integer, parameter :: link_node2 = 5
    integer, parameter :: link_xsect_type = 6
    integer, parameter :: link_xsect_wMax = 7
    integer, parameter :: link_xsect_yBot = 8
    integer, parameter :: link_q0 = 9
    integer, parameter :: link_geometry = 10
    integer, parameter :: conduit_roughness = 11
    integer, parameter :: conduit_length = 12
    integer, parameter :: link_flow = 13
    integer, parameter :: link_depth = 14
    integer, parameter :: link_volume = 15
    integer, parameter :: link_froude = 16
    integer, parameter :: link_setting = 17
    integer, parameter :: num_link_attributes = 17

    procedure(api_initialize), pointer, private :: ptr_api_initialize
    procedure(api_finalize), pointer, private :: ptr_api_finalize
    procedure(api_num_links), pointer, private :: ptr_api_num_links
    procedure(api_num_nodes), pointer, private :: ptr_api_num_nodes
    procedure(api_get_node_attribute), pointer, private :: ptr_api_get_node_attribute
    procedure(api_get_link_attribute), pointer, private :: ptr_api_get_link_attribute

contains
    subroutine initialize_api()

        integer :: ppos, num_args
        character(len=256) :: inp_file ! absolute path to .inp
        character(len=256) :: rpt_file ! absolute path to .rpt
        character(len=256) :: out_file ! absolute path to .out
        character(len = 256) :: cwd
        character(64) :: subroutine_name

        subroutine_name = 'initialize_api'

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ', subroutine_name

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

        ! TODO - change relative path
        dll%filename = trim(cwd) // os%pathsep // "libswmm5.so"

        ! Initialize API
        dll%procname = "api_initialize"
        call load_dll(os, dll, errstat, errmsg)
        call print_error(errstat, 'error: loading api_initialize')
        call c_f_procpointer(dll%procaddr, ptr_api_initialize)
        api = ptr_api_initialize(inp_file, rpt_file, out_file)
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name

    end subroutine initialize_api

    subroutine finalize_api()
        character(64) :: subroutine_name

        subroutine_name = 'finalize_api'

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ', subroutine_name

        dll%procname = "api_finalize"
        call load_dll(os, dll, errstat, errmsg )
        call print_error(errstat, 'error: loading api_finalize')
        call c_f_procpointer(dll%procaddr, ptr_api_finalize)
        call ptr_api_finalize(api)
        if (errstat /= 0) then
            call print_error(errstat, dll%procname)
            stop
        end if
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name

    end subroutine finalize_api

    function get_num_links()

        integer :: get_num_links
        character(64) :: subroutine_name

        subroutine_name = 'get_num_links'

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ', subroutine_name

        dll%procname = "api_num_links"
        call load_dll(os, dll, errstat, errmsg)
        call print_error(errstat, 'error: loading api_num_links')
        call c_f_procpointer(dll%procaddr, ptr_api_num_links)
        get_num_links = ptr_api_num_links(api)
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name

    end function get_num_links

    function get_num_nodes()

        integer :: get_num_nodes
        character(64) :: subroutine_name

        subroutine_name = 'get_num_nodes'

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ', subroutine_name

        dll%procname = "api_num_nodes"
        call load_dll(os, dll, errstat, errmsg )
        call print_error(errstat, 'error: loading api_num_nodes')
        call c_f_procpointer(dll%procaddr, ptr_api_num_nodes)
        get_num_nodes = ptr_api_num_nodes(api)
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name

    end function get_num_nodes

    function get_node_attr(node_idx, attr)

        integer :: node_idx, attr
        real :: get_node_attr
        character(64) :: subroutine_name

        subroutine_name = 'get_node_attr'

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ', subroutine_name
        print *, node_idx, attr
        if ((attr >= num_node_attributes) .or. (attr < 1)) then
            print *, "error: unexpected node attribute value"
            stop
        end if

        if ((node_idx > N_node) .or. (node_idx < 1)) then
            print *, "error: unexpected node index value"
            stop
        end if

        dll%procname = "api_get_node_attribute"
        call load_dll(os, dll, errstat, errmsg)
        call print_error(errstat, 'error: loading api_get_node_attribute')
        call c_f_procpointer(dll%procaddr, ptr_api_get_node_attribute)
        ! Fortran index starts in 1, whereas in C starts in 0
        get_node_attr = ptr_api_get_node_attribute(api, node_idx-1, attr)
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name

    end function get_node_attr

    function get_link_attr(link_idx, attr)

        integer :: link_idx, attr
        real :: get_link_attr
        character(64) :: subroutine_name

        subroutine_name = 'get_link_attr'

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ', subroutine_name

        if ((attr >= num_link_attributes) .or. (attr < 1)) then
            print *, "error: unexpected link attribute value"
            stop
        end if

        if ((link_idx > N_link) .or. (link_idx < 1)) then
            print *, "error: unexpected link index value"
            stop
        end if

        dll%procname = "api_get_link_attribute"
        call load_dll(os, dll, errstat, errmsg )
        call print_error(errstat, 'error: loading api_get_link_attribute')
        call c_f_procpointer(dll%procaddr, ptr_api_get_link_attribute)
        ! Fortran index starts in 1, whereas in C starts in 0
        get_link_attr = ptr_api_get_link_attribute(api, link_idx-1, attr)
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name

    end function get_link_attr

    function get_link_xsect_attrs(link_idx, xsect_attr)
        real :: get_link_xsect_attrs
        integer :: link_idx, xsect_attr, xsect_type
        character(64) :: subroutine_name

        subroutine_name = 'get_link_xsect_attrs'

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ', subroutine_name

        xsect_type = get_link_attr(link_idx, link_xsect_type)
        if (xsect_type == 3) then ! RECT_CLOSED
            if (xsect_attr == link_geometry) then
                get_link_xsect_attrs = lRectangular
            else if (xsect_attr == link_type) then
                get_link_xsect_attrs = lpipe
            else if (xsect_attr == link_xsect_wMax) then
                get_link_xsect_attrs = get_link_attr(link_idx, link_xsect_wMax)
            else
                get_link_xsect_attrs = nullvalueR
            end if
        else if (xsect_type == 4) then ! RECT_OPEN
            if (xsect_attr == link_geometry) then
                get_link_xsect_attrs = lRectangular
            else if (xsect_attr == link_type) then
                get_link_xsect_attrs = lchannel
            else if (xsect_attr == link_xsect_wMax) then
                get_link_xsect_attrs = get_link_attr(link_idx, link_xsect_wMax)
            else
                get_link_xsect_attrs = nullvalueR
            end if
        else if (xsect_type == 5) then ! TRAPEZOIDAL
            if (xsect_attr == link_geometry) then
                get_link_xsect_attrs = lTrapezoidal
            else if (xsect_attr == link_type) then
                get_link_xsect_attrs = lchannel
            else if (xsect_attr == link_xsect_wMax) then
                get_link_xsect_attrs = get_link_attr(link_idx, link_xsect_yBot)
            else
                get_link_xsect_attrs = nullvalueR
            end if
        else if (xsect_type == 7) then ! PARABOLIC
            if (xsect_attr == link_geometry) then
                get_link_xsect_attrs = lParabolic
            else if (xsect_attr == link_type) then
                get_link_xsect_attrs = lchannel
            else if (xsect_attr == link_xsect_wMax) then
                get_link_xsect_attrs = get_link_attr(link_idx, link_xsect_wMax)
            else
                get_link_xsect_attrs = nullvalueR
            end if
        else
            get_link_xsect_attrs = nullvalueR
        end if

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name

    end function get_link_xsect_attrs

    function get_link_slope(link_idx, node1, node2)
        integer :: link_idx, ltype, node1, node2
        real :: get_link_slope, h, sign_h, x, c_length, oneR
        character(64) :: subroutine_name

        subroutine_name = 'get_link_slope'

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ', subroutine_name

        ltype = get_link_attr(link_idx, link_type)
        oneR = 1.0

        if (ltype == 0) then ! CONDUIT
            c_length = get_link_attr(link_idx, conduit_length)
            h = get_node_attr(node1, node_invertElev) - get_node_attr(node2, node_invertElev)
            sign_h = sign(oneR, h)
            h = abs(h)
            x = sqrt(c_length**2 - h**2)
            get_link_slope = h / x
        else
            get_link_slope = 0
        end if

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name

    end function get_link_slope
end module interface

module interface
    use iso_c_binding
    use dll_mod
    use globals
    implicit none

    private

    public :: populate_tables
    ! public :: print_tables

    ! interface to C DLL
    abstract interface
        function api_initialize(f1, f2, f3, unit_system)
            use, intrinsic :: iso_c_binding
            implicit none
            character(c_char), dimension(*) :: f1, f2, f3
            integer(c_int), value :: unit_system
            type(c_ptr) :: api_initialize
        end function api_initialize

        function api_finalize(api)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
            integer (c_int) :: api_finalize
        end function api_finalize

        function api_get_node_attribute(api, k, attr)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
            integer(c_int) :: k
            integer(c_int) :: attr
            real(c_float) :: api_get_node_attribute
        end function api_get_node_attribute

        function api_get_link_attribute(api, k, attr)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
            integer(c_int) :: k
            integer(c_int) :: attr
            real(c_float) :: api_get_link_attribute
        end function api_get_link_attribute

        subroutine api_print_info(api)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: api
        end subroutine api_print_info

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

    type(os_type) :: os
    type(dll_type) :: dll
    type(c_ptr) :: api

    procedure(api_initialize), pointer :: ptr_api_initialize
    procedure(api_finalize), pointer :: ptr_api_finalize
    procedure(api_get_node_attribute), pointer :: ptr_api_get_node_attribute
    procedure(api_get_link_attribute), pointer :: ptr_api_get_link_attribute
    procedure(api_print_info), pointer :: ptr_api_print_info
    procedure(api_num_links), pointer :: ptr_api_num_links
    procedure(api_num_nodes), pointer :: ptr_api_num_nodes

    contains
        subroutine populate_tables(linkI, nodeI, linkR, nodeR, unit_system)
            integer, dimension(:,:), allocatable, target, intent(in out) :: linkI
            integer, dimension(:,:), allocatable, target, intent(in out) :: nodeI
            real, dimension(:,:), allocatable, target, intent(in out) :: linkR
            real, dimension(:,:), allocatable, target, intent(in out) :: nodeR
            integer(kind = c_int) :: unit_system
            integer :: errstat, ppos, num_args
            character(len = 1024) :: errmsg
            character(len=256) :: inp_file ! absolute path to .inp
            character(len=256) :: rpt_file ! absolute path to .rpt
            character(len=256) :: out_file ! absolute path to .out
            character(len = 256) :: cwd

            ! Initialize C API
            api = c_null_ptr

            call init_os_type(1, os)
            call init_dll(dll)

            ! Get current working directory
            call getcwd(cwd)

            !% Retrieve .inp file path from args
            num_args = command_argument_count()
            if (num_args < 1) then
                print *, "error >>> path to .inp file was not defined"
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

            dll%filename = trim(cwd) // os%pathsep // "libswmm5.so"

            dll%procname = "api_initialize"
            call load_dll(os, dll, errstat, errmsg)
            call print_error(errstat, 'error >>> loading api_initialize')
            call c_f_procpointer(dll%procaddr, ptr_api_initialize)
            api = ptr_api_initialize(inp_file, rpt_file, out_file, unit_system)

            dll%procname = "api_num_links"
            call load_dll(os, dll, errstat, errmsg)
            call print_error(errstat, 'error >>> loading api_num_links')
            call c_f_procpointer(dll%procaddr, ptr_api_num_links)
            N_link = ptr_api_num_links(api)

            dll%procname = "api_num_nodes"
            call load_dll(os, dll, errstat, errmsg )
            call print_error(errstat, 'error >>> loading api_num_nodes')
            call c_f_procpointer(dll%procaddr, ptr_api_num_nodes)
            N_node = ptr_api_num_nodes(api)

            dll%procname = "api_finalize"
            call load_dll(os, dll, errstat, errmsg )
            call print_error(errstat, 'error >>> loading api_finalize')
            call c_f_procpointer(dll%procaddr, ptr_api_finalize)
            errstat = ptr_api_finalize(api)
            if (errstat /= 0) then
                call print_error(errstat, dll%procname)
                stop
            end if

        end subroutine populate_tables
end module interface
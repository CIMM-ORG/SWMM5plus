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
        function API_initialize(api, f1, f2, f3, unit_system)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr) :: api
            character(kind = c_char), dimension(*) :: f1, f2, f3
            integer(c_int) :: unit_system, API_initialize
        end function API_initialize

        function API_finalize(api)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr) :: api
            integer (c_int) :: API_finalize
        end function API_finalize

        function API_get_node_attribute(api, k, attr)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr) :: api
            integer(c_int) :: k
            integer(c_int) :: attr
            real(c_float) :: API_get_node_attribute
        end function API_get_node_attribute

        function API_get_link_attribute(api, k, attr)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr) :: api
            integer(c_int) :: k
            integer(c_int) :: attr
            real(c_float) :: API_get_link_attribute
        end function API_get_link_attribute

        subroutine API_print_info(api)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr) :: api
        end subroutine API_print_info

        function API_num_links(api)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr) :: api
            integer(c_int) :: API_num_links
        end function API_num_links

        function API_num_nodes(api)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr) :: api
            integer(c_int) :: API_num_nodes
        end function API_num_nodes

    end interface

    type(os_type) :: os
    type(dll_type) :: dll
    type(c_ptr) :: api

    procedure(API_initialize), pointer :: ptr_API_initialize
    procedure(API_finalize), pointer :: ptr_API_finalize
    procedure(API_get_node_attribute), pointer :: ptr_API_get_node_attribute
    procedure(API_get_link_attribute), pointer :: ptr_API_get_link_attribute
    procedure(API_print_info), pointer :: ptr_API_print_info
    procedure(API_num_links), pointer :: ptr_API_num_links
    procedure(API_num_nodes), pointer :: ptr_API_num_nodes

    contains
        subroutine populate_tables(linkI, nodeI, linkR, nodeR, unit_system)
            integer, dimension(:,:), allocatable, target, intent(in out) :: linkI
            integer, dimension(:,:), allocatable, target, intent(in out) :: nodeI
            real, dimension(:,:), allocatable, target, intent(in out) :: linkR
            real, dimension(:,:), allocatable, target, intent(in out) :: nodeR
            integer :: unit_system, errstat, ppos, num_args
            character(len = 1024) :: errmsg
            character(len=1) :: path_separator
            character(len=256) :: inp_file ! absolute path to .inp
            character(len=256) :: rpt_file ! absolute path to .rpt
            character(len=256) :: out_file ! absolute path to .out
            character(len = 256) :: cwd

            call init_os_type(1, os)
            call init_dll(dll)

            ! Get current working directory and path separator
            call getcwd(cwd)
            path_separator = cwd(:1)

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

            dll%filename = trim(cwd) // path_separator // "libswmm5.so"

            dll%procname = "API_initialize"
            call load_dll(os, dll, errstat, errmsg)
            call print_error(errstat, 'error >>> loading API_initialize')
            call c_f_procpointer(dll%procaddr, ptr_API_initialize)
            errstat = ptr_API_initialize &
                ( api, &
                  trim(inp_file) // c_null_char, &
                  trim(rpt_file) // c_null_char, &
                  trim(out_file) // c_null_char, &
                  unit_system)
            if (errstat /= 0) then
                call print_error(errstat, dll%procname)
                stop
            end if

            ! dll%procname = "API_num_links"
            ! call load_dll(os, dll, errstat, errmsg)
            ! call print_error(errstat, 'error >>> loading API_num_links')
            ! call c_f_procpointer(dll%procaddr, ptr_API_num_links)
            ! N_link = ptr_API_num_links(api)

            ! dll%procname = "API_num_nodes"
            ! call load_dll(os, dll, errstat, errmsg )
            ! call print_error(errstat, 'error >>> loading API_num_nodes')
            ! call c_f_procpointer(dll%procaddr, ptr_API_num_nodes)
            ! N_node = ptr_API_num_nodes(api)

            print *, "NUM nodes", N_node, "NUM links", N_link
            print *, inp_file
            print *, rpt_file
            print *, out_file

            dll%procname = "API_finalize"
            call load_dll(os, dll, errstat, errmsg )
            call print_error(errstat, 'error >>> loading API_finalize')
            call c_f_procpointer(dll%procaddr, ptr_API_finalize)
            errstat = ptr_API_finalize(api)
            if (errstat /= 0) then
                call print_error(errstat, dll%procname)
                stop
            end if

        end subroutine populate_tables
end module interface
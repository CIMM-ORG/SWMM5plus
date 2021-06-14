module c_library

    use iso_c_binding

    implicit none

    private

    public :: c_lib_type, c_lib_load, c_lib_free

    type c_lib_type
        integer(c_intptr_t) :: fileaddr = 0
        type(c_ptr) :: fileaddrx = c_null_ptr
        type(c_funptr) :: procaddr = c_null_funptr
        character(1024) :: filename = " "
        character(1024) :: procname = " "
    end type c_lib_type

    interface
        function dlopen(filename, mode) bind(c, name="dlopen")
            use iso_c_binding
            implicit none
            type(c_ptr) :: dlopen
            character(c_char), intent(in) :: filename(*)
            integer(c_int), value :: mode
        end function

        function dlsym(handle, name) bind(c, name="dlsym")
            use iso_c_binding
            implicit none
            type(c_funptr) :: dlsym
            type(c_ptr), value :: handle
            character(c_char), intent(in) :: name(*)
        end function

        function dlclose(handle) bind(c, name="dlclose")
            use iso_c_binding
            implicit none
            integer(c_int) :: dlclose
            type(c_ptr), value :: handle
        end function
    end interface

contains

    subroutine c_lib_load(c_lib, errstat, errmsg)
        type (c_lib_type), intent(inout) :: c_lib
        integer, intent(out) :: errstat
        character(*), intent(out) :: errmsg

        errmsg = ''

        c_lib%fileaddrx = dlopen(trim(c_lib%filename) // c_null_char, 1) ! load DLL
        if( .not. c_associated(c_lib%fileaddrx) ) then
            errstat = -1
            write(errmsg, "(A, I2, A)") &
                    'The dynamic library ' // trim(c_lib%filename) // ' could not be loaded.' &
                     //' Check that the file ' // 'exists in the specified location and' &
                     //' that it is compiled for ', (c_intptr_t*8), '-bit systems.'
            return
        end if

        c_lib%procaddr = dlsym(c_lib%fileaddrx, trim(c_lib%procname) // c_null_char)
        if(.not. c_associated(c_lib%procaddr)) then
            errstat = -1
            errmsg = 'The procedure ' // trim(c_lib%procname) // ' in file ' &
                     // trim(c_lib%filename) // ' could not be loaded.'
            return
        end if

        errstat = 0
    end subroutine c_lib_load

    subroutine c_lib_free (c_lib, errstat, errmsg)
        type (c_lib_type), intent(inout) :: c_lib
        integer, intent(out) :: errstat
        character(*), intent(out) :: errmsg

        integer(c_int) :: success

        errstat = -1
        errmsg = ''

        ! close the library:
        success = dlclose( c_lib%fileaddrx )
        if ( success /= 0 ) then
            errstat = -1
            errmsg = 'The dynamic library could not be freed'
            return
        end if
    end subroutine c_lib_free
end module c_library
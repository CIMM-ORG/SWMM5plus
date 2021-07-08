module c_library

    use iso_c_binding
    use define_settings, only: setting

    implicit none

    private

    !% -------------------------------------------------------------------------------
    !% PUBLIC
    !% -------------------------------------------------------------------------------

    public :: c_lib_type
    public :: c_lib_load
    public :: c_lib_close

    !% -------------------------------------------------------------------------------
    !% PRIVATE
    !% -------------------------------------------------------------------------------

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
    !%-----------------------------------------------------------------------------
    !% Description:
    !%    Loads functions from the EPA-SWMM shared library. The name of the
    !%    function that wants to be loaded has to be defined in c_lib%procname
    !%    before executing c_lib_load.
    !%
    !%    For example:
    !%
    !%    c_lib%procname = "api_initialize"
    !%    call c_lib_load(c_lib, errstat, errmsg)
    !%
    !%-----------------------------------------------------------------------------
        type (c_lib_type), intent(inout) :: c_lib ! pointer to shared library
        integer, intent(out) :: errstat ! -1 ir there was an error, 0 if successful
        character(*), intent(out) :: errmsg
        character(64) :: subroutine_name = "c_lib_load"
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%c_library) print *, '*** enter ', subroutine_name

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

        if (setting%Debug%File%c_library) print *, '*** leave ', subroutine_name
    end subroutine c_lib_load

    subroutine c_lib_close (c_lib, errstat, errmsg)
    !%-----------------------------------------------------------------------------
    !% Description:
    !%    Closes the shared library
    !%-----------------------------------------------------------------------------
        type (c_lib_type), intent(inout) :: c_lib
        integer, intent(out) :: errstat ! -1 if error, 0 if successful
        character(*), intent(out) :: errmsg
        character(64) :: subroutine_name = "c_lib_close"
    !%-----------------------------------------------------------------------------

        if (setting%Debug%File%c_library) print *, '*** enter ', subroutine_name

        errmsg = ''
        errstat = dlclose( c_lib%fileaddrx )
        if ( errstat /= 0 ) then
            errstat = -1
            errmsg = 'The dynamic library could not be closed'
            return
        end if

        if (setting%Debug%File%c_library) print *, '*** leave ', subroutine_name
    end subroutine c_lib_close
end module c_library
!-----------------------------------------------------------------------
!module dll_module
!-----------------------------------------------------------------------
module dll_module
    use iso_c_binding
    implicit none
    private ! all by default
    public :: os_type, dll_type, load_dll, free_dll, init_os_type, init_dll, print_error
    ! general constants:
    ! the number of bits in an address (32-bit or 64-bit).
    integer, parameter :: bits_in_addr = c_intptr_t*8
    ! global error-level variables:
    integer, parameter :: errid_none = 0
    integer, parameter :: errid_info = 1
    integer, parameter :: errid_warn = 2
    integer, parameter :: errid_severe = 3
    integer, parameter :: errid_fatal = 4

    integer :: os_id

    type os_type
        character(10) :: endian
        character(len=:), allocatable :: newline
        character(len=:), allocatable :: os_desc
        character(1) :: pathsep
        character(1) :: swchar
        character(11) :: unfform
    end type os_type

    type dll_type
        integer(c_intptr_t) :: fileaddr
        type(c_ptr) :: fileaddrx
        type(c_funptr) :: procaddr
        character(1024) :: filename
        character(1024) :: procname
    end type dll_type

      ! interface to linux API
    interface
        function dlopen(filename,mode) bind(c,name="dlopen")
            ! void *dlopen(const char *filename, int mode);
            use iso_c_binding
            implicit none
            type(c_ptr) :: dlopen
            character(c_char), intent(in) :: filename(*)
            integer(c_int), value :: mode
        end function

        function dlsym(handle,name) bind(c,name="dlsym")
            ! void *dlsym(void *handle, const char *name);
            use iso_c_binding
            implicit none
            type(c_funptr) :: dlsym
            type(c_ptr), value :: handle
            character(c_char), intent(in) :: name(*)
        end function

        function dlclose(handle) bind(c,name="dlclose")
            ! int dlclose(void *handle);
            use iso_c_binding
            implicit none
            integer(c_int) :: dlclose
            type(c_ptr), value :: handle
        end function
    end interface

    contains
        !-----------------------------------------------------------------------
        !Subroutine init_dll
        !-----------------------------------------------------------------------
        subroutine init_dll(dll)
            implicit none
            type(dll_type), intent(inout) :: dll
            dll % fileaddr = 0
            dll % fileaddrx = c_null_ptr
            dll % procaddr = c_null_funptr
            dll % filename = " "
            dll % procname = " "
        end subroutine init_dll

        !-----------------------------------------------------------------------
        !Subroutine init_os_type
        !-----------------------------------------------------------------------
        subroutine init_os_type(os_id,os)
            implicit none
            integer, intent(in) :: os_id
            type(os_type), intent(inout) :: os

            select case (os_id)
            case (1) ! Linux

                os % endian = 'big_endian'
                os % newline = achar(10)
                os % os_desc = 'Linux'
                os % pathsep = '/'
                os % swchar = '-'
                os % unfform = 'unformatted'

            case (2) ! MacOS

                os % endian = 'big_endian'
                os % newline = achar(10)
                os % os_desc = 'MacOS'
                os % pathsep = '/'
                os % swchar = '-'
                os % unfform = 'unformatted'

            case default

            end select

        end subroutine init_os_type

        !-----------------------------------------------------------------------
        !Subroutine load_dll
        !-----------------------------------------------------------------------
        subroutine load_dll (os, dll, errstat, errmsg )
            ! this subroutine is used to dynamically load a dll.


            type (os_type), intent(in) :: os
            type (dll_type), intent(inout) :: dll
            integer, intent( out) :: errstat
            character(*), intent( out) :: errmsg

            integer(c_int), parameter :: rtld_lazy=1
            integer(c_int), parameter :: rtld_now=2
            integer(c_int), parameter :: rtld_global=256
            integer(c_int), parameter :: rtld_local=0

            errstat = errid_none
            errmsg = ''

            select case (os%os_desc)
            case ("Linux","MacOS")
                ! load the dll and get the file address:
                dll%fileaddrx = dlopen( trim(dll%filename)//c_null_char, rtld_lazy )
                if( .not. c_associated(dll%fileaddrx) ) then
                    errstat = errid_fatal
                    write(errmsg,'(i2)') bits_in_addr
                    errmsg = 'the dynamic library '//trim(dll%filename)//' could not be loaded. check that the file '// &
                    'exists in the specified location and that it is compiled for '//trim(errmsg)//'-bit systems.'
                    return
                end if

                ! get the procedure address:
                ! print *, trim(dll%procname)//c_null_char
                dll%procaddr = dlsym( dll%fileaddrx, trim(dll%procname)//c_null_char )
                if(.not. c_associated(dll%procaddr)) then
                    errstat = errid_fatal
                    errmsg = 'the procedure '//trim(dll%procname)//' in file '//trim(dll%filename)//' could not be loaded.'
                    return
                end if

            case ("Windows")
                errstat = errid_fatal
                errmsg = ' load_dll not implemented for '//trim(os%os_desc)

            case default
                errstat = errid_fatal
                errmsg = ' load_dll not implemented for '//trim(os%os_desc)
            end select
            return
        end subroutine load_dll

        !-----------------------------------------------------------------------
        !Subroutine free_dll
        !-----------------------------------------------------------------------
        subroutine free_dll (os, dll, errstat, errmsg )

            ! this subroutine is used to free a dynamically loaded dll
            type (os_type), intent(in) :: os
            type (dll_type), intent(inout) :: dll
            integer, intent( out) :: errstat
            character(*), intent( out) :: errmsg

            integer(c_int) :: success

            errstat = errid_none
            errmsg = ''

            select case (os%os_desc)
            case ("Linux","MacOS")

                ! close the library:
                success = dlclose( dll%fileaddrx )
                if ( success /= 0 ) then
                    errstat = errid_fatal
                    errmsg = 'the dynamic library could not be freed.'
                    return
                else
                    errstat = errid_none
                    errmsg = ''
                end if

            case ("Windows")

                errstat = errid_fatal
                errmsg = ' free_dll not implemented for '//trim(os%os_desc)

            case default
                errstat = errid_fatal
                errmsg = ' free_dll not implemented for '//trim(os%os_desc)
            end select

            return
        end subroutine free_dll

        !-----------------------------------------------------------------------
        !Subroutine print_error
        !-----------------------------------------------------------------------
        subroutine print_error (errstat, procname)
            implicit none
            integer, intent(in) :: errstat
            character(*), intent(in) :: procname
            if (errstat /= 0) then
                write(*, "(A, A, I4)") new_line(""), "Error (" // trim(procname) // "): ", errstat
            end if
        end subroutine print_error
end module dll_module

!-----------------------------------------------------------------------
!Main program
!-----------------------------------------------------------------------
program test_load_dll
    use, intrinsic :: iso_c_binding
    use dll_module
    implicit none

    ! interface to our shared lib
    abstract interface
        function swmm_close()
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int) :: swmm_close
        end function swmm_close

        function swmm_end()
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int) :: swmm_end
        end function swmm_end

        function swmm_getError(errMsg, msgLen)
            use, intrinsic :: iso_c_binding
            character(kind = c_char), dimension(*) :: errMsg
            integer(c_int), intent(in) :: msgLen
            integer(c_int) :: swmm_printInfo
        end function swmm_getError

        function swmm_getMassBalErr(runoffErr, flowErr, qualErr)
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_float) :: runoffErr, flowErr, qualErr
            integer(c_int) :: swmm_getMassBalErr
        end function swmm_getMassBalErr

        function swmm_getVersion()
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int) :: swmm_getVersion
        end function swmm_getVersion

        function swmm_getWarnings()
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int) :: swmm_getWarnings
        end function swmm_getWarnings

        function swmm_open(f1, f2, f3)
            use, intrinsic :: iso_c_binding
            implicit none
            character(kind = c_char), dimension(*) :: f1, f2, f3
            integer(c_int) :: swmm_open
        end function swmm_open

        function swmm_start(saveResults)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), intent(in) :: saveResults
            integer(c_int) :: swmm_start
        end function swmm_start

        function swmm_step(elapsedTime)
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double) :: elapsedTime
            integer(c_int) :: swmm_step
        end function swmm_step
        ! function swmm_report()

        ! end function swmm_report
        ! function swmm_run()

        ! end function swmm_run

        function swmm_printInfo(units)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), intent(in) :: units
            integer(c_int) :: swmm_printInfo
        end function swmm_printInfo
    end interface

    ! The interface runs a SWMM simulation using the C-SWMM engine
    integer :: errstat, num_args, ppos, i
    real(c_double) :: elapsedTime
    real(c_float) :: runoffErr, flowErr, qualErr
    character(len = 1024) :: errmsg
    character(len = 255) :: cwd
    character(len = 1) :: path_separator
    character(len = 256) :: inpfile, rptfile, outfile
    type(os_type) :: os
    type(dll_type) :: dll
    type(c_funptr) :: cfun
    procedure(swmm_open), pointer :: fswmm_open
    procedure(swmm_start), pointer :: fswmm_start
    procedure(swmm_step), pointer :: fswmm_step
    procedure(swmm_end), pointer :: fswmm_end
    procedure(swmm_getMassBalErr), pointer :: fswmm_getMassBalErr
    procedure(swmm_close), pointer :: fswmm_close

    ! Get current working directory and path separator
    call getcwd(cwd)
    path_separator = cwd(:1)
    call init_os_type(1,os)
    call init_dll(dll)

    ! Retrieve .inp file path from args
    num_args = command_argument_count()
    call get_command_argument(1, inpfile)
    ppos = scan(trim(inpfile), '.', back = .true.)
    if (ppos > 0) then
        rptfile = inpfile(1:ppos) // "rpt"
        outfile = inpfile(1:ppos) // "out"
    end if

    dll%filename= trim(cwd) // path_separator // trim("libswmm5.so")

    ! (1) We open the SWMM file (swmm_open)
    !   It is necessary to provide paths for the following files:
    !   - input (.inp)
    !   - report (.rpt)
    !   - output (.out)
    dll%procname = "swmm_open"
    call load_dll(os, dll, errstat, errmsg )
    call print_error(errstat, 'load_swmm_open')
    call c_f_procpointer(dll%procaddr, fswmm_open)
    errstat = fswmm_open(trim(inpfile) // c_null_char, trim(rptfile) // c_null_char, trim(outfile) // c_null_char)
    call print_error(errstat, dll%procname)

    ! (2) SWMM simulation starts (swmm_start)
    dll%procname = "swmm_start"
    call load_dll(os, dll, errstat, errmsg )
    call print_error(errstat, 'load_swmm_start')
    call c_f_procpointer(dll%procaddr, fswmm_start)
    errstat = fswmm_start(1)
    call print_error(errstat, dll%procname)

    ! (3) Run SWMM simulation step by step (swmm_step)
    dll%procname = "swmm_step"
    call load_dll(os, dll, errstat, errmsg )
    call print_error(errstat, 'load_swmm_step')
    call c_f_procpointer(dll%procaddr, fswmm_step)
    elapsedTime = 10
    do while (elapsedTime /= 0)
        errstat = fswmm_step(elapsedTime)
        call print_error(errstat, dll%procname)
    end do

    ! (4) End SWMM simulation (swmm_end)
    dll%procname = "swmm_end"
    call load_dll(os, dll, errstat, errmsg )
    call print_error(errstat, 'load_swmm_end')
    call c_f_procpointer(dll%procaddr, fswmm_end)
    errstat = fswmm_end()
    call print_error(errstat, dll%procname)

    ! (5) Get Mass Balance Error (swmm_getMassBalErr)
    dll%procname = "swmm_getMassBalErr"
    call load_dll(os, dll, errstat, errmsg )
    call print_error(errstat, 'load_swmm_getMassBalErr')
    call c_f_procpointer(dll%procaddr, fswmm_getMassBalErr)
    errstat = fswmm_getMassBalErr(runoffErr, flowErr, qualErr)
    call print_error(errstat, dll%procname)
    print *, "Run-off error: ", runoffErr, " Flow error: ", flowErr, " Qual error: ", qualErr

    ! (6) Close SWMM engine (swmm_close)
    dll%procname = "swmm_close"
    call load_dll(os, dll, errstat, errmsg )
    call print_error(errstat, 'load_swmm_close')
    call c_f_procpointer(dll%procaddr, fswmm_close)
    errstat = fswmm_close()
    call print_error(errstat, dll%procname)

    call free_dll (os, dll, errstat, errmsg)

end program test_load_dll
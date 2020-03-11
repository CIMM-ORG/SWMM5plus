module interface
    use iso_c_binding
    use dll_mod
    implicit none

    ! interface to C DLL
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

        function print_info(units)
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int), intent(in) :: units
        integer(c_int) :: print_info
        end function print_info
    end interface

    type(os_type) :: os
    type(dll_type) :: dll

    procedure(swmm_open), pointer :: fswmm_open
    procedure(swmm_start), pointer :: fswmm_start
    procedure(swmm_step), pointer :: fswmm_step
    procedure(swmm_end), pointer :: fswmm_end
    procedure(swmm_getMassBalErr), pointer :: fswmm_getMassBalErr
    procedure(swmm_close), pointer :: fswmm_close
    procedure(print_info), pointer :: fprint_info

    contains
        subroutine load_swmm_data(inpfile)
            ! The interface runs a SWMM simulation using the C-SWMM engine
            integer :: errstat, ppos, SI, US
            real(c_double) :: elapsedTime
            real(c_float) :: runoffErr, flowErr, qualErr
            character(len = 1024) :: errmsg
            character(len = 255) :: cwd
            character(len = 1) :: path_separator
            character(len = 256) :: inpfile, rptfile, outfile
            ! Define SWMM constants
            US = 0
            SI = 1
            ! Get current working directory and path separator
            call getcwd(cwd)
            path_separator = cwd(:1)
            call init_os_type(1,os)
            call init_dll(dll)

            ppos = scan(trim(inpfile), '.', back = .true.)
            if (ppos > 0) then
                rptfile = inpfile(1:ppos) // "rpt"
                outfile = inpfile(1:ppos) // "out"
            end if

            dll%filename = trim(cwd) // path_separator // "libswmm5.so"

            ! (1) SWMM file is openned (swmm_open)
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

            ! (4) Retrieve information (print_info)
            dll%procname = "print_info"
            call load_dll(os, dll, errstat, errmsg )
            call print_error(errstat, 'load_print_info')
            call c_f_procpointer(dll%procaddr, fprint_info)
            errstat = fprint_info(SI)
            call print_error(errstat, dll%procname)

            ! (5) End SWMM simulation (swmm_end)
            dll%procname = "swmm_end"
            call load_dll(os, dll, errstat, errmsg )
            call print_error(errstat, 'load_swmm_end')
            call c_f_procpointer(dll%procaddr, fswmm_end)
            errstat = fswmm_end()
            call print_error(errstat, dll%procname)

            ! (7) Close SWMM engine (swmm_close)
            dll%procname = "swmm_close"
            call load_dll(os, dll, errstat, errmsg )
            call print_error(errstat, 'load_swmm_close')
            call c_f_procpointer(dll%procaddr, fswmm_close)
            errstat = fswmm_close()
            call print_error(errstat, dll%procname)
            call free_dll (os, dll, errstat, errmsg)

        end subroutine load_swmm_data
end module interface
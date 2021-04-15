program main

    use globals
    use initialization
    use setting_definition, only: setting
    use interface

    implicit none

    integer :: i
    logical :: arg_param = .false.
    character(len=8) :: param
    character(len=256) :: arg

    ! ---  Define paths

    call getcwd(setting%Paths%project)
    setting%Paths%setting = trim(setting%Paths%project) // '/initialization/settings.json'

    ! --- Read args

    do i = 1, iargc()
        call getarg(i, arg)
        if (.not. arg_param) then
            param = arg
            if (i == 1) then
                setting%Paths%inp = arg
            elseif ((trim(arg) == "-s") .or. & ! user provides settings file
                (trim(arg) == "-t")) then  ! hard coded test case
                arg_param = .true.
            else
                write(*, *) 'The argument ' // trim(arg) // ' is unsupported'
                stop
            endif
        else
            arg_param = .false.
            if (trim(param) == '-s') then
                setting%Paths%setting = arg
            elseif (trim(param) == '-t') then
                setting%TestCase%UseTestCase = .true.
                setting%TestCase%TestName = trim(arg)
                if (trim(arg) == 'simple_channel') then
                else if (trim(arg) == 'simple_orifice') then
                else if (trim(arg) == 'simple_pipe') then
                else if (trim(arg) == 'simple_weir') then
                else if (trim(arg) == 'swashes') then
                else if (trim(arg) == 'waller_creek') then
                else if (trim(arg) == 'y_channel') then
                else if (trim(arg) == 'y_storage_channel') then
                else
                    write(*, *) 'The test case ' // trim(arg) // ' is unsupported. Please use one of the following:'
                    print *, new_line("simple_channel")
                    print *, new_line("simple_orifice")
                    print *, new_line("simple_pipe")
                    print *, new_line("simple_weir")
                    print *, new_line("swashes")
                    print *, new_line("waller_creek")
                    print *, new_line("y_channel")
                    print *, new_line("y_storage_channel")
                    stop
                endif
            endif
        endif
    enddo

    ! --- Load Settings

    call load_settings(setting%Paths%setting)

    ! --- Initialization

    call initialize_api()
    call initialize_linknode_arrays()
    call finalize_api()

end program main
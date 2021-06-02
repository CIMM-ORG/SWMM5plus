program main


    use globals
    use assign_index
    use initialization
    use setting_definition, only: setting
    use interface
    use coarray_partition
    use partitioning
    use network_define

    implicit none

    integer :: ii
    logical :: arg_param = .false.
    character(len=8) :: param
    character(len=256) :: arg

    ! --- Load Settings

    call load_settings(setting%Paths%setting)
    if (this_image() == 1) then
        call execute_command_line("if [ -d debug ]; then rm -r debug; fi && mkdir debug")
    end if

    ! --- Read args

    do ii = 1, iargc()
        call getarg(ii, arg)
        if (.not. arg_param) then
            param = arg
            if (ii == 1) then
                if (arg(:1) == '-') then
                    print *, "ERROR: it is necessary to define the path to the .inp file"
                    stop
                end if
                setting%Paths%inp = arg
            elseif ((trim(arg) == "-s") .or. & ! user provides settings file
                    (trim(arg) == "-t")) then  ! hard coded test case
                arg_param = .true.
            elseif (trim(param) == '-v') then
                setting%Verbose = .true.
            else
                write(*, *) 'The argument ' // trim(arg) // ' is unsupported'
                stop
            end if
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
                    print *, new_line('')
                    print *, "simple_channel, simple_orifice, simple_pipe"
                    print *, "simple_weir, swashes, waller_creek"
                    print *, "y_channel, y_storage_channel"
                    stop
                end if
            end if
        end if
    end do

    ! --- Initialization

    call initialize_api()
    call initialize_linknode_arrays()

    if (this_image() == 1) then
        print *, linkI(:, li_P_image)
        print *, nodeI(:, ni_P_image)
    end if

    sync all
    if (this_image() == oneI) then
       call network_initiation()
    endif
    sync all

    ! --- Finalization
    call finalize_api() ! closes link with shared library

    print*, 'End of main'

end program main
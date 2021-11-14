module utility_files

    use define_settings
    USE ifport

    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% file opening, closing, and verification
    !%-----------------------------------------------------------------------------

    private

    public :: util_file_assign_unitnumber
    public :: util_file_get_commandline
    public :: util_file_setup_input_paths_and_files
    public :: util_file_setup_output_folders

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine util_file_assign_unitnumber ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Assigns fortran unit numbers for all files in storage
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'util_file_assign_unitnumber'
        !%-----------------------------------------------------------------------------

        !% --- inp, rpt, out, setting files
        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%inp_file = setting%File%last_unit

        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%rpt_file = setting%File%last_unit

        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%out_file  = setting%File%last_unit

        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%setting_file  = setting%File%last_unit

        !% --- link and node input files
        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%links_input_file  = setting%File%last_unit

        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%nodes_input_file  = setting%File%last_unit

        !% --- debug setup link files
        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%debug_setup_linkR_file  = setting%File%last_unit

        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%debug_setup_linkI_file  = setting%File%last_unit

        !% --- debug setup node files
        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%debug_setup_nodeR_file  = setting%File%last_unit

        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%debug_setup_nodeI_file  = setting%File%last_unit

        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%debug_setup_nodeYN_file  = setting%File%last_unit

        !setting%File%last_unit = setting%File%last_unit+1
        !setting%File%UnitNumber%outputML_combined_file  = setting%File%last_unit

        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%outputML_filename_file  = setting%File%last_unit

        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%outputML_control_file  = setting%File%last_unit

        ! !% -- debug output link and node files
        ! setting%File%last_unit = setting%File%last_unit+1
        ! setting%File%UnitNumber%debug_output_linkR_file  = setting%File%last_unit

        ! setting%File%last_unit = setting%File%last_unit+1
        ! setting%File%UnitNumber%debug_output_linkI_file  = setting%File%last_unit

        ! setting%File%last_unit = setting%File%last_unit+1
        ! setting%File%UnitNumber%debug_output_nodeR_file  = setting%File%last_unit

        ! setting%File%last_unit = setting%File%last_unit+1
        ! setting%File%UnitNumber%debug_output_nodeI_file  = setting%File%last_unit

        ! setting%File%last_unit = setting%File%last_unit+1
        ! setting%File%UnitNumber%debug_output_nodeYN_file  = setting%File%last_unit

        ! setting%File%last_unit = setting%File%last_unit+1
        ! setting%File%UnitNumber%debug_output_elemR_file  = setting%File%last_unit

        ! setting%File%last_unit = setting%File%last_unit+1
        ! setting%File%UnitNumber%debug_output_faceR_file  = setting%File%last_unit

        ! setting%File%last_unit = setting%File%last_unit+1
        ! setting%File%UnitNumber%debug_output_summary_file  = setting%File%last_unit

        ! !% -- swmm5 output link and node files
        ! setting%File%last_unit = setting%File%last_unit+1
        ! setting%File%UnitNumber%swmm5_output_linkR_file  = setting%File%last_unit

        ! setting%File%last_unit = setting%File%last_unit+1
        ! setting%File%UnitNumber%swmm5_output_linkI_file  = setting%File%last_unit

        ! setting%File%last_unit = setting%File%last_unit+1
        ! setting%File%UnitNumber%swmm5_output_nodeR_file  = setting%File%last_unit

        ! setting%File%last_unit = setting%File%last_unit+1
        ! setting%File%UnitNumber%swmm5_output_nodeI_file  = setting%File%last_unit

        ! setting%File%last_unit = setting%File%last_unit+1
        ! setting%File%UnitNumber%swmm5_output_nodeYN_file  = setting%File%last_unit

    end subroutine util_file_assign_unitnumber
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_file_get_commandline ()
        !%-----------------------------------------------------------------------------
        !% Description
        !% reads the command line and stores in setting
        !%-----------------------------------------------------------------------------
        integer :: ii
        character(len=256) :: argtype, argstring
        logical :: need2arg = .false.
        character(64) :: subroutine_name = "util_file_get_commandline"
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // ' [Processor ', this_image(), ']'

        !% --- read the command line arguments (must be in pairs)
        ii = 1

        do while (ii <= iargc())

            !% --- get the first argument at command line (e.g., -i)
            call getarg(ii,argtype)

            !% --- check the arg type
            if (argtype(:1) .ne. '-') then
                if (ii==1) then
                    !% --- if no flag, then the first argument is interpreted as a filename
                    need2arg = .false.
                    argstring = argtype !% the argument is the filename
                    argtype = '-i'  !% the type is the input file
                else
                    !% --- all other arguments are in pairs beginning with a flag
                    write(*,"(A,i3,A)") 'ERROR (USER): command line argument ',ii,' is '//argtype
                    write(*,"(A)") 'Expected a flag (e.g.), -i, -p, -s as the first of a pair: -flag string'
                    stop
                end if
            else
                !% --- check if second argument is needed
                select case (argtype)
                    case ('-i')  !% input file
                        need2arg = .true.
                    case ('-l')  !% swmm library file (libswmm5.so)
                        need2arg = .true.
                    case ('-o')  ! % output path
                        need2arg = .true.
                    case ('-p')  !% project path
                        need2arg = .true.
                    case ('-s')  !% settings file
                        need2arg = .true.
                    case ('-t')  !% hard coded test cases
                        need2arg = .true.
                        print *, 'ERROR (USER/CODE): hard coded test cases (-t option) are not available or need to be revised'
                        stop 63897
                    case ('-v','-von','-voff','-w','-won','-woff')  ! single argument settings
                        need2arg = .false.
                    case default
                        write(*,"(A,i3,A)") 'ERROR (USER): unknown command line argument of '//argtype
                        stop
                end select
            end if

            if (need2arg) then
                !% --- get the second argument
                ii = ii+1
                call getarg(ii,argstring)
                if (argstring(:1) .eq. '-') then
                    write(*,"(A,i3,A)") 'ERROR (USER): command line argument ',ii,' is '//argstring
                    write(*,"(A)") 'Expected a string (e.g.), as the second of a pair: -flag string'
                    stop
                end if
            end if

            !% --- parse commmand line arguments
            !%     Note these are changed to full paths in util_file_setup_input_paths_and_files
            select case (argtype)
                case ('-i')  !% input file
                    setting%File%inp_file = trim(argstring)
                case ('-l')  !% swmm5 library folder
                    setting%File%library_folder = trim(argstring)
                case ('-o')  !% output path
                    setting%File%output_folder = trim(argstring)
                case ('-p')  !% project path
                    setting%File%project_folder = trim(argstring)
                case ('-s')  !% settings file
                    setting%File%setting_file = trim(argstring)
                case ('-t')  !% hard coded test cases
                    print *, 'hard coded test cases need to be revised'
                    stop 63897
                case ('-v','-von')  ! setting.Verbose   on
                    setting%Output%Verbose = .true.
                case ('-voff')  ! setting.Verbose  off
                    setting%Output%Verbose = .false.
                case ('-w','-won')  ! setting.Warning  on
                    setting%Output%Warning = .true.
                case ('-woff')  ! setting.Verbose  off
                    setting%Output%Warning = .false.
                case default
                    write(*,"(A,i3,A)") 'ERROR (USER): unknown command line argument of '//argtype
                    stop
            end select
            ii = ii+1
        end do

        if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine util_file_get_commandline
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_file_setup_input_paths_and_files ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Builds and checks for valid paths to the input, settings, project, and csv files
        !%-----------------------------------------------------------------------------
        integer :: ierr, ios, ireturn, i1, i2
        character(len=256) :: this_purpose
        character(len=256) :: infile_path, project_path, setting_path
        character(len=256) :: default_path, library_path, thisfile
        character(len=8) :: fext
        character(64) :: subroutine_name = "util_file_setup_input_paths_and_files"
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // ' [Processor ', this_image(), ']'

        !% --- Use the current working directory as the base folder
        ierr = getcwd(setting%File%base_folder)

        if (ierr /= 0) then
            write(*,"(A,i5)") 'ERROR (SYSTEM): getcwd() call at start returned error code', ierr
            write(*,"(A)") 'Unexpected system error, location 3799812'
            stop
        end if

        !% --- Start from the values stored in the setting structure
        !%     These may be folder/file names with no context, or complete paths.
        project_path = setting%File%project_folder
        infile_path  = setting%File%inp_file
        setting_path = setting%File%setting_file
        library_path = setting%File%library_folder

        !print *, 'libary
        !% =======================
        !% --- Library folder (for SWMM library)
        default_path = "" ! added on to base_folder
        call util_file_parse_folder_or_file_path ( &
            library_path, setting%File%base_folder, default_path, setting%File%library_folder)
        this_purpose = trim('library folder')
        ireturn = 1
        call util_file_check_if_folder_exist (setting%File%library_folder,this_purpose, ireturn)

        !print *, 'project'
        !% =======================
        !% --- Parse the project folder and path
        default_path = "" ! added on to project_folder
        call util_file_parse_folder_or_file_path ( &
             project_path, setting%File%base_folder, default_path, setting%File%project_folder)
        !% --- check to see if folder path exists
        this_purpose = trim('project folder (-p command line)')
        ireturn = 1
        call util_file_check_if_folder_exist (setting%File%project_folder,this_purpose, ireturn)

        !print *, 'input'
        !% =======================
        !% --- Parse the input file path and filename
        default_path = ""
        call util_file_parse_folder_or_file_path ( &
            infile_path, setting%File%project_folder, default_path, setting%File%inp_file)
        !% --- check to see if input file exists
        this_purpose = 'input file (-i command line)'
        ireturn = 0
        fext = '.inp'
        call util_file_check_if_file_exist ( &
            setting%File%UnitNumber%inp_file, setting%File%inp_file, &
            this_purpose, ireturn, fext)

        !% --- store the kernel of the input file name for use with output file names
        i1 = scan(setting%File%inp_file, '/', back=.true.)
        i2 = scan(setting%File%inp_file, '.inp', back=.true.)
        i1 = i1+1
        i2 = i2-4
        if (i2 > i1) then
            setting%File%input_kernel = setting%File%inp_file(i1:i2)
        else
            if (setting%CaseName%Short .ne. "") then
                setting%File%input_kernel = trim(setting%CaseName%Short)
            else
                setting%File%input_kernel = 'run'
            end if
        end if

        !print *, 'setting'
        !% =======================
        !% --- Parse the settings.json file
        !% --- HACK -- set the default path for settings file to subfolder "definitions"
        !% --- this assumes that SWMM is called from the directory with the source code.
        default_path = './definitions/settings.json'
        call util_file_parse_folder_or_file_path ( &
            setting_path, setting%File%project_folder, default_path,  setting%File%setting_file)
        !% --- check that setting file exists
        this_purpose = 'setting file (-s command line)'
        ireturn = 0
        fext = '.json'
        call util_file_check_if_file_exist ( &
            setting%File%UnitNumber%setting_file, setting%File%setting_file, &
            this_purpose, ireturn, fext)

        !print *, 'link'
        !% =======================
        !% --- links_input.csv and nodes_input.csv
        !% --- if they exist, these must be in the project directory
        !% --- if they don't exist, that will be handled in output.f90
        if (setting%Output%print_links_csv) then
            thisfile = 'links_input.csv'
            default_path = "" ! added on to project_folder
            call util_file_parse_folder_or_file_path ( &
                thisfile, setting%File%project_folder, default_path,  setting%File%links_input_file)
            this_purpose = 'link input csv file'
            ireturn = 1
            fext = '.csv'
            call util_file_check_if_file_exist ( &
                setting%File%UnitNumber%links_input_file, setting%File%links_input_file,&
                this_purpose, ireturn, fext)
            if (ireturn == 2) then
                setting%File%links_input_file_exist = .false.
            end if
        endif

        !print *, 'nodes'
        !% --- nodes_input.csv
        if (setting%Output%print_nodes_csv) then
            thisfile = 'nodes_input.csv'
            default_path = ""  ! added on to project_folder
            call util_file_parse_folder_or_file_path ( &
                thisfile, setting%File%project_folder, default_path,  setting%File%nodes_input_file)
            this_purpose = 'link input csv file'
            ireturn = 1
            fext = '.csv'
            call util_file_check_if_file_exist ( &
                setting%File%UnitNumber%nodes_input_file, setting%File%nodes_input_file, &
                this_purpose, ireturn, fext)
            if (ireturn == 2) then
                setting%File%nodes_input_file_exist = .false.
            end if
        endif

        if (setting%Debug%File%initialization) &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine util_file_setup_input_paths_and_files
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_file_setup_output_folders ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% initializes the output file path and filenames
        !%-----------------------------------------------------------------------------
        integer :: istat, ireturn, ierr
        logical :: isfolder = .false.
        character(len=256) :: cmsg, default_path, output_path, this_purpose
        character(64) :: subroutine_name = "util_file_setup_output_folders"
        !%-----------------------------------------------------------------------------

        output_path = trim(setting%File%output_folder)

        !% =======================
        !% --- Parse the first-level output folder and path
        default_path = ""
        call util_file_parse_folder_or_file_path ( &
            output_path, setting%File%project_folder, default_path, setting%File%output_folder)
        !% --- check to see if folder path exists
        this_purpose = trim('first-level output folder (-o command line)')
        ireturn = 0
        call util_file_check_if_folder_exist (setting%File%output_folder,this_purpose, ireturn)

        !% =======================
        !% --- Parse the time-stamp output folder and path
        !% --- use the input filename kernel and and '_output'
        if (setting%File%input_kernel .ne. "") then
            output_path = trim(setting%File%input_kernel)//'_output'
        else
            output_path = 'output'
        end if
        default_path = ""
        call util_file_parse_folder_or_file_path ( &
            output_path, setting%File%output_folder, default_path, setting%File%output_folder)
        !% --- check to see if folder path exists
        this_purpose = trim('time-stamp output folder (-o command line)')
        ireturn = 0
        call util_file_check_if_folder_exist (setting%File%output_folder,this_purpose, ireturn)

        !% --- use the short case name as the kernel for the time-stamp subfolder
        if (setting%CaseName%Short .ne. "") then
            setting%File%output_kernel = trim(setting%CaseName%Short)
        else
            setting%File%output_kernel = 'run'
        end if

        !% --- set up timestamp subfolder
        setting%File%output_timestamp_subfolder = trim(setting%File%output_folder) // &
            '/' // trim( setting%File%output_kernel) // '_' // trim(setting%Time%DateTimeStamp)

        !% --- check for clash if timestamp subfolder exists
        !% --- HACK --- replace with inquire()
        ierr = chdir(trim(setting%File%output_timestamp_subfolder))
        if (ierr == 0) then
            write(*,"(A)") "ERROR (operational) -- code must create a unique time-stamp output folder,..."
            write(*,"(A)") "...but folder already exists. Either delete folder or wait one minute"
            write(*,"(A)") "...to get a unique time stamp."
            write(*,"(A)") "...Folder attempted to create was..."
            write(*,"(A)") trim(setting%File%output_timestamp_subfolder)
            stop
        end if

        !% --- make directory for timestamp subfolder
        if (this_image() == 1) then
            call execute_command_line (('mkdir '//trim(setting%File%output_timestamp_subfolder)), &
                cmdstat=istat, cmdmsg=cmsg)

            if (istat /= 0) then
                write(*,"(A)") 'ERROR (user, code, or machine) -- attempting to make directory (mkdir)...'
                write(*,"(A)") '...for storing output files in directory...'
                write(*,"(A)") trim(setting%File%output_timestamp_subfolder)
                write(*,"(A)") '...system returned an error message of...'
                write(*,"(A)") trim(cmsg)
                write(*,"(A,i5)") 'The cmdstat returned was ',istat
                stop
            end if
        end if

        this_purpose = 'timestamp subfolder'
        ireturn = 0
        call util_file_check_if_folder_exist (setting%File%output_timestamp_subfolder,this_purpose, ireturn)

        !% --- make directory for timestamp/temp subfolder
        if (setting%File%output_temp_subfolder .eq. "") setting%File%output_temp_subfolder = 'temp'

        setting%File%output_temp_subfolder = trim(setting%File%output_timestamp_subfolder) // &
            '/' // trim( setting%File%output_temp_subfolder)

        if (this_image() == 1) then
            call execute_command_line (('mkdir '//trim(setting%File%output_temp_subfolder)), &
                cmdstat=istat, cmdmsg=cmsg)
            if (istat /= 0) then
                    write(*,"(A)") 'ERROR (user, code, or machine) -- attempting to make directory (mkdir)...'
                    write(*,"(A)") '...for storing output files in directory...'
                    write(*,"(A)") trim(setting%File%output_temp_subfolder)
                    write(*,"(A)") '...system returned an error message of...'
                    write(*,"(A)") trim(cmsg)
                    write(*,"(A,i5)") 'The cmdstat returned was ',istat
                    stop
            end if
        end if

        !% --- combined multi-level output file kernel in temp directory (numbers added before writing)
        if (setting%File%outputML_combinedfile_kernel .eq. "") setting%File%outputML_combinedfile_kernel = 'combined'

        setting%File%outputML_combinedfile_kernel= trim(setting%File%output_temp_subfolder) &
             // '/' //trim(setting%File%outputML_combinedfile_kernel)

        !% ---storage of output filenames (used when StoredFilenames is exceeded)
        if (trim(setting%File%outputML_filename_file) .ne. "") then
            setting%File%outputML_filename_file = trim(setting%File%output_timestamp_subfolder) &
                // '/' // trim(setting%File%outputML_filename_file)
        else
            setting%File%outputML_filename_file = trim(setting%File%output_timestamp_subfolder) &
                // '/output_filenames.txt'
        end if

        !% --- storage of control files for ML output
        if (trim(setting%File%outputML_control_file) .ne. "") then
            setting%File%outputML_control_file = trim(setting%File%output_timestamp_subfolder) &
                // '/' // trim(setting%File%outputML_control_file)
        else
            setting%File%outputML_control_file = trim(setting%File%output_timestamp_subfolder) &
                // '/output_control.unf'
        end if

        !% --- link output files
        if (setting%File%outputML_Link_kernel .eq. "") setting%File%outputML_Link_kernel = 'link'
        setting%File%outputML_Link_kernel = trim(setting%File%output_timestamp_subfolder) &
             // '/' //trim(setting%File%outputML_Link_kernel)

        !% --- node output files
             if (setting%File%outputML_Node_kernel .eq. "") setting%File%outputML_Node_kernel = 'node'
             setting%File%outputML_Node_kernel = trim(setting%File%output_timestamp_subfolder) &
                  // '/' //trim(setting%File%outputML_Node_kernel)


        !% --- setup report and output files
        setting%File%out_file = trim(setting%File%output_timestamp_subfolder) // '/' //trim(setting%File%output_kernel)//'.out'
        setting%File%rpt_file = trim(setting%File%output_timestamp_subfolder) // '/' //trim(setting%File%output_kernel)//'.rpt'

        !% --- HACK -- need error handling for all the mkdir below

        !% --- debug and swmm output folders
        if (setting%Debug%Setup) then

            !% --- setup the subfolders for links and nodes
            setting%File%debug_setup_link_folder = trim(setting%File%output_timestamp_subfolder) &
                // '/' // 'debug_setup/link'
            setting%File%debug_setup_node_folder = trim(setting%File%output_timestamp_subfolder) &
                // '/' // 'debug_setup/node'

            !% --- create directories only using image 1
            if ( this_image() == 1) then
                !% --- create the debug_setup folder
                call execute_command_line ( &
                    ('mkdir '// trim(setting%File%output_timestamp_subfolder) // '/debug_setup'), &
                    cmdstat=istat, cmdmsg=cmsg)
                !% --- create the subfolders for links and nodes
                call execute_command_line ( &
                    ('mkdir '// trim(setting%File%debug_setup_link_folder)), cmdstat=istat, cmdmsg=cmsg)
                call execute_command_line (&
                    ('mkdir '// trim(setting%File%debug_setup_node_folder)), cmdstat=istat, cmdmsg=cmsg)
            end if

        end if
            if (setting%Debug%Output .or. setting%Output%report) then

                !% --- setup the subfolders for links and nodes
                setting%File%debug_output_link_folder = trim(setting%File%output_timestamp_subfolder) &
                    // '/' // 'debug_output/link'

                setting%File%debug_output_node_folder = trim(setting%File%output_timestamp_subfolder) &
                    // '/' // 'debug_output/node'

                !% --- setup subfolders for elemR, faceR, summmary
                setting%File%debug_output_elemR_folder = trim(setting%File%output_timestamp_subfolder) &
                    // '/' // 'debug_output/elemR'

                setting%File%debug_output_faceR_folder = trim(setting%File%output_timestamp_subfolder) &
                    // '/' // 'debug_output/faceR'

                setting%File%debug_output_summary_folder = trim(setting%File%output_timestamp_subfolder) &
                    // '/' // 'debug_output/summary'

                !% --- setup swmm% subfolders for link an node
                setting%File%swmm5_output_link_folder = trim(setting%File%output_timestamp_subfolder) &
                    // '/' // 'swmm5_output/link'

                setting%File%swmm5_output_node_folder = trim(setting%File%output_timestamp_subfolder) &
                    // '/' // 'swmm5_output/node'


                !% --- create directories only using image 1
                if ( this_image() == 1) then
                    !% --- create the debug_output folder
                    call execute_command_line ( &
                        ('mkdir '// trim(setting%File%output_timestamp_subfolder) // '/debug_output'), &
                        cmdstat=istat, cmdmsg=cmsg)

                    !% --- create the subfolders for links and nodes
                    call execute_command_line ( &
                        ('mkdir '// trim(setting%File%debug_output_link_folder)), cmdstat=istat, cmdmsg=cmsg)

                    call execute_command_line (&
                        ('mkdir '// trim(setting%File%debug_output_node_folder)), cmdstat=istat, cmdmsg=cmsg)

                    !% --- create subfolders for elemR, faceR, summmary
                    call execute_command_line (&
                        ('mkdir '// trim(setting%File%debug_output_elemR_folder)), cmdstat=istat, cmdmsg=cmsg)

                    call execute_command_line (&
                        ('mkdir '// trim(setting%File%debug_output_faceR_folder)), cmdstat=istat, cmdmsg=cmsg)

                    call execute_command_line (&
                        ('mkdir '// trim(setting%File%debug_output_summary_folder)), cmdstat=istat, cmdmsg=cmsg)

                    !% --- create the swmm_output folder
                    call execute_command_line ( &
                        ('mkdir '// trim(setting%File%output_timestamp_subfolder) // '/swmm5_output'), &
                        cmdstat=istat, cmdmsg=cmsg)

                    !% --- create swmm% subfolders for link an node
                    call execute_command_line (&
                        ('mkdir '// trim(setting%File%swmm5_output_link_folder)), cmdstat=istat, cmdmsg=cmsg)

                    call execute_command_line (&
                        ('mkdir '// trim(setting%File%swmm5_output_node_folder)), cmdstat=istat, cmdmsg=cmsg)
                end if
            end if
            !if (setting%Debug%Output) then
                ! call system('mkdir debug_output/elemR')
                ! call system('mkdir debug_output/faceR')
                ! call system('mkdir debug_output/summary')
                ! !% >>> BEGIN HACK
                ! call system('mkdir debug_output/swmm5')
                ! call system('mkdir debug_output/swmm5/link')
                ! call system('mkdir debug_output/swmm5/node')
                ! !% >>> END HACK
            !end if
        ! end if

        ! !% --- debug and swmm output files
        ! if (setting%Debug%Setup) then
        !     !% --- link files
        !     setting%File%debug_setup_linkR_file = trim(setting%File%debug_setup_link_folder) // '/linkR.csv'
        !     setting%File%debug_setup_linkI_file = trim(setting%File%debug_setup_link_folder) // '/linkI.csv'
        !     !% --- node files
        !     setting%File%debug_setup_nodeR_file  = trim(setting%File%debug_setup_node_folder) // '/nodeR.csv'
        !     setting%File%debug_setup_nodeI_file  = trim(setting%File%debug_setup_node_folder) // '/nodeI.csv'
        !     setting%File%debug_setup_nodeYN_file = trim(setting%File%debug_setup_node_folder) // '/nodeYN.csv'
        ! end if
        ! if (setting%Debug%Output .or. setting%Output%report) then
        !     write(str_image, '(i5.5)') this_image()
        !     setting%File%debug_output_elemR_file = trim(setting%File%debug_output_elemR_folder) &
        !         // 'i' //trim(str_image) //'_CC_' // trim(ADJUSTL(str_link_node_idx)) &
        !         // '_' // trim(ADJUSTL(str_elem_idx))//'.csv'
        ! end if

    end subroutine util_file_setup_output_folders
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine util_file_parse_folder_or_file_path &
        (input_path, project_path, default_path, full_path)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% changes the input_path to an output path with respect to the project
        !% path and defaults
        !%-----------------------------------------------------------------------------
        integer :: ierr
        character(len=256), intent(inout) :: full_path
        character(len=256), intent(in) :: input_path, project_path, default_path
        character(64) :: subroutine_name = "util_file_parse_folder_or_file_path"
        !%-----------------------------------------------------------------------------

        !print *, subroutine_name
        !print *, 'input   :',input_path
        !print *, 'project :',project_path
        !print *, 'default :', default_path

        !% --- if no project path entered, use the current working directory
        if (project_path .eq. "") then
            ierr = getcwd(project_path)
            if (ierr /= 0) then
                write(*,"(A,i5)") 'ERROR (SYSTEM): getcwd() call at start returned error code', ierr
                write(*,"(A)") 'Unexpected system error, location 654632'
                stop
            end if
        endif

        if (input_path .eq. "") then
            !% --- if no input_path entered, then use the project path plus the default path
            if (default_path .eq. "") then
                full_path = trim(project_path)
            else
                full_path = trim(project_path) // '/' // trim(default_path)
            end if
        elseif (input_path(:1) .eq. '.') then
            !% --- if the path starts with a './' it starts from the project path
            full_path = trim(project_path)//input_path(2:)
        elseif (input_path(:1) .eq. '/') then
            !% --- if path starts with a '/' it must be a valid full path
            full_path = trim(input_path)
        else
            !% --- without . or ./ the output path must be from the project_path directory
            full_path = trim(project_path)//'/'//trim(input_path)
        end if

        !   print *, 'full path', full_path

    end subroutine util_file_parse_folder_or_file_path
!%
!%==========================================================================
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_file_check_if_folder_exist &
        (thisfolder, this_purpose, ireturn)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Checks to see if folder path exists by attempting to change directory
        !% input ireturn = 0 results in stop on failure
        !% input ireturn = 1 results will exit on failure with ireturn = 2
        !% succes is ireturn=0
        !%-----------------------------------------------------------------------------
        character(len=256), intent(in) :: thisfolder, this_purpose
        integer, intent(inout) :: ireturn
        character(len=256) :: cwd_path  !% current working directory path
        character(len=256) :: cmsg
        integer :: ierr, istat
        logical :: folder_exist
        character(64) :: subroutine_name = "util_file_check_if_folder_exist"
        !%-----------------------------------------------------------------------------

        inquire (DIRECTORY=trim(thisfolder),EXIST=folder_exist)

        if (.not. folder_exist) then
            !% --- if ireturn == 0 we don't crash on non-existence
            if (ireturn == 0) then
                if (setting%File%force_folder_creation) then
                    if (this_image() == 1) then
                        call execute_command_line (('mkdir '//trim(thisfolder)), &
                            exitstat=istat, cmdmsg=cmsg)
                        if (istat /= 0) then
                            write(*,"(A)") 'WARNING (user,system) -- attempting to make directory (mkdir)...'
                            write(*,"(A)") trim(thisfolder)
                            write(*,"(A)") '...but system returned an error message of...'
                            write(*,"(A)") trim(cmsg)
                            write(*,"(A,i5)") '...The istat returned was ',istat
                            write(*,"(A)") '...Likely problem is that base directory does not exist and cannot be created.'
                            write(*,"(A)") '...Folder purpose is: '//this_purpose
                            if (thisfolder(:1) == '/') then
                                write(*,"(A)") 'ERROR (user): the directory (see WARNING above) was an absolute directory...'
                                write(*,"(A)") '... so code must stop here.'
                                stop
                            else
                                write(*,"(A)") '...code is continuing using default directories at command line or project folder'
                            end if
                            !stop
                        end if
                    end if
                else
                    write(*,"(A,i3)") 'ERROR (user):  a required folder does not exist...'
                    write(*,"(A)") '... It is likely that the path does not exist and must be created by user.'
                    write(*,"(A)") '...Required folder entered as: '
                    write(*,"(A)") trim(thisfolder)
                    write(*,"(A)") '...Folder purpose is: '//this_purpose
                    stop
                end if
            else
                ireturn = 2
                write(*,"(A,i3)") 'ERROR (USER): a required folder does not exist...'
                write(*,"(A)") '...It is likely that the path does not exist and must be created by user.'
                write(*,"(A)") '...Required folder entered as: '
                write(*,"(A)") trim(thisfolder)
                write(*,"(A)") '...Folder purpose is: '//this_purpose
                stop
            end if
        else
            ireturn = 0
        end if


        ! call getcwd(cwd_path,ierr)
        ! if (ierr /= 0) then
        !     write(*,"(A,i3)") 'ERROR (SYSTEM): getcwd() call when checking a folder path returned error code ',ierr
        !     write(*,"(A)") 'Unexpected system error, location 98733789'
        !     stop
        ! end if

        ! call chdir(thisfolder,ierr)

        ! !% --- change directory unsuccessful
        ! if (ierr /= 0) then
        !     if (ireturn == 0) then
        !         if (setting%File%force_folder_creation) then
        !             if (this_image() == 1) then
        !                 call execute_command_line (('mkdir '//trim(thisfolder)), &
        !                     cmdstat=istat, cmdmsg=cmsg)
        !                 if (istat /= 0) then
        !                     write(*,"(A)") 'ERROR (user, code, or machine) -- attempting to make directory (mkdir)'
        !                     write(*,"(A)") 'for storing output files in directory '//trim(thisfolder)
        !                     write(*,"(A)") 'system returned an error message of...'
        !                     write(*,"(A)") trim(cmsg)
        !                     write(*,"(A,i5)") 'The cmdstat returned was ',istat
        !                 end if
        !             end if
        !         else
        !             write(*,"(A,i3)") 'ERROR (USER): chdir() to a required folder returned error code ',ierr
        !             write(*,"(A)") 'It is likely that the path does not exist and must be created by user.'
        !             write(*,"(A)") 'Required folder entered as: '//trim(thisfolder)
        !             write(*,"(A)") 'Folder purpose is: '//this_purpose
        !             stop
        !         end if
        !     else
        !         ireturn = 2
        !     end if
        ! else
        !     ! go back to current working directory
        !     call chdir(cwd_path,ierr)
        ! end if

    end subroutine util_file_check_if_folder_exist
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_file_check_if_file_exist  &
        (thisunit, thisfilename, thispurpose, ireturn, fext)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Checks to see if file exists by attempting to open and assigning new unit number
        !% ireturn = 0 will cause code to stop if file does not exist
        !% ireturn = 1 will exit with ireturn = 2 if file doess not exist
        !% success provide ireturn =0 for output
        !%-----------------------------------------------------------------------------
        character(len=256), intent(in) :: thisfilename, thispurpose
        character(len=8), intent(in) :: fext
        integer, intent(inout) :: ireturn
        integer, intent(in) :: thisunit
        integer :: ios
        logical :: file_exist
        character(64) :: subroutine_name = "util_file_check_if_file_exist"
        !%-----------------------------------------------------------------------------

        !% ---checks to see if valid file or folder (cannot distinguish between them)
        inquire (FILE=trim(thisfilename),EXIST=file_exist)

        !% --- check the file extension
        if (index(trim(thisfilename),trim(fext),BACK=.true.) == 0) then
            file_exist = .false.
        end if

        if (.not. file_exist) then
            if (ireturn == 0) then
                write(*,"(A)") 'ERROR (USER) file not found. Path or filename may be wrong'
                write(*,"(A)") 'Looking for file '//trim(thisfilename)
                write(*,"(A)") 'File purpose is '//trim(thispurpose)
                write(*,"(A)") 'File should have extension '//trim(fext)
                stop
            else
                ireturn = 2
            end if
        else
            ireturn = 0
        end if

        ! open(thisunit, &
        !     file   = trim(thisfilename),  &
        !     action = 'read', &
        !     iostat = ios)

        ! if (ios /= 0) then
        !     if (ireturn == 0) then
        !         write(*,"(A)") 'ERROR (USER) file not found. Path or filename may be wrong'
        !         write(*,"(A)") 'Looking for file '//trim(thisfilename)
        !         write(*,"(A)") 'File purpose is '//trim(thispurpose)
        !         stop
        !     else
        !         ireturn = 2
        !     end if
        ! else
        !     close(thisunit)
        ! end if

    end subroutine util_file_check_if_file_exist
!%
!%==========================================================================
!%==========================================================================
!% END OF MODULE
!%==========================================================================
!%
end module utility_files

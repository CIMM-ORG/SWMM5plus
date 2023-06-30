module utility_files
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% file opening, closing, and verification
    !%
    !%==========================================================================
    use define_settings
    use ifport
    use utility_crash, only: util_crashpoint

    implicit none

    private

    public :: util_file_assign_unitnumber
    public :: util_file_get_commandline
    public :: util_file_setup_input_paths_and_files
    public :: util_file_duplicate_input
    public :: util_file_delete_duplicate_input
    public :: util_file_setup_output_folders

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine util_file_assign_unitnumber ()
        !%------------------------------------------------------------------
        !% Description:
        !% Assigns fortran unit numbers for all files in storage
        !%------------------------------------------------------------------
        !% --- inp, rpt, out, setting files
        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%inp_file = setting%File%last_unit

        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%rpt_file = setting%File%last_unit

        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%out_file  = setting%File%last_unit

        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%setting_file  = setting%File%last_unit

        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%outputML_filename_file  = setting%File%last_unit

        setting%File%last_unit = setting%File%last_unit+1
        setting%File%UnitNumber%outputML_control_file  = setting%File%last_unit

    end subroutine util_file_assign_unitnumber
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_file_get_commandline ()
        !%------------------------------------------------------------------
        !% Description
        !% reads the command line and stores in setting
        !%------------------------------------------------------------------
        !% Declarations
            integer :: ii
            character(len=256) :: argtype, argstring
            logical :: need2arg = .false.
            character(64) :: subroutine_name = "util_file_get_commandline"
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Debug%File%initialization) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // ' [Processor ', this_image(), ']'
        !%------------------------------------------------------------------

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
                    call util_crashpoint(2298755)
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
                        call util_crashpoint(63897)
                    case ('-R','-v','-von','-voff','-w','-won','-woff')  ! single argument settings
                        need2arg = .false.
                    case default
                        write(*,"(A,i3,A)") 'ERROR (USER): unknown command line argument of '//argtype
                        call util_crashpoint(87895)
                end select
            end if

            if (need2arg) then
                !% --- get the second argument
                ii = ii+1
                call getarg(ii,argstring)
                if (argstring(:1) .eq. '-') then
                    write(*,"(A,i3,A)") 'ERROR (USER): command line argument ',ii,' is '//argstring
                    write(*,"(A)") 'Expected a string (e.g.), as the second of a pair: -flag string'
                    call util_crashpoint(7210987)
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
                    call util_crashpoint(63897)
                case ('-R')  !% Review -- stop after initialization to review
                    setting%Simulation%stopAfterInitializationYN = .true.
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
                    call util_crashpoint(3897483)
            end select
            ii = ii+1
        end do

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%initialization) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_file_get_commandline
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_file_setup_input_paths_and_files ()
        !%------------------------------------------------------------------
        !% Description:
        !% Builds and checks for valid paths to the input, settings, project, and csv files
        !%------------------------------------------------------------------
        !% Declarations
            integer :: ierr, ios, ireturn, i1, i2
            character(len=256) :: this_purpose
            character(len=256) :: infile_path, project_path, setting_path
            character(len=256) :: default_path, library_path, thisfile
            character(len=8) :: fext
            character(64) :: subroutine_name = "util_file_setup_input_paths_and_files"
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%initialization) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // ' [Processor ', this_image(), ']'
        !%------------------------------------------------------------------

        !% --- Use the current working directory as the base folder
        ierr = getcwd(setting%File%base_folder)

        if (ierr /= 0) then
            write(*,"(A,i5)") 'ERROR (SYSTEM): getcwd() call at start returned error code', ierr
            write(*,"(A)") 'Unexpected system error, location 3799812'
            call util_crashpoint(88029873)
        end if

        !% --- Start from the values stored in the setting structure
        !%     These may be folder/file names with no context, or complete paths.
        project_path = setting%File%project_folder
        infile_path  = setting%File%inp_file
        setting_path = setting%File%setting_file
        library_path = setting%File%library_folder

        !% =======================
        !% --- Library folder (for SWMM library)
        default_path = "" ! added on to base_folder
        call util_file_parse_folder_or_file_path ( &
            library_path, setting%File%base_folder, default_path, setting%File%library_folder)
        this_purpose = trim('library folder')
        ireturn = 1
        call util_file_check_if_folder_exist (setting%File%library_folder,this_purpose, ireturn)

        !% =======================
        !% --- Parse the project folder and path
        default_path = "" ! added on to project_folder
        call util_file_parse_folder_or_file_path ( &
             project_path, setting%File%base_folder, default_path, setting%File%project_folder)
        !% --- check to see if folder path exists
        this_purpose = trim('project folder (-p command line)')
        ireturn = 1
        call util_file_check_if_folder_exist (setting%File%project_folder,this_purpose, ireturn)

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

        !% =======================
        !% --- Parse the settings.json file
        if (trim(setting%File%setting_File) == "") then 
            setting%JSON_FoundFileYN = .false.
            write(*,"(A)") "************************************************************************"
            write(*,"(A)") "** setting.json file not specified (-s [filepathname] on command line **"
            write(*,"(A)") "** Running SWMM5+ with default settings and data from *.inp file only **"
            write(*,"(A)") "************************************************************************"
        else
            default_path = 'settings.json'
            call util_file_parse_folder_or_file_path ( &
                setting_path, setting%File%project_folder, default_path,  setting%File%setting_file)
            !% --- check that setting file exists
            this_purpose = 'setting file (-s command line)'
            ireturn = 0
            fext = '.json'
            call util_file_check_if_file_exist ( &
                setting%File%UnitNumber%setting_file, setting%File%setting_file, &
                this_purpose, ireturn, fext)
            if (ireturn == 0) setting%JSON_FoundFileYN = .true.
        end if 

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%initialization) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_file_setup_input_paths_and_files
!%
!%==========================================================================
!%==========================================================================
!%   
    subroutine util_file_duplicate_input ()
        !%------------------------------------------------------------------
        !% Description:
        !% Creates a duplicate input file for each image
        !% This is simpler than having a single read-in and broadcast of
        !% the input data.
        !% HACK -- need to consider better ways to handle very large input
        !% files
        !%------------------------------------------------------------------

        if (setting%File%duplicate_input_file) then
            !% --- create separate input file copies for each image (handled by image=1)
            if (this_image()==1) call util_file_input_copies (.true.,this_image())
            sync all
            !% --- reset the name of the input file for each image
            call util_file_input_copies (.false.,this_image())
        end if
        
    end subroutine util_file_duplicate_input
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_file_delete_duplicate_input ()
        !%------------------------------------------------------------------
        !% Description
        !% Deletes the duplicate input files in the tmp directory
        !% Each image will delete its own temporary input file
        !% Only applies if there is more than one image
        !%------------------------------------------------------------------

        if (setting%File%duplicate_input_file) then 
            if (num_images() > 1) then
                call execute_command_line (('rm '//setting%File%inp_file))
            end if
        end if

    end subroutine util_file_delete_duplicate_input
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_file_setup_output_folders ()
        !%------------------------------------------------------------------
        !% Description:
        !% initializes the output file path and filenames
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: istat, ireturn, ierr, kk, lc
            logical :: isfolder = .false.
            character(len=256) :: cmsg, default_path, output_path, this_purpose, temppath
        !%------------------------------------------------------------------

        output_path = trim(setting%File%output_folder)
        !% --- remove the last character of the path if it is a '/' or '\
        lc = len_trim(output_path)

        if (output_path(lc:lc) == '/') output_path(lc:lc) = ' '
        if (output_path(lc:lc) == '\') output_path(lc:lc) = ' '

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

        !% --- set up timestamp subfolder name (provisional)
        setting%File%output_timestamp_subfolder = trim(setting%File%output_folder) // &
            '/' // trim( setting%File%output_kernel) // '_' // trim(setting%Time%DateTimeStamp)

        !% --- check if provisional timestamp subfolder already exists and modify if it does exist
        if (this_image() == 1) then
            inquire (DIRECTORY=trim(setting%File%output_timestamp_subfolder),EXIST=isfolder)
            if (isfolder) then
                !% --- loop over ASCII lower case characters and try adding
                do kk=97,122
                    setting%File%output_timestamp_subfolder = trim(setting%File%output_folder) // &
                        '/' // trim(setting%File%output_kernel) // '_' // trim(setting%Time%DateTimeStamp) &
                        // '_' // char(kk)
                    inquire (DIRECTORY=trim(setting%File%output_timestamp_subfolder),EXIST=isfolder)
                    if (.not. isfolder) then
                        !% --- accept this name and leave the loop
                        exit
                    else
                        !% continue
                    end if
                end do
                !% --- if all 26 letters don't work, there's something strange going on, so stop execution.
                if (isfolder) then
                    write(*,"(A)") "ERROR (operational) -- code must create a unique time-stamp output folder,..."
                    write(*,"(A)") "...but folder already exists. Either delete folder or wait one minute"
                    write(*,"(A)") "...to get a unique time stamp."
                    write(*,"(A)") "...Last folder code attempted to create was..."
                    write(*,"(A)") trim(setting%File%output_timestamp_subfolder)
                    call util_crashpoint(983751)
                else
                    !% continue
                end if
            else
                !% continue
            end if
        else
            !% continue
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
                call util_crashpoint(736752)
            else
                !% continue
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
                    call util_crashpoint(439870)
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

    end subroutine util_file_setup_output_folders
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%   
    subroutine util_file_input_copies (createYN, imageIn)
        !%------------------------------------------------------------------
        !% Description:
        !% Handles separate EPA SWMM *.in files for each image
        !% if createYN is true, then creates a \tmp folder and creates
        !% copies of *.inp file in each folder
        !% if createYN is false then it stores the name of the imageIn
        !% input file in setting%File%inp_file
        !%------------------------------------------------------------------
        !% Declarations
            logical, intent(in) :: createYN 
            integer, intent(in) :: imageIn
            character(len=256) :: thisfolder, kernel, extension 
            character(len=256) :: newfolder, newfilename
            character(len=1)   :: divider
            integer :: ii
            logical :: doesexist
        !%------------------------------------------------------------------
        !% Preliminaries
            if (num_images() == 1) return  ! skip for a single image
        !%------------------------------------------------------------------
        !% --- find the upper level directory where the input file sits
        divider = '/'   
        call util_file_separate_kernel_and_extension    &
             (setting%File%inp_file, thisfolder, divider, extension)

        !% --- name for a temporary folder for storing input file copies
        newfolder = trim(thisfolder) // divider // 'tmp'
    
        !% --- check to see if folder exists, create it if not
        inquire(DIRECTORY=newfolder,EXIST=doesexist)
        if (.not. doesexist) then
            if (createYN) then
                call execute_command_line( ('mkdir '//trim(newfolder)) )
            else
                print *, 'CODE ERROR: need to create the tmp file before getting filenames'
                call util_crashpoint(449872)
            end if
        end if

        !% --- glue the filename back to the folder
        newfilename = trim(newfolder)//divider//trim(extension)

        !% --- split the filename extension off 
        divider = '.'   
        call util_file_separate_kernel_and_extension    &
             (newfilename, kernel, divider, extension)

        !% --- cycle through images
        do ii=1,num_images()
            call util_file_new_inp_filename_for_image &
                (newfilename, kernel, divider, extension, ii)

            if (createYN) then
                !% --- make a copy of the input file for this image
                call execute_command_line (('cp '//trim(setting%File%inp_file)//' '//trim(newfilename)))
            else
                if (ii==imageIn) then
                    !% --- store the newfilename as the input file
                    setting%File%inp_file  = trim(newfilename)
                end if
            end if
        end do

    end subroutine util_file_input_copies
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_file_separate_kernel_and_extension &
        (filename, kernel, divider, extension)
        !%------------------------------------------------------------------
        !% Description
        !% Takes the filename and finds the location of right-most "."
        !% and returns the kernel as the text to the left and the extension
        !% as the text to the right
        !%------------------------------------------------------------------
        !% Declarations
            character(len=*), intent(in) :: filename
            character(len=*), intent(inout) :: kernel, extension
            character(len=1), intent(in) :: divider
            integer :: jj
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        jj=scan(filename,divider, .true.)
        extension   = trim(filename(jj+1:))
        kernel      = filename(1:jj-1)

    end subroutine util_file_separate_kernel_and_extension
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_file_new_inp_filename_for_image &
        (newfilename, kernel, divider, extension, timage)
        !%------------------------------------------------------------------
        !% Description:
        !% creates a filename that inserts the image # after the kernel
        !% and before the extension
        !%------------------------------------------------------------------
        !% Declarations
            character(len=*), intent(inout) :: newfilename
            character(len=*), intent(in)    :: kernel, extension
            character(len=*), intent(in)    :: divider
            integer, intent(in) :: timage
        !%------------------------------------------------------------------
        if (num_images() > 99999) then
            !% --- the maximum number of images must be consistent with the formatting
            !%      in setting the newfilename
            print *, 'USER/CODE ERROR: The number of images (',num_images(),') exceeds the'
            print *, 'format size used for setting up input files. '
            print *, 'The quick fix is to limit the number of images to ',99999
            call util_crashpoint(31937)
        else
            write(newfilename,"(A,A,i5.5,A,A)") trim(kernel),'_',timage,trim(divider),trim(extension)
        end if

    end subroutine util_file_new_inp_filename_for_image
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_file_parse_folder_or_file_path &
        (input_path, project_path, default_path, full_path)
        !%------------------------------------------------------------------
        !% Description:
        !% changes the input_path to an output path with respect to the project
        !% path and defaults
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: ierr
            character(len=256), intent(inout) :: full_path
            character(len=256), intent(in) :: input_path, project_path, default_path
        !%------------------------------------------------------------------

        !% --- if no project path entered, use the current working directory
        if (project_path .eq. "") then
            ierr = getcwd(project_path)
            if (ierr /= 0) then
                write(*,"(A,i5)") 'ERROR (SYSTEM): getcwd() call at start returned error code', ierr
                write(*,"(A)") 'Unexpected system error, location 654632'
                call util_crashpoint(88133387)
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

    end subroutine util_file_parse_folder_or_file_path
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_file_check_if_folder_exist &
        (thisfolder, this_purpose, ireturn)
        !%------------------------------------------------------------------
        !% Description:
        !% Checks to see if folder path exists by attempting to change directory
        !% input ireturn = 0 results in stop on failure
        !% input ireturn = 1 results will exit on failure with ireturn = 2
        !% succes is ireturn=0
        !%------------------------------------------------------------------
        !% Declarations:
            character(len=256), intent(in) :: thisfolder, this_purpose
            integer, intent(inout) :: ireturn
            character(len=256) :: cwd_path  !% current working directory path
            character(len=256) :: cmsg
            integer :: ierr, istat
            logical :: folder_exist
        !%------------------------------------------------------------------

        inquire (DIRECTORY=trim(thisfolder),EXIST=folder_exist)

        if (.not. folder_exist) then
            !% --- if ireturn == 0 we don't crash on non-existence
            if (ireturn == 0) then
                if (setting%File%force_folder_creationYN) then
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
                                call util_crashpoint(559228)
                            else
                                write(*,"(A)") '...code is continuing using default directories at command line or project folder'
                            end if
                        end if
                    end if
                else
                    write(*,"(A,i3)") 'ERROR (user):  a required folder does not exist...'
                    write(*,"(A)") '... It is likely that the path does not exist and must be created by user.'
                    write(*,"(A)") '...Required folder entered as: '
                    write(*,"(A)") trim(thisfolder)
                    write(*,"(A)") '...Folder purpose is: '//this_purpose
                    call util_crashpoint(4422878)
                end if
            else
                ireturn = 2
                write(*,"(A,i3)") 'ERROR (USER): a required folder does not exist...'
                write(*,"(A)") '...It is likely that the path does not exist and must be created by user.'
                write(*,"(A)") '...Required folder entered as: '
                write(*,"(A)") trim(thisfolder)
                write(*,"(A)") '...Folder purpose is: '//this_purpose
                call util_crashpoint(88198722)
            end if
        else
            ireturn = 0
        end if

    end subroutine util_file_check_if_folder_exist
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_file_check_if_file_exist  &
        (thisunit, thisfilename, thispurpose, ireturn, fext)
        !%------------------------------------------------------------------
        !% Description:
        !% Checks to see if file exists by attempting to open and assigning new unit number
        !% ireturn = 0 will cause code to stop if file does not exist
        !% ireturn = 1 will exit with ireturn = 2 if file doess not exist
        !% success provide ireturn =0 for output
        !%------------------------------------------------------------------
        !% Declarations:
            character(len=256), intent(in) :: thisfilename, thispurpose
            character(len=8), intent(in) :: fext
            integer, intent(inout) :: ireturn
            integer, intent(in) :: thisunit
            integer :: ios
            logical :: file_exist
        !%------------------------------------------------------------------

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
                call util_crashpoint(32234)
            else
                ireturn = 2
            end if
        else
            ireturn = 0
        end if

    end subroutine util_file_check_if_file_exist
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module utility_files

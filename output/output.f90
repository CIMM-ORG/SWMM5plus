module output

    use define_indexes
    use define_keys
    use define_globals
    use define_settings
    use define_types
    use interface
    use utility_datetime
    !use, intrinsic :: iso_fortran_env, only: *

    implicit none

    private

    public output_read_csv_link_names
    public output_read_csv_node_names
    public output_create_link_files
    public output_create_node_files
    public output_write_link_files
    public output_write_node_files
    public output_combine_links
    public output_move_node_files
    public output_update_swmm_out

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine output_read_csv_link_names()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% reading the link input file and storing it in link_output_idx
        !%-----------------------------------------------------------------------------
        integer              :: rc, fu, pp, jj, kk, link_idx, phantom_counter
        logical              :: no_file = .false.
        logical              :: file_exists, endoffile
        character(len = 250) :: link_name
        character(len=16)    :: thispos
        character(64)        :: subroutine_name = 'output_read_csv_link_names'
        !%--------------------------------------------------------------------------
        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        !% --- abandon procedure if printout of links not needed
        link_output_idx = nullvalueI
        if (.not. setting%Output%print_links_csv) return  

        !% --- check to see if file exists
        inquire (FILE=setting%File%links_input_file, EXIST=file_exists)

        if (file_exists) then
            !% open the csv file of link names
            open(unit= setting%File%UnitNumber%links_input_file, &
                 file=trim(setting%File%links_input_file), &
                 action='read', &
                 iostat=rc)
            if (rc /= 0) then
                write (*, '(3a, i0)') 'ERROR (user): Opening file ', trim(setting%File%links_input_file), ' failed: ', rc
                stop "in " // subroutine_name
            end if 
            pp = 1 ! parent link
            endoffile = .false.
            !% ---loop through till the end of the file and save the valid links
            do while ((.not. endoffile) .and. (pp .le. setting%Output%max_links_csv))
                inquire (unit=setting%File%UnitNumber%links_input_file, position=thispos)
                if (thispos .eq. 'APPEND') then
                    endoffile = .true.
                    pp = pp-1 !% so that pp=0 indicates nothing read (empty file)
                    exit !end the do loop
                end if
                !% --- read in the link name from the csv
                read(setting%File%UnitNumber%links_input_file, *,  iostat = rc) link_name
                !% --- crash on error
                if (rc /= 0) then
                    close(setting%File%UnitNumber%links_input_file)
                    !exit
                    write(*,"(A)") 'ERROR (user): reading file ', trim(setting%File%links_input_file)
                    write(*,"(A,i5)") 'failed before end of file with error ',rc
                    stop "in " // subroutine_name
                end if

                !% --- converting link name to link idx using the interface
                link_idx = interface_find_object(object_type=API_LINK, object_name = link_name)
                !% --- crash on error
                if (link_idx == 0) then
                    write(*, "(A)") "ERROR (user): Link " // trim(link_name) // " in " // &
                        trim(setting%File%links_input_file) // " couldn't be found"
                    !exit
                    stop "in " // subroutine_name
                end if

                !% --- store index of link for output and increase index
                link_output_idx(pp) = link_idx
                pp = pp + 1

                !% checking if the link is split across processors if
                !% so then store the id of the phantom link for output
                phantom_counter = 0
                if (link%I(link_idx, li_parent_link) == link_idx) then
                    do jj = SWMM_N_link+1, N_link
                        if (link%I(jj, li_parent_link) == link_idx) then
                            link_output_idx(pp+phantom_counter+1) = jj
                            phantom_counter = phantom_counter + 1
                        end if
                    end do
                end if
                pp = pp + phantom_counter + 1
                link%I(link_idx,li_num_phantom_links) = phantom_counter

                !% --- check for stopping due to max links setting.
                if (pp .ge. setting%Output%max_links_csv) then
                    if (setting%Output%Warning) then
                        write(*,"(A)") 'WARNING: stopped reading links_input_file file due to excessive number of links'
                        write(*,"(A)") 'Filename = ',trim(setting%File%links_input_file)
                        write(*,"(A,i5)") 'Maximum links set by setting.Output.max_links_csv as: ',setting%Output%max_links_csv
                    end if
                end if
            end do
            !% --- if exited without reading (empty file) send warning
            if (pp == 0) then
                if (setting%Output%Warning) then
                    write(*,"(A)") 'WARNING: did not find an links in the links_input_file'
                    write(*,"(A)") 'Filename = ',trim(setting%File%links_input_file)
                end if
            end if
        end if

        if (.not. file_exists) then
            !% --- if links_input_file is not specified we output all the links up to the maximum allowed
            pp = 1 !% parent link
            !do link_idx = 1, SWMM_N_link
            do while ( (pp <= SWMM_N_link) .and. (pp .le. setting%Output%max_links_csv))   
                phantom_counter = 0
                link_output_idx(pp) = link_idx
                !% --- only parent links have associated phantoms
                if (link%I(link_idx, li_parent_link) == link_idx) then
                    do jj = SWMM_N_link+1, N_link
                        if (link%I(jj, li_parent_link) == link_idx) then
                            link_output_idx(pp+phantom_counter+1) = jj
                            phantom_counter = phantom_counter + 1
                        end if
                    end do
                end if
                pp = pp + phantom_counter + 1
                link%I(link_idx,li_num_phantom_links) = phantom_counter

                !% --- check for stopping due to max links setting.
                if (pp .ge. setting%Output%max_links_csv) then
                    if (setting%Output%Warning) then
                        write(*,"(A)") 'WARNING: stopped selecting links for csv output due to excessive number of links'
                        write(*,"(A,i5)") 'Maximum links set by setting.Output.max_links_csv as: ',setting%Output%max_links_csv
                    end if
                end if            
            end do            
        end if

        if (setting%Debug%File%output) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine output_read_csv_link_names
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine output_read_csv_node_names()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Reading the node input file and store nodes in note_output_idx.
        !% If file does not exist, then choose the first 1:setting%Output%max_nodes_csv 
        !% nodes for output.
        !%-----------------------------------------------------------------------------
        character(len = 250) :: node_name
        integer :: rc, fu, ii, node_idx, maxnode
        character(16) :: thispos
        logical :: file_exists, endoffile
        character(64) :: subroutine_name = 'output_read_csv_node_names'
        !%--------------------------------------------------------------------------
        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        !% --- abandon procedure if printout of nodes not needed
        node_output_idx = nullvalueI
        if (.not. setting%Output%print_nodes_csv) return    

        !% --- check to see if file exists
        inquire (FILE=setting%File%nodes_input_file, EXIST=file_exists)    

        !% --- read the file that exists
        if (file_exists) then
            !% --- open the csv file of link names
            open(unit= setting%File%UnitNumber%nodes_input_file, &
                 file=trim(setting%File%nodes_input_file), &
                 action='read', &
                 iostat=rc)
            if (rc /= 0) then
                write (*, '(3a, i0)') 'ERROR (user): Opening file ', trim(setting%File%nodes_input_file), ' failed: ', rc
                stop "in " // subroutine_name
            end if 
            ii = 1
            endoffile = .false.
            do while ((.not. endoffile) .and. (ii .le. setting%Output%max_nodes_csv))
                inquire (unit=setting%File%UnitNumber%nodes_input_file, position=thispos)
                if (thispos .eq. 'APPEND') then
                    endoffile = .true.
                    ii = ii-1 !% so that ii=0 indicates nothing read (empty file)
                    exit !end the do loop
                end if  
                !% --- read in nodes              
                read(setting%File%UnitNumber%nodes_input_file, *, iostat = rc) node_name
                if (rc /= 0) then
                    node_output_idx(ii:) = nullvalueI
                    close(setting%File%UnitNumber%nodes_input_file)
                    !exit
                    write(*,"(A)") 'ERROR (user): reading file ', trim(setting%File%nodes_input_file)
                    write(*,"(A,i5)") 'failed before end of file with error ',rc
                    stop "in " // subroutine_name
                end if
                !% --- converting node name to node idx using the interface
                node_idx = interface_find_object(object_type=API_NODE, object_name = node_name)
                if (node_idx == 0) then
                    write(*, "(A)") "Node " // trim(node_name) // " in " // &
                    trim(setting%File%nodes_input_file) // " couldn't be found"
                    stop "in " // subroutine_name
                end if
                node_output_idx(ii) = node_idx
                ii = ii + 1
            end do
            !% --- if exited without reading (empty file) send warning
            if (ii == 0) then
                if (setting%Output%Warning) then
                    write(*,"(A)") 'WARNING: did not find an nodes in the nodes_input_file'
                    write(*,"(A)") 'Filename = ',trim(setting%File%nodes_input_file)
                end if
            end if        
        end if    

        !% --- store 1:setting%Output%max_nodes_csv nodes for output when file doesn't exist
        if (.not. file_exists) then
            !% --- Output all nodes (up to max) if CSV file
            maxnode = max( SWMM_N_node, setting%Output%max_nodes_csv)
            node_output_idx = (/ (ii, ii =1, maxnode)/)
            !% --- send warning if output is truncated
            if (SWMM_N_node > maxnode) then
                write(*,"(A)") 'WARNING: stopped selecting nodes for csv output due to excessive number of nodes'
                write(*,"(A,i5)") 'Maximum nodes set by setting.Output.max_nodes_csv as: ',setting%Output%max_links_csv
            end if
        end if

        if (setting%Debug%File%output) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

    end subroutine output_read_csv_node_names
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine output_create_link_files
        !%-----------------------------------------------------------------------------
        !% Description:
        !%   Creation of link files and header for the files
        !%-----------------------------------------------------------------------------
        integer :: ii,fu, open_status, temp_link_idx
        character(len = 250) :: file_name
        character(len = 100) :: link_name
        character(len = 5)   :: str_image
        character(len = 10)  :: str_idx
        character(64) :: subroutine_name = 'output_create_link_files'
        !%--------------------------------------------------------------------------
        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        write(str_image, '(i5.5)') this_image()

        do ii=1, size(link%P%have_output)

            !% HACK needed to ensure that a valid index is provided
            if (link%P%have_output(ii) > 0) then

                !% check if the link is a phantom link and if so find original link name and open correct file for the correct processor
                !% otherwise open file in the usual format of "link name_imageID.csv"
                if (link%P%have_output(ii) > size(link%names(:))) then
                    write(str_idx, '(i5.5)') link%P%have_output(ii)
                    temp_link_idx = link%P%have_output(ii)
                    !file_name = "debug_output/link/"//trim(link%names(link%I(temp_link_idx,li_parent_link))%str) &
                    !    //"_"//trim(str_image)//"_"//trim(str_idx)//".csv"
                    file_name = trim(setting%File%debug_output_link_folder) &
                        //trim(link%names(link%I(temp_link_idx,li_parent_link))%str) &
                        //"_"//trim(str_image)//"_"//trim(str_idx)//".csv"    
                else
                    !file_name = "debug_output/link/"//trim(link%names(link%P%have_output(ii))%str) &
                    !    //"_"//trim(str_image)//".csv"
                    file_name = trim(setting%File%debug_output_link_folder) &
                        //trim(link%names(link%P%have_output(ii))%str) &
                        //"_"//trim(str_image)//".csv"
                end if

                open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
                    form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (*, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                    stop "in " // subroutine_name
                end if

                !% Write the header of the file, set end for next write and then close file
                write(fu, *) "Timestamp,Time_In_Secs,flowrate"
                endfile(fu)
                close(fu)
            end if    
        end do

        if (setting%Debug%File%output) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine output_create_link_files
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine output_create_node_files
        !%-----------------------------------------------------------------------------
        !% Description:
        !%   Creation of node files and header for the files
        !%-----------------------------------------------------------------------------
        integer :: ii,fu, open_status
        character(len = 250) :: file_name
        character(len = 100) :: node_name
        character(len = 5)   :: str_image
        character(64) :: subroutine_name = 'output_create_node_files'
        !%--------------------------------------------------------------------------
        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        !% Get current image as a string
        write(str_image, '(i5.5)') this_image()

        do ii=1, size(node%P%have_output) 

            !% HACK -- needed because node%P%have_output(ii)= 0 is assigned somehwere
            if (node%P%have_output(ii) > 0 ) then

                !% Open the node file
                !file_name = "debug_output/node/"//trim(node%names(node%P%have_output(ii))%str) &
                !    //"_"//trim(str_image)//".csv"
                file_name = trim(setting%File%debug_output_node_folder) &
                    //trim(node%names(node%P%have_output(ii))%str) &
                    //"_"//trim(str_image)//".csv"   

                open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
                    form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (*, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                    stop "in " // subroutine_name
                end if

                !% Write the header, this endfile and close the file
                write(fu, *) "Timestamp,Time_In_Secs,Head"
                endfile(fu)
                close(fu)
            end if    
        end do

        if (setting%Debug%File%output) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine output_create_node_files
!%
!%==========================================================================
!%==========================================================================
!%  
    subroutine output_write_link_files
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Called at the report time step to calculate the flowrate in the link and 
        !% write it to the file
        !%-----------------------------------------------------------------------------
        integer :: ii, fu, open_status, yr, mnth, dy, hr, min, sec
        integer :: start_elem, end_elem, temp_link_idx
        real(8) :: time_secs, time_epoch, avg_flowrate
        character(len = 250) :: file_name
        character(len = 100) :: link_name
        character(len = 5)   :: str_image
        character(len = 10)  :: str_idx
        character(64) :: subroutine_name = 'output_write_link_files'
        !%--------------------------------------------------------------------------
        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        write(str_image, '(i5.5)') this_image()
        time_secs = setting%Time%Now
        time_epoch = util_datetime_secs_to_epoch(time_secs)
        call util_datetime_decodedate(time_epoch, yr, mnth, dy)
        call util_datetime_decodetime(time_epoch, hr, min, sec)

        do ii=1, size(link%P%have_output)

            if (link%P%have_output(ii) > 0 ) then

                !% store the store the location of the start and end elem for easier reading
                start_elem = link%I(link%P%have_output(ii),li_first_elem_idx)
                end_elem = link%I(link%P%have_output(ii),li_last_elem_idx)

                !% calculate average flowrate by summing up the elems and diving about the number of elems
                avg_flowrate = sum(elemR(start_elem:end_elem,er_Flowrate))/(end_elem-start_elem)

                if (start_elem == end_elem) then
                    avg_flowrate = elemR(start_elem,er_flowrate)
                end if
                !% check if the link is a phantom link and if so find original link name and open correct file for the correct processor
                !% otherwise open file in the usual format of "link name_imageID.csv"
                if (link%P%have_output(ii) > size(link%names(:))) then
                    write(str_idx, '(i5.5)') link%P%have_output(ii)
                    temp_link_idx = link%P%have_output(ii)
                    
                    !file_name = "debug_output/link/"//trim(link%names(link%I(temp_link_idx,li_parent_link))%str) &
                    !    //"_"//trim(str_image)//"_"//trim(str_idx)//".csv"
                    file_name = trim(setting%File%debug_output_link_folder) &
                        //trim(link%names(link%I(temp_link_idx,li_parent_link))%str) &
                        //"_"//trim(str_image)//"_"//trim(str_idx)//".csv"
                else
                    !file_name = "debug_output/link/"//trim(link%names(link%P%have_output(ii))%str) &
                    !    //"_"//trim(str_image)//".csv"
                    file_name = trim(setting%File%debug_output_link_folder) &
                        //trim(link%names(link%P%have_output(ii))%str) &
                        //"_"//trim(str_image)//".csv"
                end if

                open(newunit=fu, file = file_name, status = 'old',access = 'append', &
                    form   = 'formatted', action = 'write', iostat = open_status)

                !% writing timestamped output to file for average flowrate across the link
                write(fu,fmt='(i4, 2(a,i2.2))',advance = 'no') yr,"_",mnth,"_",dy
                write(fu,fmt = '(A)',advance = 'no') '_'
                write(fu,fmt='(2(i2.2,a), i2.2)',advance = 'no') hr,":",min,":",sec
                write(fu,'(A)', advance = 'no') ','
                write(fu, '(F0.16)', advance = 'no') time_secs
                write(fu,'(A)', advance = 'no') ','
                write(fu, '(*(G0.6 : ","))') avg_flowrate

                !% set the end of the file for next write and close file
                endfile(fu)
                close(fu)

            end if
        end do

        if (setting%Debug%File%output) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine output_write_link_files
!%
!%==========================================================================
!%==========================================================================
!%  
    subroutine output_write_node_files
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Called at the report time step to calculate the head in the node and 
        !% write it to the file
        !%-----------------------------------------------------------------------------
        integer :: ii, fu, open_status, yr, mnth, dy, hr, min, sec
        integer :: temp_node_idx
        real(8) :: time_secs, time_epoch, avg_head
        character(len = 250) :: file_name
        character(len = 100) :: link_name
        character(len = 5)   :: str_image
        character(64) :: subroutine_name = 'output_write_node_files'

        !%--------------------------------------------------------------------------
        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        !% converter image ID to string, as well as get current time
        write(str_image, '(i5.5)') this_image()
        time_secs = setting%Time%Now
        time_epoch = util_datetime_secs_to_epoch(time_secs)
        call util_datetime_decodedate(time_epoch, yr, mnth, dy)
        call util_datetime_decodetime(time_epoch, hr, min, sec)

        !% loop through nodes we have to output
        do ii=1, size(node%P%have_output)

            if (node%P%have_output(ii) > 0 ) then

                !% open node file
                !file_name = "debug_output/node/"//trim(node%names(node%P%have_output(ii))%str) &
                !    //"_"//trim(str_image)//".csv"
                file_name = trim(setting%File%debug_output_node_folder) &
                    //trim(node%names(node%P%have_output(ii))%str) &
                    //"_"//trim(str_image)//".csv"

                open(newunit=fu, file = file_name, status = 'old',access = 'append', &
                    form   = 'formatted', action = 'write', iostat = open_status)

                if (open_status /= 0) then
                    write (*, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                    stop "in " // subroutine_name
                end if

                !% temp value for easier to read code
                temp_node_idx = node%P%have_output(ii)

                !% check if the node is a BC type node or a junction type node, then write appropriate data
                if (node%I(temp_node_idx,ni_node_type) == nBCup .or. node%I(temp_node_idx,ni_node_type) == nBCdn) then
                    write(fu,fmt='(i4, 2(a,i2.2))',advance = 'no') yr,"/",mnth,"/",dy
                    write(fu,fmt = '(A)',advance = 'no') ' '
                    write(fu,fmt='(2(i2.2,a), i2.2)',advance = 'no') hr,":",min,":",sec
                    write(fu,'(A)', advance = 'no') ','
                    write(fu, '(F0.16)', advance = 'no') time_secs
                    write(fu,'(A)', advance = 'no') ','
                    write(fu, '(*(G0.6 : ","))') faceR(node%I(temp_node_idx,ni_elemface_idx),fr_Head_d)

                else if (node%I(temp_node_idx,ni_node_type) == nJ2 .or. node%I(temp_node_idx,ni_node_type) == nJm) then
                    write(fu,fmt='(i4, 2(a,i2.2))',advance = 'no') yr,"/",mnth,"/",dy
                    write(fu,fmt = '(A)',advance = 'no') ' '
                    write(fu,fmt='(2(i2.2,a), i2.2)',advance = 'no') hr,":",min,":",sec
                    write(fu,'(A)', advance = 'no') ','
                    write(fu, '(F0.16)', advance = 'no') time_secs
                    write(fu,'(A)', advance = 'no') ','
                    write(fu, '(*(G0.6 : ","))') elemR(node%I(temp_node_idx,ni_elemface_idx),er_Head)
                !% if the node type is neither a BC or junction type then print the warning
                else
                    print *, "WARNING: node selected is neither BCup, BCdn, nJ2 or nJM node, no output will be written"
                    print *, "temp_node_idx",temp_node_idx
                    print *, "node%I(temp_node_idx,ni_node_type)", node%I(temp_node_idx,ni_node_type)
                end if

                !% call endfile and close the file
                endfile(fu)
                close(fu)
            end if    
        end do

        if (setting%Debug%File%output) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

    end subroutine output_write_node_files
!%
!%==========================================================================
!%==========================================================================
!%  
    subroutine output_combine_links
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Combines phantom links into single output links
        !%-----------------------------------------------------------------------------
        integer :: ii, jj, pp, rc, open_status, N_parents, N_phantoms
        integer :: temp_link_idx, temp_phantom_link ,link_output_idx_length
        integer :: start_elem, end_elem,num_elems
        integer, allocatable :: file_idx(:)
        real(8)              :: time_secs, time_epoch, flowrate
        real(8), pointer     :: avg_flowrate(:)
        real(8), allocatable :: full_length(:), tt(:)
        logical, allocatable :: first_iteration(:)
        logical, allocatable :: its_over(:)
        character(len = 250) :: parent_file_name, phantom_file_name
        character(len = 250) :: final_file_name
        character(len = 100) :: link_name
        character(len = 5)   :: str_image
        character(len = 10)  :: str_idx
        character(len = 19)  :: str_time
        character(64) :: subroutine_name = 'output_write_link_files'
        !%-----------------------------------------------------------------------------

        link_output_idx_length = count(link_output_idx(:) /= nullvalueI)
        N_phantoms = sum(link%I(:, li_num_phantom_links))
        N_parents = link_output_idx_length - N_phantoms
        
        !print *, link_output_idx_length
        !print *, N_phantoms
        !print *, N_parents
        

        allocate(first_iteration(N_parents))
        allocate(full_length(N_parents))
        allocate(its_over(N_parents))
        allocate(tt(N_parents))
        allocate(file_idx(link_output_idx_length+N_parents))

        avg_flowrate => link%R(:, lr_flowrate)

        first_iteration(:) = .true.
        full_length(:) = zeroR
        avg_flowrate(:) = zeroR
        its_over(:) = .false.
        tt(:) = 1
        file_idx = 0       

        !% only execute if there are actually links to output
        if (sum(link_output_idx(:)) > 0) then 

            do while (any(.not. its_over))
                ii = 1
                pp = 1
                do while (pp <= N_parents)
                    !% Set some starting values
                    !% we use temp_link_idx as the index of the parent link we are getting the output for
                    temp_link_idx = link_output_idx(ii)

                    if (first_iteration(pp)) then
                        !% define parent filename
                        write(str_image, '(i5.5)') link%I(temp_link_idx,li_P_image)

                        !parent_file_name = "debug_output/link/"// &
                        !    trim(link%names(link%I(temp_link_idx,li_parent_link))%str) &
                        !        //"_"//trim(str_image)//".csv"
                        parent_file_name = trim(setting%File%debug_output_link_folder)// &
                                trim(link%names(link%I(temp_link_idx,li_parent_link))%str) &
                                    //"_"//trim(str_image)//".csv"

                        !final_file_name = "swmm5_output/link/"//trim(link%names(temp_link_idx)%str)//".csv"
                        final_file_name = trim(setting%File%swmm5_output_link_folder) &
                            //trim(link%names(temp_link_idx)%str)//".csv"

                        !% Open parent file
                        open(newunit=file_idx(ii), action='read', file=parent_file_name, iostat=rc)
                        if (rc /= 0) then
                            write (*, '(3a, i0)') 'Opening file "', trim(parent_file_name), '" failed: ', rc
                            stop "in " // subroutine_name
                        end if
                        read (file_idx(ii), *, iostat=rc) str_time ! advance one line (skip header)

                        ! start calculating the full_length of the link for the phantom links
                        full_length(pp) = full_length(pp) + link%R(temp_link_idx,lr_AdjustedLength)

                        !Now we open the Final file for the link output which is in a different location and just the name of the link
                        !We also write the header
                        open(newunit=file_idx(link_output_idx_length+pp), &
                            file = final_file_name, status = 'replace',access = 'sequential', &
                            form = 'formatted', action = 'write', iostat = open_status)
                        if (open_status /= 0) then
                            write (*, '(3a, i0)') 'Opening file "', trim(Final_File_NAME), '" failed: ', open_status
                            stop "in " // subroutine_name
                        end if
                        write(file_idx(link_output_idx_length+pp), *) "Timestamp,Time_In_Secs,flowrate"

                        !Now we check if the parent link has phantoms related to it and loop through those values
                        !While looping we open each of the phantom link files and calculate the full_length
                        if (link%I(temp_link_idx, li_num_phantom_links) > 0) then
                            do jj = 1, link%I(temp_link_idx, li_num_phantom_links)
                                temp_phantom_link = link_output_idx(ii+jj)
                                write(str_image, '(i5.5)') link%I(temp_phantom_link,li_P_image)
                                write(str_idx, '(i5.5)')   temp_phantom_link

                                !phantom_file_name = "debug_output/link/"// &
                                !    trim(link%names(link%I(temp_phantom_link,li_parent_link))%str) &
                                !   //"_"//trim(str_image)//"_"//trim(str_idx)//".csv"
                                phantom_file_name =  trim(setting%File%debug_output_link_folder)// &
                                    trim(link%names(link%I(temp_phantom_link,li_parent_link))%str) &
                                    //"_"//trim(str_image)//"_"//trim(str_idx)//".csv"
                                full_length(pp) = full_length(pp) + link%R(temp_phantom_link,lr_AdjustedLength)

                                open(newunit=file_idx(ii+jj), action='read', file=phantom_file_name, iostat=rc)
                                if (rc /= 0) then
                                    write (*, '(3a, i0)') 'Opening file "', trim(phantom_file_name), '" failed: ', rc
                                    stop "in " // subroutine_name
                                end if
                                read (file_idx(ii+jj), *, iostat=rc) str_time ! advance one line (skip header)
                            end do
                            tt(pp) = tt(pp) + 1
                        else
                            !If the link doesn't have any phantom links related to it we close the parent and final file and simply rename and move the parent file to the final file location.
                            !Then increment one.
                            close(file_idx(ii))
                            close(file_idx(link_output_idx_length+pp))
                            call rename(trim(parent_file_name), trim(final_file_name))
                            its_over(pp) = .true.
                        end if
                        first_iteration(pp) = .false.
                    else
                        if (.not. its_over(pp)) then ! here we will only encounter parents with phantom links
                            ! Now that we have the parent file, the final file and the phantom files open
                            ! we can start recombing the flowrates and writing the final file
                            read (file_idx(ii), *, iostat=rc) str_time, time_secs, flowrate
                            if (rc /= 0) then
                                its_over(pp) = .true.
                                ii = ii + link%I(temp_link_idx, li_num_phantom_links) + 1
                                pp = pp + 1
                                cycle
                            end if

                            avg_flowrate(pp) = avg_flowrate(pp) + &
                                (flowrate * link%R(temp_link_idx,lr_AdjustedLength) / full_length(pp))

                            do jj = 1, link%I(temp_link_idx, li_num_phantom_links)
                                temp_phantom_link = link_output_idx(ii+jj)
                                read (file_idx(ii+jj), *, iostat=rc) str_time, time_secs, flowrate
                                avg_flowrate(pp) = avg_flowrate(pp) + &
                                    (flowrate * link%R(temp_phantom_link,lr_AdjustedLength) / full_length(pp))
                            end do

                            write(file_idx(link_output_idx_length+pp), '(A)',     advance = 'no') str_time
                            write(file_idx(link_output_idx_length+pp), '(A)',     advance = 'no') ','
                            write(file_idx(link_output_idx_length+pp), '(F0.16)', advance = 'no') time_secs
                            write(file_idx(link_output_idx_length+pp), '(A)',     advance = 'no') ','
                            write(file_idx(link_output_idx_length+pp), '(*(G0.6 : ","))') avg_flowrate(pp)

                            avg_flowrate(pp) = 0
                            tt(pp) = tt(pp) + 1
                        end if
                    end if
                    !% Now we increment based off of how many phantom links there where related to the parent link that we combined the output
                    ii = ii + link%I(temp_link_idx, li_num_phantom_links) + 1
                    pp = pp + 1
                end do
            end do
            do ii = 1, size(file_idx)
                close(file_idx(ii))
            end do
        end if

        deallocate(first_iteration)
        deallocate(full_length)
        deallocate(its_over)
        deallocate(file_idx)

    end subroutine output_combine_links
!%
!%==========================================================================
!%==========================================================================
!%  
    subroutine output_move_node_files
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Moves the node files to the swmm5 folder
        !%-----------------------------------------------------------------------------
        integer :: ii, fu, open_status
        character(len = 250) :: file_name, file_name_new
        character(len = 100) :: node_name
        character(len = 5)   :: str_image
        character(24) :: timestamp
        character(64) :: subroutine_name = 'output_move_node_files'
        !%--------------------------------------------------------------------------
        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        !% Get current image as a string
        write(str_image, '(i5.5)') this_image()

        do ii=1, size(node%P%have_output)

            if (node%P%have_output(ii) > 0) then

                !% Open the node file
                !file_name = "debug_output/node/"//trim(node%names(node%P%have_output(ii))%str) &
                !    //"_"//trim(str_image)//".csv"
                file_name = trim(setting%File%debug_output_node_folder) &
                    //trim(node%names(node%P%have_output(ii))%str) &
                    //"_"//trim(str_image)//".csv"

                !file_name_new = "swmm5_output/node/"//trim(node%names(node%P%have_output(ii))%str)//".csv"
                file_name_new = trim(setting%File%swmm5_output_node_folder) &
                    //trim(node%names(node%P%have_output(ii))%str)//".csv"  

                call rename(file_name, file_name_new)
            end if
        end do

        if (setting%Debug%File%output) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"

    end subroutine output_move_node_files
!%
!%==========================================================================
!%==========================================================================
!%  
    subroutine output_update_swmm_out()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------
        character(len = 250) :: fname
        character(len = 24) :: timestamp
        logical :: wrote_all_links = .false.
        logical :: wrote_all_nodes = .false.
        integer :: ii, rc, node_idx, link_idx
        integer, allocatable :: fus_nodes(:), fus_links(:)
        real(8) :: node_head, node_result
        real(8) :: link_flowrate, link_result
        real(8) :: timesecs
        character(64) :: subroutine_name = 'output_update_swmm_out'
        !%--------------------------------------------------------------------------
        if (setting%Debug%File%output) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        if (this_image() == 1) then
            allocate(fus_nodes(size(node%P%have_output)))
            allocate(fus_links(size(link%P%have_output)))
            fus_links = nullvalueI
            fus_nodes = nullvalueI

            do while(.not. (wrote_all_links .and. wrote_all_nodes))
  
                if (size(node%P%have_output) < 1) then
                    wrote_all_nodes = .true.
                else
                    do ii=1, size(node%P%have_output)
                        !print *, 'node ',ii,  size(node%P%have_output),  node%P%have_output(ii)
                        node_idx = node%P%have_output(ii)    
                        if (node_idx > 0) then
                            if (fus_nodes(ii) == nullvalueI) then
                                !% open files to process .out
                                !fname = "swmm5_output/node/"//trim(node%names(node_idx)%str)//".csv"
                                fname = trim(setting%File%swmm5_output_node_folder) &
                                    //trim(node%names(node_idx)%str)//".csv"
                                open(action='read', file=trim(fname), iostat=rc, newunit=fus_nodes(ii))
                                read(fus_nodes(ii), *, iostat = rc) timestamp
                            end if

                            read(fus_nodes(ii), "(A,2F10.8)", iostat = rc) timestamp, timesecs, node_head
                            if (rc /= 0) then
                                wrote_all_nodes = .true.
                                close(fus_nodes(ii))
                                !% Write line of .out
                                !exit
                            else   
                                wrote_all_nodes = .false.
                                node_result = node_head - node%R(node_idx,nr_Zbottom)  
                                !% stage values in .out
                                call interface_update_nodeResult(node_idx, api_output_node_depth, node_result)
                            end if 
                        else
                            if (ii .eq. size(node%P%have_output)) then
                                wrote_all_nodes = .true.
                            end if    
                        end if    
                        if (ii .lt. size(node%P%have_output)) then
                            wrote_all_nodes = .false.
                        end if
                    end do
                end if    

                if (size(link%P%have_output) < 1) then
                    wrote_all_links = .true.
                else
                    do ii=1, size(link%P%have_output)
                        !print *, 'link ',ii,  size(link%P%have_output),  link%P%have_output(ii)
                        link_idx = link%P%have_output(ii)
                        if (link_idx > 0) then
                            if (fus_links(ii) == nullvalueI) then
                                !% open files to process .out
                                !fname = "swmm5_output/link/"//trim(link%names(link_idx)%str)//".csv"
                                fname =  trim(setting%File%swmm5_output_link_folder) &
                                    //trim(link%names(link_idx)%str)//".csv"
                                open(action='read', file=trim(fname), iostat=rc, newunit=fus_links(ii))
                                read(fus_links(ii), *, iostat = rc) timestamp
                            end if

                            read(fus_links(ii), "(A,2F10.8)", iostat = rc) timestamp, timesecs, link_flowrate
                            if (rc /= 0) then
                                wrote_all_links = .true.
                                close(fus_links(ii))
                                !% Write line of .out
                                !exit
                            else
                                wrote_all_links = .false.    
                                link_result = link_flowrate
                                !% stage values in .out
                                call interface_update_linkResult(link_idx, api_output_link_flow, link_result)
                            end if
                        else
                            if (ii .eq. size(link%P%have_output)) then
                                wrote_all_links = .true.
                            end if        
                        end if      
                        if (ii .lt. size(link%P%have_output)) then
                            wrote_all_links = .false.
                        end if
                    end do
                end if   
                call interface_write_output_line(timesecs)
            end do

            deallocate(fus_links)
            deallocate(fus_nodes)

            ! do ii = 1, size(link%P%have_output)
            !     link_idx = link%P%have_output(ii)
            !     call interface_export_link_results(link_idx)
            ! end do
        end if

        if (setting%Debug%File%output) &
        write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
        
    end subroutine output_update_swmm_out
!%
!%==========================================================================
!% END MODULE   
!%==========================================================================
!%      
end module output 

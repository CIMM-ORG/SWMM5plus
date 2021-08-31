module output

    use define_indexes
    use define_keys
    use define_globals
    use define_settings
    use define_types
    use interface
    use utility_datetime
    use, intrinsic :: iso_fortran_env, only: error_unit

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

contains

    !% subroutine for reading the link input file and storing it in link_output_idx

    subroutine output_read_csv_link_names()

        integer              :: rc, fu, pp, jj, kk, link_idx, phantom_counter
        logical              :: no_file = .false.
        character(len = 250) :: link_name
        character(64)        :: subroutine_name = 'output_read_csv_link_names'

        !%--------------------------------------------------------------------------

        if (setting%Debug%File%output) print *, '*** enter ', this_image(), subroutine_name

        if (trim(setting%Output%links_file) == "") no_file = .true.
        link_output_idx = nullvalueI
        rc = -1

        if (.not. no_file) then
            !% open csv file of link names
            open(action='read', file=trim(setting%Output%links_file), iostat=rc, newunit=fu)
            if (rc /= 0) then
                write (error_unit, '(3a, i0)') 'Opening file "', trim(setting%Output%links_file), '" failed: ', rc
                stop "in " // subroutine_name
            end if
        end if

        if (no_file) then
            !% if links_file is not specified we output all the links
            pp = 1 !% parent link
            do link_idx = 1, SWMM_N_link
                phantom_counter = 0
                link_output_idx(pp) = link_idx
                !% only parent links have associated phantoms
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
            end do
        else
            pp = 1 ! parent link
            !% loop through till the end of the file and save the valid links
            do
                !% read in the link name from the csv
                read(fu, *, iostat = rc) link_name
                if (rc /= 0) then
                    close(fu)
                    exit
                end if

                !% converting link name to link idx using the interface
                link_idx = interface_find_object(object_type=API_LINK, object_name = link_name)
                if (link_idx == 0) then
                    write(error_unit, "(A)") "Link " // trim(link_name) // " in " // &
                        trim(setting%Output%links_file) // " couldn't be found"
                    exit
                end if

                !% store index of link for output and increase index
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
            end do
        end if
        if (setting%Debug%File%output) print *, '*** leave ', this_image(),subroutine_name
    end subroutine output_read_csv_link_names

    subroutine output_read_csv_node_names()
        character(len = 250) :: node_name
        integer :: rc, fu, ii, node_idx
        character(64) :: subroutine_name = 'output_read_csv_node_names'

        !%--------------------------------------------------------------------------
        if (setting%Debug%File%output) print *, '*** enter ', this_image(), subroutine_name

        if (trim(setting%Output%nodes_file) == "") then
            !% Output all nodes if user does not specify CSV file
            node_output_idx = (/ (ii, ii =1, SWMM_N_node)/)
        else
            ii = 1
            open(action='read', file=trim(setting%Output%nodes_file), iostat=rc, newunit=fu)
            do
                read(fu, *, iostat = rc) node_name
                if (rc /= 0) then
                    node_output_idx(ii:) = nullvalueI
                    close(fu)
                    exit
                end if

                node_idx = interface_find_object(object_type=API_NODE, object_name = node_name)
                if (node_idx == 0) then
                    write(error_unit, "(A)") "Node " // trim(node_name) // " in " // &
                    trim(setting%Output%nodes_file) // " couldn't be found"
                    stop "in " // subroutine_name
                end if
                node_output_idx(ii) = node_idx
                ii = ii + 1
            end do
        end if

        print *, node_output_idx, "node_output_idx"
        stop "in " // subroutine_name
        if (setting%Debug%File%output) print *, '*** leave ', this_image(),subroutine_name

    end subroutine output_read_csv_node_names


    !% Creation of link files and header for the files
    subroutine output_create_link_files

        integer :: ii,fu, open_status, temp_link_idx
        character(len = 250) :: file_name
        character(len = 100) :: link_name
        character(len = 5)   :: str_image
        character(len = 10)  :: str_idx
        character(64) :: subroutine_name = 'output_create_link_files'

        !%--------------------------------------------------------------------------
        if (setting%Debug%File%output) print *, '*** enter ', this_image(), subroutine_name


        write(str_image, '(i5.5)') this_image()

        do ii=1, size(link%P%have_output)

            !% check if the link is a phantom link and if so find original link name and open correct file for the correct processor
            !% otherwise open file in the usual format of "link name_imageID.csv"
            if (link%P%have_output(ii) > size(link%names(:))) then
                write(str_idx, '(i5.5)') link%P%have_output(ii)
                temp_link_idx = link%P%have_output(ii)
                file_name = "debug_output/link/"//trim(link%names(link%I(temp_link_idx,li_parent_link))%str) &
                    //"_"//trim(str_image)//"_"//trim(str_idx)//".csv"
            else
                file_name = "debug_output/link/"//trim(link%names(link%P%have_output(ii))%str) &
                    //"_"//trim(str_image)//".csv"
            end if

            open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
                form   = 'formatted', action = 'write', iostat = open_status)

            if (open_status /= 0) then
                write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                stop "in " // subroutine_name
            end if

            !% Write the header of the file, set end for next write and then close file
            write(fu, *) "Timestamp,Time_In_Secs,flowrate"
            endfile(fu)
            close(fu)
        end do

        if (setting%Debug%File%output) print *, '*** leave ', this_image(),subroutine_name
    end subroutine output_create_link_files

    subroutine output_create_node_files
        integer :: ii,fu, open_status
        character(len = 250) :: file_name
        character(len = 100) :: node_name
        character(len = 5)   :: str_image
        character(64) :: subroutine_name = 'output_create_node_files'
        !%--------------------------------------------------------------------------
        if (setting%Debug%File%output) print *, '*** enter ', this_image(),subroutine_name

        !% Get current image as a string
        write(str_image, '(i5.5)') this_image()

        do ii=1, size(node%P%have_output)

            !% Open the node file
            file_name = "debug_output/node/"//trim(node%names(node%P%have_output(ii))%str) &
                //"_"//trim(str_image)//".csv"
            open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
                form   = 'formatted', action = 'write', iostat = open_status)

            if (open_status /= 0) then
                write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
                stop "in " // subroutine_name
            end if

            !% Write the header, this endfile and close the file
            write(fu, *) "Timestamp,Time_In_Secs,Head"
            endfile(fu)
            close(fu)
        end do

        if (setting%Debug%File%output) print *, '*** leave ', this_image(),subroutine_name
    end subroutine output_create_node_files


    !% This will be called at the report time step to calculate the flowrate in the link and write it to the file
    subroutine output_write_link_files

        integer :: ii, fu, open_status, yr, mnth, dy, hr, min, sec
        integer :: start_elem, end_elem, temp_link_idx
        real(8) :: time_secs, time_epoch, avg_flowrate
        character(len = 250) :: file_name
        character(len = 100) :: link_name
        character(len = 5)   :: str_image
        character(len = 10)  :: str_idx
        character(64) :: subroutine_name = 'output_write_link_files'

        !%--------------------------------------------------------------------------
        if (setting%Debug%File%output) print *, '*** enter ', this_image(),subroutine_name

        write(str_image, '(i5.5)') this_image()
        time_secs = setting%Time%Now
        time_epoch = util_datetime_secs_to_epoch(time_secs)
        call util_datetime_decodedate(time_epoch, yr, mnth, dy)
        call util_datetime_decodetime(time_epoch, hr, min, sec)

        do ii=1, size(link%P%have_output)

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
                file_name = "debug_output/link/"//trim(link%names(link%I(temp_link_idx,li_parent_link))%str) &
                    //"_"//trim(str_image)//"_"//trim(str_idx)//".csv"
            else
                file_name = "debug_output/link/"//trim(link%names(link%P%have_output(ii))%str) &
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

        end do

        if (setting%Debug%File%output) print *, '*** leave ', this_image(),subroutine_name
    end subroutine output_write_link_files

    subroutine output_write_node_files

        integer :: ii, fu, open_status, yr, mnth, dy, hr, min, sec
        integer :: temp_node_idx
        real(8) :: time_secs, time_epoch, avg_head
        character(len = 250) :: file_name
        character(len = 100) :: link_name
        character(len = 5)   :: str_image
        character(64) :: subroutine_name = 'output_write_node_files'

        !%--------------------------------------------------------------------------
        if (setting%Debug%File%output) print *, '*** enter ', this_image(),subroutine_name

        !% converter image ID to string, as well as get current time
        write(str_image, '(i5.5)') this_image()
        time_secs = setting%Time%Now
        time_epoch = util_datetime_secs_to_epoch(time_secs)
        call util_datetime_decodedate(time_epoch, yr, mnth, dy)
        call util_datetime_decodetime(time_epoch, hr, min, sec)

        !% loop through nodes we have to output
        do ii=1, size(node%P%have_output)

            !% open node file
            file_name = "debug_output/node/"//trim(node%names(node%P%have_output(ii))%str) &
                //"_"//trim(str_image)//".csv"

            open(newunit=fu, file = file_name, status = 'old',access = 'append', &
                form   = 'formatted', action = 'write', iostat = open_status)

            if (open_status /= 0) then
                write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
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
        end do

        if (setting%Debug%File%output) print *, '*** leave ', this_image(),subroutine_name

    end subroutine output_write_node_files

    subroutine output_combine_links

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

        link_output_idx_length = count(link_output_idx(:) /= nullvalueI)
        N_phantoms = sum(link%I(:, li_num_phantom_links))
        N_parents = link_output_idx_length - N_phantoms

        allocate(first_iteration(N_parents))
        allocate(full_length(N_parents))
        allocate(its_over(N_parents))
        allocate(tt(N_parents))
        allocate(file_idx(link_output_idx_length+N_parents))
        avg_flowrate => link%R(:, lr_flowrate)

        first_iteration(:) = .true.
        full_length(:) = 0
        avg_flowrate(:) = 0
        its_over(:) = .false.
        tt(:) = 1
        file_idx = 0

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
                    parent_file_name = "debug_output/link/"// &
                        trim(link%names(link%I(temp_link_idx,li_parent_link))%str) &
                            //"_"//trim(str_image)//".csv"
                    final_file_name = "swmm5_output/link/"//trim(link%names(temp_link_idx)%str)//".csv"

                    !% Open parent file
                    open(newunit=file_idx(ii), action='read', file=parent_file_name, iostat=rc)
                    if (rc /= 0) then
                        write (error_unit, '(3a, i0)') 'Opening file "', trim(parent_file_name), '" failed: ', rc
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
                        write (error_unit, '(3a, i0)') 'Opening file "', trim(Final_File_NAME), '" failed: ', open_status
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
                            phantom_file_name = "debug_output/link/"// &
                                trim(link%names(link%I(temp_phantom_link,li_parent_link))%str) &
                                //"_"//trim(str_image)//"_"//trim(str_idx)//".csv"

                            full_length(pp) = full_length(pp) + link%R(temp_phantom_link,lr_AdjustedLength)

                            open(newunit=file_idx(ii+jj), action='read', file=phantom_file_name, iostat=rc)
                            if (rc /= 0) then
                                write (error_unit, '(3a, i0)') 'Opening file "', trim(phantom_file_name), '" failed: ', rc
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

                        !% Stage entry for .out
                        call inteface_update_linkResult(pp, api_output_link_flow, real(avg_flowrate(pp),8))
                        avg_flowrate(pp) = 0
                        tt(pp) = tt(pp) + 1
                    end if
                end if
                !% Now we increment based off of how many phantom links there where related to the parent link that we combined the output
                ii = ii + link%I(temp_link_idx, li_num_phantom_links) + 1
                pp = pp + 1
            end do
            !% Write line of .out
            call interface_write_output_line(time_secs)
        end do

        do ii = 1, size(file_idx)
            close(file_idx(ii))
        end do

        deallocate(first_iteration)
        deallocate(full_length)
        deallocate(its_over)
        deallocate(file_idx)
    end subroutine output_combine_links

    subroutine output_move_node_files
        integer :: ii,fu, open_status
        character(len = 250) :: file_name, file_name_new
        character(len = 100) :: node_name
        character(len = 5)   :: str_image
        character(64) :: subroutine_name = 'output_move_node_files'
        !%--------------------------------------------------------------------------
        if (setting%Debug%File%output) print *, '*** enter ', this_image(),subroutine_name

        !% Get current image as a string
        write(str_image, '(i5.5)') this_image()

        do ii=1, size(node%P%have_output)

            !% Open the node file
            file_name = "debug_output/node/"//trim(node%names(node%P%have_output(ii))%str) &
                //"_"//trim(str_image)//".csv"
            file_name_new = "swmm5_output/node/"//trim(node%names(node%P%have_output(ii))%str)//".csv"
            call rename(file_name, file_name_new)

        end do

        if (setting%Debug%File%output) print *, '*** leave ', this_image(),subroutine_name

    end subroutine output_move_node_files
end module output


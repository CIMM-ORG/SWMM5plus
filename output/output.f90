Module output

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
    public output_create_link_files
    public output_write_link_files
    
contains
    
    !subroutine for reading the link input file and storing it in link_output_idx
    
    subroutine output_read_csv_link_names(file_name)
        
        character(len = *), intent(in) :: file_name
        character(len = 250) :: link_name
        integer :: rc, fu, ii, jj, link_temp_idx, temp_node_idx, phantom_link_idx
        integer, allocatable :: array_test(:)

        character(64) :: subroutine_name = 'output_read_csv_link_names'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%output) print *, '*** enter ', this_image(),subroutine_name

        !open csv file of link names
        open (action='read', file=file_name, iostat=rc, newunit=fu)
        ii = 1

        if(rc /= 0) then
            write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', rc
        end if

        !read the first line which is just the titles of the columns 
        read(fu, *, iostat = rc) link_name

        !loop through till the end of the file and save the valid links 
        do
            !read in the link name from the csv
            read(fu, *, iostat = rc) link_name
            if(rc /= 0) then
                exit
            end if

            !%converting link name to link idx using the interface
            link_temp_idx = interface_find_object(object_type=3, object_name = link_name)

            !% if it is an invalid link found while reading skip the loop and read the next line
            if(link_temp_idx == 0) then
                cycle 
            end if

            !store index of link for output and increase index 
            link_output_idx(ii) = link_temp_idx
            ii = ii + 1

            !checking if the link is spit across processors if so then store the id of the phantom link for output
            
            if(link%I(link_temp_idx,li_Mnode_d) >= (N_node)) then
                temp_node_idx = link%I(link_temp_idx,li_Mnode_d)
                link_output_idx(ii) = node%I(temp_node_idx,ni_Mlink_d1)
                ii = ii + 1
                
            end if
            
        end do
        N_link_output = ii - 1
        !set the rest of the array to null
        link_output_idx(ii:N_link) = nullvalueI

        if (setting%Debug%File%output) print *, '*** leave ', this_image(),subroutine_name
        
    end subroutine output_read_csv_link_names

    !Creation of link files and header for the files
    subroutine output_create_link_files

        integer :: ii,fu, open_status, temp_node_idx
        character(len = 250) :: file_name
        character(len = 100) :: link_name
        character(len = 4)   :: str_image
        character(64) :: subroutine_name = 'output_create_link_files'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%output) print *, '*** enter ', this_image(),subroutine_name
        
        
        write(str_image, '(i1)') this_image()
        
        do ii=1, size(link%P%have_output)

            !check if the link is a phantom link and if so find original link name and open correct file for the correct processor
            !otherwise open file in the usual format of "link name_imageID.csv"
            if(link%P%have_output(ii) > size(link%names(:))) then
                temp_node_idx = link%I(link%P%have_output(ii),li_Mnode_u)

                file_name = "debug_output/link/"//trim(link%names(node%I(temp_node_idx,ni_Mlink_u1))%str) &
                            //"_"//trim(str_image)//".csv"
            else
                file_name = "debug_output/link/"//trim(link%names(link%P%have_output(ii))%str)//"_"//trim(str_image)//".csv"
            end if
            
            open(newunit=fu, file = file_name, status = 'replace',access = 'sequential', &
                form   = 'formatted', action = 'write', iostat = open_status)

            if (open_status /= 0) then
                write (error_unit, '(3a, i0)') 'Opening file "', trim(FILE_NAME), '" failed: ', open_status
            end if

            !Write the header of the file, set end for next write and then close file
            
            write(fu, *) "Timestamp,Time_In_Secs,flowrate"
            endfile(fu)
            close(fu)

            
        end do

        if (setting%Debug%File%output) print *, '*** leave ', this_image(),subroutine_name
    end subroutine output_create_link_files


    
    !This will be called at the report time step to calculate the flowrate in the link and write it to the file
    subroutine output_write_link_files

        integer :: ii, fu, open_status, yr, mnth, dy, hr, min, sec
        integer :: start_elem, end_elem, temp_node_idx
        real(8) :: time_secs, time_epoch, avg_flowrate
        character(len = 250) :: file_name
        character(len = 100) :: link_name
        character(len = 4)   :: str_image
        character(64) :: subroutine_name = 'output_write_link_files'

        !--------------------------------------------------------------------------
        if (setting%Debug%File%output) print *, '*** enter ', this_image(),subroutine_name

        write(str_image, '(i1)') this_image()
        time_secs = setting%Time%Hydraulics%timeNow
        time_epoch = util_datetime_secs_to_epoch(time_secs)
        call util_datetime_decodedate(time_epoch, yr, mnth, dy)
        call util_datetime_decodetime(time_epoch, hr, min, sec)

        do ii=1, size(link%P%have_output)         

            !store the store the location of the start and end elem for easier reading
            start_elem = link%I(link%P%have_output(ii),li_first_elem_idx)
            end_elem = link%I(link%P%have_output(ii),li_last_elem_idx)

            !calculate average flowrate by summing up the elems and diving about the number of elems
            avg_flowrate = sum(elemR(start_elem:end_elem,er_Flowrate))/(end_elem-start_elem)

            !check if the link is a phantom link and if so find original link name and open correct file for the correct processor
            !otherwise open file in the usual format of "link name_imageID.csv"
            if(link%P%have_output(ii) > size(link%names(:))) then
                temp_node_idx = link%I(link%P%have_output(ii),li_Mnode_u)

                file_name = "debug_output/link/"//trim(link%names(node%I(temp_node_idx,ni_Mlink_u1))%str) &
                    //"_"//trim(str_image)//".csv"
            else
                file_name = "debug_output/link/"//trim(link%names(link%P%have_output(ii))%str)//"_"//trim(str_image)//".csv"
            end if
           
            
            open(newunit=fu, file = file_name, status = 'old',access = 'append', &
                form   = 'formatted', action = 'write', iostat = open_status)

            !writing timestamped output to file for average flowrate across the link
            
            write(fu,fmt='(i4, 2(a,i2.2))',advance = 'no') yr,"/",mnth,"/",dy
            write(fu,fmt = '(A)',advance = 'no') ' '
            write(fu,fmt='(2(i2.2,a), i2.2)',advance = 'no') hr,":",min,":",sec
            write(fu,'(A)', advance = 'no') ', '
            write(fu, '(F32.16)', advance = 'no') time_secs
            write(fu,'(A)', advance = 'no') ', '
            write(fu, '(*(G0.6 : ","))') avg_flowrate
            


            !set the end of the file for next write and close file
            endfile(fu)
            close(fu)

        end do

        if (setting%Debug%File%output) print *, '*** leave ', this_image(),subroutine_name
    end subroutine output_write_link_files

end module output

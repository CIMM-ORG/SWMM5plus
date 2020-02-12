!==========================================================================
!
! 2019-11-11 ==> contributed by Sazzad Sharior
!==========================================================================
module postProcessing


    implicit none
    private

    public  :: get_specific_link_data

    integer :: debuglevel = 0
    integer :: debuglevelall = 0
    integer :: nullvalueR = -998877.00

contains
!
!
!==========================================================================
!==========================================================================
!
!
!==========================================================================
!==========================================================================
!
subroutine get_specific_link_data &
    (iunit, n_cells, n_links, n_linkItems, max_linkItems, n_timeSteps, &
    time_steps, data_idx, length_idx, link_lengths, link_data, &
    specific_link, specific_linkData)

    character(64) :: subroutine_name = 'get_specific_link_data'

    integer, intent(in)                             :: iunit
    integer                                         :: istat
    integer, intent(inout)                          :: n_cells
    integer, intent(inout)                          :: n_links
    integer, intent(inout)                          :: max_linkItems
    integer, intent(inout)                          :: n_timeSteps
    integer, intent(in)                             :: specific_link
    integer, dimension(:), allocatable,intent(inout):: time_steps, data_idx
    integer, dimension(:), allocatable,intent(inout):: n_linkItems
    integer, dimension(:), allocatable,intent(inout):: length_idx 
    real, dimension(:,:), allocatable, intent(inout):: link_data
    real, dimension(:), allocatable, intent(inout)  :: link_lengths        
    real, dimension(:,:), allocatable, intent(out)  :: specific_linkData

    integer :: allocation_status
    character(len=99) :: emsg

    integer :: ii, jj, kk

!--------------------------------------------------------------------------
if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

    
    call get_time_steps &
        (iunit, n_cells, n_links, n_timeSteps, time_steps)
    call get_link_items &
        (iunit, n_cells, n_links, specific_link, n_linkItems, max_linkItems)
    call get_data_index &
        (iunit, max_linkItems, length_idx, data_idx, n_cells,n_timeSteps)
    call get_link_lengths &
        (iunit, n_links, max_linkItems, length_idx, specific_link, link_lengths)    
    call get_all_link_data &
        (iunit, n_cells, max_linkItems, data_idx, link_data, specific_link, n_links)
    
    ! All the link data is saved here in a single array
    allocate(specific_linkData(n_timeSteps, max_linkItems))
    specific_linkData(:,:) = nullvalueR
    
    specific_linkData(:,:) = link_data(:,:)

if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
end subroutine get_specific_link_data
!
!==========================================================================
!==========================================================================
!
! This function get the specified data of the link items
subroutine get_all_link_data &
    (iunit, n_cells, max_linkItems, data_idx, link_data, specific_link, n_links)

    character(64) :: subroutine_name = 'get_all_link_data'

    integer, intent(in)                             :: iunit, n_cells, specific_link
    integer, intent(in)                             :: max_linkItems, n_links
    integer, dimension(:), intent(in)               :: data_idx     
    integer                                         :: istat, nn
    real, dimension(:), allocatable                 :: link_data_temp
    real, dimension(:,:), allocatable, intent(out)  :: link_data

    integer :: allocation_status
    character(len=99) :: emsg

    integer :: ii, jj, kk, ll

!--------------------------------------------------------------------------
if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
    nn = n_cells/n_links

    allocate(link_data_temp (max_linkItems))
    allocate(link_data(nn, max_linkItems))

    link_data_temp(:) = nullvalueR
    link_data(:,:) = nullvalueR

    jj = 1
    kk = specific_link
    ll = 1
    rewind(iunit)
    ! This is the same algorithm as the 'get_link_lengths'
    do ii = 1, data_idx(n_cells)
        jj = data_idx(kk)
        if (ii .lt. jj) then
            read(iunit, *)
        else
            read(iunit, *)link_data_temp
            !T1 Link1: |Element1Data Element2Data ... ... ...|
            !   Link2: |Element1Data Element2Data ... ... ...|
            !   Link3: |Element1Data Element2Data ... ... ...|
            !T2 Link1: |Element1Data Element2Data ... ... ...|
            !   Link2: |Element1Data Element2Data ... ... ...|
            !   Link3: |Element1Data Element2Data ... ... ...|
            !.................................................
            link_data(ll,:) = link_data_temp
            if (kk .le. n_cells) then
                kk = kk + 3
                ll = ll + 1
                if (kk .gt. n_cells) Then
                    exit
                end if
            else
                exit
            endif
        end if
    end do
    rewind(iunit)

if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
end subroutine get_all_link_data
!
!==========================================================================
!==========================================================================
! This function get the lenghts of the link items
subroutine get_link_lengths &
    (iunit, n_links, max_linkItems, length_idx, specific_link, link_lengths)

    character(64) :: subroutine_name = 'get_link_lengths'

    integer, intent(in)                             :: iunit, n_links, specific_link
    integer, intent(in)                             :: max_linkItems
    integer, dimension(:), intent(in)               :: length_idx   
    integer                                         :: istat
    integer, dimension(:), allocatable              :: length_idx_short
    real, dimension(:), allocatable                 :: link_lengths_temp
    real, dimension(:), allocatable, intent(out)    :: link_lengths

    integer :: allocation_status
    character(len=99) :: emsg

    integer :: ii, jj, kk

!--------------------------------------------------------------------------
if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
   
    allocate(link_lengths_temp (max_linkItems))
    allocate(link_lengths(max_linkItems))
    allocate(length_idx_short(n_links))

    link_lengths_temp(:) = nullvalueR
    link_lengths(:) = nullvalueR
    ! As the lengths are repeated, only the first 'n_links' indexes are used
    length_idx_short(:) = length_idx(1:n_links)

    jj = 1
    kk = specific_link
    rewind(iunit)
    ! The first do loop is ran till the last index of saved data
    do ii = 1, length_idx_short(n_links)
        jj = length_idx_short(kk)
        if (ii .lt. jj) then
            read(iunit, *)
        else
            ! Reads through the line untill it reaches the index of 
            !saved data
            ! Then it reads the data in that index and saves in a 
            !temp array
            read(iunit, *)link_lengths
            ! Saves the link lenght in an array
            ! [Link1Lenghts ... ... ...]
            exit

        end if
    end do
    rewind(iunit)

if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
end subroutine get_link_lengths
!
!==========================================================================
!==========================================================================
!
! This function gets the index of saved link data
subroutine get_data_index &
    (iunit, max_linkItems, length_idx, data_idx, n_cells,n_timeSteps)

    character(64) :: subroutine_name = 'get_data_index'

    integer, intent(in)                             :: iunit, max_linkItems, n_timeSteps, n_cells
    integer, dimension(:), allocatable              :: n_linkItems
    character(len=512)                              :: tmp
    integer                                         :: istat
    integer, dimension(:), allocatable              :: temp_idx  
    integer, intent(out), dimension(:), allocatable :: length_idx, data_idx

    integer :: allocation_status
    character(len=99) :: emsg

    integer :: ii, jj

!--------------------------------------------------------------------------
if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
    
    allocate(temp_idx(n_cells))
    allocate(length_idx(n_cells))
    allocate(data_idx(n_cells))

    temp_idx(:) = nullvalueR
    length_idx(:) = nullvalueR
    data_idx(:) = nullvalueR

    ii = 0
    jj = 0
    rewind(iunit)
    do
        read(iunit, fmt = '(A)', iostat=istat) tmp
        if (is_iostat_end(istat)) exit
            ii = ii + 1
        ! Reads the text file line by line to find '=rows_this_link_X_data'
        ! All the data are saved after '=rows_this_link_X_data'
        if (tmp(14:36) == '=rows_this_link_X_data') then
            jj = jj + 1
        ! Find the line number of '=rows_this_link_X_data'
            temp_idx(jj) = ii
        end if
    end do
    rewind(iunit)
    ! The line after '=rows_this_link_X_data' is link lenghts. 
    length_idx = temp_idx + 1
    ! The 2nd line after '=rows_this_link_X_data' is link data. 
    data_idx = temp_idx + 2

if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
end subroutine get_data_index
!
!==========================================================================
!==========================================================================
!
! This function gets the number of elements of a link
subroutine get_link_items &
    (iunit, n_cells, n_links, specific_link, n_linkItems, max_linkItems)

    character(64) :: subroutine_name = 'get_link_items'

    integer,intent(in)                              :: iunit, specific_link
    integer, intent(in out)                         :: n_cells, n_links 
    character(len=512)                              :: tmp
    integer                                         :: istat
    integer, dimension(:), allocatable              :: n_linkItems_tmp  
    integer, intent(out), dimension(:), allocatable :: n_linkItems 
    integer, intent(out)                            :: max_linkItems

    integer :: allocation_status
    character(len=99) :: emsg

    integer :: ii

!--------------------------------------------------------------------------
if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

    n_cells = read_number_of_cells(iunit)
    n_links = read_number_of_links(iunit)

    allocate(n_linkItems_tmp(n_cells))
    allocate(n_linkItems(n_links))

    n_linkItems_tmp(:) = nullvalueR
    n_linkItems(:) = nullvalueR

    ii = 0
    rewind(iunit)
    do
        read(iunit, fmt = '(A)', iostat=istat) tmp
        if (is_iostat_end(istat)) exit
        ! Reads through the text file line by line and saves the item each
        ! link in 'tmp'
        if (tmp(14:29) == '=items_this_link') then
            ii = ii + 1
        ! Saves the tmp value in a temporary array
        ! [LinkItem1 LinkItem2 LinkItem3 LinkItem1 LinkItem2 LinkItem3 ..] 
            read(tmp(1:12), '(I12)') n_linkItems_tmp(ii)
        end if
    end do
    ! Takes the first 'n_links' values of n_linkItems_tmp array
    ! [LinkItem1 LinkItem2 LinkItem3]
    n_linkItems(:) = n_linkItems_tmp(1:n_links)
    max_linkItems = n_linkItems_tmp(specific_link)
    rewind(iunit)

if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
end subroutine get_link_items
!
!==========================================================================
!==========================================================================
!
! This function gets the time step array from the output file
subroutine get_time_steps &
    (iunit, n_cells, n_links, n_timeSteps, time_steps)

    character(64) :: subroutine_name = 'get_time_steps'

    integer,intent(in)                              :: iunit
    integer, intent(in out)                         :: n_cells, n_links 
    character(len=512)                              :: tmp
    integer                                         :: istat
    integer, dimension(:), allocatable              :: time_steps_tmp
    integer, intent(out)                            :: n_timeSteps
    integer, intent(out), dimension(:), allocatable :: time_steps

    integer :: allocation_status
    character(len=99) :: emsg

    integer :: ii, jj, kk

!--------------------------------------------------------------------------
if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
    
    n_cells = read_number_of_cells(iunit)
    n_links = read_number_of_links(iunit)

    allocate(time_steps_tmp(n_cells))
    n_timeSteps = n_cells/n_links
    allocate(time_steps(n_timeSteps))
    time_steps_tmp(:) = nullvalueR
    time_steps(:) = nullvalueR

    ii = 0
    rewind(iunit)
    ! This do loop creates a temporary array of the time steps
    ! [T1 T1 T1 T2 T2 T2 T3 T3 T3 .... ...]
    do
        read(iunit, fmt = '(A)', iostat=istat) tmp
        if (is_iostat_end(istat)) exit

        if (tmp(14:23) == '=this step') then
            ii = ii + 1
            read(tmp(1:12), '(I12)') time_steps_tmp(ii)
        end if
    end do

    ! This do loop trims the temporary array to [T1 T2 T3 ... ... ...]
    kk = 0
    do jj = 1,n_timeSteps
        kk = kk + n_links
        time_steps(jj) = time_steps_tmp(kk)
    end do
    rewind(iunit)

if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
end subroutine get_time_steps
!
!==========================================================================
!==========================================================================
!
! This function read the number of links in the output file
function read_number_of_links (iunit) result(n_links)

    character(64) :: subroutine_name = 'read_number_of_links'

    integer, intent(in)                     :: iunit   
    integer                                 :: n_links 
    character(len=512)                      :: tmp
    integer                                 :: istat

!--------------------------------------------------------------------------
if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

    rewind(iunit)
    do
        read(iunit, fmt = '(A)', iostat=istat) tmp
        if (is_iostat_end(istat)) exit
        ! Reads through the text file line by line and sets the last link 
        ! index as the n_links 
        if (tmp(14:29) == '=this_link_index') then
        read(tmp(1:12), '(I12)') n_links
        end if
    end do
    rewind(iunit)

if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
end function read_number_of_links
!
!==========================================================================
!==========================================================================
!
! This function read the number of cells in the output file
!--------------------------------------------------------------------------
!           0 =this step
!           1 =this_link_index
!          33 =items_this_link
!           2 =rows_this_link_X_data
! -107.444443      -53.7222214      -1.07444441     .....   .....
!  2.99999976       2.99999976       2.99999976     .....   ..... 
!--------------------------------------------------------------------------
! This is an example of a cell in the output file
!
function read_number_of_cells (iunit) result(n_cells)

    character(64) :: subroutine_name = 'read_number_of_cells'

    integer, intent(in)                     :: iunit   
    integer                                 :: n_cells 
    character(len=512)                      :: tmp
    integer                                 :: istat

!--------------------------------------------------------------------------
if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

    rewind(iunit)
    ! Initializing the number of cells 
    n_cells = 0
    do
        ! Reads the text file line by line and saves the line in 'tmp'
        read(iunit, fmt = '(A)', iostat=istat) tmp
        if (is_iostat_end(istat)) exit
        ! '=this step' in the text file is used as the starting point 
        ! of a cell
        if (tmp(14:23) == '=this step') then
        n_cells = n_cells + 1
        end if
    end do
    rewind(iunit)

if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
end function read_number_of_cells

end module postProcessing
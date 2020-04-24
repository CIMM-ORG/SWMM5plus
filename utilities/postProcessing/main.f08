!==========================================================================
!
! 2019-11-11 ==> contributed by Sazzad Sharior
!==========================================================================
program Plot

use postProcessing

    integer                                 :: iunit = 21, specific_linkElement, plot_type, j, k
    integer                                 :: istat
    integer                                 :: n_cells
    integer                                 :: n_links
    integer                                 :: max_linkItems
    integer                                 :: n_timeSteps
    integer                                 :: specific_link
    integer, dimension(:), allocatable      :: time_steps, data_idx, n_linkItems
    integer, dimension(:), allocatable      :: length_idx 
    real,    dimension(:,:), allocatable    :: link_data
    real,    dimension(:), allocatable      :: link_lengths
    real,    dimension(:,:), allocatable    :: link_long_data
    real,    dimension(:),   allocatable    :: link_long_lengths, z_bottoms        
    
    real,    dimension(:,:), allocatable       :: specific_linkData
    real,    dimension(:)  , allocatable       :: xx,yy


open(newunit=iunit, file='/home/saz/SWMM/SWMMengine/OutputThreaded/out_eta_20200423_1719.txt', status='OLD')
specific_link = 1
call get_specific_link_data &
    (iunit, n_cells, n_links, n_linkItems, max_linkItems, n_timeSteps, &
    time_steps, data_idx, length_idx, link_lengths, link_data, &
    specific_link, specific_linkData, link_long_data, link_long_lengths, &
    z_bottoms)



plot_type = 3
specific_linkElement = 5
! 1 = time series, 2 = longitudinal, 3 = hardcoded weir long
if (plot_type .eq. 1)  then
    xx = time_steps
    yy = specific_linkData(:,specific_linkElement)
! print*, time_steps
! print*, yy
    open (unit = 7, action = 'write', file = 'data.txt', status = 'replace')
    do i = 1,n_timeSteps
    write(7,*)xx(i), yy(i)
    end do
    call system('gnuplot -p plot.plt')
    close(7,status='delete')

elseif (plot_type .eq. 2)  then
    xx = link_lengths
    yy = specific_linkData(1,:)
    ! print*, xx
    ! print*, yy
    ! print*, size(link_lengths)
    open (unit = 7, action = 'write', file = 'data.txt', status = 'replace')
    do i = 1,size(link_lengths)
    write(7,*)xx(i), yy(i)
    ! print*, xx
    ! print*, yy
    end do
    call system('gnuplot -p plot_long.plt')
    close(7,status='delete')

elseif (plot_type .eq. 3)  then !this is hard coded for simple weir case
    do j = 1,n_timeSteps
        xx = link_long_lengths
        yy = link_long_data(j,:)
        ! yy = link_long_data(j,:) + z_bottoms !this is only for depth plot
        open (unit = 7, action = 'write', file = 'data.txt')
        do i = 1,size(link_long_lengths)
            write(7,*)xx(i), yy(i)
        end do
    end do
    call system('gnuplot -p plot_long.plt')
    close(7,status='delete')
else
    print*,'Invalid print type'
endif

end program Plot

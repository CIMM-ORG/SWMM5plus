!==========================================================================
!
! 2019-11-11 ==> contributed by Sazzad Sharior
!==========================================================================
program Plot

use postProcessing

    integer                                 :: iunit = 21, specific_linkElement
    integer                                 :: istat
    integer                                 :: n_cells
    integer                                 :: n_links
    integer                                 :: max_linkItems
    integer                                 :: n_timeSteps
    integer                                 :: specific_link
    integer, dimension(:), allocatable      :: time_steps, data_idx, n_linkItems
    integer, dimension(:), allocatable      :: length_idx 
    real,    dimension(:,:), allocatable    :: link_data
    real,    dimension(:,:), allocatable    :: link_lengths        
    
    real,    dimension(:,:), allocatable       :: specific_linkData
    real,    dimension(:)  , allocatable       :: xx,yy


open(newunit=iunit, file='/home/saz/SWMM/SWMMengine/OutputThreaded/out_eta_20200207_1540.txt', status='OLD')

specific_link = 2
call get_specific_link_data &
    (iunit, n_cells, n_links, n_linkItems, max_linkItems, n_timeSteps, &
    time_steps, data_idx, length_idx, link_lengths, link_data, &
    specific_link, specific_linkData)

specific_linkElement = 5
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


end program Plot

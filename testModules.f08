program a1

use read_width_depth

integer :: unit = 11
integer :: n_rows_in_file_node = 0
integer :: max_number_of_pairs = 0
integer, dimension(:), allocatable :: ID
 integer, dimension(:), allocatable :: numberPairs
 
 real, dimension(:), allocatable :: ManningsN
 real, dimension(:), allocatable :: Length
 real, dimension(:), allocatable :: zBottom
 real, dimension(:), allocatable :: xDistance
 real, dimension(:), allocatable :: Breadth
 
 real, dimension(:,:,:), allocatable :: widthDepthData
 
 character(len=:), allocatable :: cellType(:)

open(newunit=unit, file='WLR_WidthDepthList.txt', status='OLD')
n_rows_in_file_node = read_number_of_cells(unit)
max_number_of_pairs = read_max_number_of_pairs(unit)
print*, n_rows_in_file_node
print*, max_number_of_pairs

call read_widthdepth_pairs &
    (unit, ID, numberPairs, ManningsN, Length, zBottom, xDistance, &
     Breadth, widthDepthData, cellType)

print*, cellType


end program a1

! program a1
! 
! use read_width_depth
! 
! integer :: unit = 11
! integer :: n_rows_in_file_node = 0
! integer :: max_number_of_pairs = 0
! 
! ! real, dimension(:), allocatable :: pair
! ! character, parameter :: sep = ' '
! ! character(len=256) :: str = "2.216 0.030968"
! 
! open(newunit=unit, file='WLR_WidthDepthList.txt', status='OLD')
! n_rows_in_file_node = read_number_of_cells(unit)
! max_number_of_pairs = read_max_number_of_pairs(unit)
! print*, n_rows_in_file_node
! print*, max_number_of_pairs
! 
! call read_widthdepth_pairs (iunit)
! !call split_read_int(trim(str), sep, pair)
! ! print *,pair(1)
! ! print*, pair(2)
! 
! 
! end program a1

program a1

 use read_width_depth
 use setting_definition
 use case_waller_creek

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
 real, dimension(:), allocatable :: faceZBottom
 
 real, dimension(:,:,:), allocatable :: widthDepthData
 
 integer, dimension(:),     allocatable :: newID
 integer, dimension(:),     allocatable :: newNumberPairs
 real,    dimension(:),     allocatable :: newManningsN
 real,    dimension(:),     allocatable :: newLength
 real,    dimension(:),     allocatable :: newZBottom
 real,    dimension(:),     allocatable :: newXDistance
 real,    dimension(:),     allocatable :: newBreadth
 real,    dimension(:,:,:), allocatable :: newWidthDepthData
 
 character(len=:), allocatable :: cellType(:)

 open(newunit=unit, file='WLR_WidthDepthList.txt', status='OLD')
 n_rows_in_file_node = read_number_of_cells(unit)
 max_number_of_pairs = read_max_number_of_pairs(unit)
 print*, n_rows_in_file_node
 print*, max_number_of_pairs

 call read_widthdepth_pairs &
     (unit, ID, numberPairs, ManningsN, Length, zBottom, xDistance, &
      Breadth, widthDepthData, cellType)

 call nonmonotonic_subdivide &
      (ID, numberPairs, ManningsN, Length, zBottom, xDistance,             &
      Breadth, widthDepthData, cellType, n_rows_in_file_node, faceZBottom, &
      max_number_of_pairs, newID, newNumberPairs, newManningsN, newLength, &
      newZBottom, newXDistance, newBreadth, newWidthDepthData,             &
      subdivide_length_check)

 

end program a1

! program a1
! 
! use read_width_depth
! use setting_definition
! 
! print*, setting%Method%AdjustWidthDepth%cellSizeTarget
! 
! end program a1

! program a1
! 
! use read_width_depth
! use array_index
! use case_waller_creek
! 
!  integer :: unit = 11
!  integer :: n_rows_in_file_node = 0
!  integer :: max_number_of_pairs = 0
!  integer, dimension(:), allocatable :: ID
!  integer, dimension(:), allocatable :: numberPairs
!  
!  real, dimension(:), allocatable :: ManningsN
!  real, dimension(:), allocatable :: Length
!  real, dimension(:), allocatable :: zBottom
!  real, dimension(:), allocatable :: xDistance
!  real, dimension(:), allocatable :: Breadth
!  
!  real, target, dimension(:,:,:), allocatable :: widthDepthData
!  
!  character(len=:), allocatable :: cellType(:)
!  
!  real, dimension(:,:), allocatable :: dWidth
!  real, dimension(:,:), allocatable :: dDepth
!  
!  real, pointer :: up2W(:,:)
!  
!  integer :: nfix,ii,depth,width
! 
!  open(newunit=unit, file='WLR_WidthDepthList.txt', status='OLD')
!  n_rows_in_file_node = read_number_of_cells(unit)
!  max_number_of_pairs = read_max_number_of_pairs(unit)
!  print*, n_rows_in_file_node
!  print*, max_number_of_pairs
! 
!  call read_widthdepth_pairs &
!      (unit, ID, numberPairs, ManningsN, Length, zBottom, xDistance, &
!       Breadth, widthDepthData, cellType)
!       
!  allocate(dWidth(size(widthDepthData,1),size(widthDepthData,2)))
!  dWidth(:,:) = 0.0
!  allocate(dDepth(size(widthDepthData,1),size(widthDepthData,2)))
!  dDepth(:,:) = 0.0
! 
!  nfix = 0
!  do i=1, size(numberPairs)
!     do j=1, numberPairs(i)
!         !compute the difference in width across each level
!         dWidth(i,1:j-1) = widthDepthData(i,2:j,1) - widthDepthData(i,1:j-1,1)
!         !compute difference in depth across eacg level
!         dDepth(i,1:j-1) = widthDepthData(i,2:j,2) - widthDepthData(i,1:j-1,2)
!     enddo
!  enddo
!  
!  nfix = nfix + count(dWidth < 0.0 .or. dDepth < 0.0)
!  
!  print*, nfix
!  
!  ii=1
!  depth = 2
!  print*, widthDepthData(ii,numberPairs(ii)+1,depth)
!  
!  width = 1
!  up2W => widthDepthData(:, :, width)
!  
!  print*, up2W(2,2)
!  
!  print*, sign(-12,5)
! 
! 
! 
! end program a1

! program a1
! 
! use read_width_depth
! 
! integer :: unit = 11
! integer :: n_rows_in_file_node = 0
! integer :: max_number_of_pairs = 0
! integer, dimension(:), allocatable :: ID
!  integer, dimension(:), allocatable :: numberPairs
!  
!  real, dimension(:), allocatable :: ManningsN
!  real, dimension(:), allocatable :: Length
!  real, dimension(:), allocatable :: zBottom
!  real, dimension(:), allocatable :: xDistance
!  real, dimension(:), allocatable :: Breadth
!  
!  real, dimension(:,:,:), allocatable :: widthDepthData
!  
!  character(len=:), allocatable :: cellType(:)
!  
!  real, dimension(:), allocatable :: aaaa
! allocate(aaaa(5))
! 
! open(newunit=unit, file='WLR_WidthDepthList.txt', status='OLD')
! n_rows_in_file_node = read_number_of_cells(unit)
! max_number_of_pairs = read_max_number_of_pairs(unit)
! print*, n_rows_in_file_node
! print*, max_number_of_pairs
! 
! call read_widthdepth_pairs &
!     (unit, ID, numberPairs, ManningsN, Length, zBottom, xDistance, &
!      Breadth, widthDepthData, cellType)
! 
! !print*, cellType
! 
! aaaa(:) = 0.0
! 
! print*, size(widthDepthData,1)
! print*, size(widthDepthData,2)
! print*, size(widthDepthData,3)
! 
! aaaa(2:4) = widthDepthData(1,1:3,2) - widthDepthData(1,1:3,1)
! print*, "======"
! 
! print*, aaaa
! 
! print*, "======"
! 
! 
! print*, widthDepthData(1,1:5,1)
! 
! print*, "======"
! print*, widthDepthData(1,1:5,2)
! 
! 
! end program a1

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

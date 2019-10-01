program a1

use read_width_depth

integer :: unit = 11
integer :: n_rows_in_file_node = 0
integer :: max_number_of_pairs = 0

open(newunit=unit, file='WLR_WidthDepthList.txt', status='OLD')
n_rows_in_file_node = read_number_of_cells(unit)
max_number_of_pairs = read_max_number_of_pairs(unit)
print*, n_rows_in_file_node
print*, max_number_of_pairs


end program a1
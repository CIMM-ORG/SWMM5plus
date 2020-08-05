module tables

    use linked_list
    implicit none

    !%  diagnostic%Volume
    type table
        character(64) :: id
        integer :: curve_type ! table type
        integer, dimension(2) :: shape
        type(list) :: data
    end type table

! contains

!     subroutine add_time_series(id, type, swmm, X, Y, filepath)
!         character(64) :: id
!         integer :: type
!         logical :: swmm
!         real, dimension(*) :: X
!         real, dimension(*) :: Y
!         character(256) :: filepath

!     end subroutine add_time_series
end module tables
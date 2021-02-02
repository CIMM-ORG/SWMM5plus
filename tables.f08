module tables

    use dynamic_array

    implicit none

    integer, private :: debuglevel = 0

    ! Table types
    integer, parameter :: tseries_table = 1
    integer, parameter :: curve_table = 2

    type real_table_XY
        integer :: id
        integer :: table_type
        integer :: len
        type(real_array) :: x
        type(real_array) :: y
    end type real_table_XY

contains

    subroutine tables_add_entry(table, x, y)
        type(real_table_XY), intent(inout) :: table
        real, intent(in) :: x
        real, intent(in) :: y
        character(64) :: subroutine_name

        subroutine_name = 'tables_add_entry'

        if (debuglevel > 0) print *, '*** enter ', subroutine_name

        call dyna_real_append(table%x, x)
        call dyna_real_append(table%y, y)

        table%len = table%x%len
        if (debuglevel > 0)  print *, '*** leave ', subroutine_name
    end subroutine tables_add_entry

    subroutine free_table(table)
        type(real_table_XY), intent(inout) :: table
        character(64) :: subroutine_name

        subroutine_name = 'free_table'

        if (debuglevel > 0) print *, '*** enter ', subroutine_name

        call free_real_array(table%x)
        call free_real_array(table%y)

        if (debuglevel > 0)  print *, '*** leave ', subroutine_name
    end subroutine free_table

end module tables
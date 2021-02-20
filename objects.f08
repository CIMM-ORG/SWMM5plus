module objects

    use errors
    use dynamic_array
    use type_definitions

    implicit none

    public

    integer, private :: debuglevel = 0

    ! ----------------------------------------------
    type(real_table), allocatable :: all_tseries(:)
    type(pattern), allocatable :: all_patterns(:)
    type(totalInflow), allocatable :: total_inflows(:)
    ! ----------------------------------------------
end module objects
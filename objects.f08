module objects

    use errors
    use dynamic_array
    use type_definitions

    implicit none

    public

    integer, private :: debuglevel = 0

    ! ----------------------------------------------
    type(tseries), allocatable :: all_tseries(:)
    type(pattern), allocatable :: all_patterns(:)
    type(real_table), allocatable :: total_inflows(:)
    ! ----------------------------------------------
    type(extInflow), allocatable :: ext_inflows(:)
    type(dwfInflow), allocatable :: dwf_inflows(:)

end module objects
module tables

    use dynamic_array
    use errors
    use objects
    use globals

    implicit none

    integer, private :: debuglevel = 0

    ! Table types
    integer, parameter :: tseries_table = 1
    integer, parameter :: curve_table = 2
    integer, parameter :: tinflow = 3

    ! Interpolation Types
    integer :: INTERPOLATION_LINEAR

contains

    function new_real_table(ttype, dim)
        integer, intent(in) :: ttype
        integer, intent(in) :: dim
        type(real_table) :: new_real_table
        new_real_table%table_type = ttype
        allocate(new_real_table%tsize(dim))
        new_real_table%tsize(:) = 0
        allocate(new_real_table%data(dim))
    end function new_real_table

    subroutine tables_add_entry(table, entry, axis)
        type(real_table), intent(inout) :: table
        real, intent(in) :: entry(:)
        integer, optional, intent(in) :: axis

        integer :: n, i
        character(64) :: subroutine_name  = 'tables_add_entry'
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** enter ', subroutine_name

        n = size(entry)
        if (present(axis)) then
            if ((axis > 0) .and. (axis <= n)) then
                do i = 1, n
                    call dyna_real_append(table%data(axis), entry(i))
                end do
                table%tsize(axis) = table%tsize(axis) + n
                return
            endif
        endif

        do i = 1, n
            call dyna_real_append(table%data(i), entry(i))
            table%tsize(i) = table%tsize(i) + 1
        end do

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name
    end subroutine tables_add_entry

    subroutine free_table(table)
        type(real_table), intent(inout) :: table
        integer :: i
        character(64) :: subroutine_name

        subroutine_name = 'free_table'

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** enter ', subroutine_name
        deallocate(table%tsize)
        do i = 1, table%dim
            call free_real_array(table%data(i))
        end do
        deallocate(table%data)
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name
    end subroutine free_table

    ! Interpolation
    ! ---------------------------------------------

    function find_next_xy_between(x1, x2, y1, y2, resolution_type) result(xy)
        ! x1, x2 are in days
        real(8), intent(in) :: x1, x2, y1, y2
        integer, intent(in) :: resolution_type
        real(8), dimension(2) :: xy

        xy(1) = datetime_get_next_time(x1, resolution_type)
        if (xy(1) >= x2) then
            xy(1) = x2
            xy(2) = y2
        else
            xy(2) = interpolate(xy(1), x1, x2, y1, y2, INTERPOLATION_LINEAR)
        endif
    end function find_next_xy_between

    subroutine table_resample(tablexy, resolution_type)
        type(real_table), intent(inout) :: tablexy
        integer, intent(in) :: resolution_type
        real(8), allocatable :: x(:)
        real(8), allocatable :: y(:)
        real(8) :: x1, x2, y1, y2
        real(8) :: xy(2) = -1
        integer :: i, tsize

        tsize = tablexy%table%tsize(1)

        tablexy%table%data(1)%len = 0
        tablexy%table%data(2)%len = 0

        allocate(x(tsize))
        allocate(y(tsize))

        x(:) = tablexy%table%data(1)%array(1:tsize)
        y(:) = tablexy%table%data(2)%array(1:tsize)

        do i = 1, size(x)-1
            x1 = x(i)
            x2 = x(i+1)
            y1 = y(i)
            y2 = y(i+1)

            call tables_add_entry(tablexy%table, x(i), y(i))
            do while (xy(1) <= x(i+1))
                xy = find_next_xy_between(x1, x2, y1, y2, resolution_type)
                if (xy(1) >= x(i+1)) exit
                x1 = xy(1)
                y1 = xy(2)
                call tables_add_entry(tablexy%table, x1, y1)
            enddo
        enddo

        call tables_add_entry(tablexy%table, x(size(x)), y(size(y)))
        deallocate(x)
        deallocate(y)
    end subroutine table_resample

    function tables_find_time(table, t) result(idx)
        type(real_table), intent(in) :: table
        real(8), intent(in) :: t ! time

        integer :: idx
        integer :: n
        real :: tmp
        character(64) :: subroutine_name
        real, allocatable :: x(:)

        subroutine_name = 'tables_find_time'

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** enter ', subroutine_name
        if (table%dim .ne. 2) then
            print *, MSG_FEATURE_NOT_COMPATIBLE
            stop
        endif

        n = table%tsize(1) / 2
        x = table%data(1)%array
        ! performs binary search in tseries

        do while (.true.)
            tmp = x(n)
            if (t < tmp) then
                if (n/2 > 0) then
                    n = n / 2
                else
                    if (n == 0) then
                        idx = n
                        exit
                    else
                        if (abs(x(n) - t) < abs(x(n-1) - t)) then
                            idx = n
                            exit
                        else
                            idx = n-1
                            exit
                        end if
                    end if
                end if
            else
                if (n/2 > 0) then
                    n = n + (n/2)
                else
                    if (n == table%tsize(1)) then
                        idx = n
                        exit
                    else
                        if (abs(x(n) - t) < abs(x(n+1) - t)) then
                            idx = n
                            exit
                        else
                            idx = n+1
                            exit
                        end if
                    end if
                end if
            end if
        end do

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name

    end function tables_find_time

    function interpolate(x, x1, x2, y1, y2, interpolation_type)
        real(8), intent(in) :: x, x1, x2, y1, y2
        integer, intent(in) :: interpolation_type

        real :: interpolate

        if (interpolation_type == INTERPOLATION_LINEAR) then
            interpolate = (y2-y1)*(x-x1)/(x2-x1) + y1
        else
            print *, MSG_FEATURE_NOT_COMPATIBLE
            stop
        endif
    end function
end module tables
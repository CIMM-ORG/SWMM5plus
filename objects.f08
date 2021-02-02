module objects

    use errors
    use tables

    implicit none

    public

    type tseries
        type(real_table_XY) :: table
    end type tseries

    type pattern
        integer :: type
        integer :: count
        real, dimension(24) :: factor
    end type pattern

    integer, private :: debuglevel = 0
contains

    function find_time_tseries(val, t)
        real, intent(in) :: val
        type(tseries), intent(in) :: t
        integer :: n, find_time_tseries
        real :: tmp
        character(64) :: subroutine_name

        subroutine_name = 'find_time_tseries'

        if (debuglevel > 0) print *, '*** enter ', subroutine_name

        n = t%table%len / 2

        ! performs binary search in tseries

        do while (.true.)
            tmp = t%table%x%array(n)
            if (val < tmp) then
                if (n/2 > 0) then
                    n = n / 2
                else
                    if (n == 0) then
                        find_time_tseries = n
                        exit
                    else
                        if (abs(t%table%x%array(n) - val) < abs(t%table%x%array(n-1) - val)) then
                            find_time_tseries = n
                            exit
                        else
                            find_time_tseries = n-1
                            exit
                        end if
                    end if
                end if
            else
                if (n/2 > 0) then
                    n = n + (n/2)
                else
                    if (n == t%table%len) then
                        find_time_tseries = n
                        exit
                    else
                        if (abs(t%table%x%array(n) - val) < abs(t%table%x%array(n+1) - val)) then
                            find_time_tseries = n
                            exit
                        else
                            find_time_tseries = n+1
                            exit
                        end if
                    end if
                end if
            end if
        end do

        if (debuglevel > 0)  print *, '*** leave ', subroutine_name

    end function find_time_tseries

end module objects
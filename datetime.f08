module datetime

    implicit none

    integer :: dayspermonth(2,12) = &
        reshape((/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, & ! normal years
        31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/), (/2,12/)) ! leap years

    integer :: datedelta = 693594
    integer :: secsperday = 86400

    contains

    function isleapyear(year)
        integer, intent(in) :: year
        integer :: isleapyear
        if ((mod(year,4) == 0) .and. ((mod(year,100) .ne. 0) .or. (mod(year,400) == 0))) then
            isleapyear = 2
        else
            isleapyear = 1
        end if
    end function isleapyear

    subroutine divmod(n, d, result, remainder)
        integer, intent(in) :: n, d
        integer, intent(inout) :: result, remainder

        if (d == 0) then
            result = 0
            remainder = 0
        else
            result = n/d
            remainder = n - d*result
        endif
    end subroutine

    function datetime_encodedate(year, month, day)
        integer, intent(in) :: year, month, day
        integer :: i, j, dday
        real(4) :: datetime_encodedate

        i = isleapyear(year)
        dday = day
        if ((year >=1) .and. (year <= 9999) .and. (month >= 1) &
            .and. (month <= 12) .and. (day >= 1) .and. (day <= dayspermonth(i,month))) then
            do j = 1, month
                dday = dday + dayspermonth(i,j)
            enddo
            i = year
            datetime_encodedate = i*365 + i/4 - i/100 + i/400 + dday - datedelta
        endif
        datetime_encodedate = -datedelta
    end function datetime_encodedate

    function datetime_encodetime(hour, minute, second)
        integer, intent(in) :: hour, minute, second
        real(4) :: datetime_encodetime, s
        if ((hour >= 0) .and. (minute >= 0) .and. (second >= 0)) then
            s = (hour * 3600 + minute * 60 + second)
            datetime_encodetime = s/secsperday
        endif
        datetime_encodetime = 0
    end function datetime_encodetime

    subroutine datetime_decodedate(date, year, month, day)
        real(4), intent(in) :: date
        integer, intent(inout) :: year, month, day
        integer :: d1, d4, d100, d400
        integer :: y, m, d, i, k, t

        d1 = 365
        d4 = d1 * 4 + 1
        d100 = d4 * 25 - 1
        d400 = d100 * 4 + 1

        t = int(floor(date)) + datedelta
        if (t <= 0) then
            year = 0
            month = 1
            day = 1
        else
            t = t - 1
            y = 1
            do while (t >= d400)
                t = t - d400
                y = y + 400
            enddo
            call divmod(t, d100, i, d)
            if (i == 4) then
                i = i - 1
                d = d + d100
            endif
            y = y + i*100
            call divmod(d, d4, i, d)
            y = y + i*4
            call divmod(d, d1, i, d)
            if (i == 4) then
                i = i - 1
                d = d + d1
            endif
            y = y + i
            k = isleapyear(y)
            m = 1
            do while (.true.)
                i = dayspermonth(k, m)
                if (d < i) exit
                d = d - i
                m = m + 1
            enddo
            year = y
            month = m
            day = d + 1
        endif
    end subroutine

    subroutine datetime_decodetime(time, h, m, s)
        real(4), intent(in) :: time
        integer, intent(inout) :: h, m, s
        integer :: secs, mins
        real(4) :: fracday

        fracday = (time - floor(time)) * secsperday
        secs = int(floor(fracday + 0.5))
        if (secs >= secsperday) secs = 86399
        call divmod(secs, 60, mins, s)
        call divmod(mins, 60, h, m)
        if (h > 23) h = 0
    end subroutine datetime_decodetime

    function datetime_dayofweek(date)
        real(4), intent(in) :: date
        integer :: t, datetime_dayofweek
        t = int(floor(date)) + datedelta
        datetime_dayofweek = mod(t, 7) + 1
    end function

    function datetime_get_next_month(date1)
        real(4), intent(in) :: date1
        real(4) :: datetime_get_next_month
        real(4) :: days_til_next_day, days_til_next_month
        integer :: yy, mm, dd, i

        call datetime_decodedate(date1, yy, mm, dd)
        i = isleapyear(yy)
        days_til_next_day = 1 - floor(date1)
        days_til_next_month = dayspermonth(i, mm) - dd + days_til_next_day

        datetime_get_next_month = date1 + days_til_next_month
    end function datetime_get_next_month

    function datetime_get_next_day(date1)
        real(4), intent(in) :: date1
        real(4) :: datetime_get_next_day
        datetime_get_next_day = ceiling(date1)
    end function datetime_get_next_day

    function datetime_get_next_hour(date1)
        real(4), intent(in) :: date1
        real(4) :: datetime_get_next_hour
        datetime_get_next_hour = ceiling(date1 * 24) / 24
    end function datetime_get_next_hour

    recursive function datetime_get_next_weekendday_hour(date1) result(next)
        ! sun = 1, ..., sat = 7
        real(4), intent(in) :: date1
        real(4) :: next
        real(4) :: days_til_weekendday, d1
        integer :: dayofweek

        dayofweek = datetime_dayofweek(date1)
        if ((dayofweek > 1) .and. (dayofweek < 7)) then
            days_til_weekendday = 7 - dayofweek - floor(date1)
            next = date1 + days_til_weekendday
        else
            d1 = datetime_get_next_hour(date1)
            dayofweek = datetime_dayofweek(d1)
            if ((dayofweek .ne. 1) .and. (dayofweek .ne. 7)) then
                next = datetime_get_next_weekendday_hour(d1)
            endif
        endif
    end function datetime_get_next_weekendday_hour

end module datetime
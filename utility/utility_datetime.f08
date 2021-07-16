module utility_datetime

    use define_api_keys, only: api_daily, &
                               api_hourly, &
                               api_weekend, &
                               api_monthly
    implicit none

    integer :: dayspermonth(12,2) = &
        reshape((/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, & ! normal years
        31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/), (/12,2/)) ! leap years

    integer, parameter :: datedelta = 693594
    integer, parameter :: secsperday = 86400

    contains

    function util_datetime_get_next_time(date_in_days, resolution_type)
        real(8), intent(in) :: date_in_days
        integer, intent(in) :: resolution_type
        real(8) :: util_datetime_get_next_time

        if (resolution_type == api_daily) then
            util_datetime_get_next_time = util_datetime_get_next_day(date_in_days)
        else if (resolution_type == api_hourly) then
            util_datetime_get_next_time = util_datetime_get_next_hour(date_in_days)
        else if (resolution_type == api_monthly) then
            util_datetime_get_next_time = util_datetime_get_next_month(date_in_days)
        else if (resolution_type == api_weekend) then
            util_datetime_get_next_time = util_datetime_get_next_weekendday_hour(date_in_days)
        else
            print *, "Resolution type not supported, use"
            print *, "(1) monthly, (2) daily, (3) hourly, (4) weekend"
            stop
        endif
    end function util_datetime_get_next_time

    function util_datetime_days_to_secs(date_in_days, start_date_in_days)
        real(8), intent(in) :: date_in_days
        real(8), intent(in) :: start_date_in_days
        real(8) :: util_datetime_days_to_secs
        util_datetime_days_to_secs = (date_in_days - start_date_in_days) * real(secsperday)
    end function util_datetime_days_to_secs

    function util_datetime_secs_to_days(date_in_secs, start_date_in_days)
        real(8), intent(in) :: date_in_secs
        real(8), intent(in) :: start_date_in_days
        real(8) :: util_datetime_secs_to_days
        util_datetime_secs_to_days = date_in_secs/real(secsperday) + start_date_in_days
    end function util_datetime_secs_to_days

    function util_datetime_isleapyear(year)
        integer, intent(in) :: year
        integer :: util_datetime_isleapyear
        if ((mod(year,4) == 0) .and. ((mod(year,100) .ne. 0) .or. (mod(year,400) == 0))) then
            util_datetime_isleapyear = 2
        else
            util_datetime_isleapyear = 1
        end if
    end function util_datetime_isleapyear

    function util_datetime_encodedate(year, month, day)
        integer, intent(in) :: year, month, day
        integer :: i, j, dday
        real(8) :: util_datetime_encodedate

        i = util_datetime_isleapyear(year)
        dday = day
        if ((year >=1) .and. (year <= 9999) .and. (month >= 1) &
            .and. (month <= 12) .and. (day >= 1) .and. (day <= dayspermonth(month,i))) then
            do j = 1, month-1
                dday = dday + dayspermonth(j,i)
            enddo
            i = year-1
            util_datetime_encodedate = i*365 + i/4 - i/100 + i/400 + dday - datedelta
            return
        endif
        util_datetime_encodedate = -datedelta
    end function util_datetime_encodedate

    function util_datetime_encodetime(hour, minute, second)
        integer, intent(in) :: hour, minute, second
        real(8) :: util_datetime_encodetime, s
        if ((hour >= 0) .and. (minute >= 0) .and. (second >= 0)) then
            s = (hour * 3600 + minute * 60 + second)
            util_datetime_encodetime = s/real(secsperday)
            return
        endif
        util_datetime_encodetime = 0
    end function util_datetime_encodetime

    subroutine util_datetime_decodedate(date_in_days, year, month, day)
        !-----------------------------------------------------------------------------
        ! Description:
        !
        !
        ! Method:
        !    
        !-----------------------------------------------------------------------------
        real(8), intent(in) :: date_in_days
        integer, intent(inout) :: year, month, day
        integer :: d1, d4, d100, d400
        integer :: y, m, d, i, k, t

        d1 = 365
        d4 = d1 * 4 + 1
        d100 = d4 * 25 - 1
        d400 = d100 * 4 + 1

        t = int(floor(date_in_days)) + datedelta
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
            call util_datetime_divmod(t, d100, i, d)
            if (i == 4) then
                i = i - 1
                d = d + d100
            endif
            y = y + i*100
            call util_datetime_divmod(d, d4, i, d)
            y = y + i*4
            call util_datetime_divmod(d, d1, i, d)
            if (i == 4) then
                i = i - 1
                d = d + d1
            endif
            y = y + i
            k = util_datetime_isleapyear(y)
            m = 1
            do while (.true.)
                i = dayspermonth(m,k)
                if (d < i) exit
                d = d - i
                m = m + 1
            enddo
            year = y
            month = m
            day = d + 1
        endif
    end subroutine util_datetime_decodedate

    subroutine util_datetime_decodetime(time_in_days, h, m, s)
        !-----------------------------------------------------------------------------
        ! Description:
        !
        !
        ! Method:
        !    
        !-----------------------------------------------------------------------------
        real(8), intent(in) :: time_in_days
        integer, intent(inout) :: h, m, s
        integer :: secs, mins
        real(8) :: fracday

        fracday = (time_in_days - floor(time_in_days)) * real(secsperday)
        secs = int(floor(fracday + 0.5))
        if (secs >= real(secsperday)) secs = 86399
        call util_datetime_divmod(secs, 60, mins, s)
        call util_datetime_divmod(mins, 60, h, m)
        if (h > 23) h = 0
    end subroutine util_datetime_decodetime

    function util_datetime_dayofweek(date_in_days)
        real(8), intent(in) :: date_in_days
        integer :: t, util_datetime_dayofweek
        t = floor(date_in_days) + datedelta
        util_datetime_dayofweek = mod(t, 7)+1
    end function

    function util_datetime_get_next_month(date_in_days)
        real(8), intent(in) :: date_in_days
        real(8) :: util_datetime_get_next_month
        real(8) :: elapsed_days, days_til_next_month
        integer :: yy, mm, dd, i

        call util_datetime_decodedate(date_in_days, yy, mm, dd)
        i = util_datetime_isleapyear(yy)
        elapsed_days = real(date_in_days - int(date_in_days))
        days_til_next_month = real(dayspermonth(mm,i)) - real(dd)
        util_datetime_get_next_month = date_in_days + days_til_next_month + 1 - elapsed_days
    end function util_datetime_get_next_month

    function util_datetime_get_next_day(date_in_days)
        real(8), intent(in) :: date_in_days
        real(8) :: util_datetime_get_next_day
        if (date_in_days - int(date_in_days) > 0) then
            util_datetime_get_next_day = real(ceiling(date_in_days))
        else
            util_datetime_get_next_day = date_in_days + real(1.0)
        endif
    end function util_datetime_get_next_day

    function util_datetime_get_next_hour(date_in_days)
        real(8), intent(in) :: date_in_days
        real(8) :: util_datetime_get_next_hour
        real(8) :: n24 = 24.0
        real(8) :: n1 = 1.0
        util_datetime_get_next_hour = (int(date_in_days*n24) + n1) / n24
    end function util_datetime_get_next_hour

    recursive function util_datetime_get_next_weekendday_hour(date_in_days) result(next)
        ! sun = 1, ..., sat = 7
        real(8), intent(in) :: date_in_days
        real(8) :: next
        real(8) :: days_til_weekendday
        integer :: dayofweek, h, m, s

        dayofweek = util_datetime_dayofweek(date_in_days)
        if ((dayofweek > 1) .and. (dayofweek < 7)) then
            days_til_weekendday = real(7) - real(dayofweek) - real(date_in_days - floor(date_in_days))
            next = date_in_days + days_til_weekendday
        else
            next = util_datetime_get_next_hour(date_in_days)
            dayofweek = util_datetime_dayofweek(next)
            call util_datetime_decodetime(next, h, m, s)
            if ((dayofweek .ne. 1) .and. (dayofweek .ne. 7)) then
                if ((h == 0) .and. (m == 0) .and. (s == 0) .and. (dayofweek == 2)) return
                next = util_datetime_get_next_weekendday_hour(next)
            endif
        endif
    end function util_datetime_get_next_weekendday_hour

    subroutine util_datetime_divmod(n, d, result, remainder)
        !-----------------------------------------------------------------------------
	    ! Description:
	    !
        !
        ! Method:
        !    
        !-----------------------------------------------------------------------------
        integer, intent(in) :: n, d
        integer, intent(inout) :: result, remainder

        if (d == 0) then
            result = 0
            remainder = 0
        else
            result = n/d
            remainder = n - d*result
        endif
    end subroutine util_datetime_divmod

end module utility_datetime
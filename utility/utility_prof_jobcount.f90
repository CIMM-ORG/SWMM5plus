module utility_prof_jobcount

    use utility_prof_array
    use define_types
    type job
       type(f_array) :: time_stamps
       real :: end_time
       integer :: id = -1
    end type job

contains

    subroutine util_restart(this)
        type(job), intent(inout) :: this
        real :: start_time

        call cpu_time(start_time)
        if ( mod(this%time_stamps%len, 2) /= 0 ) then
            this%time_stamps%arr(this%time_stamps%len) = start_time
        else
            call util_prof_append(this%time_stamps, start_time)
        end if

    end subroutine util_restart

    subroutine util_finish(this)
        type(job), intent(inout) :: this
        if (mod(this%time_stamps%len, 2) /= 0) then
            call cpu_time(this%end_time)
            call util_prof_append(this%time_stamps, this%end_time)
        end if
    end subroutine util_finish

    function duration(this) result(d)
        type(job), intent(inout) :: this
        if (this%time_stamps%len >= 2) then
            d = this%time_stamps%arr(this%time_stamps%len) - this%time_stamps%arr(this%time_stamps%len-1)
        else
            d = 0
        end if
    end function duration

end module utility_prof_jobcount
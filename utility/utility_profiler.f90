module job_count

    use util_prof_array
    type job
        type(f_array) :: time_stamps
        real :: end_time
        integer :: id = -1
    end type job

contains

    subroutine restart(this)
        type(job), intent(inout) :: this
        real :: start_time

        call cpu_time(start_time)
        if ( mod(this%time_stamps%len, 2) /= 0 ) then
            this%time_stamps%arr(this%time_stamps%len) = start_time
        else
            call append(this%time_stamps, start_time)
        end if

    end subroutine restart

    subroutine finish(this)
        type(job), intent(inout) :: this
        if (mod(this%time_stamps%len, 2) /= 0) then
            call cpu_time(this%end_time)
            call append(this%time_stamps, this%end_time)
        end if
    end subroutine finish

    function duration(this) result(d)
        type(job), intent(inout) :: this
        if (this%time_stamps%len >= 2) then
            d = this%time_stamps%arr(this%time_stamps%len) - this%time_stamps%arr(this%time_stamps%len-1)
        else
            d = 0
        end if
    end function duration

end module job_count

module util_profiler
    use job_count
    implicit none

    type wall_clk
        type(job), allocatable :: jobs(:)
        integer :: max_num_jobs = 0
        integer :: num_jobs = 0
    end type wall_clk

contains

    subroutine tic(this, id)
        type(wall_clk), intent(inout) :: this
        type(job), allocatable :: resized_arr(:)
        integer :: id, i

        ! Takes care of dynamic allocation of jobs
        if (this%max_num_jobs == 0) then
            allocate(this%jobs(10))
            this%max_num_jobs = 10
            do i = 1, 10
                call init_job(this, i)
            end do
        else if (this%num_jobs == this%max_num_jobs) then
            allocate(resized_arr(this%max_num_jobs * 2))
            resized_arr(1:this%max_num_jobs) = this%jobs(1:this%max_num_jobs)
            this%max_num_jobs = this%max_num_jobs * 2
            call free_jobs(this)
            this%jobs = resized_arr
            do i = this%num_jobs, this%max_num_jobs
                call init_job(this, i)
            end do
        end if

        if (this%jobs(id)%id == -1) then
            this%jobs(id)%id = id
            this%num_jobs = this%num_jobs + 1
        end if

        call restart(this%jobs(id))
    end subroutine tic

    subroutine toc(this, id)
        type(wall_clk), intent(inout):: this
        integer, intent(in) :: id
        call finish(this%jobs(id))
    end subroutine toc

    subroutine free_jobs(this)
        type(wall_clk), intent(inout) :: this
        integer :: i
        if (this%num_jobs > 0) then
            do i = 1, this%num_jobs
                call free_arr(this%jobs(i)%time_stamps)
            end do
            deallocate(this%jobs)
        end if
    end subroutine free_jobs

    subroutine init_job(wclk, id)
        type(wall_clk), intent(inout) :: wclk
        integer, intent(in) :: id
        type(job) :: new_job
        wclk%jobs(id) = new_job
    end subroutine init_job

end module util_profiler


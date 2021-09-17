
module utility_profiler
    use utility_prof_jobcount
    use define_types
    implicit none

    type wall_clk
       type(job), allocatable :: jobs(:)
       integer :: max_num_jobs = 0
       integer :: num_jobs = 0
    end type wall_clk

contains

    subroutine util_tic(this, id)
        type(wall_clk), intent(inout) :: this
        type(job), allocatable :: resized_arr(:)
        integer :: id, i

        ! Takes care of dynamic allocation of jobs
        if (this%max_num_jobs == 0) then
            allocate(this%jobs(10))
            this%max_num_jobs = 10
            do i = 1, 10
                call util_init_job(this, i)
            end do
        else if (this%num_jobs == this%max_num_jobs) then
            allocate(resized_arr(this%max_num_jobs * 2))
            resized_arr(1:this%max_num_jobs) = this%jobs(1:this%max_num_jobs)
            this%max_num_jobs = this%max_num_jobs * 2
            call util_free_jobs(this)
            this%jobs = resized_arr
            do i = this%num_jobs, this%max_num_jobs
                call util_init_job(this, i)
            end do
        end if

        if (this%jobs(id)%id == -1) then
            this%jobs(id)%id = id
            this%num_jobs = this%num_jobs + 1
        end if
        call util_restart(this%jobs(id))
    end subroutine util_tic

    subroutine util_toc(this, id)
        type(wall_clk), intent(inout):: this
        integer, intent(in) :: id
        call util_finish(this%jobs(id))
    end subroutine util_toc

    subroutine util_free_jobs(this)
        type(wall_clk), intent(inout) :: this
        integer :: i
        if (this%num_jobs > 0) then
            do i = 1, this%num_jobs
                call util_free_arr(this%jobs(i)%time_stamps)
            end do
            deallocate(this%jobs)
        end if
    end subroutine util_free_jobs

    subroutine util_init_job(wclk, id)
        type(wall_clk), intent(inout) :: wclk
        integer, intent(in) :: id
        type(job) :: new_job
        wclk%jobs(id) = new_job
    end subroutine util_init_job

end module utility_profiler


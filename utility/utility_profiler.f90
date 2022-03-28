
module utility_profiler
    !use utility_prof_jobcount
    !use define_types
    use define_indexes
    use define_globals

    implicit none
!-----------------------------------------------------------------------------
! Description:
! Computes cumulative elapsed CPU time on each processor for procedures with
! indexes in the pf_#subroutine name# enumerated in define_indexes
!-----------------------------------------------------------------------------  
    private
    public :: util_profiler_start
    public :: util_profiler_stop
    public :: util_profiler_print_summary 

contains
    !%
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    subroutine util_profiler_start (this_procedure)
        !%--------------------------------------------------------------------------
        ! starts timing for the current processor and procedure  
        !%--------------------------------------------------------------------------
        integer, intent(in) :: this_procedure
        !%--------------------------------------------------------------------------    
        call cpu_time (profiler_data(pfr_thisstart,this_procedure))
        
    end subroutine util_profiler_start
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine util_profiler_stop (this_procedure)
        !%--------------------------------------------------------------------------
        ! stops timing for the current processor and procedure and records accumulation  
        !%--------------------------------------------------------------------------
        integer, intent(in) :: this_procedure
        !%--------------------------------------------------------------------------    
        call cpu_time (profiler_data(pfr_thisend,this_procedure))
        profiler_data(pfr_cumulative,this_procedure) =        &
              profiler_data(pfr_cumulative,this_procedure)    &
            + profiler_data(pfr_thisend   ,this_procedure)    &
            - profiler_data(pfr_thisstart ,this_procedure)

    end subroutine util_profiler_stop
    !%
    !%==========================================================================
    !%==========================================================================
    !%
    subroutine util_profiler_print_summary ()
        !%--------------------------------------------------------------------------
        !% provides command-line listing of profiling
        !%-------------------------------------------------------------------------- 
        integer :: ii, mp,kk  
        !%-------------------------------------------------------------------------- 
        do mp = 1, num_images()
            if (this_image() == mp) then
                print *, '-------------- image ',mp ,' profile data------------------'
                do kk =1,3
                    print *, 'Level ',kk,' procedures'
                    do ii =1,Ncol_pf
                        if (profiler_procedure_level(ii) == kk) then
                            print *, profiler_data(pfr_cumulative,ii),  trim(profiler_procedure_name(ii))
                        endif
                    end do
                end do
                print *, 'unclassified levels'
                do ii =1,Ncol_pf
                    if ((profiler_procedure_level(ii) < 1) .or. (profiler_procedure_level(ii) > 3)) then
                        print *, profiler_data(pfr_cumulative,ii),  &
                                 trim(profiler_procedure_name(ii)), &
                                 profiler_procedure_level(ii)
                    end if    
                end do
            end if
        end do
    end subroutine util_profiler_print_summary
    !%
    !%==========================================================================
    !%==========================================================================
    !%        
    ! subroutine util_tic(this, id)
    !     type(wall_clk), intent(inout) :: this
    !     type(job), allocatable :: resized_arr(:)
    !     integer :: id, i

    !     ! Takes care of dynamic allocation of jobs
    !     if (this%max_num_jobs == 0) then
    !         allocate(this%jobs(10))
    !         this%max_num_jobs = 10
    !         do i = 1, 10
    !             call util_init_job(this, i)
    !         end do
    !     else if (this%num_jobs == this%max_num_jobs) then
    !         allocate(resized_arr(this%max_num_jobs * 2))
    !         resized_arr(1:this%max_num_jobs) = this%jobs(1:this%max_num_jobs)
    !         this%max_num_jobs = this%max_num_jobs * 2
    !         call util_free_jobs(this)
    !         this%jobs = resized_arr
    !         do i = this%num_jobs, this%max_num_jobs
    !             call util_init_job(this, i)
    !         end do
    !     end if

    !     if (this%jobs(id)%id == -1) then
    !         this%jobs(id)%id = id
    !         this%num_jobs = this%num_jobs + 1
    !     end if
    !     call util_restart(this%jobs(id))
    ! end subroutine util_tic

    ! subroutine util_toc(this, id)
    !     type(wall_clk), intent(inout):: this
    !     integer, intent(in) :: id
    !     call util_finish(this%jobs(id))
    ! end subroutine util_toc

    ! subroutine util_free_jobs(this)
    !     type(wall_clk), intent(inout) :: this
    !     integer :: i
    !     if (this%num_jobs > 0) then
    !         do i = 1, this%num_jobs
    !             call util_free_arr(this%jobs(i)%time_stamps)
    !         end do
    !         deallocate(this%jobs)
    !     end if
    ! end subroutine util_free_jobs

    ! subroutine util_init_job(wclk, id)
    !     type(wall_clk), intent(inout) :: wclk
    !     integer, intent(in) :: id
    !     type(job) :: new_job
    !     wclk%jobs(id) = new_job
    ! end subroutine util_init_job

end module utility_profiler


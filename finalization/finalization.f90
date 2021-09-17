module finalization

    use define_settings, only: setting
    use interface
    use utility_deallocate
    use utility_profiler
    use utility_prof_jobcount
    implicit none

contains

    subroutine finalize_all()
        character(64) :: subroutine_name = 'finalize_all'
        type(wall_clk) :: timer

        real(8) :: start, intermediate, finish
        call cpu_time(start)

        if (setting%Debug%File%finalization) print *, '*** enter ', this_image(), subroutine_name
        if (setting%Profile%File%finalization) call util_tic(timer, 4)
        !------------------------------------------------------------------------------

        call interface_finalize()
        call util_deallocate_network_data()

        if (setting%Profile%File%finalization) then
            call util_toc(timer, 4)
            print *, '** time', this_image(),subroutine_name, ' = ', duration(timer%jobs(4))
            ! call util_free_jobs(timer)
        end if

    end subroutine finalize_all


end module finalization
module finalization

    use define_settings, only: setting
    use define_globals, only: timer
    use interface
    use utility_deallocate
    use utility_profiler
    use utility_prof_jobcount
    use output

    implicit none

contains

    subroutine finalize_all()
        character(64) :: subroutine_name = 'finalize_all'
        !------------------------------------------------------------------------------

        if (setting%Debug%File%finalization) print *, '*** enter ', this_image(), subroutine_name

        if (setting%Profile%File%finalization) call util_tic(timer, 4)
        if ((this_image() == 1) .and. &
            (setting%Output%report .or. setting%Debug%Output)) then
            call output_combine_links()
        end if

        sync all

        if (setting%Output%report .or. setting%Debug%Output) call output_move_node_files
        if (setting%Output%report .or. setting%Debug%Output) call output_update_swmm_out
        call interface_finalize()
        call util_deallocate_network_data()

        if (setting%Profile%File%finalization) then
            call util_toc(timer, 4)
            print *, '** time', this_image(),subroutine_name, ' = ', duration(timer%jobs(4))
            call util_free_jobs(timer)
        end if

        if ((this_image() == 1) .and. (.not. setting%Debug%Output)) then
            call system('mkdir -p debug_output')
            call system('rm -r debug_output')
        end if
    end subroutine finalize_all


end module finalization

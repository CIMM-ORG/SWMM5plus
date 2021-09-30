module finalization

    use define_settings, only: setting
    use define_globals, only: timer
    use interface
    use utility_deallocate
    use utility_profiler
    !use utility_prof_jobcount
    use output

    implicit none

contains

    subroutine finalize_all()
        character(64) :: subroutine_name = 'finalize_all'
        !------------------------------------------------------------------------------

        if (setting%Debug%File%finalization) print *, '*** enter ', this_image(), subroutine_name

        call util_profiler_print_summary()

        if ((this_image() == 1) .and. &
            (setting%Output%report .or. setting%Debug%Output)) then
            call output_combine_links()
        end if

        
        sync all

        if (setting%Output%report .or. setting%Debug%Output) call output_move_node_files
        if (setting%Output%report .or. setting%Debug%Output) call output_update_swmm_out
        call interface_finalize()
        
        call util_deallocate_network_data()

        !stop 9709
        
        !% we don't want an old debug_output directory to confuse things, so
        !% here we first create the directory (if it doesn't exist) and
        !% then remove the new directory or (recursively) remove the entire old directory
        if ((this_image() == 1) .and. (.not. setting%Debug%Output)) then
            call system('mkdir -p debug_output')
            call system('rm -r debug_output')
        end if
        
    end subroutine finalize_all


end module finalization

module finalization

    use interface
    use utility_deallocate
    use output
    use define_settings, only: setting

    implicit none

contains

    subroutine finalize_all()
        if ((this_image() == 1) .and. &
            (setting%Output%report .or. setting%Debug%Output)) then
            call output_combine_links()
        end if
        sync all
        if (setting%Output%report .or. setting%Debug%Output) call output_move_node_files
        call interface_finalize()
        call util_deallocate_network_data()

        if ((this_image() == 1) .and. (.not. setting%Debug%Output)) then
            call system('mkdir -p debug_output')
            call system('rm -r debug_output')
        end if
    end subroutine finalize_all


end module finalization

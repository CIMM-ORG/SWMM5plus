module finalization

    use define_settings, only: setting
    use define_globals, only: timer
    use interface
    use utility_deallocate
    use utility_profiler
    !use utility_prof_jobcount
    use output

    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% 
    !%-----------------------------------------------------------------------------

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine finalize_toplevel()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------
        logical :: isLastStep
        character(64) :: subroutine_name = 'finalize_toplevel'
        !------------------------------------------------------------------------------
        if (setting%Debug%File%finalization) print *, '*** enter ', this_image(), subroutine_name

        if (setting%Output%Verbose) &
            write(*,"(2A,i5,A)") new_line(" "), 'finalize [Processor ', this_image(), "]"

        !% finalize the profiler and print times
        call util_profiler_print_summary()

        !% write a final combined multi-level files
        call outputML_store_data (.true.)
        
        !% write the control file for the stored mult-level files
        call outputML_write_control_file ()

        ! !brh20211006 if ((this_image() == 1) .and. &
        ! !brh20211006     (setting%Output%report .or. setting%Debug%Output)) then
        ! !brh20211006     call outputD_combine_links()
        ! !brh20211006 end if

        sync all

        call outputML_convert_elements_to_linknode_and_write ()

        !brh20211006 if (setting%Output%report .or. setting%Debug%Output) call outputD_move_node_files
        !brh20211006 if (setting%Output%report .or. setting%Debug%Output) call outputD_update_swmm_out

        ! call interface_finalize()

        ! call util_deallocate_network_data()
        
        ! !% brh 20211005 -- commented this. since debug_output is now in a time-stamped
        ! !% subdirectory we can leave it for user to delete
        ! ! !% we don't want an old debug_output directory to confuse things, so
        ! ! !% here we first create the directory (if it doesn't exist) and
        ! ! !% then remove the new directory or (recursively) remove the entire old directory
        ! ! if ((this_image() == 1) .and. (.not. setting%Debug%Output)) then
        ! !     call system('mkdir -p debug_output')
        ! !     call system('rm -r debug_output')
        ! ! end if

        ! sync all

        ! call cpu_time(setting%Time%CPU%EpochFinishSeconds)

        write(*, "(A,i5,A,G0.6,A)") &
            new_line(" ") // 'Processor ', this_image(), " | Simulation Time = ", &
            (setting%Time%CPU%EpochFinishSeconds - setting%Time%CPU%EpochStartSeconds), " [s]"
        write(*,"(A)") '========================= SWMM5+ finished =================================='
        write(*,"(A)") ''    
        
    end subroutine finalize_toplevel
!%
!%==========================================================================
!% END MODULE
!%==========================================================================
!%
end module finalization

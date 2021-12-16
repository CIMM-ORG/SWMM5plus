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
        !%-------------------------------------------------------------------
        !% Description:
        !% Final output and cleanup
        !%-------------------------------------------------------------------
        !% Declarations:
            logical :: isLastStep
            character(64) :: subroutine_name = 'finalize_toplevel'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%finalization) print *, '*** enter ', this_image(), subroutine_name
            if (setting%Output%Verbose) &
                write(*,"(2A,i5,A)") new_line(" "), 'finalize [Processor ', this_image(), "]"
        !--------------------------------------------------------------------
        !% finalize the profiler and print times
        if (setting%Profile%useYN) call util_profiler_print_summary()

        if (setting%Simulation%useHydraulics) then !% brh20211208 -- only if N_link > 0  
            if ((setting%Output%Report%provideYN) .and. &
                (.not. setting%Output%Report%suppress_MultiLevel_Output)) then    

                !% write a final combined multi-level files
                call outputML_store_data (.true.)

                !% write the control file for the stored mult-level files
                call outputML_write_control_file ()

                sync all

                call outputML_convert_elements_to_linknode_and_write ()
            end if

        end if !% brh20211208    

        !% HACK brh20211214 -- segmentation fault from interface_finalize()
        !% when running YJ_subcritical_all_constant_inflow.inp testcase
        !call interface_finalize()

        call util_deallocate_network_data()

        sync all

        call cpu_time(setting%Time%CPU%EpochFinishSeconds)

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

module finalization

    use define_settings, only: setting
    use define_globals
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
            real(8) :: looptime
            logical :: isLastStep
            character(16) :: timevalue
            character(64) :: subroutine_name = 'finalize_toplevel'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%finalization) print *, '*** enter ', this_image(), subroutine_name
            if (setting%Output%Verbose) &
                write(*,"(2A,i5,A)") new_line(" "), 'finalize [Processor ', this_image(), "]"
        !--------------------------------------------------------------------
        !% finalize the profiler and print times
        if (setting%Profile%useYN) call util_profiler_print_summary()

        if (setting%Simulation%useHydraulics) then !% otherwise N_Link = 0 causes crash
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

        !% --- shut down EPA SWMM-C and delete the API
        print *, 'calling interface finalize'
        call interface_finalize()

        !% --- close all the allocated data
        call util_deallocate_network_data()

        sync all

        call cpu_time(setting%Time%CPU%EpochFinishSeconds)

        looptime = real(setting%Time%Real%ClockLoopEnd - setting%Time%Real%ClockLoopStart,kind=8)
        looptime = looptime / real(setting%Time%Real%ClockCountRate,8)
        if (looptime <= 90.0) then
            timevalue = 'seconds'
        else
            looptime = looptime / sixtyR
            timevalue = 'minutes'
            if (looptime > 90.0) then 
                looptime = looptime / sixtyR
                timevalue = 'hours'
                if (looptime > 48.0) then
                    looptime = looptime / twentyfourR
                    timevalue = 'days'
                end if
            end if
        end if



        !--------------------------------------------------------------------
        !% Closing
        write(*, "(A,i5,A,G0.6,A)") &
            new_line(" ") // 'Processor ', this_image(), " | Simulation Time = ", &
            (setting%Time%CPU%EpochFinishSeconds - setting%Time%CPU%EpochStartSeconds), " [s]"
        write(*,"(A,F9.2,A,A)") ' Wall-clock time spent in time-marching : ',looptime, ' ',trim(timevalue)    
        write(*,"(A)") '========================= SWMM5+ finished =================================='
        write(*,"(A)") ''

    end subroutine finalize_toplevel
!%
!%==========================================================================
!% END MODULE
!%==========================================================================
!%
end module finalization

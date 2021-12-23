module finalization

    use define_settings, only: setting
    use define_globals
    use interface
    use utility_datetime
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
            integer(kind=8) :: crate, cmax, cval
            real(8) :: total_time,timemarch_time, hydraulics_time
            real(8) :: hydrology_time, loopoutput_time, initialization_time
            real(8) :: lastoutput_time
            logical :: isLastStep
            character(8) :: total_units, timemarch_units, hydraulics_units
            character(8) :: hydrology_units, loopoutput_units, initialization_units
            character(8) :: lastoutput_units
            character(64) :: subroutine_name = 'finalize_toplevel'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%finalization) print *, '*** enter ', this_image(), subroutine_name
            if (setting%Output%Verbose) &
                write(*,"(2A,i5,A)") new_line(" "), 'finalize [Processor ', this_image(), "]"
        !--------------------------------------------------------------------
        !% finalize the profiler and print times
        if (setting%Profile%useYN) call util_profiler_print_summary()

        if (this_image()==1) then
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%FinalOutputStart = cval
        end if 

        if (setting%Simulation%useHydraulics) then !% otherwise N_Link = 0 causes crash
            if ((setting%Output%Report%provideYN) .and. &
                (.not. setting%Output%Report%suppress_MultiLevel_Output)) then    

                write(*,*) 'beginning final conversion of output (this can be slow)...'
                !% write a final combined multi-level files
                call outputML_store_data (.true.)

                !% write the control file for the stored mult-level files
                call outputML_write_control_file ()

                sync all

                call outputML_convert_elements_to_linknode_and_write ()
            end if

        end if !% brh20211208    

        !% --- shut down EPA SWMM-C and delete the API
        !print *, 'calling interface finalize'
        call interface_finalize()

        !% --- close all the allocated data
        call util_deallocate_network_data()

        sync all

        call cpu_time(setting%Time%CPU%EpochFinishSeconds)

        if (this_image()==1) then
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%End = cval

            total_time = real(setting%Time%WallClock%End - setting%Time%WallClock%Start,kind=8)
            total_time = total_time / real(setting%Time%WallClock%CountRate,kind=8)
            call util_datetime_display_time (total_time, total_units)

            !% finalize the timemarch time counter for display
            timemarch_time = real(setting%Time%WallClock%TimeMarchEnd &
                             - setting%Time%WallClock%TimeMarchStart,kind=8)
            timemarch_time = timemarch_time / real(setting%Time%WallClock%CountRate,kind=8)
            call util_datetime_display_time (timemarch_time, timemarch_units)

            !% hydraulics time
            hydraulics_time = real(setting%Time%WallClock%HydraulicsCumulative,kind=8) &
                            / real(setting%Time%WallClock%CountRate,kind=8)
            call util_datetime_display_time (hydraulics_time, hydraulics_units)

            !% hydrology time
            hydrology_time = real(setting%Time%WallClock%HydrologyCumulative,kind=8) &
                            / real(setting%Time%WallClock%CountRate,kind=8)           
            call util_datetime_display_time (hydrology_time, hydrology_units)                     
            
            !% finalize the output time counter for display
            loopoutput_time = real(setting%Time%WallClock%LoopOutputCumulative,kind=8) &
                            / real(setting%Time%WallClock%CountRate,kind=8)
            call util_datetime_display_time (loopoutput_time, loopoutput_units)   

            initialization_time = real(setting%Time%WallClock%InitializationEnd &
                                     - setting%Time%WallClock%Start,kind=8)
            initialization_time = initialization_time / real(setting%Time%WallClock%CountRate,kind=8)
            call util_datetime_display_time (initialization_time, initialization_units)          
            
            lastoutput_time = real(setting%Time%WallClock%End &
                                -  setting%Time%WallClock%FinalOutputStart)
            lastoutput_time = lastoutput_time / real(setting%Time%WallClock%CountRate,kind=8)                    
            call util_datetime_display_time (lastoutput_time, lastoutput_units)
        end if


        !--------------------------------------------------------------------
        !% Closing
        write(*, "(A,i5,A,G0.6,A)") &
            new_line(" ") // 'Processor ', this_image(), " | CPU Time = ", &
            (setting%Time%CPU%EpochFinishSeconds - setting%Time%CPU%EpochStartSeconds), " [s]"
        if (this_image() == 1) then
            write(*,*) ' '
            write(*,"(A,F9.2,A,A)") ' Wall-clock time in total                : ',total_time, ' ',trim(total_units)
            write(*,"(A,F9.2,A,A)") ' Wall-clock time spent in initialization : ',initialization_time, ' ',trim(initialization_units)
            write(*,"(A,F9.2,A,A)") ' Wall-clock time spent in time-marching  : ',timemarch_time, ' ',trim(timemarch_units)
            write(*,"(A,F9.2,A,A)") ' Wall-clock time spent in loop output    : ',loopoutput_time, ' ',trim(loopoutput_units)
            write(*,"(A,F9.2,A,A)") ' Wall-clock time spent in hydrology      : ',hydrology_time, ' ',trim(hydrology_units)
            write(*,"(A,F9.2,A,A)") ' Wall-clock time spent in hydraulics     : ',hydraulics_time, ' ',trim(hydraulics_units)
            write(*,"(A,F9.2,A,A)") ' Wall-clock time spent in final output   : ',lastoutput_time, ' ',trim(lastoutput_units)    
            write(*,*) ' '
            write(*,"(A,i9)") ' Number of serial writes to file     = ',setting%Output%NumberOfWriteSteps
            write(*,"(A,i9)") ' Total number of time levels written = ',setting%Output%NumberOfTimeLevelSaved
            write(*,"(A)") '========================= SWMM5+ finished =================================='
            write(*,"(A)") ''
        end if

    end subroutine finalize_toplevel
!%
!%==========================================================================
!% END MODULE
!%==========================================================================
!%
end module finalization

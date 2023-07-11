module finalization
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Closes and finishes a SWMM5+ run
    !%
    !%==========================================================================
    use define_settings, only: setting
    use define_globals
    use define_indexes
    use interface_
    use utility
    use utility_datetime
    use utility_deallocate
    use utility_profiler
    use output
    use utility_crash
    use utility_files, only: util_file_delete_duplicate_input

    implicit none

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
            integer :: ii
            real(8) :: total_time,timemarch_time, hydraulics_time
            real(8) :: hydrology_time, loopoutput_time, initialization_time
            real(8) :: lastoutput_time, shared_time, volume_nonconservation
            real(8) :: timemarch_seconds, shared_seconds, partition_time
            logical :: isLastStep
            character(8) :: total_units, timemarch_units, hydraulics_units
            character(8) :: hydrology_units, loopoutput_units, initialization_units
            character(8) :: lastoutput_units, shared_units, partition_units
            character(64) :: subroutine_name = 'finalize_toplevel'
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%finalization) print *, '*** enter ', this_image(), subroutine_name
        !--------------------------------------------------------------------

        !% --- finalize the profiler and print times
        if (setting%Profile%useYN) call util_profiler_print_summary()

        !% --- compute total volume non-conservation
        call util_total_volume_conservation (volume_nonconservation)

        sync all
        if (this_image()==1) then
            !% --- start the wall clock for output
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%FinalOutputStart = cval
        end if 

        !% --- combine and save all the output
        if (setting%Simulation%useHydraulics) then !% otherwise N_Link = 0 causes crash
            if ((setting%Output%Report%provideYN) .and. &
                (.not. setting%Output%Report%suppress_MultiLevel_Output)) then    

                if (this_image() == 1) write(*,"(A,i5)") '... beginning final conversion of output (this can be slow)...'
                !% --- write a final combined multi-level files
                call outputML_store_data (.true.)

                !% --- write the control file for the stored mult-level files
                call outputML_write_control_file ()

                !% --- if output static elem output is given, go through and combine those values into one output array
                if(N_Out_static_TypeElem > 0) then 
                    call outputML_combine_static_elem_data ()
                end if
                sync all
 
                call outputML_convert_elements_to_linknode_and_write ()
            end if
        end if  
        call util_crashstop(31903)

        !% --- shut down EPA SWMM-C and delete the API
        call interface_finalize()

        sync all

        !% --- delete the duplicate input files
        ! call util_file_delete_duplicate_input ()

        !% --- stop the CPU time clock
        call cpu_time(setting%Time%CPU%EpochFinishSeconds)

        if (this_image()==1) then
            !% --- stop the wall clock
            call system_clock(count=cval,count_rate=crate,count_max=cmax)
            setting%Time%WallClock%End = cval

            !% --- compute total time
            total_time = real(setting%Time%WallClock%End - setting%Time%WallClock%Start,kind=8)
            total_time = total_time / real(setting%Time%WallClock%CountRate,kind=8)
            call util_datetime_display_time (total_time, total_units)

            !% --- finalize the timemarch time counter for display
            timemarch_time = real(setting%Time%WallClock%TimeMarchEnd &
                             - setting%Time%WallClock%TimeMarchStart,kind=8)
            timemarch_time = timemarch_time / real(setting%Time%WallClock%CountRate,kind=8)
            timemarch_seconds = timemarch_time
            call util_datetime_display_time (timemarch_time, timemarch_units)

            !% --- hydraulics time
            hydraulics_time = real(setting%Time%WallClock%HydraulicsCumulative,kind=8) &
                            / real(setting%Time%WallClock%CountRate,kind=8)
            call util_datetime_display_time (hydraulics_time, hydraulics_units)

            !% --- hydrology time
            hydrology_time = real(setting%Time%WallClock%HydrologyCumulative,kind=8) &
                            / real(setting%Time%WallClock%CountRate,kind=8)           
            call util_datetime_display_time (hydrology_time, hydrology_units)      
            
            !% --- time spent in shared communication across processors
            shared_time = real(setting%Time%WallClock%SharedCumulative,kind=8) &
                            / real(setting%Time%WallClock%CountRate,kind=8)  
            shared_seconds = shared_time        
            call util_datetime_display_time (shared_time, shared_units)   
            
            !% --- output processing during time loop
            loopoutput_time = real(setting%Time%WallClock%LoopOutputCumulative,kind=8) &
                            / real(setting%Time%WallClock%CountRate,kind=8)
            call util_datetime_display_time (loopoutput_time, loopoutput_units)   

            !% --- time spent in initialization
            initialization_time = real(setting%Time%WallClock%InitializationEnd &
                                     - setting%Time%WallClock%Start,kind=8)
            initialization_time = initialization_time / real(setting%Time%WallClock%CountRate,kind=8)
            call util_datetime_display_time (initialization_time, initialization_units) 

            !% --- time spent in partitioning
            partition_time = real(setting%Time%WallClock%PartitionEnd &
                                     - setting%Time%WallClock%Start,kind=8)
                                     partition_time = partition_time / real(setting%Time%WallClock%CountRate,kind=8)
            call util_datetime_display_time (partition_time, partition_units) 
            
            !% --- time spent in the final output
            lastoutput_time = real(setting%Time%WallClock%End &
                                -  setting%Time%WallClock%FinalOutputStart)
            lastoutput_time = lastoutput_time / real(setting%Time%WallClock%CountRate,kind=8)                    
            call util_datetime_display_time (lastoutput_time, lastoutput_units)
        end if

        sync all
        !--------------------------------------------------------------------
        !% Closing
        if (setting%Output%Verbose) then
            do ii =1,num_images()
                if (ii == this_image()) then
                    if (ii==1) print *, ' '
                    write(*, "(A,i5,A,G0.6,A)") ' Processor ', this_image(), " | CPU Time = ", &
                        (setting%Time%CPU%EpochFinishSeconds - setting%Time%CPU%EpochStartSeconds), " [s]"
                end if
                sync all
            end do
            sync all

            if (this_image() == 1) then
                write(*,*) ' '
                write(*,"(A,e12.3)")' Total volume non-conservation (m^3)              = ',volume_nonconservation
                write(*,*) ' '
                write(*,"(A,i9)")  ' Number of serial writes to output files          = ',setting%Output%NumberOfWriteSteps
                write(*,"(A,i9)")  ' Total number of time levels written              = ',setting%Output%NumberOfTimeLevelSaved
                write(*,"(A,i9)")  ' Total number of finite-volume elements written   = ',sum(N_OutElem(:))
                write(*,"(A,i9)")  ' Total number of finite-volume faces written      = ',sum(N_OutFace(:))
                write(*,"(A,2i9)") ' Total number of SWMM links, nodes in system      = ',setting%SWMMinput%N_link, setting%SWMMinput%N_node
                write(*,"(A,2i9)") ' Total number of finite-volume elements in system = ',sum(N_elem)
                write(*,"(A,2i9)") ' Total number of time steps                       = ',setting%Time%step
                write(*,"(A,e12.3)") ' FLOP scale (time steps x elements)               = ',&
                    real(setting%Time%step,8) * real(sum(N_elem),8)
                write(*,"(A,e12.3)") ' FLOP scale / processor                           = ', &
                    real(setting%Time%step,8) * real(sum(N_elem),8) / real(num_images(),8)
                write(*,"(A,e12.3)") ' FLOP scale / (processor x timemarching wall-clock total time) = ', &
                    real(setting%Time%step,8) * real(sum(N_elem),8) / (real(num_images(),8) * timemarch_seconds)
                write(*,"(A,e12.3)") ' FLOP scale / (processor x (timemarch - shared) wall-clock total time) = ', &
                    real(setting%Time%step,8) * real(sum(N_elem),8) / (real(num_images(),8) * (timemarch_seconds - shared_seconds) )

                write(*,*) ' '
                write(*,"(A,F9.2,A,A)") ' Wall-clock time in total                : ',total_time, ' ',trim(total_units)
                write(*,"(A,F9.2,A,A)") ' Wall-clock time spent in initialization : ',initialization_time, ' ',trim(initialization_units)
                write(*,"(A,F9.2,A,A)") ' Wall-clock time spent in partitioning   : ',partition_time, ' ',trim(partition_units)
                write(*,"(A,F9.2,A,A)") ' Wall-clock time spent in time-marching  : ',timemarch_time, ' ',trim(timemarch_units)
                write(*,"(A,F9.2,A,A)") ' Wall-clock time spent in loop output    : ',loopoutput_time, ' ',trim(loopoutput_units)
                write(*,"(A,F9.2,A,A)") ' Wall-clock time spent in hydrology      : ',hydrology_time, ' ',trim(hydrology_units)
                write(*,"(A,F9.2,A,A)") ' Wall-clock time spent in hydraulics     : ',hydraulics_time, ' ',trim(hydraulics_units)
                write(*,"(A,F9.2,A,A)") ' Wall-clock time spent in communication  : ',shared_time, ' ',trim(shared_units)
                write(*,"(A,F9.2,A,A)") ' Wall-clock time spent in final output   : ',lastoutput_time, ' ',trim(lastoutput_units)    

                ! print *, ' '
                ! print *,  real(setting%Time%WallClock%SharedCumulative,kind=8) / real(setting%Time%WallClock%CountRate,kind=8)
                ! print *,  real(setting%Time%WallClock%SharedCumulative_A,kind=8) / real(setting%Time%WallClock%CountRate,kind=8)
                ! print *,  real(setting%Time%WallClock%SharedCumulative_B,kind=8) / real(setting%Time%WallClock%CountRate,kind=8)
                ! print *,  real(setting%Time%WallClock%SharedCumulative_C,kind=8) / real(setting%Time%WallClock%CountRate,kind=8)
            end if



        end if

        !% --- print the output folder (used for profile animation)
        if ((setting%Output%Verbose) .and. (N_profiles>zeroI)) then
            write(*,"(A)") ''
            write(*,"(A)") 'For profile_animate.py, the latest output (-o arg) is stored at... '
            write(*,"(A)") trim(setting%File%output_timestamp_subfolder)
            write(*,"(A)") ''
            write(*,"(A)") 'Copy the following for animation of first profile'
            write(*,"(A)") 'python profile_animate.py -o '//trim(setting%File%output_timestamp_subfolder)//' -s '//trim(output_profile_names(1))
        end if

        !% --- print valid profile names
        if ((setting%Output%Verbose) .and. (N_profiles>zeroI)) then
            write(*,"(A)") ''
            write(*,"(A)") 'Valid profile_animate.py names (-s arg) are...'
            do ii=1,N_profiles 
                write(*,"(A)") trim(output_profile_names(ii))
            end do
        end if

        !% --- close all the allocated data 
        !%     this is located here because it is often a seg fault when creating new code
        call util_deallocate_network_data()

        if (setting%Output%Verbose) then
            if (this_image() == 1) then
                write(*,"(A)") ''
                write(*,"(A)") '========================= SWMM5+ finished =================================='
                write(*,"(A)") ''
            end if
        end if

    end subroutine finalize_toplevel
!%
!%==========================================================================
!% END MODULE
!%==========================================================================
!%
end module finalization

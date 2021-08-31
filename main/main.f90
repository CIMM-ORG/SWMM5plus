program main

    use define_globals
    use define_indexes
    use initialization
    use interface
    use timeloop
    use finalization
    use define_settings, only: setting

    implicit none

    !real :: start_time, end_time

    ! --- store CPU clock start time
    call cpu_time(setting%Time%CPU%EpochStartSeconds)
    ! --- store Real time
    setting%Time%Real%EpochStartSeconds = time()

    ! --- Initialization
    if (setting%Verbose) &
        write(*,"(2A,i5.5,A)") new_line(" "), 'begin initialization [Image ', this_image(), "] ..."
    call initialize_all()

    ! --- Time Loop
    if (setting%Verbose) &
        write(*,"(2A,i5.5,A)") new_line(" "), 'begin timeloop [Image ', this_image(), "]"
    setting%Time%Real%EpochTimeLoopStartSeconds = time()
    call timeloop_toplevel()

    if (setting%Verbose) &
        write(*,"(2A,i5.5,A)") new_line(" "), 'finalize [Image ', this_image(), "]"
    ! --- Finalization
    call finalize_all()

    sync all

    call cpu_time(setting%Time%CPU%EpochFinishSeconds)

    print *, new_line("")
    write(*, "(A,i2,A,G0.6,A)") 'Processor ', this_image(), " | Simulation Time = ", &
        (setting%Time%CPU%EpochFinishSeconds - setting%Time%CPU%EpochStartSeconds), " [s]"

end program main

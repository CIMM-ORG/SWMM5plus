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

    if (setting%Verbose) print *, 'begin initialization...'
    ! --- Initialization
    call initialize_all()

    if (setting%Verbose) print *, 'begin timeloop'
    setting%Time%Real%EpochTimeLoopStartSeconds = time()
    ! --- Time Loop
    call timeloop_toplevel()

    if (setting%Verbose) print *, 'finalize'
    ! --- Finalization
    call finalize_all()

    sync all

    call cpu_time(setting%Time%CPU%EpochFinishSeconds)

    print *, new_line("")
    write(*, "(A,i2,A,G0.6,A)") 'Processor ', this_image(), " | Simulation Time = ", &
        (setting%Time%CPU%EpochFinishSeconds - setting%Time%CPU%EpochStartSeconds), " [s]"

end program main
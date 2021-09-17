program main

    use define_globals
    use define_indexes
    use initialization
    use interface
    use timeloop
    use finalization
    use define_settings, only: setting

    implicit none

    real :: start_time, end_time

    ! --- Clock the simulation time
    call cpu_time(start_time)
    ! --- Initialization

    call initialize_all()

    ! --- Time Loop
    call timeloop_toplevel()
    
    ! --- Finalization
    call finalize_all()

    sync all

    call cpu_time(end_time)

    print *, new_line("")
    write(*, "(A,i2,A,G0.6,A)") 'Processor ', this_image(), " | Simulation Time = ", (end_time - start_time), " [s]"

end program main
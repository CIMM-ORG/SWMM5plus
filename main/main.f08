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

    print*, 'finished simulation'

    sync all 

    call cpu_time(end_time)
    print*, 'Simulation Time = ', end_time - start_time

end program main
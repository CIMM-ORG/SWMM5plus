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
    call cpu_time(start_time)
    ! --- Initialization
    call initialize_all()

    ! --- Time Loop
    call timeloop_toplevel()

    ! --- Finalization
    call finalize_all()
    sync all 
    call cpu_time(end_time)
    print*, 'finished simulation'
    print*, 'Time = ', end_time - start_time

end program main
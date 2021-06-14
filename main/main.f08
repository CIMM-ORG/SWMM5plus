program main

    use define_globals
    use define_indexes
    use initialization
    use interface
    use timeloop
    use finalization
    use define_settings, only: setting

    implicit none

    ! --- Initialization
    call initialize_all()

    ! --- Time Loop
     call timeloop_toplevel()

    ! --- Finalization
    call finalize_all()

end program main
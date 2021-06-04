program main


    use define_globals
    use define_indexes
    use initialization
    use interface
    use finalization
    use define_settings, only: setting

    implicit none

    ! --- Initialization
    call initialize_all()

    ! --- Time Loop

    ! --- Finalization
    call finalize_all()

end program main
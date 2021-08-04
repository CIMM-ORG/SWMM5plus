module finalization

    use interface
    use utility_deallocate
    implicit none

contains

    subroutine finalize_all()
        call interface_finalize()
        call util_deallocate_network_data()
    end subroutine finalize_all


end module finalization
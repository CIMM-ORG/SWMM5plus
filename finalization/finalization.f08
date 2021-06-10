module finalization

    use interface
    use utility_deallocate
    implicit none

contains

    subroutine finalize_all()
        call interface_finalize()
    end subroutine finalize_all

    subroutine release_memory_all()
        call util_deallocate_network_data()
        call util_deallocate_api_data()
    end subroutine release_memory_all

end module finalization
module finalization

    use interface

    implicit none

contains

    subroutine finalize_all()
        call interface_finalize()
    end subroutine finalize_all

end module finalization
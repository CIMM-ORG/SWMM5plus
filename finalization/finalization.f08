module finalization

    use interface

    implicit none

contains

    subroutine finalize_all()
        call finalize_api()
    end subroutine finalize_all

end module finalization
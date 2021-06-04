module finalization

    use interface

    implicit none

contains

    subroutine finalize_all()
        call final_interface()
    end subroutine finalize_all

end module finalization
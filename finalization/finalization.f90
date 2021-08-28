module finalization

    use interface
    use utility_deallocate
    use output
    implicit none

contains

    subroutine finalize_all()
        if(this_image() == 1) then
            call output_combine_links()
        end if
        sync all 
        call interface_finalize()
        call util_deallocate_network_data()
    end subroutine finalize_all


end module finalization

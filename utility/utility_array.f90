module utility_array
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Allocation of run-time arrays. 
    !% Test module for deterimning partitioning of BIPquick with link and node
    !% arrays
    !%==========================================================================
    use utility_allocate
    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use utility_crash, only: util_crashpoint

    implicit none

    public :: util_image_number_calculation

    contains
!%
!%=========================================================================
!% PUBLIC
!%=========================================================================
!%
    subroutine util_image_number_calculation(nimgs_assign, unique_imagenum)
        !%------------------------------------------------------------------
        !% Description
        !% Get the unique list of images from partitioning modules
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(inout), allocatable :: unique_imagenum(:)
            integer, intent(inout) :: nimgs_assign
            integer, allocatable :: img_arr(:), unique(:)
            integer :: ii=0, min_val, max_val
            character(64) :: subroutine_name = 'util_image_number_calculation'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%utility_array) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------

        allocate(img_arr(size(link%I,1)))
        allocate( unique(size(link%I,1)))

        img_arr = link%I(:,li_P_image) ! original information from BIPquick -- image numbers

        min_val = minval(img_arr) - 1
        max_val = maxval(img_arr, mask=img_arr<nullValueI)

        do while (min_val .lt. max_val)
            ii = ii+1
            min_val = minval(img_arr, mask=img_arr>min_val)
            unique(ii) = min_val
        end do

        allocate(unique_imagenum(ii), source = unique(1:ii)) ! The list of image number from BIPquick

        nimgs_assign = size(unique_imagenum,1) ! The number of images assigned by BIPquick

        if ( nimgs_assign /= num_images() ) then
            write(*,"(A,i5,A)") "in subroutine " // trim(subroutine_name) // " [Processor ", this_image(), "]"
            write(*,"(A,2i5)") "There is a mismatch between the assigned images and num_images", nimgs_assign, num_images()
            call util_crashpoint(49703)
        end if

        !%-----------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%utility_array)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_image_number_calculation
!%
!%=========================================================================
!% END OF MODULE
!%=========================================================================
end module utility_array

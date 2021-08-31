! this is a test module for figuring out how to distribute the BIPquick output to Coarray
!
! BIPquick information is merged to link%I and node%I (for now)
!
!
!==========================================================================
!
module utility_array
    use utility_allocate
    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting

    implicit none

    public :: util_image_number_calculation

    contains
    !
    !==========================================================================
    ! PUBLIC
    !==========================================================================
    !
    subroutine util_image_number_calculation(nimgs_assign, unique_imagenum)
        ! Get the unique list of images from partitioning modules

        integer, intent(inout), allocatable :: unique_imagenum(:)
        integer, intent(inout) :: nimgs_assign
        integer, allocatable :: img_arr(:), unique(:)
        integer :: ii=0, min_val, max_val
        character(64) :: subroutine_name = 'util_image_number_calculation'

        if (setting%Debug%File%utility_array) print *, '*** enter ', this_image(),subroutine_name

        allocate(img_arr(size(link%I,1)))
        allocate(unique(size(link%I,1)))

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
            print*, "There is a mismatch between the assigned images and num_images", nimgs_assign, num_images()
            stop
        end if

        if (setting%Debug%File%utility_array)  print *, '*** leave ', this_image(),subroutine_name
    end subroutine util_image_number_calculation
    !
    !==========================================================================
    ! END OF MODULE
    !==========================================================================
    !
end module utility_array

! this is a test module for figuring out how to distribute the BIPquick output to Coarray
! 
! BIPquick output will be two arrays: Node(Num_N,2), Link(Num_elem,2)
!
! 
!==========================================================================
! 

module coarray_bipquick_distribute
    
    use allocate_storage
    use globals
    use array_index
    use data_keys
    use setting_definition
    
    implicit none 

    private

    public :: coarray_distribution



contains

subroutine get_image_num(nimgs)
    ! simply return the number of images
    integer, intent(inout) :: nimgs
    nimgs = num_images()
end subroutine get_image_num


subroutine calculate_array_lenth(s)

end subroutine calculate_array_lenth

subroutine coarray_distribution(link, node)

end subroutine coarray_distribution


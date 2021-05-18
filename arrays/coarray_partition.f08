! this is a test module for figuring out how to distribute the BIPquick output to Coarray
! 
! BIPquick information is merged to linkI and nodeI (for now)
!
! 
!==========================================================================
! 
module coarray_partition
    use allocate_storage
    use globals
    use assign_index
    use data_keys
    use setting_definition, only: setting
    
    implicit none

    public :: coarray_length_calculation

    contains

    subroutine image_number_calculation(nimgs_assign, unique_imagenum)
        ! Get the unique list of images from partitioning modules

        integer, intent(inout), allocatable :: unique_imagenum(:)
        integer, intent(inout) :: nimgs_assign
        integer, allocatable :: img_arr(:), unique(:)
        integer :: ii=0, min_val, max_val
        character(64) :: subroutine_name = 'image_number_calculation'

        if (setting%Debug%File%coarray_bipquick) print *, '*** enter ',subroutine_name

        allocate(img_arr(size(linkI,1)))
        allocate(unique(size(linkI,1)))

        img_arr = linkI(:,li_BQ_image) ! original information from BIPquick -- image numbers

        min_val = minval(img_arr) - 1
        max_val = maxval(img_arr)

        do while (min_val .lt. max_val)
            ii = ii+1 
            min_val = minval(img_arr, mask=img_arr>min_val)
            unique(ii) = min_val
        enddo

        allocate(unique_imagenum(ii), source = unique(1:ii)) ! The list of image number from BIPquick
        
        nimgs_assign = size(unique_imagenum,1) ! The number of images assigned by BIPquick
        
        if (setting%Debug%File%coarray_bipquick)  print *, '*** leave ',subroutine_name
    end subroutine image_number_calculation

    
    subroutine coarray_length_calculation()
        ! for coarray length determination
        integer :: nimgs_assign
        integer, allocatable :: unique_imagenum(:)
        integer :: ii, jj, kk, idx, counter, elem_counter=0, face_counter=0
        integer, allocatable :: temp_elem_N(:), temp_face_N(:)
        integer, allocatable :: node_index(:), link_index(:), temp_arr(:)
        character(64) :: subroutine_name = 'array_length_calculation'
        
        if (setting%Debug%File%coarray_bipquick) print *, '*** enter ',subroutine_name

        call image_number_calculation(nimgs_assign, unique_imagenum)

        allocate(temp_elem_N(size(unique_imagenum,1)))
        allocate(temp_face_N(size(unique_imagenum,1)))

        do ii=1, size(unique_imagenum,1)
            node_index = PACK([(counter, counter=1,size(nodeI,1))], nodeI(:, ni_BQ_image) .eq. unique_imagenum(ii))
            link_index = PACK([(counter, counter=1,size(linkI,1))], linkI(:, li_BQ_image) .eq. unique_imagenum(ii))
            ! create corresponding indices for node and link in this image
            
            ! The number of elements and faces is actually decided by the junctions
            ! So we will calculate the number of junction and base on different scenarios to decided
            ! how many elem/face are assigned to each image

            do jj = 1, size(node_index,1)
                idx = node_index(jj)
                if ( nodeI(idx,ni_node_type) .eq. nJm) then 
                    
                    elem_counter = elem_counter + J_elem_add ! 6 adjacent small elems + junction itself
                    face_counter = face_counter + J_face_add ! 6 faces for each junction
                endif
            enddo
            elem_counter = elem_counter + sum(linkI(link_index, li_N_element))
            face_counter = face_counter + sum(linkI(link_index, li_N_element)) +1

            temp_elem_N(ii) = elem_counter
            temp_face_N(ii) = face_counter

            elem_counter = 0 ! reset the counter
            face_counter = 0
        enddo        
        
        max_caf_elem_N = maxval(temp_elem_N)
        max_caf_face_N = maxval(temp_face_N) ! assign the max value 

        if (setting%Debug%File%coarray_bipquick)  print *, '*** leave ',subroutine_name

    end subroutine coarray_length_calculation



end module coarray_partition
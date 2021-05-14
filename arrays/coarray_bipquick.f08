! this is a test module for figuring out how to distribute the BIPquick output to Coarray
! 
! BIPquick information is merged to linkI and nodeI (for now)
!
! 
!==========================================================================
! 
module coarray_bipquick
    use allocate_storage
    use globals
    use array_index
    use data_keys
    use setting_definition
    
    implicit none

    integer, private :: debuglevel=0

    public :: coarray_length_calculation
    contains

    subroutine get_num_of_images()
        ! A subroutine for returing the number of involved images
        !integer, intent(inout) :: nimgs
        character(64) :: subroutine_name = 'get_num_of_images'

        if (debuglevel > 0) print *, '*** enter ',subroutine_name

        nimgs = num_images()

        if (debuglevel > 0)  print *, '*** leave ',subroutine_name
    end subroutine get_num_of_images


    subroutine image_number_calculation(nimgs_assign, unique_imagenum)
        ! Although the BIPquick output assign the image number based on numerical order, this function is another safety to make sure 
        ! we use the right number of image processors

        integer, intent(inout), allocatable :: unique_imagenum(:)
        integer, intent(inout) :: nimgs_assign
        integer, allocatable :: img_arr(:), unique(:)
        integer :: ii=0, min_val, max_val
        character(64) :: subroutine_name = 'image_number_calculation'

        if (debuglevel > 0) print *, '*** enter ',subroutine_name

        allocate(img_arr(size(linkI,1)))
        allocate(unique(size(linkI,1)))

        img_arr = linkI(:,li_BQ_image) ! original information from BIPquick - image numbers

        min_val = minval(img_arr) - 1
        max_val = maxval(img_arr)

        do while (min_val .lt. max_val)
            ii = ii+1 
            min_val = minval(img_arr, mask=img_arr>min_val)
            unique(ii) = min_val
        enddo

        allocate(unique_imagenum(ii), source = unique(1:ii)) ! The list of image number from BIPquick
        
        nimgs_assign = size(unique_imagenum,1) ! The number of images assigned by BIPquick
        
        if (debuglevel > 0)  print *, '*** leave ',subroutine_name
    end subroutine image_number_calculation

    
    subroutine coarray_length_calculation()
        ! for coarray length determination
        integer :: nimgs_assign
        integer, allocatable :: unique_imagenum(:)
        integer :: ii, jj, kk, idx, counter, elem_counter=0, face_counter=0
        integer, allocatable :: temp_elem_N(:), temp_face_N(:)
        integer, allocatable :: node_index(:), link_index(:), temp_arr(:)
        character(64) :: subroutine_name = 'array_length_calculation'
        
        if (debuglevel > 0) print *, '*** enter ',subroutine_name

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
                if ( nodeI(idx,ni_node_type) .eq. nJm) then !.and. ( nodeI(idx, ni_BQ_edge) .eq. 1 ) ) then
                    ! This needs work - we loop through all up/dn links to check if 
                    ! any of them is assigned to the same image
                    !
                    !do kk = ni_Mlink_u1, ni_Mlink_d3 ! loop through all 6 up/dn stream links
                    !    if ( any(linkI(link_index, li_idx) == nodeI(idx, kk) )) then
                    !        elem_counter = elem_counter+1
                    !        ! if any up/dn link is adjacent to the junction and in this image
                    !        ! elem + 1 (because we need to split a section to junction) 
                    !    endif
                    !enddo
                    elem_counter = elem_counter + 7 ! Junction itself
                    face_counter = face_counter + 6 ! 6 faces for each junction
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

        if (debuglevel > 0)  print *, '*** leave ',subroutine_name

    end subroutine coarray_length_calculation



!    subroutine coarray_allocation(nimgs_assign, coarray_length)
!        ! allocate the coarray to images
!
!        integer, intent(in) :: nimgs_assign, coarray_length
!        character(64) :: subroutine_name = 'coarray_allocation'
!
!        if (debuglevel > 0) print *, '*** enter ',subroutine_name
!
!        allocate(Link_c(coarray_length,2)[*]) ![1:nimgs_assign])
!        allocate(Node_c(coarray_length,2)[*]) ![1:nimgs_assign])
!        sync all 
!
!        if (debuglevel > 0)  print *, '*** leave ',subroutine_name
!
!    end subroutine coarray_allocation
!
!    subroutine coarray_data_assignment(arr)
!        integer, intent(in) :: arr(:,:)
!        integer :: nimgs = 0, nimgs_assign = 0
!        integer :: coarray_length = 0
!        integer, allocatable :: unique_imagenum(:), temp_arr(:)
!        integer :: ii
!        character(64) :: subroutine_name = 'coarray_data_assignment'
!
!        if (debuglevel > 0) print *, '*** enter ',subroutine_name
!
!
!        call get_num_of_images(nimgs) ! compute the available image number
!        call image_number_calculation(arr, nimgs_assign, unique_imagenum) ! get the images assigned by BIPquick
!
!        if (nimgs_assign .gt. nimgs) then 
!            print *, "Assigned image number is more than available image"
!            stop
!        endif
!        ! check the number of available images and make sure it is always >=1
!        ! May need refinement later
!        if ( (nimgs .lt. 1) .or. (nimgs_assign .lt. 1) ) then
!            print *, "Available number of images = ", nimgs
!            print *, "Assigned image number = ", nimgs_assign
!            stop
!        endif
!
!        ! get the max length so that we can know the memory space we should assign to coarray
!        call array_length_calculation(arr, nimgs_assign, unique_imagenum, coarray_length)
!
!        ! allocate memory
!        call coarray_allocation(nimgs_assign, coarray_length)
!
!        if (this_image() .eq. 1) then
!            do ii=1,size(unique_imagenum,1) !if BIPquick follows the numerical order (starts from 1) to assign the image number, then it should match system's default image number 
!                temp_arr = PACK(arr(:,1), arr(:,2) .eq. unique_imagenum(ii))
!                if (coarray_length .ge. size(temp_arr,1)) then
!                    Link_c(:,1)[ii] = temp_arr
!                else if (coarray_length .lt. size(temp_arr,1)) then
!                    print *, "Not enough allocated space on image =", ii, " for accommodating partitioned array from BIPquick:"
!                    print *, temp_arr
!                    stop
!                endif
!            enddo
!        endif
!        sync all
!
!        if (debuglevel > 0)  print *, '*** leave ',subroutine_name
!
!   end subroutine coarray_data_assignment

end module coarray_bipquick
! this is a test module for figuring out how to distribute the BIPquick output to Coarray
!
! BIPquick information is merged to linkI and nodeI (for now)
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

    public :: coarray_length_calculation

    contains

    subroutine image_number_calculation(nimgs_assign, unique_imagenum)
        ! Get the unique list of images from partitioning modules

        integer, intent(inout), allocatable :: unique_imagenum(:)
        integer, intent(inout) :: nimgs_assign
        integer, allocatable :: img_arr(:), unique(:)
        integer :: ii=0, min_val, max_val
        character(64) :: subroutine_name = 'image_number_calculation'

        if (setting%Debug%File%utility_array) print *, '*** enter ',subroutine_name

        allocate(img_arr(size(linkI,1)))
        allocate(unique(size(linkI,1)))

        img_arr = linkI(:,li_P_image) ! original information from BIPquick -- image numbers

        min_val = minval(img_arr) - 1
        max_val = maxval(img_arr)

        do while (min_val .lt. max_val)
            ii = ii+1
            min_val = minval(img_arr, mask=img_arr>min_val)
            unique(ii) = min_val
        enddo

        allocate(unique_imagenum(ii), source = unique(1:ii)) ! The list of image number from BIPquick

        nimgs_assign = size(unique_imagenum,1) ! The number of images assigned by BIPquick

        if (setting%Debug%File%utility_array)  print *, '*** leave ',subroutine_name
    end subroutine image_number_calculation


    subroutine coarray_length_calculation()
        ! for coarray length determination
        integer :: nimgs_assign
        integer, allocatable :: unique_imagenum(:)
        integer :: ii, jj, kk, idx, counter, elem_counter=0, face_counter=0, junction_counter=0
        integer, allocatable :: node_index(:), link_index(:), temp_arr(:)
        character(64) :: subroutine_name = 'array_length_calculation'

        if (setting%Debug%File%utility_array) print *, '*** enter ',subroutine_name

        call image_number_calculation(nimgs_assign, unique_imagenum)

        allocate(N_elem(size(unique_imagenum,1)))
        allocate(N_face(size(unique_imagenum,1)))

        do ii=1, size(unique_imagenum,1)
            node_index = PACK([(counter, counter=1,size(nodeI,1))], nodeI(:, ni_P_image) .eq. unique_imagenum(ii))
            link_index = PACK([(counter, counter=1,size(linkI,1))], linkI(:, li_P_image) .eq. unique_imagenum(ii))
            ! create corresponding indices for node and link in this image

            ! The number of elements and faces is actually decided by the junctions
            ! So we will calculate the number of junction and base on different scenarios to decided
            ! how many elem/face are assigned to each image
            junction_counter = count(nodeI(node_index, ni_node_type) == nJm)

            !% first calculate the number of nodes in each partition, assign elems/faces for junctions
            elem_counter = elem_counter + J_elem_add * junction_counter
            face_counter = face_counter + J_face_add * junction_counter

            !% loop through the links in the parition ->
            do jj = 1, size(node_index,1)
                idx = node_index(jj)
                if ( nodeI(idx, ni_node_type) .eq. nJm ) then
                    face_counter = face_counter + nodeI(idx,ni_N_link_u) ! need face for closing the upstream links
                elseif (nodeI(idx, ni_node_type) .eq. nBCdn) then
                    face_counter = face_counter + 1 !% downstream BC face
                endif

                !% if the cut at the normal junction (1up 1 down) -> make a face space
                if ( ( nodeI(idx, ni_P_is_boundary) .eq. 1 ) .and. ( nodeI(idx, ni_node_type) .ne. nJm) ) then
                    face_counter = face_counter + 1
                endif
            enddo

            !% this loop is for handling duplicated faces, we check the up/dn node of a link
            !% if the node is edge, but does not belong to this image -> make a face space for
            !% duplicating the face
            do jj = 1, size(link_index,1)
                idx = link_index(jj)
                !% check upstream node first
                if ( ( nodeI(linkI(idx, li_Mnode_u), ni_P_is_boundary) .eq. 1) .and. &
                    ( nodeI(linkI(idx, li_Mnode_u), ni_P_image) .ne. ii) ) then
                    face_counter = face_counter +1
                endif
                ! then downstream node
                if ( ( nodeI(linkI(idx, li_Mnode_d), ni_P_is_boundary) .eq. 1) .and. &
                    ( nodeI(linkI(idx, li_Mnode_d), ni_P_image) .ne. ii) ) then
                    face_counter = face_counter +1
                endif

            enddo

            elem_counter = elem_counter + sum(linkI(link_index, li_N_element))
            face_counter = face_counter + sum(linkI(link_index, li_N_element))

            N_elem(ii) = elem_counter
            N_face(ii) = face_counter

            elem_counter = 0 ! reset the counter
            face_counter = 0
        enddo

        max_caf_elem_N = maxval(N_elem)
        max_caf_face_N = maxval(N_face) ! assign the max value

        if (setting%Debug%File%utility_array) then
            do ii = 1, size(unique_imagenum,1)
                print*, 'Image => ', ii
                print*, 'Elements expected ', N_elem(ii)
                print*, 'Faces expected    ', N_face(ii)
            end do
        endif

        if (setting%Debug%File%utility_array)  print *, '*** leave ',subroutine_name

    end subroutine coarray_length_calculation



end module utility_array

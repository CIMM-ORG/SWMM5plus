module utility_allocate

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use interface

    ! use utility, only: utility_check_allocation

    implicit none

!-----------------------------------------------------------------------------
!
! Description:
!   This has all the large array allocation. The only other allocation
!   should be in boundary conditions (module bc)
!   All the top-level storage should be allocated in this module
!
!-----------------------------------------------------------------------------

    private

    ! utility_allocate constants
    integer           :: allocation_status
    character(len=99) ::              emsg

    ! public members
    public :: util_allocate_linknode
    public :: util_allocate_partitioning_arrays
    public :: util_allocate_elemX_faceX
    public :: util_allocate_columns
    public :: util_allocate_bc

contains
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_linknode()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   Allocates the link and node storage used for the coarse representation
    !   of the network connectivity
    !
    ! Method:
    !   The tables node%I, link%I, node%R, link%R, node%YN, link%YN, are allocated
    !   These are defined in globals.f08). Every time memory is allocated, the
    !   util_allocate_check functionality (from utility.f08) is used to
    !   determine wheter or not there was an error during the allocation.
    !
    !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'util_allocate_linknode'
        integer       :: additional_rows = 0
        integer       :: ii, obj_name_len

    !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% If BIPquick is being used for Partitioning, include additional rows to the link-node arrays
        if (setting%Partitioning%PartitioningMethod == BQuick) then
            additional_rows = num_images() - 1
        end if

        allocate(node%I(N_node + additional_rows, Ncol_nodeI), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        node%I(:,:) = nullvalueI

        allocate(link%I(N_link + additional_rows, Ncol_linkI), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        link%I(:,:) = nullvalueI

        allocate(node%R(N_node + additional_rows, Ncol_nodeR), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        node%R(:,:) = nullvalueR

        allocate(link%R(N_link + additional_rows, Ncol_linkR), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        link%R(:,:) = nullvalueR

        allocate(node%YN(N_node + additional_rows, Ncol_nodeYN), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        node%YN(:,:) = nullvalueL

        allocate(link%YN(N_link + additional_rows, Ncol_linkYN), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        link%YN(:,:) = nullvalueL

        !% |
        !% | Only names of objects present in EPA-SWMM are stored
        !% v

        allocate(node%Names(N_node), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)

        do ii = 1, N_node
            obj_name_len = interface_get_obj_name_len(ii, API_NODE)
            allocate(character(obj_name_len) :: node%Names(ii)%str, stat=allocation_status, errmsg=emsg)
            call util_allocate_check(allocation_status, emsg)
        end do

        allocate(link%Names(N_link), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)

        do ii = 1, N_link
            obj_name_len = interface_get_obj_name_len(ii, API_LINK)
            allocate(character(obj_name_len) :: link%Names(ii)%str, stat=allocation_status, errmsg=emsg)
            call util_allocate_check(allocation_status, emsg)
        end do

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name
    end subroutine util_allocate_linknode
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_partitioning_arrays()
        allocate(adjacent_links(max_branch_per_node))
        allocate(elem_per_image(num_images()))
        allocate(image_full(num_images()))

        !% If BIPquick is being used for Partitioning, allocate additional arrays
        if (setting%Partitioning%PartitioningMethod == BQuick) then
            allocate(B_nodeI(size(nodeI,1), max_us_branch_per_node))
            allocate(B_nodeR(size(nodeR,1), twoI))
            allocate(totalweight_visited_nodes(size(nodeI, oneI)))
            allocate(partitioned_nodes(size(nodeI, oneI)))
            allocate(partitioned_links(size(linkI, oneI)))
            allocate(weight_range(size(linkI, oneI), twoI))
            allocate(accounted_for_links(size(linkI, oneI)))
        end if
    end subroutine util_allocate_partitioning_arrays
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_deallocate_partitioning_arrays()
        deallocate(adjacent_links)
        deallocate(elem_per_image)
        deallocate(image_full)

                !% If BIPquick is being used for Partitioning, allocate additional arrays
        if (setting%Partitioning%PartitioningMethod == BQuick) then
            deallocate(B_nodeI)
            deallocate(B_nodeR)
            deallocate(totalweight_visited_nodes)
            deallocate(partitioned_nodes)
            deallocate(partitioned_links)
            deallocate(weight_range)
            deallocate(accounted_for_links)
        end if

    end subroutine util_deallocate_partitioning_arrays
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_elemX_faceX ()
        ! the max_caf_elem and max_caf_face are the maximum length of the coarray
        ! across all employed images
        ! ==========================
        ! This will be excuted at parallel level
        ! ==========================
        integer :: ii
        integer, pointer :: ncol
        character(64) :: subroutine_name = 'util_allocate_elemX_faceX'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !==== elem allocation ====
        ncol => Ncol_elemR ! the maxmiumu number of columns
        allocate(elemR(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        elemR(:,:) = nullvalueR

        ncol => Ncol_elemI
        allocate(elemI(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        elemI(:,:) = nullvalueI

        ncol => Ncol_elemYN
        allocate(elemYN(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        elemYN(:,:) = nullvalueL

        ncol => Ncol_elemP
        allocate(elemP(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        elemP(:,:) = nullvalueI

        ncol => Ncol_elemPGalltm
        allocate(elemPGalltm(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        elemPGalltm(:,:) = nullvalueI

        ncol => Ncol_elemPGetm
        allocate(elemPGetm(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        elemPGetm(:,:) = nullvalueI

        ncol => Ncol_elemPGac
        allocate(elemPGac(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        elemPGac(:,:) = nullvalueI

        ncol => Ncol_elemSI
        allocate(elemSI(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        elemSI(:,:) = nullvalueI

        ncol => Ncol_elemSR
        allocate(elemSR(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        elemSR(:,:) = nullvalueR

        ncol => Ncol_elemSGR
        allocate(elemSGR(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        elemSGR(:,:) = nullvalueR

        !==== face allocation ====
        ncol => Ncol_faceR
        allocate(faceR(max_caf_face_N, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        faceR(:,:) = nullvalueR

        ncol=> Ncol_faceI
        allocate(faceI(max_caf_face_N, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        faceI(:,:) = nullvalueI

        ncol=> Ncol_faceYN
        allocate(faceYN(max_caf_face_N, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        faceYN(:,:) = nullvalueL

        ncol=> Ncol_faceP
        allocate(faceP(max_caf_face_N, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        faceP(:,:) = nullvalueI

        allocate(facePS(max_caf_face_N, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        facePS(:,:) = nullvalueI

        ncol=> Ncol_faceM
        allocate(faceM(max_caf_face_N, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg)
        faceM(:,:) = nullvalueL

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name


    end subroutine util_allocate_elemX_faceX
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_columns()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   All the enumerated variables can not be used as pointers. Thus the
        !   variables are stored in col_elemX(:) arrays that is a target
        !
        !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'util_allocate_columns'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% allocation of the col_elemX and npack_elemX
        call util_allocate_col_elemI
        call util_allocate_col_elemP
        call util_allocate_col_elemPGalltm
        call util_allocate_col_elemPGac
        call util_allocate_col_elemPGetm
        call util_allocate_col_elemR
        call util_allocate_col_elemSI
        call util_allocate_col_elemSR
        call util_allocate_col_elemSGR
        call util_allocate_col_elemWDI
        call util_allocate_col_elemWDR
        call util_allocate_col_elemYN
        call util_allocate_col_faceI
        call util_allocate_col_faceM
        call util_allocate_col_faceP
        call util_allocate_col_facePS
        call util_allocate_col_faceR
        call util_allocate_col_faceYN

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_columns
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_elemI()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemI is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated ei_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii, jj
        character(64)       :: subroutine_name = 'util_allocate_col_elemI'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        ncol => Ncol_elemI

        !% allocate an array for storing the column
        allocate( col_elemI(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemI(:) = [(ii,ii=1,ncol)]

        !%--------------------------------------------------------------
        !% the code below is a quick print check to see if
        !% the coarray have been set up properly
        ! if (this_image() == 1) then
        !     do jj = 1, num_images()
        !         print*, jj, 'image no'
        !         print*, col_elemI(:)[jj], 'col_elemI(:)[jj]'
        !     enddo
        ! endif
        ! print*, 'press return to continue'
        ! read(*,*)
        !%--------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_elemI
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_elemP()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemP is a vector of the columns in the elemP arrays
        !   that correspond to the enumerated ep_... array_index parameter
        !
        !   the npack_elemP(:) vector contains the number of packed elements
        !   for a given column.
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemP'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemP

        !% allocate an array for storing the size of each packed type
        allocate(npack_elemP(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% allocate an array for storing the column of each packed type
        allocate( col_elemP(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemP(:) = [(ii,ii=1,ncol)]

        !% zero the number of packed items (to be defined in the packing)
        npack_elemP(:) = 0

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_elemP
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_elemPGalltm()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemPGalltm is a vector of the columns in the elemPGalltm array
        !   that correspond to the enumerated ePG_... array_index parameters
        !
        !   the npack_elemPGalltm(:) vector contains the number of packed elements
        !   for a given column.
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemPGalltm'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemPGalltm !% whatever the last item in the enumerator

        !% allocate an array for storing the size of each packed type
        allocate( npack_elemPGalltm(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% allocate an array for storing the enum type of each column
        allocate( col_elemPGalltm(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemPGalltm(:) = [(ii,ii=1,ncol)]

        !% zero the number of packed items (to be defined in the packing)
        npack_elemPGalltm(:) = 0

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_elemPGalltm
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_elemPGac()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemPGac is a vector of the columns in the elemPGac arrays
        !   that correspond to the enumerated epg_... array_index parameters
        !
        !   the npack_elemPGac(:) vector contains the number of packed elements
        !   for a given column.
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemPGac'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemPGac!% whatever the last item in the enumerator

        !% allocate an array for storing the size of each packed type
        allocate( npack_elemPGac(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% allocate an array for storing the enum type of each column
        allocate( col_elemPGac(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemPGac(:) = [(ii,ii=1,ncol)]

        !% zero the number of packed items (to be defined in the packing)
        npack_elemPGac(:) = 0

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_elemPGac
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_elemPGetm()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemPGetm is a vector of the columns in the elemPGetm arrays
        !   that correspond to the enumerated epg_... array_index parameters
        !
        !   the npack_elemPGetm(:) vector contains the number of packed elements
        !   for a given column.
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemPGetm'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemPGetm   !% whatever the last item in the enumerator

        !% allocate an array for storing the size of each packed type
        allocate( npack_elemPGetm(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% allocate an array for storing the enum type of each column
        allocate( col_elemPGetm(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemPGetm(:) = [(ii,ii=1,ncol)]

        !% zero the number of packed items (to be defined in the packing)
        npack_elemPGetm(:) = 0

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_elemPGetm
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_elemR()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemR is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated er_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemR'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemR

        !% allocate an array for storing the column
        allocate( col_elemR(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemR(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_elemR
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_elemSI()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemSI is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated eSI_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemSI'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemSI

        !% allocate an array for storing the column
        allocate( col_elemSI(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemSI(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_elemSI
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_elemSR()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemSR is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated eSR_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemSR'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemSR

        !% allocate an array for storing the column
        allocate( col_elemSR(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemSR(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_elemSR
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_elemSGR()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemSGR is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated eSGR_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemSGR'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemSGR

        !% allocate an array for storing the column
        allocate( col_elemSGR(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemSGR(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_elemSGR
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_elemWDI()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemWDI is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated eWDI_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemWDI'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemWDI

        !% allocate an array for storing the column
        allocate( col_elemWDI(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemWDI(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_elemWDI
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_elemWDR()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemWDR is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated eWDR_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemWDI'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemWDR

        !% allocate an array for storing the column
        allocate( col_elemWDR(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemWDR(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_elemWDR
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_elemYN()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemYN is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated eYN_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemYN'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_elemYN

        !% allocate an array for storing the column
        allocate( col_elemYN(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_elemYN(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_elemYN
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_faceI()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_faceI is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated fi_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_faceI'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_faceI

        !% allocate an array for storing the column
        allocate( col_faceI(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_faceI(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_faceI
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_faceM()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_faceM is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated fM_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_faceM'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_faceM

        !% allocate an array for storing the column
        allocate( col_faceM(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_faceM(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_faceM
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_faceP()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   packed arrays for faces
        !   the col_faceP is a vector of the columns in the faceP arrays
        !   that correspond to the enumerated fp_... array_index parameters
        !
        !   the npack_faceP(:) vector contains the number of packed elements
        !   for a given column.
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_faceP'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_faceP

        !% allocate an array for storing the size of each packed type
        allocate( npack_faceP(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% allocate an array for storing the column of each packed type
        allocate( col_faceP(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_faceP(:) = [(ii,ii=1,ncol)]

        !% zero the number of packed items (to be defined in the packing)
        npack_faceP(:) = 0

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_faceP
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_facePS()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   packed arrays for the shared (internal boundary)faces
        !   the col_facePS is a vector of the columns in the facePS arrays
        !   that correspond to the enumerated fp_... array_index parameters
        !   col_facePS has the same number of columns as col_faceP because
        !   all the packs for internal faces are needed for shared faces as
        !   well.
        !
        !   the npack_facePS(:) vector contains the number of packed elements
        !   for a given column.
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_facePS'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_faceP

        !% allocate an array for storing the size of each packed type
        allocate( npack_facePS(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% allocate an array for storing the column of each packed type
        allocate( col_facePS(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_facePS(:) = [(ii,ii=1,ncol)]

        !% zero the number of packed items (to be defined in the packing)
        npack_facePS(:) = 0

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_facePS
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_faceR()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_faceR is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated fr_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_faceR'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_faceR  !global

        !% allocate an array of column indexes that can be used as targets of pointers
        allocate( col_faceR(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_faceR(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_faceR
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_col_faceYN()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_faceYN is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated fYN_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_faceYN'

        !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        !% define the maximum number of columns as
        ncol => Ncol_faceYN  !global

        !% allocate an array of column indexes that can be used as targets of pointers
        allocate( col_faceYN(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg)

        !% this array can be used as a pointer target in defining masks
        col_faceYN(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_col_faceYN
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_bc()
    !
    ! allocate storage for boundary conditions.
    !
    !-----------------------------------------------------------------------------

        character(64)      :: subroutine_name = 'util_allocate_bc'
        integer            :: ii, allocation_status, bc_node
        character(len=99)  :: emsg

    !-----------------------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) print *, '*** enter ',subroutine_name

        if (setting%BC%BCSlots < 2) then
            print *, "Error: the number of BCSlots has to be greater than 2"
            stop
        end if

        if (N_headBC > 0) then
            allocate(BC%headI(N_headBC, N_headI), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg)

            allocate(BC%headR_timeseries(N_headBC, setting%BC%BCSlots, N_headR), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg)

            allocate(BC%headIdx(N_headBC), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg)

            allocate(BC%headRI(N_headBC), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg)
        end if

        if (N_flowBC > 0) then
            allocate(BC%flowI(N_flowBC, N_flowI), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg)

            allocate(BC%flowR_timeseries(N_flowBC, setting%BC%BCSlots, N_flowR), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg)

            allocate(BC%flowIdx(N_flowBC), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg)

            allocate(BC%flowRI(N_flowBC), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg)
        end if

        if (setting%Debug%File%utility_allocate) print *, '*** leave ',subroutine_name
    end subroutine util_allocate_bc
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine util_allocate_check(allocation_status, emsg)
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   Checks allocation status and stops if there is an error
        !
        !-----------------------------------------------------------------------------

            integer,           intent(in   ) :: allocation_status
            character(len=*),  intent(in   ) :: emsg

            character(64):: subroutine_name = 'util_allocate_check'

        !-----------------------------------------------------------------------------

            if (setting%Debug%File%utility) print *, '*** enter ',subroutine_name

            if (allocation_status > 0) then
                print *, trim(emsg)
                stop
            end if

            if (setting%Debug%File%utility) print *, '*** leave ',subroutine_name

    end subroutine util_allocate_check

end module utility_allocate

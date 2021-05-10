module allocate_storage

    use array_index
    use globals
    use utility, only: utility_check_allocation

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

    ! allocate_storage constants
    integer           :: allocation_status
    integer           ::        debuglevel = 0
    character(len=99) ::              emsg

    ! public members
    public :: allocate_linknode_storage

contains

    subroutine allocate_linknode_storage ()

    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   Allocates the link and node storage used for the coarse representation
    !   of the network connectivity
    !
    ! Method:
    !   The tables nodeI, linkI, nodeR, linkR, nodeYN, linkYN, are allocated
    !   These are defined in globals.f08). Every time memory is allocated, the
    !   utility_check_allocation functionality (from utility.f08) is used to
    !   determine wheter or not there was an error during the allocation.
    !
    !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'allocate_linknode_storage'

    !-----------------------------------------------------------------------------

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        allocate(nodeI(N_node, ni_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        nodeI(:,:) = nullvalueI

        allocate(linkI(N_link, li_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        linkI(:,:) = nullvalueI

        allocate(B_nodeI(N_node, 3), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        B_nodeI(:,:) = nullvalueI

        allocate(B_linkI(N_link, 2), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        B_linkI(:,:) = nullvalueI

        allocate(nodeR(N_node, nr_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        nodeR(:,:) = nullvalueR

        allocate(linkR(N_link, lr_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        linkR(:,:) = nullvalueR

        allocate(nodeYN(N_node, nYN_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        nodeYN(:,:) = nullvalueL

        allocate(linkYN(N_link, lYN_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        linkYN(:,:) = nullvalueL

        allocate(nodeName(N_node), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)

        allocate(linkName(N_link), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
    end subroutine allocate_linknode_storage


    subroutine coarray_storage_allocation(max_caf_elem_N, max_caf_face_N)
        ! the max_caf_elem and max_caf_face are the maximum length of the coarray 
        ! across all employed images
        integer, intent(in) :: max_caf_elem_N, max_caf_face_N
        integer, pointer :: ncol
        !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'coarray_storage_allocation'
    
        !-----------------------------------------------------------------------------

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


        ncol => Ncol_elemR ! the maxmiumu number of columns
        allocate(elemR_caf(max_caf_elem_N, ncol)[*], stat=allocation_status, errmsg=emsg)
        
        ncol => Ncol_elemI
        allocate(elemI_caf(max_caf_elem_N, ncol)[*], stat=allocation_status, errmsg=emsg)

        ncol => Ncol_elemYN
    end subroutine coarray_storage_allocation
end module allocate_storage

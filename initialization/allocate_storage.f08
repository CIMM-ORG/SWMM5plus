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


    subroutine coarray_storage_allocation()
        ! the max_caf_elem and max_caf_face are the maximum length of the coarray 
        ! across all employed images
        ! ==========================
        ! This will be excuted at parallel level
        ! ==========================

        character(64) :: subroutine_name = 'coarray_storage_allocation'
    
        !-----------------------------------------------------------------------------

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !==== elem allocation ====
        !ncol => Ncol_elemR ! the maxmiumu number of columns
        allocate(elemR_caf(max_caf_elem_N, Ncol_elemR)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemR_caf(:,:) = nullvalueR

        !ncol => Ncol_elemI
        allocate(elemI_caf(max_caf_elem_N, Ncol_elemI)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemI_caf(:,:) = nullvalueI

        !ncol => Ncol_elemYN
        allocate(elemYN_caf(max_caf_elem_N, Ncol_elemYN)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemYN_caf(:,:) = nullvalueL

        !ncol => Ncol_elemP
        allocate(elemP_caf(max_caf_elem_N, Ncol_elemP)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemP_caf(:,:) = nullvalueI

        ! Ncol_elemPG not defined yet. Comment out
        !ncol => Ncol_elemPG
        !allocate(elemPG_caf(max_caf_elem_N, Ncol_elemPG)[*], stat=allocation_status, errmsg=emsg)
        !call utility_check_allocation(allocation_status, emsg)
        !elemPG_caf(:,:)[:] = nullvalueI

        !ncol => Ncol_elemSI
        allocate(elemSI_caf(max_caf_elem_N, Ncol_elemSI)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemSI_caf(:,:) = nullvalueI

        !ncol => Ncol_elemSR
        allocate(elemSR_caf(max_caf_elem_N, Ncol_elemSR)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemSR_caf(:,:) = nullvalueR

        !==== face allocation ====
        !ncol => Ncol_faceR 
        allocate(faceR_caf(max_caf_face_N, Ncol_faceR )[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        faceR_caf(:,:) = nullvalueR

        !ncol=> Ncol_faceI
        allocate(faceI_caf(max_caf_face_N, Ncol_faceI)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        faceI_caf(:,:) = nullvalueI

        !ncol=> Ncol_faceYN
        allocate(faceYN_caf(max_caf_face_N, Ncol_faceYN)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        faceYN_caf(:,:) = nullvalueL


        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name


    end subroutine coarray_storage_allocation
end module allocate_storage

module allocate_storage

    use array_index
    use globals
    use utility, only: utility_check_allocation
    use setting_definition, only: setting

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

        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        allocate(nodeI(N_node, ni_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        nodeI(:,:) = nullvalueI

        allocate(linkI(N_link, li_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        linkI(:,:) = nullvalueI

        allocate(P_nodeI(N_node, 3), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        P_nodeI(:,:) = nullvalueI

        allocate(P_linkI(N_link, 2), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        P_linkI(:,:) = nullvalueI

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

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name
    end subroutine allocate_linknode_storage


    subroutine coarray_storage_allocation()
        ! the max_caf_elem and max_caf_face are the maximum length of the coarray 
        ! across all employed images
        ! ==========================
        ! This will be excuted at parallel level
        ! ==========================

        character(64) :: subroutine_name = 'coarray_storage_allocation'
    
        !-----------------------------------------------------------------------------

        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

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


        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name


    end subroutine coarray_storage_allocation
    ! subroutine allocate_bc()
    ! !
    ! ! allocate storage for boundary conditions.
    ! !
    ! !-----------------------------------------------------------------------------

    !     character(64)      :: subroutine_name = 'allocate_bc'
    !     integer            :: ii
    !     integer            :: allocation_status
    !     character(len=99)  :: emsg

    ! !-----------------------------------------------------------------------------

    !     if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

    !     !% the Upstream and Downstream bc structure
    !     allocate(bcdataUp(N_BCupstream), stat=allocation_status, errmsg=emsg)
    !     call utility_check_allocation (allocation_status, emsg)

    !     allocate(bcdataDn(N_BCdnstream), stat=allocation_status, errmsg=emsg)
    !     call utility_check_allocation (allocation_status, emsg)

    !     !% the downstream arrays - HACK default downstream is elevation
    !     do ii=1,N_BCdnstream
    !         bcdataDn(ii)%idx = ii
    !         bcdataDn(ii)%Updn          = bc_updn_downstream
    !         bcdataDn(ii)%Category      = bc_category_elevation
    !         bcdataDn(ii)%NodeID         = nullvalueI
    !         bcdataDn(ii)%FaceID         = nullvalueI
    !         bcdataDn(ii)%ElemGhostID    = nullvalueI
    !         bcdataDn(ii)%ElemInsideID   = nullvalueI
    !         bcdataDn(ii)%ThisValue      = nullvalueR
    !         bcdataDn(ii)%ThisTime       = nullvalueR
    !         bcdataDn(ii)%ThisFlowrate   = nullvalueR
    !         allocate(bcdataUp(ii)%TimeArray(setting%Constant%BCSlots))
    !         allocate(bcdataUp(ii)%ValueArray(setting%Constant%BCSlots))
    !     end do

    !     !% the upstream arrays = HACK default upstream is flowrate
    !     do ii=1,N_BCupstream
    !         bcdataUp(ii)%idx = ii
    !         bcdataUp(ii)%Updn       = bc_updn_upstream
    !         bcdataUp(ii)%category   = bc_category_inflowrate
    !         bcdataUp(ii)%NodeID         = nullvalueI
    !         bcdataUp(ii)%FaceID         = nullvalueI
    !         bcdataUp(ii)%ElemGhostID    = nullvalueI
    !         bcdataUp(ii)%ElemInsideID   = nullvalueI
    !         bcdataUp(ii)%ThisValue      = nullvalueR
    !         bcdataUp(ii)%ThisTime       = nullvalueR
    !         bcdataUp(ii)%ThisFlowrate   = nullvalueR
    !         allocate(bcdataUp(ii)%TimeArray(setting%Constant%BCSlots))
    !         allocate(bcdataUp(ii)%ValueArray(setting%Constant%BCSlots))
    !     end do

    !     if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name
    ! end subroutine allocate_bc

end module allocate_storage

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
    public :: allocate_coarray_storage

contains

    subroutine allocate_all()
        call allocate_linknode_storage()
    end subroutine allocate_all

    subroutine deallocate_all()
        ! from allocate_linknode_storage
        deallocate(nodeI)
        deallocate(linkI)
        deallocate(P_nodeI)
        deallocate(P_linkI)
        deallocate(nodeR)
        deallocate(linkR)
        deallocate(nodeYN)
        deallocate(linkYN)
    end subroutine deallocate_all

    subroutine deallocate_all_temporal()
        ! deallocate(totalInflows%val)
        ! deallocate(tempInflows%val)
    end subroutine deallocate_all_temporal

    subroutine allocate_linknode_storage()
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
        integer       :: additional_rows = 0

    !-----------------------------------------------------------------------------

        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% If BIPquick is being used for Partitioning, include additional rows to the link-node arrays
        if (setting%Partitioning%PartitioningMethod == BQuick) then
            additional_rows = setting%Partitioning%N_Image - 1
        end if

        allocate(nodeI(N_node + additional_rows, ni_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        nodeI(:,:) = nullvalueI

        allocate(linkI(N_link + additional_rows, li_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        linkI(:,:) = nullvalueI

        allocate(P_nodeI(N_node + additional_rows, 3), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        P_nodeI(:,:) = nullvalueI

        allocate(P_linkI(N_link + additional_rows, 2), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        P_linkI(:,:) = nullvalueI

        allocate(nodeR(N_node + additional_rows, nr_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        nodeR(:,:) = nullvalueR

        allocate(linkR(N_link + additional_rows, lr_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        linkR(:,:) = nullvalueR

        allocate(nodeYN(N_node + additional_rows, nYN_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        nodeYN(:,:) = nullvalueL

        allocate(linkYN(N_link + additional_rows, lYN_idx_max), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        linkYN(:,:) = nullvalueL

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name
    end subroutine allocate_linknode_storage
    !
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    !
    subroutine allocate_coarray_storage()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   the max_caf_elem and max_caf_face are the maximum length of the coarray
    !   across all employed images
    !   ==========================
    !   This will be excuted at parallel level
    !   ==========================
    !
    ! Method:
    !
    !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'allocate_coarray_storage'

    !-----------------------------------------------------------------------------

        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        ! --- element arrays allocation ---

        allocate(elemR(max_caf_elem_N, Ncol_elemR)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemR(:,:) = nullvalueR

        allocate(elemI(max_caf_elem_N, Ncol_elemI)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemI(:,:) = nullvalueI

        allocate(elemYN(max_caf_elem_N, Ncol_elemYN)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemYN(:,:) = nullvalueL

        allocate(elemP(max_caf_elem_N, Ncol_elemP)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemP(:,:) = nullvalueI

        ! Ncol_elemPG not defined yet. Comment out
        !ncol => Ncol_elemPG
        !allocate(elemPG_caf(max_caf_elem_N, Ncol_elemPG)[*], stat=allocation_status, errmsg=emsg)
        !call utility_check_allocation(allocation_status, emsg)
        !elemPG_caf(:,:)[:] = nullvalueI

        allocate(elemSI(max_caf_elem_N, Ncol_elemSI)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemSI(:,:) = nullvalueI

        allocate(elemSR(max_caf_elem_N, Ncol_elemSR)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        elemSR(:,:) = nullvalueR

        !==== face allocation ====
        allocate(faceR(max_caf_face_N, Ncol_faceR )[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        faceR(:,:) = nullvalueR

        allocate(faceI(max_caf_face_N, Ncol_faceI)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        faceI(:,:) = nullvalueI

        allocate(faceYN(max_caf_face_N, Ncol_faceYN)[*], stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)
        faceYN(:,:) = nullvalueL

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name


    end subroutine allocate_coarray_storage

    subroutine allocate_temporal_arrays()
        type(LimiterArraySizeType), pointer :: a_size => setting%Limiter%ArraySize

        allocate(tempInflows%val(N_inflow, num_nodeInflow_attributes, a_size%TemporalBC))
        ! 3D array - (node, (time, values))
        allocate(totalInflows%val(N_inflow, 2, a_size%TemporalBC))
    end subroutine allocate_temporal_arrays

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

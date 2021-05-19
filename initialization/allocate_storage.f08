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
        integer       :: additional_rows = 0

    !-----------------------------------------------------------------------------

        if (setting%Debug%File%allocate_storage) print *, '*** enter ',subroutine_name

        !% If BIPquick is being used for Partitioning, include additional rows to the link-node arrays
        if (setting%Partitioning%PartitioningMethod == BIPquick) then
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

        allocate(nodeName(N_node + additional_rows), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)

        allocate(linkName(N_link + additional_rows), stat=allocation_status, errmsg=emsg)
        call utility_check_allocation(allocation_status, emsg)

        if (setting%Debug%File%allocate_storage) print *, '*** leave ',subroutine_name
    end subroutine allocate_linknode_storage

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

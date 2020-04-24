! module allocate_storage
!
! This has all the large array allocation.
! The only other allocation should be in boundary conditions (module bc)
!
!==========================================================================
!
 module allocate_storage
!
! all the top-level storage should be allocated in this module
!
    use array_index
    use globals
    use utility

    implicit none
    private

    public :: allocate_linknode_storage
    public :: allocate_data_storage

    integer, private            :: allocation_status
    character(len=99), private  :: emsg

    integer, private    :: debuglevel = 0

 contains
!
!============================================================================
!============================================================================
!
 subroutine allocate_linknode_storage &
    (linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName)
!
! allocates the link and node storage used for the coarse representation
! of the network connectivity
!
 character(64) :: subroutine_name = 'allocate_linknode_storage'

 integer,   dimension(:,:), allocatable, target, intent(out)    :: linkI
 integer,   dimension(:,:), allocatable, target, intent(out)    :: nodeI

 real,      dimension(:,:), allocatable, target, intent(out)    :: linkR
 real,      dimension(:,:), allocatable, target, intent(out)    :: nodeR

 logical,   dimension(:,:), allocatable, target, intent(out)    :: linkYN
 logical,   dimension(:,:), allocatable, target, intent(out)    :: nodeYN

 type(string), dimension(:), allocatable, target, intent(out)   :: linkName
 type(string), dimension(:), allocatable, target, intent(out)   :: nodeName


!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

 allocate( nodeI(N_node, ni_idx_max), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 nodeI(:,:) = nullvalueI

 allocate( linkI(N_link, li_idx_max), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 linkI(:,:) = nullvalueI

 allocate( nodeR(N_node, nr_idx_max), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 nodeR(:,:) = nullvalueR

 allocate( linkR(N_link, lr_idx_max), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 linkR(:,:) = nullvalueR

 allocate( nodeYN(N_node, nYN_idx_max), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 nodeYN(:,:) = nullvalueL

 allocate( linkYN(N_link, lYN_idx_max), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 linkYN(:,:) = nullvalueL

 allocate( nodeName(N_node), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)

 allocate( linkName(N_link), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine allocate_linknode_storage
!
!==========================================================================
!==========================================================================
!
 subroutine allocate_data_storage &
    (elem2R, elemMR, faceR, elem2I, elemMI, faceI, elem2YN, elemMYN, faceYN, elem2Name, elemMname, faceName)
!
! allocation of the element and face arrays for the high-resolution network
! All storage includes an extra dummy space (e.g. N_elem+1) to allow for
! where statements that remain in bounds during mapping.
!
 character(64) :: subroutine_name = 'allocate_data_storage'

 real,       dimension(:,:), allocatable, target, intent(out)    :: elem2R       ! real data for elements
 integer,    dimension(:,:), allocatable, target, intent(out)    :: elem2I       ! integer data for elements
 logical,    dimension(:,:), allocatable, target, intent(out)    :: elem2YN      ! logical data for elements

 real,       dimension(:,:), allocatable, target, intent(out)    :: elemMR       ! real data for elements
 integer,    dimension(:,:), allocatable, target, intent(out)    :: elemMI       ! integer data for elements
 logical,    dimension(:,:), allocatable, target, intent(out)    :: elemMYN      ! logical data for elements

 real,       dimension(:,:), allocatable, target, intent(out)    :: faceR       ! real data for faces
 integer,    dimension(:,:), allocatable, target, intent(out)    :: faceI       ! integer data for faces
 logical,    dimension(:,:), allocatable, target, intent(out)    :: faceYN      ! logical data for faces

 type(string), dimension(:), allocatable, target, intent(out)    :: elemMName, elem2Name, faceName

!--------------------------------------------------------------------------
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

 allocate( elem2YN(first_elem2_index:first_elem2_index+N_elem2, e2YN_idx_max), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 elem2YN(:,:) = .false.

 allocate( elemMYN(first_elemM_index:first_elemM_index+N_elemM, eMYN_idx_max), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 elemMYN(:,:) = .false.

 allocate( elem2R(first_elem2_index:first_elem2_index+N_elem2, e2r_idx_max), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 elem2R(:,:) = 0.0

 allocate( elemMR(first_elemM_index:first_elemM_index+N_elemM, eMr_idx_max), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 elemMR(:,:) = 0.0

 allocate( elem2I(first_elem2_index:first_elem2_index+N_elem2, e2i_idx_max), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 elem2I(:,:) = 0

 allocate( elemMI(first_elemM_index:first_elemM_index+N_elemM, eMi_idx_max), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 elemMI(:,:) = 0



 allocate( faceYN(first_face_index:first_face_index+N_face, fYN_idx_max), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 faceYN(:,:) = .false.

 allocate( faceR(first_face_index:first_face_index+N_face, fr_idx_max), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 faceR(:,:) = 0.0

 allocate( faceI(first_face_index:first_face_index+N_face, fi_idx_max), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 faceI(:,:) = 0


 allocate( elem2Name(first_elem2_index:first_elem2_index+N_elem2), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)

 allocate( elemMName(first_elemM_index:first_elemM_index+N_elemM), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)

 allocate( faceName(first_face_index:first_face_index+N_face), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)

 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ',subroutine_name
 end subroutine allocate_data_storage
!
!==========================================================================
! END OF MODULE allocate_storage
!==========================================================================
 end module allocate_storage

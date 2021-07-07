! module define_globals
!
! These are global parameters and variables that are not associated with
! the 2D data storage arrays or the setting structure.
!
! The globals for data storage are found in modules array_index and data_keys
! The globals for setting are found in module setting_definition
!
!==========================================================================
module define_globals

    use define_types
    use define_api_keys

    implicit none

    public

    ! rm gr24269 - there is an elementh_length defined in setting%Discretization%NominalElemLength
    ! real(8), parameter :: element_length = 10.0 ! This is a temporary element length

    !% ===================================================================================
    !% ARRAYS
    !% ===================================================================================

    !% Number of maximum branches for a junction
    !integer, parameter :: max_us_branch_per_node = 3
    !integer, parameter :: max_ds_branch_per_node = 3
    integer, parameter :: max_branch_per_node = 6
    real(8) :: branchsign(max_branch_per_node)

    !% Main Arrays - Allocated in allocate_storage.f08
    !%  links are the building blocks from SWMM link-node formulation
    type(LinkArray), target :: link

    !%  nodes are the building blocks from the SWMM link-node formulation
    type(NodeArray), target :: node

    !% Boundary conditions
    type(BCArray), target :: BC

    !%  columns of element and face arrays
    integer, allocatable, target :: col_elemI(:)[:]                                !% columns of elemI array
    integer, allocatable, target :: col_elemP(:)[:],       npack_elemP(:)[:]       !% columns and number of packs for elemP array
    integer, allocatable, target :: col_elemPGalltm(:)[:], npack_elemPGalltm(:)[:] !% columns and number of packs for elemPG array for all tm
    integer, allocatable, target :: col_elemPGac(:)[:],    npack_elemPGac(:)[:]    !% columns and number of packs for elemPG array for ac tm
    integer, allocatable, target :: col_elemPGetm(:)[:],   npack_elemPGetm(:)[:]   !% columns and number of packs for elemPG array for etm
    integer, allocatable, target :: col_elemR(:)[:]                                !% columns of elemR array
    integer, allocatable, target :: col_elemSI(:)[:]                               !% columns of elemSI array
    integer, allocatable, target :: col_elemSR(:)[:]                               !% columns of elemSR array
    integer, allocatable, target :: col_elemSGR(:)[:]                              !% columns of elemSGR array
    integer, allocatable, target :: col_elemWDI(:)[:]                              !% columns of elemWDI array
    integer, allocatable, target :: col_elemWDR(:)[:]                              !% columns of elemWDR array
    integer, allocatable, target :: col_elemYN(:)[:]                               !% columns of elemYN array
    integer, allocatable, target :: col_faceI(:)[:]                                !% columns of faceI array
    integer, allocatable, target :: col_faceM(:)[:]                                !% columns of faceM array
    integer, allocatable, target :: col_faceP(:)[:],       npack_faceP(:)[:]       !% columns and number of packs for faceP array
    integer, allocatable, target :: col_facePS(:)[:],      npack_facePS(:)[:]      !% columns and number of packs for facePS array
    integer, allocatable, target :: col_faceR(:)[:]                                !% columns of faceR array
    integer, allocatable, target :: col_faceYN(:)[:]                               !% columns of faceYN array

    !%  number of dummy row in elem arrays (do not change)
    integer, parameter :: N_dummy_elem = 1

    !%  vector of number of elements and faces across images
    integer, allocatable, target :: N_elem(:)
    integer, allocatable, target :: N_face(:)
    integer, allocatable, target :: N_unique_face(:)

    !%  elems in coarray
    real(8), allocatable, target :: elemR(:,:)[:]       !% coarray for elements
    integer, allocatable, target :: elemI(:,:)[:]       !% coarray for element Interger
    logical, allocatable, target :: elemYN(:,:)[:]      !% coarray for element logical
    integer, allocatable, target :: elemP(:,:)[:]       !% coarray for element pack array

    integer, allocatable, target :: elemPGalltm(:,:)[:] !% coarray for element pack geometry array
    integer, allocatable, target :: elemPGac(:,:)[:]    !% coarray for element pack geometry array
    integer, allocatable, target :: elemPGetm(:,:)[:]   !% coarray for element pack geometry array
    integer, allocatable, target :: elemSI(:,:)[:]      !% coarray for special element Integer
    real(8), allocatable, target :: elemSR(:,:)[:]      !% coarray for special elemen Real
    real(8), allocatable, target :: elemSGR(:,:)[:]     !% coarray for special element geometry Real

    !%  faces in coarray
    real(8), allocatable, target :: faceR(:,:)[:]       !% coarray for faces real data
    integer, allocatable, target :: faceI(:,:)[:]       !% coarray for faces integer data
    logical, allocatable, target :: faceYN(:,:)[:]      !% coarray for faces logical data
    integer, allocatable, target :: faceP(:,:)[:]       !% coarray for faces pack array
    integer, allocatable, target :: facePS(:,:)[:]      !% coarray for shared faces pack array
    logical, allocatable, target :: faceM(:,:)[:]       !% coarray for faces mask array

    !% BIPquick Arrays
    integer, allocatable :: B_nodeI(:,:)
    real(8), allocatable :: B_nodeR(:,:)

    !% Partitioning Module Allocatables - Allocated and Deallocated in execute_partitioning.f08
    integer, allocatable :: adjacent_links(:)
    integer, allocatable :: elem_per_image(:)
    logical, allocatable :: image_full(:)

    !% ===================================================================================
    !% CONSTANTS
    !% ===================================================================================

    !% note that nullvalueI < 0 is required
    integer, parameter :: nullvalueI = 998877
    real(8), parameter :: nullvalueR = 9.98877e16
    logical, parameter :: nullvalueL = .false.
    real(8), parameter :: negoneR = -1.0
    real(8), parameter :: zeroR = 0.0
    real(8), parameter :: oneR = 1.0
    real(8), parameter :: twoR = 2.0
    real(8), parameter :: threeR = 3.0
    real(8), parameter :: fourR = 4.0
    real(8), parameter :: fiveR = 5.0
    real(8), parameter :: sixR = 6.0
    real(8), parameter :: eightR = 8.0
    real(8), parameter :: tenR = 10.0
    real(8), parameter :: pi = 4.d0*datan(1.d0)

    real(8), parameter :: oneeighthR = oneR / eightR
    real(8), parameter :: onesixthR = oneR / sixR
    real(8), parameter :: onefourthR = oneR / fourR
    real(8), parameter :: onethirdR = oneR / threeR
    real(8), parameter :: onehalfR = oneR / twoR
    real(8), parameter :: twothirdR = twoR / threeR
    real(8), parameter :: threefourthR = threeR / fourR
    real(8), parameter :: threehalfR = threeR / twoR
    real(8), parameter :: fourthirdsR = fourR / threeR

    integer, parameter :: zeroI = 0
    integer, parameter :: oneI = 1
    integer, parameter :: twoI = 2
    integer, parameter :: threeI = 3
    integer, parameter :: fourI = 4
    integer, parameter :: fiveI = 5
    integer, parameter :: sixI = 6

    !% Number of objects
    integer :: N_link
    integer :: N_node
    integer :: N_headBC
    integer :: N_flowBC
    integer :: N_nBCup
    integer :: N_nBCdn
    integer :: N_nJm
    integer :: N_nStorage
    integer :: N_nJ2
    integer :: N_diag
    integer :: N_ac
    integer :: N_etm

    !% Number of API parameters
    integer, parameter :: N_api_node_attributes = api_node_overflow
    integer, parameter :: N_api_link_attributes = api_conduit_length
    integer, parameter :: N_api_link_xsect_attributes = api_link_xsect_yBot - N_api_link_attributes
    integer, parameter :: N_api_total_link_attributes = N_api_link_attributes + N_api_link_xsect_attributes

    !% Coarray variables
    integer :: max_caf_elem_N    ! size of all elem array in coarray
    integer :: max_caf_face_N    ! size of all face array in coarray

    !% assign status parameters for nodes
    integer, parameter :: nUnassigned = 998877
    integer, parameter :: nAssigned   = 1
    integer, parameter :: nDeferred   = -1

    ! assign status parameters for links
    integer, parameter :: lUnassigned = 998877
    integer, parameter :: lAssigned   = 1
    integer, parameter :: lDeferred   = -1

    ! Constants for Junction
    integer, parameter :: J_elem_add = 7 ! Supplement elements for junction
    integer, parameter :: J_face_add = 6 ! Supplement faces for junction

    ! default number of elements for different node types
    integer, parameter :: N_elem_nJ2      = 0 ! 2-link nodes are assigned to a single face
    integer, parameter :: N_elem_nJm      = 7 ! M-link nodes are assigned a maximum of 7 elements
    integer, parameter :: N_elem_nStorage = 1 ! Storage nodes are assigned to 1 element
    integer, parameter :: N_elem_nBCdn    = 0 ! Downstream BC nodes are assigned to 0 element (only a face)
    integer, parameter :: N_elem_nBCup    = 0 ! Upstream BC nodes are assigned to 0 element (only a face)

    ! useful shortcuts
    !% NOTE: don't use setting%... structure in define_globals to prevent linking problems
    !rm 20210607 brh real(8), pointer :: dt => setting%time%dt  !% need different Hydrology and Hydraulics dt
    !rm 20210610 brh real(8), pointer :: grav => setting%constant%gravity
    !rm 20210629 griano integer, parameter :: debuglevelall = 0 ! set to 1 to get print of subroutine calls

    !% 20210607 brh Moved these from globals and put in Discretization
    !real(8), pointer :: elem_nominal_length => setting%Discretization%NominalElemLength
    !real(8), pointer :: elem_shorten_cof => setting%Discretization%LinkShortingFactor

    integer, parameter :: NoAdjust       = 1   !% no link length adjustment has done
    integer, parameter :: OneSideAdjust  = 2   !% one sided link length adjustment has done
    integer, parameter :: BothSideAdjust = 3   !% both sided link length adjustment has done

    ! default for edge and non-edge node
    integer, parameter :: nonEdgeNode = 0 ! Upstream BC nodes are assigned to 1 element
    integer, parameter :: EdgeNode    = 1 ! Edge node of a partition

    ! defaults from initial depth type in links
    integer, parameter :: Uniform           = 1
    integer, parameter :: LinearlyVarying   = 2
    integer, parameter :: ExponentialDecay  = 3

end module define_globals

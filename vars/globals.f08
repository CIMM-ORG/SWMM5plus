! module globals
!
! These are global parameters and variables that are not associated with
! the 2D data storage arrays or the setting structure.
!
! The globals for data storage are found in modules array_index and data_keys
! The globals for setting are found in module setting_definition
!
!==========================================================================
module globals

    use setting_definition
    use type_definitions

    implicit none

    public

    ! Main Arrays
    !%  links are the building blocks from SWMM link-node formulation
    real(8), dimension(:,:), allocatable, target :: linkR ! real data for links
    integer, dimension(:,:), allocatable, target :: linkI ! integer data for links
    integer, dimension(:,:), allocatable, target :: P_linkI ! Partitioning output for links
    logical, dimension(:,:), allocatable, target :: linkYN ! logical data for links

    !%  nodes are the building blocks from teh SWMM link-node formulation
    real(8), dimension(:,:), allocatable, target :: nodeR ! real data for nodes
    integer, dimension(:,:), allocatable, target :: nodeI ! integer data for nodes
    integer, dimension(:,:), allocatable, target :: P_nodeI ! Partitioning output for nodes
    logical, dimension(:,:), allocatable, target :: nodeYN ! logical data for nodes

    !% element coarrays
    real(8), allocatable :: elemR(:,:)[:]   ! coarray for elements
    integer, allocatable :: elemI(:,:)[:]    ! coarray for element Interger
    logical, allocatable :: elemYN(:,:)[:]   ! coarray for element logical
    integer, allocatable :: elemP(:,:)[:]    ! coarray for element pack array
    !integer, allocatable, target :: elemPG(:,:)[:]   ! coarray for element pack geometry array   [NOTE] elemPG not defined yet
    integer, allocatable :: elemSI(:,:)[:]   ! coarray for special element Integer
    real(8), allocatable :: elemSR(:,:)[:]   ! coarray for special elemen Real
    real(8), allocatable :: elemSGR(:,:)[:]  ! coarray for special element geometry Real

    !% face coarrays
    real(8), allocatable, target :: faceR(:,:)[:]    ! coarray for faces real data
    integer, allocatable, target :: faceI(:,:)[:]    ! coarray for faces integer data
    logical, allocatable, target :: faceYN(:,:)[:]   ! coarray for faces logical data
    integer, allocatable, target :: faceP(:,:)[:]    ! coarray for faces pack array


    type(string), dimension(:), allocatable, target :: nodeName ! array of character strings

    ! note that nullvalueI < 0 is required
    integer, parameter :: nullvalueI = 998877
    real(8), parameter :: nullvalueR = 9.98877e16
    logical, parameter :: nullvalueL = .false.
    real(8), parameter :: negoneR = -1.0
    real(8), parameter :: zeroR = 0.0
    real(8), parameter :: oneR = 1.0
    real(8), parameter :: twoR = 2.0
    real(8), parameter :: threeR = 3.0
    real(8), parameter :: fourR = 4.0
    real(8), parameter :: sixR = 6.0
    real(8), parameter :: eightR = 8.0
    real(8), parameter :: tenR = 10.0
    real(8), parameter :: pi = 4.d0*datan(1.d0)

    real(8), parameter :: oneeighthR = oneR / eightR
    real(8), parameter :: onefourthR = oneR / fourR
    real(8), parameter :: onethirdR = oneR / threeR
    real(8), parameter :: onehalfR = oneR / twoR
    real(8), parameter :: twothirdR = twoR / threeR
    real(8), parameter :: threefourthR = threeR / fourR

    integer, parameter :: zeroI = 0
    integer, parameter :: oneI = 1
    integer, parameter :: twoI = 2
    integer, parameter :: threeI = 3

    ! Number of objects
    integer :: N_link
    integer :: N_node
    integer :: N_curve
    integer :: N_tseries
    integer :: N_pattern
    integer :: N_inflow
    integer :: N_BCupstream
    integer :: N_BCdnstream
    integer :: N_nBCup
    integer :: N_nBCdn
    integer :: N_nJm
    integer :: N_nStorage
    integer :: N_nJ2

    ! Coarray variables
    integer :: max_caf_elem_N ! size of all elem array in coarray
    integer :: max_caf_face_N ! size of all face array in coarray

    ! Constants for Junction
    integer :: J_elem_add = 7 ! Supplement elements for junction
    integer :: J_face_add = 6 ! Supplement faces for junction

    ! useful shortcuts
    real(8), pointer :: dt => setting%time%dt
    real(8), pointer :: grav => setting%constant%gravity
    integer, parameter :: debuglevelall = 0 ! set to 1 to get print of subroutine calls
    real(8), pointer :: elem_nominal_length => setting%Discretization%NominalElemLength
    real(8), pointer :: elem_branch_factor => setting%Discretization%BranchFactor

    ! 3D Circular Queues
    type(Array3D), target :: tempInflows
    type(Array3D), target :: totalInflows

    ! Boundary Conditions

end module globals
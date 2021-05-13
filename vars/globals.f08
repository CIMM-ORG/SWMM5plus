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

    real(8), parameter :: element_length = 10.0 ! This is a temperory element length

    ! Main Arrays
    !%  links are the building blocks from SWMM link-node formulation
    real(8), dimension(:,:), allocatable, target :: linkR ! real data for links
    integer, dimension(:,:), allocatable, target :: linkI ! integer data for links
    integer, dimension(:,:), allocatable, target :: B_linkI ! BIPquick output for links
    logical, dimension(:,:), allocatable, target :: linkYN ! logical data for links

    type(string), dimension(:), allocatable, target :: linkName ! array of character strings

    !%  nodes are the building blocks from teh SWMM link-node formulation
    real(8), dimension(:,:), allocatable, target :: nodeR ! real data for nodes
    integer, dimension(:,:), allocatable, target :: nodeI ! integer data for nodes
    integer, dimension(:,:), allocatable, target :: B_nodeI ! BIPquick output for nodes
    logical, dimension(:,:), allocatable, target :: nodeYN ! logical data for nodes

    !% elems in coarray
    real(8), allocatable :: elemR_caf(:,:)[:]   ! coarray for elements
    integer, allocatable :: elemI_caf(:,:)[:]    ! coarray for element Interger
    logical, allocatable :: elemYN_caf(:,:)[:]   ! coarray for element logical
    integer, allocatable :: elemP_caf(:,:)[:]    ! coarray for element pack array
    !integer, allocatable, target :: elemPG_caf(:,:)[:]   ! coarray for element pack geometry array   [NOTE] elemPG not defined yet
    integer, allocatable :: elemSI_caf(:,:)[:]   ! coarray for special element Integer
    real(8), allocatable :: elemSR_caf(:,:)[:]   ! coarray for special elemen Real
    real(8), allocatable :: elemSGR_caf(:,:)[:]  ! coarray for special element geometry Real

    !% faces in coarray
    real(8), allocatable, target :: faceR_caf(:,:)[:]    ! coarray for faces real data
    integer, allocatable, target :: faceI_caf(:,:)[:]    ! coarray for faces integer data
    logical, allocatable, target :: faceYN_caf(:,:)[:]   ! coarray for faces logical data
    integer, allocatable, target :: faceP_caf(:,:)[:]    ! coarray for faces pack array
    

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

    ! Coarray variables
    integer :: max_caf_elem_N ! size of all elem array in coarray
    integer :: max_caf_face_N ! size of all face array in coarray

    ! images information
    integer :: nimgs
   

    ! useful shortcuts
    real(8), pointer :: dt => setting%time%dt
    real(8), pointer :: grav => setting%constant%gravity
    integer, parameter :: debuglevelall = 0 ! set to 1 to get print of subroutine calls

    ! Tables
    type(real_table), allocatable :: all_tseries(:)
    type(pattern), allocatable :: all_patterns(:)
    type(totalInflow), allocatable :: total_inflows(:)

end module globals
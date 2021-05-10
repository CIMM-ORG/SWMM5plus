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
    logical, dimension(:,:), allocatable, target :: linkYN ! logical data for links

    type(string), dimension(:), allocatable, target :: linkName ! array of character strings

    !%  nodes are the building blocks from teh SWMM link-node formulation
    real(8), dimension(:,:), allocatable, target :: nodeR ! real data for nodes
    integer, dimension(:,:), allocatable, target :: nodeI ! integer data for nodes
    logical, dimension(:,:), allocatable, target :: nodeYN ! logical data for nodes

    !% elems in coarray
    real(8), dimension(:,:)[:], allocatable, target :: elemR_caf ! coarray for elements
    integer, dimension(:,:)[:], allocatable, target :: elemI_caf ! coarray for element Interger
    logical, dimension(:,:)[:], allocatable, target :: elemYN_caf ! coarray for element logical

    !% faces in coarray
    real(8), dimension(:,:)[:]. allocatable, target :: faceR_caf ! coarray for faces real data
    integer , dimension(:,:)[:], allocatable, target :: faceI_caf ! coarray for faces integer data
    logical, dimension(:,:)[:], allocatable, target :: faceYN_caf ! coarray for faces logical data
    

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

    ! useful shortcuts
    real(8), pointer :: dt => setting%time%dt
    real(8), pointer :: grav => setting%constant%gravity
    integer, parameter :: debuglevelall = 0 ! set to 1 to get print of subroutine calls

    ! Tables
    type(real_table), allocatable :: all_tseries(:)
    type(pattern), allocatable :: all_patterns(:)
    type(totalInflow), allocatable :: total_inflows(:)

end module globals
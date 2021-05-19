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
    integer, dimension(:,:), allocatable, target :: P_linkI ! Partitioning output for links
    logical, dimension(:,:), allocatable, target :: linkYN ! logical data for links

    type(string), dimension(:), allocatable, target :: linkName ! array of character strings

    !%  nodes are the building blocks from teh SWMM link-node formulation
    real(8), dimension(:,:), allocatable, target :: nodeR ! real data for nodes
    integer, dimension(:,:), allocatable, target :: nodeI ! integer data for nodes
    integer, dimension(:,:), allocatable, target :: P_nodeI ! Partitioning output for nodes
    logical, dimension(:,:), allocatable, target :: nodeYN ! logical data for nodes

    !%  columns of element and face arrays
    integer, dimension(:), allocatable, target :: col_elemI[:]                                  ! columns of elemI array
    integer, dimension(:), allocatable, target :: col_elemP[:], npack_elemP[:]                  ! columns and number of packs for elemP array
    integer, dimension(:), allocatable, target :: col_elemPGalltm[:], npack_elemPGalltm[:]      ! columns and number of packs for elemPG array for all tm
    integer, dimension(:), allocatable, target :: col_elemPGac[:], npack_elemPGac[:]            ! columns and number of packs for elemPG array for ac tm
    integer, dimension(:), allocatable, target :: col_elemPGetm[:], npack_elemPGetm[:]          ! columns and number of packs for elemPG array for etm
    integer, dimension(:), allocatable, target :: col_elemR[:]                                  ! columns of elemR array
    integer, dimension(:), allocatable, target :: col_elemSI[:]                                 ! columns of elemSI array
    integer, dimension(:), allocatable, target :: col_elemSR[:]                                 ! columns of elemSR array
    integer, dimension(:), allocatable, target :: col_elemSGR[:]                                ! columns of elemSGR array
    integer, dimension(:), allocatable, target :: col_elemWDI[:]                                ! columns of elemWDI array
    integer, dimension(:), allocatable, target :: col_elemWDR[:]                                ! columns of elemWDR array
    integer, dimension(:), allocatable, target :: col_elemYN[:]                                 ! columns of elemYN array
    integer, dimension(:), allocatable, target :: col_faceI[:]                                  ! columns of faceI array
    integer, dimension(:), allocatable, target :: col_faceM[:]                                  ! columns of faceM array
    integer, dimension(:), allocatable, target :: col_faceP[:], npack_faceP[:]                  ! columns and number of packs for faceP array
    integer, dimension(:), allocatable, target :: col_faceR[:]                                  ! columns of faceR array
    integer, dimension(:), allocatable, target :: col_faceYN[:]                                 ! columns of faceYN array

    !%  vector of number of elements and faces across images
    integer, dimension(:), allocatable, target :: N_elem
    integer, dimension(:), allocatable, target :: N_face 

    !%  elems in coarray
    real(8), allocatable, target :: elemR(:,:)[:]    ! coarray for elements
    integer, allocatable, target :: elemI(:,:)[:]    ! coarray for element Interger
    logical, allocatable, target :: elemYN(:,:)[:]   ! coarray for element logical
    integer, allocatable, target :: elemP(:,:)[:]    ! coarray for element pack array
    integer, allocatable, target :: elemPG(:,:)[:]   ! coarray for element pack geometry array   [NOTE] elemPG not defined yet
    integer, allocatable, target :: elemSI(:,:)[:]   ! coarray for special element Integer
    real(8), allocatable, target :: elemSR(:,:)[:]   ! coarray for special elemen Real
    real(8), allocatable, target :: elemSGR(:,:)[:]  ! coarray for special element geometry Real

    !%  faces in coarray
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

    ! images control
    integer :: image
    ! Number of objects
    integer :: N_link
    integer :: N_node
    integer :: N_curve
    integer :: N_tseries
    integer :: N_pattern
    integer :: N_BCupstream
    integer :: N_BCdnstream

    ! Coarray variables
    integer :: max_caf_elem_N ! size of all elem array in coarray
    integer :: max_caf_face_N ! size of all face array in coarray


    ! Constants for Junction
    integer, target :: J_elem_add = 7 ! Supplement elements for junction
    integer, target :: J_face_add = 6 ! Supplement faces    for junction
   

    ! useful shortcuts
    real(8), pointer :: dt => setting%time%dt
    real(8), pointer :: grav => setting%constant%gravity
    integer, parameter :: debuglevelall = 0 ! set to 1 to get print of subroutine calls
    real(8), pointer :: elem_shorten_cof => setting%ElementLengthAdjust%LinkShortingFactor

    ! Tables
    type(real_table), allocatable :: all_tseries(:)
    type(pattern), allocatable :: all_patterns(:)
    type(totalInflow), allocatable :: total_inflows(:,:,:)

    ! Boundary Conditions
    real(8), allocatable :: bcdataDn
    real(8), allocatable :: bcdataUp

end module globals
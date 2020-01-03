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

    character(len=99) :: casename

    ! note that nullvalueI < 0 is required
    integer, parameter :: nullvalueI = -998877
    real,    parameter :: nullvalueR = -9.98877e16
    logical, parameter :: nullvalueL = .false.
    real,    parameter :: zeroR      = 0.0
    real,    parameter :: oneR       = 1.0
    real,    parameter :: twoR       = 2.0
    real,    parameter :: threeR     = 3.0
    real,    parameter :: fourR      = 4.0
    real,    parameter :: six        = 6.0
    real,    parameter :: eightR     = 8.0
    real,    parameter :: tenR       = 10.0
    real,    parameter :: pi         = 4.d0*datan(1.d0)


    real,    parameter :: oneeighthR    = oneR   / eightR
    real,    parameter :: onefourthR   = oneR   / fourR
    real,    parameter :: onethirdR    = oneR   / threeR
    real,    parameter :: onehalfR     = oneR   / twoR
    real,    parameter :: twothirdR    = twoR   / threeR
    real,    parameter :: threefourthR = threeR / fourR

    integer, parameter :: zeroI      = 0
    integer, parameter :: oneI       = 1

    integer :: N_link
    integer :: N_node
    integer :: N_elem2
    integer :: N_elemM
    integer :: N_face
    integer :: N_BCupstream
    integer :: N_BCdnstream
    integer :: dummy_face_index
    integer :: dummy_elem2_index
    integer :: dummy_elemM_index


    integer :: next_e2i_temparray = 1
    integer :: next_e2r_temparray = 1
    integer :: next_e2YN_temparray = 1

    integer :: next_eMi_temparray = 1
    integer :: next_eMr_temparray = 1
    integer :: next_eMYN_temparray = 1

    integer :: next_fi_temparray = 1
    integer :: next_fr_temparray = 1
    integer :: next_fYN_temparray = 1

    integer :: outputfile_next_unitnumber = 10 ! used for fileopening

    ! useful shortcuts
    real, pointer :: dt   => setting%time%dt
    real, pointer :: grav => setting%constant%gravity

    integer :: debugcounter = 0

    integer, parameter :: debuglevelall = 0 ! set to 1 to get print of subroutine calls

    integer :: idummyarray(1)


!==========================================================================
! END OF MODULE globals
!==========================================================================
 end module globals

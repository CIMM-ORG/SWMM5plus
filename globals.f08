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
    integer,    parameter :: nullvalueI = -998877
    real(8),    parameter :: nullvalueR = -9.98877e16
    logical,    parameter :: nullvalueL = .false.
    real(8),    parameter :: negoneR    = -1.0
    real(8),    parameter :: zeroR      = 0.0
    real(8),    parameter :: oneR       = 1.0
    real(8),    parameter :: twoR       = 2.0
    real(8),    parameter :: threeR     = 3.0
    real(8),    parameter :: fourR      = 4.0
    real(8),    parameter :: sixR       = 6.0

    real(8),    parameter :: eightR     = 8.0
    real(8),    parameter :: tenR       = 10.0
    real(8),    parameter :: pi         = 4.d0*datan(1.d0)


    real(8),    parameter :: oneeighthR    = oneR   / eightR
    real(8),    parameter :: onefourthR   = oneR   / fourR
    real(8),    parameter :: onethirdR    = oneR   / threeR
    real(8),    parameter :: onehalfR     = oneR   / twoR
    real(8),    parameter :: twothirdR    = twoR   / threeR
    real(8),    parameter :: threefourthR = threeR / fourR

    integer, parameter :: zeroI      = 0
    integer, parameter :: oneI       = 1
    integer, parameter :: twoI       = 2
    integer, parameter :: threeI     = 3

    integer :: N_link
    integer :: N_node
    integer :: N_elem2
    integer :: N_elemM
    integer :: N_face
    integer :: N_BCupstream
    integer :: N_BCdnstream
    integer :: N_Gates
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
    real(8), pointer :: dt   => setting%time%dt
    real(8), pointer :: grav => setting%constant%gravity

    type(graph) :: swmm_graph
    integer :: debugcounter = 0

    integer, parameter :: debuglevelall = 0 ! set to 1 to get print of subroutine calls

    integer :: idummyarray(1)

    !==========================================================================
    ! END OF MODULE globals
    !==========================================================================
end module globals

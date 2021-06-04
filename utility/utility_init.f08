module utility_init

    use define_indexes
    use define_globals
    use define_settings, only: setting
    use utility, only: utility_check_allocation

    implicit none

!-----------------------------------------------------------------------------
!
! Description:
!   This module contains all of the utility functions that are used in the 
!   initialization arm of the code. 
!
!-----------------------------------------------------------------------------

    private

    ! utility_allocate constants
    integer           :: allocation_status
    character(len=99) ::              emsg




contains
